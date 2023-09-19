/*
 *  Copyright 2022-2023 LAIKA. Authored by Mitch Prater.
 *
 *  Licensed under the Apache License Version 2.0 http://apache.org/licenses/LICENSE-2.0,
 *  or the MIT license http://opensource.org/licenses/MIT, at your option.
 *
 *  This program may not be copied, modified, or distributed except according to those terms.
 */

#include "RixRNG.h"
#include "RixBxdf.h"
#include "RixBxdfLobe.h"
#include "RixShading.h"
#include "RixShadingUtils.h"
#include "RixPredefinedStrings.hpp"
#include "bxdf/PxrSurfaceOpacity.h"

// Provides a non-const pointer to a const data set so it can be written to.
#define NON_CONST_PTR(type,ptr) const_cast< type >( static_cast< const type >( ptr ))


/*
==================================================================
The RixBxdfFactory class: the bxdf shader plugin.
This contains data and methods for defining the interactions of
the plugin with the renderer and integrator, including evaluation
of its user parameters and their shading network connections.
==================================================================
*/
class BxdfFactory : public RixBxdfFactory
{
  private:

    // Any general-purpose Rix Interface handles we might
    // need that are accessible within the Init() method.
    RixMessages*  rixMsg;

    // Data struct that contains plugin parameter values and other
    // data that can be computed once for each unique invocation
    // (a.k.a. instance) of the plugin.
    struct pluginInstanceData
    {
        // Stores the bitfield value that will be returned by
        // GetInstanceHints(), which informs the integrator what
        // types of presence/opacity/interior computations this
        // bxdf plugin performs.
        int  poiHints;

        // Frees the struct's memory.
        static void Delete( void* data )
        {
            auto  pluginData = static_cast< pluginInstanceData* >( data );
            delete  pluginData;
        }
    };

    // Default parameter values (from the .args file).
    const float  def_Gain = 1.0f;
    const RtColorRGB  def_Color = RtColorRGB( 0.5f );
    const float  def_Direction = 0.5f;
    const float  def_Dispersion = 0.5f;
    const float  def_Presence = 1.0f;
    const int    def_Socket = 0;

    // pTable indices.
    enum pTableEnries
    {
        in_Gain = 0,
        in_Color,
        in_Direction,
        in_Dispersion,
        in_Presence,
        in_Socket
    };

  public:

    BxdfFactory() {}
    ~BxdfFactory() {}

    // The parameter table specifies the parameter name and type of any
    // parameters the plugin wants to make use of. The parameters and
    // their UI are defined in the plugin's corresponding .args file.
    const RixSCParamInfo* GetParamTable()
    {
        static RixSCParamInfo pTable[] =
        {
            RixSCParamInfo( RtUString( "Gain" ), k_RixSCFloat, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Color" ), k_RixSCColor, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Direction" ), k_RixSCFloat, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Dispersion" ), k_RixSCFloat, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Presence" ), k_RixSCFloat, k_RixSCPresenceInput ),
            RixSCParamInfo( RtUString( "Socket" ), k_RixSCInteger ),
            RixSCParamInfo() // Ends the table.
        };

        return &pTable[0];
    }

    // Get any general-purpose RixInterface handles we need.
    int Init( RixContext& ctx, const RtUString plugPathName )
    {
        rixMsg = static_cast< RixMessages* >( ctx.GetRixInterface( k_RixMessages ));
        return !rixMsg ? -1 : 0;
    }
    // Since Init() didn't allocate any memory, Finalize() is a no-op.
    void Finalize( RixContext& ctx ) {}

    // Sets this bxdf's characteristic Light Path Expression variables at RenderBegin.
    void Synchronize( RixContext&, RixSCSyncMsg, const RixParameterList* );

    // CreateInstanceData() is called once for each unique set of plugin parameter values
    // (a.k.a. a plugin instance) and is used to store those values and other data that
    // can be computed once.
    void CreateInstanceData
    (
        RixContext&             ctx,
        const RtUString         plugNodeName,
        const RixParameterList* pList,
        InstanceData*           instanceData
    )
    {
        // Initialize the InstanceData struct.
        instanceData->datalen = 0;
        instanceData->data = NULL;
        instanceData->freefunc = NULL;

        // Allocate new memory for this plugin's per-instance data.
        auto pluginData = new BxdfFactory::pluginInstanceData;
        if( !pluginData ) return;

        // Set this plugin's (presence/opacity/interior) hints
        // based on how the Presence parameter is being set.
        // The emun InstanceHints entries are defined in RixBxdf.h
        RixSCType            type;
        RixSCConnectionInfo  cinfo;

        pList->GetParamInfo( in_Presence, &type, &cinfo );
        switch( cinfo )
        {
            case k_RixSCDefaultValue:
            {
                pluginData->poiHints = k_TriviallyOpaque;
            }
            case k_RixSCParameterListValue:
            {
                float  Presence( def_Presence );
                pList->EvalParam( in_Presence, -1, &Presence );
                pluginData->poiHints = ( Presence != def_Presence ) ? k_ComputesPresence : k_TriviallyOpaque;
            }
            case k_RixSCNetworkValue:
            {
                // Presumes any connection to Presence contains values < 1.
                pluginData->poiHints = k_ComputesPresence;
            }
        }

        // Set the InstanceData struct members to access this plugin's new per-instance data.
        instanceData->datalen = sizeof( *pluginData );
        instanceData->data = static_cast< void* >( pluginData );
        instanceData->freefunc = BxdfFactory::pluginInstanceData::Delete;
    }

    // GetInstanceHints() provides information to the integrator
    // about this bxdf's presence/opacity/interior computations.
    int GetInstanceHints( void* data ) const
    {
        if( data ) return static_cast< BxdfFactory::pluginInstanceData* >( data )->poiHints;
        else return k_TriviallyOpaque;
    }

    //  Unused. Using CreateInstanceData() instead.
    void SynchronizeInstanceData( RixContext&, const RtUString, const RixParameterList*, const uint32_t, InstanceData* ) {}

    // No useful comments about this in RixBxdf.h, and all examples just return 1.
    float GetIndexOfRefraction( void* data ) const { return 1.0f; }

    // Only used by volume bxdf's that have time-varying data. See PxrVolume.cpp.
    void  RegisterTemporalVolumeParams( void*, std::vector< int >& ) const {}

    //------------------------------------------------
    // Begin/End methods provide the integrator with
    // access to various characteristics of the bxdf.
    //------------------------------------------------

    // BeginScatter() returns a RixBxdf object that encapsulates this bxdf's light-scattering behavior.
    // EndScatter() is called by RixBxdf::Release() to release the object created by BeingScatter().
    RixBxdf* BeginScatter( const RixShadingContext*, const RixBXLobeTraits&, RixSCShadingMode, void*, void* );
    void     EndScatter( RixBxdf* );

    // BeginOpacity() returns a RixOpacity object that encapsulates this bxdf's presence and/or opacity behavior.
    // EndOpacity() is called by RixOpacity::Release() to release the object created by BeginOpacity().
    RixOpacity* BeginOpacity( const RixShadingContext*, RixSCShadingMode, void*, void* );
    void        EndOpacity( RixOpacity* );

    // These functional blocks are not used in this bxdf.
    RixVolumeIntegrator* BeginInterior( const RixShadingContext*, RixSCShadingMode, void*, void* ) { return NULL; }
    void                 EndInterior( RixVolumeIntegrator* ) {}
    RixVolumeIntegrator* BeginSubsurface( const RixShadingContext*, RixSCShadingMode, void*, void* ) { return NULL; }
    void                 EndSubsurface( RixVolumeIntegrator* ) {}
    RixPostLighting*     BeginPostLighting( const RixShadingContext*, RixSCShadingMode, void*, void* ) { return NULL; }
    void                 EndPostLighting( RixPostLighting* ) {}
};


/*
=============================================================================
These static global variables contain the Light Path Expression (LPE) traits
of this bxdf's response(s). Since these values must be communicated across
the RixBxdfFactory and RixBxdf class boundaries, they are declared as static
global variables of this plugin.
=============================================================================
*/
// Declare one RixBXLobeSampled and RixBXLobeTraits pair of variables per response.
// This data type contains information about one response only.
static RixBXLobeSampled  sg_PraterScatter_LS;
// This data type can contain any number of responses. Used to create/define sets of responses.
static RixBXLobeTraits   sg_PraterScatter_LT;

// Synchronize() sets the Lobe variables.
// This can't be done statically since RixBXLookupLobeByName() requires the
// Light Path Expression system which is not available until k_RixSCRenderBegin.
void BxdfFactory::Synchronize( RixContext& ctx, RixSCSyncMsg syncMsg, const RixParameterList* pList )
{
    if( syncMsg != k_RixSCRenderBegin ) return;

    // Query the Light Path Expression (LPE) entry for the "PraterScatter"
    // response (specified in the rendermn.ini file) and set its traits.
    // Each response will require its own static global variables and
    // RixBXLookupLobeByName() and RixBXLobeTraits() calls to set them.
    // rendermn.ini entry: /prman/lpe/specular6  PraterScatter
    // LPE: color lpe:CS6.*[<L.>O]
    sg_PraterScatter_LS = RixBXLookupLobeByName( ctx,
                        false, // not discrete ⇒ samples over a solid angle.
                        true,  // specular ⇒ not diffuse.
                        true,  // reflected scattering ⇒ not transmitted.
                        false, // not a user response ⇒ standard (specular or diffuse).
                        0,     // response "id" number. not used by prman or LPE.
                        "PraterScatter" // The response's Name (used in rendermn.ini).
                        );

    sg_PraterScatter_LT = RixBXLobeTraits( sg_PraterScatter_LS );
}


/* 
=============================================================
The RixBxdf class defines this bxdf's response interactions
with individual rays as its part of the integration process:
this plugin's bxdf response closure. Its GenerateSamples(),
EvaluateSamples(), EvaluateSamplesAtIndex(), and EmitLocal()
methods define those interactions.
=============================================================
*/
class BxdfClosure : public RixBxdf
{
  private:

    // Contains the intersection of the response (lobes) the
    // integrator wants with the lobes this bxdf can produce.
    RixBXLobeTraits   bxdfLobes;

    // Shading context data needed to compute this bxdf.
    // These will consist of numPts values: one for each shaded point.
    // Note: a shading context is also known as a rendering "grid".
    const int  numPts;
    const float*      Gain;
    const RtColorRGB* Color;
    const float*      Direction;
    const float*      Dispersion;
    const RtNormal3*  Ng;
    const RtVector3*  Vn;

  public:

    // The RixBxdf object constructor.
    // This is instantiated by the RixBxdfFactory::BeginScatter()
    // method, which also provides it with any data it will need.
    BxdfClosure
    (
        RixBxdfFactory*          bFac,
        const RixShadingContext* sCtx,
        const RixBXLobeTraits&   lobesWanted, // by the integrator, per bxdf closure.

        // Parameters containing data needed to compute this bxdf's response(s).
        const int         _numPts,
        const float*      _Gain,
        const RtColorRGB* _Color,
        const float*      _Direction,
        const float*      _Dispersion
    ):
        // Initializes the protected RixBxdf class members
        // 'shadingCtx' (shading context) to sCtx, and 'bxdfFactory' to bFac.
        RixBxdf( sCtx, bFac ),

        // Initialize this BxdfClosure's bxdfLobes member variable
        // to the responses the integrator wants from it.
        bxdfLobes( lobesWanted ),

        // Initialize this BxdfClosure's parameter data pointers.
        numPts( _numPts ),
        Gain( _Gain ),
        Color( _Color ),
        Direction( _Direction ),
        Dispersion( _Dispersion )
    {
        // Intersect the response(s) wanted by the integrator (lobesWanted)
        // with the response(s) this bxdf produces. The result defines the set
        // of responses we need to compute in the BxdfClosure::*Sample() methods.
        bxdfLobes &= sg_PraterScatter_LT; // Additional response sg_*_LT values are | together.

        // Save some shading context data in the BxdfClosure's member
        // variables that we'll need later in its *Sample() methods.
        sCtx->GetBuiltinVar( RixShadingContext::k_Ngn, &Ng ); // wg
        sCtx->GetBuiltinVar( RixShadingContext::k_Vn,  &Vn ); // wo
    }
    // Destructor.
    ~BxdfClosure() {}

    // Provides from what direction(s) (a.k.a. domains) around the shaded
    // point rays can come to which this bxdf might react. Ensures that
    // the integrator doesn't bother sending it rays from other directions.
    RixBXEvaluateDomain GetEvaluateDomain() { return k_RixBXOutsideReflect; }

    // Returns a bit field containing all the responses (lobes) this bxdf
    // is capable of producing intersected with those the integrator wants.
    void GetAggregateLobeTraits( RixBXLobeTraits *t ) { *t = bxdfLobes; }

    // If this bxdf has a MaterialIor or Albedo property, return that here.
    // Otherwise, do what all example plugins do and return none: k_RixSCInvalidDetail.
    RixSCDetail GetProperty( BxdfProperty, const void** ) const { return k_RixSCInvalidDetail; }

    // Inline the implementation of this plugin's BxdfClosure methods.
    // These define this bxdf's response interactions with the integrator.
    #include "bxdf/PraterScatterSampling.inl"
};


/*
===============================================================================
The renderer's operating model is that the RixBxdf is a closure with functions
GenerateSample(), EvaluateSample(), EvaluateSamplesAtIndex(), and EmitLocal().
The integrator invokes BeginScatter() to obtain a BxdfClosure object.
The RixBxdfFactory should stash any necessary state within the RixBxdf object,
and consider that the RixBxdf lifetime is under control of the integrator. 
Such state includes any needed parameter values or built-in variables.
===============================================================================
*/
RixBxdf* BxdfFactory::BeginScatter
(
    const RixShadingContext* sCtx,
    const RixBXLobeTraits&   lobesWanted, // by the integrator, per bxdf closure.
    RixSCShadingMode         shadingMode,
    void*                    parentData, // See RixBxdf.h
    void*                    data // The per-instance data.
)
{
    int  numPts = sCtx->numPts;

    // This plugin's per-instance data.
    // auto  pluginData = static_cast< BxdfFactory::pluginInstanceData* >( data );

    // Evaluate the Socket parameter to trigger connected node evaluation.
    const int*  Socket;
    sCtx->EvalParam( in_Socket, -1, &Socket );

    // Evaluate (the potentially varying) parameters.
    const float*      Gain;
    const RtColorRGB* Color;
    const float*      Direction;
    const float*      Dispersion;

    sCtx->EvalParam( in_Gain, -1, &Gain, &def_Gain, true );
    sCtx->EvalParam( in_Color, -1, &Color, &def_Color, true );
    sCtx->EvalParam( in_Direction, -1, &Direction, &def_Direction, true );
    sCtx->EvalParam( in_Dispersion, -1, &Dispersion, &def_Dispersion, true );

    // Convert Direction and Dispersion parameter values (in-place) to g and f.
    float*  writeDirection = NON_CONST_PTR( float*, Direction );
    float*  writeDispersion = NON_CONST_PTR( float*, Dispersion );
    for( int i=0; i < numPts; i++ )
    {
        writeDirection[i] = 0.99f * Direction[i]; // -1 < g < +1
        writeDispersion[i] = 1.0f - Dispersion[i]; // 0 ≤ f ≤ 1
    }


    // Create a shading context memory pool.
    RixShadingContext::Allocator  pool( sCtx );

    // Allocate memory for this shader's RixBxdf object.
    void*  mem = pool.AllocForBxdf< BxdfClosure >(1);

    // Create an instance of this shader's bxdf closure
    // and pass any necessary data to it.
    BxdfClosure*  bxdf = new (mem) BxdfClosure( this, sCtx, lobesWanted,
                                        numPts,
                                        Gain,
                                        Color,
                                        Direction,
                                        Dispersion
                                        );

    return bxdf;
}
// Releases the RixBxdf object created by BeingScatter().
void BxdfFactory::EndScatter( RixBxdf* ) {}


/*
==============================================
Returns a RixOpacity object that encapsulates
this bxdf's presence and/or opacity behavior.
==============================================
*/
RixOpacity* BxdfFactory::BeginOpacity
(
    const RixShadingContext* sCtx,
    RixSCShadingMode         shadingMode,
    void*                    parentData, // See RixBxdf.h
    void*                    data // The per-instance data.
)
{
    // This plugin's per-instance data.
    auto  pluginData = static_cast< BxdfFactory::pluginInstanceData* >( data );

    // Evaluate the Presence parameter and determine whether it is varying or uniform.
    const float*  Presence = NULL;
    bool  uniformPresence = false;

    if( pluginData->poiHints & k_ComputesPresence )
    {
        RixSCDetail  detailPresence = sCtx->EvalParam( in_Presence, -1, &Presence, &def_Presence );
        uniformPresence = k_RixSCUniform == detailPresence;
    }

    if( Presence )
    {
        RixShadingContext::Allocator  pool( sCtx );
        void*  mem = pool.AllocForBxdf< PxrSurfaceOpacity >( 1 );
        return  new (mem) PxrSurfaceOpacity( sCtx, this, Presence, NULL, uniformPresence );
    }
    else return NULL;
}
// Releases the RixOpacity object created by BeingOpacity().
void BxdfFactory::EndOpacity( RixOpacity* ) {}


/*
==============================================
Entrypoints to this plugin from the renderer.
==============================================
*/
extern "C" PRMANEXPORT RixBxdfFactory* CreateRixBxdfFactory( RtConstString )
{
    return new BxdfFactory();
}

extern "C" PRMANEXPORT void DestroyRixBxdfFactory( RixBxdfFactory* bFac )
{
    delete static_cast< BxdfFactory* >( bFac );
}

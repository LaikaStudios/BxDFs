/*
 *  Copyright 2022 LAIKA. Authored by Mitch J Prater.
 *
 *  Licensed under the Apache License Version 2.0 http://apache.org/licenses/LICENSE-2.0,
 *  or the MIT license http://opensource.org/licenses/MIT, at your option.
 *
 *  This program may not be copied, modified, or distributed except according to those terms.
 */
/*  Implements:
 *
 *      2003 Koenderink-Pont "asperity" response.
 *      The secret of velvety skin.
 *      Koenderink, Jan and Pont, Sylvia
 *      Machine Vision and Applications, Vol. 14, p. 260-268 (2003)
 *      http://dx.doi.org/10.1007/s00138-002-0089-7
 */

#include "RixRNG.h"
#include "RixBxdf.h"
#include "RixBxdfLobe.h"
#include "RixShading.h"
#include "RixShadingUtils.h"
#include "RixPredefinedStrings.hpp"


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
    friend class BxdfClosure;

  private:

    // Data struct that contains plugin parameter values and other
    // data that can be computed once for each unique invocation
    // (a.k.a. instance) of the plugin.
    struct pluginInstanceData
    {
        // Orientation parameter connection info.
        RixSCConnectionInfo  cinfoOrientation;

        // Frees the struct's memory.
        static void Delete( void* data )
        {
            auto  pluginData = static_cast< pluginInstanceData* >( data );
            delete  pluginData;
        }
    };

    // Default parameter values (from the .args file).
    const float      def_Gain = 1.0f;
    const RtColorRGB def_Color = RtColorRGB( 0.5f );
    const RtNormal3  def_Orientation = RtNormal3( 0.0f );
    const float      def_Length = 1.0f;

    // pTable indices.
    enum pTableEntries
    {
        in_Gain = 0,
        in_Color,
        in_Orientation,
        in_Length
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
            RixSCParamInfo( RtUString( "Orientation" ), k_RixSCNormal, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Length" ), k_RixSCFloat, k_RixSCScatterInput ),
            RixSCParamInfo() // Ends the table.
        };

        return &pTable[0];
    }

    // Sets this bxdf's characteristic Light Path Expression variables at RenderBegin.
    void Synchronize( RixContext&, RixSCSyncMsg, const RixParameterList* );

    // CreateInstanceData() is called once for each unique set of plugin parameter values
    // (a.k.a. a plugin instance) and is used to store those values and other data that
    // can be computed once per instance.
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

        // Get the Orientation parameter's connection info.
        RixSCType  type;
        pList->GetParamInfo( in_Orientation, &type, &(pluginData->cinfoOrientation) );

        // Set the InstanceData struct members to access this plugin's new per-instance data.
        instanceData->datalen = sizeof( *pluginData );
        instanceData->data = static_cast< void* >( pluginData );
        instanceData->freefunc = BxdfFactory::pluginInstanceData::Delete;
    }

    // GetInstanceHints() provides information to the integrator
    // about this bxdf's presence/opacity/interior computations.
    int GetInstanceHints( void* data ) const { return k_TriviallyOpaque; }

    // No useful comments about this in RixBxdf.h, and all examples just return 1.
    float GetIndexOfRefraction( void* data ) const { return 1.0f; }

    // Unused here.
    int  Init( RixContext&, const RtUString ) { return 0; }
    void Finalize( RixContext& ctx ) {}
    void SynchronizeInstanceData( RixContext&, const RtUString, const RixParameterList*, const uint32_t, InstanceData* ) {}
    void RegisterTemporalVolumeParams( void*, std::vector< int >& ) const {}

    //------------------------------------------------
    // Begin/End methods provide the integrator with
    // access to various characteristics of the bxdf.
    //------------------------------------------------

    // BeginScatter() returns a RixBxdf object that encapsulates this bxdf's light-scattering behavior.
    // EndScatter() is called by RixBxdf::Release() to release the object created by BeginScatter().
    RixBxdf* BeginScatter( const RixShadingContext*, const RixBXLobeTraits&, RixSCShadingMode, void*, void* );
    void     EndScatter( RixBxdf* );

    // These functional blocks are not used in this bxdf.
    RixOpacity*          BeginOpacity( const RixShadingContext*, RixSCShadingMode, void*, void* ) { return NULL; }
    void                 EndOpacity( RixOpacity* ) {}
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
static RixBXLobeSampled  sg_KPAsperity_LS;
// This data type can contain any number of responses. Used to create/define sets of responses.
static RixBXLobeTraits   sg_KPAsperity_LT;

// Synchronize() sets the Lobe variables.
// This can't be done statically since RixBXLookupLobeByName() requires the
// Light Path Expression system which is not available until k_RixSCRenderBegin.
void BxdfFactory::Synchronize( RixContext& ctx, RixSCSyncMsg syncMsg, const RixParameterList* pList )
{
    if( syncMsg != k_RixSCRenderBegin ) return;

    // Query the Light Path Expression (LPE) entry for the "KPAsperity"
    // response (specified in the rendermn.ini file) and set its traits.
    // Each response will require its own static global variables and
    // RixBXLookupLobeByName() and RixBXLobeTraits() calls to set them.
    // rendermn.ini entry: /prman/lpe/specular6  KPAsperity
    // LPE: color lpe:CS6.*[<L.>O]
    sg_KPAsperity_LS = RixBXLookupLobeByName( ctx,
                        false, // not discrete ⇒ samples over a solid angle.
                        true,  // specular ⇒ not diffuse.
                        true,  // reflected scattering ⇒ not transmitted.
                        false, // not a user response ⇒ standard (specular or diffuse).
                        0,     // response "id" number. not used by prman or LPE.
                        "KPAsperity" // The response's Name (used in rendermn.ini).
                        );

    sg_KPAsperity_LT = RixBXLobeTraits( sg_KPAsperity_LS );
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
    const int         numPts;
    const float*      Gain;
    const RtColorRGB* Color;
    const RtNormal3*  Orientation;
    const float*      Length;
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
        const RtNormal3*  _Orientation,
        const float*      _Length,
        const RtVector3*  _Vn
    ):
        // Initializes the protected RixBxdf class members
        // 'shadingCtx' (shading context) to sCtx, and 'bxdfFactory' to bFac.
        RixBxdf( sCtx, bFac ),

        // Initialize this BxdfClosure's bxdfLobes member variable
        // to the responses the integrator wants from it.
        bxdfLobes( lobesWanted ),

        // Initialize this BxdfClosure's member variables.
        numPts( _numPts ),
        Gain( _Gain ),
        Color( _Color ),
        Orientation( _Orientation ),
        Length( _Length ),
        Vn( _Vn )
    {
        // Intersect the response(s) wanted by the integrator (lobesWanted)
        // with the response(s) this bxdf produces. The result defines the set
        // of responses we need to compute in the BxdfClosure::*Sample() methods.
        bxdfLobes &= sg_KPAsperity_LT; // Additional response sg_*_LT values are | together.
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
    #include "bxdf/KPAsperitySampling.inl"
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
    const int  numPts = sCtx->numPts;

    // This plugin's per-instance data.
    auto  pluginData = static_cast< BxdfFactory::pluginInstanceData* >( data );

    // Evaluate (the potentially varying) parameters.
    const float*      Gain;
    const RtColorRGB* Color;
    const RtNormal3*  Orientation;
    const float*      Length;

    sCtx->EvalParam( in_Gain, -1, &Gain, &def_Gain, true );
    sCtx->EvalParam( in_Color, -1, &Color, &def_Color, true );
    sCtx->EvalParam( in_Orientation, -1, &Orientation, &def_Orientation, true );
    sCtx->EvalParam( in_Length, -1, &Length, &def_Length, true );

    // Get some shading context variable pointers.
    const RtNormal3*  Nn;
    const RtNormal3*  Ng;
    const RtVector3*  Vn;
    sCtx->GetBuiltinVar( RixShadingContext::k_Nn,  &Nn );
    sCtx->GetBuiltinVar( RixShadingContext::k_Ngn, &Ng ); // wg
    sCtx->GetBuiltinVar( RixShadingContext::k_Vn,  &Vn ); // wo

    // Orientation -> Nn if it hasn't been set to anything.
    if( k_RixSCDefaultValue == pluginData->cinfoOrientation )
    {
        Orientation = Nn;
    }
    // Otherwise, be sure its value is normalized.
    else
    {
        RtNormal3*  writeOrientation = const_cast< RtNormal3* >( Orientation );
        for( int i=0; i < numPts; ++i ) writeOrientation[i].Normalize();
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
                                        Orientation,
                                        Length,
                                        Vn
                                        );

    return bxdf;
}
// Releases the RixBxdf object created by BeginScatter().
void BxdfFactory::EndScatter( RixBxdf* ) {}


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

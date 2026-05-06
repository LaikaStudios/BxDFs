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
 *      2022 Prater "scatter" response.
 *      Simulates the presence of a scattering medium in a boundary
 *      layer of the surface, like dust or a random jumble of fibers.
 *      Originally developed in the 1990's, but has undergone continuous
 *      revision since then. Unpublished.
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

    // Default parameter values (from the .args file).
    const float      def_Gain = 1.0f;
    const RtColorRGB def_Color = RtColorRGB( 0.5f );
    const float      def_Direction = 0.5f;
    const float      def_Dispersion = 0.5f;

    // pTable indices.
    enum pTableEntries
    {
        in_Gain = 0,
        in_Color,
        in_Direction,
        in_Dispersion
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
            RixSCParamInfo() // Ends the table.
        };

        return &pTable[0];
    }

    // Sets this bxdf's characteristic Light Path Expression variables at RenderBegin.
    void Synchronize( RixContext&, RixSCSyncMsg, const RixParameterList* );

    // GetInstanceHints() provides information to the integrator
    // about this bxdf's presence/opacity/interior computations.
    int GetInstanceHints( void* data ) const { return k_TriviallyOpaque; }

    // No useful comments about this in RixBxdf.h, and all examples just return 1.
    float GetIndexOfRefraction( void* data ) const { return 1.0f; }

    // Unused here.
    int  Init( RixContext&, const RtUString ) { return 0; }
    void Finalize( RixContext& ctx ) {}
    void CreateInstanceData( RixContext&, const RtUString, const RixParameterList*, InstanceData* ) {}
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
    const int         numPts;
    const float*      Gain;
    const RtColorRGB* Color;
    const float*      Direction;
    const float*      Dispersion;
    const RtNormal3*  Nn;
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
        const float*      _Dispersion,
        const RtNormal3*  _Nn,
        const RtNormal3*  _Ng,
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
        Direction( _Direction ),
        Dispersion( _Dispersion ),
        Nn( _Nn ),
        Ng( _Ng ),
        Vn( _Vn )
    {
        // Intersect the response(s) wanted by the integrator (lobesWanted)
        // with the response(s) this bxdf produces. The result defines the set
        // of responses we need to compute in the BxdfClosure::*Sample() methods.
        bxdfLobes &= sg_PraterScatter_LT; // Additional response sg_*_LT values are | together.
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
    const int  numPts = sCtx->numPts;

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
    float*  writeDirection = const_cast< float* >( Direction );
    float*  writeDispersion = const_cast< float* >( Dispersion );
    for( int i=0; i < numPts; ++i )
    {
        writeDirection[i] = 0.99f * Direction[i]; // -1 < g < +1
        writeDispersion[i] = 1.0f - Dispersion[i]; // 0 ≤ f ≤ 1
    }

    // Get some shading context variable pointers.
    const RtNormal3*  Nn;
    const RtNormal3*  Ng;
    const RtVector3*  Vn;
    sCtx->GetBuiltinVar( RixShadingContext::k_Nn,  &Nn ); // wn
    sCtx->GetBuiltinVar( RixShadingContext::k_Ngn, &Ng ); // wg
    sCtx->GetBuiltinVar( RixShadingContext::k_Vn,  &Vn ); // wo

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
                                            Dispersion,
                                            Nn,
                                            Ng,
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

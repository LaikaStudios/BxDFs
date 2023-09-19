/*
 *  Copyright 2022-2023 LAIKA. Authored by Mitch Prater.
 *
 *  Licensed under the Apache License Version 2.0 http://apache.org/licenses/LICENSE-2.0,
 *  or the MIT license http://opensource.org/licenses/MIT, at your option.
 *
 *  This program may not be copied, modified, or distributed except according to those terms.
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>

#include "RixRNG.h"
#include "RixBxdf.h"
#include "RixBxdfLobe.h"
#include "RixShading.h"
#include "RixShadingUtils.h"
#include "RixPredefinedStrings.hpp"
#include "bxdf/PxrSurfaceOpacity.h"

// Provides a non-const pointer to a const data set so it can be modified in place.
#define NON_CONST_PTR(type,ptr) const_cast< type >( static_cast< const type >( ptr ))


//===================================================================
//  2022 Prater "fuzz" response function.
//  Simulates the presence of a surface boundary layer consisting
//  of fibers oriented perpendicularly to the surface: e.g. velvet.
//  Originally developed in the 1990's, but has undergone continuous
//  revision since then. Unpublished.
//===================================================================

/*
 *  Notation (vectors originate at the surface):
 *      wi - incident direction: the "light" direction. a.k.a. "incoming".
 *      wo - observer direction: the "view" direction. a.k.a. "outgoing".
 *      wg - the geometric (modeled surface) normal.
 *      ws - surface shading normal. possibly "bumped" relative to wg.
 *      wn - response computation normal. generally equivalent to ws.
 *      Cs - the substance's characteristic (un-lit) coloration (spectrum).
 *      Cr - the response color & magnitude (spectrum*intensity).
 *      W  - the response weight.
 *      fPdf - forward pdf: probability of light moving from wi toward wo.
 *      rPdf - reverse pdf: probability of light moving from wo toward wi.
 */

PRMAN_INLINE
float PraterFuzzResponse
(
    const float  g, // Direction: -1 < g < +1
    const float  f, // 1-Dispersion: 0 ≤ f ≤ 1
    const float  cos_theta, // Θ = wi ∠ wo
    const float  cos_wnwi,
    const float  cos_wnwo
)
{
    // Forward/backward scattering weight.
    // t = ray/unit-circle intersection distance.
    // ray direction d = ( cos_theta, sin_theta )
    // ray origin p = ( g, 0 ); 0 < g ⇒ forward scattering.
    const float  dot_pp = g*g;
    const float  dot_pd = cos_theta*g;
    const float  t = std::sqrt( dot_pd*dot_pd - dot_pp + 1.0f ) - dot_pd;

    // Asperity.
    const float  asperity_wi = 1.0f - std::abs( cos_wnwi );
    const float  asperity_wo = 1.0f - std::abs( cos_wnwo );

    const float  asperity_broad = (asperity_wi + asperity_wo) * 0.5f;
    const float  asperity_narrow = asperity_wi * asperity_wo;
    const float  asperity = f*f*(asperity_narrow - asperity_broad) + asperity_broad;

    // Response function.
    const float  e = 0.5f + 3.5f*f; // lerp( 0.5, 4.0, f );
    const float  r = std::pow( t * asperity, e );

    return r;
}


// Uniformly and regularly spaced spherical sample directions
// used for (brute force) numerical response integration, which
// is then used to normalize the response energy.
#include "bxdf/sampleDir.h"


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
        // Stores the bitfield value that will be returned by
        // GetInstanceHints(), which informs the integrator what
        // types of presence/opacity/interior computations this
        // bxdf plugin performs.
        int  poiHints;

        // Orientation param connection info.
        RixSCConnectionInfo  cinfoOrientation;

        // ½ degree increment wn ∠ wi response integral values in the
        // [0,π/2] range for the current g and f.
        float  g;
        float  f;
        std::array< float, 181 >  responseIntegral;

        // Compute the response integration values at every ½ degree
        // of wn ∠ wi in [0,π/2] based on the current g and f.
        void SetNormalization( const float _g, const float _f )
        {
            g = _g;
            f = _f;

            const RtVector3  wn = RtVector3( 0.0, 0.0, 1.0 );

            for( int halfDeg = 0; halfDeg < responseIntegral.size(); halfDeg++ ) // [0,π/2]
            {
                const float  cos_wnwi = std::cos( halfDeg * (F_PI / 360.0f) );
                const RtVector3  wi = RtVector3( 0.0, cos_wnwi, std::sqrt( 1.0 - cos_wnwi*cos_wnwi ));

                double  sum = 0.0;
                const int  sampleNum = sizeof( sampleDir ) / sizeof(float) / 3; // Number of sample vectors.
                for( int s = 0; s < sampleNum; s++ )
                {
                    RtVector3  wo = RtVector3( sampleDir[s][0], sampleDir[s][1], sampleDir[s][2] );
                    const float  cos_wnwo = wo.Dot( wn );
                    const float  cos_theta = wo.Dot( wi );

                    sum += PraterFuzzResponse( g, f, cos_theta, cos_wnwi, cos_wnwo );
                }
                sum *= F_FOURPI / sampleNum; // Apply sample solid angle scale outside the loop.

                responseIntegral[ halfDeg ] = sum;
            }
        }

        // Interpolate the response integral data for the given
        // cos_wnwi and return the response normalization value.
        // Note: response is symmetric about cos_wnwi = 0.
        float GetNormalization( const float cos_wnwi )
        {
            const float halfDeg = std::acos( std::abs( cos_wnwi )) * (360.0f / F_PI);
            const int   i0 = std::floor( halfDeg );
            const int   i1 = std::ceil( halfDeg );
            const float t = halfDeg - i0;
            return  1.0f / RixMix( responseIntegral[i0], responseIntegral[i1], t );
        }

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
    const float      def_Direction = -0.5f;
    const float      def_Dispersion = 0.5f;
    const float      def_Presence = 1.0f;
    const int        def_Socket = 0;

    // pTable indices.
    enum pTableEnries
    {
        in_Gain = 0,
        in_Color,
        in_Orientation,
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
            RixSCParamInfo( RtUString( "Orientation" ), k_RixSCNormal, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Direction" ), k_RixSCFloat, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Dispersion" ), k_RixSCFloat, k_RixSCScatterInput ),
            RixSCParamInfo( RtUString( "Presence" ), k_RixSCFloat, k_RixSCPresenceInput ),
            RixSCParamInfo( RtUString( "Socket" ), k_RixSCInteger ),
            RixSCParamInfo() // Ends the table.
        };

        return &pTable[0];
    }

    // Get any general-purpose RixInterface handles we need.
    int Init( RixContext&, const RtUString ) { return 0; }

    // Since Init() didn't allocate any memory, Finalize() is a no-op.
    void Finalize( RixContext& ctx ) {}

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

        // Convert Direction and Dispersion parameter values to g and f
        // and use those to compute the response normalization data.
        float  Direction( def_Direction );
        float  Dispersion( def_Dispersion );
        pList->EvalParam( in_Direction, -1, &Direction );
        pList->EvalParam( in_Dispersion, -1, &Dispersion );

        const float  g = 0.99f * Direction; // -1 < g < +1
        const float  f = 1.0f - Dispersion; // 0 ≤ f ≤ 1
        pluginData->SetNormalization( g, f );

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

        // Get the Orientation parameter's connection info.
        pList->GetParamInfo( in_Orientation, &type, &(pluginData->cinfoOrientation) );

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
static RixBXLobeSampled  sg_PraterFuzz_LS;
// This data type can contain any number of responses. Used to create/define sets of responses.
static RixBXLobeTraits   sg_PraterFuzz_LT;

// Synchronize() sets the Lobe variables.
// This can't be done statically since RixBXLookupLobeByName() requires the
// Light Path Expression system which is not available until k_RixSCRenderBegin.
void BxdfFactory::Synchronize( RixContext& ctx, RixSCSyncMsg syncMsg, const RixParameterList* pList )
{
    if( syncMsg != k_RixSCRenderBegin ) return;

    // Query the Light Path Expression (LPE) entry for the "PraterFuzz"
    // response (specified in the rendermn.ini file) and set its traits.
    // Each response will require its own static global variables and
    // RixBXLookupLobeByName() and RixBXLobeTraits() calls to set them.
    // rendermn.ini entry: /prman/lpe/specular6  PraterFuzz
    // LPE: color lpe:CS6.*[<L.>O]
    sg_PraterFuzz_LS = RixBXLookupLobeByName( ctx,
                        false, // not discrete ⇒ samples over a solid angle.
                        true,  // specular ⇒ not diffuse.
                        true,  // reflected scattering ⇒ not transmitted.
                        false, // not a user response ⇒ standard (specular or diffuse).
                        0,     // response "id" number. not used by prman or LPE.
                        "PraterFuzz" // The response's Name (used in rendermn.ini).
                        );

    sg_PraterFuzz_LT = RixBXLobeTraits( sg_PraterFuzz_LS );
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

    // Shading instance data needed to compute this bxdf.
    BxdfFactory::pluginInstanceData*  pluginData; // InstanceData->data.

    // Shading context data needed to compute this bxdf.
    // These will consist of numPts values: one for each shaded point.
    // Note: a shading context is also known as a rendering "grid".
    const int         numPts;
    const float*      Gain;
    const RtColorRGB* Color;
    const RtNormal3*  Orientation;
    const RtNormal3*  Ng;
    const RtVector3*  Vn;

    void PraterFuzzPdf
    (
        const float  cos_theta, // Θ = wi ∠ wo
        const float  cos_wnwi,
        const float  cos_wnwo,
        // Results:
        float&  fPdf,
        float&  rPdf,
        float&  W
    )
    {
        W = PraterFuzzResponse( pluginData->g, pluginData->f, cos_theta, cos_wnwi, cos_wnwo )
          * pluginData->GetNormalization( cos_wnwi );
        
        // Volume scattering, so we don't scale by projected solid angle.
        // fPdf = W * cos_wgwi;
        // rPdf = W * cos_wgwo;
        rPdf = fPdf = W;
    }

    void PraterFuzz
    (
        const float  cos_theta,
        const float  cos_wnwi,
        const float  cos_wnwo,
        const RtColorRGB  Cs,
        // Results:
        float&  fPdf,
        float&  rPdf,
        RtColorRGB& Cr
    )
    {
        float  W;
        PraterFuzzPdf( cos_theta, cos_wnwi, cos_wnwo, fPdf, rPdf, W );
        Cr = Cs * W;
    }

    bool Generate
    (
        const RtNormal3  wg,
        const RtVector3  wn,
        const RtVector3  wo,
        const RtFloat2   xi,
        const RtColorRGB Cs,
        // Results:
        RtVector3&  wi,
        float&  fPdf,
        float&  rPdf,
        RtColorRGB& Cr
    )
    {
        // Test the observer visibility.
        if( wg.Dot(wo) < 0.00001f ) return false;

        // Generate an incident sample direction (wi).
        RtVector3  wt, wb;
        wg.CreateOrthonormalBasis( wt, wb );

        // Toroidal distribution produced by squaring xi.y.
        float  dummy;
        const RtFloat2  newXi = RtFloat2( xi[0], xi[1]*xi[1] );
        RixUniformDirectionalDistribution( newXi, wg, wt, wb, wi, dummy );

        // Test the incident sample visibility: reject those below the horizon.
        // No rejection testing necessary, since only reflection hemisphere
        // samples are generated.
        // if( wg.Dot(wi) < 0.00001f ) return false;

        const float  cos_theta = wi.Dot(wo);
        const float  cos_wnwi = wn.Dot(wi);
        const float  cos_wnwo = wn.Dot(wo);

        PraterFuzz( cos_theta, cos_wnwi, cos_wnwo, Cs, fPdf, rPdf, Cr );
        return true;
    }

    bool Evaluate
    (
        const RtNormal3  wg,
        const RtVector3  wn,
        const RtVector3  wo,
        const RtVector3  wi,
        const RtColorRGB Cs,
        // Results:
        float&  fPdf,
        float&  rPdf,
        RtColorRGB& Cr
    )
    {
        // Test the observer and incident visibility.
        if( wg.Dot(wo) < 0.00001f ) return false;
        if( wg.Dot(wi) < 0.00001f ) return false;

        const float  cos_theta = wi.Dot(wo);
        const float  cos_wnwi = wn.Dot(wi);
        const float  cos_wnwo = wn.Dot(wo);

        PraterFuzz( cos_theta, cos_wnwi, cos_wnwo, Cs, fPdf, rPdf, Cr );
        return true;
    }

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
        BxdfFactory::pluginInstanceData* _pluginData, // InstanceData->data.
        const int         _numPts,
        const float*      _Gain,
        const RtColorRGB* _Color,
        const RtNormal3*  _Orientation,
        const RtNormal3*  _Ng,
        const RtVector3*  _Vn
    ):
        // Initializes the protected RixBxdf class members
        // 'shadingCtx' (shading context) to sCtx, and 'bxdfFactory' to bFac.
        RixBxdf( sCtx, bFac ),

        // Initialize this BxdfClosure's bxdfLobes member variable
        // to the responses the integrator wants from it.
        bxdfLobes( lobesWanted ),

        // Initialize this BxdfClosure's other member variables.
        pluginData( _pluginData ),
        numPts( _numPts ),
        Gain( _Gain ),
        Color( _Color ),
        Orientation( _Orientation ),
        Ng( _Ng ),
        Vn( _Vn )
    {
        // Intersect the response(s) wanted by the integrator (lobesWanted)
        // with the response(s) this bxdf produces. The result defines the set
        // of responses we need to compute in the BxdfClosure::*Sample() methods.
        bxdfLobes &= sg_PraterFuzz_LT; // Additional response sg_*_LT values are | together.
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
    // #include "bxdf/PraterFuzzSampling.inl"

    /*
    ================================================================
    GenerateSample() provides the integrator with a shading context
    set of samples generated from this bxdf's response lobe(s).
    ================================================================
    */
    void GenerateSample
    (
        RixBXTransportTrait     transportTrait, // Direct, indirect, or both bit field.
        const RixBXLobeTraits*  lobesWanted, // by the integrator, per shading context point.
        RixRNG*                 rixRng, // handle to the random number generator.
        // Generated results:
        RixBXLobeSampled*       lobeGenerated, // which lobe type was generated.
        RtVector3*              wi, // incoming sample direction per shading context point.
        RixBXLobeWeights&       lobeWeights, // response weights.
        float*                  fPdf, // forward pdf.
        float*                  rPdf, // reverse pdf.
        RtColorRGB*             compTrans // compositing transparency (PRMan "Oi")?
    )
    {
        RtFloat2*  xi = static_cast< RtFloat2* >( RixAlloca( numPts*sizeof(RtFloat2) ));
        rixRng->DrawSamples2D( xi );

        const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

        RtColorRGB*  PraterFuzzWeight = NULL;

        for( int i=0; i < numPts; i++ )
        {
            lobeGenerated[i].SetValid( false );

            const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

            const bool  doPraterFuzz = ( lobesToConsider & sg_PraterFuzz_LT ).HasAny();

            if( doPraterFuzz )
            {
                const RtColorRGB Cs = Color[i]*Gain[i];
                const RtNormal3  wn = Orientation[i];
                const RtNormal3  wg = Ng[i];
                const RtVector3  wo = Vn[i];

                if( !PraterFuzzWeight ) PraterFuzzWeight = lobeWeights.AddActiveLobe( sg_PraterFuzz_LS );

                if( Generate( wg, wn, wo, xi[i], Cs, wi[i], fPdf[i], rPdf[i], PraterFuzzWeight[i] ))
                {
                    lobeGenerated[i] = sg_PraterFuzz_LS;
                }
            }
        }
    }

    /*
    ================================================================
    EvaluateSample() provides the integrator with a shading context
    set of samples evaluated using this bxdf's response lobe(s).
    ================================================================
    */
    void EvaluateSample
    (
        RixBXTransportTrait     transportTrait, // Direct, indirect, or both bit field.
        const RixBXLobeTraits*  lobesWanted, // by the integrator, per shading context point.
        RixRNG*                 rixRng, // handle to the random number generator.
        RixBXLobeTraits*        lobesEvaluated, // Returned value.
        const RtVector3*        wi, // incoming sample direction per shading context point.
        // Evaluated results:
        RixBXLobeWeights&       lobeWeights, // sample weight.
        float*                  fPdf, // forward pdf.
        float*                  rPdf  // reverse pdf.
    )
    {
        const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

        RtColorRGB*  PraterFuzzWeight = NULL;

        for( int i=0; i < numPts; i++ )
        {
            lobesEvaluated[i].SetNone();

            const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

            const bool  doPraterFuzz = ( lobesToConsider & sg_PraterFuzz_LT ).HasAny();

            if( doPraterFuzz )
            {
                const RtColorRGB Cs = Color[i]*Gain[i];
                const RtNormal3  wn = Orientation[i];
                const RtNormal3  wg = Ng[i];
                const RtVector3  wo = Vn[i];

                if( !PraterFuzzWeight ) PraterFuzzWeight = lobeWeights.AddActiveLobe( sg_PraterFuzz_LS );

                if( Evaluate( wg, wn, wo, wi[i], Cs, fPdf[i], rPdf[i], PraterFuzzWeight[i] ))
                {
                    lobesEvaluated[i] |= sg_PraterFuzz_LT;
                }
            }
        }
    }

    /*
    =============================================================
    Like EvaluateSample(), but does multiple evaluations of this
    bxdf's response lobe(s) at a single shading context point.
    =============================================================
    */
    void EvaluateSamplesAtIndex
    (
        RixBXTransportTrait     transportTrait, // Direct, indirect, or both bit field.
        const RixBXLobeTraits&  lobesWanted, // by the integrator, at the scIndex point.
        RixRNG*                 rixRng, // handle to the random number generator.
        int                     scIndex, // shading context point to evaluate.
        int                     nSamps, // number of wi samples to evaluate.
        RixBXLobeTraits*        lobesEvaluated, // Returned value.
        const RtVector3*        wi, // nSamps incoming sample directions.
        // Evaluated results:
        RixBXLobeWeights&       lobeWeights, // sample weight.
        float*                  fPdf, // forward pdf.
        float*                  rPdf  // reverse pdf.
    )
    {
        const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();
        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted;

        const bool  doPraterFuzz = ( lobesToConsider & sg_PraterFuzz_LT ).HasAny();

        RtColorRGB*  PraterFuzzWeight = NULL;
        if( doPraterFuzz ) PraterFuzzWeight = lobeWeights.AddActiveLobe( sg_PraterFuzz_LS );

        for( int i=0; i < nSamps; i++ )
        {
            lobesEvaluated[i].SetNone();

            if( doPraterFuzz )
            {
                const RtColorRGB Cs = Color[scIndex]*Gain[scIndex];
                const RtNormal3  wn = Orientation[scIndex];
                const RtNormal3  wg = Ng[scIndex];
                const RtVector3  wo = Vn[scIndex];

                if( Evaluate( wg, wn, wo, wi[i], Cs, fPdf[i], rPdf[i], PraterFuzzWeight[i] ))
                {
                    lobesEvaluated[i] |= sg_PraterFuzz_LT;
                }
            }
        }
    }

    /*
    ======================================================================
    EmitLocal() produces a shading context set of this bxdf's baked
    (pre-integrated) illumination responses and/or light emission results.
    ======================================================================
    */
    bool EmitLocal( RtColorRGB* ) { return false; } // None in this case.
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
    auto  pluginData = static_cast< BxdfFactory::pluginInstanceData* >( data );

    // Get some shading context variable pointers.
    const RtNormal3*  Nn;
    const RtNormal3*  Ng;
    const RtVector3*  Vn;
    sCtx->GetBuiltinVar( RixShadingContext::k_Nn,  &Nn ); // wn
    sCtx->GetBuiltinVar( RixShadingContext::k_Ngn, &Ng ); // wg
    sCtx->GetBuiltinVar( RixShadingContext::k_Vn,  &Vn ); // wo

    // Evaluate the Socket parameter to trigger connected node evaluation.
    const int*  Socket;
    sCtx->EvalParam( in_Socket, -1, &Socket );

    // Evaluate (the potentially varying) parameters.
    const float*      Gain;
    const RtColorRGB* Color;
    const RtNormal3*  Orientation;

    sCtx->EvalParam( in_Gain, -1, &Gain, &def_Gain, true );
    sCtx->EvalParam( in_Color, -1, &Color, &def_Color, true );
    sCtx->EvalParam( in_Orientation, -1, &Orientation, &def_Orientation, true );

    // Point Orientation to Nn if it hasn't been set to anything.
    if( k_RixSCDefaultValue == pluginData->cinfoOrientation )
    {
        Orientation = Nn;
    }
    // Otherwise, be sure its value is normalized.
    else
    {
        RtNormal3*  writeOrientation = NON_CONST_PTR( RtNormal3*, Orientation );
        for( int i=0; i < numPts; ++i ) writeOrientation[i].Normalize();
    }

    // Create a shading context memory pool.
    RixShadingContext::Allocator  pool( sCtx );

    // Allocate memory for this shader's RixBxdf object.
    void*  mem = pool.AllocForBxdf< BxdfClosure >(1);

    // Create an instance of this shader's bxdf closure
    // and pass any necessary data to it.
    BxdfClosure*  bxdf = new (mem) BxdfClosure( this, sCtx, lobesWanted,
                                            pluginData,
                                            numPts,
                                            Gain,
                                            Color,
                                            Orientation,
                                            Ng,
                                            Vn
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

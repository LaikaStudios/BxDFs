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

#include "prmanapi.h"
#include "RiTypesHelper.h"

#include "RixBxdf.h"
#include "RixBxdfLobe.h"
#include "RixShading.h"
#include "RixShadingUtils.h"

#include <cmath>

/*
 *  Defines the response sample methods of the RixBxdf object:
 *      GenerateSample()
 *      EvaluateSample()
 *      EvaluateSamplesAtIndex()
 *      EmitLocal() - optional
 *
 *  These methods implement the response(s) produced by this plugin.
 * 
 *  Notation (vectors originate at the surface):
 *      wi - incident direction: the "light" direction. a.k.a. "incoming".
 *      wo - observer direction: the "view" direction. a.k.a. "outgoing".
 *      wg - the geometric (modeled + displaced surface) normal.
 *      wn - surface shading normal. possibly "bumped" relative to wg.
 *      Cs - the substance's characteristic (un-lit) coloration (albedo).
 *      Cr - the response color & magnitude (albedo*intensity).
 *      w  - the response weight (intensity).
 *      fPdf - forward pdf: probability of light moving from wi toward wo.
 *      rPdf - reverse pdf: probability of light moving from wo toward wi.
 *      Unless otherwise noted:
 *      θ - spherical coordinate polar angle from +z (the surface normal).
 *      φ - spherical coordinate azimuth angle around z.
 */

inline
float KPAsperityResponse
(
    const float          lm, // ∆/λ : fiber length / scattering mean free path.
    const RtVector3      wn, // n
    const RtVector3      wo, // u
    const RtVector3      wi, // v
    const float    cos_wnwo, // u·n
    const float    cos_wnwi, // v·n
    const float    cos_wowi  // p(-u·v)
)
{
    const RtVector3  u_add_v = wo + wi; // u+v
    const float  u_add_v_dot_n = u_add_v.Dot(wn); // (u+v)·n

    const float  dot_uv = cos_wowi; // p(-u·v)
    const float  dot_un = cos_wnwo; // u·n
    const float  dot_vn = cos_wnwi; // v·n

    const float  r = dot_uv*( dot_vn / u_add_v_dot_n )*( 1.0f - std::exp( -lm * u_add_v_dot_n/( dot_un*dot_vn )));

    return r;
}

inline
void KPAsperity
(
    const RtColorRGB     Cs,
    const float          lm, // ∆/λ : fiber length / scattering mean free path.
    const RtVector3      wn, // n
    const RtVector3      wo, // u
    const RtVector3      wi, // v
    const float    cos_wnwo, // u·n
    const float    cos_wnwi, // v·n
    const float    cos_wowi, // p(-u·v)

    // Results:
    float&  fPdf,
    float&  rPdf,
    RtColorRGB&  Cr
)
{
    // Response.
    float  w = KPAsperityResponse( lm, wn, wo, wi, cos_wnwo, cos_wnwi, cos_wowi );
    Cr = Cs * w;

    // Pdf. This must match the sampling distribution, not the response function.
    // This is because Generate() and Evaluate() must report the same pdf for a
    // given sample direction to ensure their results are consistent with each other.
    // But regardless of what samples are used for integration, the response itself
    // is always evaluated using the bxdf's normalized response function.
    //
    // The sample generation produces a distribution whose pdf is 1/(4π√cosθ).
    fPdf = F_INVFOURPI / std::sqrt( std::max( 0.000001f, cos_wnwi ));
    rPdf = F_INVFOURPI / std::sqrt( std::max( 0.000001f, cos_wnwo ));
}


//===================================================================
// Generate() and Evaluate() connect the response functions above to
// the API methods below. They also perform any necessary visibility
// or other sanity checks and can bypass a sample if necessary.
//===================================================================
inline
bool Generate
(
    const RtColorRGB  Cs,
    const float       lm, // ∆/λ : fiber length / scattering mean free path.
    const RtNormal3   wn, // n
    const RtVector3   wo, // u
    const RtFloat2    xi,

    // Results:
    RtVector3&  wi,
    float&  fPdf,
    float&  rPdf,
    RtColorRGB&  Cr
)
{
    const float  cos_wnwo = wn.Dot(wo); // u·n

    // Generate an incident sample direction (wi).
    // Stub-in isotropic uniform samples for now.
    RtVector3  wt, wb;
    wn.CreateOrthonormalBasis( wt, wb );

    // Squaring xi[1] gives more importance to samples near the
    // horizon and produces a distribution with a pdf of 1/(4π√cosθ).
    const RtFloat2  newXi = RtFloat2( xi[0], xi[1]*xi[1] );
    float  cos_wnwi; // v·n
    RixUniformDirectionalDistribution( newXi, wn, wt, wb, wi, cos_wnwi );

    const float  cos_wowi = std::max( 0.0f, -(wo.Dot(wi)) ); // p(-u·v)

    KPAsperity( Cs, lm, wn, wo, wi, cos_wnwo, cos_wnwi, cos_wowi, fPdf, rPdf, Cr );
    return true;
}

inline
bool Evaluate
(
    const RtColorRGB  Cs,
    const float       lm, // ∆/λ : fiber length / scattering mean free path.
    const RtNormal3   wn, // n
    const RtVector3   wo, // u
    const RtVector3   wi, // v

    // Results:
    float&  fPdf,
    float&  rPdf,
    RtColorRGB&  Cr
)
{
    const float  cos_wnwo = wn.Dot(wo); // u·n
    const float  cos_wnwi = wn.Dot(wi); // v·n
    const float  cos_wowi = std::max( 0.0f, -(wo.Dot(wi)) ); // p(-u·v)

    KPAsperity( Cs, lm, wn, wo, wi, cos_wnwo, cos_wnwi, cos_wowi, fPdf, rPdf, Cr );
    return true;
}


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

    RtColorRGB*  KPAsperityWeight = NULL;

    for( int i=0; i < numPts; ++i )
    {
        lobeGenerated[i].SetValid( false );

        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        const bool  doKPAsperity = ( lobesToConsider & sg_KPAsperity_LT ).HasAny();

        if( doKPAsperity )
        {
            // Create the sample weight array for this response, if needed.
            if( !KPAsperityWeight ) KPAsperityWeight = lobeWeights.AddActiveLobe( sg_KPAsperity_LS );

            // Initialize any data needed to compute the response.
            const RtColorRGB Cs = Color[i]*Gain[i];
            const float      lm = Length[i];
            const RtNormal3  wn = Orientation[i];
            const RtVector3  wo = Vn[i];

            if( Generate( Cs, lm, wn, wo, xi[i], wi[i], fPdf[i], rPdf[i], KPAsperityWeight[i] ))
            {
                lobeGenerated[i] = sg_KPAsperity_LS;
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

    RtColorRGB*  KPAsperityWeight = NULL;

    for( int i=0; i < numPts; ++i )
    {
        lobesEvaluated[i].SetNone();

        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        const bool  doKPAsperity = ( lobesToConsider & sg_KPAsperity_LT ).HasAny();
        if( doKPAsperity )
        {
            // Initialize any data needed to compute the response.
            const RtColorRGB Cs = Color[i]*Gain[i];
            const float      lm = Length[i];
            const RtNormal3  wn = Orientation[i];
            const RtVector3  wo = Vn[i];

            if( !KPAsperityWeight ) KPAsperityWeight = lobeWeights.AddActiveLobe( sg_KPAsperity_LS );

            if( Evaluate( Cs, lm, wn, wo, wi[i], fPdf[i], rPdf[i], KPAsperityWeight[i] ))
            {
                lobesEvaluated[i] = sg_KPAsperity_LT;
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

    RtColorRGB*  KPAsperityWeight = NULL;

    for( int i=0; i < nSamps; ++i )
    {
        lobesEvaluated[i].SetNone();


        const bool  doKPAsperity = ( lobesToConsider & sg_KPAsperity_LT ).HasAny();
        if( doKPAsperity )
        {
            if( !KPAsperityWeight ) KPAsperityWeight = lobeWeights.AddActiveLobe( sg_KPAsperity_LS );

            const RtColorRGB Cs = Color[scIndex]*Gain[scIndex];
            const float      lm = Length[scIndex];
            const RtNormal3  wn = Orientation[scIndex];
            const RtVector3  wo = Vn[scIndex];

            if( Evaluate( Cs, lm, wn, wo, wi[i], fPdf[i], rPdf[i], KPAsperityWeight[i] ))
            {
                lobesEvaluated[i] = sg_KPAsperity_LT;
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

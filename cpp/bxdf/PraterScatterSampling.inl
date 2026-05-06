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
float PraterScatterResponse
(
    const float  g, // Direction: -1 < g < +1
    const float  f, // 1-Dispersion: 0 ≤ f ≤ 1
    const float  cos_theta // Phase angle Θ = wi ∠ wo
)
{
    // Forward/backward scattering weight.
    // t = ray/unit-circle intersection distance.
    // ray direction d = ( cos_theta, sin_theta )
    // ray origin p = ( g, 0 ); 0 < g ⇒ forward scattering.
    const float  dot_pp = g*g;
    const float  dot_pd = cos_theta*g;
    const float  t = std::sqrt( dot_pd*dot_pd - dot_pp + 1.0f ) - dot_pd;

    // Response function.
    const float  e = 0.5 + 3.5*f; // RixMix( 0.5, 4.0, f );
    const float  r = std::pow( t, e );

    return r;
}

// Compute the response normalization for the entire (spherical domain)
// response lobe for the given parameters g and f.
inline
float PraterScatterNormalize
(
    const float  g,
    const float  f
)
{
    const float  p00 =  0.1255e+2;
    const float  p10 = -0.2098e+0;
    const float  p01 =  0.6372e+0;
    const float  p20 =  0.1247e+1;
    const float  p11 = -0.7389e+0;
    const float  p02 = -0.2639e+1;
    const float  p30 = -0.3036e+1;
    const float  p21 = -0.2737e+2;
    const float  p12 =  0.1714e+2;
    const float  p03 =  0.1369e+1;
    const float  p40 =  0.6518e+1;
    const float  p31 =  0.3178e+1;
    const float  p22 =  0.4763e+2;
    const float  p13 = -0.3391e+2;
    const float  p04 =  0.3654e+1;
    const float  p50 = -0.7246e+1;
    const float  p41 =  0.1571e+2;
    const float  p32 = -0.2700e+2;
    const float  p23 =  0.2575e+1;
    const float  p14 =  0.1531e+2;
    const float  p05 = -0.2945e+1;

    // Polynomial represents ½ the integral surface, which is symmetric
    // about g=0: ±g produces the same spherical domain integral value.
    const float  g1 = std::abs( g );
    const float  g2 = g1*g1;
    const float  g3 = g1*g2;
    const float  g4 = g1*g3;
    const float  g5 = g1*g4;
    const float  f2 = f*f;
    const float  f3 = f*f2;
    const float  f4 = f*f3;
    const float  f5 = f*f4;

    const float  integral = p00 + p10*g1 + p01*f + p20*g2 + p11*g1*f + p02*f2 + p30*g3 + p21*g2*f 
                            + p12*g1*f2 + p03*f3 + p40*g4 + p31*g3*f + p22*g2*f2 
                            + p13*g1*f3 + p04*f4 + p50*g5 + p41*g4*f + p32*g3*f2 
                            + p23*g2*f3 + p14*g1*f4 + p05*f5;

    return 1.0 / integral;
}

inline
void PraterScatter
(
    const RtColorRGB  Cs,
    const float  g,
    const float  f,
    const float  cos_theta, // Phase angle Θ = wi ∠ wo

    // Results:
    float&  fPdf,
    float&  rPdf,
    RtColorRGB&  Cr
)
{
    // Response.
    float  w = PraterScatterResponse( g, f, cos_theta ) * PraterScatterNormalize( g, f );
    Cr = Cs * w;

    // Pdf. This must match the sampling distribution, not the response function.
    // This is because Generate() and Evaluate() must report the same pdf for a
    // given sample direction to ensure their results are consistent with each other.
    // But regardless of what samples are used for integration, the response itself
    // is always evaluated using the bxdf's normalized response function.
    //
    // Since the generated samples are uniformly distributed over the hemisphere,
    // each sample has the same probability, which is equal to the inverse of the
    // hemisphere's area: 1/(2π).
    rPdf = fPdf = F_INVTWOPI;
}


//===================================================================
// Generate() and Evaluate() connect the response functions above to
// the API methods below. They also perform any necessary visibility
// or other sanity checks and can bypass a sample if necessary.
//===================================================================
inline
bool Generate
(
    const RtColorRGB Cs,
    const float  g,
    const float  f,
    const RtNormal3  wn,
    const RtNormal3  wg,
    const RtVector3  wo,
    const RtFloat2   xi,

    // Results:
    RtVector3&  wi,
    float&  fPdf,
    float&  rPdf,
    RtColorRGB&  Cr
)
{
    // Test the observer visibility.
    //if( wg.Dot(wo) < 0.00001f ) return false;

    // Generate an incident sample direction (wi).
    // Stub-in isotropic uniform samples for now.
    RtVector3  wt, wb;
    wn.CreateOrthonormalBasis( wt, wb );

    float  dummy;
    RixUniformDirectionalDistribution( xi, wn, wt, wb, wi, dummy );

    // Test the incident sample visibility: reject those below the horizon.
    //if( wg.Dot(wi) < 0.00001f ) return false;

    const float  cos_theta = wi.Dot(wo);

    PraterScatter( Cs, g, f, cos_theta, fPdf, rPdf, Cr );
    return true;
}

inline
bool Evaluate
(
    const RtColorRGB Cs,
    const float  g,
    const float  f,
    const RtNormal3  wn,
    const RtNormal3  wg,
    const RtVector3  wo,
    const RtVector3  wi,

    // Results:
    float&  fPdf,
    float&  rPdf,
    RtColorRGB&  Cr
)
{
    // Test the observer and incident visibility.
    //if( wg.Dot(wo) < 0.00001f ) return false;
    //if( wg.Dot(wi) < 0.00001f ) return false;

    const float  cos_theta = wi.Dot(wo);
    
    PraterScatter( Cs, g, f, cos_theta, fPdf, rPdf, Cr );
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

    RtColorRGB*  PraterScatterWeight = NULL;

    for( int i=0; i < numPts; ++i )
    {
        lobeGenerated[i].SetValid( false );

        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        const bool  doPraterScatter = ( lobesToConsider & sg_PraterScatter_LT ).HasAny();

        if( doPraterScatter )
        {
            const RtColorRGB Cs = Color[i]*Gain[i];
            const float  g = Direction[i];
            const float  f = Dispersion[i];
            const RtNormal3  wn = Nn[i];
            const RtNormal3  wg = Ng[i];
            const RtVector3  wo = Vn[i];

            if( !PraterScatterWeight ) PraterScatterWeight = lobeWeights.AddActiveLobe( sg_PraterScatter_LS );

            if( Generate( Cs, g, f, wn, wg, wo, xi[i], wi[i], fPdf[i], rPdf[i], PraterScatterWeight[i] ))
            {
                lobeGenerated[i] = sg_PraterScatter_LS;
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

    RtColorRGB*  PraterScatterWeight = NULL;

    for( int i=0; i < numPts; ++i )
    {
        lobesEvaluated[i].SetNone();

        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        const bool  doPraterScatter = ( lobesToConsider & sg_PraterScatter_LT ).HasAny();

        if( doPraterScatter )
        {
            const RtColorRGB Cs = Color[i]*Gain[i];
            const float  g = Direction[i];
            const float  f = Dispersion[i];
            const RtNormal3  wn = Nn[i];
            const RtNormal3  wg = Ng[i];
            const RtVector3  wo = Vn[i];

            if( !PraterScatterWeight ) PraterScatterWeight = lobeWeights.AddActiveLobe( sg_PraterScatter_LS );

            if( Evaluate( Cs, g, f, wn, wg, wo, wi[i], fPdf[i], rPdf[i], PraterScatterWeight[i] ))
            {
                lobesEvaluated[i] |= sg_PraterScatter_LT;
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

    const bool  doPraterScatter = ( lobesToConsider & sg_PraterScatter_LT ).HasAny();

    RtColorRGB*  PraterScatterWeight = NULL;
    if( doPraterScatter ) PraterScatterWeight = lobeWeights.AddActiveLobe( sg_PraterScatter_LS );

    for( int i=0; i < nSamps; ++i )
    {
        lobesEvaluated[i].SetNone();

        if( doPraterScatter )
        {
            const RtColorRGB Cs = Color[scIndex]*Gain[scIndex];
            const float  g = Direction[scIndex];
            const float  f = Dispersion[scIndex];
            const RtNormal3  wn = Nn[scIndex];
            const RtNormal3  wg = Ng[scIndex];
            const RtVector3  wo = Vn[scIndex];

            if( Evaluate( Cs, g, f, wn, wg, wo, wi[i], fPdf[i], rPdf[i], PraterScatterWeight[i] ))
            {
                lobesEvaluated[i] |= sg_PraterScatter_LT;
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

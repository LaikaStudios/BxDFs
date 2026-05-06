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
 *      1941 Henyey-Greenstein "scatter" response.
 *      Henyey, Louis G. and Greenstein, Jesse L.
 *      Astrophysical Journal, Vol. 93, p. 70-83 (1941)
 *      http://dx.doi.org/10.1086/144246
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
void HenyeyGreensteinPdf
(
    const float  g, // (-1,1); |g| ⇒ 1 is not very usable.
    const float  cos_wowi, // cos(θ); phase angle Θ = wo ∠ wi

    // Results:
    float&  fPdf,
    float&  rPdf,
    float&  w
)
{
    const float  one_plus_g2 = 1.0f + g*g;
    const float  one_minus_g2 = 1.0f - g*g;

    // wo (Vn) is inverted with respect to the original HG phase function definition.
    // Using -cos_wowi here so that positive g produces forward scattering as intended.
    // Also, settting γ = Ι = 1: spherical albedo = diffuse intensity = 1.
    const float  pdf = ( one_minus_g2*F_INVFOURPI ) / std::pow( one_plus_g2 + 2.0f*g*cos_wowi, 1.5f );

    rPdf = fPdf = w = pdf;
}

inline
float HenyeyGreensteinInvCdf
(
    const float  g, // g ≠ 0.
    const float  xi // [0,1)
)
{
    // Inverse CDF for HG cos(Θ); phase angle Θ = wo ∠ wi.
    const float  two_g = 2.0f*g;
    const float  g2 = g*g;

    const float  tmp = ( 1.0f - g2 ) / ( 1.0f - g + two_g*xi );
    const float  cos_theta = (( 1.0f + g2 ) - tmp*tmp ) / two_g;

    // wo (Vn) is inverted with respect to the original HG phase function definition.
    // Return -cos_theta so that positive g produces forward scattering as intended.
    return -cos_theta;
}

inline
bool HenyeyGreenstein
(
    const RtColorRGB  Cs,
    const float  g,
    const float  cos_wowi,

    // Results:
    float&  fPdf,
    float&  rPdf,
    RtColorRGB& Cr
)
{
    // Test the incident view factor (observer view factor tested prior to calling).
    //if( cos_wnwi < 0.00001f ) return false;

    float  w;
    HenyeyGreensteinPdf( g, cos_wowi, fPdf, rPdf, w );
    Cr = Cs * w;

    return true;
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
    const float     g,
    const RtNormal3  wn,
    const RtVector3  wo,
    const RtFloat2   xi,
    // Results:
    RtVector3&  wi,
    float&  fPdf,
    float&  rPdf,
    RtColorRGB& Cr
)
{
    // Test the observer view factor.
    const float  cos_wnwo = wn.Dot(wo);
    if( cos_wnwo < 0.00001f ) return false;

    // Generate an incident sample direction (wi) and its view factor.
    float  cos_theta;

    if( std::abs(g) < 0.001f ) // Isotropic.
    {
        RixUniformDirectionalDistribution( xi, wn, wi, cos_theta );
    }
    else // HG Inverse CDF.
    {
        cos_theta = HenyeyGreensteinInvCdf( g, xi.x );
        const float  sin_theta = std::sqrt( 1.0f - cos_theta*cos_theta );

        const float  phi = xi.y * F_TWOPI;
        float  sin_phi, cos_phi;
        RixSinCos( phi, &sin_phi, &cos_phi );

        RtVector3  wt, wb;
        wo.CreateOrthonormalBasis( wt, wb );

        // Samples where wn∙wi < 0 will be rejected by HenyeyGreensteinEval().
        wi = sin_theta*sin_phi*wt + sin_theta*cos_phi*wb + cos_theta*wo;
    }

    const float  cos_wowi = wo.Dot(wi);

    return HenyeyGreenstein( Cs, g, cos_wowi, fPdf, rPdf, Cr );
}

inline
bool Evaluate
(
    const RtColorRGB  Cs,
    const float  g,
    const RtNormal3  wn,
    const RtVector3  wo,
    const RtVector3  wi,
    // Results:
    float&  fPdf,
    float&  rPdf,
    RtColorRGB& Cr
)
{
    // Test the observer view factor.
    const float  cos_wnwo = wn.Dot(wo);
    if( cos_wnwo < 0.00001f ) return false;

    const float  cos_wowi = wo.Dot(wi);

    return HenyeyGreenstein( Cs, g, cos_wowi, fPdf, rPdf, Cr );
}


/*
================================================================
GenerateSample() provides the integrator with a shading context
set of samples generated from this bxdf's response lobe(s).
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
    // Get all the responses this bxdf closure instance might need to generate.
    const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

    // The response weight arrays: declare one for each response.
    // These will point to either a diffuse, specular, or user LPE classified
    // reservoir created inside the RixBXLobeWeights (lobeWeights) struct.
    // This mechanism uses the response's sg_*_LS value to determine which
    // reservoir to put the response in. See RixBxdfLobe.h
    RtColorRGB*  HGScatterWeight = NULL;

    // For each point in the shading context, evaluate its response.
    for( int i=0; i < numPts; ++i )
    {
        // Initialize this shaded point to no generated response.
        lobesEvaluated[i].SetNone();

        // Determine whether we need to evaluate a sample direction for this i by
        // intersecting the responses we're supposed to produce (bxdfLobes: which
        // is based on what the renderer wants & what this shader produces), with
        // what responses the integrator wants at this shading point i.
        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        // Do we need to compute the "HGScatter" response for this shading point?
        const bool  doHGScatter = ( lobesToConsider & sg_HGScatter_LT ).HasAny();
        if( doHGScatter )
        {
            // Initialize any data needed to compute the response.
            const RtColorRGB Cs = Color[i]*Gain[i];
            const float       g = Direction[i];
            const RtNormal3  wn = Nn[i];
            const RtVector3  wo = Vn[i];

            // Create the sample weight array for this response, if needed.
            if( !HGScatterWeight ) HGScatterWeight = lobeWeights.AddActiveLobe( sg_HGScatter_LS );

            // Evaluate the response sample.
            if( Evaluate( Cs, g, wn, wo, wi[i], fPdf[i], rPdf[i], HGScatterWeight[i] ))
            {
                // If successful, set the i sample to the response we evaluated.
                lobesEvaluated[i] = sg_HGScatter_LT;
            }
        }
    }
}


//
//  Like EvaluateSample(), but does multiple evaluations
//  of this bxdf's response at a single shading point.
//
void EvaluateSamplesAtIndex
(
    RixBXTransportTrait     transportTrait, // Direct, indirect, or both bit field.
    const RixBXLobeTraits&  lobesWanted, // by the integrator, per bxdf closure instance.
    RixRNG*                 rixRng, // handle to the random number generator.
    int                     scIndex, // shading context point to evaluate.
    int                     nSamps, // number of samples to evaluate.
    RixBXLobeTraits*        lobesEvaluated, // Returned value.
    const RtVector3*        wi, // nSamps incoming sample directions.
    // Evaluated results:
    RixBXLobeWeights&       lobeWeights, // sample weight.
    float*                  fPdf, // forward pdf.
    float*                  rPdf  // reverse pdf.
)
{
    // Get all the responses this bxdf closure instance might need to generate.
    // This value was set by the BxdfClosure constructor.
    const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

    // The response weight arrays: declare one for each response.
    // These will point to either a diffuse, specular, or user LPE classified
    // reservoir created inside the RixBXLobeWeights (lobeWeights) struct.
    // This mechanism uses the response's sg_*_LS value to determine which
    // reservoir to put the response in. See RixBxdfLobe.h
    RtColorRGB*  HGScatterWeight = NULL;

    // For the indexed sample in the shading context, evaluate the response nSamps times.
    for( int i=0; i < nSamps; ++i )
    {
        // Initialize this shaded point to an invalid generated lobe.
        lobesEvaluated[i].SetNone();

        // Determine whether we need to evaluate a sample direction for this i by
        // intersecting the responses we're supposed to produce (bxdfLobes: which
        // is based on what the renderer wants & what this shader produces), with
        // what responses the integrator wants at this shading point i.
        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted;

        // Do we need to compute the "HGScatter" response for this point?
        const bool  doHGScatter = ( lobesToConsider & sg_HGScatter_LT ).HasAny();
        if( doHGScatter )
        {
            // Create the sample weight array for this response, if needed.
            if( !HGScatterWeight ) HGScatterWeight = lobeWeights.AddActiveLobe( sg_HGScatter_LS );

            // Initialize any data needed to compute the response.
            const RtColorRGB Cs = Color[scIndex]*Gain[scIndex];
            const float       g = Direction[i];
            const RtNormal3  wn = Nn[scIndex];
            const RtVector3  wo = Vn[scIndex];

            // Evaluate the response sample.
            // Scale G: |g| ⇒ 1 is not very usable.
            if( Evaluate( Cs, g, wn, wo, wi[i], fPdf[i], rPdf[i], HGScatterWeight[i] ))
            {
                // If successful, set the i sample to the response we evaluated.
                lobesEvaluated[i] = sg_HGScatter_LT;
            }
        }
    }
}


//
//  Provides the integrator with a shading context set of samples
//  that are generated from this bxdf's response lobe(s).
//
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
    // Generate random number pairs for use in computing sample directions.
    RtFloat2*  xi = static_cast< RtFloat2* >( RixAlloca( numPts*sizeof(RtFloat2) ));
    rixRng->DrawSamples2D( xi );

    // Get all the responses this bxdf closure instance might need to generate.
    // This value was set by the BxdfClosure constructor.
    const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

    // The response weight arrays: declare one for each response.
    // These will point to either a diffuse, specular, or user LPE classified
    // reservoir created inside the RixBXLobeWeights (lobeWeights) struct.
    // This mechanism uses the response's sg_*_LS value to determine which
    // reservoir to put the response in. See RixBxdfLobe.h
    RtColorRGB*  HGScatterWeight = NULL;

    // For each point in the shading context, generate a sample (wi) if necessary.
    for( int i=0; i < numPts; ++i )
    {
        // Initialize this shaded point to an invalid generated response.
        lobeGenerated[i].SetValid( false );

        // Determine whether we need to generate a sample direction for this i by
        // intersecting the responses we're supposed to produce (bxdfLobes: which
        // is based on what the renderer wants & what this shader produces), with
        // the responses the integrator wants at this shading point i.
        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        // Do we need to compute the "HGScatter" response for this point?
        const bool  doHGScatter = ( lobesToConsider & sg_HGScatter_LT ).HasAny();
        if( doHGScatter )
        {
            // Create the sample weight array for this response, if needed.
            if( !HGScatterWeight ) HGScatterWeight = lobeWeights.AddActiveLobe( sg_HGScatter_LS );

            // Initialize any data needed to compute the response.
            const RtColorRGB Cs = Color[i]*Gain[i];
            const float       g = Direction[i];
            const RtNormal3  wn = Nn[i];
            const RtVector3  wo = Vn[i];

            // Generate the response sample.
            // Scale G: |g| ⇒ 1 is not very usable.
            if( Generate( Cs, g, wn, wo, xi[i], wi[i], fPdf[i], rPdf[i], HGScatterWeight[i] ))
            {
                // If successful, set the i sample to the response we generated.
                lobeGenerated[i] = sg_HGScatter_LS;
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

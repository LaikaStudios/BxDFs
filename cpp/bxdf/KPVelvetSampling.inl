/*
 *  Copyright 2022-2023 LAIKA. Authored by Mitch Prater.
 *
 *  Licensed under the Apache License Version 2.0 http://apache.org/licenses/LICENSE-2.0,
 *  or the MIT license http://opensource.org/licenses/MIT, at your option.
 *
 *  This program may not be copied, modified, or distributed except according to those terms.
 */
#include "prmanapi.h"
#include "RiTypesHelper.h"

#include "RixBxdf.h"
#include "RixBxdfLobe.h"
#include "RixShading.h"
#include "RixShadingUtils.h"


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
 *      wg - the geometric (modeled surface) normal.
 *      ws - surface shading normal. possibly "bumped" relative to wg.
 *      wn - response computation normal. generally equivalent to ws.
 *      Cs - the substance's characteristic (un-lit) coloration (spectrum).
 *      Cr - the response color & magnitude (spectrum*intensity).
 *      W  - the response weight.
 *      fPdf - forward pdf: probability of light moving from wi toward wo.
 *      rPdf - reverse pdf: probability of light moving from wo toward wi.
 */


//===============================================================
//  Koenderink-Pont asperity response.
//  Machine Vision and Applications, Vol. 14, p. 260-268 (2003)
//  http://dx.doi.org/10.1007/s00138-002-0089-7
//===============================================================
PRMAN_INLINE
void KoenderinkPontPdf
(
    const RtFloat        lm, // ∆/λ - fiber length / scattering mean free path.
    const RtVector3      wn, // n
    const RtVector3      wo, // u
    const RtVector3      wi, // v
    const RtFloat  cos_wnwo, // u·n
    const RtFloat  cos_wnwi, // v·n
    const RtFloat  cos_wowi, // p(-u·v)
    // Results:
    RtFloat&  fPdf,
    RtFloat&  rPdf,
    RtFloat&  W
)
{
    const RtVector3  u_add_v = wo + wi; // u+v
    const float  u_add_v_dot_n = u_add_v.Dot(wn); // (u+v)·n

    const float  dot_uv = cos_wowi; // p(-u·v)
    const float  dot_un = cos_wnwo; // u·n
    const float  dot_vn = cos_wnwi; // v·n

    const float  pdf = dot_uv*( dot_vn / u_add_v_dot_n )*( 1.0f - std::exp( -lm * u_add_v_dot_n/( dot_un*dot_vn )));

    rPdf = fPdf = W = pdf;
}

PRMAN_INLINE
bool KoenderinkPontEval
(
    const RtColorRGB     Cs,
    const RtFloat        lm, // ∆/λ - fiber length / scattering mean free path.
    const RtVector3      wn, // n
    const RtVector3      wo, // u
    const RtVector3      wi, // v
    const RtFloat  cos_wnwo, // u·n
    const RtFloat  cos_wnwi, // v·n
    // Results:
    RtFloat&  fPdf,
    RtFloat&  rPdf,
    RtColorRGB& Cr
)
{
    // Test the incident view factor (observer view factor tested prior to calling).
    if( cos_wnwi < 0.00001f ) return false;

    // Eliminates forward scattering.
    const float  cos_wowi = std::max( 0.0f, -(wo.Dot(wi)) ); // p(-u·v)

    RtFloat  W;
    KoenderinkPontPdf( lm, wn, wo, wi, cos_wnwo, cos_wnwi, cos_wowi, fPdf, rPdf, W );
    Cr = Cs * W;

    return true;
}

PRMAN_INLINE
bool KPEvaluate
(
    const RtColorRGB  Cs,
    const RtFloat     lm, // ∆/λ - fiber length / scattering mean free path.
    const RtNormal3   wn, // n
    const RtVector3   wo, // u
    const RtVector3   wi, // v
    // Results:
    RtFloat&  fPdf,
    RtFloat&  rPdf,
    RtColorRGB& Cr
)
{
    // Test the observer view factor.
    const float  cos_wnwo = wn.Dot(wo);  // u·n
    if( cos_wnwo < 0.00001f ) return false;

    const float  cos_wnwi = wn.Dot(wi); // v·n

    return KoenderinkPontEval( Cs, lm, wn, wo, wi, cos_wnwo, cos_wnwi, fPdf, rPdf, Cr );
}

PRMAN_INLINE
bool KPGenerate
(
    const RtColorRGB  Cs,
    const RtFloat     lm, // ∆/λ - fiber length / scattering mean free path.
    const RtNormal3   wn, // n
    const RtVector3   wo, // u
    const RtFloat2    xi,
    // Results:
    RtVector3&  wi,
    RtFloat&  fPdf,
    RtFloat&  rPdf,
    RtColorRGB& Cr
)
{
    // Test the observer view factor.
    const float  cos_wnwo = wn.Dot(wo);
    if( cos_wnwo < 0.00001f ) return false;

    // Generate an incident sample direction (wi) and its view factor.
    // XXX stub-in uniform samples for now.
    RtVector3  wt, wb;
    wn.CreateOrthonormalBasis( wt, wb );

    RtFloat  cos_wnwi;
    RixUniformDirectionalDistribution( xi, wn, wt, wb, wi, cos_wnwi );

    return KoenderinkPontEval( Cs, lm, wn, wo, wi, cos_wnwo, cos_wnwi, fPdf, rPdf, Cr );
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
    RtFloat*                fPdf, // forward pdf.
    RtFloat*                rPdf, // reverse pdf.
    RtColorRGB*             compTrans // compositing transparency (PRMan "Oi")?
)
{
    RtFloat2*  xi = static_cast< RtFloat2* >( RixAlloca( numPts*sizeof(RtFloat2) ));
    rixRng->DrawSamples2D( xi );

    const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

    RtColorRGB*  KPVelvetWeight = NULL;

    for( int i=0; i < numPts; i++ )
    {
        lobeGenerated[i].SetValid( false );

        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        const bool  doKPVelvet = ( lobesToConsider & sg_KPVelvet_LT ).HasAny();
        if( doKPVelvet )
        {
            // Create the sample weight array for this response, if needed.
            if( !KPVelvetWeight ) KPVelvetWeight = lobeWeights.AddActiveLobe( sg_KPVelvet_LS );

            // Initialize any data needed to compute the response.
            const RtColorRGB Cs = Color[i]*Gain[i];
            const RtFloat    lm = Length[i];
            const RtNormal3  wn = Orientation[i];
            const RtVector3  wo = Vn[i];

            if( KPGenerate( Cs, lm, wn, wo, xi[i], wi[i], fPdf[i], rPdf[i], KPVelvetWeight[i] ))
            {
                lobeGenerated[i] = sg_KPVelvet_LS;
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
    RtFloat*                fPdf, // forward pdf.
    RtFloat*                rPdf  // reverse pdf.
)
{
    const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

    RtColorRGB*  KPVelvetWeight = NULL;

    for( int i=0; i < numPts; i++ )
    {
        lobesEvaluated[i].SetNone();

        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted[i];

        const bool  doKPVelvet = ( lobesToConsider & sg_KPVelvet_LT ).HasAny();
        if( doKPVelvet )
        {
            // Initialize any data needed to compute the response.
            const RtColorRGB Cs = Color[i]*Gain[i];
            const RtFloat    lm = Length[i];
            const RtNormal3  wn = Orientation[i];
            const RtVector3  wo = Vn[i];

            if( !KPVelvetWeight ) KPVelvetWeight = lobeWeights.AddActiveLobe( sg_KPVelvet_LS );

            if( KPEvaluate( Cs, lm, wn, wo, wi[i], fPdf[i], rPdf[i], KPVelvetWeight[i] ))
            {
                lobesEvaluated[i] = sg_KPVelvet_LT;
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
    RtFloat*                fPdf, // forward pdf.
    RtFloat*                rPdf  // reverse pdf.
)
{
    const RixBXLobeTraits  bxdfLobes = GetAllLobeTraits();

    RtColorRGB*  KPVelvetWeight = NULL;

    for( int i=0; i < nSamps; i++ )
    {
        lobesEvaluated[i].SetNone();

        const RixBXLobeTraits  lobesToConsider = bxdfLobes & lobesWanted;

        const bool  doKPVelvet = ( lobesToConsider & sg_KPVelvet_LT ).HasAny();
        if( doKPVelvet )
        {
            if( !KPVelvetWeight ) KPVelvetWeight = lobeWeights.AddActiveLobe( sg_KPVelvet_LS );

            const RtColorRGB Cs = Color[scIndex]*Gain[scIndex];
            const RtFloat    lm = Length[scIndex];
            const RtNormal3  wn = Orientation[scIndex];
            const RtVector3  wo = Vn[scIndex];

            if( KPEvaluate( Cs, lm, wn, wo, wi[i], fPdf[i], rPdf[i], KPVelvetWeight[i] ))
            {
                lobesEvaluated[i] = sg_KPVelvet_LT;
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

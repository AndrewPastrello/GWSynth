#pragma once

#include "dpc_common.hpp"
#include "gamma.hpp"

#ifndef LAL_E
#define LAL_E 2.718281828459045235360287471352662498
#endif

/**
 * This function calculates the EOB A function which using the pre-computed
 * coefficients which should already have been calculated.
 */
static REAL8 XLALCalculateEOBA(
    const REAL8 r, /**<< Orbital separation (in units of total mass M) */
    EOBACoefficients
        *coeffs /**<< Pre-computed coefficients for the A function */
) {
  REAL8 r2, r3, r4, r5;
  REAL8 NA, DA;

  /* Note that this function uses pre-computed coefficients,
   * and assumes they have been calculated. Since this is a static function,
   * so only used here, I assume it is okay to neglect error checking
   */

  r2 = r * r;
  r3 = r2 * r;
  r4 = r2 * r2;
  r5 = r4 * r;

  NA = r4 * coeffs->n4 + r5 * coeffs->n5;

  DA = coeffs->d0 + r * coeffs->d1 + r2 * coeffs->d2 + r3 * coeffs->d3 +
       r4 * coeffs->d4 + r5 * coeffs->d5;

  return NA / DA;
}

/**
 * Function to calculate the EOB effective Hamiltonian for the
 * given values of the dynamical variables. The coefficients in the
 * A potential function should already have been computed.
 * Note that the pr used here is the tortoise co-ordinate.
 */
static REAL8 XLALEffectiveHamiltonian(
    const REAL8 eta,          /**<< Symmetric mass ratio */
    const REAL8 r,            /**<< Orbital separation */
    const REAL8 pr,           /**<< Tortoise co-ordinate */
    const REAL8 pp,           /**<< Momentum pphi */
    EOBACoefficients *aCoeffs /**<< Pre-computed coefficients in A function
                               */
) {
  /* The pr used in here is the tortoise co-ordinate */
  REAL8 r2, pr2, pp2, z3, eoba;

  r2 = r * r;
  pr2 = pr * pr;
  pp2 = pp * pp;

  eoba = XLALCalculateEOBA(r, aCoeffs);
  z3 = 2. * (4. - 3. * eta) * eta;
  return sycl::sqrt(pr2 + eoba * (1. + pp2 / r2 + z3 * pr2 * pr2 / r2));
}

/**
 * Calculated the derivative of the EOB A function with respect to
 * r, using the pre-computed A coefficients
 */
static REAL8 XLALCalculateEOBdAdr(
    const REAL8 r, /**<< Orbital separation (in units of total mass M) */
    EOBACoefficients
        *coeffs /**<< Pre-computed coefficients for the A function */
) {
  REAL8 r2, r3, r4, r5;

  REAL8 NA, DA, dNA, dDA, dA;

  r2 = r * r;
  r3 = r2 * r;
  r4 = r2 * r2;
  r5 = r4 * r;

  NA = r4 * coeffs->n4 + r5 * coeffs->n5;

  DA = coeffs->d0 + r * coeffs->d1 + r2 * coeffs->d2 + r3 * coeffs->d3 +
       r4 * coeffs->d4 + r5 * coeffs->d5;

  dNA = 4. * coeffs->n4 * r3 + 5. * coeffs->n5 * r4;

  dDA = coeffs->d1 + 2. * coeffs->d2 * r + 3. * coeffs->d3 * r2 +
        4. * coeffs->d4 * r3 + 5. * coeffs->d5 * r4;

  dA = dNA * DA - dDA * NA;

  return dA / (DA * DA);
}

/**
 * Computes the non-Keplerian correction to the velocity as determined from the
 * frequency obtained assuming a circular orbit. In the early stages of the
 * evolution, this should be a number close to 1.
 */
static REAL8 nonKeplerianCoefficient(
    const REAL8 *values,     /**<< Dynamics r, phi, pr, pphi */
    const REAL8 eta,         /**<< Symmetric mass ratio */
    EOBACoefficients *coeffs /**<< Pre-computed A coefficients */
) {
  REAL8 r = values[0];
  REAL8 pphi = values[3];

  REAL8 A = XLALCalculateEOBA(r, coeffs);
  REAL8 dA = XLALCalculateEOBdAdr(r, coeffs);

  return 2. *
         (1. +
          2. * eta * (-1. + sycl::sqrt((1. + pphi * pphi / (r * r)) * A))) /
         (r * r * dA);
}

static REAL8 XLALAssociatedLegendreXIsZero(const int l, const int m) {
  REAL8 legendre;

  /* we will switch on the values of m and n */
  switch (l) {
    case 1:
      switch (m) {
        case 1:
          legendre = -1.;
          break;
      }
      break;
    case 2:
      switch (m) {
        case 2:
          legendre = 3.;
          break;
        case 1:
          legendre = 0.;
          break;
      }
      break;
    case 3:
      switch (m) {
        case 3:
          legendre = -15.;
          break;
        case 2:
          legendre = 0.;
          break;
        case 1:
          legendre = 1.5;
          break;
      }
      break;
    case 4:
      switch (m) {
        case 4:
          legendre = 105.;
          break;
        case 3:
          legendre = 0.;
          break;
        case 2:
          legendre = -7.5;
          break;
        case 1:
          legendre = 0;
          break;
      }
      break;
    case 5:
      switch (m) {
        case 5:
          legendre = -945.;
          break;
        case 4:
          legendre = 0.;
          break;
        case 3:
          legendre = 52.5;
          break;
        case 2:
          legendre = 0;
          break;
        case 1:
          legendre = -1.875;
          break;
      }
      break;
    case 6:
      switch (m) {
        case 6:
          legendre = 10395.;
          break;
        case 5:
          legendre = 0.;
          break;
        case 4:
          legendre = -472.5;
          break;
        case 3:
          legendre = 0;
          break;
        case 2:
          legendre = 13.125;
          break;
        case 1:
          legendre = 0;
          break;
      }
      break;
    case 7:
      switch (m) {
        case 7:
          legendre = -135135.;
          break;
        case 6:
          legendre = 0.;
          break;
        case 5:
          legendre = 5197.5;
          break;
        case 4:
          legendre = 0.;
          break;
        case 3:
          legendre = -118.125;
          break;
        case 2:
          legendre = 0.;
          break;
        case 1:
          legendre = 2.1875;
          break;
      }
      break;
    case 8:
      switch (m) {
        case 8:
          legendre = 2027025.;
          break;
        case 7:
          legendre = 0.;
          break;
        case 6:
          legendre = -67567.5;
          break;
        case 5:
          legendre = 0.;
          break;
        case 4:
          legendre = 1299.375;
          break;
        case 3:
          legendre = 0.;
          break;
        case 2:
          legendre = -19.6875;
          break;
        case 1:
          legendre = 0.;
          break;
      }
      break;
  }

  legendre *= sycl::sqrt((REAL8)(2 * l + 1) * factorial(l - m) /
                         (4. * LAL_PI * factorial(l + m)));

  return legendre;
}

/* In the calculation of the Newtonian multipole, we only use
 * the spherical harmonic with theta set to pi/2. Since this
 * is always the case, we can use this information to use a
 * faster version of the spherical harmonic code
 */

static std::complex<REAL8> XLALScalarSphHarmThetaPiBy2(INT4 l, INT4 m,
                                                       REAL8 phi) {
  REAL8 legendre;
  INT4 absM = sycl::abs(m);

  /* For some reason GSL will not take negative m */
  /* We will have to use the relation between sph harmonics of +ve and -ve m */
  legendre = XLALAssociatedLegendreXIsZero(l, absM);

  /* Compute the values for the spherical harmonic */
  // std::complex<REAL8> y =
  //     crect(legendre * sycl::cos(m * phi), legendre * sycl::sin(m * phi));

  std::complex<REAL8> y = {legendre * sycl::cos(m * phi),
                           legendre * sycl::sin(m * phi)};

  /* If m is negative, perform some jiggery-pokery */
  if (m < 0 && absM % 2 == 1) {
    y = -(y);
  }

  return y;
}

/**
 * This function calculates the Newtonian multipole part of the
 * factorized waveform for the EOBNRv2 model. This is defined in Eq. 4.
 */
UNUSED static std::complex<REAL8> XLALSimIMREOBCalculateNewtonianMultipole(
    REAL8 x,          /**<< Dimensionless parameter \f$\equiv v^2\f$ */
    UNUSED REAL8 r,   /**<< Orbital separation (units of total mass M) */
    REAL8 phi,        /**<< Orbital phase (in radians) */
    UINT4 l,          /**<< Mode l */
    INT4 m,           /**<< Mode m */
    EOBParams *params /**<< Pre-computed coefficients, parameters, etc. */
) {
  INT4 epsilon = (l + m) % 2;

  /* Calculate the necessary Ylm */
  std::complex<REAL8> y = XLALScalarSphHarmThetaPiBy2(l - epsilon, -m, phi);

  /* Special treatment for (2,1) and (4,4) modes, defined in Eq. 17ab of
   * PRD84:124052 2011 */
  std::complex<REAL8> multipole;
  if ((l == 4 && m == 4) || (l == 2 && m == 1)) {
    multipole = params->prefixes->values[l][m] *
                sycl::pow(x, (REAL8)(l + epsilon) / 2.0 - 1.0) / r;
  } else {
    multipole = params->prefixes->values[l][m] *
                sycl::pow(x, (REAL8)(l + epsilon) / 2.0);
  }
  multipole *= y;

  return multipole;
}

/**
 * Computes the factorized waveform according to the prescription
 * given in Pan et al, arXiv:1106.1021v1 [gr-qc], for a given
 * mode l,m, for the given values of the dynamics at that point.
 * The function returns XLAL_SUCCESS if everything works out properly,
 * otherwise XLAL_FAILURE will be returned.
 */
UNUSED static std::complex<REAL8> XLALSimIMREOBGetFactorizedWaveform(
    const REAL8 *values, /**<< Vector containing dynamics r, phi, pr,
                              pphi for a given point */
    const REAL8 v,       /**<< Velocity (in geometric units) */
    const INT4 l,        /**<< Mode l */
    const INT4 m,        /**<< Mode m */
    EOBParams
        *params /**<< Structure containing pre-computed coefficients, etc. */
) {
  /* Status of function calls */
  INT4 i;

  REAL8 eta;
  REAL8 r, pr, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;
  REAL8 Hreal, Heff, Slm, deltalm, rholm, rholmPwrl;

  REAL8 z2;

  /* Non-Keplerian velocity */
  REAL8 vPhi;

  /* Pre-computed coefficients */
  FacWaveformCoeffs *hCoeffs = params->hCoeffs;

  if (sycl::abs(m) > (INT4)l) {
    // XLAL_ERROR(XLAL_EINVAL);
  }

  eta = params->eta;

  /* Check our eta was sensible */
  if (eta > 0.25) {
    // XLALPrintError("Eta seems to be > 0.25 - this isn't allowed!\n");
    // XLAL_ERROR(XLAL_EINVAL);
  } else if (eta == 0.25 && m % 2) {
    /* If m is odd and dM = 0, hLM will be zero */
    return {0, 0};
  }

  r = values[0];
  pr = values[2];
  pp = values[3];

  Heff = XLALEffectiveHamiltonian(eta, r, pr, pp, params->aCoeffs);
  Hreal = sycl::sqrt(1.0 + 2.0 * eta * (Heff - 1.0));
  v2 = v * v;
  Omega = v2 * v;
  vh3 = Hreal * Omega;
  vh = sycl::cbrt(vh3);
  eulerlogxabs = LAL_GAMMA + sycl::log(2.0 * (REAL8)m * v);

  /* Calculate the non-Keplerian velocity */
  /* given by Eq. (18) of Pan et al, PRD84, 124052(2011) */
  /* psi given by Eq. (19) of Pan et al, PRD84, 124052(2011) */
  /* Assign temporarily to vPhi */
  vPhi = nonKeplerianCoefficient(values, eta, params->aCoeffs);
  /* Assign rOmega value temporarily to vPhi */
  vPhi = r * sycl::cbrt(vPhi);
  /* Assign rOmega * Omega to vPhi */
  vPhi *= Omega;

  /* Calculate the newtonian multipole */
  std::complex<REAL8> hNewton = XLALSimIMREOBCalculateNewtonianMultipole(
      vPhi * vPhi, vPhi / Omega, values[1], (UINT4)l, m, params);

  /* Calculate the source term */
  if (((l + m) % 2) == 0) {
    Slm = Heff;
  } else {
    Slm = v * pp;
  }

  /* Calculate the Tail term */
  k = m * Omega;
  hathatk = Hreal * k;

  auto [lnr1, arg1] = lngamma_lanczos_complex(l + 1.0, -2.0 * hathatk);

  z2 = factorial(l);

  // std::complex<REAL8> Tlm =
  //     std::exp((lnr1 + LAL_PI * hathatk) +
  //              I * (arg1 + 2.0 * hathatk * sycl::log(4.0 * k /
  //              sycl::sqrt(LAL_E))));

  std::complex<REAL8> Tlm = complex_exp(
      (lnr1 + LAL_PI * hathatk) +
      I * (arg1 + 2.0 * hathatk * sycl::log(4.0 * k / sycl::sqrt(LAL_E))));

  Tlm /= z2;

  /* Calculate the residue phase and amplitude terms */
  switch (l) {
    case 2:
      switch (sycl::abs(m)) {
        case 2:
          deltalm = vh3 * (hCoeffs->delta22vh3 +
                           vh3 * (hCoeffs->delta22vh6 +
                                  vh * vh * (hCoeffs->delta22vh9 * vh))) +
                    hCoeffs->delta22v5 * v * v2 * v2 +
                    hCoeffs->delta22v8 * v2 * v2 * v2 * v2;
          rholm =
              1. + v2 * (hCoeffs->rho22v2 +
                         v * (hCoeffs->rho22v3 +
                              v * (hCoeffs->rho22v4 +
                                   v * (hCoeffs->rho22v5 +
                                        v * (hCoeffs->rho22v6 +
                                             hCoeffs->rho22v6l * eulerlogxabs +
                                             v * (hCoeffs->rho22v7 +
                                                  v * (hCoeffs->rho22v8 +
                                                       hCoeffs->rho22v8l *
                                                           eulerlogxabs +
                                                       (hCoeffs->rho22v10 +
                                                        hCoeffs->rho22v10l *
                                                            eulerlogxabs) *
                                                           v2)))))));
          break;
        case 1:
          deltalm = vh3 * (hCoeffs->delta21vh3 +
                           vh3 * (hCoeffs->delta21vh6 +
                                  vh * (hCoeffs->delta21vh7 +
                                        (hCoeffs->delta21vh9) * vh * vh))) +
                    hCoeffs->delta21v5 * v * v2 * v2 +
                    hCoeffs->delta21v7 * v2 * v2 * v2 * v;
          rholm =
              1. +
              v * (hCoeffs->rho21v1 +
                   v * (hCoeffs->rho21v2 +
                        v * (hCoeffs->rho21v3 +
                             v * (hCoeffs->rho21v4 +
                                  v * (hCoeffs->rho21v5 +
                                       v * (hCoeffs->rho21v6 +
                                            hCoeffs->rho21v6l * eulerlogxabs +
                                            v * (hCoeffs->rho21v7 +
                                                 hCoeffs->rho21v7l *
                                                     eulerlogxabs +
                                                 v * (hCoeffs->rho21v8 +
                                                      hCoeffs->rho21v8l *
                                                          eulerlogxabs +
                                                      (hCoeffs->rho21v10 +
                                                       hCoeffs->rho21v10l *
                                                           eulerlogxabs) *
                                                          v2))))))));
          break;
        default:
          //   XLAL_ERROR(XLAL_EINVAL);
          break;
      }
      break;
    case 3:
      switch (m) {
        case 3:
          deltalm =
              vh3 * (hCoeffs->delta33vh3 +
                     vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
              hCoeffs->delta33v5 * v * v2 * v2 +
              hCoeffs->delta33v7 * v2 * v2 * v2 * v;
          rholm =
              1. +
              v2 *
                  (hCoeffs->rho33v2 +
                   v * (hCoeffs->rho33v3 +
                        v * (hCoeffs->rho33v4 +
                             v * (hCoeffs->rho33v5 +
                                  v * (hCoeffs->rho33v6 +
                                       hCoeffs->rho33v6l * eulerlogxabs +
                                       v * (hCoeffs->rho33v7 +
                                            (hCoeffs->rho33v8 +
                                             hCoeffs->rho33v8l * eulerlogxabs) *
                                                v))))));
          break;
        case 2:
          deltalm =
              vh3 *
              (hCoeffs->delta32vh3 +
               vh * (hCoeffs->delta32vh4 +
                     vh * vh *
                         (hCoeffs->delta32vh6 + hCoeffs->delta32vh9 * vh3)));
          rholm =
              1. +
              v * (hCoeffs->rho32v +
                   v * (hCoeffs->rho32v2 +
                        v * (hCoeffs->rho32v3 +
                             v * (hCoeffs->rho32v4 +
                                  v * (hCoeffs->rho32v5 +
                                       v * (hCoeffs->rho32v6 +
                                            hCoeffs->rho32v6l * eulerlogxabs +
                                            (hCoeffs->rho32v8 +
                                             hCoeffs->rho32v8l * eulerlogxabs) *
                                                v2))))));
          break;
        case 1:
          deltalm = vh3 * (hCoeffs->delta31vh3 +
                           vh3 * (hCoeffs->delta31vh6 +
                                  vh * (hCoeffs->delta31vh7 +
                                        hCoeffs->delta31vh9 * vh * vh))) +
                    hCoeffs->delta31v5 * v * v2 * v2;
          rholm =
              1. +
              v2 *
                  (hCoeffs->rho31v2 +
                   v * (hCoeffs->rho31v3 +
                        v * (hCoeffs->rho31v4 +
                             v * (hCoeffs->rho31v5 +
                                  v * (hCoeffs->rho31v6 +
                                       hCoeffs->rho31v6l * eulerlogxabs +
                                       v * (hCoeffs->rho31v7 +
                                            (hCoeffs->rho31v8 +
                                             hCoeffs->rho31v8l * eulerlogxabs) *
                                                v))))));
          break;
        default:
          //   XLAL_ERROR(XLAL_EINVAL);
          break;
      }
      break;
    case 4:
      switch (m) {
        case 4:
          deltalm = vh3 * (hCoeffs->delta44vh3 + hCoeffs->delta44vh6 * vh3) +
                    hCoeffs->delta44v5 * v2 * v2 * v;
          rholm = 1. + v2 * (hCoeffs->rho44v2 +
                             v * (hCoeffs->rho44v3 +
                                  v * (hCoeffs->rho44v4 +
                                       v * (hCoeffs->rho44v5 +
                                            (hCoeffs->rho44v6 +
                                             hCoeffs->rho44v6l * eulerlogxabs) *
                                                v))));
          break;
        case 3:
          deltalm =
              vh3 *
              (hCoeffs->delta43vh3 +
               vh * (hCoeffs->delta43vh4 + hCoeffs->delta43vh6 * vh * vh));
          rholm = 1. + v * (hCoeffs->rho43v +
                            v * (hCoeffs->rho43v2 +
                                 v2 * (hCoeffs->rho43v4 +
                                       v * (hCoeffs->rho43v5 +
                                            (hCoeffs->rho43v6 +
                                             hCoeffs->rho43v6l * eulerlogxabs) *
                                                v))));
          break;
        case 2:
          deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
          rholm = 1. + v2 * (hCoeffs->rho42v2 +
                             v * (hCoeffs->rho42v3 +
                                  v * (hCoeffs->rho42v4 +
                                       v * (hCoeffs->rho42v5 +
                                            (hCoeffs->rho42v6 +
                                             hCoeffs->rho42v6l * eulerlogxabs) *
                                                v))));
          break;
        case 1:
          deltalm =
              vh3 *
              (hCoeffs->delta41vh3 +
               vh * (hCoeffs->delta41vh4 + hCoeffs->delta41vh6 * vh * vh));
          rholm = 1. + v * (hCoeffs->rho41v +
                            v * (hCoeffs->rho41v2 +
                                 v2 * (hCoeffs->rho41v4 +
                                       v * (hCoeffs->rho41v5 +
                                            (hCoeffs->rho41v6 +
                                             hCoeffs->rho41v6l * eulerlogxabs) *
                                                v))));
          break;
        default:
          //   XLAL_ERROR(XLAL_EINVAL);
          break;
      }
      break;
    case 5:
      switch (m) {
        case 5:
          deltalm =
              hCoeffs->delta55vh3 * vh3 + hCoeffs->delta55v5 * v2 * v2 * v;
          rholm =
              1. +
              v2 * (hCoeffs->rho55v2 +
                    v * (hCoeffs->rho55v3 +
                         v * (hCoeffs->rho55v4 +
                              v * (hCoeffs->rho55v5 + hCoeffs->rho55v6 * v))));
          break;
        case 4:
          deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
          rholm = 1. + v2 * (hCoeffs->rho54v2 +
                             v * (hCoeffs->rho54v3 + hCoeffs->rho54v4 * v));
          break;
        case 3:
          deltalm = hCoeffs->delta53vh3 * vh3;
          rholm =
              1. + v2 * (hCoeffs->rho53v2 +
                         v * (hCoeffs->rho53v3 +
                              v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
          break;
        case 2:
          deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
          rholm = 1. + v2 * (hCoeffs->rho52v2 +
                             v * (hCoeffs->rho52v3 + hCoeffs->rho52v4 * v));
          break;
        case 1:
          deltalm = hCoeffs->delta51vh3 * vh3;
          rholm =
              1. + v2 * (hCoeffs->rho51v2 +
                         v * (hCoeffs->rho51v3 +
                              v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)));
          break;
        default:
          //   XLAL_ERROR(XLAL_EINVAL);
          break;
      }
      break;
    case 6:
      switch (m) {
        case 6:
          deltalm = hCoeffs->delta66vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho66v2 +
                             v * (hCoeffs->rho66v3 + hCoeffs->rho66v4 * v));
          break;
        case 5:
          deltalm = hCoeffs->delta65vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
          break;
        case 4:
          deltalm = hCoeffs->delta64vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho64v2 +
                             v * (hCoeffs->rho64v3 + hCoeffs->rho64v4 * v));
          break;
        case 3:
          deltalm = hCoeffs->delta63vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
          break;
        case 2:
          deltalm = hCoeffs->delta62vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho62v2 +
                             v * (hCoeffs->rho62v3 + hCoeffs->rho62v4 * v));
          break;
        case 1:
          deltalm = hCoeffs->delta61vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
          break;
        default:
          //   XLAL_ERROR(XLAL_EINVAL);
          break;
      }
      break;
    case 7:
      switch (m) {
        case 7:
          deltalm = hCoeffs->delta77vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
          break;
        case 6:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho76v2 * v2;
          break;
        case 5:
          deltalm = hCoeffs->delta75vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
          break;
        case 4:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho74v2 * v2;
          break;
        case 3:
          deltalm = hCoeffs->delta73vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
          break;
        case 2:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho72v2 * v2;
          break;
        case 1:
          deltalm = hCoeffs->delta71vh3 * vh3;
          rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
          break;
        default:
          //   XLAL_ERROR(XLAL_EINVAL);
          break;
      }
      break;
    case 8:
      switch (m) {
        case 8:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho88v2 * v2;
          break;
        case 7:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho87v2 * v2;
          break;
        case 6:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho86v2 * v2;
          break;
        case 5:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho85v2 * v2;
          break;
        case 4:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho84v2 * v2;
          break;
        case 3:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho83v2 * v2;
          break;
        case 2:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho82v2 * v2;
          break;
        case 1:
          deltalm = 0.0;
          rholm = 1. + hCoeffs->rho81v2 * v2;
          break;
        default:
          //   XLAL_ERROR(XLAL_EINVAL);
          break;
      }
      break;
    default:
      //   XLAL_ERROR(XLAL_EINVAL);
      break;
  }

  /* Raise rholm to the lth power */
  rholmPwrl = 1.0;
  i = l;
  while (i--) {
    rholmPwrl *= rholm;
  }

  // std::complex<REAL8> hlm = Tlm * std::exp(I * deltalm) * Slm * rholmPwrl;
  std::complex<REAL8> hlm = Tlm * complex_exp(I * deltalm) * Slm * rholmPwrl;

  hlm *= hNewton;

  return hlm;
}

/**
 * This function calculates the non-quasicircular correction to apply to
 * the waveform. The form of this correction can be found in Pan et al,
 * PRD84, 124052(2011) [arXiv:1106.1021], Eq.(22), and also in the DCC
 document
 * T1100433. Note that when calling this function, the NQC coefficients
 should
 * already have been pre-computed.
 */
UNUSED static std::complex<REAL8> XLALSimIMREOBNonQCCorrection(
    const REAL8 *values,
    /**<< Dynamics r, phi, pr, pphi */
    const REAL8 omega, /**<< Angular frequency */
    EOBNonQCCoeffs *coeffs
    /**<< NQC coefficients */
) {
  REAL8 rOmega, rOmegaSq;
  REAL8 r, p, sqrtR;

  REAL8 mag, phase;

  r = values[0];
  p = values[2];

  sqrtR = sycl::sqrt(r);

  rOmega = r * omega;
  rOmegaSq = rOmega * rOmega;
  /*printf("a1 = %.16e, a2 = %.16e, a3 = %.16e, a3S = %.16e, a4 = %.16e, a5 =
  %.16e\n",coeffs->a1,coeffs->a2,coeffs->a3,coeffs->a3S,
  coeffs->a4,coeffs->a5); printf("b1 = %.16e, b2 = %.16e, b3 = %.16e, b4 =
  %.16e\n",coeffs->b1,coeffs->b2,coeffs->b3,coeffs->b4);*/
  /* In EOBNRv2, coeffs->a3S, coeffs->a4 and coeffs->a5 are set to zero */
  /* through XLALSimIMREOBGetCalibratedNQCCoeffs() */
  /* and XLALSimIMREOBCalculateNQCCoefficients() */
  mag = 1. + (p * p / rOmegaSq) *
                 (coeffs->a1 + coeffs->a2 / r +
                  (coeffs->a3 + coeffs->a3S) / (r * sqrtR) +
                  coeffs->a4 / (r * r) + coeffs->a5 / (r * r * sqrtR));
  // printf("NQC INFO mag = %.16e, r = %.16e, p = %.16e\n",mag,r,p);
  phase =
      coeffs->b1 * p / rOmega +
      p * p * p / rOmega * (coeffs->b2 + coeffs->b3 / sqrtR + coeffs->b4 / r);

  std::complex<REAL8> nqc = mag * sycl::cos(phase);
  nqc += I * mag * sycl::sin(phase);
  /*printf("r = %.16e, pr = %.16e, omega = %.16e\n",r,p,omega);
  printf("NQC mag = %.16e, arg = %.16e\n",mag,phase);*/
  return nqc;
}

/**
 * This function calculates the factorized flux in the EOB dynamics for
 * the EOBNR (and potentially subsequent) models. The flux function
 * is found in Phys.Rev.D79:064004,2009.
 */
static REAL8 XLALSimIMREOBFactorizedFlux(
    const REAL8 *values, /**<< Dynamics r, phi, pr, pphi */
    const REAL8 omega,   /**<< Angular frequency omega */
    EOBParams *ak,       /**<< Structure containing pre-computed parameters
                          */
    const INT4 lMax      /**<< Maximum l to include when calculating flux
      (between 2
                         and 8) */
)

{
  REAL8 flux = 0.0;
  REAL8 v;
  REAL8 omegaSq;
  INT4 l, m;

  EOBNonQCCoeffs *nqcCoeffs;

  // #ifndef LAL_NDEBUG
  //   if (!values || !ak) {
  //     XLAL_ERROR_REAL8(XLAL_EFAULT);
  //   }
  // #endif

  //   if (lMax < 2) {
  //     XLAL_ERROR_REAL8(XLAL_EINVAL);
  //   }

  nqcCoeffs = ak->nqcCoeffs;

  /* Omegs is the derivative of phi */
  omegaSq = omega * omega;

  v = sycl::cbrt(omega);

  /* We need to apply the NQC for the (2,2) mode */
  /* To avoid having an if statement in the loop we will */
  /* deal with (2,2) and (2,1) separately */
  /* (2,2) */
  l = 2;
  m = 2;

  std::complex<REAL8> hNQC =
      XLALSimIMREOBNonQCCorrection(values, omega, nqcCoeffs);

  std::complex<REAL8> hLM =
      XLALSimIMREOBGetFactorizedWaveform(values, v, l, m, ak);

  /* For the 2,2 mode, we apply NQC correction to the flux */
  hLM *= hNQC;

  flux = (REAL8)(m * m) * omegaSq *
         (std::real(hLM) * std::real(hLM) + std::imag(hLM) * std::imag(hLM));

  /* (2,1) */
  l = 2;
  m = 1;

  hLM = XLALSimIMREOBGetFactorizedWaveform(values, v, l, m, ak);

  flux += (REAL8)(m * m) * omegaSq *
          (std::real(hLM) * std::real(hLM) + std::imag(hLM) * std::imag(hLM));

  /* All other modes */
  for (l = 3; l <= lMax; l++) {
    /*INT4 minM = l-3;
    if ( minM < 1 )
      minM = 1;*/

    for (m = 1; m <= l; m++) {
      hLM = XLALSimIMREOBGetFactorizedWaveform(values, v, l, m, ak);

      flux +=
          (REAL8)(m * m) * omegaSq *
          (std::real(hLM) * std::real(hLM) + std::imag(hLM) * std::imag(hLM));
    }
  }

  return flux * LAL_1_PI / 8.0;
}

/**
 * Calculate the EOB D function.
 */
static REAL8 XLALCalculateEOBD(
    REAL8 r,  /**<< Orbital separation (in units of total mass M) */
    REAL8 eta /**<< Symmetric mass ratio */
) {
  REAL8 u, u2, u3;

  u = 1. / r;
  u2 = u * u;
  u3 = u2 * u;

  return 1. / (1. + 6. * eta * u2 + 2. * eta * (26. - 3. * eta) * u3);
}

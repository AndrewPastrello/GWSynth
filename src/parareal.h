#pragma once

#include <math.h>

#include "LALSimIMREOBNRv2.h"

#ifdef __cplusplus
#include "derivative.hpp"
#include "dpc_common.hpp"
#define pr_sqrt sycl::sqrt
#else
#include "LALSimIMREOBFactorizedFlux.c"
#define pr_sqrt sqrt
#endif

static const int nn = 4;

#ifdef __cplusplus
extern "C" {
#endif
int parareal(EOBParams *params, const REAL8 y_init[nn], const REAL8 t_init,
             REAL8 t_end, const REAL8 dt, const REAL8 pr_tol,
             const int n_slices, REAL8Array **y_out);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#define pr_max std::max
#else
#define pr_max fmax
#endif
/*
 * This function computes the derivatives of the EOB Hamiltonian w.r.t. the
 * dynamical variables, and therefore the derivatives of the dynamical variables
 * w.r.t. time. As such this gets called in the Runge-Kutta integration of the
 * orbit.
 */
static int LALHCapDerivativesP4PN(
    double UNUSED t, /**<< Current time (GSL requires it to be a parameter, but
                        it's irrelevant) */
    const REAL8 values[], /**<< The dynamics r, phi, pr, pphi */
    // REAL8 values[],
    // const REAL8 r, const REAL8 phi, const REAL8 pr, const REAL8 pphi,
    REAL8
        dvalues[], /**<< The derivatives dr/dt, dphi/dt. dpr/dt and dpphi/dt */
    void *funcParams /**<< Structure containing all the necessary parameters */
) {
  EOBParams *params = NULL;

  /* Max l to sum up to in the factorized flux */
  const INT4 lMax = 8;

  REAL8 eta, omega;

  REAL8 r2, p2, p4, q2;
  REAL8 u, u2, u3;
  REAL8 A, AoverSqrtD, dAdr, Heff, Hreal;
  REAL8 HeffHreal;
  REAL8 flux;
  REAL8 z3;

  params = (EOBParams *)funcParams;

  eta = params->eta;

  z3 = 2. * (4. - 3. * eta) * eta;

  REAL8 r, phi, pr, pphi;
  r = values[0];
  phi = values[1];
  pr = values[2];
  pphi = values[3];

  u = 1.0 / r;
  u2 = u * u;
  u3 = u2 * u;
  r2 = r * r;
  p2 = pr * pr;
  p4 = p2 * p2;
  q2 = pphi * pphi;

  A = XLALCalculateEOBA(r, params->aCoeffs);
  dAdr = XLALCalculateEOBdAdr(r, params->aCoeffs);
  AoverSqrtD = A / pr_sqrt(XLALCalculateEOBD(r, eta));

  /* Note that Hreal as given here is missing a factor of 1/eta */
  /* This is because it only enters into the derivatives in     */
  /* the combination eta*Hreal*Heff, so the eta would get       */
  /* cancelled out anyway. */

  Heff = XLALEffectiveHamiltonian(eta, r, pr, pphi, params->aCoeffs);
  Hreal = pr_sqrt(1. + 2. * eta * (Heff - 1.));

  HeffHreal = Heff * Hreal;

  /* rDot */
  dvalues[0] = AoverSqrtD * u2 * pr * (r2 + 2. * p2 * z3 * A) / HeffHreal;

  /* sDot */
  dvalues[1] = omega = pphi * A * u2 / HeffHreal;

  /* Note that the only field of dvalues used in the flux is dvalues->data[1] */
  /* which we have already calculated. */
  flux = XLALSimIMREOBFactorizedFlux(values, omega, params, lMax);

  /* pDot */
  dvalues[2] = 0.5 * AoverSqrtD * u3 *
                   (2.0 * (q2 + p4 * z3) * A - r * (q2 + r2 + p4 * z3) * dAdr) /
                   HeffHreal -
               (pr / pphi) * (flux / (eta * omega));

  /* qDot */
  dvalues[3] = -flux / (eta * omega);

  return GSL_SUCCESS;
}

/*
 * Function which will determine whether to stop the evolution for the initial,
 * user-requested sample rate. We stop in this case when we have reached the
 * peak orbital frequency.
 */
static int XLALFirstStoppingCondition(
    double UNUSED t,              /**<< Current time (required by GSL) */
    const double UNUSED values[], /**<< Current dynamics (required by GSL) */
    double dvalues[],             /**<< Derivatives of dynamics w.r.t. time */
    void *funcParams /**<< Structure containing necessary parameters */
) {
  EOBParams *params = (EOBParams *)funcParams;
  double omega = dvalues[1];

  if (omega < params->omega) {
    return 1;
  }

  params->omega = omega;
  return GSL_SUCCESS;
}

/*
 * Function which will determine whether to stop the evolution for the high
 * sample rate. In this case, the data obtained will be used to attach the
 * ringdown, so to make sure we won't be interpolating data too near the final
 * points, we push this integration as far as is feasible.
 */
static UNUSED int XLALHighSRStoppingCondition(
    double UNUSED t,        /**<< Current time (required by GSL) */
    const double values[],  /**<< Current dynamics  */
    double dvalues[],       /**<< Derivatives of dynamics w.r.t. time */
    void UNUSED *funcParams /**<< Structure containing necessary parameters */
) {
  EOBParams *params = (EOBParams *)funcParams;
  REAL8 rstop;
  if (params->eta > 0.1) {
    rstop = 1.25 - params->eta;
  } else {
    rstop = 2.1 - 10.0 * params->eta;
  }

  if (values[0] <= rstop || isnan(dvalues[3]) || isnan(dvalues[2]) ||
      isnan(dvalues[1]) || isnan(dvalues[0])) {
    return 1;
  }
  return 0;
}

static UNUSED int rk42(void *params, const REAL8 y_init[nn], const REAL8 t_init,
                       const REAL8 t_end, REAL8 dt, const REAL8 tol,
                       const int stop_on_t, REAL8 *t_out, REAL8 *y_out,
                       int *steps_out) {
  REAL8 t = t_init;

  REAL8 y[nn];
  REAL8 dydt[nn];
  REAL8 y2[nn];
  REAL8 y4[nn];
  REAL8 k1[nn];
  REAL8 k2[nn];
  REAL8 k3[nn];
  REAL8 k4[nn];
  REAL8 y_temp[nn];

  for (int i = 0; i < nn; i++) y[i] = y_init[i];

  t_out[0] = t_init;
  for (int i = 0; i < nn; i++) {
    y_out[i] = y_init[i];
  }

  int steps = 0;
  int stopped_on_condition = false;

  REAL8 m = 1.0;
  do {
    LALHCapDerivativesP4PN(0.0, y, k1, params);
    while (true) {
      dt *= 0.9 * m;

      if (stop_on_t && t + dt > t_end) dt = t_end - t;

      for (int i = 0; i < nn; i++) {
        y_temp[i] = y[i] + 0.5 * dt * k1[i];
      }

      LALHCapDerivativesP4PN(0.0, y_temp, k2, params);

      for (int i = 0; i < nn; i++) {
        y_temp[i] = y[i] + 0.5 * dt * k2[i];
      }

      LALHCapDerivativesP4PN(0.0, y_temp, k3, params);

      for (int i = 0; i < nn; i++) {
        y2[i] = y[i] + dt * k3[i];
      }

      LALHCapDerivativesP4PN(0.0, y2, k4, params);

      for (int i = 0; i < nn; i++) {
        dydt[i] = 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
      }

      for (int i = 0; i < nn; i++) {
        y4[i] = y[i] + dt * dydt[i];
      }

      REAL8 max_rerr = 0;
      for (int i = 0; i < nn; i++) {
        max_rerr = pr_max(max_rerr, fabs((y4[i] - y2[i]) / y4[i]));
      }

      m = pr_sqrt(tol / max_rerr);

      if (m >= 1.0) break;
    }
    for (int i = 0; i < nn; i++) y[i] = y4[i];

    t += dt;

    steps++;

    t_out[steps] = t;
    for (int i = 0; i < nn; i++) {
      y_out[steps * nn + i] = y[i];
    }

    if (stop_on_t && t >= t_end) break;
    stopped_on_condition =
        !(XLALFirstStoppingCondition(0.0, y, dydt, params) == GSL_SUCCESS);
  } while (!stopped_on_condition);
  *steps_out = steps;

  return stopped_on_condition;
}
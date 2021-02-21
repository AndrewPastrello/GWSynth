#include "parareal.h"
//
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <lal/AVFactories.h>
#include <lal/XLALGSL.h>

#include "dpc_common.hpp"

EOBParams copy(EOBParams *a, sycl::queue q) {
  EOBParams b = *a;
  b.aCoeffs = sycl::malloc<EOBACoefficients>(1, q, sycl::usm::alloc::shared);
  b.hCoeffs = sycl::malloc<FacWaveformCoeffs>(1, q, sycl::usm::alloc::shared);
  b.nqcCoeffs = sycl::malloc<EOBNonQCCoeffs>(1, q, sycl::usm::alloc::shared);
  b.prefixes =
      sycl::malloc<NewtonMultipolePrefixes>(1, q, sycl::usm::alloc::shared);

  q.memcpy(b.aCoeffs, a->aCoeffs, sizeof(EOBACoefficients));
  q.memcpy(b.hCoeffs, a->hCoeffs, sizeof(FacWaveformCoeffs));
  q.memcpy(b.nqcCoeffs, a->nqcCoeffs, sizeof(EOBNonQCCoeffs));
  q.memcpy(b.prefixes, a->prefixes, sizeof(NewtonMultipolePrefixes));
  q.wait();
  return b;
}

void free(EOBParams a, sycl::queue q) {
  sycl::free(a.aCoeffs, q);
  sycl::free(a.hCoeffs, q);
  sycl::free(a.nqcCoeffs, q);
  sycl::free(a.prefixes, q);
  q.wait();
}

extern "C" {
int parareal(EOBParams *params_init, const REAL8 y_init[nn], const REAL8 t_init,
             REAL8 t_end, const REAL8 dt, const REAL8 pr_tol,
             const int n_slices, REAL8Array **y_out) {
  // Phase cancels with tol 1e-7
  const REAL8 coarse_tol = 1e-5;
  const REAL8 fine_tol = 1e-7;

  sycl::queue q(sycl::cpu_selector{});

  auto params = copy(params_init, q);

  REAL8 coarse_dt = (t_end - t_init) / n_slices;

  int bufferlength = (int)((t_end - t_init) / dt) + 2;
  REAL8 *buf =
      sycl::malloc<REAL8>((nn + 1) * bufferlength, q, sycl::usm::alloc::shared);
  REAL8 *t_buf = buf;
  REAL8 *y_buf = buf + bufferlength;

  REAL8 *coarse_buf =
      sycl::malloc<REAL8>((nn + 1) * bufferlength, q, sycl::usm::alloc::shared);
  REAL8 *coarse_t_buf = coarse_buf;
  REAL8 *coarse_y_buf = coarse_buf + bufferlength;

  // Fine solve
  REAL8 *y_initf =
      sycl::malloc<REAL8>(n_slices * nn, q, sycl::usm::alloc::shared);
  REAL8 *y_initc =
      sycl::malloc<REAL8>(n_slices * nn, q, sycl::usm::alloc::shared);
  REAL8 *y_correct =
      sycl::malloc<REAL8>(n_slices * nn, q, sycl::usm::alloc::shared);
  REAL8 *y_fine0 =
      sycl::malloc<REAL8>(n_slices * nn, q, sycl::usm::alloc::shared);
  REAL8 *y_coarse0 =
      sycl::malloc<REAL8>(n_slices * nn, q, sycl::usm::alloc::shared);
  REAL8 *y_coarse1 =
      sycl::malloc<REAL8>(n_slices * nn, q, sycl::usm::alloc::shared);

  REAL8 *t_initf = sycl::malloc<REAL8>(n_slices, q, sycl::usm::alloc::shared);
  REAL8 *t_endf = sycl::malloc<REAL8>(n_slices, q, sycl::usm::alloc::shared);

  EOBParams coarse_params = params;

  EOBParams *fine_params =
      sycl::malloc<EOBParams>(n_slices, q, sycl::usm::alloc::shared);
  int *fine_steps = sycl::malloc<int>(n_slices, q, sycl::usm::alloc::shared);
  q.wait();

initialize:
  for (int i = 0; i < nn; i++) y_initc[i] = y_init[i];
  for (int k = 0; k < nn; k++) y_initf[k] = y_init[k];
  // TODO: Load balance fine step by reinitializing t_initf uniformly from
  // coarse_t_buf location, not uniform time
  for (int i = 0; i < n_slices; i++) {
    t_initf[i] = t_init + i * coarse_dt;
    t_endf[i] = t_initf[i] + coarse_dt;
    fine_params[i] = params;
  }

  // Initial coarse solve
  for (int i = 0; i < n_slices; i++) {
    int steps;
    int stop_on_time = i == n_slices - 1 ? false : true;
    int stopped_on_condition =
        rk42(&coarse_params, &y_initc[i * nn], t_initf[i], t_endf[i], dt,
             coarse_tol, stop_on_time, coarse_t_buf, coarse_y_buf, &steps);
    //  Reinitialize if the guessed t_end is not in the last slice
    if (stopped_on_condition && i < n_slices - 1) {
      t_end = coarse_t_buf[steps];
      coarse_dt = (t_end - t_init) / n_slices;
      coarse_params.omega = params.omega;
      goto initialize;
    }
    for (int k = 0; k < nn; k++) {
      y_coarse0[i * nn + k] = coarse_y_buf[steps * nn + k];
      y_coarse1[i * nn + k] = coarse_y_buf[steps * nn + k];
      if (i < n_slices - 1) {
        y_initc[(i + 1) * nn + k] = y_coarse1[i * nn + k];
        y_initf[(i + 1) * nn + k] = y_coarse1[i * nn + k];
      }
    }
  }
  coarse_params.omega = params.omega;

  // Parareal converges in max n_slices iterations
  for (int iteration = 0; iteration < n_slices; iteration++) {
    // Fine solve
    q.submit([&](auto &cgh) {
       cgh.parallel_for(sycl::range<1>(n_slices - iteration), [=](auto ind) {
         auto i = ind + iteration;
         int buf_offset = i * (bufferlength / n_slices);
         // Last thread stops on stopping condition, others stop on time
         int stop_on_time = i == n_slices - 1 ? false : true;
         rk42(&fine_params[i], &y_initf[i * nn], t_initf[i], t_endf[i], dt,
              fine_tol, stop_on_time, &t_buf[buf_offset],
              &y_buf[buf_offset * nn], &fine_steps[i]);
         for (int k = 0; k < nn; k++) {
           y_fine0[i * nn + k] = y_buf[(buf_offset + fine_steps[i]) * nn + k];
         }
       });
     }).wait();

    // Additional slice solved after each iteration
    if (iteration < n_slices - 1) {
      for (int k = 0; k < nn; k++) {
        y_initc[(iteration + 1) * nn + k] = y_fine0[iteration * nn + k];
        y_initf[(iteration + 1) * nn + k] = y_fine0[iteration * nn + k];
      }
    }
    // Coarse solve
    for (int i = iteration + 1; i < n_slices; i++) {
      int stop_on_time = i == n_slices - 1 ? false : true;
      int steps;
      rk42(&coarse_params, &y_initc[i * nn], t_initf[i], t_endf[i], dt,
           coarse_tol, stop_on_time, coarse_t_buf, coarse_y_buf, &steps);
      for (int k = 0; k < nn; k++) {
        y_coarse1[i * nn + k] = coarse_y_buf[steps * nn + k];
        if (i < n_slices - 1) {
          y_initc[(i + 1) * nn + k] = y_coarse1[i * nn + k];
        }
      }
    }

    // Correction
    for (int k = 0; k < nn; k++) {
      y_correct[iteration * nn + k] = y_fine0[iteration * nn + k];
    }
    for (int i = iteration + 1; i < n_slices; i++) {
      for (int k = 0; k < nn; k++) {
        int ind = i * nn + k;
        y_correct[ind] = y_coarse1[ind] + y_fine0[ind] - y_coarse0[ind];
      }
    }

    REAL8 max_rerr = 0;
    for (int i = iteration + 2; i < n_slices; i++) {
      // Last point of previous time slice
      REAL8 *y_new = y_correct + (i - 1) * nn;
      REAL8 *y_old = y_initc + i * nn;
      for (int k = 0; k < nn; k++) {
        REAL8 rerr = std::abs((y_old[k] - y_new[k]) / y_new[k]);
        max_rerr = std::max(max_rerr, rerr);
      }
    }
    if (max_rerr < pr_tol) break;

    // Reinitialize
    for (int i = iteration + 2; i < n_slices; i++) {
      for (int k = 0; k < nn; k++) {
        y_initc[i * nn + k] = y_correct[(i - 1) * nn + k];
        y_initf[i * nn + k] = y_coarse1[(i - 1) * nn + k];
      }
    }
    for (int i = iteration + 1; i < n_slices; i++) {
      fine_params[i].omega = params.omega;
    }
    REAL8 *tmp = y_coarse0;
    y_coarse0 = y_coarse1;
    y_coarse1 = tmp;
    coarse_params.omega = params.omega;
  }

  // Compact bufs
  int integrationlength = fine_steps[0] + 1;
  int t_ind = integrationlength;
  int y_ind = integrationlength * nn;
  for (int j = 1; j < n_slices; j++) {
    integrationlength += fine_steps[j];
    int buf_offset = j * (bufferlength / n_slices);
    // Skip initial value
    for (int i = 1; i <= fine_steps[j]; i++) {
      t_buf[t_ind++] = t_buf[buf_offset + i];
      for (int k = 0; k < nn; k++) {
        y_buf[y_ind++] = y_buf[(buf_offset + i) * nn + k];
      }
    }
  }

  REAL8 *y_flat = sycl::malloc<REAL8>((nn + 1) * integrationlength, q,
                                      sycl::usm::alloc::shared);
  for (int k = 0; k < nn; k++) {
    for (int i = 0; i < integrationlength; i++) {
      y_flat[integrationlength * k + i] = y_buf[i * nn + k];
    }
  }

  t_end = t_buf[integrationlength - 1];

  int outputlen = (int)(t_end / dt) + 1;
  REAL8Array *output = XLALCreateREAL8ArrayL(2, nn + 1, outputlen);

  /* needed for the final interpolation for unequal timesteps*/
  gsl_spline *interp = gsl_spline_alloc(gsl_interp_cspline, integrationlength);
  gsl_interp_accel *accel = gsl_interp_accel_alloc();

  for (int j = 0; j < outputlen; j++) output->data[j] = t_init + dt * j;

  /* interpolate! */
  for (int i = 1; i <= nn; i++) {
    gsl_spline_init(interp, t_buf, &y_flat[integrationlength * (i - 1)],
                    integrationlength);

    for (int j = 0; j < outputlen; j++) {
      gsl_spline_eval_e(interp, output->data[j], accel,
                        &(output->data[outputlen * i + j]));
    }
  }

  if (interp) XLAL_CALLGSL(gsl_spline_free(interp));
  if (accel) XLAL_CALLGSL(gsl_interp_accel_free(accel));

  if (buf) sycl::free(buf, q);
  if (coarse_buf) sycl::free(coarse_buf, q);
  if (y_flat) sycl::free(y_flat, q);
  if (fine_steps) sycl::free(fine_steps, q);
  if (t_initf) sycl::free(t_initf, q);
  if (t_endf) sycl::free(t_endf, q);
  if (y_initf) sycl::free(y_initf, q);
  if (y_initc) sycl::free(y_initc, q);
  if (y_correct) sycl::free(y_correct, q);
  if (y_fine0) sycl::free(y_fine0, q);
  if (y_coarse0) sycl::free(y_coarse0, q);
  if (y_coarse1) sycl::free(y_coarse1, q);
  if (fine_params) sycl::free(fine_params, q);
  free(params, q);

  *y_out = output;
  return outputlen;
}
}
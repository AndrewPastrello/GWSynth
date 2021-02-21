#include "parareal.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <lal/AVFactories.h>
#include <lal/XLALGSL.h>
#include <math.h>

int parareal(EOBParams *params, const REAL8 y_init[nn], const REAL8 t_init,
             REAL8 t_end, const REAL8 dt, const REAL8 pr_tol, int n_slices,
             REAL8Array **y_out) {
  // Phase cancels with tol 1e-7
  const REAL8 coarse_tol = 1e-5;
  const REAL8 fine_tol = 1e-7;

  REAL8 coarse_dt = (t_end - t_init) / n_slices;

  int bufferlength = (int)((t_end - t_init) / dt) + 2;
  REAL8 *buf = (REAL8 *)malloc((nn + 1) * sizeof(REAL8) * bufferlength);
  REAL8 *t_buf = buf;
  REAL8 *y_buf = buf + bufferlength;

  REAL8 *coarse_buf = (REAL8 *)malloc((nn + 1) * sizeof(REAL8) * bufferlength);
  REAL8 *coarse_t_buf = coarse_buf;
  REAL8 *coarse_y_buf = coarse_buf + bufferlength;

  // Fine solve
  REAL8 *y_initf = (REAL8 *)malloc(sizeof(REAL8) * n_slices * nn);
  REAL8 *y_initc = (REAL8 *)malloc(sizeof(REAL8) * n_slices * nn);
  REAL8 *y_correct = (REAL8 *)malloc(sizeof(REAL8) * n_slices * nn);
  REAL8 *y_fine0 = (REAL8 *)malloc(sizeof(REAL8) * n_slices * nn);
  REAL8 *y_coarse0 = (REAL8 *)malloc(sizeof(REAL8) * n_slices * nn);
  REAL8 *y_coarse1 = (REAL8 *)malloc(sizeof(REAL8) * n_slices * nn);

  REAL8 *t_initf = (REAL8 *)malloc(sizeof(REAL8) * n_slices);
  REAL8 *t_endf = (REAL8 *)malloc(sizeof(REAL8) * n_slices);

  EOBParams coarse_params;
  memcpy(&coarse_params, params, sizeof(*params));
  EOBParams *fine_params = (EOBParams *)malloc(sizeof(EOBParams) * n_slices);
  int *fine_steps = (int *)malloc(sizeof(int) * n_slices);

initialize:
  for (int i = 0; i < nn; i++) y_initc[i] = y_init[i];
  for (int k = 0; k < nn; k++) y_initf[k] = y_init[k];
  for (int i = 0; i < n_slices; i++) {
    t_initf[i] = t_init + i * coarse_dt;
    t_endf[i] = t_initf[i] + coarse_dt;
    memcpy(&fine_params[i], params, sizeof(EOBParams));
  }

  // Initial coarse solve
  for (int i = 0; i < n_slices; i++) {
    int steps;
    int stop_on_time = i == n_slices - 1 ? false : true;
    int stopped_on_condition =
        rk42(&coarse_params, &y_initc[i * nn], t_initf[i], t_endf[i], dt,
             coarse_tol, stop_on_time, coarse_t_buf, coarse_y_buf, &steps);
    if (stopped_on_condition && i < n_slices - 1) {
      t_end = coarse_t_buf[steps];
      coarse_dt = (t_end - t_init) / n_slices;
      memcpy(&coarse_params, params, sizeof(*params));
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
  memcpy(&coarse_params, params, sizeof(*params));

  // Parareal converges in max n_slices iterations
  for (int iteration = 0; iteration < n_slices; iteration++) {
    // Fine solve
#pragma omp parallel for
    for (int i = iteration; i < n_slices; i++) {
      int buf_offset = i * (bufferlength / n_slices);
      int stop_on_time = i == n_slices - 1 ? false : true;
      rk42(&fine_params[i], &y_initf[i * nn], t_initf[i], t_endf[i], dt,
           fine_tol, stop_on_time, &t_buf[buf_offset], &y_buf[buf_offset * nn],
           &fine_steps[i]);
      int ind = (buf_offset + fine_steps[i]) * nn;
      for (int k = 0; k < nn; k++) {
        y_fine0[i * nn + k] = y_buf[ind + k];
      }
    }
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
      REAL8 *y_new = y_correct + (i - 1) * nn;
      REAL8 *y_old = y_initc + i * nn;
      for (int k = 0; k < nn; k++) {
        REAL8 rerr = fabs((y_old[k] - y_new[k]) / y_new[k]);
        max_rerr = fmax(max_rerr, rerr);
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
      memcpy(&fine_params[i], params, sizeof(EOBParams));
    }
    REAL8 *tmp = y_coarse0;
    y_coarse0 = y_coarse1;
    y_coarse1 = tmp;
    memcpy(&coarse_params, params, sizeof(*params));
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

  REAL8 *y_flat = (REAL8 *)malloc((nn + 1) * sizeof(REAL8) * integrationlength);
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

  if (buf) free(buf);
  if (coarse_buf) free(coarse_buf);
  if (y_flat) free(y_flat);
  if (fine_steps) free(fine_steps);
  if (t_initf) free(t_initf);
  if (t_endf) free(t_endf);
  if (y_initf) free(y_initf);
  if (y_initc) free(y_initc);
  if (y_correct) free(y_correct);
  if (y_fine0) free(y_fine0);
  if (y_coarse0) free(y_coarse0);
  if (y_coarse1) free(y_coarse1);
  if (fine_params) free(fine_params);

  *y_out = output;
  return outputlen;
}

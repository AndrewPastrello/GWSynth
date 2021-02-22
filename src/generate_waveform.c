

#include <math.h>

#include "LALSimIMREOBNRv2.h"
#include "gwsynth.h"

#ifdef __cplusplus
extern "C" {
#endif
int generate_waveform(float* buf, const int buf_size, const int n_channels,
                      const int sample_rate, const float mass1,
                      const float mass2, const float inclination,
                      const float freq_min, const int n_slices) {
  REAL8TimeSeries* hplus = NULL;
  REAL8TimeSeries* hcross = NULL;
  const REAL8 phiRef = 0;
  const REAL8 deltaT = 1.0 / sample_rate;
  const REAL8 distance = 1e6;

  XLALSimIMREOBNRv2Generator(&hplus, &hcross, NULL, phiRef, deltaT,
                             mass1 * LAL_MSUN_SI, mass2 * LAL_MSUN_SI, freq_min,
                             distance, inclination, 1, n_slices);

  int length = hplus->data->length;
  double max = 0;
  for (int i = 0; i < length; i++) {
    double x = fabs(hplus->data->data[i]);
    max = x > max ? x : max;
  }
  if (n_channels == 2) {
    for (int i = 0; i < length; i++) {
      double y = fabs(hcross->data->data[i]);
      max = y > max ? y : max;
    }
  }
  for (int i = 0; i < length; i++) {
    buf[i] = (float)(hplus->data->data[i] / max);
  }
  if (n_channels == 2) {
    for (int i = 0; i < length; i++) {
      buf[i + length] = (float)(hcross->data->data[i] / max);
    }
  }

  if (hplus) XLALDestroyREAL8TimeSeries(hplus);
  if (hcross) XLALDestroyREAL8TimeSeries(hcross);
  return length;
}
#ifdef __cplusplus
}
#endif

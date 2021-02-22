#pragma once

#ifdef __cplusplus
extern "C" {
#endif
// Generates a binary black hole inspiral-merger-ringdown gravitational waveform
// using the EOBNRv2HM approximant
int generate_waveform(float* buf, const int buf_size, const int n_channels,
                      const int sample_rate, const float mass1,
                      const float mass2, const float inclination,
                      const float freq_min, const int n_slices);
#ifdef __cplusplus
}
#endif
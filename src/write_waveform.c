
#include <math.h>
#include <stdlib.h>

#include "gwsynth.h"
#include "tinywav/tinywav.c"

int write_wav(const char* filename, float* buf, const int length,
              const int n_channels, const int sample_rate) {
  // 1 to output hplus in mono, 2 to output hplus and hcross in stereo
  const int NUM_CHANNELS = 1;

  TinyWav tw;
  tinywav_open_write(
      &tw, NUM_CHANNELS, sample_rate,
      TW_FLOAT32,  // the output samples will be 32-bit floats. TW_INT16 is also
                   // supported
      TW_INLINE,   // the samples will be presented inlined in a single buffer.
                   // Other options include TW_INTERLEAVED and TW_SPLIT
      filename     // the output path
  );

  tinywav_write_f(&tw, buf, length);
  tinywav_close_write(&tw);

  return 0;
}

int main(int argc, char* argv[]) {
  const char* usage =
      "Generate a binary black hole inspiral-merger-ringdown "
      "gravitational waveform using the EOBNRv2 approximant from the "
      "lalsimulation library\n\n"
      "The following options can be given (will assume a default value if "
      "omitted):\n"
      "--m1 M1                    Mass of the 1st object in solar masses "
      "(default 10)\n"
      "--m2 M2                    Mass of the 2nd object in solar masses "
      "(default 1.5)\n"
      "--inclination IOTA         Angle in radians between line of sight (N) "
      "and \n"
      "                           orbital angular momentum (L) at the "
      "reference\n"
      "                           (default 0, face on)\n"
      "--filename FNAME           Output to file FNAME (default waveform.wav)\n"

      "--sample-rate SRATE        Sampling rate in Hz (default 96000)\n"
      "--n-channels NCHANNEL      The number of audio channels 1=mono, "
      "2=stereo (default 1)"
      "--n-slices NSLICE          The number of time slices (threads) to be "
      "used\n"
      "                           in the parareal algorithm (default "
      "32)\n";

  const char* default_filename = "waveform.wav";
  const char* filename = default_filename;

  float m1 = 10;
  float m2 = 1.5;
  float inclination = 0;
  float sample_rate = 96000;
  int n_channels = 1;
  int n_slices = 32;

  for (int i = 1; i < argc; ++i) {
    if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
      printf("%s", usage);
      exit(0);
    } else if (strcmp(argv[i], "--m1") == 0) {
      m1 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--m2") == 0) {
      m2 = atof(argv[++i]);
    } else if (strcmp(argv[i], "--inclination") == 0) {
      inclination = atof(argv[++i]);
    } else if (strcmp(argv[i], "--filename") == 0) {
      filename = argv[i + 1];
    } else if (strcmp(argv[i], "--sample-rate") == 0) {
      sample_rate = atof(argv[++i]);
    } else if (strcmp(argv[i], "--sample-rate") == 0) {
      n_channels = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--n-slices") == 0) {
      n_slices = atoi(argv[++i]);
    }
  }

  if (m2 > m1) {
    float tmp = m1;
    m1 = m2;
    m2 = tmp;
  }

  m1 = fmin(m1, 200);
  m1 = fmax(m1, 10);
  m2 = fmin(m2, 150);
  m2 = fmax(m2, 1.5);
  m2 = fmin(m2, 0.6 * m1);
  m2 = fmax(m2, 0.15 * m1);

  float freq_min = 10;
  freq_min = -1.0 / 3.0 * (m1 + m2) + 45;
  freq_min = fmin(freq_min, 45);
  freq_min = fmax(freq_min, 10);

  const int buf_length = n_channels * 1000000;
  float* buf = (float*)malloc(sizeof(float) * buf_length);

  int length = generate_waveform(buf, buf_length, n_channels, sample_rate, m1,
                                 m2, inclination, freq_min, n_slices);

  write_wav(filename, buf, length, n_channels, sample_rate);
  printf("Wrote %s\n", filename);
  free(buf);
}

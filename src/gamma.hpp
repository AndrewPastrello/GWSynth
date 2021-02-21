#pragma once

#include <cmath>
#include <complex>
#include <tuple>

#include "dpc_common.hpp"

#define LogRootTwoPi_ 0.9189385332046727418

constexpr double factorial(int n) {
  double r = 1;
  while (1 < n) r *= n--;
  return r;
}

auto complex_exp(std::complex<REAL8> exponent) {
  auto rx = std::real(exponent);
  auto iy = std::imag(exponent);

  return sycl::exp(rx) * (sycl::cos(iy) + I * sycl::sin(iy));
}

auto complex_log(const double zr, const double zi) {
  /* CHECK_POINTER(lnr) */
  /* CHECK_POINTER(theta) */

  const double ax = sycl::abs(zr);
  const double ay = sycl::abs(zi);
  const double min = sycl::min(ax, ay);
  const double max = sycl::max(ax, ay);
  auto lnr = sycl::log(max) + 0.5 * sycl::log(1.0 + (min / max) * (min / max));
  auto theta = sycl::atan2(zi, zr);
  return std::tuple{lnr, theta};
}

/* coefficients for gamma=7, kmax=8  Lanczos method */
constexpr double lanczos_7_c[9] = {0.99999999999980993227684700473478,
                                   676.520368121885098567009190444019,
                                   -1259.13921672240287047156078755283,
                                   771.3234287776530788486528258894,
                                   -176.61502916214059906584551354,
                                   12.507343278686904814458936853,
                                   -0.13857109526572011689554707,
                                   9.984369578019570859563e-6,
                                   1.50563273514931155834e-7};

/* complex version of Lanczos method; this is not safe for export
 * since it becomes bad in the left half-plane
 */
auto lngamma_lanczos_complex(double zr, double zi) {
  int k;

  double Ag_r, Ag_i;
  double yi_tmp_val, yi_tmp_err;

  zr -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag_r = lanczos_7_c[0];
  Ag_i = 0.0;
  for (k = 1; k <= 8; k++) {
    double R = zr + k;
    double I = zi;
    double a = lanczos_7_c[k] / (R * R + I * I);
    Ag_r += a * R;
    Ag_i -= a * I;
  }

  auto [log1_r, log1_i] = complex_log(zr + 7.5, zi);
  auto [logAg_r, logAg_i] = complex_log(Ag_r, Ag_i);

  /* (z+0.5)*log(z+7.5) - (z+7.5) + LogRootTwoPi_ + log(Ag(z)) */
  auto yr =
      (zr + 0.5) * log1_r - zi * log1_i - (zr + 7.5) + LogRootTwoPi_ + logAg_r;
  auto yi = zi * log1_r + (zr + 0.5) * log1_i - zi + logAg_i;

  return std::tuple{yr, yi};
}
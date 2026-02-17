#include "jet_tools/include/z_jet.h"
#include "jet_tools/include/math_helpers.h"

#include <cmath>
#include <limits>

namespace jet_tools {

std::optional<double> compute_z_inv(const P4& p_lab, const P4& q_lab,
                                    const P4& pjet_lab) {
  const double denom = minkowski_dot(p_lab, q_lab);
  if (!std::isfinite(denom) || std::abs(denom) <= 0.0)
    return std::nullopt;
  const double z_inv = minkowski_dot(p_lab, pjet_lab) / denom;
  if (!std::isfinite(z_inv))
    return std::nullopt;
  return z_inv;
}

std::optional<double> compute_z_breit(double E_B, double pz_B, double Q,
                                      double proton_sign) {
  if (!std::isfinite(Q) || Q <= 0.0)
    return std::nullopt;
  const double z_breit = (E_B + proton_sign * pz_B) / Q;
  if (!std::isfinite(z_breit))
    return std::nullopt;
  return z_breit;
}

}  // namespace jet_tools

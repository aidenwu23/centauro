#pragma once

#include <optional>

#include <Math/Vector4D.h>

namespace jet_tools {

using P4 = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;

std::optional<double> compute_z_inv(const P4& p_lab, const P4& q_lab,
                                    const P4& pjet_lab);

std::optional<double> compute_z_breit(double E_B, double pz_B, double Q,
                                      double proton_sign);

}  // namespace jet_tools

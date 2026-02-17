#pragma once

#include <limits>

#include <Math/Boost.h>
#include <Math/Rotation3D.h>
#include <Math/Vector4D.h>

#include "jet_tools/include/types.h"

namespace jet_tools {

enum class BreitXSource {
  Derived,
  Input,
};

bool transform_breit(const P4 &p_i, const P4 &e_i, const P4 &e_f,
                     ROOT::Math::Boost &out_boost,
                     ROOT::Math::Rotation3D &out_rot, BreitXSource x_source,
                     double x_input = std::numeric_limits<double>::quiet_NaN(),
                     double *Q_out = nullptr);

P4 breit_to_lab(const P4 &p_breit, const ROOT::Math::Boost &boost,
                const ROOT::Math::Rotation3D &rot);

P4 lab_to_breit(const P4 &p_lab, const ROOT::Math::Boost &boost,
                const ROOT::Math::Rotation3D &rot);

SimpleJet breit_to_lab(const SimpleJet &jet_breit,
                       const ROOT::Math::Boost &boost,
                       const ROOT::Math::Rotation3D &rot);

SimpleJet lab_to_breit(const SimpleJet &jet_lab,
                       const ROOT::Math::Boost &boost,
                       const ROOT::Math::Rotation3D &rot);

} // namespace jet_tools

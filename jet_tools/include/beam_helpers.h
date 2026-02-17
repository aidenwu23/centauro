#pragma once

#include <array>

#include <Math/Vector3D.h>

#include "jet_tools/include/types.h"

namespace edm4hep {
class MCParticleCollection;
}

namespace jet_tools {

P4 round_beam(const ROOT::Math::XYZVector &p_in, double mass,
              const std::array<double, 3> &nominal_pz,
              double crossing_angle);

bool find_beam(const edm4hep::MCParticleCollection &mc, int pdg, P4 &out);

}  // namespace jet_tools

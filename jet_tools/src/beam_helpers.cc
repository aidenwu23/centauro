#include "jet_tools/include/beam_helpers.h"

#include <cmath>

#include <edm4hep/MCParticleCollection.h>

namespace jet_tools {

P4 round_beam(const ROOT::Math::XYZVector &p_in, double mass,
              const std::array<double, 3> &nominal_pz,
              double crossing_angle) {
  double best_pz = nominal_pz.empty() ? 0.0 : nominal_pz.front();
  if (std::isfinite(p_in.Z()) && !nominal_pz.empty()) {
    double best_delta = std::abs(p_in.Z() - best_pz);
    for (const double pz : nominal_pz) {
      const double delta = std::abs(p_in.Z() - pz);
      if (delta < best_delta) {
        best_delta = delta;
        best_pz = pz;
      }
    }
    const double denom = std::abs(best_pz);
    if (denom > 0.0 && (best_delta / denom) > 0.1) {
      best_pz = p_in.Z();
    }
  } else if (std::isfinite(p_in.Z())) {
    best_pz = p_in.Z();
  }

  const double px = best_pz * std::sin(crossing_angle);
  const double pz = best_pz * std::cos(crossing_angle);
  const double py = 0.0;
  const double energy = std::sqrt(std::max(0.0, mass * mass + px * px + py * py + pz * pz));
  P4 out;
  out.SetPxPyPzE(px, py, pz, energy);
  return out;
}

bool find_beam(const edm4hep::MCParticleCollection &mc, int pdg, P4 &out) {
  constexpr double kElectronMass = 0.000511;
  constexpr double kProtonMass = 0.938272;

  for (const auto &p : mc) {
    if (p.getGeneratorStatus() == 4 && p.getPDG() == pdg) {
      const auto m = p.getMomentum();
      double e = p.getEnergy();
      if (!std::isfinite(e) || e <= 0.0) {
        const double mass = (std::abs(pdg) == 11) ? kElectronMass : kProtonMass;
        e = std::sqrt(m.x * m.x + m.y * m.y + m.z * m.z + mass * mass);
      }
      out.SetPxPyPzE(m.x, m.y, m.z, e);
      return true;
    }
  }

  for (const auto &p : mc) {
    if (p.getPDG() == pdg) {
      const auto m = p.getMomentum();
      double e = p.getEnergy();
      if (!std::isfinite(e) || e <= 0.0) {
        const double mass = (std::abs(pdg) == 11) ? kElectronMass : kProtonMass;
        e = std::sqrt(m.x * m.x + m.y * m.y + m.z * m.z + mass * mass);
      }
      out.SetPxPyPzE(m.x, m.y, m.z, e);
      return true;
    }
  }

  return false;
}

}  // namespace jet_tools

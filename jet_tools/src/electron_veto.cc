#include "jet_tools/include/electron_veto.h"
#include "jet_tools/include/math_helpers.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_set>
#include <vector>

#include <edm4eic/InclusiveKinematicsCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <podio/Frame.h>

namespace jet_tools {

static bool has_collection(const podio::Frame &frame, const std::string &name) {
  for (const auto &n : frame.getAvailableCollections()) {
    if (n == name)
      return true;
  }
  return false;
}

static std::vector<int> sorted_indices(const std::unordered_set<int> &indices) {
  std::vector<int> out(indices.begin(), indices.end());
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

int find_scattered_electron(const edm4eic::ReconstructedParticleCollection &particles) {
  auto pick_impl = [&](bool require_negative) {
    int index = -1;
    double best_pT = -1.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
      const auto &p = particles[i];
      if (std::abs(p.getPDG()) != 11) {
        continue;
      }
      if (require_negative && p.getCharge() >= 0) {
        continue;
      }
      const auto mom = p.getMomentum();
      const double pT = std::hypot(mom.x, mom.y);
      if (pT > best_pT) {
        best_pT = pT;
        index = static_cast<int>(i);
      }
    }
    return index;
  };

  int index = pick_impl(true);
  if (index < 0) {
    index = pick_impl(false);
  }
  return index;
}

TruthScatPickResult find_scattered_electron_truth(
    const edm4eic::ReconstructedParticleCollection &truth_particles) {
  TruthScatPickResult result;
  result.index = find_scattered_electron(truth_particles);
  return result;
}

TruthScatPickResult find_scattered_electron_truth_breit(
    const edm4eic::ReconstructedParticleCollection &lab_particles,
    const edm4eic::ReconstructedParticleCollection &breit_particles) {
  const int lab_index = find_scattered_electron(lab_particles);
  if (lab_index >= 0 && lab_particles.size() == breit_particles.size() &&
      static_cast<std::size_t>(lab_index) < breit_particles.size()) {
    TruthScatPickResult mapped;
    mapped.index = lab_index;
    return mapped;
  }
  TruthScatPickResult fallback;
  fallback.index = find_scattered_electron(breit_particles);
  return fallback;
}

ElectronVetoResult find_scattered_particles_reco(
    const podio::Frame &frame,
    const edm4eic::ReconstructedParticleCollection &reco_particles,
    const JetToolsCollections &collections, const ElectronVetoCuts &cuts) {
  ElectronVetoResult res;
  std::unordered_set<int> drop;

  if (!has_collection(frame, collections.electron_kinematics)) {
    return res;
  }
  const auto &kin = frame.get<edm4eic::InclusiveKinematicsCollection>(collections.electron_kinematics);
  if (kin.empty()) {
    return res;
  }

  const auto scat = kin[0].getScat();
  if (!scat.isAvailable()) {
    return res;
  }

  const int scat_index = scat.getObjectID().index;
  res.scat_index = scat_index;
  if (scat_index < 0 || static_cast<std::size_t>(scat_index) >= reco_particles.size()) {
    return res;
  }

  const auto mom_s = scat.getMomentum();
  const double pT_s = std::hypot(mom_s.x, mom_s.y);
  if (pT_s <= 0.0) {
    return res;
  }

  const double eta_s = std::asinh(mom_s.z / pT_s);
  const double phi_s = std::atan2(mom_s.y, mom_s.x);

  const auto &rp = reco_particles[static_cast<std::size_t>(scat_index)];
  const auto rm = rp.getMomentum();
  const double pT_r = std::hypot(rm.x, rm.y);
  const bool mom_ok = pT_r > 0.0 && std::abs(pT_r - pT_s) / pT_s < cuts.mom_rel_tol;
  const bool charge_ok = (scat.getCharge() == 0) || (rp.getCharge() == scat.getCharge());
  if (!mom_ok || !charge_ok) {
    return res;
  }

  drop.insert(scat_index);
  if (cuts.scat_cone_dr > 0.0) {
    for (std::size_t i = 0; i < reco_particles.size(); ++i) {
      const auto &cand = reco_particles[i];
      const auto mom = cand.getMomentum();
      const double pT = std::hypot(mom.x, mom.y);
      if (pT <= 0.0) {
        continue;
      }

      const double eta = std::asinh(mom.z / pT);
      const double phi = std::atan2(mom.y, mom.x);
      const double dphi = delta_phi(phi, phi_s);
      const double deta = eta - eta_s;
      const double dR = std::sqrt(dphi * dphi + deta * deta);
      if (dR > cuts.scat_cone_dr) {
        continue;
      }
      drop.insert(static_cast<int>(i));
    }
  }

  res.veto_indices = sorted_indices(drop);
  return res;
}

ElectronVetoResult find_scattered_particles_reco_breit(
    const podio::Frame &frame,
    const edm4eic::ReconstructedParticleCollection &reco_breit,
    const JetToolsCollections &collections, const ElectronVetoCuts &cuts) {
  ElectronVetoResult res;
  if (!has_collection(frame, collections.reco_lab_collection))
    return res;

  const auto &reco_lab = frame.get<edm4eic::ReconstructedParticleCollection>(collections.reco_lab_collection);
  const ElectronVetoResult veto_lab =
      find_scattered_particles_reco(frame, reco_lab, collections, cuts);

  // Mirror lab veto when lab and breit collections align; otherwise use direct breit veto.
  if (reco_lab.size() != reco_breit.size()) {
    return find_scattered_particles_reco(frame, reco_breit, collections, cuts);
  }
  if (veto_lab.veto_indices.empty())
    return res;

  std::unordered_set<int> drop;
  for (int index : veto_lab.veto_indices) {
    if (index < 0 || static_cast<std::size_t>(index) >= reco_breit.size())
      continue;
    drop.insert(index);
  }

  res.veto_indices = sorted_indices(drop);
  if (!res.veto_indices.empty()) {
    res.scat_index = veto_lab.scat_index;
  }
  return res;
}

int find_scattered_electron_index(
    const podio::Frame &frame,
    const edm4eic::ReconstructedParticleCollection &particles, bool use_reco,
    const JetToolsCollections &collections, const ElectronVetoCuts &cuts) {
  if (!use_reco)
    return find_scattered_electron(particles);
  const auto veto = find_scattered_particles_reco(frame, particles, collections, cuts);
  if (veto.scat_index >= 0)
    return veto.scat_index;
  if (!veto.veto_indices.empty())
    return veto.veto_indices.front();
  return -1;
}

} // namespace jet_tools

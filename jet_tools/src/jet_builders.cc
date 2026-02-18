#include "jet_tools/include/jet_builders.h"

#include <cmath>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

#include <edm4eic/ReconstructedParticleCollection.h>
#include <fastjet/ClusterSequence.hh>

namespace {

jet_tools::SimpleJet fill_one_simplejet(
    const fastjet::PseudoJet &jet,
    const edm4eic::ReconstructedParticleCollection *particles) {
  jet_tools::SimpleJet simple_jet{
      jet.perp(),
      jet.px(),
      jet.py(),
      jet.pz(),
      jet.E(),
      jet.eta(),
      jet.phi_std(),
      {},
      {},
      0,
      0,
      std::numeric_limits<double>::quiet_NaN()};

  const auto constituents = jet.constituents();
  simple_jet.n_reco_constituents = constituents.size();
  simple_jet.n_constituents = simple_jet.n_reco_constituents;

  simple_jet.constituent_indices.reserve(constituents.size());
  simple_jet.constituent_E.reserve(constituents.size());

  for (const auto &constituent : constituents) {
    const int user_index = constituent.user_index();
    double energy = constituent.E();
    if (particles && user_index >= 0 && static_cast<std::size_t>(user_index) < particles->size()) {
      const auto &particle = particles->at(static_cast<std::size_t>(user_index));
      energy = particle.getEnergy();
    }
    if (!std::isfinite(energy) || energy <= 0.0) {
      continue;
    }

    simple_jet.constituent_indices.push_back(user_index);
    simple_jet.constituent_E.push_back(energy);
  }
  return simple_jet;
}

} // namespace

namespace jet_tools {

void build_skip_indices(
    const edm4eic::ReconstructedParticleCollection &particles,
    const ConstituentGateCuts &cuts,
    const std::unordered_set<int> *base_skip_indices, // From the scattered-electron list
    std::unordered_set<int> &skip_out) {
  skip_out.clear();

  // Set the scattered electron indices as the bare minimum.
  if (base_skip_indices) {
    skip_out = *base_skip_indices;
  }

  // Loop through particles and add on to the scattered electron indices.
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const auto &particle = particles[i];

    // Neutrino cut.
    const int ap = std::abs(particle.getPDG());
    if (ap == 12 || ap == 14 || ap == 16) {
      skip_out.insert(static_cast<int>(i));
      continue;
    }

    // Particle pT cut.
    const auto mom = particle.getMomentum();
    const double pT = std::hypot(mom.x, mom.y);
    if (!std::isfinite(pT) || pT <= 0.0 || pT < cuts.min_cst_pT ||
        pT > cuts.max_cst_pT) {
      skip_out.insert(static_cast<int>(i));
      continue;
    }

    // Particle eta cut.
    const double eta = std::asinh(mom.z / pT);
    if (!std::isfinite(eta) || std::abs(eta) > cuts.max_cst_eta) {
      skip_out.insert(static_cast<int>(i));
      continue;
    }

    // If valid energy, add to list of particles to skip before clustering.
    const double energy = particle.getEnergy();
    if (!std::isfinite(energy) || energy <= 0.0) {
      skip_out.insert(static_cast<int>(i));
    }
  }
}

std::vector<fastjet::PseudoJet> build_input_pseudojets(
    const edm4eic::ReconstructedParticleCollection &particles,
    const std::unordered_set<int> *skip_indices) {
  std::vector<fastjet::PseudoJet> inputs;
  inputs.reserve(particles.size());

  for (std::size_t index = 0; index < particles.size(); ++index) {
    // If part of a valid skip list, skip.
    if (skip_indices &&
        skip_indices->count(static_cast<int>(index)) > 0) {
      continue;
    }

    // If not a particle to skip, grab particle momenta.
    const auto &particle = particles[index];
    const auto mom = particle.getMomentum();
    const double px = mom.x;
    const double py = mom.y;
    const double pz = mom.z;
    const double energy = particle.getEnergy();
    if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(energy)) {
      continue;
    }

    // Build a FastJet input and keep the source particle index.
    fastjet::PseudoJet pseudojet(px, py, pz, energy);
    pseudojet.set_user_index(static_cast<int>(index));
    inputs.push_back(std::move(pseudojet));
  }

  return inputs;
}

std::vector<SimpleJet> build_simplejet_from_pseudojets(
    const std::vector<fastjet::PseudoJet> &inputs,
    const fastjet::JetDefinition &jet_def,
    const edm4eic::ReconstructedParticleCollection *particles, double max_eta,
    double min_pT) {
  if (inputs.empty()) {
    return {};
  }

  fastjet::ClusterSequence sequence(inputs, jet_def);

  // Use inclusive jets to apply the momenta cut.
  const auto jets = sequence.inclusive_jets(min_pT);
  std::vector<SimpleJet> out;
  out.reserve(jets.size());

  // Apply the jet eta cut.
  for (const auto &jet : jets) {
    if (std::abs(jet.eta()) > max_eta) {
      continue;
    }
    out.push_back(fill_one_simplejet(jet, particles));
  }

  return out;
}

} // namespace jet_tools

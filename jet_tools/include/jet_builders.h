#pragma once

#include <fastjet/PseudoJet.hh>

#include <limits>
#include <unordered_set>
#include <vector>

#include "jet_tools/include/types.h"

namespace edm4eic {
class ReconstructedParticleCollection;
}  // namespace edm4eic

namespace fastjet {
class JetDefinition;
}  // namespace fastjet

namespace jet_tools {

struct ConstituentGateCuts {
  double min_cst_pT = 0.0;
  double max_cst_pT = std::numeric_limits<double>::infinity();
  double max_cst_eta = std::numeric_limits<double>::infinity();
};

void build_skip_indices(
    const edm4eic::ReconstructedParticleCollection &particles,
    const ConstituentGateCuts &cuts,
    const std::unordered_set<int> *base_skip_indices,
    std::unordered_set<int> &skip_out);

std::vector<fastjet::PseudoJet> build_input_pseudojets(
    const edm4eic::ReconstructedParticleCollection &particles,
    const std::unordered_set<int> *skip_indices = nullptr);

std::vector<SimpleJet> build_simplejet_from_pseudojets(
    const std::vector<fastjet::PseudoJet> &inputs,
    const fastjet::JetDefinition &jet_def,
    const edm4eic::ReconstructedParticleCollection *particles, double max_eta,
    double min_pT);

} // namespace jet_tools

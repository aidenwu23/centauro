#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <podio/Frame.h>
#include <edm4eic/ReconstructedParticleCollection.h>

namespace jet_tools {

struct JetToolsCollections {
  std::string electron_kinematics = "InclusiveKinematicsElectron";
  std::string reco_lab_collection = "ReconstructedParticles";
  std::string scat_electron_truth = "ScatteredElectronsTruth";
};

struct ElectronVetoCuts {
  double max_dpT_rel = 0.3;
  double max_dr = 0.1;
  double mom_rel_tol = 1e-3;
  double scat_cone_dr = 0.35;
};

struct ElectronVetoResult {
  std::vector<int> veto_indices; // sorted, unique
  int scat_index = -1; // selected scattered-electron anchor index
};

struct TruthScatPickResult {
  int index = -1;
};

// Pick highest-pT scattered-electron candidate by PDG only.
int find_scattered_electron(
    const edm4eic::ReconstructedParticleCollection &particles);

// Pick the index of the highest-pT truth electron.
TruthScatPickResult find_scattered_electron_truth(
    const edm4eic::ReconstructedParticleCollection &truth_particles);

// Pick truth scattered electron in Breit collection by lab->Breit index mapping when aligned.
TruthScatPickResult find_scattered_electron_truth_breit(
    const edm4eic::ReconstructedParticleCollection &lab_particles,
    const edm4eic::ReconstructedParticleCollection &breit_particles);

// Pick scattered-electron index for reco or truth collections.
int find_scattered_electron_index(
    const podio::Frame &frame,
    const edm4eic::ReconstructedParticleCollection &particles, bool use_reco,
    const JetToolsCollections &collections = {},
    const ElectronVetoCuts &cuts = {});

// Find reconstructed-particle indices belonging to the scattered electron in the lab frame.
ElectronVetoResult find_scattered_particles_reco(
    const podio::Frame &frame,
    const edm4eic::ReconstructedParticleCollection &reco_particles,
    const JetToolsCollections &collections = {},
    const ElectronVetoCuts &cuts = {});

// Find reconstructed-particle indices belonging to the scattered electron in the Breit frame.
ElectronVetoResult find_scattered_particles_reco_breit(
    const podio::Frame &frame,
    const edm4eic::ReconstructedParticleCollection &reco_breit,
    const JetToolsCollections &collections = {},
    const ElectronVetoCuts &cuts = {});

} // namespace jet_tools

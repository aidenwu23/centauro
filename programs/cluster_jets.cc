/*
Example:
Runs clustering with jet radius = 0.6, eta < 4, and jet cuts applied in whatever frame the input
collection is. 

./run.sh programs/cluster_jets.cc \
  -i data/events/filtered/merged_0-999.root \
  -o data/jets/raw/jets_0-999_lab_cuts.root --threshold-frame lab

./run.sh programs/cluster_jets.cc \
  -i data/events/filtered/merged_0-999.root \
  -o data/jets/raw/jets_0-999_native_cuts.root --threshold-frame native \

Runs clustering with jet radius = 0.6, Breit eta < 4, breit pT > 5, and jet cuts applied in whatever 
frame the input collection is. 

./run.sh programs/cluster_jets.cc \
  -i data/events/filtered/merged_0-999.root \
  -o data/jets/raw/jets_0-999_lab_cuts.root \
  --threshold-frame lab \
  --min-jet-pT 0 --max-eta 6 --min-cst-pT 0.2 --max-cst-eta 4 
*/
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <array>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include <Math/Boost.h>
#include <Math/Rotation3D.h>
#include <Math/Vector3D.h>
#include <TFile.h>
#include <TTree.h>
#include <fastjet/JetDefinition.hh>
#include <fastjet/contrib/Centauro.hh>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/InclusiveKinematicsCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4hep/MCParticleCollection.h>

#include "jet_tools/include/beam_helpers.h"
#include "jet_tools/include/electron_veto.h"
#include "jet_tools/include/event_progress.h"
#include "jet_tools/include/jet_builders.h"
#include "jet_tools/include/root_io.h"
#include "jet_tools/include/transform_breit.h"
#include "jet_tools/include/types.h"
#include "jet_tools/include/z_jet.h"

namespace {

using P4 = jet_tools::P4;
using Vec3 = ROOT::Math::XYZVector;

constexpr double kElectronMass = 0.000511;
constexpr double kProtonMass = 0.938272;
constexpr double kCrossingAngle = -0.025;
const std::array<double, 3> kElectronNominalPz{-5.0, -10.0, -18.0};
const std::array<double, 3> kProtonNominalPz{41.0, 100.0, 275.0};

enum class ThresholdFrame {
  Native,
  Lab,
  Breit,
};

enum class JetFrame {
  Lab,
  Breit,
};

struct EventKinematics {
  P4 p_in;
  P4 q_lab;
  double Q2 = std::numeric_limits<double>::quiet_NaN();
  ROOT::Math::Boost boost;
  ROOT::Math::Rotation3D rot;
  bool have_transform = false;
  bool valid = false;
};

struct Args {
  std::string input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  double R = 0.6;
  double max_eta = 4.0;
  double min_jet_pT = 5.0;
  double min_cst_pT = 0.0;
  double max_cst_pT = std::numeric_limits<double>::infinity();
  double max_cst_eta = std::numeric_limits<double>::infinity();
  bool drop_scat_electron = true;
  ThresholdFrame threshold_frame = ThresholdFrame::Native;
};

//=========================================================================================
// CLI helpers
//=========================================================================================
void usage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " -i INPUT.root -o OUTPUT.root [--max-events N] [-R VAL]\n"
            << "       [--max-eta VAL] [--min-jet-pT VAL] [--min-cst-pT VAL]\n"
            << "       [--max-cst-pT VAL] [--max-cst-eta VAL] [--keep-scat-electron]\n"
            << "       [--threshold-frame native|lab|breit]\n";
}

bool parse_threshold_frame(std::string_view value, ThresholdFrame &out) {
  if (value == "native") {
    out = ThresholdFrame::Native;
    return true;
  }
  if (value == "lab") {
    out = ThresholdFrame::Lab;
    return true;
  }
  if (value == "breit") {
    out = ThresholdFrame::Breit;
    return true;
  }
  return false;
}

Args parse_args(int argc, char *argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    std::string_view v(argv[i]);
    if ((v == "-i" || v == "--input") && i + 1 < argc) {
      args.input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--max-events" && i + 1 < argc) {
      const long long n = std::stoll(argv[++i]);
      if (n < 0) {
        std::cerr << "[cluster_jets] max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else if (v == "-R" && i + 1 < argc) {
      args.R = std::atof(argv[++i]);
    } else if (v == "--max-eta" && i + 1 < argc) {
      args.max_eta = std::atof(argv[++i]);
    } else if (v == "--min-jet-pT" && i + 1 < argc) {
      args.min_jet_pT = std::atof(argv[++i]);
    } else if (v == "--min-cst-pT" && i + 1 < argc) {
      args.min_cst_pT = std::atof(argv[++i]);
    } else if (v == "--max-cst-pT" && i + 1 < argc) {
      args.max_cst_pT = std::atof(argv[++i]);
    } else if (v == "--max-cst-eta" && i + 1 < argc) {
      args.max_cst_eta = std::atof(argv[++i]);
    } else if (v == "--keep-scat-electron") {
      args.drop_scat_electron = false;
    } else if (v == "--threshold-frame" && i + 1 < argc) {
      const std::string_view mode(argv[++i]);
      if (!parse_threshold_frame(mode, args.threshold_frame)) {
        std::cerr << "[cluster_jets] threshold-frame must be native, lab, or breit\n";
        std::exit(1);
      }
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (!std::isfinite(args.R) || args.R <= 0.0) {
    std::cerr << "[cluster_jets] R must be finite and > 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.max_eta) || args.max_eta <= 0.0) {
    std::cerr << "[cluster_jets] max-eta must be finite and > 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.min_jet_pT) || args.min_jet_pT < 0.0) {
    std::cerr << "[cluster_jets] min-jet-pT must be finite and >= 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.min_cst_pT) || args.min_cst_pT < 0.0) {
    std::cerr << "[cluster_jets] min-cst-pT must be finite and >= 0\n";
    std::exit(1);
  }
  if (args.max_cst_pT <= 0.0) {
    std::cerr << "[cluster_jets] max-cst-pT must be > 0\n";
    std::exit(1);
  }
  if (args.max_cst_eta <= 0.0) {
    std::cerr << "[cluster_jets] max-cst-eta must be > 0\n";
    std::exit(1);
  }
  return args;
}

bool has_collection(const podio::Frame &frame, const std::string &name) {
  for (const auto &n : frame.getAvailableCollections()) {
    if (n == name) {
      return true;
    }
  }
  return false;
}

bool build_event_kinematics(const podio::Frame &frame, bool use_reco,
                            const jet_tools::JetToolsCollections &collections,
                            const jet_tools::ElectronVetoCuts &cuts,
                            EventKinematics &out) {
  out = EventKinematics{};
  if (!has_collection(frame, "MCParticles")) {
    return false;
  }
  const auto &mc = frame.get<edm4hep::MCParticleCollection>("MCParticles");
  if (mc.empty()) {
    return false;
  }

  P4 p_in_raw;
  P4 k_in_raw;
  if (!jet_tools::find_beam(mc, 2212, p_in_raw) ||
      !jet_tools::find_beam(mc, 11, k_in_raw)) {
    return false;
  }

  const Vec3 p_vec(p_in_raw.Px(), p_in_raw.Py(), p_in_raw.Pz());
  const Vec3 k_vec(k_in_raw.Px(), k_in_raw.Py(), k_in_raw.Pz());
  const P4 p_in = jet_tools::round_beam(p_vec, kProtonMass, kProtonNominalPz, kCrossingAngle);
  const P4 k_in = jet_tools::round_beam(k_vec, kElectronMass, kElectronNominalPz, 0.0);

  const std::string particles_name = use_reco ? collections.reco_lab_collection : "GeneratedParticles";
  if (!has_collection(frame, particles_name)) {
    return false;
  }
  const auto &particles = frame.get<edm4eic::ReconstructedParticleCollection>(particles_name);
  if (particles.empty()) {
    return false;
  }

  const int scat_idx = jet_tools::find_scattered_electron_index(
      frame, particles, use_reco, collections, cuts);
  if (scat_idx < 0 || static_cast<std::size_t>(scat_idx) >= particles.size()) {
    return false;
  }

  const auto &scat = particles[static_cast<std::size_t>(scat_idx)];
  const auto m = scat.getMomentum();
  double e_f = scat.getEnergy();
  if (!std::isfinite(e_f) || e_f <= 0.0) {
    e_f = std::sqrt(m.x * m.x + m.y * m.y + m.z * m.z +
                    kElectronMass * kElectronMass);
  }
  const P4 k_out(m.x, m.y, m.z, e_f);
  const P4 q_lab = k_in - k_out;
  const double Q2 = -q_lab.M2();
  if (!std::isfinite(Q2) || Q2 <= 0.0) {
    return false;
  }

  jet_tools::BreitXSource x_source = jet_tools::BreitXSource::Derived;
  double x_input = std::numeric_limits<double>::quiet_NaN();
  const std::string kin_name = use_reco ? collections.electron_kinematics : "InclusiveKinematicsTruth";
  if (has_collection(frame, kin_name)) {
    const auto &kin =
        frame.get<edm4eic::InclusiveKinematicsCollection>(kin_name);
    if (!kin.empty()) {
      const double x = kin[0].getX();
      if (std::isfinite(x) && x > 0.0) {
        x_source = jet_tools::BreitXSource::Input;
        x_input = x;
      }
    }
  }

  out.p_in = p_in;
  out.q_lab = q_lab;
  out.Q2 = Q2;
  out.have_transform = jet_tools::transform_breit(p_in, k_in, k_out, out.boost, out.rot, x_source, x_input);
  out.valid = true;
  return true;
}

std::unordered_set<int> to_set(const std::vector<int> &indices) {
  std::unordered_set<int> out;
  out.reserve(indices.size());
  for (int idx : indices) {
    out.insert(idx);
  }
  return out;
}

bool pass_gates_jet(const jet_tools::SimpleJet &jet, double min_jet_pT,
                    double max_eta) {
  if (!std::isfinite(jet.pT) || jet.pT < min_jet_pT) {
    return false;
  }
  if (!std::isfinite(jet.eta) || std::abs(jet.eta) > max_eta) {
    return false;
  }
  return true;
}

// Make sure that the pT and eta cuts either happen in the collection's native frame,
// entirely in the lab frame, or entirely in the breit frame.
bool pass_threshold_frame(const jet_tools::SimpleJet &jet, JetFrame jet_frame,
                          ThresholdFrame threshold_frame, double min_jet_pT,
                          double max_eta, bool have_transform,
                          const ROOT::Math::Boost &boost,
                          const ROOT::Math::Rotation3D &rot) {
  if (threshold_frame == ThresholdFrame::Native) {
    return pass_gates_jet(jet, min_jet_pT, max_eta);
  }

  // Lab frame cuts.
  if (threshold_frame == ThresholdFrame::Lab) {
    if (jet_frame == JetFrame::Lab) {
      return pass_gates_jet(jet, min_jet_pT, max_eta);
    }
    if (!have_transform) {
      return false;
    }
    // Jet is in Breit, so transform to lab for gate checks.
    const auto jet_lab = jet_tools::breit_to_lab(jet, boost, rot);
    return pass_gates_jet(jet_lab, min_jet_pT, max_eta);
  }

  // Breit frame cuts.
  if (jet_frame == JetFrame::Breit) {
    return pass_gates_jet(jet, min_jet_pT, max_eta);
  }
  if (!have_transform) {
    return false;
  }
  // Jet is in lab, so transform to Breit for gate checks.
  const auto jet_breit = jet_tools::lab_to_breit(jet, boost, rot);
  // Return true/false.
  return pass_gates_jet(jet_breit, min_jet_pT, max_eta);
}

std::vector<jet_tools::SimpleJet> filter_threshold_frame(
    const std::vector<jet_tools::SimpleJet> &jets, JetFrame jet_frame,
    ThresholdFrame threshold_frame, double min_jet_pT, double max_eta,
    bool have_transform, const ROOT::Math::Boost &boost,
    const ROOT::Math::Rotation3D &rot) {
  std::vector<jet_tools::SimpleJet> out;
  out.reserve(jets.size());
  for (const auto &jet : jets) {
    if (pass_threshold_frame(jet, jet_frame, threshold_frame, min_jet_pT,
                             max_eta, have_transform, boost, rot)) {
      out.push_back(jet);
    }
  }
  return out;
}

void attach_z_inv(std::vector<jet_tools::SimpleJet> &jets, JetFrame jet_frame,
                  const EventKinematics &kin,
                  const edm4eic::ReconstructedParticleCollection *lab_particles) {
  (void)jet_frame;
  for (auto &jet : jets) {
    jet.z_inv = std::numeric_limits<double>::quiet_NaN();
    if (!kin.valid || !lab_particles || jet.constituents.empty()) {
      continue;
    }

    double sum_e = 0.0;
    double sum_px = 0.0;
    double sum_py = 0.0;
    double sum_pz = 0.0;
    bool have_constituent = false;
    for (int index : jet.constituents) {
      if (index < 0) {
        continue;
      }
      const std::size_t ui = static_cast<std::size_t>(index);
      if (ui >= lab_particles->size()) {
        continue;
      }
      const auto &particle = (*lab_particles)[ui];
      const auto mom = particle.getMomentum();
      if (!std::isfinite(mom.x) || !std::isfinite(mom.y) || !std::isfinite(mom.z)) {
        continue;
      }
      const double e_in = particle.getEnergy();
      const double e_use = (std::isfinite(e_in) && e_in > 0.0)
                               ? e_in
                               : std::sqrt(mom.x * mom.x + mom.y * mom.y + mom.z * mom.z);
      sum_e += e_use;
      sum_px += mom.x;
      sum_py += mom.y;
      sum_pz += mom.z;
      have_constituent = true;
    }

    if (!have_constituent || !std::isfinite(sum_e) || sum_e <= 0.0) {
      continue;
    }
    const P4 p_jet_lab(sum_px, sum_py, sum_pz, sum_e);

    const auto z_inv = jet_tools::compute_z_inv(kin.p_in, kin.q_lab, p_jet_lab);
    if (z_inv.has_value() && std::isfinite(*z_inv)) {
      jet.z_inv = *z_inv;
    }
  }
}

struct JetTreeWriter {
  explicit JetTreeWriter(const char *name, const char *title) : tree(name, title) {
    tree.Branch("event_id", &event_id);
    tree.Branch("pT", &pT);
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("E", &E);
    tree.Branch("eta", &eta);
    tree.Branch("phi", &phi);
    tree.Branch("Q2", &Q2);
    tree.Branch("z_inv", &z_inv);
    tree.Branch("n_constituents", &n_constituents);
    tree.Branch("n_reco_constituents", &n_reco_constituents);
    tree.Branch("constituents", &constituents);
    tree.Branch("constituent_px", &constituent_px);
    tree.Branch("constituent_py", &constituent_py);
    tree.Branch("constituent_pz", &constituent_pz);
    tree.Branch("constituent_E", &constituent_E);
  }

  void fill(ULong64_t event, const jet_tools::SimpleJet &jet, double q2) {
    event_id = event;
    pT = jet.pT;
    px = jet.px;
    py = jet.py;
    pz = jet.pz;
    E = jet.E;
    eta = jet.eta;
    phi = jet.phi;
    Q2 = q2;
    z_inv = jet.z_inv;
    n_constituents = static_cast<int>(jet.n_constituents);
    n_reco_constituents = static_cast<int>(jet.n_reco_constituents);
    constituents = jet.constituents;
    constituent_px = jet.constituent_px;
    constituent_py = jet.constituent_py;
    constituent_pz = jet.constituent_pz;
    constituent_E = jet.constituent_E;
    tree.Fill();
  }

  void fill_many(ULong64_t event, const std::vector<jet_tools::SimpleJet> &jets,
                 double q2) {
    for (const auto &jet : jets) {
      fill(event, jet, q2);
    }
  }

  TTree tree;
  ULong64_t event_id = 0;
  double pT = 0.0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double E = 0.0;
  double eta = 0.0;
  double phi = 0.0;
  double Q2 = std::numeric_limits<double>::quiet_NaN();
  double z_inv = std::numeric_limits<double>::quiet_NaN();
  int n_constituents = 0;
  int n_reco_constituents = 0;
  std::vector<int> constituents;
  std::vector<double> constituent_px;
  std::vector<double> constituent_py;
  std::vector<double> constituent_pz;
  std::vector<double> constituent_E;
};

} // namespace

int main(int argc, char *argv[]) {
  const Args args = parse_args(argc, argv);
  jet_tools::ensure_parent(args.output);

  // Grab events.
  podio::ROOTReader reader;
  reader.openFile(args.input);
  const std::size_t entries_total = reader.getEntries("events");
  if (entries_total == 0) {
    std::cerr << "[cluster_jets] no events in input\n";
    return 1;
  }
  const std::size_t entries = std::min(entries_total, args.max_events);

  // Jet definitions.
  fastjet::JetDefinition jetdef_antikt(fastjet::antikt_algorithm, args.R, fastjet::E_scheme);
  fastjet::JetDefinition jetdef_centauro(new fastjet::contrib::CentauroPlugin(args.R));
  jetdef_centauro.set_recombination_scheme(fastjet::E_scheme);

  const jet_tools::JetToolsCollections collections;
  jet_tools::ElectronVetoCuts veto_cuts;

  // Disable native jet-level pT/eta gates when threshold_frame != native.
  // This defers jet thresholding to filter_threshold_frame so all jet collections
  // are evaluated in one consistent frame.
  const double native_max_eta = (args.threshold_frame == ThresholdFrame::Native)
                                    ? args.max_eta : std::numeric_limits<double>::infinity();
  const double native_min_jet_pT = (args.threshold_frame == ThresholdFrame::Native)
                                       ? args.min_jet_pT : 0.0;

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[cluster_jets] failed to open output " << args.output << "\n";
    return 1;
  }

  // Initialize SimpleJet --> ROOT tree writers.
  JetTreeWriter t_lab_antikt_truth("LabFrameAntiktTruth", "lab antikt truth jets");
  JetTreeWriter t_lab_cent_truth("LabFrameCentauroTruth", "lab centauro truth jets");
  JetTreeWriter t_breit_antikt_truth("BreitFrameAntiktTruth", "breit antikt truth jets");
  JetTreeWriter t_breit_cent_truth("BreitFrameCentauroTruth", "breit centauro truth jets");
  JetTreeWriter t_lab_antikt_reco("LabFrameAntiktReco", "lab antikt reco jets");
  JetTreeWriter t_lab_cent_reco("LabFrameCentauroReco", "lab centauro reco jets");
  JetTreeWriter t_breit_antikt_reco("BreitFrameAntiktReco", "breit antikt reco jets");
  JetTreeWriter t_breit_cent_reco("BreitFrameCentauroReco", "breit centauro reco jets");
  jet_tools::ProgressTicker progress;

  // Event loop.
  for (std::size_t i = 0; i < entries; ++i) {
    if (progress.should_report(i)) {
      const std::string message = std::string("[cluster_jets] event ") +
                                  std::to_string(i) + "/" +
                                  std::to_string(entries);
      progress.report(message);
    }

    auto data = reader.readEntry("events", i);
    if (!data) {
      continue;
    }
    podio::Frame frame(std::move(data));

    // Grab event id.
    const ULong64_t event_id = static_cast<ULong64_t>(i);

    EventKinematics truth_kin;
    EventKinematics reco_kin;
    build_event_kinematics(frame, false, collections, veto_cuts, truth_kin);
    build_event_kinematics(frame, true, collections, veto_cuts, reco_kin);

    //=========================================================================================
    // Grab inputs.
    //=========================================================================================
    // Check for collections.
    const bool have_truth_lab = has_collection(frame, "GeneratedParticles");
    const bool have_truth_breit = has_collection(frame, "GeneratedBreitFrameParticles");
    const bool have_reco_lab = has_collection(frame, "ReconstructedParticles");
    const bool have_reco_breit = has_collection(frame, "ReconstructedBreitFrameParticles");

    const edm4eic::ReconstructedParticleCollection *truth_lab_particles = nullptr;
    const edm4eic::ReconstructedParticleCollection *truth_breit_particles = nullptr;
    const edm4eic::ReconstructedParticleCollection *reco_lab_particles = nullptr;
    const edm4eic::ReconstructedParticleCollection *reco_breit_particles = nullptr;

    // If collections exist, grab the particles specifically for this event per collection.
    if (have_truth_lab) {
      truth_lab_particles = &frame.get<edm4eic::ReconstructedParticleCollection>("GeneratedParticles");
    }
    if (have_truth_breit) {
      truth_breit_particles = &frame.get<edm4eic::ReconstructedParticleCollection>("GeneratedBreitFrameParticles");
    }
    if (have_reco_lab) {
      reco_lab_particles = &frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
    }
    if (have_reco_breit) {
      reco_breit_particles = &frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedBreitFrameParticles");
    }

    //=========================================================================================
    // Build list of scattered electron indices to skip.
    //=========================================================================================
    std::unordered_set<int> skip_truth_lab;
    std::unordered_set<int> skip_truth_breit;
    std::unordered_set<int> skip_reco_lab;
    std::unordered_set<int> skip_reco_breit;
    std::unordered_set<int> skip_truth_lab_all;
    std::unordered_set<int> skip_truth_breit_all;
    std::unordered_set<int> skip_reco_lab_all;
    std::unordered_set<int> skip_reco_breit_all;

    if (args.drop_scat_electron) {
      // For GeneratedParticles
      if (truth_lab_particles) {
        const auto pick = jet_tools::find_scattered_electron_truth(*truth_lab_particles);
        if (pick.index >= 0 && static_cast<std::size_t>(pick.index) < truth_lab_particles->size()) {
          skip_truth_lab.insert(pick.index);
        }
      }

      // For GeneratedBreitFrameParticles
      if (truth_breit_particles) {
        jet_tools::TruthScatPickResult pick;
        // Use the lab particle indices assuming the Breit collection lines up.
        if (truth_lab_particles) {
          pick = jet_tools::find_scattered_electron_truth_breit(*truth_lab_particles, *truth_breit_particles);
        } else {
          pick = jet_tools::find_scattered_electron_truth(*truth_breit_particles);
        }
        if (pick.index >= 0 && static_cast<std::size_t>(pick.index) < truth_breit_particles->size()) {
          skip_truth_breit.insert(pick.index);
        }
      }

      // For ReconstructedParticles.
      if (reco_lab_particles) {
        const auto veto = jet_tools::find_scattered_particles_reco(frame, *reco_lab_particles, collections, 
            veto_cuts);
        skip_reco_lab = to_set(veto.veto_indices);
      }

      // For ReconstructedBreitFrameParticles.
      if (reco_breit_particles) {
        const auto veto = jet_tools::find_scattered_particles_reco_breit(frame, *reco_breit_particles, 
            collections, veto_cuts);
        skip_reco_breit = to_set(veto.veto_indices);
      }
    }

    const jet_tools::ConstituentGateCuts cst_cuts{
        args.min_cst_pT, args.max_cst_pT, args.max_cst_eta};

    if (args.threshold_frame == ThresholdFrame::Native) {
      // Each collection builds the skip indices from themselves.
      if (truth_lab_particles) {
        jet_tools::build_skip_indices(*truth_lab_particles, cst_cuts, &skip_truth_lab, skip_truth_lab_all);
      }
      if (truth_breit_particles) {
        jet_tools::build_skip_indices(*truth_breit_particles, cst_cuts, &skip_truth_breit, skip_truth_breit_all);
      }
      if (reco_lab_particles) {
        jet_tools::build_skip_indices(*reco_lab_particles, cst_cuts, &skip_reco_lab, skip_reco_lab_all);
      }
      if (reco_breit_particles) {
        jet_tools::build_skip_indices(*reco_breit_particles, cst_cuts, &skip_reco_breit, skip_reco_breit_all);
      }
    } else {
      // Build skip indices in the chosen threshold source frame, then reuse them in
      // the target frame collections; this relies on 1-to-1 lab/breit index mapping.

      // Example: for lab-frame cuts on breit collections, build skip indices from lab particles and
      // apply the same indices to breit particles.
      const edm4eic::ReconstructedParticleCollection *truth_source_particles = nullptr;
      const edm4eic::ReconstructedParticleCollection *reco_source_particles = nullptr;

      // Scattered-electron skip indices in the selected source frame.
      const std::unordered_set<int> *truth_electron_skip = nullptr;
      const std::unordered_set<int> *reco_electron_skip = nullptr;
      if (args.threshold_frame == ThresholdFrame::Lab) {
        // Use lab frame particles if lab frame cuts.
        truth_source_particles = truth_lab_particles;
        reco_source_particles = reco_lab_particles;

        // Use scattered electron indices as the base list, and then append additional indices later on.
        // (mostly particles that fail the configured pT, eta, and neutrino cuts).
        truth_electron_skip = &skip_truth_lab;
        reco_electron_skip = &skip_reco_lab;
      } else {
        // Breit particles for breit frame cuts.
        truth_source_particles = truth_breit_particles;
        reco_source_particles = reco_breit_particles;
        
        truth_electron_skip = &skip_truth_breit;
        reco_electron_skip = &skip_reco_breit;
      }

      if (truth_lab_particles && truth_breit_particles &&
          truth_lab_particles->size() != truth_breit_particles->size()) {
        std::cerr << "[cluster_jets] truth lab/breit size mismatch at event "
                  << event_id << "\n";
        return 1;
      }
      if (reco_lab_particles && reco_breit_particles &&
          reco_lab_particles->size() != reco_breit_particles->size()) {
        std::cerr << "[cluster_jets] reco lab/breit size mismatch at event "
                  << event_id << "\n";
        return 1;
      }

      // Append additional indices to the base electron list.
      if (truth_source_particles) {
        jet_tools::build_skip_indices(*truth_source_particles, cst_cuts, truth_electron_skip, skip_truth_lab_all);
        skip_truth_breit_all = skip_truth_lab_all;
      } else {
        truth_lab_particles = nullptr;
        truth_breit_particles = nullptr;
      }

      if (reco_source_particles) {
        jet_tools::build_skip_indices(*reco_source_particles, cst_cuts, reco_electron_skip, skip_reco_lab_all);
        skip_reco_breit_all = skip_reco_lab_all;
      } else {
        reco_lab_particles = nullptr;
        reco_breit_particles = nullptr;
      }
    }

    auto cluster_and_write =
        [&](const edm4eic::ReconstructedParticleCollection *particles,
            const std::unordered_set<int> &skip_indices, JetFrame jet_frame,
            const EventKinematics &kin,
            const edm4eic::ReconstructedParticleCollection *lab_particles,
            JetTreeWriter &antikt_writer, JetTreeWriter &cent_writer) {
          if (!particles) {
            return;
          }

          // At this point skip indices should include the electron list and the additional particle list.
          const std::unordered_set<int> *skip_ptr = skip_indices.empty() ? nullptr : &skip_indices;

          const auto inputs = jet_tools::build_input_pseudojets(*particles, skip_ptr);

          auto antikt = jet_tools::build_simplejet_from_pseudojets(
              inputs, jetdef_antikt, particles, native_max_eta, native_min_jet_pT);
              
          auto centauro = jet_tools::build_simplejet_from_pseudojets(
              inputs, jetdef_centauro, particles, native_max_eta, native_min_jet_pT);

          // Native mode already applies jet pT/eta in the builder, so no need for an additional jet-level
          // filter.
          if (args.threshold_frame != ThresholdFrame::Native) {
            antikt = filter_threshold_frame(
                antikt, jet_frame, args.threshold_frame, args.min_jet_pT,
                args.max_eta, kin.have_transform, kin.boost, kin.rot);
            centauro = filter_threshold_frame(
                centauro, jet_frame, args.threshold_frame, args.min_jet_pT,
                args.max_eta, kin.have_transform, kin.boost, kin.rot);
          }

          attach_z_inv(antikt, jet_frame, kin, lab_particles);
          attach_z_inv(centauro, jet_frame, kin, lab_particles);
          antikt_writer.fill_many(event_id, antikt, kin.Q2);
          cent_writer.fill_many(event_id, centauro, kin.Q2);
        };

    cluster_and_write(truth_lab_particles, skip_truth_lab_all, JetFrame::Lab, truth_kin, truth_lab_particles, 
                      t_lab_antikt_truth, t_lab_cent_truth);
    cluster_and_write(truth_breit_particles, skip_truth_breit_all, JetFrame::Breit, truth_kin, truth_lab_particles,
                      t_breit_antikt_truth, t_breit_cent_truth);
    cluster_and_write(reco_lab_particles, skip_reco_lab_all, JetFrame::Lab, reco_kin, reco_lab_particles, 
                      t_lab_antikt_reco, t_lab_cent_reco);
    cluster_and_write(reco_breit_particles, skip_reco_breit_all, JetFrame::Breit, reco_kin, reco_lab_particles,
                      t_breit_antikt_reco, t_breit_cent_reco);
  }

  progress.finish(std::string("[cluster_jets] event ") + std::to_string(entries) +
                  "/" + std::to_string(entries));

  // Write out the trees.
  t_lab_antikt_truth.tree.Write();
  t_lab_cent_truth.tree.Write();
  t_breit_antikt_truth.tree.Write();
  t_breit_cent_truth.tree.Write();
  t_lab_antikt_reco.tree.Write();
  t_lab_cent_reco.tree.Write();
  t_breit_antikt_reco.tree.Write();
  t_breit_cent_reco.tree.Write();

  fout.Close();
  return 0;
}

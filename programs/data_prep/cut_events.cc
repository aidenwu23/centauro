// Apply truth-based DIS cuts (Q2, y) in addition to Breit frame cuts in order to an EDM4hep/EDM4eic 
// file and write a skim.
/*
Run:
./run.sh programs/data_prep/cut_events.cc -i data/events/filtered/merged_pre.root -o data/events/filtered/merged.root --require-breit
./run.sh programs/data_prep/cut_events.cc -i data/events/filtered/merged_pre.root --list

Always require breit.

=============================================================================
Log
=============================================================================
Basic kinematic cuts such as Q^2 and y are done in this file.

Why these cuts matter (physics view)
- We restrict to DIS events in the perturbative regime (Q^2 >= 100 GeV^2, 0.05 < y < 0.7) so jet
  structure and boosts to the Breit frame are well-defined and comparable to the analysis spec.
- Here y is the DIS inelasticity (energy transferred from the lepton to the hadronic system). Very low
  y gives poor energy resolution and unstable boosts; very high y suffers from acceptance/radiative
  effects. The 0.05-0.7 window keeps kinematics well measured.
- Truth-level DIS kinematics ensure the sample corresponds to the intended process before reconstruction
  effects; this gate keeps backgrounds and low-virtuality tails out of the calibration chain.
- Consistent Breit boosts require a valid scattered electron; events with missing inputs or failing
  Breit invariants/residuals are dropped to avoid corrupting jet kinematics.

What this step does
- Reads EDM4hep/eic frames, applies truth-based DIS gates, and optionally enforces the Breit-frame
  sanity checks; writes a skimmed ROOT for downstream jet finding and calibration.
- Tracks how many events fail each gate (Q^2, y, missing truth, missing Breit inputs, invariants) for
  QA and sanity checks.

*/

#include <edm4eic/InclusiveKinematicsCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4hep/MCParticleCollection.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/ROOTWriter.h>

#include <Math/GenVector/Boost.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzE4D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/Vector3D.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "jet_tools/include/beam_helpers.h"
#include "jet_tools/include/math_helpers.h"

namespace {

struct Args {
  std::string input;
  std::string output;
  bool list_only = false;
  double min_q2 = 100.0;
  double y_min = 0.05;
  double y_max = 0.7;
  bool require_breit = true;
  bool require_min_breit = false;
};

struct Counters {
  std::size_t processed = 0;
  std::size_t written = 0;
  std::size_t fail_q2 = 0;
  std::size_t fail_y = 0;
  std::size_t missing_truth = 0;
  std::size_t miss_breit_inputs = 0;
  std::size_t fail_breit_finite = 0;
  std::size_t fail_breit_resid = 0;
  std::size_t fail_breit_inv = 0;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root [-o OUTPUT.root] [--list]\n"
            << "       [--min-q2 GEVC2] [--y-min VAL] [--y-max VAL]\n"
            << "       [--require-breit] [--require-min-breit] [--no-breit]\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    std::string_view v(argv[i]);
    if ((v == "-i" || v == "--input") && i + 1 < argc) {
      args.input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--list") {
      args.list_only = true;
    } else if (v == "--min-q2" && i + 1 < argc) {
      args.min_q2 = std::atof(argv[++i]);
    } else if (v == "--y-min" && i + 1 < argc) {
      args.y_min = std::atof(argv[++i]);
    } else if (v == "--y-max" && i + 1 < argc) {
      args.y_max = std::atof(argv[++i]);
    } else if (v == "--require-breit") {
      args.require_breit = true;
      args.require_min_breit = false;
    } else if (v == "--require-min-breit") {
      args.require_min_breit = true;
      args.require_breit = false;
    } else if (v == "--no-breit") {
      args.require_breit = false;
      args.require_min_breit = false;
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }
  if (args.input.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (!args.list_only && args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  return args;
}

// Collect names of all available collections in the frame.
std::unordered_set<std::string> collect_names(const podio::Frame& frame) {
  std::unordered_set<std::string> names;
  for (const auto& n : frame.getAvailableCollections()) {
    names.emplace(n);
  }
  return names;
}

// Truth-based DIS gate on Q2 and y using InclusiveKinematicsTruth.
bool passes_kinematics(const podio::Frame& frame, const std::unordered_set<std::string>& names,
                       const Args& args, Counters& counts) {
  if (!names.count("InclusiveKinematicsTruth")) {
    ++counts.missing_truth;
    return false;
  }
  const auto& kin = frame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsTruth");
  if (kin.empty()) {
    ++counts.missing_truth;
    return false;
  }
  const auto q2 = kin[0].getQ2();
  const auto y = kin[0].getY();
  if (!(q2 > args.min_q2)) {
    ++counts.fail_q2;
    return false;
  }
  if (!(y >= args.y_min && y <= args.y_max)) {
    ++counts.fail_y;
    return false;
  }
  return true;
}

using P4 = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using Vec3 = ROOT::Math::XYZVector;
using jet_tools::round_beam;
using jet_tools::minkowski_dot;

// Snap a beam momentum onto nominal rails, rotate by the crossing angle, and build a four-vector.

edm4hep::MCParticle find_beam(const edm4hep::MCParticleCollection& mc,
                              const std::vector<int>& pdgs) {
  for (const auto& p : mc) {
    if (p.getGeneratorStatus() == 4) {
      for (int code : pdgs) {
        if (p.getPDG() == code) return p;
      }
    }
  }
  return edm4hep::MCParticle();
}

bool passes_breit(const podio::Frame& frame, const std::unordered_set<std::string>& names,
                  Counters& counts) {
  // Thresholds mirrored from macros/check_breit.C
  constexpr double kBjorkenMax = 5e-4;
  constexpr double kQ2Max = 1e-2;          // GeV^2
  constexpr double kPhotonETol = 1e-3;     // GeV
  constexpr double kPhotonTTol = 1e-3;     // GeV
  constexpr double kProtonPTTol = 1e-3;    // GeV
  constexpr double kProtonPzMax = 1e-2;    // GeV
  constexpr double me = 0.000511;          // GeV
  constexpr double mp = 0.938272;          // GeV
  constexpr double crossing = -0.025;      // rad
  const std::array<double, 3> eNomPz{-5.0, -10.0, -18.0};
  const std::array<double, 3> hNomPz{41.0, 100.0, 275.0};

  if (!names.count("InclusiveKinematicsElectron") || !names.count("MCParticles")) {
    ++counts.miss_breit_inputs;
    return false;
  }

  const auto& kin = frame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsElectron");
  const auto& mc = frame.get<edm4hep::MCParticleCollection>("MCParticles");
  if (kin.empty() || mc.empty()) {
    ++counts.miss_breit_inputs;
    return false;
  }

  const auto& k = kin[0];
  const double x = k.getX();
  const double Q2 = k.getQ2();
  auto scat = k.getScat();
  if (!scat.isAvailable()) {
    ++counts.miss_breit_inputs;
    return false;
  }

  const auto eInMC = find_beam(mc, {11});
  const auto pInMC = find_beam(mc, {2212, 2112});
  if (!eInMC.isAvailable() || !pInMC.isAvailable()) {
    ++counts.miss_breit_inputs;
    return false;
  }

  // Build incoming/outgoing lepton and proton four-vectors.
  const auto& eMomIn = eInMC.getMomentum();
  const auto& pMomIn = pInMC.getMomentum();
  P4 e_i = round_beam(Vec3(eMomIn.x, eMomIn.y, eMomIn.z), me, eNomPz, 0.0);
  P4 p_i = round_beam(Vec3(pMomIn.x, pMomIn.y, pMomIn.z), mp, hNomPz, crossing);
  P4 e_f;
  const auto& scatMom = scat.getMomentum();
  e_f.SetPxPyPzE(scatMom.x, scatMom.y, scatMom.z, scat.getEnergy());

  P4 q = e_i - e_f;
  const double Q2_ref_candidate = -q.M2();
  const double denom_ref = 2.0 * minkowski_dot(p_i, q);
  if (!(Q2_ref_candidate > 0.0) || std::abs(denom_ref) < 1e-12) {
    ++counts.fail_breit_inv;
    return false;
  }
  const double x_ref = Q2_ref_candidate / denom_ref;
  const double Q2_ref = Q2_ref_candidate;

  const auto P3 = p_i.Vect();
  const auto q3 = q.Vect();
  const double denom = (2.0 * x * p_i.E() + q.E());
  if (std::abs(denom) < 1e-9) {
    ++counts.fail_breit_inv;
    return false;
  }

  ROOT::Math::Boost boost(-(2.0 * x * P3 + q3) / denom);
  P4 p_i_b = boost * p_i;
  P4 e_i_b = boost * e_i;
  P4 e_f_b = boost * e_f;
  P4 q_b = boost * q;

  // Rotate so the virtual photon sits on -z.
  auto zhat = -q_b.Vect().Unit();
  auto yhat = e_i_b.Vect().Cross(e_f_b.Vect());
  const double ymag = yhat.R();
  if (ymag > 0) {
    yhat = yhat * (1.0 / ymag);
  } else {
    auto tmp = Vec3(0, 1, 0);
    if (std::abs(zhat.Dot(tmp)) > 0.9) tmp = Vec3(1, 0, 0);
    yhat = (tmp.Cross(zhat)).Unit();
  }
  auto xhat = yhat.Cross(zhat);
  ROOT::Math::Rotation3D Rinv(xhat, yhat, zhat);
  auto R = Rinv.Inverse();

  p_i_b = R * p_i_b;
  e_i_b = R * e_i_b;
  e_f_b = R * e_f_b;
  q_b = R * q_b;

  // Residual gate
  const double qE = std::abs(q_b.E());
  const double qTx = q_b.Px(), qTy = q_b.Py();
  const double qT = std::sqrt(qTx * qTx + qTy * qTy);
  const double sqrtQ2 = Q2 > 0.0 ? std::sqrt(Q2) : 0.0;
  const double qz_expected = -sqrtQ2;
  const double qz_dev = std::abs(q_b.Pz() - qz_expected);
  const double pT = std::sqrt(p_i_b.Px() * p_i_b.Px() + p_i_b.Py() * p_i_b.Py());
  const double pz_expected = (sqrtQ2 > 0.0 && x > 0.0) ? sqrtQ2 / (2.0 * x) : 0.0;
  const double pz_dev = std::abs(p_i_b.Pz() - pz_expected);

  const bool breit_ok = (qE <= kPhotonETol && qT <= kPhotonTTol &&
                         pT <= kProtonPTTol && pz_dev <= kProtonPzMax);
  if (!breit_ok) {
    ++counts.fail_breit_resid;
    return false;
  }

  // Invariant gate (x, Q2) comparing reconstructed to Breit recomputation.
  P4 q_b_actual = e_i_b - e_f_b;
  const double Q2_breit_candidate = -q_b_actual.M2();
  const double denom_breit = 2.0 * minkowski_dot(p_i_b, q_b_actual);
  if (!(Q2_breit_candidate > 0.0) || std::abs(denom_breit) < 1e-12) {
    ++counts.fail_breit_inv;
    return false;
  }
  const double x_breit = Q2_breit_candidate / denom_breit;
  const double Q2_breit = Q2_breit_candidate;
  const double dx = std::abs(x_ref - x_breit);
  const double dQ2 = std::abs(Q2_ref - Q2_breit);
  const bool inv_ok = (dx <= kBjorkenMax && dQ2 <= kQ2Max);
  if (!inv_ok) {
    ++counts.fail_breit_inv;
    return false;
  }

  return true;
}

bool passes_breit_min(const podio::Frame& frame, const std::unordered_set<std::string>& names,
                      Counters& counts) {
  constexpr double me = 0.000511;          // GeV
  constexpr double mp = 0.938272;          // GeV
  constexpr double crossing = -0.025;      // rad
  const std::array<double, 3> eNomPz{-5.0, -10.0, -18.0};
  const std::array<double, 3> hNomPz{41.0, 100.0, 275.0};

  if (!names.count("InclusiveKinematicsElectron") || !names.count("MCParticles")) {
    ++counts.miss_breit_inputs;
    return false;
  }

  const auto& kin = frame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsElectron");
  const auto& mc = frame.get<edm4hep::MCParticleCollection>("MCParticles");
  if (kin.empty() || mc.empty()) {
    ++counts.miss_breit_inputs;
    return false;
  }

  const auto& k = kin[0];
  const double x = k.getX();
  const double Q2 = k.getQ2();
  auto scat = k.getScat();
  if (!scat.isAvailable()) {
    ++counts.miss_breit_inputs;
    return false;
  }
  if (!std::isfinite(x) || !std::isfinite(Q2)) {
    ++counts.fail_breit_finite;
    return false;
  }

  const auto eInMC = find_beam(mc, {11});
  const auto pInMC = find_beam(mc, {2212, 2112});
  if (!eInMC.isAvailable() || !pInMC.isAvailable()) {
    ++counts.miss_breit_inputs;
    return false;
  }

  const auto& eMomIn = eInMC.getMomentum();
  const auto& pMomIn = pInMC.getMomentum();
  P4 e_i = round_beam(Vec3(eMomIn.x, eMomIn.y, eMomIn.z), me, eNomPz, 0.0);
  P4 p_i = round_beam(Vec3(pMomIn.x, pMomIn.y, pMomIn.z), mp, hNomPz, crossing);
  const auto& scatMom = scat.getMomentum();
  P4 e_f;
  e_f.SetPxPyPzE(scatMom.x, scatMom.y, scatMom.z, scat.getEnergy());

  if (!std::isfinite(e_i.E()) || !std::isfinite(p_i.E()) || !std::isfinite(e_f.E())) {
    ++counts.fail_breit_finite;
    return false;
  }

  P4 q = e_i - e_f;
  const double Q2_ref_candidate = -q.M2();
  const double denom_ref = 2.0 * minkowski_dot(p_i, q);
  if (!std::isfinite(Q2_ref_candidate) || !std::isfinite(denom_ref) ||
      std::abs(denom_ref) < 1e-12) {
    ++counts.fail_breit_finite;
    return false;
  }
  const double x_ref = Q2_ref_candidate / denom_ref;
  if (!std::isfinite(x_ref)) {
    ++counts.fail_breit_finite;
    return false;
  }

  const auto P3 = p_i.Vect();
  const auto q3 = q.Vect();
  const double denom = (2.0 * x * p_i.E() + q.E());
  if (!std::isfinite(denom) || std::abs(denom) < 1e-9) {
    ++counts.fail_breit_finite;
    return false;
  }

  ROOT::Math::Boost boost(-(2.0 * x * P3 + q3) / denom);
  P4 p_i_b = boost * p_i;
  P4 e_i_b = boost * e_i;
  P4 e_f_b = boost * e_f;
  P4 q_b = boost * q;

  auto zhat = -q_b.Vect().Unit();
  auto yhat = e_i_b.Vect().Cross(e_f_b.Vect());
  const double ymag = yhat.R();
  if (ymag > 0) {
    yhat = yhat * (1.0 / ymag);
  } else {
    auto tmp = Vec3(0, 1, 0);
    if (std::abs(zhat.Dot(tmp)) > 0.9) tmp = Vec3(1, 0, 0);
    yhat = (tmp.Cross(zhat)).Unit();
  }
  auto xhat = yhat.Cross(zhat);
  ROOT::Math::Rotation3D Rinv(xhat, yhat, zhat);
  auto R = Rinv.Inverse();

  p_i_b = R * p_i_b;
  e_i_b = R * e_i_b;
  e_f_b = R * e_f_b;
  q_b = R * q_b;

  const double qE = q_b.E();
  const double qTx = q_b.Px();
  const double qTy = q_b.Py();
  const double pT = std::sqrt(p_i_b.Px() * p_i_b.Px() + p_i_b.Py() * p_i_b.Py());
  const double pz = p_i_b.Pz();
  if (!std::isfinite(qE) || !std::isfinite(qTx) || !std::isfinite(qTy) ||
      !std::isfinite(pT) || !std::isfinite(pz)) {
    ++counts.fail_breit_finite;
    return false;
  }

  return true;
}

}  // namespace

int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);

  // Open the input file.
  podio::ROOTReader reader;
  try {
    reader.openFile(args.input);
  } catch (const std::exception& ex) {
    std::cerr << "[cut_events] Failed to open input: " << ex.what() << "\n";
    return 1;
  }

  const auto entries = reader.getEntries("events");
  if (entries == 0) {
    std::cerr << "[cut_events] No events in input.\n";
    return 1;
  }

  // List mode: print collections and exit.
  if (args.list_only) {
    auto data = reader.readEntry("events", 0);
    if (!data) {
      std::cerr << "[cut_events] Failed to read first event for listing.\n";
      return 1;
    }
    podio::Frame frame(std::move(data));
    std::cout << "[cut_events] Collections in input:\n";
    for (const auto& name : frame.getAvailableCollections()) {
      std::cout << "  " << name << "\n";
    }
    return 0;
  }

  podio::ROOTWriter writer(args.output);
  Counters counts;

  // Event loop with truth-based cuts.
  for (std::size_t i = 0; i < entries; ++i) {
    ++counts.processed;
    auto data = reader.readEntry("events", i);
    if (!data) continue;
    podio::Frame frame(std::move(data));
    const auto names = collect_names(frame);

    // Truth DIS gate first, optional Breit gate second.
    if (!passes_kinematics(frame, names, args, counts)) continue;
    if (args.require_breit && !passes_breit(frame, names, counts)) continue;
    if (args.require_min_breit && !passes_breit_min(frame, names, counts)) continue;

    writer.writeFrame(frame, "events");
    ++counts.written;
  }

  writer.finish();

  std::cout << "[cut_events] Processed " << counts.processed << " events\n";
  std::cout << "[cut_events] Wrote " << counts.written << " events to " << args.output << "\n";
  std::cout << "[cut_events] Fail Q2: " << counts.fail_q2 << "\n";
  std::cout << "[cut_events] Fail y: " << counts.fail_y << "\n";
  std::cout << "[cut_events] Missing truth kinematics: " << counts.missing_truth << "\n";
  if (args.require_breit || args.require_min_breit) {
    std::cout << "[cut_events] Missing Breit inputs: " << counts.miss_breit_inputs << "\n";
    std::cout << "[cut_events] Fail Breit finite: " << counts.fail_breit_finite << "\n";
    if (args.require_breit) {
      std::cout << "[cut_events] Fail Breit residuals: " << counts.fail_breit_resid << "\n";
      std::cout << "[cut_events] Fail Breit invariants: " << counts.fail_breit_inv << "\n";
    }
  }

  return 0;
}

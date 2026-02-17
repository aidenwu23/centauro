/*
Run:
./run.sh programs/same_frame/performance.cc \
  -m data/jets/matched/same_frame_matches.root \
  -j data/jets/raw/jets_0-999_native_cuts.root \
  -o data/graphs/performance.root \
  --abs-eta-max 3 --pT-min 5
*/

#include <Rtypes.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>

#include "jet_tools/include/plot_helpers.h"
#include "jet_tools/include/root_io.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace {

//=====================================================================================================
// Types And Data Models.
//=====================================================================================================
struct Args {
  std::string matches_input;
  std::string jets_input;
  std::string output;
  std::size_t max_events = std::numeric_limits<std::size_t>::max();
  double max_eta = 0.5;
  double min_pT = 5.0;
  bool use_lab = false;
};

struct AlgData {
  std::vector<double> pT_truth;   // Truth pT in Breit frame.
  std::vector<double> rel_resid;  // (reco_pT - truth_pT) / truth_pT
  std::vector<double> resp_raw;   // pT_reco/pT_truth
  std::vector<double> dphi;       // |phi_reco - phi_truth|
  std::vector<double> deta;       // |eta_reco - eta_truth|
  std::vector<double> dR;         // sqrt(dphi^2 + deta^2)
};

using EventJets = jet_tools::EventJets;

//=====================================================================================================
// CLI Helpers.
//=====================================================================================================
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " -m MATCHES.root -j JETS.root -o OUTPUT.root [--max-events N]\n"
            << "       [--abs-eta-max VAL] [--pT-min VAL] [--use-lab]\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view v(argv[i]);
    if ((v == "-m" || v == "--matches-input" || v == "-i" || v == "--input") &&
        i + 1 < argc) {
      args.matches_input = argv[++i];
    } else if ((v == "-j" || v == "--jets-input") && i + 1 < argc) {
      args.jets_input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--max-events" && i + 1 < argc) {
      const long long n = std::stoll(argv[++i]);
      if (n < 0) {
        std::cerr << "[performance] max-events must be >= 0\n";
        std::exit(1);
      }
      args.max_events = static_cast<std::size_t>(n);
    } else if (v == "--abs-eta-max" && i + 1 < argc) {
      args.max_eta = std::atof(argv[++i]);
    } else if (v == "--pT-min" && i + 1 < argc) {
      args.min_pT = std::atof(argv[++i]);
    } else if (v == "--use-lab") {
      args.use_lab = true;
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.matches_input.empty() || args.jets_input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  if (std::isfinite(args.max_eta) && args.max_eta <= 0.0) {
    std::cerr << "[performance] abs-eta-max must be > 0\n";
    std::exit(1);
  }
  if (!std::isfinite(args.min_pT) || args.min_pT < 0.0) {
    std::cerr << "[performance] pT-min must be finite and >= 0\n";
    std::exit(1);
  }
  return args;
}

//=====================================================================================================
// Math And Utility Helpers.
//=====================================================================================================
double delta_phi(double phi1, double phi2) {
  constexpr double pi = 3.14159265358979323846;
  if (!std::isfinite(phi1) || !std::isfinite(phi2)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double dphi = phi1 - phi2;
  while (dphi > pi) {
    dphi -= 2.0 * pi;
  }
  while (dphi <= -pi) {
    dphi += 2.0 * pi;
  }
  return dphi;
}

std::vector<double> finite_copy(const std::vector<double>& values) {
  std::vector<double> out;
  out.reserve(values.size());
  for (const double value : values) {
    if (std::isfinite(value)) {
      out.push_back(value);
    }
  }
  return out;
}

double mean_finite(const std::vector<double>& values) {
  const std::vector<double> finite = finite_copy(values);
  if (finite.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double sum = 0.0;
  for (const double value : finite) {
    sum += value;
  }
  return sum / static_cast<double>(finite.size());
}

//=====================================================================================================
// Statistical Helpers.
//=====================================================================================================
double quantile_sorted(const std::vector<double>& sorted_values, double q) {
  if (sorted_values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (q <= 0.0) {
    return sorted_values.front();
  }
  if (q >= 1.0) {
    return sorted_values.back();
  }
  const double pos = q * static_cast<double>(sorted_values.size() - 1U);
  const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
  if (lo == hi) {
    return sorted_values[lo];
  }
  const double frac = pos - static_cast<double>(lo);
  return (1.0 - frac) * sorted_values[lo] + frac * sorted_values[hi];
}

double effSigma68(std::vector<double> values) {
  values = finite_copy(values);
  if (values.size() < 2U) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::sort(values.begin(), values.end());

  // Total length and 68% length (rounded up).
  const std::size_t n = values.size();
  const std::size_t n68 = static_cast<std::size_t>(std::ceil(0.68 * static_cast<double>(n)));
  if (n68 < 2U || n68 > n) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Find the shortest interval that contains 68% of the data.
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i + n68 - 1U < n; ++i) {
    const double width = values[i + n68 - 1U] - values[i];
    if (width < best) {
      best = width;
    }
  }
  if (!std::isfinite(best)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Return half-width.
  return 0.5 * best;
}

double effSigma68_err_bootstrap(const std::vector<double>& values, std::uint32_t seed,
                                int n_samples = 300) {
  return jet_tools::bootstrap_err_stat(values, seed, &effSigma68, n_samples);
}

int find_bin(const std::vector<double>& edges, double value) {
  for (std::size_t i = 0; i + 1U < edges.size(); ++i) {
    const bool is_last = (i + 2U == edges.size());
    if (value >= edges[i] && (value < edges[i + 1U] || (is_last && value <= edges[i + 1U]))) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

bool pass_jet_cuts(const jet_tools::SimpleJet& jet, const Args& args) {
  if (!std::isfinite(jet.pT) || jet.pT < args.min_pT) {
    return false;
  }
  if (!std::isfinite(jet.eta)) {
    return false;
  }
  if (std::isfinite(args.max_eta) && std::abs(jet.eta) > args.max_eta) {
    return false;
  }
  return true;
}

void write_canvas(TDirectory* dir, TCanvas& canvas) {
  if (!dir) {
    return;
  }
  dir->cd();
  canvas.Write();
}

//=====================================================================================================
// Core Data Collection.
//=====================================================================================================
void collect_alg_data(const std::string& label, const EventJets& truth_jets,
                      const EventJets& reco_jets,
                      const std::vector<jet_tools::TruthRecoMatchRow>& matches,
                      const Args& args,
                      AlgData& out) {
  // Create and reserve structure to compute metrics later.
  out = AlgData{};
  out.pT_truth.reserve(matches.size());
  out.rel_resid.reserve(matches.size());
  out.resp_raw.reserve(matches.size());
  out.dphi.reserve(matches.size());
  out.deta.reserve(matches.size());
  out.dR.reserve(matches.size());

  std::size_t last_chars = 0;
  for (std::size_t i = 0; i < matches.size(); ++i) {
    if (i % 100 == 0) {
      const std::string message = std::string("[performance] collect ") + label + " " +
                                  std::to_string(i) + "/" +
                                  std::to_string(matches.size());
      jet_tools::print_progress(message, last_chars);
    }
    const auto& match = matches[i];
    if (match.truth_index < 0 || match.reco_index < 0) {
      continue;
    }

    const auto truth_iterator = truth_jets.find(match.event);
    const auto reco_iterator = reco_jets.find(match.event);
    if (truth_iterator == truth_jets.end() || reco_iterator == reco_jets.end()) {
      continue;
    }

    const std::size_t truth_index = static_cast<std::size_t>(match.truth_index);
    const std::size_t reco_index = static_cast<std::size_t>(match.reco_index);
    if (truth_index >= truth_iterator->second.size() || reco_index >= reco_iterator->second.size()) {
      continue;
    }

    const auto& truth = truth_iterator->second[truth_index];
    const auto& reco = reco_iterator->second[reco_index];
    if (!pass_jet_cuts(truth, args)) {
      continue;
    }
    if (!std::isfinite(truth.pT) || truth.pT <= 0.0 || !std::isfinite(reco.pT)) {
      continue;
    }
    if (!std::isfinite(truth.eta) || !std::isfinite(reco.eta) || !std::isfinite(truth.phi) || 
        !std::isfinite(reco.phi)) {
      continue;
    }

    const double response = reco.pT / truth.pT;
    const double rel_residual = (reco.pT - truth.pT) / truth.pT;
    if (!std::isfinite(response) || !std::isfinite(rel_residual)) {
      continue;
    }

    const double dphi = delta_phi(reco.phi, truth.phi);
    const double deta = reco.eta - truth.eta;
    double pair_dR = match.dR;
    if (!std::isfinite(pair_dR) && std::isfinite(dphi) && std::isfinite(deta)) {
      pair_dR = std::hypot(dphi, deta);
    }

    out.pT_truth.push_back(truth.pT);
    out.rel_resid.push_back(rel_residual);
    out.resp_raw.push_back(response);
    out.dphi.push_back(dphi);
    out.deta.push_back(deta);
    out.dR.push_back(pair_dR);
  }
  const std::string done_message = std::string("[performance] collect ") + label + " " +
                                   std::to_string(matches.size()) + "/" +
                                   std::to_string(matches.size());
  jet_tools::print_progress(done_message, last_chars);
  std::cout << "\n";
}


//=====================================================================================================
// Summary And Output Writers.
//=====================================================================================================
void write_outputs(const Args& args, const AlgData (&data)[2], TFile& fout,
                   TDirectory* hists_dir, TDirectory* raw_dir) {
  const double max_pT_plot = 120.0;
  const std::vector<double> pT_edges = {5.0, 7.0, 8.0, 12.0, 20.0, 60.0};

  //===================================================================================================
  // Histogram Definitions.
  //===================================================================================================
  std::unique_ptr<TH2D> hrel_resid_pT[2];
  std::unique_ptr<TH2D> hresp_raw_pT[2];
  std::unique_ptr<TH1D> hrel_resid_bias[2];
  std::unique_ptr<TH1D> hrel_resid_effsigma68[2];
  std::unique_ptr<TH1D> hresp_mean_raw_pT[2];
  std::unique_ptr<TH1D> hdphi_bias[2];
  std::unique_ptr<TH1D> hdphi_effsigma68[2];
  std::unique_ptr<TH1D> hdeta_bias[2];
  std::unique_ptr<TH1D> hdeta_effsigma68[2];
  std::unique_ptr<TH1D> hdR_median[2];

  //===================================================================================================
  // Histogram Filling And Binning.
  //===================================================================================================
  for (int alg = 0; alg < 2; ++alg) {
    const char* alg_name = (alg == 0) ? "antikt" : "centauro";
    const char* alg_label = (alg == 0) ? "anti-k_{t}" : "Centauro";

    double rel_min = std::numeric_limits<double>::infinity();
    double rel_max = -std::numeric_limits<double>::infinity();
    for (double value : data[alg].rel_resid) {
      if (!std::isfinite(value)) {
        continue;
      }
      rel_min = std::min(rel_min, value);
      rel_max = std::max(rel_max, value);
    }
    if (!std::isfinite(rel_min) || !std::isfinite(rel_max)) {
      rel_min = -1.0;
      rel_max = 1.0;
    }
    if (rel_max <= rel_min) {
      rel_min -= 1.0;
      rel_max += 1.0;
    }
    const double rel_half = 0.5 * (rel_max - rel_min);
    const double rel_center = 0.5 * (rel_max + rel_min);
    const double rel_span = (rel_half > 0.0) ? (1.1 * rel_half) : 1.0;

    hrel_resid_pT[alg] = std::make_unique<TH2D>(
        std::string("rel_resid_pT_vs_pT_truth_").append(alg_name).c_str(),
        (std::string("Relative residual vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        60, args.min_pT, max_pT_plot, 120,
        rel_center - rel_span, rel_center + rel_span);
    hrel_resid_pT[alg]->SetDirectory(nullptr);
    hrel_resid_pT[alg]->SetStats(false);
    hrel_resid_pT[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hrel_resid_pT[alg]->GetYaxis()->SetTitle("(p_{T,Breit}^{raw}-p_{T,Breit}^{truth})/p_{T,Breit}^{truth}");

    hresp_raw_pT[alg] = std::make_unique<TH2D>(
        std::string("resp_raw_pT_vs_pT_truth_").append(alg_name).c_str(),
        (std::string("Response vs truth p_{T,Breit} (raw, ") + alg_label + ")").c_str(),
        60, args.min_pT, max_pT_plot, 120, 0.0, 2.0);
    hresp_raw_pT[alg]->SetDirectory(nullptr);
    hresp_raw_pT[alg]->SetStats(false);
    hresp_raw_pT[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hresp_raw_pT[alg]->GetYaxis()->SetTitle("p_{T,Breit}^{raw}/p_{T,Breit}^{truth}");

    hrel_resid_bias[alg] = std::make_unique<TH1D>(
        std::string("bias_rel_resid_pT_").append(alg_name).c_str(),
        (std::string("Mean relative p_{T,Breit} residual vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hrel_resid_bias[alg]->SetDirectory(nullptr);
    hrel_resid_bias[alg]->SetStats(false);
    hrel_resid_bias[alg]->SetMarkerStyle(20);
    hrel_resid_bias[alg]->SetMarkerSize(1.0);
    hrel_resid_bias[alg]->SetLineWidth(2);

    hrel_resid_effsigma68[alg] = std::make_unique<TH1D>(
        std::string("effSigma68_rel_resid_pT_").append(alg_name).c_str(),
        (std::string("Core relative p_{T,Breit} resolution vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hrel_resid_effsigma68[alg]->SetDirectory(nullptr);
    hrel_resid_effsigma68[alg]->SetStats(false);
    hrel_resid_effsigma68[alg]->SetMarkerStyle(20);
    hrel_resid_effsigma68[alg]->SetMarkerSize(1.0);
    hrel_resid_effsigma68[alg]->SetLineWidth(2);

    hresp_mean_raw_pT[alg] = std::make_unique<TH1D>(
        std::string("mean_resp_raw_pT_").append(alg_name).c_str(),
        (std::string("Response vs truth p_{T,Breit} (raw, ") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hresp_mean_raw_pT[alg]->SetDirectory(nullptr);
    hresp_mean_raw_pT[alg]->SetStats(false);
    hresp_mean_raw_pT[alg]->SetMarkerStyle(20);
    hresp_mean_raw_pT[alg]->SetMarkerSize(1.0);
    hresp_mean_raw_pT[alg]->SetLineWidth(2);

    hdphi_bias[alg] = std::make_unique<TH1D>(
        std::string("bias_dphi_pT_").append(alg_name).c_str(),
        (std::string("Mean #Delta#phi_{Breit} vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hdphi_bias[alg]->SetDirectory(nullptr);
    hdphi_bias[alg]->SetStats(false);
    hdphi_bias[alg]->SetMarkerStyle(20);
    hdphi_bias[alg]->SetMarkerSize(1.0);
    hdphi_bias[alg]->SetLineWidth(2);

    hdphi_effsigma68[alg] = std::make_unique<TH1D>(
        std::string("effSigma68_dphi_pT_").append(alg_name).c_str(),
        (std::string("#Delta#phi_{Breit} core width vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hdphi_effsigma68[alg]->SetDirectory(nullptr);
    hdphi_effsigma68[alg]->SetStats(false);
    hdphi_effsigma68[alg]->SetMarkerStyle(20);
    hdphi_effsigma68[alg]->SetMarkerSize(1.0);
    hdphi_effsigma68[alg]->SetLineWidth(2);

    hdeta_bias[alg] = std::make_unique<TH1D>(
        std::string("bias_deta_pT_").append(alg_name).c_str(),
        (std::string("Mean #Delta#eta_{Breit} vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hdeta_bias[alg]->SetDirectory(nullptr);
    hdeta_bias[alg]->SetStats(false);
    hdeta_bias[alg]->SetMarkerStyle(20);
    hdeta_bias[alg]->SetMarkerSize(1.0);
    hdeta_bias[alg]->SetLineWidth(2);

    hdeta_effsigma68[alg] = std::make_unique<TH1D>(
        std::string("effSigma68_deta_pT_").append(alg_name).c_str(),
        (std::string("#Delta#eta_{Breit} core width vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hdeta_effsigma68[alg]->SetDirectory(nullptr);
    hdeta_effsigma68[alg]->SetStats(false);
    hdeta_effsigma68[alg]->SetMarkerStyle(20);
    hdeta_effsigma68[alg]->SetMarkerSize(1.0);
    hdeta_effsigma68[alg]->SetLineWidth(2);

    hdR_median[alg] = std::make_unique<TH1D>(
        std::string("median_dR_pT_").append(alg_name).c_str(),
        (std::string("Median #DeltaR_{Breit} vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
        static_cast<int>(pT_edges.size() - 1U), pT_edges.data());
    hdR_median[alg]->SetDirectory(nullptr);
    hdR_median[alg]->SetStats(false);
    hdR_median[alg]->SetMarkerStyle(20);
    hdR_median[alg]->SetMarkerSize(1.0);
    hdR_median[alg]->SetLineWidth(2);

    std::vector<std::vector<double>> rel_resid_bins(pT_edges.size() - 1U);
    std::vector<std::vector<double>> resp_raw_bins(pT_edges.size() - 1U);
    std::vector<std::vector<double>> dphi_bins(pT_edges.size() - 1U);
    std::vector<std::vector<double>> deta_bins(pT_edges.size() - 1U);
    std::vector<std::vector<double>> dR_bins(pT_edges.size() - 1U);

    const std::size_t n_points = data[alg].pT_truth.size();
    for (std::size_t i = 0; i < n_points; ++i) {
      const double pT_truth = data[alg].pT_truth[i];
      const double rel_resid = data[alg].rel_resid[i];
      const double resp_raw = data[alg].resp_raw[i];
      const double dphi = data[alg].dphi[i];
      const double deta = data[alg].deta[i];
      const double dR = data[alg].dR[i];

      if (std::isfinite(rel_resid)) {
        hrel_resid_pT[alg]->Fill(pT_truth, rel_resid);
      }

      if (std::isfinite(resp_raw)) {
        hresp_raw_pT[alg]->Fill(pT_truth, resp_raw);
      }

      const int pT_bin = find_bin(pT_edges, pT_truth);
      if (pT_bin < 0) {
        continue;
      }
      const std::size_t bin = static_cast<std::size_t>(pT_bin);

      if (std::isfinite(rel_resid)) {
        rel_resid_bins[bin].push_back(rel_resid);
      }

      if (std::isfinite(resp_raw)) {
        resp_raw_bins[bin].push_back(resp_raw);
      }

      if (std::isfinite(dphi)) {
        dphi_bins[bin].push_back(dphi);
      }

      if (std::isfinite(deta)) {
        deta_bins[bin].push_back(deta);
      }

      if (std::isfinite(dR)) {
        dR_bins[bin].push_back(dR);
      }
    }

    const std::uint32_t seed = 424242U + static_cast<std::uint32_t>(alg * 1337);

    for (std::size_t bin = 0; bin + 1U < pT_edges.size(); ++bin) {
      const int root_bin = static_cast<int>(bin + 1U);

      const double mean_rel = mean_finite(rel_resid_bins[bin]);
      hrel_resid_bias[alg]->SetBinContent(root_bin, std::isfinite(mean_rel) ? mean_rel : 0.0);
      hrel_resid_bias[alg]->SetBinError(
          root_bin,
          jet_tools::bootstrap_mean_err(rel_resid_bins[bin], seed + static_cast<std::uint32_t>(bin)));

      double eff_rel = 0.0;
      double eff_rel_err = 0.0;
      if (finite_copy(rel_resid_bins[bin]).size() >= 2U) {
        eff_rel = effSigma68(rel_resid_bins[bin]);
        eff_rel_err = effSigma68_err_bootstrap(
            rel_resid_bins[bin], seed + static_cast<std::uint32_t>(bin));
      }
      hrel_resid_effsigma68[alg]->SetBinContent(root_bin, std::isfinite(eff_rel) ? eff_rel : 0.0);
      hrel_resid_effsigma68[alg]->SetBinError(root_bin, eff_rel_err);

      const double mean_resp = mean_finite(resp_raw_bins[bin]);
      hresp_mean_raw_pT[alg]->SetBinContent(root_bin, std::isfinite(mean_resp) ? mean_resp : 0.0);
      hresp_mean_raw_pT[alg]->SetBinError(root_bin, jet_tools::bootstrap_mean_err(resp_raw_bins[bin], 
          seed + static_cast<std::uint32_t>(bin + 5000U)));

      const double mean_dphi = mean_finite(dphi_bins[bin]);
      hdphi_bias[alg]->SetBinContent(root_bin, std::isfinite(mean_dphi) ? mean_dphi : 0.0);
      hdphi_bias[alg]->SetBinError(root_bin, jet_tools::bootstrap_mean_err(dphi_bins[bin], 
          seed + static_cast<std::uint32_t>(bin + 10000U)));

      double eff_dphi = 0.0;
      double eff_dphi_err = 0.0;
      if (finite_copy(dphi_bins[bin]).size() >= 2U) {
        eff_dphi = effSigma68(dphi_bins[bin]);
        eff_dphi_err = effSigma68_err_bootstrap(
            dphi_bins[bin], seed + static_cast<std::uint32_t>(bin + 1000U));
      }
      hdphi_effsigma68[alg]->SetBinContent(root_bin, std::isfinite(eff_dphi) ? eff_dphi : 0.0);
      hdphi_effsigma68[alg]->SetBinError(root_bin, eff_dphi_err);

      const double mean_deta = mean_finite(deta_bins[bin]);
      hdeta_bias[alg]->SetBinContent(root_bin, std::isfinite(mean_deta) ? mean_deta : 0.0);
      hdeta_bias[alg]->SetBinError( root_bin,
          jet_tools::bootstrap_mean_err(deta_bins[bin], seed + static_cast<std::uint32_t>(bin + 20000U)));

      double eff_deta = 0.0;
      double eff_deta_err = 0.0;
      if (finite_copy(deta_bins[bin]).size() >= 2U) {
        eff_deta = effSigma68(deta_bins[bin]);
        eff_deta_err = effSigma68_err_bootstrap(
            deta_bins[bin], seed + static_cast<std::uint32_t>(bin + 2000U));
      }
      hdeta_effsigma68[alg]->SetBinContent(root_bin, std::isfinite(eff_deta) ? eff_deta : 0.0);
      hdeta_effsigma68[alg]->SetBinError(root_bin, eff_deta_err);

      double med_dR = 0.0;
      if (!dR_bins[bin].empty()) {
        auto sorted = finite_copy(dR_bins[bin]);
        if (!sorted.empty()) {
          std::sort(sorted.begin(), sorted.end());
          med_dR = quantile_sorted(sorted, 0.5);
        }
      }
      hdR_median[alg]->SetBinContent(root_bin, std::isfinite(med_dR) ? med_dR : 0.0);
      hdR_median[alg]->SetBinError(root_bin, 0.0);
    }
  }

  //===================================================================================================
  // Summary Tree.
  //===================================================================================================
  fout.cd();
  TTree summary("summary", "paper summary");
  int s_alg = 0;
  double s_mean_resp_raw = 0.0;
  double s_effsigma68_rel_resid = 0.0;

  summary.Branch("alg", &s_alg);
  summary.Branch("mean_resp_raw", &s_mean_resp_raw);
  summary.Branch("effSigma68_rel_resid", &s_effsigma68_rel_resid);
  for (int alg = 0; alg < 2; ++alg) {
    s_alg = alg;
    s_mean_resp_raw = mean_finite(data[alg].resp_raw);
    s_effsigma68_rel_resid = effSigma68(data[alg].rel_resid);
    summary.Fill();
  }
  summary.Write();

  //===================================================================================================
  // Histogram Writes.
  //===================================================================================================
  if (hists_dir) {
    hists_dir->cd();
  } else {
    fout.cd();
  }
  for (int alg = 0; alg < 2; ++alg) {
    hrel_resid_bias[alg]->Write();
    hresp_mean_raw_pT[alg]->Write();
    hdphi_bias[alg]->Write();
    hdeta_bias[alg]->Write();
    hrel_resid_effsigma68[alg]->Write();
    hdphi_effsigma68[alg]->Write();
    hdeta_effsigma68[alg]->Write();
    hdR_median[alg]->Write();
    hrel_resid_pT[alg]->Write();
    hresp_raw_pT[alg]->Write();
  }

  //===================================================================================================
  // Canvas Writes.
  //===================================================================================================
  fout.cd();

  for (int alg = 0; alg < 2; ++alg) {
    const char* alg_name = (alg == 0) ? "antikt" : "centauro";
    const char* alg_label = (alg == 0) ? "anti-k_{t}" : "Centauro";

    TCanvas c_rel_resid_bias(std::string("c_bias_rel_resid_pT_").append(alg_name).c_str(),
                             (std::string("Mean relative p_{T,Breit} residual vs truth p_{T,Breit} (") +
                              alg_label + ")")
                                 .c_str(),
                             800, 600);
    hrel_resid_bias[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hrel_resid_bias[alg]->GetYaxis()->SetTitle(
        "mean((p_{T,Breit}^{raw}-p_{T,Breit}^{truth})/p_{T,Breit}^{truth})");
    hrel_resid_bias[alg]->Draw("E1");
    write_canvas(raw_dir, c_rel_resid_bias);

    TCanvas c_resp_raw_mean(std::string("c_mean_resp_raw_pT_").append(alg_name).c_str(),
                            (std::string("Response vs truth p_{T,Breit} (raw, ") + alg_label + ")").c_str(),
                            800, 600);
    hresp_mean_raw_pT[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hresp_mean_raw_pT[alg]->GetYaxis()->SetTitle("mean(p_{T,Breit}^{raw}/p_{T,Breit}^{truth})");
    hresp_mean_raw_pT[alg]->Draw("E1");
    write_canvas(raw_dir, c_resp_raw_mean);

    TCanvas c_dphi_bias(std::string("c_bias_dphi_pT_").append(alg_name).c_str(),
                        (std::string("Mean #Delta#phi_{Breit} vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
                        800, 600);
    hdphi_bias[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hdphi_bias[alg]->GetYaxis()->SetTitle("mean(#Delta#phi_{Breit})");
    hdphi_bias[alg]->Draw("E1");
    write_canvas(raw_dir, c_dphi_bias);

    TCanvas c_deta_bias(std::string("c_bias_deta_pT_").append(alg_name).c_str(),
                        (std::string("Mean #Delta#eta_{Breit} vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
                        800, 600);
    hdeta_bias[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hdeta_bias[alg]->GetYaxis()->SetTitle("mean(#Delta#eta_{Breit})");
    hdeta_bias[alg]->Draw("E1");
    write_canvas(raw_dir, c_deta_bias);

    TCanvas c_rel_resid_rms(std::string("c_effSigma68_rel_resid_pT_").append(alg_name).c_str(),
                            (std::string("Core relative p_{T,Breit} resolution vs truth p_{T,Breit} (") +
                             alg_label + ")")
                                .c_str(),
                            800, 600);
    hrel_resid_effsigma68[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hrel_resid_effsigma68[alg]->GetYaxis()->SetTitle(
        "#sigma_{68}^{eff}((p_{T,Breit}^{raw}-p_{T,Breit}^{truth})/p_{T,Breit}^{truth})");
    hrel_resid_effsigma68[alg]->Draw("E1");
    write_canvas(raw_dir, c_rel_resid_rms);

    TCanvas c_dphi_rms(std::string("c_effSigma68_dphi_pT_").append(alg_name).c_str(),
                       (std::string("#Delta#phi_{Breit} core width vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
                       800, 600);
    hdphi_effsigma68[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hdphi_effsigma68[alg]->GetYaxis()->SetTitle("#sigma_{68}^{eff}(#Delta#phi_{Breit})");
    hdphi_effsigma68[alg]->Draw("E1");
    write_canvas(raw_dir, c_dphi_rms);

    TCanvas c_deta_rms(std::string("c_effSigma68_deta_pT_").append(alg_name).c_str(),
                       (std::string("#Delta#eta_{Breit} core width vs truth p_{T,Breit} (") + alg_label + ")").c_str(),
                       800, 600);
    hdeta_effsigma68[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hdeta_effsigma68[alg]->GetYaxis()->SetTitle("#sigma_{68}^{eff}(#Delta#eta_{Breit})");
    hdeta_effsigma68[alg]->Draw("E1");
    write_canvas(raw_dir, c_deta_rms);

    TCanvas c_dR_med(std::string("c_median_dR_pT_").append(alg_name).c_str(),
                     (std::string("#DeltaR_{Breit} vs Truth p_{T,Breit} (Median, ") + alg_label + ")").c_str(),
                     800, 600);
    hdR_median[alg]->GetXaxis()->SetTitle("p_{T,Breit}^{truth} [GeV]");
    hdR_median[alg]->GetYaxis()->SetTitle("median(#DeltaR_{Breit})");
    hdR_median[alg]->Draw("E1");
    write_canvas(raw_dir, c_dR_med);
  }

  TCanvas c_rel_resid_pT_both("c_rel_resid_pT_both", "Relative residual vs truth p_{T,Breit} (side-by-side)",
                              1200, 600);
  c_rel_resid_pT_both.Divide(2, 1);
  if (auto* pad = c_rel_resid_pT_both.cd(1)) {
    pad->SetRightMargin(0.12);
    hrel_resid_pT[0]->Draw("COLZ");
  }
  if (auto* pad = c_rel_resid_pT_both.cd(2)) {
    pad->SetRightMargin(0.12);
    hrel_resid_pT[1]->Draw("COLZ");
  }
  write_canvas(raw_dir, c_rel_resid_pT_both);

  TCanvas c_resp_raw_pT_both("c_resp_raw_pT_both", "Response vs truth p_{T,Breit} (raw, side-by-side)",
                             1200, 600);
  c_resp_raw_pT_both.Divide(2, 1);
  if (auto* pad = c_resp_raw_pT_both.cd(1)) {
    pad->SetRightMargin(0.12);
    hresp_raw_pT[0]->Draw("COLZ");
  }
  if (auto* pad = c_resp_raw_pT_both.cd(2)) {
    pad->SetRightMargin(0.12);
    hresp_raw_pT[1]->Draw("COLZ");
  }
  write_canvas(raw_dir, c_resp_raw_pT_both);
}

}  // namespace

//=====================================================================================================
// Program Entry Point.
//=====================================================================================================
int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);
  gStyle->SetOptStat(1110);
  jet_tools::ensure_parent(args.output);

  TFile f_matches(args.matches_input.c_str(), "READ");
  if (!f_matches.IsOpen()) {
    std::cerr << "[performance] failed to open matches input " << args.matches_input << "\n";
    return 1;
  }

  TFile f_jets(args.jets_input.c_str(), "READ");
  if (!f_jets.IsOpen()) {
    std::cerr << "[performance] failed to open jets input " << args.jets_input << "\n";
    return 1;
  }

  const char* antikt_truth_name = args.use_lab ? "LabFrameAntiktTruth" : "BreitFrameAntiktTruth";
  const char* antikt_reco_name = args.use_lab ? "LabFrameAntiktReco" : "BreitFrameAntiktReco";
  const char* centauro_truth_name = args.use_lab ? "LabFrameCentauroTruth" : "BreitFrameCentauroTruth";
  const char* centauro_reco_name = args.use_lab ? "LabFrameCentauroReco" : "BreitFrameCentauroReco";

  auto* t_antikt_truth = jet_tools::get_required_tree(f_jets, antikt_truth_name, "performance");
  auto* t_antikt_reco = jet_tools::get_required_tree(f_jets, antikt_reco_name, "performance");
  auto* t_centauro_truth = jet_tools::get_required_tree(f_jets, centauro_truth_name, "performance");
  auto* t_centauro_reco = jet_tools::get_required_tree(f_jets, centauro_reco_name, "performance");

  auto* t_antikt_matches = jet_tools::get_required_tree(f_matches, "antikt_matches", "performance");
  auto* t_centauro_matches = jet_tools::get_required_tree(f_matches, "centauro_matches", "performance");

  EventJets antikt_truth;
  EventJets antikt_reco;
  EventJets centauro_truth;
  EventJets centauro_reco;
  jet_tools::read_jet_tree(*t_antikt_truth, args.max_events, antikt_truth, "performance");
  jet_tools::read_jet_tree(*t_antikt_reco, args.max_events, antikt_reco, "performance");
  jet_tools::read_jet_tree(*t_centauro_truth, args.max_events, centauro_truth, "performance");
  jet_tools::read_jet_tree(*t_centauro_reco, args.max_events, centauro_reco, "performance");

  const std::vector<jet_tools::TruthRecoMatchRow> antikt_matches =
      jet_tools::read_match_tree(*t_antikt_matches, args.max_events, "performance");
  const std::vector<jet_tools::TruthRecoMatchRow> centauro_matches =
      jet_tools::read_match_tree(*t_centauro_matches, args.max_events, "performance");

  AlgData data[2];
  collect_alg_data("antikt", antikt_truth, antikt_reco, antikt_matches, args, data[0]);
  collect_alg_data("centauro", centauro_truth, centauro_reco, centauro_matches, args, data[1]);

  if (data[0].pT_truth.empty() && data[1].pT_truth.empty()) {
    std::cerr << "[performance] no matched jets after selections\n";
    return 1;
  }

  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[performance] failed to open output " << args.output << "\n";
    return 1;
  }

  TDirectory* hists_dir = fout.mkdir("hists");
  TDirectory* raw_dir = fout.mkdir("raw");

  write_outputs(args, data, fout, hists_dir, raw_dir);

  fout.Close();
  f_matches.Close();
  f_jets.Close();
  std::cout << "[performance] wrote " << args.output << "\n";
  return 0;
}

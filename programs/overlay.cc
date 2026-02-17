// Generate overlay canvases from a ROOT file.
/*
Run:
./run.sh programs/overlay.cc -i data/graphs/performance.root -o data/graphs/performance_overlays.root
./run.sh programs/overlay.cc -i data/graphs/controls.root -o data/graphs/controls_overlays.root --no-errors --same-algorithm

*/

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TList.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>

#include "jet_tools/include/plot_helpers.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

struct Args {
  std::string input;
  std::string output;
  bool no_errors = false;
  bool overlay_variants = false;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " -i INPUT.root -o OUTPUT.root [--no-errors]\n"
            << "       [--same-algorithm]\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    std::string_view v(argv[i]);
    if ((v == "-i" || v == "--input") && i + 1 < argc) {
      args.input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else if (v == "--no-errors") {
      args.no_errors = true;
    } else if (v == "--same-algorithm") {
      args.overlay_variants = true;
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }
  if (args.input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  return args;
}

void ensure_parent(const std::string& path) {
  std::filesystem::path p(path);
  if (p.has_parent_path()) {
    std::filesystem::create_directories(p.parent_path());
  }
}

std::unique_ptr<TH1> clone_hist(const TH1& src) {
  auto* clone = static_cast<TH1*>(src.Clone());
  clone->SetDirectory(nullptr);
  return std::unique_ptr<TH1>(clone);
}

std::unique_ptr<TCanvas> clone_canvas(const TCanvas& src) {
  auto* clone = static_cast<TCanvas*>(src.Clone());
  return std::unique_ptr<TCanvas>(clone);
}

void style_hist(TH1& h, int color, int marker) {
  h.SetLineColor(color);
  h.SetMarkerColor(color);
  h.SetMarkerStyle(marker);
  h.SetLineWidth(2);
}

TGraphErrors build_graph(const TH1& h, double x_offset, double x_err_frac) {
  TGraphErrors gr;
  gr.SetName(std::string(h.GetName()).append("_gr").c_str());
  int p = 0;
  for (int i = 1; i <= h.GetNbinsX(); ++i) {
    const double x = h.GetXaxis()->GetBinCenter(i) + x_offset;
    const double y = h.GetBinContent(i);
    const double ey = h.GetBinError(i);
    const double ex = x_err_frac * h.GetXaxis()->GetBinWidth(i);
    gr.SetPoint(p, x, y);
    gr.SetPointError(p, ex, ey);
    ++p;
  }
  return gr;
}

void style_graph(TGraphErrors& g, int color, int marker) {
  g.SetLineColor(color);
  g.SetMarkerColor(color);
  g.SetMarkerStyle(marker);
  g.SetLineWidth(2);
}

jet_tools::AxisRange hist_range(const TH1* a, const TH1* b,
                                double fallback_min, double fallback_max) {
  double min_v = std::numeric_limits<double>::infinity();
  double max_v = -std::numeric_limits<double>::infinity();
  auto scan = [&](const TH1* h) {
    if (!h) return;
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
      const double y = h->GetBinContent(bin);
      const double ey = h->GetBinError(bin);
      if (!std::isfinite(y) || !std::isfinite(ey)) continue;
      min_v = std::min(min_v, y - ey);
      max_v = std::max(max_v, y + ey);
    }
  };
  scan(a);
  scan(b);
  if (!std::isfinite(min_v) || !std::isfinite(max_v)) {
    return {fallback_min, fallback_max};
  }
  if (min_v == max_v) {
    const double pad = std::abs(min_v) > 0.0 ? 0.1 * std::abs(min_v) : 0.1;
    return {min_v - pad, max_v + pad};
  }
  const double pad = 0.1 * (max_v - min_v);
  return {min_v - pad, max_v + pad};
}

std::string strip_once(std::string s, const std::string& needle) {
  const std::size_t pos = s.find(needle);
  if (pos == std::string::npos) return s;
  s.erase(pos, needle.size());
  return s;
}

std::string sanitize_title(std::string title) {
  title = strip_once(std::move(title), " (Raw, anti-k_{t})");
  title = strip_once(std::move(title), " (Corr, anti-k_{t})");
  title = strip_once(std::move(title), " (Calib, anti-k_{t})");
  title = strip_once(std::move(title), " (Raw, anti-k_{T})");
  title = strip_once(std::move(title), " (Corr, anti-k_{T})");
  title = strip_once(std::move(title), " (Calib, anti-k_{T})");
  title = strip_once(std::move(title), " (anti-k_{t})");
  title = strip_once(std::move(title), " (anti-k_{T})");
  title = strip_once(std::move(title), " (Centauro)");
  title = strip_once(std::move(title), ", anti-k_{t}");
  title = strip_once(std::move(title), ", anti-k_{T}");
  title = strip_once(std::move(title), ", Centauro");
  return title;
}

std::string replace_all(std::string s, const std::string& from,
                        const std::string& to) {
  if (from.empty())
    return s;
  std::size_t pos = 0;
  while ((pos = s.find(from, pos)) != std::string::npos) {
    s.replace(pos, from.size(), to);
    pos += to.size();
  }
  return s;
}

std::string normalize_title(std::string title) {
  title = sanitize_title(std::move(title));

  const std::vector<std::pair<std::string, std::string>> map = {
      {"Effective sigma (68%) vs Truth p_{T,Breit} (Relative Residual)",
       "Core relative p_{T,Breit} resolution vs truth p_{T,Breit}"},
      {"Effective sigma (68%) vs Truth p_{T,Breit} (Relative Residual, Calibrated)",
       "Core relative p_{T,Breit} resolution vs truth p_{T,Breit} (calibrated)"},
      {"Bias vs Truth p_{T,Breit} (Relative Residual)",
       "Mean relative p_{T,Breit} residual vs truth p_{T,Breit}"},
      {"Bias vs Truth p_{T,Breit} (Relative Residual, Calibrated)",
       "Mean relative p_{T,Breit} residual vs truth p_{T,Breit} (calibrated)"},
      {"#DeltaR_{Breit} vs Truth p_{T,Breit} (Median)",
       "Median #DeltaR_{Breit} vs truth p_{T,Breit}"},
      {"Response vs Truth p_{T,Breit}", "Response vs truth p_{T,Breit}"},
      {"Response vs Truth p_{T,Breit} (Raw)",
       "Raw response vs truth p_{T,Breit}"},
      {"Response vs Truth p_{T,Breit} (Calibrated)",
       "Response vs truth p_{T,Breit} (calibrated)"},
      {"Response vs Truth p_{T,Breit} (Rho-corrected)",
       "Response vs truth p_{T,Breit} (rho-corrected)"},
      {"Tail rate |Relative residual|>0.2 vs Truth p_{T,Breit} (Calibrated)",
       "Tail rate |relative residual|>0.2 vs truth p_{T,Breit} (calibrated)"},
      {"|Relative residual| p95 vs Truth p_{T,Breit} (Calibrated)",
       "|relative residual| p95 vs truth p_{T,Breit} (calibrated)"},
      {"|Relative residual| p99 vs Truth p_{T,Breit} (Calibrated)",
       "|relative residual| p99 vs truth p_{T,Breit} (calibrated)"},
      {"Reco anti-k_{t} jets per event (Breit)", "Reco jets per event"},
      {"Reco Centauro jets per event (Breit)", "Reco jets per event"},
      {"Reco anti-k_{t} jets per event (Breit, baseline selection)", "Reco jets per event"},
      {"Reco Centauro jets per event (Breit, baseline selection)", "Reco jets per event"},
      {"Anti-k_{t} constituent #eta_{Breit} (truth)", "Truth constituent #eta_{Breit}"},
      {"Centauro constituent #eta_{Breit} (truth)", "Truth constituent #eta_{Breit}"},
      {"Anti-k_{t} constituent #eta_{Breit} (reco)", "Reco constituent #eta_{Breit}"},
      {"Centauro constituent #eta_{Breit} (reco)", "Reco constituent #eta_{Breit}"},
      {"Anti-k_{t} constituent #eta_{Breit} (truth, normalized)", "Truth constituent #eta_{Breit} normalized"},
      {"Centauro constituent #eta_{Breit} (truth, normalized)", "Truth constituent #eta_{Breit} normalized"},
      {"Anti-k_{t} constituent #eta_{Breit} (reco, normalized)", "Reco constituent #eta_{Breit} normalized"},
      {"Centauro constituent #eta_{Breit} (reco, normalized)", "Reco constituent #eta_{Breit} normalized"},
  };
  for (const auto& [from, to] : map) {
    if (title == from)
      return to;
  }

  title = replace_all(std::move(title), "vs Truth", "vs truth");
  return title;
}

enum class Stage {
  Raw,
  Corr,
  Calib,
};

Stage stage_from_base(const std::string& base) {
  if (base.find("calib") != std::string::npos) return Stage::Calib;
  if (base.find("corr") != std::string::npos) return Stage::Corr;
  return Stage::Raw;
}

TDirectory* ensure_dir(TDirectory* root, const std::string& rel_path) {
  if (!root) return nullptr;
  if (rel_path.empty()) return root;

  TDirectory* cur = root;
  std::size_t start = 0;
  while (start < rel_path.size()) {
    const auto end = rel_path.find('/', start);
    const std::string part = rel_path.substr(
        start, end == std::string::npos ? rel_path.size() - start : end - start);
    if (!part.empty()) {
      TDirectory* next = cur->GetDirectory(part.c_str());
      if (!next) next = cur->mkdir(part.c_str());
      cur = next;
    }
    if (end == std::string::npos) break;
    start = end + 1;
  }
  return cur;
}

bool same_binning_1d(const TH1& a, const TH1& b) {
  if (a.GetNbinsX() != b.GetNbinsX()) return false;
  const auto* ax = a.GetXaxis();
  const auto* bx = b.GetXaxis();
  if (!ax || !bx) return false;
  for (int bin = 1; bin <= a.GetNbinsX(); ++bin) {
    const double a_lo = ax->GetBinLowEdge(bin);
    const double a_hi = ax->GetBinUpEdge(bin);
    const double b_lo = bx->GetBinLowEdge(bin);
    const double b_hi = bx->GetBinUpEdge(bin);
    if (a_lo != b_lo || a_hi != b_hi) return false;
  }
  return true;
}

std::string base_name_from_alg(const std::string& name, std::string* alg_out) {
  const std::string a = "_antikt";
  const std::string c = "_centauro";
  const std::string a_label = "_anti-k_{t}";
  const std::string a_label_caps = "_anti-k_{T}";
  const std::string c_label = "_Centauro";
  const std::string a_mid = "_antikt_";
  const std::string c_mid = "_centauro_";
  const std::string a_mid_label = "_anti-k_{t}_";
  const std::string a_mid_label_caps = "_anti-k_{T}_";
  const std::string c_mid_label = "_Centauro_";
  const std::string a_beg = "antikt_";
  const std::string c_beg = "centauro_";
  const std::string a_beg_label = "anti-k_{t}_";
  const std::string a_beg_label_caps = "anti-k_{T}_";
  const std::string c_beg_label = "Centauro_";
  if (name.size() > a.size() && name.compare(name.size() - a.size(), a.size(), a) == 0) {
    if (alg_out) *alg_out = "antikt";
    return name.substr(0, name.size() - a.size());
  }
  if (name.size() > c.size() && name.compare(name.size() - c.size(), c.size(), c) == 0) {
    if (alg_out) *alg_out = "centauro";
    return name.substr(0, name.size() - c.size());
  }
  if (name.size() > a_label.size() &&
      name.compare(name.size() - a_label.size(), a_label.size(), a_label) == 0) {
    if (alg_out) *alg_out = "antikt";
    return name.substr(0, name.size() - a_label.size());
  }
  if (name.size() > a_label_caps.size() &&
      name.compare(name.size() - a_label_caps.size(), a_label_caps.size(), a_label_caps) == 0) {
    if (alg_out) *alg_out = "antikt";
    return name.substr(0, name.size() - a_label_caps.size());
  }
  if (name.size() > c_label.size() &&
      name.compare(name.size() - c_label.size(), c_label.size(), c_label) == 0) {
    if (alg_out) *alg_out = "centauro";
    return name.substr(0, name.size() - c_label.size());
  }
  if (name.size() > a_beg.size() && name.compare(0, a_beg.size(), a_beg) == 0) {
    if (alg_out) *alg_out = "antikt";
    return name.substr(a_beg.size());
  }
  if (name.size() > c_beg.size() && name.compare(0, c_beg.size(), c_beg) == 0) {
    if (alg_out) *alg_out = "centauro";
    return name.substr(c_beg.size());
  }
  if (name.size() > a_beg_label.size() && name.compare(0, a_beg_label.size(), a_beg_label) == 0) {
    if (alg_out) *alg_out = "antikt";
    return name.substr(a_beg_label.size());
  }
  if (name.size() > a_beg_label_caps.size() &&
      name.compare(0, a_beg_label_caps.size(), a_beg_label_caps) == 0) {
    if (alg_out) *alg_out = "antikt";
    return name.substr(a_beg_label_caps.size());
  }
  if (name.size() > c_beg_label.size() && name.compare(0, c_beg_label.size(), c_beg_label) == 0) {
    if (alg_out) *alg_out = "centauro";
    return name.substr(c_beg_label.size());
  }
  const bool has_a_mid = name.find(a_mid) != std::string::npos ||
                         name.find(a_mid_label) != std::string::npos ||
                         name.find(a_mid_label_caps) != std::string::npos;
  const bool has_c_mid = name.find(c_mid) != std::string::npos ||
                         name.find(c_mid_label) != std::string::npos;
  if (has_a_mid == has_c_mid)
    return {};
  const std::string& tag =
      name.find(a_mid) != std::string::npos
          ? a_mid
          : (name.find(a_mid_label) != std::string::npos
                 ? a_mid_label
                 : (name.find(a_mid_label_caps) != std::string::npos
                        ? a_mid_label_caps
                        : (name.find(c_mid) != std::string::npos ? c_mid : c_mid_label)));
  if (alg_out) *alg_out = has_a_mid ? "antikt" : "centauro";
  std::string base = name;
  const std::size_t pos = base.find(tag);
  if (pos == std::string::npos)
    return {};
  base.erase(pos, tag.size());
  return base;
}

bool strip_variant_suffix(std::string& base, std::string& variant) {
  const std::string raw = "_raw";
  const std::string groomed = "_groomed";
  if (base.size() > raw.size() &&
      base.compare(base.size() - raw.size(), raw.size(), raw) == 0) {
    variant = "raw";
    base.erase(base.size() - raw.size());
    return true;
  }
  if (base.size() > groomed.size() &&
      base.compare(base.size() - groomed.size(), groomed.size(), groomed) == 0) {
    variant = "groomed";
    base.erase(base.size() - groomed.size());
    return true;
  }
  return false;
}

bool extract_variant_key(const std::string& name, std::string& base,
                         std::string& alg, std::string& variant) {
  base = base_name_from_alg(name, &alg);
  if (base.empty())
    return false;
  return strip_variant_suffix(base, variant);
}

std::string frame_label_from_base(const std::string& base) {
  if (base.find("breit") != std::string::npos ||
      base.find("pT_B") != std::string::npos) {
    return "p_{T,B}";
  }
  return "p_{T,lab}";
}

TH1* pick_hist_from_canvas(TCanvas& c) {
  TList* prims = c.GetListOfPrimitives();
  if (!prims)
    return nullptr;
  for (TObject* obj : *prims) {
    if (!obj)
      continue;
    if (obj->InheritsFrom(TH1::Class()))
      return static_cast<TH1*>(obj);
  }
  return nullptr;
}

struct AxisTitles {
  std::string title;
  std::string x;
  std::string y;
};

using AxisTitleMap = std::map<std::string, AxisTitles>;

void apply_axis_titles(const AxisTitles& t, TH1& h) {
  if (!t.x.empty())
    h.GetXaxis()->SetTitle(t.x.c_str());
  if (!t.y.empty())
    h.GetYaxis()->SetTitle(t.y.c_str());
}

void write_overlay(TDirectory& out_dir, const std::string& base,
                   std::unique_ptr<TH1> h_a, std::unique_ptr<TH1> h_c,
                   const AxisTitleMap* axis_titles, bool no_errors) {
  if (!h_a || !h_c) return;
  if (!same_binning_1d(*h_a, *h_c)) return;

  std::string title =
      normalize_title(h_a->GetTitle() ? h_a->GetTitle() : base);
  if (title.empty() || title == h_a->GetName() || title == h_c->GetName() ||
      title.find("_antikt") != std::string::npos ||
      title.find("_centauro") != std::string::npos) {
    title = "Overlay " + base;
  }
  if (base.find("mean_resp_raw_pT") != std::string::npos &&
      title == "Response vs truth p_{T,Breit}") {
    title = "Raw response vs truth p_{T,Breit}";
  }
  const std::string canvas_name = "c_overlay_" + base;

  style_hist(*h_a, kBlue + 1, 20);
  style_hist(*h_c, kRed + 1, 21);
  h_a->SetStats(false);
  h_c->SetStats(false);
  h_a->SetTitle(title.c_str());
  h_c->SetTitle(title.c_str());
  if (axis_titles) {
    auto it = axis_titles->find(base);
    if (it != axis_titles->end()) {
      if (std::string(h_a->GetXaxis()->GetTitle()).empty() ||
          std::string(h_a->GetYaxis()->GetTitle()).empty()) {
        apply_axis_titles(it->second, *h_a);
        apply_axis_titles(it->second, *h_c);
      }
    }
  }
  if (base.find("mean_resp_raw_pT") != std::string::npos) {
    title = "Raw response vs truth p_{T,Breit}";
    h_a->SetTitle(title.c_str());
    h_c->SetTitle(title.c_str());
  }

  const auto y_range = hist_range(h_a.get(), h_c.get(), -1.0, 1.0);
  double y_min = y_range.min;
  double y_max = y_range.max;
  if (y_max - y_min < 1e-3) {
    y_min = 0.95;
    y_max = 1.05;
  }
  out_dir.cd();
  TCanvas c(canvas_name.c_str(), title.c_str(), 900, 600);
  const double x_min = h_a->GetXaxis()->GetXmin();
  const double x_max = h_a->GetXaxis()->GetXmax();
  TH1D frame(("axes_overlay_" + base).c_str(), title.c_str(), 100, x_min, x_max);
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.SetMinimum(y_min);
  frame.SetMaximum(y_max);
  frame.GetXaxis()->SetTitle(h_a->GetXaxis()->GetTitle());
  frame.GetYaxis()->SetTitle(h_a->GetYaxis()->GetTitle());
  frame.Draw("HIST");
  h_a->SetLineColor(kBlue + 1);
  h_c->SetLineColor(kRed + 1);
  h_a->SetLineWidth(2);
  h_c->SetLineWidth(2);
  h_a->Draw("HIST SAME");
  h_c->Draw("HIST SAME");
  std::unique_ptr<TGraphErrors> gr_a;
  std::unique_ptr<TGraphErrors> gr_c;
  if (!no_errors) {
    bool has_errors = false;
    for (int bin = 1; bin <= h_a->GetNbinsX(); ++bin) {
      if (h_a->GetBinError(bin) > 0.0 || h_c->GetBinError(bin) > 0.0) {
        has_errors = true;
        break;
      }
    }
    if (has_errors) {
      const double bw = h_a->GetNbinsX() > 0 ? h_a->GetXaxis()->GetBinWidth(1) : 0.0;
      const double dx = 0.15 * bw;
      const double x_err_frac = 0.0;
      gr_a = std::make_unique<TGraphErrors>(build_graph(*h_a, -dx, x_err_frac));
      gr_c = std::make_unique<TGraphErrors>(build_graph(*h_c, dx, x_err_frac));
      style_graph(*gr_a, kBlue + 1, 20);
      style_graph(*gr_c, kRed + 1, 21);
      gr_a->SetMarkerSize(0.9);
      gr_c->SetMarkerSize(0.9);
      gr_a->Draw("E1P SAME");
      gr_c->Draw("E1P SAME");
    }
  }

  TLegend leg(0.72, 0.74, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.035);
  leg.SetEntrySeparation(0.12);
  leg.AddEntry(h_a.get(), "anti-k_{t}", "L");
  leg.AddEntry(h_c.get(), "Centauro", "L");
  leg.Draw();
  c.Write();
}

void write_overlay_variant(TDirectory& out_dir, const std::string& base,
                           const std::string& label_a,
                           const std::string& label_b,
                           const std::string& title_override,
                           std::unique_ptr<TH1> h_a,
                           std::unique_ptr<TH1> h_b,
                           const AxisTitleMap* axis_titles, bool no_errors) {
  if (!h_a || !h_b) return;
  if (!same_binning_1d(*h_a, *h_b)) return;

  std::string title = title_override;
  if (title.empty()) {
    title = normalize_title(h_a->GetTitle() ? h_a->GetTitle() : base);
    if (title.empty() || title == h_a->GetName() || title == h_b->GetName()) {
      title = "Overlay " + base;
    }
  }
  const std::string canvas_name = "c_overlay_" + base;

  style_hist(*h_a, kBlue + 1, 20);
  style_hist(*h_b, kRed + 1, 21);
  h_a->SetStats(false);
  h_b->SetStats(false);
  h_a->SetTitle(title.c_str());
  h_b->SetTitle(title.c_str());
  if (axis_titles) {
    auto it = axis_titles->find(base);
    if (it != axis_titles->end()) {
      apply_axis_titles(it->second, *h_a);
      apply_axis_titles(it->second, *h_b);
    }
  }

  const auto y_range = hist_range(h_a.get(), h_b.get(), -1.0, 1.0);
  double y_min = y_range.min;
  double y_max = y_range.max;
  if (y_max - y_min < 1e-3) {
    y_min = 0.95;
    y_max = 1.05;
  }
  out_dir.cd();
  TCanvas c(canvas_name.c_str(), title.c_str(), 900, 600);
  const double x_min = h_a->GetXaxis()->GetXmin();
  const double x_max = h_a->GetXaxis()->GetXmax();
  TH1D frame(("axes_overlay_" + base).c_str(), title.c_str(), 100, x_min, x_max);
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.SetMinimum(y_min);
  frame.SetMaximum(y_max);
  frame.GetXaxis()->SetTitle(h_a->GetXaxis()->GetTitle());
  frame.GetYaxis()->SetTitle(h_a->GetYaxis()->GetTitle());
  frame.Draw("HIST");
  h_a->SetLineColor(kBlue + 1);
  h_b->SetLineColor(kRed + 1);
  h_a->SetLineWidth(2);
  h_b->SetLineWidth(2);
  h_a->Draw("HIST SAME");
  h_b->Draw("HIST SAME");
  std::unique_ptr<TGraphErrors> gr_a;
  std::unique_ptr<TGraphErrors> gr_b;
  if (!no_errors) {
    bool has_errors = false;
    for (int bin = 1; bin <= h_a->GetNbinsX(); ++bin) {
      if (h_a->GetBinError(bin) > 0.0 || h_b->GetBinError(bin) > 0.0) {
        has_errors = true;
        break;
      }
    }
    if (has_errors) {
      const double bw = h_a->GetNbinsX() > 0 ? h_a->GetXaxis()->GetBinWidth(1) : 0.0;
      const double dx = 0.15 * bw;
      const double x_err_frac = 0.0;
      gr_a = std::make_unique<TGraphErrors>(build_graph(*h_a, -dx, x_err_frac));
      gr_b = std::make_unique<TGraphErrors>(build_graph(*h_b, dx, x_err_frac));
      style_graph(*gr_a, kBlue + 1, 20);
      style_graph(*gr_b, kRed + 1, 21);
      gr_a->SetMarkerSize(0.9);
      gr_b->SetMarkerSize(0.9);
      gr_a->Draw("E1P SAME");
      gr_b->Draw("E1P SAME");
    }
  }

  TLegend leg(0.72, 0.74, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.035);
  leg.SetEntrySeparation(0.12);
  leg.AddEntry(h_a.get(), label_a.c_str(), "L");
  leg.AddEntry(h_b.get(), label_b.c_str(), "L");
  leg.Draw();
  c.Write();
}

void style_legend(TLegend& leg) {
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.035);
  leg.SetEntrySeparation(0.12);
  leg.SetX1NDC(0.72);
  leg.SetY1NDC(0.74);
  leg.SetX2NDC(0.88);
  leg.SetY2NDC(0.88);
}

void write_canvas_copy(TDirectory& out_dir, std::unique_ptr<TCanvas> c) {
  if (!c)
    return;
  c->SetTitle(normalize_title(c->GetTitle() ? c->GetTitle() : "").c_str());
  if (TList* prims = c->GetListOfPrimitives()) {
    for (TObject* obj : *prims) {
      if (!obj)
        continue;
      if (obj->InheritsFrom(TLegend::Class())) {
        style_legend(*static_cast<TLegend*>(obj));
        continue;
      }
      if (!obj->InheritsFrom(TH1::Class()))
        continue;
      auto* h = static_cast<TH1*>(obj);
      h->SetTitle(normalize_title(h->GetTitle() ? h->GetTitle() : "").c_str());
    }
  }
  out_dir.cd();
  c->Write();
}

void harvest_axis_titles(TDirectory& in_dir, AxisTitleMap& axis_titles,
                         bool overlay_variants) {
  TIter next(in_dir.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next())) {
    const char* name = key->GetName();
    if (!name) continue;
    std::unique_ptr<TObject> obj(key->ReadObj());
    if (!obj) continue;
    if (obj->InheritsFrom(TDirectory::Class())) {
      auto* sub_in = dynamic_cast<TDirectory*>(obj.get());
      if (sub_in)
        harvest_axis_titles(*sub_in, axis_titles, overlay_variants);
      continue;
    }
    if (!obj->InheritsFrom(TCanvas::Class()))
      continue;
    auto* src = dynamic_cast<TCanvas*>(obj.get());
    if (!src)
      continue;
    if (TList* prims = src->GetListOfPrimitives()) {
      for (TObject* prim : *prims) {
        if (!prim || !prim->InheritsFrom(TH1::Class()))
          continue;
        auto* h = static_cast<TH1*>(prim);
        std::string alg;
        std::string base;
        if (overlay_variants) {
          std::string variant;
          if (!extract_variant_key(h->GetName(), base, alg, variant))
            continue;
        } else {
          base = base_name_from_alg(h->GetName(), &alg);
          if (base.empty())
            continue;
        }
        AxisTitles& t = axis_titles[base];
        if (t.title.empty() && h->GetTitle())
          t.title = normalize_title(h->GetTitle());
        const char* xt = h->GetXaxis() ? h->GetXaxis()->GetTitle() : "";
        const char* yt = h->GetYaxis() ? h->GetYaxis()->GetTitle() : "";
        if (t.x.empty() && xt)
          t.x = xt;
        if (t.y.empty() && yt)
          t.y = yt;
      }
    }
  }
}

void process_dir(TDirectory& in_dir, TDirectory& out_dir,
                 const AxisTitleMap* axis_titles, bool no_errors,
                 bool overlay_variants) {
  struct Pair {
    std::unique_ptr<TH1> antikt;
    std::unique_ptr<TH1> centauro;
  };

  std::map<std::string, Pair> pairs;

  struct VariantPair {
    std::unique_ptr<TH1> raw;
    std::unique_ptr<TH1> groomed;
    std::string alg;
  };

  std::map<std::string, VariantPair> variant_pairs;

  TIter next(in_dir.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next())) {
    const char* name = key->GetName();
    if (!name) continue;

    std::unique_ptr<TObject> obj(key->ReadObj());
    if (!obj) continue;

    if (obj->InheritsFrom(TDirectory::Class())) {
      auto* sub_in = dynamic_cast<TDirectory*>(obj.get());
      if (!sub_in) continue;
      auto* sub_out = out_dir.mkdir(name);
      if (!sub_out) continue;
      process_dir(*sub_in, *sub_out, axis_titles, no_errors, overlay_variants);
      continue;
    }

    if (obj->InheritsFrom(TCanvas::Class())) {
      auto* src = dynamic_cast<TCanvas*>(obj.get());
      if (!src)
        continue;
      bool added = false;
      if (TList* prims = src->GetListOfPrimitives()) {
        for (TObject* prim : *prims) {
          if (!prim || !prim->InheritsFrom(TH1::Class()))
            continue;
          auto* h_src = static_cast<TH1*>(prim);
          std::string alg;
          std::string base;
          std::string variant;
          if (overlay_variants) {
            if (!extract_variant_key(h_src->GetName(), base, alg, variant))
              continue;
          } else {
            base = base_name_from_alg(h_src->GetName(), &alg);
            if (base.empty())
              continue;
          }
          auto h = clone_hist(*h_src);
          if (!h)
            continue;
          if (overlay_variants) {
            const std::string key = base + "_" + alg;
            auto& slot = variant_pairs[key];
            slot.alg = alg;
            if (variant == "raw") {
              slot.raw = std::move(h);
              added = true;
            } else if (variant == "groomed") {
              slot.groomed = std::move(h);
              added = true;
            }
          } else {
            auto& slot = pairs[base];
            if (alg == "antikt") {
              slot.antikt = std::move(h);
              added = true;
            } else if (alg == "centauro") {
              slot.centauro = std::move(h);
              added = true;
            }
          }
        }
      }
      if (!added) {
        write_canvas_copy(out_dir, clone_canvas(*src));
      }
      continue;
    }

    if (!obj->InheritsFrom(TH1::Class())) continue;
    auto* src = dynamic_cast<TH1*>(obj.get());
    if (!src) continue;
    if (src->GetDimension() != 1) continue;

    std::string alg;
    std::string base;
    std::string variant;
    if (overlay_variants) {
      if (!extract_variant_key(name, base, alg, variant))
        continue;
    } else {
      base = base_name_from_alg(name, &alg);
      if (base.empty()) continue;
    }

    auto h = clone_hist(*src);
    if (!h) continue;

    if (overlay_variants) {
      const std::string key = base + "_" + alg;
      auto& slot = variant_pairs[key];
      slot.alg = alg;
      if (variant == "raw") {
        slot.raw = std::move(h);
      } else if (variant == "groomed") {
        slot.groomed = std::move(h);
      }
    } else {
      auto& slot = pairs[base];
      if (alg == "antikt") {
        slot.antikt = std::move(h);
      } else if (alg == "centauro") {
        slot.centauro = std::move(h);
      }
    }
  }

  if (overlay_variants) {
    for (auto& [key, pair] : variant_pairs) {
      if (!pair.raw || !pair.groomed) continue;
      const std::string alg_label =
          (pair.alg == "antikt") ? "anti-k_{t}" : "Centauro";
      const std::string frame_label = frame_label_from_base(key);
      const std::string title =
          alg_label + " N_{const} vs " + frame_label;
      write_overlay_variant(out_dir, key, "Ungroomed", "Soft Drop", title,
                            std::move(pair.raw), std::move(pair.groomed),
                            axis_titles, no_errors);
    }
  } else {
    for (auto& [base, pair] : pairs) {
      if (!pair.antikt || !pair.centauro) continue;
      write_overlay(out_dir, base, std::move(pair.antikt), std::move(pair.centauro),
                    axis_titles, no_errors);
    }
  }
}

void process_hists_tree(TDirectory& in_dir, TDirectory* raw_root, TDirectory* corr_root,
                        TDirectory* calib_root, const std::string& rel_path,
                        const AxisTitleMap* axis_titles, bool no_errors,
                        bool overlay_variants) {
  struct Pair {
    std::unique_ptr<TH1> antikt;
    std::unique_ptr<TH1> centauro;
  };

  std::map<std::string, Pair> pairs;

  struct VariantPair {
    std::unique_ptr<TH1> raw;
    std::unique_ptr<TH1> groomed;
    std::string alg;
  };

  std::map<std::string, VariantPair> variant_pairs;

  TIter next(in_dir.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next())) {
    const char* name = key->GetName();
    if (!name) continue;

    std::unique_ptr<TObject> obj(key->ReadObj());
    if (!obj) continue;

    if (obj->InheritsFrom(TDirectory::Class())) {
      auto* sub = dynamic_cast<TDirectory*>(obj.get());
      if (!sub) continue;
      const std::string next_rel = rel_path.empty() ? std::string(name)
                                                    : (rel_path + "/" + std::string(name));
      process_hists_tree(*sub, raw_root, corr_root, calib_root, next_rel,
                         axis_titles, no_errors, overlay_variants);
      continue;
    }

    if (!obj->InheritsFrom(TH1::Class())) continue;
    auto* src = dynamic_cast<TH1*>(obj.get());
    if (!src) continue;
    if (src->GetDimension() != 1) continue;

    std::string alg;
    std::string base;
    std::string variant;
    if (overlay_variants) {
      if (!extract_variant_key(name, base, alg, variant))
        continue;
    } else {
      base = base_name_from_alg(name, &alg);
      if (base.empty()) continue;
    }

    auto h = clone_hist(*src);
    if (!h) continue;

    if (overlay_variants) {
      const std::string key = base + "_" + alg;
      auto& slot = variant_pairs[key];
      slot.alg = alg;
      if (variant == "raw") {
        slot.raw = std::move(h);
      } else if (variant == "groomed") {
        slot.groomed = std::move(h);
      }
    } else {
      auto& slot = pairs[base];
      if (alg == "antikt") {
        slot.antikt = std::move(h);
      } else if (alg == "centauro") {
        slot.centauro = std::move(h);
      }
    }
  }

  if (overlay_variants) {
    for (auto& [key, pair] : variant_pairs) {
      if (!pair.raw || !pair.groomed) continue;

      TDirectory* out_root = raw_root;
      const Stage stage = stage_from_base(key);
      if (stage == Stage::Corr) out_root = corr_root;
      if (stage == Stage::Calib) out_root = calib_root;
      if (!out_root) continue;

      TDirectory* out_dir = ensure_dir(out_root, rel_path);
      if (!out_dir) continue;
      const std::string alg_label =
          (pair.alg == "antikt") ? "anti-k_{t}" : "Centauro";
      const std::string frame_label = frame_label_from_base(key);
      const std::string title =
          alg_label + " N_{const} vs " + frame_label;
      write_overlay_variant(*out_dir, key, "Ungroomed", "Soft Drop", title,
                            std::move(pair.raw), std::move(pair.groomed),
                            axis_titles, no_errors);
    }
  } else {
    for (auto& [base, pair] : pairs) {
      if (!pair.antikt || !pair.centauro) continue;

      TDirectory* out_root = raw_root;
      const Stage stage = stage_from_base(base);
      if (stage == Stage::Corr) out_root = corr_root;
      if (stage == Stage::Calib) out_root = calib_root;
      if (!out_root) continue;

      TDirectory* out_dir = ensure_dir(out_root, rel_path);
      if (!out_dir) continue;
      write_overlay(*out_dir, base, std::move(pair.antikt), std::move(pair.centauro),
                    axis_titles, no_errors);
    }
  }
}

void process_stage_dir(TDirectory& stage_in, TDirectory* stage_out_root,
                       const AxisTitleMap* axis_titles, bool no_errors,
                       bool overlay_variants) {
  if (!stage_out_root)
    return;
  process_dir(stage_in, *stage_out_root, axis_titles, no_errors, overlay_variants);
}

}  // namespace

int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);

  TFile fin(args.input.c_str(), "READ");
  if (!fin.IsOpen()) {
    std::cerr << "[overlay] Failed to open input " << args.input << "\n";
    return 1;
  }

  ensure_parent(args.output);
  TFile fout(args.output.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[overlay] Failed to open output " << args.output << "\n";
    return 1;
  }

  gStyle->SetOptStat(0);

  TDirectory* in_raw = fin.GetDirectory("raw");
  TDirectory* in_corr = fin.GetDirectory("corr");
  TDirectory* in_calib = fin.GetDirectory("calib");
  TDirectory* hists = fin.GetDirectory("hists");
  const bool has_stage_dirs = (in_raw || in_corr || in_calib);
  TDirectory* raw_root = nullptr;
  TDirectory* corr_root = nullptr;
  TDirectory* calib_root = nullptr;
  if (has_stage_dirs) {
    if (in_raw)
      raw_root = fout.mkdir("raw");
    if (in_corr)
      corr_root = fout.mkdir("corr");
    if (in_calib)
      calib_root = fout.mkdir("calib");
  } else {
    raw_root = &fout;
    corr_root = &fout;
    calib_root = &fout;
  }
  AxisTitleMap axis_titles;
  if (in_raw)
    harvest_axis_titles(*in_raw, axis_titles, args.overlay_variants);
  if (in_corr)
    harvest_axis_titles(*in_corr, axis_titles, args.overlay_variants);
  if (in_calib)
    harvest_axis_titles(*in_calib, axis_titles, args.overlay_variants);

  if (hists) {
    process_hists_tree(*hists, raw_root, corr_root, calib_root, "", &axis_titles,
                       args.no_errors, args.overlay_variants);
  } else if (in_raw || in_corr || in_calib) {
    if (in_raw)
      process_stage_dir(*in_raw, raw_root, &axis_titles, args.no_errors,
                        args.overlay_variants);
    if (in_corr)
      process_stage_dir(*in_corr, corr_root, &axis_titles, args.no_errors,
                        args.overlay_variants);
    if (in_calib)
      process_stage_dir(*in_calib, calib_root, &axis_titles, args.no_errors,
                        args.overlay_variants);
  } else {
    process_dir(fin, fout, &axis_titles, args.no_errors, args.overlay_variants);
  }

  fout.Close();
  std::cout << "\n[overlay] Wrote overlays to " << args.output << "\n";
  return 0;
}

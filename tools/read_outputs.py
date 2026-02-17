#!/usr/bin/env python3
import argparse
import math

import ROOT


def eff_sigma68(values):
    vals = [v for v in values if math.isfinite(v)]
    if not vals:
        return float("nan")
    vals.sort()
    n = len(vals)
    n68 = int(math.ceil(0.68 * n))
    if n68 < 2 or n68 > n:
        return float("nan")
    best = None
    for i in range(0, n - n68 + 1):
        width = vals[i + n68 - 1] - vals[i]
        if best is None or width < best:
            best = width
    if best is None or not math.isfinite(best):
        return float("nan")
    return 0.5 * best


def open_root(path, label):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Failed to open {label}: {path}")
    return f


def mean_hist(f, name):
    h = f.Get(name)
    if not h:
        raise RuntimeError(f"Missing histogram: {name}")
    return h.GetMean()


def print_header(title):
    print(f"\n== {title} ==")


def fmt(v, nd=4):
    if not math.isfinite(v):
        return "nan"
    return f"{v:.{nd}f}"


def main():
    parser = argparse.ArgumentParser(description="Read ROOT outputs for data_only.tex tables.")
    parser.add_argument("--calib", default="data/jets/matched/calib_merged_0-999.root")
    parser.add_argument("--performance", default="data/graphs/performance_0-999.root")
    parser.add_argument("--controls", default="data/graphs/controls_0-999.root")
    parser.add_argument("--frag-truth", default="data/graphs/frag_shapes_0-999_truth.root")
    parser.add_argument("--frag-reco", default="data/graphs/frag_shapes_0-999_reco.root")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)

    f_calib = open_root(args.calib, "calibration matches")
    f_perf = open_root(args.performance, "performance")
    f_ctrl = open_root(args.controls, "controls")
    f_frag_truth = open_root(args.frag_truth, "frag shapes truth")
    f_frag_reco = open_root(args.frag_reco, "frag shapes reco")

    # Integrated response/resolution (5-20 GeV, |eta|<=3)
    t = f_calib.Get("matches")
    if not t:
        raise RuntimeError("Missing matches tree in calibration file.")
    vals = {0: {"resp": [], "resid": [], "resp_calib": [], "resid_calib": []},
            1: {"resp": [], "resid": [], "resp_calib": [], "resid_calib": []}}
    for entry in t:
        if abs(entry.eta_truth) > 3:
            continue
        if entry.pT_truth < 5 or entry.pT_truth >= 20:
            continue
        if entry.alg < 0 or entry.alg > 1:
            continue
        if not math.isfinite(entry.pT_reco_raw):
            continue
        if not (math.isfinite(entry.phi_reco) and math.isfinite(entry.phi_truth)):
            continue
        if not (math.isfinite(entry.eta_reco) and math.isfinite(entry.eta_truth)):
            continue
        if entry.pT_truth <= 0:
            continue
        resp = entry.pT_reco_raw / entry.pT_truth
        resid = (entry.pT_reco_raw - entry.pT_truth) / entry.pT_truth
        resp_calib = entry.resp_calib
        if not math.isfinite(resp_calib):
            resp_calib = entry.pT_calib / entry.pT_truth if entry.pT_truth > 0 else float("nan")
        resid_calib = resp_calib - 1.0 if math.isfinite(resp_calib) else float("nan")
        vals[int(entry.alg)]["resp"].append(resp)
        vals[int(entry.alg)]["resid"].append(resid)
        vals[int(entry.alg)]["resp_calib"].append(resp_calib)
        vals[int(entry.alg)]["resid_calib"].append(resid_calib)

    print_header("Integrated Response/Resolution (5-20 GeV)")
    for alg in (0, 1):
        resp_vals = vals[alg]["resp"]
        resid_vals = vals[alg]["resid"]
        mean_resp = sum(resp_vals) / len(resp_vals) if resp_vals else float("nan")
        eff = eff_sigma68(resid_vals)
        label = "anti-k_t" if alg == 0 else "Centauro"
        print(f"{label}: mean_resp_raw={fmt(mean_resp)} effSigma68={fmt(eff)}")

    # Performance binned
    h_resp = {
        0: f_perf.Get("hists/mean_resp_raw_pT_antikt"),
        1: f_perf.Get("hists/mean_resp_raw_pT_centauro"),
    }
    h_bias = {
        0: f_perf.Get("hists/bias_rel_resid_pT_antikt"),
        1: f_perf.Get("hists/bias_rel_resid_pT_centauro"),
    }
    h_eff = {
        0: f_perf.Get("hists/effSigma68_rel_resid_pT_antikt"),
        1: f_perf.Get("hists/effSigma68_rel_resid_pT_centauro"),
    }
    h_med_dr = {
        0: f_perf.Get("hists/median_dR_pT_antikt"),
        1: f_perf.Get("hists/median_dR_pT_centauro"),
    }
    h_dphi = {
        0: f_perf.Get("hists/effSigma68_dphi_pT_antikt"),
        1: f_perf.Get("hists/effSigma68_dphi_pT_centauro"),
    }
    h_deta = {
        0: f_perf.Get("hists/effSigma68_deta_pT_antikt"),
        1: f_perf.Get("hists/effSigma68_deta_pT_centauro"),
    }

    print_header("Binned Response/Resolution (pT bins 1-4)")
    for b in range(1, 5):
        lo = h_resp[0].GetXaxis().GetBinLowEdge(b)
        hi = h_resp[0].GetXaxis().GetBinUpEdge(b)
        print(
            f"{int(lo)}-{int(hi)} "
            f"resp={fmt(h_resp[0].GetBinContent(b))}+-{fmt(h_resp[0].GetBinError(b))},"
            f"{fmt(h_resp[1].GetBinContent(b))}+-{fmt(h_resp[1].GetBinError(b))} "
            f"bias={fmt(h_bias[0].GetBinContent(b))}+-{fmt(h_bias[0].GetBinError(b))},"
            f"{fmt(h_bias[1].GetBinContent(b))}+-{fmt(h_bias[1].GetBinError(b))} "
            f"effSigma68={fmt(h_eff[0].GetBinContent(b))}+-{fmt(h_eff[0].GetBinError(b))},"
            f"{fmt(h_eff[1].GetBinContent(b))}+-{fmt(h_eff[1].GetBinError(b))}"
        )

    print_header("Binned Angular Metrics (pT bins 1-4)")
    for b in range(1, 5):
        lo = h_med_dr[0].GetXaxis().GetBinLowEdge(b)
        hi = h_med_dr[0].GetXaxis().GetBinUpEdge(b)
        print(
            f"{int(lo)}-{int(hi)} "
            f"median_dR={fmt(h_med_dr[0].GetBinContent(b))},{fmt(h_med_dr[1].GetBinContent(b))} "
            f"dphi={fmt(h_dphi[0].GetBinContent(b))}+-{fmt(h_dphi[0].GetBinError(b))},"
            f"{fmt(h_dphi[1].GetBinContent(b))}+-{fmt(h_dphi[1].GetBinError(b))} "
            f"deta={fmt(h_deta[0].GetBinContent(b))}+-{fmt(h_deta[0].GetBinError(b))},"
            f"{fmt(h_deta[1].GetBinContent(b))}+-{fmt(h_deta[1].GetBinError(b))}"
        )

    # Yield summary
    print_header("Yield Summary")
    truth_jets_event = {
        "antikt": mean_hist(f_ctrl, "hists/njets_event_truth_antikt"),
        "cent": mean_hist(f_ctrl, "hists/njets_event_truth_centauro"),
    }
    truth_pT_mean = {
        "antikt": mean_hist(f_ctrl, "hists/truth_jet_pT_antikt"),
        "cent": mean_hist(f_ctrl, "hists/truth_jet_pT_centauro"),
    }
    truth_eta_mean = {
        "antikt": mean_hist(f_ctrl, "hists/truth_jet_eta_breit_antikt"),
        "cent": mean_hist(f_ctrl, "hists/truth_jet_eta_breit_centauro"),
    }
    reco_jets_event = {
        "antikt": mean_hist(f_ctrl, "hists/njets_event_reco_antikt"),
        "cent": mean_hist(f_ctrl, "hists/njets_event_reco_centauro"),
    }
    reco_unmatched_event = {
        "antikt": mean_hist(f_ctrl, "hists/njets_event_reco_unmatched_antikt"),
        "cent": mean_hist(f_ctrl, "hists/njets_event_reco_unmatched_centauro"),
    }
    reco_pT_mean = {
        "antikt": mean_hist(f_ctrl, "hists/reco_jet_pT_antikt"),
        "cent": mean_hist(f_ctrl, "hists/reco_jet_pT_centauro"),
    }
    reco_eta_mean = {
        "antikt": mean_hist(f_ctrl, "hists/reco_jet_eta_antikt"),
        "cent": mean_hist(f_ctrl, "hists/reco_jet_eta_centauro"),
    }
    print(
        "truth_jets_event "
        f"{fmt(truth_jets_event['antikt'], 3)} {fmt(truth_jets_event['cent'], 3)}"
    )
    print(
        "truth_pT_mean "
        f"{fmt(truth_pT_mean['antikt'], 3)} {fmt(truth_pT_mean['cent'], 3)}"
    )
    print(
        "truth_eta_mean "
        f"{fmt(truth_eta_mean['antikt'], 3)} {fmt(truth_eta_mean['cent'], 3)}"
    )
    print(
        "reco_jets_event "
        f"{fmt(reco_jets_event['antikt'], 3)} {fmt(reco_jets_event['cent'], 3)}"
    )
    print(
        "reco_unmatched_event "
        f"{fmt(reco_unmatched_event['antikt'], 3)} {fmt(reco_unmatched_event['cent'], 3)}"
    )
    print(
        "reco_pT_mean "
        f"{fmt(reco_pT_mean['antikt'], 3)} {fmt(reco_pT_mean['cent'], 3)}"
    )
    print(
        "reco_eta_mean "
        f"{fmt(reco_eta_mean['antikt'], 3)} {fmt(reco_eta_mean['cent'], 3)}"
    )

    # Efficiency vs pT and eta
    print_header("Efficiency vs pT (bins 1-4)")
    h_eff_pT = {
        0: f_ctrl.Get("hists/eff_match_pT_antikt"),
        1: f_ctrl.Get("hists/eff_match_pT_centauro"),
    }
    for b in range(1, 5):
        lo = h_eff_pT[0].GetXaxis().GetBinLowEdge(b)
        hi = h_eff_pT[0].GetXaxis().GetBinUpEdge(b)
        print(
            f"{int(lo)}-{int(hi)} "
            f"{fmt(h_eff_pT[0].GetBinContent(b))}+-{fmt(h_eff_pT[0].GetBinError(b))},"
            f"{fmt(h_eff_pT[1].GetBinContent(b))}+-{fmt(h_eff_pT[1].GetBinError(b))}"
        )

    print_header("Efficiency vs eta (all bins)")
    h_eff_eta = {
        0: f_ctrl.Get("hists/eff_match_eta_antikt"),
        1: f_ctrl.Get("hists/eff_match_eta_centauro"),
    }
    nb = h_eff_eta[0].GetNbinsX()
    for b in range(1, nb + 1):
        lo = h_eff_eta[0].GetXaxis().GetBinLowEdge(b)
        hi = h_eff_eta[0].GetXaxis().GetBinUpEdge(b)
        print(
            f"[{lo:.1f},{hi:.1f}] "
            f"{fmt(h_eff_eta[0].GetBinContent(b))}+-{fmt(h_eff_eta[0].GetBinError(b))},"
            f"{fmt(h_eff_eta[1].GetBinContent(b))}+-{fmt(h_eff_eta[1].GetBinError(b))}"
        )

    # Fragmentation summary
    def frag_summary(f, tag):
        h_n = f.Get(f"jet_obs/n_const_{tag}")
        h_dz = f.Get(f"inclusive/all/Dz_{tag}")
        h_psi = f.Get(f"inclusive/all/psi_{tag}")
        if not h_n or not h_dz or not h_psi:
            raise RuntimeError(f"Missing frag shapes objects for {tag}")
        n_const = h_n.GetMean()
        bmin = h_dz.FindBin(0.010001)
        bmax = h_dz.FindBin(0.999999)
        dz_int = h_dz.Integral(bmin, bmax, "width")
        bpsi = h_psi.FindBin(0.6 - 1e-6)
        psi_val = h_psi.GetBinContent(bpsi)
        return n_const, dz_int, psi_val

    print_header("Fragmentation Summary")
    for label, f in [("truth", f_frag_truth), ("reco", f_frag_reco)]:
        antikt = frag_summary(f, "antikt")
        cent = frag_summary(f, "centauro")
        print(
            f"{label} anti-k_t {fmt(antikt[0])} {fmt(antikt[1])} {fmt(antikt[2])}"
        )
        print(
            f"{label} Centauro {fmt(cent[0])} {fmt(cent[1])} {fmt(cent[2])}"
        )

if __name__ == "__main__":
    main()

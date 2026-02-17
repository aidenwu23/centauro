#include "jet_tools/include/calibration.h"
#include "jet_tools/include/math_helpers.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

namespace jet_tools {
namespace {

void coarsen_bins(std::vector<std::vector<double>> &bin_resp,
                  std::vector<double> &scales, std::size_t N_min) {
  const std::size_t n = bin_resp.size();
  scales.assign(n, 1.0);
  if (n == 0)
    return;

  std::vector<bool> active(n, true);
  std::size_t total_samples = 0;
  for (const auto &v : bin_resp)
    total_samples += v.size();
  if (total_samples == 0)
    return;

  // Merge sparse neighboring bins until each active bin passes the N_min threshold.
  auto count = [&](std::size_t i) {
    return active[i] ? bin_resp[i].size() : 0UL;
  };
  auto nearest_active_left = [&](std::size_t index) {
    for (std::size_t i = index; i-- > 0;) {
      if (active[i])
        return std::optional<std::size_t>(i);
    }
    return std::optional<std::size_t>();
  };
  auto nearest_active_right = [&](std::size_t index) {
    for (std::size_t i = index + 1; i < n; ++i) {
      if (active[i])
        return std::optional<std::size_t>(i);
    }
    return std::optional<std::size_t>();
  };

  while (true) {
    std::size_t index = n;
    for (std::size_t i = 0; i < n; ++i) {
      if (active[i] && count(i) < N_min) {
        index = i;
        break;
      }
    }
    if (index == n)
      break;
    auto l = nearest_active_left(index);
    auto r = nearest_active_right(index);
    if (!l && !r)
      break;
    std::size_t nb = n;
    const std::size_t lc = l ? count(*l) : 0;
    const std::size_t rc = r ? count(*r) : 0;
    if (!l) {
      nb = *r;
    } else if (!r) {
      nb = *l;
    } else {
      nb = (rc > lc) ? *r : *l;
    }
    if (nb == n)
      break;
    bin_resp[nb].insert(bin_resp[nb].end(), bin_resp[index].begin(),
                        bin_resp[index].end());
    bin_resp[index].clear();
    active[index] = false;
  }

  std::vector<double> med(n, std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 0; i < n; ++i) {
    if (!active[i] || bin_resp[i].empty())
      continue;
    med[i] = median_value(bin_resp[i]);
  }

  for (std::size_t i = 0; i < n; ++i) {
    double best_val = 1.0;
    std::size_t best_dist = std::numeric_limits<std::size_t>::max();
    for (std::size_t j = 0; j < n; ++j) {
      if (!active[j])
        continue;
      if (!std::isfinite(med[j]) || med[j] <= 0.0)
        continue;
      const std::size_t dist = (j > i) ? (j - i) : (i - j);
      if (dist < best_dist) {
        best_dist = dist;
        best_val = 1.0 / med[j];
        if (best_dist == 0)
          break;
      }
    }
    scales[i] = best_val;
  }
}

std::vector<double> default_y_edges() {
  return {0.05, 0.2, 0.35, 0.5, 0.7};
}

std::vector<double> coarsen_edges(const std::vector<double> &edges,
                                  std::size_t max_bins) {
  if (edges.size() < 2)
    return edges;
  const std::size_t bins = edges.size() - 1;
  if (bins <= max_bins)
    return edges;
  const std::size_t step =
      static_cast<std::size_t>(std::ceil(static_cast<double>(bins) /
                                         static_cast<double>(max_bins)));
  std::vector<double> out;
  out.reserve(max_bins + 1);
  for (std::size_t i = 0; i < bins; i += step) {
    out.push_back(edges[i]);
  }
  if (out.empty() || out.back() != edges.back())
    out.push_back(edges.back());
  return out;
}

std::size_t bin_index_clamped(const std::vector<double> &edges, double x) {
  if (edges.size() < 2)
    return 0;
  const std::size_t n = edges.size() - 1;
  if (x < edges.front())
    return 0;
  if (x >= edges.back())
    return n - 1;
  for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
    if (x >= edges[b] && x < edges[b + 1])
      return b;
  }
  return n - 1;
}

}  // namespace

double interpolate_bins(const std::vector<double> &edges,
                        const std::vector<double> &vals, double x) {
  if (edges.size() < 2 || vals.empty())
    return 1.0;
  for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
    if (x >= edges[b] && x < edges[b + 1]) {
      return vals[b];
    }
  }
  if (x < edges.front())
    return vals.front();
  return vals.back();
}

CalibMaps derive_calibration_maps(const CalibConfig &cfg,
                                  const std::vector<MatchRecord> &records) {
  CalibMaps maps;
  maps.eta_edges = cfg.eta_bins;
  maps.pT_edges = cfg.pT_bins;
  maps.pT_dis_edges = coarsen_edges(cfg.pT_bins, 6);
  maps.y_edges = default_y_edges();
  const std::size_t N_min = cfg.N_min;
  double pT_min = 5.0;
  double pT_max = 60.0;
  if (cfg.pT_bins.size() >= 2 && cfg.pT_bins.front() < cfg.pT_bins.back()) {
    pT_min = cfg.pT_bins.front();
    pT_max = cfg.pT_bins.back();
  }

  for (int alg = 0; alg < 2; ++alg) {
    maps.crel[alg].resize(maps.eta_edges.size() - 1, 1.0);
    std::vector<std::vector<double>> bin_resp(maps.eta_edges.size() - 1);
    for (const auto &r : records) {
      if (!r.is_train)
        continue;
      if (r.alg != alg)
        continue;
      if (r.pT_truth < pT_min)
        continue;
      if (r.pT_reco_raw < pT_min || r.pT_reco_raw > pT_max)
        continue;
      const double eta = r.eta_reco;
      for (std::size_t b = 0; b + 1 < maps.eta_edges.size(); ++b) {
        if (eta >= maps.eta_edges[b] && eta < maps.eta_edges[b + 1]) {
          const double resp = r.pT_reco_raw / r.pT_truth;
          if (std::isfinite(resp) && resp > 0.0)
            bin_resp[b].push_back(resp);
          break;
        }
      }
    }
    coarsen_bins(bin_resp, maps.crel[alg], N_min);
  }

  for (int alg = 0; alg < 2; ++alg) {
    maps.cabs[alg].resize(maps.pT_edges.size() - 1, 1.0);
    std::vector<std::vector<double>> bin_resp(maps.pT_edges.size() - 1);
    for (const auto &r : records) {
      if (!r.is_train)
        continue;
      if (r.alg != alg)
        continue;
      if (r.pT_truth < pT_min)
        continue;
      if (std::abs(r.eta_reco) >= cfg.max_abs_eta)
        continue;

      double crel = 1.0;
      for (std::size_t b = 0; b + 1 < maps.eta_edges.size(); ++b) {
        if (r.eta_reco >= maps.eta_edges[b] && r.eta_reco < maps.eta_edges[b + 1]) {
          crel = maps.crel[alg][b];
          break;
        }
      }
      const double pT_for_abs = r.pT_reco_raw * crel;
      const double resp = pT_for_abs / r.pT_truth;

      for (std::size_t b = 0; b + 1 < maps.pT_edges.size(); ++b) {
        if (pT_for_abs >= maps.pT_edges[b] && pT_for_abs < maps.pT_edges[b + 1]) {
          if (std::isfinite(resp) && resp > 0.0)
            bin_resp[b].push_back(resp);
          break;
        }
      }
    }
    coarsen_bins(bin_resp, maps.cabs[alg], N_min);
  }

  const std::size_t npt = maps.pT_dis_edges.size() > 1
                              ? maps.pT_dis_edges.size() - 1
                              : 0;
  const std::size_t ny = maps.y_edges.size() > 1 ? maps.y_edges.size() - 1 : 0;
  if (npt > 0 && ny > 0) {
    for (int alg = 0; alg < 2; ++alg) {
      maps.cdis_pt_y[alg].assign(npt * ny, 1.0);
      std::vector<std::vector<double>> bin_resp(npt * ny);
      std::vector<std::size_t> counts(npt * ny, 0);

      for (const auto &r : records) {
        if (!r.is_train)
          continue;
        if (r.alg != alg)
          continue;
        if (r.pT_truth < pT_min)
          continue;
        if (std::abs(r.eta_reco) >= cfg.max_abs_eta)
          continue;
        if (!std::isfinite(r.y))
          continue;

        const double crel =
            interpolate_bins(maps.eta_edges, maps.crel[alg], r.eta_reco);
        const double pT_after_eta = r.pT_reco_raw * crel;
        if (pT_after_eta < pT_min || pT_after_eta > pT_max)
          continue;
        const double cabs =
            interpolate_bins(maps.pT_edges, maps.cabs[alg], pT_after_eta);
        const double pT0 = pT_after_eta * cabs;
        const double resp = pT0 / r.pT_truth;
        if (!std::isfinite(resp) || resp <= 0.0)
          continue;

        const std::size_t ipt = bin_index_clamped(maps.pT_dis_edges, pT_after_eta);
        const std::size_t iy = bin_index_clamped(maps.y_edges, r.y);
        const std::size_t index = ipt * ny + iy;
        if (index >= bin_resp.size())
          continue;
        bin_resp[index].push_back(resp);
        counts[index] += 1;
      }

      std::vector<double> c_raw(npt * ny, 1.0);
      std::vector<bool> active(npt * ny, false);
      for (std::size_t i = 0; i < bin_resp.size(); ++i) {
        if (bin_resp[i].size() < N_min)
          continue;
        const double med = median_value(bin_resp[i]);
        if (!std::isfinite(med) || med <= 0.0)
          continue;
        c_raw[i] = 1.0 / med;
        active[i] = true;
      }

      for (std::size_t ipt = 0; ipt < npt; ++ipt) {
        std::vector<double> c_fill(ny, 1.0);
        bool has_active = false;
        for (std::size_t iy = 0; iy < ny; ++iy) {
          const std::size_t idx = ipt * ny + iy;
          if (!active[idx])
            continue;
          c_fill[iy] = c_raw[idx];
          has_active = true;
        }
        if (has_active) {
          for (std::size_t iy = 0; iy < ny; ++iy) {
            const std::size_t idx = ipt * ny + iy;
            if (active[idx])
              continue;
            std::size_t best = ny;
            std::size_t best_dist = std::numeric_limits<std::size_t>::max();
            for (std::size_t j = 0; j < ny; ++j) {
              const std::size_t jdx = ipt * ny + j;
              if (!active[jdx])
                continue;
              const std::size_t dist = (j > iy) ? (j - iy) : (iy - j);
              if (dist < best_dist) {
                best_dist = dist;
                best = j;
                if (best_dist == 0)
                  break;
              }
            }
            if (best < ny)
              c_fill[iy] = c_raw[ipt * ny + best];
          }
        }
        // Normalize per pT bin so weighted mean Cdis is 1.
        double wsum = 0.0;
        double csum = 0.0;
        for (std::size_t iy = 0; iy < ny; ++iy) {
          const std::size_t idx = ipt * ny + iy;
          const double w = static_cast<double>(counts[idx]);
          wsum += w;
          csum += w * c_fill[iy];
        }
        const double norm = (wsum > 0.0) ? (csum / wsum) : 1.0;
        for (std::size_t iy = 0; iy < ny; ++iy) {
          const std::size_t idx = ipt * ny + iy;
          maps.cdis_pt_y[alg][idx] = c_fill[iy] / norm;
        }
      }
    }
  }

  // Reclose pT after the DIS stage to keep inclusive scale stable.
  for (int alg = 0; alg < 2; ++alg) {
    maps.cabs_final[alg].resize(maps.pT_edges.size() - 1, 1.0);
    std::vector<std::vector<double>> bin_resp(maps.pT_edges.size() - 1);
    for (const auto &r : records) {
      if (!r.is_train)
        continue;
      if (r.alg != alg)
        continue;
      if (r.pT_truth < pT_min)
        continue;
      if (std::abs(r.eta_reco) >= cfg.max_abs_eta)
        continue;

      const double crel =
          interpolate_bins(maps.eta_edges, maps.crel[alg], r.eta_reco);
      const double pT_after_eta = r.pT_reco_raw * crel;
      if (pT_after_eta < pT_min || pT_after_eta > pT_max)
        continue;
      const double cabs =
          interpolate_bins(maps.pT_edges, maps.cabs[alg], pT_after_eta);
      const double pT0 = pT_after_eta * cabs;
      double cdis = 1.0;
      if (!maps.cdis_pt_y[alg].empty() && npt > 0 && ny > 0) {
        const std::size_t ipt = bin_index_clamped(maps.pT_dis_edges, pT_after_eta);
        const std::size_t iy = bin_index_clamped(maps.y_edges, r.y);
        const std::size_t idx = ipt * ny + iy;
        if (idx < maps.cdis_pt_y[alg].size())
          cdis = maps.cdis_pt_y[alg][idx];
      }
      const double pT1 = pT0 * cdis;
      const double resp = pT1 / r.pT_truth;
      if (!std::isfinite(resp) || resp <= 0.0)
        continue;

      for (std::size_t b = 0; b + 1 < maps.pT_edges.size(); ++b) {
        if (pT_after_eta >= maps.pT_edges[b] &&
            pT_after_eta < maps.pT_edges[b + 1]) {
          bin_resp[b].push_back(resp);
          break;
        }
      }
    }
    coarsen_bins(bin_resp, maps.cabs_final[alg], N_min);
  }

  if (cfg.renorm) {
    for (int alg = 0; alg < 2; ++alg) {
      double sum_resp = 0.0;
      std::size_t n_resp = 0;
      for (const auto &r : records) {
        if (!r.is_train)
          continue;
        if (r.alg != alg)
          continue;
        if (std::abs(r.eta_reco) >= cfg.max_abs_eta)
          continue;
        if (r.pT_truth <= 0.0)
          continue;

        const double crel =
            interpolate_bins(maps.eta_edges, maps.crel[alg], r.eta_reco);
        const double pT_after_eta = r.pT_reco_raw * crel;
        if (pT_after_eta < pT_min || pT_after_eta > pT_max)
          continue;
        if (r.pT_truth < pT_min)
          continue;

        const double cabs =
            interpolate_bins(maps.pT_edges, maps.cabs[alg], pT_after_eta);
        double cdis = 1.0;
        if (!maps.cdis_pt_y[alg].empty() && npt > 0 && ny > 0) {
          const std::size_t ipt =
              bin_index_clamped(maps.pT_dis_edges, pT_after_eta);
          const std::size_t iy = bin_index_clamped(maps.y_edges, r.y);
          const std::size_t idx = ipt * ny + iy;
          if (idx < maps.cdis_pt_y[alg].size())
            cdis = maps.cdis_pt_y[alg][idx];
        }
        const double cfinal =
            interpolate_bins(maps.pT_edges, maps.cabs_final[alg], pT_after_eta);
        const double resp =
            (r.pT_reco_raw * crel * cabs * cdis * cfinal) / r.pT_truth;
        if (!std::isfinite(resp) || resp <= 0.0)
          continue;
        sum_resp += resp;
        ++n_resp;
      }
      if (n_resp == 0)
        continue;

      const double mean_resp = sum_resp / static_cast<double>(n_resp);
      if (!std::isfinite(mean_resp) || mean_resp <= 0.0)
        continue;
      const double scale = 1.0 / mean_resp;
      for (auto &v : maps.cabs_final[alg])
        v *= scale;
      std::cout << "[jet_calib] Applied global scale " << scale
                << " to C_abs_final for alg " << alg
                << " to enforce unit mean response\n";
    }
  }

  return maps;
}

}  // namespace jet_tools

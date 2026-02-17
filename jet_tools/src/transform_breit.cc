#include "jet_tools/include/transform_breit.h"

#include <cmath>
#include <limits>

#include <Math/Vector3D.h>

#include "jet_tools/include/math_helpers.h"

namespace jet_tools {

using Vec3 = ROOT::Math::XYZVector;

namespace {

// Compute pseudorapidity with guards against bad inputs.
double safe_eta(double px, double py, double pz) {
  const double pT = std::hypot(px, py);
  if (!std::isfinite(pT) || pT <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double eta = std::asinh(pz / pT);
  return std::isfinite(eta) ? eta : std::numeric_limits<double>::quiet_NaN();
}

// Compute azimuth with guards against bad inputs.
double safe_phi(double px, double py) {
  const double phi = std::atan2(py, px);
  return std::isfinite(phi) ? phi : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace

// On success, returns true and fills the output boost and rotation (by reference).
// Breit frame definition:
// 1. Boost such that 2xP + q = 0 (momentum of quark + momentum of boson = 0).
// 2. Rotate such that:
//    - z-axis is along -q (proton direction).
//    - y-axis is perpendicular to the lepton scattering plane.
//    - x-axis completes the right-handed system.
bool transform_breit(const P4 &p_i, const P4 &e_i, const P4 &e_f,
                     ROOT::Math::Boost &out_boost,
                     ROOT::Math::Rotation3D &out_rot, BreitXSource x_source,
                     double x_input, double *Q_out) {
  // Build the exchanged 4-momentum q with the virtual boson.
  const P4 q = e_i - e_f;
  // Compute Q2.
  const double Q2 = -q.M2();
  if (!std::isfinite(Q2) || Q2 <= 0.0)
    return false;

  // Determine Bjorken x.
  double bjorken_x = std::numeric_limits<double>::quiet_NaN();
  // If provided a Bjorken x input, use the input.
  if (x_source == BreitXSource::Input) {
    bjorken_x = x_input;
  } else {
    // Else, derive Bjorken x using invariants.
    const double denom = 2.0 * minkowski_dot(p_i, q);
    if (!std::isfinite(denom) || std::abs(denom) < 1e-12)
      return false;
    bjorken_x = Q2 / denom;
  }
  if (!std::isfinite(bjorken_x) || bjorken_x <= 0.0)
    return false;

  //---------------------------------------------------------------------------
  // Main Breit boost
  //---------------------------------------------------------------------------
  const auto P3 = p_i.Vect(); // incoming hadron 3-momentum
  const auto q3 = q.Vect();   // exchanged boson 3-momentum

  // Time component.
  const double denom = (2.0 * bjorken_x * p_i.E() + q.E());
  if (!std::isfinite(denom) || std::abs(denom) < 1e-9)
    return false;

  // Build the boost to bring V-vec to zero in the boosted frame.
  out_boost = ROOT::Math::Boost(-(2.0 * bjorken_x * P3 + q3) * (1.0 / denom));

  // Apply the boost to the incoming lepton, outgoing lepton, and the momentum transfer q.
  const P4 e_i_b = out_boost * e_i;
  const P4 e_f_b = out_boost * e_f;
  const P4 q_b = out_boost * q;

  // Define Breit coordinate axis.
  const auto qhat = q_b.Vect(); // extract boosted q three-vector.
  if (qhat.R() <= 0.0)
    return false;
  // z-axis: opposite direction of boosted q.
  const auto zhat = -qhat.Unit();
  // y-axis: normal to the boosted lepton plane (e_in cross e_out)
  auto yhat = (e_i_b.Vect().Cross(e_f_b.Vect()));
  const double ymag = yhat.R();
  if (std::isfinite(ymag) && ymag > 0.0) {
    // Normalize the plane if well-defined.
    yhat = yhat * (1.0 / ymag);
  } else {
    // Fallback: choose a temporary vector not parallel to zhat, 
    // then take tmp cross zhat to get a stable y direction.
    Vec3 tmp(0.0, 1.0, 0.0);
    // If tmp is too aligned with zhat, swap tmp to another axis to avoid a near-zero vector.
    if (std::abs(zhat.Dot(tmp)) > 0.9)
      tmp = Vec3(1.0, 0.0, 0.0);
    // Build a perpendicular unit vector for yhat.
    yhat = (tmp.Cross(zhat)).Unit();
  }
  
  // x-axis: complete right-handed system (x = y cross z)
  auto xhat = yhat.Cross(zhat);
  if (xhat.R() <= 0.0)
    return false;
  // Normalize x.
  xhat = xhat.Unit();
  // Recompute y to enforce exact orthonormality (y = z cross x)
  yhat = zhat.Cross(xhat).Unit();

  // Build rotation.
  const ROOT::Math::Rotation3D Rinv(xhat, yhat, zhat);
  // Store the rotation that takes vectors into the Breit axes.
  out_rot = Rinv.Inverse();
  if (Q_out)
    *Q_out = std::sqrt(Q2);
  return true;
}

// Boost from breit back to lab.
P4 breit_to_lab(const P4 &p_breit, const ROOT::Math::Boost &boost,
                const ROOT::Math::Rotation3D &rot) {
  return boost.Inverse() * (rot.Inverse() * p_breit);
}

// Boost from lab to breit.
P4 lab_to_breit(const P4 &p_lab, const ROOT::Math::Boost &boost,
                const ROOT::Math::Rotation3D &rot) {
  return rot * (boost * p_lab);
}

// Boost a Breit frame jet back to lab frame.
SimpleJet breit_to_lab(const SimpleJet &jet_breit, const ROOT::Math::Boost &boost,
                       const ROOT::Math::Rotation3D &rot) {
  // 4-vectors.
  const P4 p_breit(jet_breit.px, jet_breit.py, jet_breit.pz, jet_breit.E);
  const P4 p_lab = breit_to_lab(p_breit, boost, rot);

  // Make output jet.
  SimpleJet out = jet_breit;

  // Update 4-vectors and compute other kinematics.
  out.px = p_lab.Px();
  out.py = p_lab.Py();
  out.pz = p_lab.Pz();
  out.E = p_lab.E();
  out.pT = std::hypot(out.px, out.py);
  out.eta = safe_eta(out.px, out.py, out.pz);
  out.phi = safe_phi(out.px, out.py);
  return out;
}

// Boost a jet from lab frame to Breit frame.
SimpleJet lab_to_breit(const SimpleJet &jet_lab, const ROOT::Math::Boost &boost,
                       const ROOT::Math::Rotation3D &rot) {
  // 4-vectors.
  const P4 p_lab(jet_lab.px, jet_lab.py, jet_lab.pz, jet_lab.E);
  const P4 p_breit = lab_to_breit(p_lab, boost, rot);

  // Make output jet.
  SimpleJet out = jet_lab;

  // Update 4-vectors and compute other kinematics.
  out.px = p_breit.Px();
  out.py = p_breit.Py();
  out.pz = p_breit.Pz();
  out.E = p_breit.E();
  out.pT = std::hypot(out.px, out.py);
  out.eta = safe_eta(out.px, out.py, out.pz);
  out.phi = safe_phi(out.px, out.py);
  return out;
}

} // namespace jet_tools

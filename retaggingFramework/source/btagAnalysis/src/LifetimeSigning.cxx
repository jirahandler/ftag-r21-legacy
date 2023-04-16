#include "../btagAnalysis/LifetimeSigning.hh"
#include "../btagAnalysis/BTagTrackAccessors.hh"

// the logic here is copied from https://goo.gl/iWLv5T

namespace {
  Amg::Vector3D vec3(const std::vector<float>& vec) {
    assert(vec.size() == 3);
    return Eigen::Vector3f(vec.data()).cast<double>();
  }
}

double get_3d_lifetime_sign(const SG::AuxElement& track,
                            const Amg::Vector3D& jet) {
  BTagTrackAccessors acc;
  auto trk_mom = vec3(acc.momentum(track));
  auto trk_dis = vec3(acc.displacement(track));
  double sign = (jet.cross(trk_mom)).dot(trk_mom.cross(-trk_dis));
  return sign >= 0 ? 1 : -1;
}
double get_2d_lifetime_sign(const SG::AuxElement& track,
                            const Amg::Vector3D& jet) {
  BTagTrackAccessors acc;
  double d0  = acc.d0(track);
  double phi = vec3(acc.momentum(track)).phi();
  double vs = std::sin( jet.phi() - phi )*d0;
  return vs >= 0 ? 1 : -1;
}

double get_z_lifetime_sign(const SG::AuxElement& track,
                           const Amg::Vector3D& jet) {
  BTagTrackAccessors acc;
  double z0 = acc.z0(track);
  double zs = (jet.eta() - vec3(acc.momentum(track)).eta()) * z0;
  return zs > 0 ? 1 : -1;
}

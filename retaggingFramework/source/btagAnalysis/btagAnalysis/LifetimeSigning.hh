#ifndef LIFETIME_SIGNING_HH
#define LIFETIME_SIGNING_HH

#include "AthContainers/AuxElement.h"
#include "GeoPrimitives/GeoPrimitives.h"

// gets the lifetime sign for a track (assumed decorated by the
// BTagTrackAugmenter) and the jet direction.
double get_3d_lifetime_sign(const SG::AuxElement& track,
                            const Amg::Vector3D& jet);
double get_2d_lifetime_sign(const SG::AuxElement& track,
                            const Amg::Vector3D& jet);
double get_z_lifetime_sign(const SG::AuxElement& track,
                           const Amg::Vector3D& jet);

#endif

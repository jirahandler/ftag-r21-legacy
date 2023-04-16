#ifndef BTAG_TRACK_ACCESSORS_HH
#define BTAG_TRACK_ACCESSORS_HH

#include "AthContainers/AuxElement.h"
// #include "GeoPrimitives/GeoPrimitives.h"

struct BTagTrackAccessors {
  SG::AuxElement::ConstAccessor< float > d0;
  SG::AuxElement::ConstAccessor< float > z0;
  SG::AuxElement::ConstAccessor< float > d0_sigma;
  SG::AuxElement::ConstAccessor< float > z0_sigma;

  SG::AuxElement::ConstAccessor< std::vector<float> > displacement;
  SG::AuxElement::ConstAccessor< std::vector<float> > momentum;

  BTagTrackAccessors():
    /*d0("btagIp_d0"),
    z0("btagIp_z0SinTheta"),
    d0_sigma("btagIp_d0Uncertainty"),
    z0_sigma("btagIp_z0SinThetaUncertainty"),
    displacement("btagIp_trackDisplacement"),
    momentum("btagIp_trackMomentum")*/
    // name change -- A.X.
    d0("btag_ip_d0"),
    z0("btag_ip_z0"),
    d0_sigma("btag_ip_d0_sigma"),
    z0_sigma("btag_ip_z0_sigma"),
    displacement("btag_track_displacement"),
    momentum("btag_track_momentum")
    {
    }

};

#endif

#ifndef TRACK_BRANCH_BUFFER_HH
#define TRACK_BRANCH_BUFFER_HH

#include <vector>

struct TrackBranchBuffer {
  std::vector<int> *ntrk;

  std::vector<std::vector<float> > *pt;
  std::vector<std::vector<float> > *eta;
  std::vector<std::vector<float> > *theta;
  std::vector<std::vector<float> > *phi;
  std::vector<std::vector<float> > *qoverp;

  std::vector<std::vector<float> > *charge;

  std::vector<std::vector<float> > *chi2;
  std::vector<std::vector<float> > *ndf;

  // note: these don't use the EDM name, which is unnessisarily confusing
  // but we do it for backward compatibitlity.
  std::vector<std::vector<int> > *nNextToInnHits;
  std::vector<std::vector<int> > *nInnHits;
  std::vector<std::vector<int> > *nBLHits; // soo this will be deprecated
  std::vector<std::vector<int> > *nsharedBLHits;
  std::vector<std::vector<int> > *nsplitBLHits;
  std::vector<std::vector<int> > *nPixHits;
  std::vector<std::vector<int> > *nPixHoles;
  std::vector<std::vector<int> > *nsharedPixHits;
  std::vector<std::vector<int> > *nsplitPixHits;
  std::vector<std::vector<int> > *nSCTHits;
  std::vector<std::vector<int> > *nSCTHoles;
  std::vector<std::vector<int> > *nsharedSCTHits;

  // note: where possible (and when adding new branches) it's better to
  // use the EDM name as we do here.
  std::vector<std::vector<int> > *expectBLayerHit;
  std::vector<std::vector<int> > *expectInnermostPixelLayerHit;
  std::vector<std::vector<int> > *expectNextToInnermostPixelLayerHit;

  // actual d0 variables (not lifetime-signed)
  std::vector<std::vector<float> > *d0;
  std::vector<std::vector<float> > *z0;

  std::vector<std::vector<float> > *ip3d_d0;
  std::vector<std::vector<float> > *ip3d_d02D;
  std::vector<std::vector<float> > *ip3d_z0;
  std::vector<std::vector<float> > *ip3d_d0sig;
  std::vector<std::vector<float> > *ip3d_z0sig;

  std::vector<int>* v_jet_btag_ntrk;
  std::vector<int>* v_jet_ip3d_ntrk;
  std::vector<std::vector<int> > *v_jet_trk_algo         ;
  std::vector<std::vector<float> > *v_jet_trk_vtx_X      ;
  std::vector<std::vector<float> > *v_jet_trk_vtx_Y      ;
  std::vector<std::vector<float> > *v_jet_trk_vtx_Z      ;
  std::vector<std::vector<int> >* v_jet_trk_IP3D_grade   ;
  std::vector<std::vector<float> >* v_jet_trk_IP3D_llr   ;
  std::vector<std::vector<int> > *v_jet_trk_pdg_id       ;
  std::vector<std::vector<int> > *v_jet_trk_barcode      ;
  std::vector<std::vector<int> > *v_jet_trk_parent_pdgid ;

};


#endif // TRACK_BRANCH_BUFFER_HH

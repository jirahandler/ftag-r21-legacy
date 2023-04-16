#ifndef SOFTMUON_BRANCH_BUFFER_HH
#define SOFTMUON_BRANCH_BUFFER_HH

#include <vector>

struct SoftMuonBranchBuffer {

  std::vector<double> *v_jet_mu_smt;
  std::vector<float> *v_jet_mu_dR;
  std::vector<float> *v_jet_mu_pTrel;
  std::vector<float> *v_jet_mu_qOverPratio;
  std::vector<float> *v_jet_mu_mombalsignif;
  std::vector<float> *v_jet_mu_scatneighsignif;
  std::vector<float> *v_jet_mu_VtxTyp;
  std::vector<float> *v_jet_mu_pt;
  std::vector<float> *v_jet_mu_eta;
  std::vector<float> *v_jet_mu_phi;
  std::vector<float> *v_jet_mu_d0;
  std::vector<float> *v_jet_mu_z0;
  std::vector<float> *v_jet_mu_parent_pdgid;
  std::vector<float> *v_jet_mu_ID_qOverP_var;
  std::vector<float> *v_jet_mu_muonType;



};

#endif // SOFTMUON_BRANCH_BUFFER_HH

#ifndef JETPROPERTIES_BRANCH_BUFFER_HH
#define JETPROPERTIES_BRANCH_BUFFER_HH

#include <vector>

struct JetPropertiesBranchBuffer {

  std::vector<float> *v_jet_pt;
  std::vector<float> *v_jet_eta;
  std::vector<float> *v_jet_phi;
  std::vector<float> *v_jet_E;
  std::vector<float> *v_jet_pt_orig;
  std::vector<float> *v_jet_eta_orig;
  std::vector<float> *v_jet_phi_orig;
  std::vector<float> *v_jet_E_orig;
  std::vector<int> 	 *v_jet_LabDr_HadF;
  std::vector<int>   *v_jet_DoubleHadLabel;
  std::vector<float> *v_jet_JVT;

  std::vector<float> *v_jet_m;
  std::vector<float> *v_jet_nConst;
  std::vector<float> *v_jet_dRiso;
  std::vector<int>   *v_jet_truthMatch;
  std::vector<int>   *v_jet_isPU;
  std::vector<int>   *v_jet_aliveAfterOR;
  std::vector<int>   *v_jet_aliveAfterORmu;
  std::vector<int>   *v_jet_isBadMedium;
  std::vector<float> *v_jet_truthPt;
  std::vector<float> *v_jet_dRminToB;
  std::vector<float> *v_jet_dRminToC;
  std::vector<float> *v_jet_dRminToT;

  std::vector<int>   *v_jet_isLLPDecayProd;
  std::vector<float> *v_jet_dRtoLLPDecayProd;
  std::vector<int>   *v_jet_truthLLPJetLabel;
  std::vector<float> *v_jet_truthLLP_Decay_x;
  std::vector<float> *v_jet_truthLLP_Decay_y;
  std::vector<float> *v_jet_truthLLP_Decay_z;

};

#endif // JETPROPERTIES_BRANCH_BUFFER_HH

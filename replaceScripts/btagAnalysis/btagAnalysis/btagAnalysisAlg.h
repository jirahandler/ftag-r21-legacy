#ifndef BTAGANALYSIS_BTAGANALYSISALG_H
#define BTAGANALYSIS_BTAGANALYSISALG_H

// Additions from Dan: classes to manage branches
#include "JetPropertiesBranches.hh"
#include "TaggerScoreBranches.hh"
#include "ClusterBranches.hh"
#include "SubjetBranches.hh"
#include "TrackCovBranches.hh"
#include "TrackBranches.hh"
#include "BHadronBranches.hh"
#include "KShortBranches.hh"
#include "KShortRecoBranches.hh"
#include "ImpactParameterBranches.hh"
#include "JetFitterBranches.hh"
#include "SubstructureMomentBranches.hh"
#include "UnclusteredVertexBranches.hh"
#include "SVBranches.hh"
#include "SoftMuonBranches.hh"
#include "BTagTrackAccessors.hh"
#include "ExampleBranches.hh"

#include "AthenaBaseComps/AthHistogramAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "TFile.h"
#include "TTree.h"

#ifndef __MAKECINT__
#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/TruthParticle.h"
#endif // not __MAKECINT__

#include "TrigDecisionTool/TrigDecisionTool.h"

// forward declarations
class IJetSelector;
class IJetCalibrationTool;
namespace InDet { class IInDetTrackSelectionTool; }
namespace CP {
  class ITrackVertexAssociationTool;
  class IPileupReweightingTool;
}

namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}

class IJetUpdateJvt;

class ArbitraryJetBranches;

class btagAnalysisAlg: public ::AthHistogramAlgorithm {
 public:
  btagAnalysisAlg( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~btagAnalysisAlg();

  virtual StatusCode  initialize();
  virtual StatusCode  initializeTools();
  virtual StatusCode  execute();
  virtual StatusCode  finalize();

private:
  TTree* m_tree;
  std::string m_stream;

  std::string m_jetCollectionName;
  std::string m_track_associator;

  float m_jetPtCut;


  //Flags
  bool m_calibrateJets;
  bool m_cleanJets;
  bool m_clean_parent_jet;
  bool m_doJVT;

  bool m_EventInfo;
  bool m_retriveTruthJets;
  bool m_JetProperties;
  bool m_TaggerScoresInfo;
  bool m_exampleBranchInfo;
  bool m_ImpactParameterInfo;
  bool m_SVInfo;
  bool m_JetFitterInfo;
  bool m_SoftMuoninfo;
  bool m_bHadInfo;
  bool m_bHadExtraInfo;
  bool m_kshortInfo;
  bool m_TrackInfo;
  bool m_TrackCovariance;

  bool m_showDebug;
  bool m_access_btag_object;
  int m_n_required_si_hits;

  //default value handling
  std::map<std::string, double > m_defaultDict;
  bool m_replaceDefaults;

  /// tool handle for jet cleaning tool
  ToolHandle< IJetSelector > m_jetCleaningTool;

  /// tool handle for jet calibration tool
  ToolHandle< IJetCalibrationTool > m_jetCalibrationTool;

  //JVT update tool
  ToolHandle<IJetUpdateJvt> m_jvt;

  //PileUp Reweighting
  ToolHandle<CP::IPileupReweightingTool> m_PUtool;



  //Event Info
  int runnumber;
  int eventnumber;
  int mcchannel;
  float mcweight;
  float mu;
  int Act_mu;
  float PV_x;
  float PV_y;
  float PV_z;
  float truth_PV_x;
  float truth_PV_y;
  float truth_PV_z;
  int njets;

  //Branches
  ExampleBranches m_example_branches;
  JetPropertiesBranches m_jet_properties_branches;
  TaggerScoreBranches m_tagger_scores_branches;
  ImpactParameterBranches m_impact_parameter_branches;
  JetFitterBranches m_jetfitter_branches;
  SoftMuonBranches m_softmuon_branches;
  KShortBranches m_kshort_branches;
  KShortRecoBranches m_kshortreco_branches;
  BHadronBranches m_bhadron_branches;
  // subjet dumper
  std::map<std::string, std::string> m_subjet_collections;
  std::vector<std::pair<std::string, SubjetBranches*> > m_subjet_branches;

  // SVx branches
  // first is the prefix in the ntuple, second is edm name
  std::map<std::string, std::string> m_svx_collections;
  std::vector<SVBranches*> m_svx_branches;
  TrackBranches m_track_branches;
  TrackCovBranches m_track_cov_branches;

  ArbitraryJetBranches* m_arb_branches;
  std::vector<std::string> m_arb_double_names;
  std::vector<std::string> m_arb_float_vec_names;

  // A.X.
  int nmymc;
  std::vector<float> *mymc_pt;
  std::vector<float> *mymc_eta;
  std::vector<float> *mymc_phi;
  std::vector<int> *mymc_pdgId;
  std::vector<float> *mymc_decayVtx_x;
  std::vector<float> *mymc_decayVtx_y;
  std::vector<float> *mymc_decayVtx_z;
  std::vector<int> *mymc_ix1;
  std::vector<int> *mymc_ix2;
  int nmymc1;
  std::vector<float> *mymc1_pt;
  std::vector<float> *mymc1_eta;
  std::vector<float> *mymc1_phi;
  std::vector<int> *mymc1_pdgId;



  bool isFromWZ( const xAOD::TruthParticle* particle );
  const xAOD::Jet* GetParentJet(const xAOD::Jet* Jet, std::string Keyname);

};

#endif //> !BTAGANALYSIS_BTAGANALYSISALG_H

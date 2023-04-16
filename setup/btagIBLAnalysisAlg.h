#ifndef BTAGIBLANALYSIS_BTAGIBLANALYSISALG_H
#define BTAGIBLANALYSIS_BTAGIBLANALYSISALG_H 1

// Additions from Dan: classes to manage branches
#include "ClusterBranches.hh"
#include "SubjetBranches.hh"
#include "TrackCovBranches.hh"
#include "TrackBranches.hh"
#include "BHadronBranches.hh"
#include "KShortBranches.hh"
#include "KShortRecoBranches.hh"
#include "JetFitterBranches.hh"
#include "SubstructureMomentBranches.hh"
#include "UnclusteredVertexBranches.hh"
#include "SVBranches.hh"
#include "BTagTrackAccessors.hh"

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
//namespace Trk  { class ITrackToVertexIPEstimator; }

namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}

class IGoodRunsListSelectionTool;
class IJetUpdateJvt;

class ArbitraryJetBranches;

enum TAGGERALGO{ IP2D=0,
		 IP3D,
		 SV0,
		 SV1,
		 JF };

// moved to BHadronBranches
// enum TRKORIGIN{ PUFAKE=-1,
// 		FROMB,
// 		FROMC,
// 		FRAG,
// 		GEANT };

enum LFCalibType{
  nothing = 0x0,
  hasKShort = 0x1,
  hasLambda = 0x2,
  hasConversion = 0x4,
  hasHadMatInt = 0x8,
  hasFake = 0x10 };

class btagIBLAnalysisAlg: public ::AthHistogramAlgorithm {
 public:
  btagIBLAnalysisAlg( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~btagIBLAnalysisAlg();

  virtual StatusCode  initialize();
  virtual StatusCode  execute();
  virtual StatusCode  finalize();

private:
  TFile* output;
  TTree* tree;
  std::string m_stream;

  // general event info
  int runnumber;
  int eventnumber;
  int mcchannel;
  double mcweight;
  int lbn;
  int coreFlag;
  int larError;
  int tileError;
  int npv;
  double mu;
  int Act_mu;
  double PV_x;
  double PV_y;
  double PV_z;
  double truth_PV_x;
  double truth_PV_y;
  double truth_PV_z;
  float  truth_LeadJet_pt;

  bool* v_L1trigger;
  std::vector<std::string> v_L1triggerNames;

  // jet info
  int njets;
  int nbjets;
  int nbjets_q;
  int nbjets_HadI;
  int nbjets_HadF;
  std::vector<float> *v_jet_pt;
  std::vector<float> *v_jet_eta;
  std::vector<float> *v_jet_pt_orig;
  std::vector<float> *v_jet_eta_orig;
  std::vector<float> *v_jet_phi_orig;
  std::vector<float> *v_jet_E_orig;
  std::vector<float> *v_jet_sumtrkS_pt;
  std::vector<float> *v_jet_sumtrkV_pt;
  std::vector<float> *v_jet_sumtrkV_eta;
  std::vector<float> *v_jet_sumtrkV_phi;
  std::vector<int>   *v_jet_sumtrk_ntrk;
  std::vector<float> *v_jet_phi;
  std::vector<float> *v_jet_E;
  std::vector<float> *v_jet_m;
  std::vector<int> *v_jet_nConst;
  std::vector<int> *v_jet_truthflav;

  std::vector<int> *v_jet_GhostL_q;
  std::vector<int> *v_jet_GhostL_HadI;
  std::vector<int> *v_jet_GhostL_HadF;
  std::vector<int> *v_jet_LabDr_HadF;
  std::vector<int> *v_jet_DoubleHadLabel;
  std::vector<int> *v_jet_aliveAfterOR;
  std::vector<int> *v_jet_aliveAfterORmu;
  std::vector<int> *v_jet_truthMatch;
  std::vector<int> *v_jet_isPU;
  std::vector<int> *v_jet_isBadMedium;
  std::vector<float> *v_jet_truthPt;
  std::vector<float> *v_jet_dRiso;
  std::vector<float> *v_jet_JVT;
  std::vector<float> *v_jet_JVF;

  std::vector<float> *v_jet_dRminToB;
  std::vector<float> *v_jet_dRminToC;
  std::vector<float> *v_jet_dRminToT;

  // IP2D
  std::vector<float> *v_jet_ip2d_pb;
  std::vector<float> *v_jet_ip2d_pc;
  std::vector<float> *v_jet_ip2d_pu;
  std::vector<float> *v_jet_ip2d_llr;

  // IP3D
  std::vector<float> *v_jet_ip3d_pb;
  std::vector<float> *v_jet_ip3d_pc;
  std::vector<float> *v_jet_ip3d_pu;
  std::vector<float> *v_jet_ip3d_llr;

  // SV0
  // std::vector<std::vector<float> > *v_jet_sv0_vtxx;
  // std::vector<std::vector<float> > *v_jet_sv0_vtxy;
  // std::vector<std::vector<float> > *v_jet_sv0_vtxz;

  // SV1
  std::vector<float> *v_jet_sv1_pb;
  std::vector<float> *v_jet_sv1_pc;
  std::vector<float> *v_jet_sv1_pu;
  std::vector<float> *v_jet_sv1_llr;
  // std::vector<std::vector<float> > *v_jet_sv1_vtxx;
  // std::vector<std::vector<float> > *v_jet_sv1_vtxy;
  // std::vector<std::vector<float> > *v_jet_sv1_vtxz;



  // Other
  std::vector<double> *v_jet_sv1ip3d;
  std::vector<double> *v_jet_mv1;
  std::vector<double> *v_jet_mv1c;
  std::vector<double> *v_jet_mv2c00;
  std::vector<double> *v_jet_mv2c10;
  std::vector<double> *v_jet_mv2c10mu;
  std::vector<double> *v_jet_mv2c10rnn;
  std::vector<double> *v_jet_mv2c20;
  std::vector<double> *v_jet_mv2c100;
  std::vector<double> *v_jet_mv2cl100;
  std::vector<double> *v_jet_mv2m_pb;
  std::vector<double> *v_jet_mv2m_pc;
  std::vector<double> *v_jet_mv2m_pu;
  std::vector<double> *v_jet_mvb;

  std::vector<float> *v_jet_dl1_pb;
  std::vector<float> *v_jet_dl1_pc;
  std::vector<float> *v_jet_dl1_pu;

  std::vector<float> *v_jet_dl1mu_pb;
  std::vector<float> *v_jet_dl1mu_pc;
  std::vector<float> *v_jet_dl1mu_pu;

  std::vector<float> *v_jet_dl1rnn_pb;
  std::vector<float> *v_jet_dl1rnn_pc;
  std::vector<float> *v_jet_dl1rnn_pu;

  std::vector<float> *v_jet_mv2c20flip;
  std::vector<float> *v_jet_mv2c10flip;

  //MSV
  std::vector<double> *v_jet_multisvbb1;
  std::vector<double> *v_jet_multisvbb2;
  std::vector<int> *v_jet_msv_N2Tpair;
  std::vector<float> *v_jet_msv_energyTrkInJet;
  std::vector<int> *v_jet_msv_nvsec;
  std::vector<float> *v_jet_msv_normdist;
  std::vector<std::vector<float> > *v_jet_msv_vtx_cov0;
  std::vector<std::vector<float> > *v_jet_msv_vtx_cov1;
  std::vector<std::vector<float> > *v_jet_msv_vtx_cov2;
  std::vector<std::vector<float> > *v_jet_msv_vtx_cov3;
  std::vector<std::vector<float> > *v_jet_msv_vtx_cov4;
  std::vector<std::vector<float> > *v_jet_msv_vtx_cov5;
  std::vector<std::vector<float> > *v_jet_msv_vtx_mass;
  std::vector<std::vector<float> > *v_jet_msv_vtx_efrc;
  std::vector<std::vector<float> > *v_jet_msv_vtx_ntrk;
  std::vector<std::vector<float> > *v_jet_msv_vtx_pt;
  std::vector<std::vector<float> > *v_jet_msv_vtx_eta;
  std::vector<std::vector<float> > *v_jet_msv_vtx_phi;
  std::vector<std::vector<float> > *v_jet_msv_vtx_dls;
  std::vector<std::vector<float> > *v_jet_msv_vtx_x;
  std::vector<std::vector<float> > *v_jet_msv_vtx_y;
  std::vector<std::vector<float> > *v_jet_msv_vtx_z;
  std::vector<std::vector<float> > *v_jet_msv_vtx_chi;
  std::vector<std::vector<float> > *v_jet_msv_vtx_ndf;

  // Exktbb
  std::vector<double> *v_jet_ExKtbb_Hbb_DoubleMV2c20;
  std::vector<double> *v_jet_ExKtbb_Hbb_SingleMV2c20;
  std::vector<double> *v_jet_ExKtbb_Hbb_MV2Only;
  std::vector<double> *v_jet_ExKtbb_Hbb_MV2andJFDRSig;
  std::vector<double> *v_jet_ExKtbb_Hbb_MV2andTopos;

  // track info
  std::vector<int>   *v_jet_btag_ntrk;
  std::vector<int>   *v_jet_LFCalibType;

  std::vector<std::vector<float> > *v_jet_trk_dr;
  std::vector<std::vector<int> > *v_jet_trk_assoc_msv;

  std::vector<std::vector<int> > *v_jet_trk_algo;
  //std::vector<std::vector<int> > *v_jet_trk_orig; moved to BHadronBranches
  std::vector<std::vector<int> > *v_jet_trk_pdg_id;
  std::vector<std::vector<int> > *v_jet_trk_barcode;
  std::vector<std::vector<int> > *v_jet_trk_parent_pdgid;
  std::vector<std::vector<int> > *v_jet_trk_is_tracking_cp_loose;
  std::vector<std::vector<float> > *v_jet_trk_vtx_dx;
  std::vector<std::vector<float> > *v_jet_trk_vtx_dy;
  std::vector<std::vector<float> > *v_jet_trk_vtx_X;
  std::vector<std::vector<float> > *v_jet_trk_vtx_Y;
  std::vector<std::vector<float> > *v_jet_trk_vtx_Z;

  std::vector<std::vector<float> > *v_jet_trk_d0_truth;
  std::vector<std::vector<float> > *v_jet_trk_z0_truth;

  std::vector<std::vector<int> > *v_jet_trk_IP3D_grade;
  std::vector<std::vector<float> > *v_jet_trk_IP3D_d0;
  std::vector<std::vector<float> > *v_jet_trk_IP3D_d02D;
  std::vector<std::vector<float> > *v_jet_trk_IP3D_z0;
  std::vector<std::vector<float> > *v_jet_trk_IP3D_d0sig;
  std::vector<std::vector<float> > *v_jet_trk_IP3D_z0sig;
  std::vector<std::vector<float> > *v_jet_trk_IP3D_llr;

  //std::vector<std::vector<int> > *v_jet_trk_jf_Vertex; moved to JetFitterBranches

  // those are just quick accessors
  std::vector<int>   *v_jet_sv1_ntrk;
  std::vector<int>   *v_jet_ip3d_ntrk;
  std::vector<int>   *v_jet_jf_ntrk;

  // MVb variables
  std::vector<float> *v_jet_width;
  std::vector<int>   *v_jet_n_trk_sigd0cut;
  std::vector<float> *v_jet_trk3_d0sig;
  std::vector<float> *v_jet_trk3_z0sig;
  std::vector<float> *v_jet_sv_scaled_efc;
  std::vector<float> *v_jet_jf_scaled_efc;

  // additions by Andrea
  std::vector<double> *v_jet_mu_smt;
  std::vector<float> *v_jet_mu_truthflav;
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
  std::vector<float> *v_jet_mu_assJet_pt;
  // additions by nikola
  std::vector<int> *v_jet_mu_fatjet_nMu;
  std::vector<float> *v_jet_mu_fatjet_pTmax_pT;
  std::vector<float> *v_jet_mu_fatjet_pTmax_pTrel;
  std::vector<float> *v_jet_mu_fatjet_pTmax_pTrelFrac;

  void clearvectors();

#ifndef __MAKECINT__
  const xAOD::TruthParticle*  truthParticle(const xAOD::TrackParticle *trkPart) const;
  void GetParentTracks(const xAOD::TruthParticle* part,
		       std::vector<const xAOD::TruthParticle*> &tracksFromB,
		       std::vector<const xAOD::TruthParticle*> &tracksFromC,
		       bool isfromC, std::string indent);
  bool decorateTruth(const xAOD::TruthParticle & particle);
  int parent_classify(const xAOD::TruthParticle * theParticle);
  //bool particleInCollection( const xAOD::TrackParticle *trkPart,
  //std::vector< ElementLink< xAOD::TrackParticleContainer > > trkColl);
  void fillGhostTracks(const xAOD::Jet& jet,
                       const xAOD::Vertex& vx);


#endif // not __MAKECINT__


  /// from outside
  bool m_reduceInfo;               // if set to true is allows to run over xAOD and not crashing when info are missing
  bool m_essentialInfo;            // basically as slim as possible ntuple which will only allow to make efficiency plots
  bool m_subjetInfo;
  bool m_dumpCaloInfo;
  bool m_dumpTrackCovariance;
  bool m_dumpGATracks;
  bool m_doMSV;                    // if set to true it includes variables from multi SV tagger
  bool m_rel20;                    // if set to true code works for rel20, if set to false it will work for rel19
  bool m_SMT;
  bool m_bHadronInfo;
  bool m_kShortInfo;
  bool m_kShortRecoInfo;

  std::string m_jetCollectionName; // name of the jet collection to work with
  float m_jetPtCut;                // pT cut to apply
  bool m_calibrateJets;
  bool m_cleanJets;
  bool m_clean_parent_jet;
  std::string m_track_associator;

  bool m_saveTrkAlgInfo;           // Force to save available information about tracks entering into the algorithms
  bool m_saveJetLFCalibType;  // Save jet clarissification according to LFCalibTyep flag

  std::string m_triggerLogic;

  // additions by Dan: branch collections

  // SVx branches
  // first is the prefix in the ntuple, second is edm name
  std::map<std::string, std::string> m_svx_collections;
  std::vector<SVBranches*> m_svx_branches;

  // B-hadron quantities
  BHadronBranches m_bhadron_branches;

  // K0 short quantities
  KShortBranches m_kshort_branches;
  KShortRecoBranches m_kshortreco_branches;

  JetFitterBranches m_jetfitter_branches;

  // cluster dumper
  ClusterBranches m_cluster_branches;
  SubstructureMomentBranches m_substructure_moment_branches;
  // subjet dumper
  // first string is the name we call it, second is the EDM name
  std::map<std::string, std::string> m_subjet_collections;
  std::vector<std::pair<std::string, SubjetBranches*> > m_subjet_branches;
  // track dumper
  TrackBranches m_track_branches;
  TrackCovBranches m_track_cov_branches;
  TrackBranches m_ga_track_branches;
  TrackCovBranches m_ga_track_cov_branches;
  // must be initialized after the constructor
  ArbitraryJetBranches* m_arb_branches;
  std::vector<std::string> m_arb_double_names;
  std::vector<std::string> m_arb_float_vec_names;

  /// tool handle for jet cleaning tool
  ToolHandle< IJetSelector > m_jetCleaningTool;

  /// tool handle for jet calibration tool
  ToolHandle< IJetCalibrationTool > m_jetCalibrationTool;

  /** InDetTrackSelectorTool (temporary: to be moved to a separate Tool) */
  ToolHandle< InDet::IInDetTrackSelectionTool > m_InDetTrackSelectorTool;
  ToolHandle< InDet::IInDetTrackSelectionTool > m_CPTrackingLooseLabel;

  /** TrackVertex associator (temporary: to be moved to a separate Tool) */
  ToolHandle< CP::ITrackVertexAssociationTool > m_TightTrackVertexAssociationTool;

  /** GP: Tool for the estimation of the IPs to the Vertex */
  //ToolHandle< Trk::ITrackToVertexIPEstimator > m_trackToVertexIPEstimator;

  ToolHandle< Trig::TrigDecisionTool > m_tdt;

  std::string m_GRLname;
  ToolHandle<IGoodRunsListSelectionTool> m_GRLSelectionTool;

  ToolHandle<IJetUpdateJvt> m_jvt;

  ToolHandle<CP::IPileupReweightingTool> m_PUtool;

  BTagTrackAccessors m_track_accessors;

  // determine whether particle is B hadron or not
  bool isBHadron(int pdgid);

  // compute dR between two objects
  float deltaR(float eta1, float eta2, float phi1, float phi2);

  const xAOD::Jet* GetParentJet(const xAOD::Jet* Jet, std::string Keyname);

 int getTrackOrigin(const xAOD::TrackParticle *tmpTrk,
	            std::vector<const xAOD::TruthParticle*> tracksFromB,
		    std::vector<const xAOD::TruthParticle*> tracksFromC,
		    std::vector<const xAOD::TruthParticle*> tracksFromCc,
		    std::vector<const xAOD::TruthParticle*> tracksFromB1,
		    std::vector<const xAOD::TruthParticle*> tracksFromB2,
		    std::vector<const xAOD::TruthParticle*> tracksFromC1,
		    std::vector<const xAOD::TruthParticle*> tracksFromC2,
		    std::vector<const xAOD::TruthParticle*> tracksFromCNotFromB1,
		    std::vector<const xAOD::TruthParticle*> tracksFromCNotFromB2);
  
 int getLFCalibType(const xAOD::Jet  * inJet, const xAOD::Vertex* primVtx_ptr);
};

#endif //> !BTAGIBLANALYSIS_BTAGIBLANALYSISALG_H

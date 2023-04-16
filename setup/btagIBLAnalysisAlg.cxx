///////////////////////////////////////
// btagIBLAnalysisAlg.cxx
///////////////////////////////////////
// Author(s): marx@cern.ch
// Description: athena-based code to process xAODs
///////////////////////////////////////
#include "track_to_vertex_associators.hh"
#include "LifetimeSigning.hh"

#include <utility>
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ServiceHandle.h"

#include "btagIBLAnalysisAlg.h"
#include "GAFlavourLabel.h"
#include "ArbitraryJetBranches.hh"

#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "ParticleJetTools/JetFlavourInfo.h"
#include "xAODBTagging/SecVtxHelper.h"

#include "JetInterface/IJetSelector.h"
#include "JetCalibTools/IJetCalibrationTool.h"
#include "xAODEventShape/EventShape.h"

#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"
#include "TrackVertexAssociationTool/ITrackVertexAssociationTool.h"
//#include "TrkVertexFitterInterfaces/ITrackToVertexIPEstimator.h"
#include "PileupReweighting/PileupReweightingTool.h"

//#include "xAODTrigger/JetRoIContainer.h"
//#include "xAODTrigger/MuonRoIContainer.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "GoodRunsLists/IGoodRunsListSelectionTool.h"
#include "JetInterface/IJetUpdateJvt.h"

// some tracking mumbo jumbo
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

using xAOD::IParticle;

// DG: I'd like to make this NAN, but maybe someone is attached to
// the current convention.
const int MISSING_VALUE = -999;
// this is the key we use to keep track of how many primary vertices
// we found if it's missing we didn't run the BTagVertexAugmenter.
const std::string VX_COUNT_KEY = "BTaggingNumberOfPrimaryVertices";
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool particleInCollection( const xAOD::TrackParticle *trkPart, std::vector< ElementLink< xAOD::TrackParticleContainer > > trkColl ) {
  for (unsigned int iT = 0; iT < trkColl.size(); iT++) {
    if (trkPart == *(trkColl.at(iT))) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool xaodJetPtSorting(const xAOD::Jet *jet1, const xAOD::Jet *jet2) {
  return jet1->pt() > jet2->pt();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isFromWZ( const xAOD::TruthParticle* particle ) {
  if ( particle==0 ) return false;

  if ( fabs(particle->pdgId())!= 11 && fabs(particle->pdgId())!= 13) return false;

  const xAOD::TruthVertex* prodvtx = particle->prodVtx();
  if ( prodvtx==0 ) return false;

  if (  prodvtx->nIncomingParticles()==0 ) return false;

  if (  prodvtx->nIncomingParticles()>1 ) {
    int nCharge=0;
    int nNeutral=0;
    for(unsigned j = 0; j < prodvtx->nIncomingParticles(); j++){
      if ( fabs( prodvtx->incomingParticle(j)->pdgId() )==11 || fabs( prodvtx->incomingParticle(j)->pdgId() )==13 ) nCharge++;
      if ( fabs( prodvtx->incomingParticle(j)->pdgId() )==12 || fabs( prodvtx->incomingParticle(j)->pdgId() )==14 ) nNeutral++;
    }
    if ( nCharge>1 ) return true;
    if ( nCharge+nNeutral>1 ) return true;
    return false;
  }
  int absPDG=fabs(prodvtx->incomingParticle(0)->pdgId());
  if ( absPDG==15) return false;
  else if ( absPDG==24 || absPDG==23 ) return true;
  return isFromWZ( prodvtx->incomingParticle(0) );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
btagIBLAnalysisAlg::btagIBLAnalysisAlg( const std::string& name, ISvcLocator *pSvcLocator ) :
  AthHistogramAlgorithm(name, pSvcLocator),
  m_stream("BTAGSTREAM"),
  m_dumpCaloInfo(false),
  m_dumpTrackCovariance(false),
  m_dumpGATracks(false),
  m_doMSV(false),
  m_svx_collections(),
  m_bhadron_branches(),
  m_kshort_branches(),
  m_kshortreco_branches(),
  m_jetfitter_branches(),
  m_cluster_branches(),
  m_substructure_moment_branches(),
  m_subjet_collections(),
  m_track_branches(),
  m_track_cov_branches(),
  m_ga_track_branches(),
  m_ga_track_cov_branches(),
  m_arb_branches(0),
  m_jetCleaningTool("JetCleaningTool/JetCleaningTool", this),
  m_jetCalibrationTool(""),
  m_InDetTrackSelectorTool(""),
  m_CPTrackingLooseLabel(""),
  m_TightTrackVertexAssociationTool(""),
  m_tdt("Trig::TrigDecisionTool/TrigDecisionTool"),
  m_GRLSelectionTool("GoodRunsListSelectionTool/GoodRunsListSelectionTool", this),
  m_jvt("")
{
  m_triggerLogic="";
  declareProperty( "Stream", m_stream );

  declareProperty( "JetCleaningTool", m_jetCleaningTool );
  declareProperty( "JetCalibrationTool", m_jetCalibrationTool );

  declareProperty( "InDetTrackSelectionTool", m_InDetTrackSelectorTool );
  declareProperty( "CPTrackingLooseLabel", m_CPTrackingLooseLabel);
  declareProperty( "TrackVertexAssociationTool", m_TightTrackVertexAssociationTool );
  //declareProperty( "TrackToVertexIPEstimator", m_trackToVertexIPEstimator );
  declareProperty( "JVTtool", m_jvt );

  declareProperty( "EssentialInfo", m_essentialInfo =true );
  declareProperty( "ReduceInfo"   , m_reduceInfo=false );
  declareProperty( "SaveTrkAlgInfo"   , m_saveTrkAlgInfo=false );
  declareProperty( "SaveJetLFCalibType"   , m_saveJetLFCalibType=false );
  declareProperty( "SubjetInfo"   , m_subjetInfo=false );
  declareProperty( "DumpTrackCovariance"   , m_dumpTrackCovariance=false );
  declareProperty( "DumpGATracks"   , m_dumpGATracks=false );
  declareProperty( "Rel20", m_rel20 = false );
  declareProperty( "DoMSV", m_doMSV = false );
  declareProperty( "doSMT", m_SMT = false );
  declareProperty( "bHadronBranches", m_bHadronInfo = true );
  declareProperty( "kShortBranches", m_kShortInfo = false );
  declareProperty( "kShortRecoBranches", m_kShortRecoInfo = false );
  declareProperty( "CalibrateJets", m_calibrateJets = true );
  declareProperty( "CleanJets", m_cleanJets = true );
  declareProperty( "CleanParentJet", m_clean_parent_jet = false );
  declareProperty( "TrackAssociator",
                   m_track_associator = "BTagTrackToJetAssociator");

  declareProperty( "GRLname", m_GRLname = "" );
  declareProperty( "JetCollectionName", m_jetCollectionName = "AntiKt4LCTopoJets" );
  declareProperty( "JetPtCut", m_jetPtCut = 20.e3 );

  declareProperty( "TriggerLogic", m_triggerLogic );
  declareProperty( "TriggerDecisionTool", m_tdt);

  declareProperty( "DumpCaloInfo", m_dumpCaloInfo);
  declareProperty( "ArbitraryDoubleBranches", m_arb_double_names);
  declareProperty( "ArbitraryFloatVectorBranches", m_arb_float_vec_names);

  declareProperty( "subjetCollections", m_subjet_collections);
  declareProperty( "svxCollections", m_svx_collections);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
btagIBLAnalysisAlg::~btagIBLAnalysisAlg() {
  delete m_arb_branches;
  // FIXME: we're leaking memory with all the vectors, that we never
  // delete, but I suppose there are bigger issues with this code.
  for (auto coll_br: m_subjet_branches) {
    delete coll_br.second;
    coll_br.second = 0;
  }
  for (auto coll: m_svx_branches) {
    delete coll;
    coll = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
StatusCode btagIBLAnalysisAlg::initialize() {
  ATH_MSG_INFO ("Initializing " << name() << "...");

  if (m_essentialInfo) m_reduceInfo=true;

  // Register histograms
  //ATH_CHECK( book( TH1F("hist_Lxy_denom", "Lxy", 200, 0.0, 100.0) ) );

  // Register output tree
  ServiceHandle<ITHistSvc> histSvc("THistSvc",name());
  CHECK( histSvc.retrieve() );

  if (m_GRLname != "") ATH_CHECK(m_GRLSelectionTool.retrieve());

  ATH_MSG_INFO(m_jetCollectionName);

  tree = new TTree( ("bTag_" + m_jetCollectionName).c_str(),
		    ("bTag"  + m_jetCollectionName).c_str() );
  ATH_MSG_INFO ("VALERIO: registering tree in stream: " << m_stream);
  CHECK( histSvc->regTree("/" + m_stream + "/tree_" + m_jetCollectionName, tree) );

  // Retrieve the jet cleaning tool
  CHECK( m_jetCleaningTool.retrieve() );
  CHECK( m_jvt.retrieve() );

  // Retrieve the jet calibration tool
  m_jetCalibrationTool.setTypeAndName("JetCalibrationTool/BTagDumpAlg_" + m_jetCollectionName + "_JCalib");
  if (m_calibrateJets) CHECK( m_jetCalibrationTool.retrieve() );

  // retrieve other special tools
  if (m_InDetTrackSelectorTool.retrieve().isFailure())  {
    ATH_MSG_FATAL("#BTAG# Failed to retrieve tool " << m_InDetTrackSelectorTool);
    return StatusCode::FAILURE;
  }
  else {
    ATH_MSG_DEBUG("#BTAG# Retrieved tool " << m_InDetTrackSelectorTool);
  }
  if (!m_CPTrackingLooseLabel.empty()) {
    if (m_CPTrackingLooseLabel.retrieve().isFailure()) {
      ATH_MSG_FATAL("#BTAG# Failed to retrieve tool " << m_CPTrackingLooseLabel);
      return StatusCode::FAILURE;
    } else {
      ATH_MSG_DEBUG("#BTAG# Retrieved tool " << m_CPTrackingLooseLabel);
    }
  } else {
    ATH_MSG_WARNING("#BTAG# tracking cp loose tool not specified, "
                    "will not fill that branch");
  }
  if (m_TightTrackVertexAssociationTool.retrieve().isFailure())  {
    ATH_MSG_FATAL("#BTAG# Failed to retrieve tool " <<  m_TightTrackVertexAssociationTool);
    return StatusCode::FAILURE;
  }
  //if (m_trackToVertexIPEstimator.retrieve().isFailure()) {
  //  ATH_MSG_FATAL("#BTAG# Failed to retrieve tool " << m_trackToVertexIPEstimator);
  //  return StatusCode::FAILURE;
  //}
  //else {
  //  ATH_MSG_DEBUG("#BTAG# Retrieved tool " << m_trackToVertexIPEstimator);
  //}

  m_PUtool.setTypeAndName("CP::PileupReweightingTool/prw");
  CHECK( m_PUtool.retrieve() );

  ATH_CHECK(m_tdt.retrieve());

  m_jetfitter_branches.set_tree(*tree);

  m_bhadron_branches.set_tree(*tree, !m_reduceInfo && m_bHadronInfo );

  if(m_kShortInfo){  m_kshort_branches.set_tree(*tree); }
  if(m_kShortRecoInfo){  m_kshortreco_branches.set_tree(*tree); }

  // addition from Dan: create cluster branches
  if (m_dumpCaloInfo) {
    m_cluster_branches.set_tree(*tree);
    m_substructure_moment_branches.set_tree(*tree);
  }
  for (const auto& prefix_name_collection_name: m_subjet_collections) {
    const auto& prefix = prefix_name_collection_name.first;
    const auto& collection = prefix_name_collection_name.second;
    m_subjet_branches.emplace_back(collection, new SubjetBranches);
    m_subjet_branches.back().second->set_tree(*tree, prefix);
  }
  for (const auto& prefix_name_edm_name: m_svx_collections) {
    const auto& ntuple_prefix = prefix_name_edm_name.first;
    const auto& edm_name = prefix_name_edm_name.second;
    m_svx_branches.emplace_back(new SVBranches(edm_name));
    m_svx_branches.back()->set_tree(*tree, ntuple_prefix);
  }
  if (m_dumpTrackCovariance && !m_reduceInfo) {
    m_track_cov_branches.set_tree(*tree, "jet_trk_");
  }
  if (m_dumpGATracks) {
    m_ga_track_branches.set_tree(*tree, "jet_ga_trk_");
    if (m_dumpTrackCovariance) {
      m_ga_track_cov_branches.set_tree(*tree, "jet_ga_trk_");
    }
  }
  if (m_arb_double_names.size() + m_arb_float_vec_names.size() > 0) {
    m_arb_branches = new ArbitraryJetBranches(m_arb_double_names,
                                              m_arb_float_vec_names);
    m_arb_branches->set_tree(*tree, "jet_");
  }

  // Setup branches
  v_jet_pt = new std::vector<float>(); //v_jet_pt->reserve(15);
  v_jet_eta = new std::vector<float>(); //v_jet_eta->reserve(15);
  v_jet_phi = new std::vector<float>(); //v_jet_phi->reserve(15);
  v_jet_pt_orig =new std::vector<float>();
  v_jet_eta_orig = new std::vector<float>();
  v_jet_phi_orig = new std::vector<float>();
  v_jet_E_orig = new std::vector<float>();
  v_jet_sumtrkS_pt = new std::vector<float>();
  v_jet_sumtrkV_pt = new std::vector<float>();
  v_jet_sumtrkV_phi = new std::vector<float>();
  v_jet_sumtrkV_eta = new std::vector<float>();
  v_jet_sumtrk_ntrk = new std::vector<int>();
  v_jet_E = new std::vector<float>(); //v_jet_E->reserve(15);
  v_jet_m = new std::vector<float>(); //v_jet_m->reserve(15);
  v_jet_nConst     =new std::vector<int>();
  v_jet_truthflav = new std::vector<int>();
  v_jet_GhostL_q = new std::vector<int>();
  v_jet_GhostL_HadI = new std::vector<int>();
  v_jet_GhostL_HadF = new std::vector<int>();
  v_jet_LabDr_HadF = new std::vector<int>();
  v_jet_DoubleHadLabel = new std::vector<int>();
  v_jet_aliveAfterOR = new std::vector<int>();
  v_jet_aliveAfterORmu =new std::vector<int>();
  v_jet_truthMatch = new std::vector<int>();
  v_jet_isPU = new std::vector<int>();
  v_jet_isBadMedium = new std::vector<int>();
  v_jet_truthPt = new std::vector<float>();
  v_jet_dRiso = new std::vector<float>();
  v_jet_JVT = new std::vector<float>();
  v_jet_JVF = new std::vector<float>();
  v_jet_dRminToB = new std::vector<float>();
  v_jet_dRminToC = new std::vector<float>();
  v_jet_dRminToT = new std::vector<float>();

  v_jet_ip2d_pb = new std::vector<float>();
  v_jet_ip2d_pc = new std::vector<float>();
  v_jet_ip2d_pu = new std::vector<float>();
  v_jet_ip2d_llr = new std::vector<float>();
  v_jet_ip3d_pb = new std::vector<float>();
  v_jet_ip3d_pc = new std::vector<float>();
  v_jet_ip3d_pu = new std::vector<float>();
  v_jet_ip3d_llr = new std::vector<float>();
  v_jet_sv1_pb = new std::vector<float>();
  v_jet_sv1_pc = new std::vector<float>();
  v_jet_sv1_pu = new std::vector<float>();
  v_jet_sv1_llr = new std::vector<float>();

  v_jet_dl1_pb=new std::vector<float>();
  v_jet_dl1_pc=new std::vector<float>();
  v_jet_dl1_pu=new std::vector<float>();

  v_jet_dl1mu_pb=new std::vector<float>();
  v_jet_dl1mu_pc=new std::vector<float>();
  v_jet_dl1mu_pu=new std::vector<float>();

  v_jet_dl1rnn_pb=new std::vector<float>();
  v_jet_dl1rnn_pc=new std::vector<float>();
  v_jet_dl1rnn_pu=new std::vector<float>();

  v_jet_sv1ip3d = new std::vector<double>();
  v_jet_mv1 = new std::vector<double>();
  v_jet_mv1c = new std::vector<double>();
  v_jet_mv2c00 = new std::vector<double>();
  v_jet_mv2c10 = new std::vector<double>();
  v_jet_mv2c10mu = new std::vector<double>();
  v_jet_mv2c10rnn = new std::vector<double>();
  v_jet_mv2c20 = new std::vector<double>();
  v_jet_mv2c100 = new std::vector<double>();
  v_jet_mv2cl100 = new std::vector<double>();
  v_jet_mv2m_pu = new std::vector<double>();
  v_jet_mv2m_pc = new std::vector<double>();
  v_jet_mv2m_pb = new std::vector<double>();
  v_jet_mvb = new std::vector<double>();

  v_jet_mv2c10flip = new std::vector<float>();
  v_jet_mv2c20flip = new std::vector<float>();

  v_jet_multisvbb1 = new std::vector<double>();
  v_jet_multisvbb2 = new std::vector<double>();
  v_jet_msv_N2Tpair = new std::vector<int>();
  v_jet_msv_energyTrkInJet = new std::vector<float>();
  v_jet_msv_nvsec = new std::vector<int>();
  v_jet_msv_normdist = new std::vector<float>();
  v_jet_msv_vtx_cov0 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov1 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov2 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov3 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov4 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov5 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_mass = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_efrc = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_ntrk = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_pt = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_eta = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_phi = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_dls = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_x = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_y = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_z = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_chi = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_ndf = new std::vector<std::vector<float> >();

  v_jet_ExKtbb_Hbb_DoubleMV2c20 = new std::vector<double>();
  v_jet_ExKtbb_Hbb_SingleMV2c20 = new std::vector<double>();
  v_jet_ExKtbb_Hbb_MV2Only = new std::vector<double>();
  v_jet_ExKtbb_Hbb_MV2andJFDRSig = new std::vector<double>();
  v_jet_ExKtbb_Hbb_MV2andTopos = new std::vector<double>();

  v_jet_mv2c10flip = new std::vector<float>();
  v_jet_mv2c20flip = new std::vector<float>();

  v_jet_multisvbb1 = new std::vector<double>();
  v_jet_multisvbb2 = new std::vector<double>();
  v_jet_msv_N2Tpair = new std::vector<int>();
  v_jet_msv_energyTrkInJet = new std::vector<float>();
  v_jet_msv_nvsec = new std::vector<int>();
  v_jet_msv_normdist = new std::vector<float>();
  v_jet_msv_vtx_cov0 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov1 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov2 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov3 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov4 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_cov5 = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_mass = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_efrc = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_ntrk = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_pt = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_eta = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_phi = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_dls = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_x = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_y = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_z = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_chi = new std::vector<std::vector<float> >();
  v_jet_msv_vtx_ndf = new std::vector<std::vector<float> >();

  v_jet_ExKtbb_Hbb_DoubleMV2c20 = new std::vector<double>();
  v_jet_ExKtbb_Hbb_SingleMV2c20 = new std::vector<double>();
  v_jet_ExKtbb_Hbb_MV2Only = new std::vector<double>();
  v_jet_ExKtbb_Hbb_MV2andJFDRSig = new std::vector<double>();
  v_jet_ExKtbb_Hbb_MV2andTopos = new std::vector<double>();

  v_jet_LFCalibType = new std::vector<int>();
  v_jet_btag_ntrk = new std::vector<int>();
  v_jet_trk_dr = new std::vector<std::vector<float> >();
  v_jet_trk_assoc_msv = new std::vector<std::vector<int> >();   // mod nikola
  v_jet_trk_algo = new std::vector<std::vector<int> >();
  // v_jet_trk_orig = new std::vector<std::vector<int> >(); moved to BHadronBranches
  v_jet_trk_is_tracking_cp_loose = new std::vector<std::vector<int> >();

  v_jet_trk_vtx_X = new std::vector<std::vector<float> >();
  v_jet_trk_vtx_Y = new std::vector<std::vector<float> >();
  v_jet_trk_vtx_Z = new std::vector<std::vector<float> >();
  v_jet_trk_vtx_dx = new std::vector<std::vector<float> >();
  v_jet_trk_vtx_dy = new std::vector<std::vector<float> >();

  // v_jet_trk_d0sig = new std::vector<std::vector<float> >();
  // v_jet_trk_z0sig = new std::vector<std::vector<float> >();
  v_jet_trk_d0_truth = new std::vector<std::vector<float> >();
  v_jet_trk_z0_truth = new std::vector<std::vector<float> >();

  v_jet_trk_IP3D_grade = new std::vector<std::vector<int> >();
  v_jet_trk_IP3D_d0 = new std::vector<std::vector<float> >();
  v_jet_trk_IP3D_d02D = new std::vector<std::vector<float> >();
  v_jet_trk_IP3D_z0 = new std::vector<std::vector<float> >();
  v_jet_trk_IP3D_d0sig = new std::vector<std::vector<float> >();
  v_jet_trk_IP3D_z0sig = new std::vector<std::vector<float> >();
  v_jet_trk_IP3D_llr = new std::vector<std::vector<float> >();

  v_jet_trk_pdg_id = new std::vector<std::vector<int> >();
  v_jet_trk_barcode = new std::vector<std::vector<int> >();
  v_jet_trk_parent_pdgid = new std::vector<std::vector<int> >();

  // those are just quick accessors
  v_jet_sv1_ntrk = new std::vector<int>();
  v_jet_ip3d_ntrk = new std::vector<int>();
  v_jet_jf_ntrk = new std::vector<int>();

  // MVb variables
  v_jet_width  = new std::vector<float>();
  v_jet_n_trk_sigd0cut  = new std::vector<int>();
  v_jet_trk3_d0sig  = new std::vector<float>();
  v_jet_trk3_z0sig  = new std::vector<float>();
  v_jet_sv_scaled_efc  = new std::vector<float>();
  v_jet_jf_scaled_efc  = new std::vector<float>();

  // additions by andrea
  v_jet_mu_smt = new std::vector<double>();
  v_jet_mu_assJet_pt = new std::vector<float>();
  v_jet_mu_truthflav = new std::vector<float>();
  v_jet_mu_dR = new std::vector<float>();
  v_jet_mu_pTrel = new std::vector<float>();;
  v_jet_mu_qOverPratio = new std::vector<float>();
  v_jet_mu_mombalsignif = new std::vector<float>();
  v_jet_mu_scatneighsignif = new std::vector<float>();
  v_jet_mu_VtxTyp = new std::vector<float>();
  v_jet_mu_pt = new std::vector<float>();
  v_jet_mu_eta = new std::vector<float>();
  v_jet_mu_phi = new std::vector<float>();
  v_jet_mu_d0 = new std::vector<float>();
  v_jet_mu_z0 = new std::vector<float>();
  v_jet_mu_parent_pdgid = new std::vector<float>();
  v_jet_mu_ID_qOverP_var = new std::vector<float>();
  v_jet_mu_muonType = new std::vector<float>();
  // additions by nikola
  v_jet_mu_fatjet_nMu = new std::vector<int>();
  v_jet_mu_fatjet_pTmax_pT = new std::vector<float>();
  v_jet_mu_fatjet_pTmax_pTrel = new std::vector<float>();
  v_jet_mu_fatjet_pTmax_pTrelFrac = new std::vector<float>();

  tree->Branch("runnb", &runnumber);
  tree->Branch("eventnb", &eventnumber);
  tree->Branch("mcchan", &mcchannel);
  tree->Branch("mcwg", &mcweight);
  tree->Branch("lbn", &lbn);
  tree->Branch("coreFlag", &coreFlag);
  tree->Branch("larError", &larError);
  tree->Branch("tileError", &tileError);
  tree->Branch("nPV", &npv);
  tree->Branch("avgmu", &mu);
  tree->Branch("actmu",&Act_mu);
  tree->Branch("PVx", &PV_x);
  tree->Branch("PVy", &PV_y);
  tree->Branch("PVz", &PV_z);
  tree->Branch("truth_PVx", &truth_PV_x);
  tree->Branch("truth_PVy", &truth_PV_y);
  tree->Branch("truth_PVz", &truth_PV_z);
  tree->Branch("truth_LeadJet_pt", &truth_LeadJet_pt);

  tree->Branch("njets", &njets);
  if (!m_essentialInfo) tree->Branch("nbjets", &nbjets);
  if (!m_essentialInfo) tree->Branch("nbjets_q", &nbjets_q);
  if (!m_essentialInfo) tree->Branch("nbjets_HadI", &nbjets_HadI);
  if (!m_essentialInfo) tree->Branch("nbjets_HadF", &nbjets_HadF);
  tree->Branch("jet_pt", &v_jet_pt);
  tree->Branch("jet_eta", &v_jet_eta);
  tree->Branch("jet_phi", &v_jet_phi);
  tree->Branch("jet_pt_orig", &v_jet_pt_orig);
  tree->Branch("jet_eta_orig", &v_jet_eta_orig);
  if (!m_essentialInfo) tree->Branch("jet_phi_orig", &v_jet_phi_orig);
  if (!m_essentialInfo) tree->Branch("jet_E_orig", &v_jet_E_orig);
  if (!m_essentialInfo) tree->Branch("jet_sumtrkS_pt", &v_jet_sumtrkS_pt);
  tree->Branch("jet_sumtrkV_pt", &v_jet_sumtrkV_pt);
  tree->Branch("jet_sumtrkV_eta", &v_jet_sumtrkV_eta);
  if (!m_essentialInfo) tree->Branch("jet_sumtrkV_phi", &v_jet_sumtrkV_phi);
  tree->Branch("jet_sumtrk_ntrk", &v_jet_sumtrk_ntrk);
  tree->Branch("jet_E", &v_jet_E);
  if (!m_essentialInfo) tree->Branch("jet_m", &v_jet_m);
  tree->Branch("jet_nConst",&v_jet_nConst);
  if (!m_essentialInfo) tree->Branch("jet_truthflav", &v_jet_truthflav);
  if (!m_essentialInfo) tree->Branch("jet_GhostL_q", &v_jet_GhostL_q);
  if (!m_essentialInfo) tree->Branch("jet_GhostL_HadI", &v_jet_GhostL_HadI);
  if (!m_essentialInfo) tree->Branch("jet_GhostL_HadF", &v_jet_GhostL_HadF);
  tree->Branch("jet_LabDr_HadF", &v_jet_LabDr_HadF);
  tree->Branch("jet_DoubleHadLabel", &v_jet_DoubleHadLabel);
  tree->Branch("jet_aliveAfterOR", &v_jet_aliveAfterOR);
  tree->Branch("jet_aliveAfterORmu" ,&v_jet_aliveAfterORmu);
  tree->Branch("jet_truthMatch", &v_jet_truthMatch);
  tree->Branch("jet_isPU", &v_jet_isPU);
  tree->Branch("jet_isBadMedium", &v_jet_isBadMedium);
  tree->Branch("jet_truthPt", &v_jet_truthPt);
  tree->Branch("jet_dRiso", &v_jet_dRiso);
  tree->Branch("jet_JVT", &v_jet_JVT);
  if (!m_essentialInfo) tree->Branch("jet_JVF", &v_jet_JVF);
  tree->Branch("jet_dRminToB", &v_jet_dRminToB);
  tree->Branch("jet_dRminToC", &v_jet_dRminToC);
  tree->Branch("jet_dRminToT", &v_jet_dRminToT);

  if (!m_essentialInfo) tree->Branch("jet_ip2d_pb", &v_jet_ip2d_pb);
  if (!m_essentialInfo) tree->Branch("jet_ip2d_pc", &v_jet_ip2d_pc);
  if (!m_essentialInfo) tree->Branch("jet_ip2d_pu", &v_jet_ip2d_pu);
  tree->Branch("jet_ip2d_llr", &v_jet_ip2d_llr);

  tree->Branch("jet_ip3d_pb", &v_jet_ip3d_pb);
  tree->Branch("jet_ip3d_pc", &v_jet_ip3d_pc);
  tree->Branch("jet_ip3d_pu", &v_jet_ip3d_pu);
  tree->Branch("jet_ip3d_llr", &v_jet_ip3d_llr);

  if (!m_essentialInfo) tree->Branch("jet_sv1_pb", &v_jet_sv1_pb);
  if (!m_essentialInfo) tree->Branch("jet_sv1_pc", &v_jet_sv1_pc);
  if (!m_essentialInfo) tree->Branch("jet_sv1_pu", &v_jet_sv1_pu);
  tree->Branch("jet_sv1_llr", &v_jet_sv1_llr);

  tree->Branch("jet_dl1_pb",&v_jet_dl1_pb);
  tree->Branch("jet_dl1_pc",&v_jet_dl1_pc);
  tree->Branch("jet_dl1_pu",&v_jet_dl1_pu);

  tree->Branch("jet_dl1mu_pb",&v_jet_dl1mu_pb);
  tree->Branch("jet_dl1mu_pc",&v_jet_dl1mu_pc);
  tree->Branch("jet_dl1mu_pu",&v_jet_dl1mu_pu);

  tree->Branch("jet_dl1rnn_pb",&v_jet_dl1rnn_pb);
  tree->Branch("jet_dl1rnn_pc",&v_jet_dl1rnn_pc);
  tree->Branch("jet_dl1rnn_pu",&v_jet_dl1rnn_pu);

  if (!m_essentialInfo) tree->Branch("jet_sv1ip3d", &v_jet_sv1ip3d);
  if (!m_essentialInfo) tree->Branch("jet_mv1", &v_jet_mv1);
  if (!m_essentialInfo) tree->Branch("jet_mv1c", &v_jet_mv1c);
  tree->Branch("jet_mv2c00", &v_jet_mv2c00);
  tree->Branch("jet_mv2c10", &v_jet_mv2c10);
  tree->Branch("jet_mv2c10mu", &v_jet_mv2c10mu);
  tree->Branch("jet_mv2c10rnn", &v_jet_mv2c10rnn);
  tree->Branch("jet_mv2c20", &v_jet_mv2c20);
  tree->Branch("jet_mv2c100", &v_jet_mv2c100);
  tree->Branch("jet_mv2cl100", &v_jet_mv2cl100);
  tree->Branch("jet_mv2m_pu", &v_jet_mv2m_pu);
  tree->Branch("jet_mv2m_pc", &v_jet_mv2m_pc);
  tree->Branch("jet_mv2m_pb", &v_jet_mv2m_pb);
  tree->Branch("jet_mvb", &v_jet_mvb);

  tree->Branch("jet_mv2c10Flip", &v_jet_mv2c10flip);
  tree->Branch("jet_mv2c20Flip", &v_jet_mv2c20flip);

  if (!m_reduceInfo && m_doMSV) {
    tree->Branch("jet_multisvbb1", &v_jet_multisvbb1);
    tree->Branch("jet_multisvbb2", &v_jet_multisvbb2);
    tree->Branch("jet_msv_N2Tpair", &v_jet_msv_N2Tpair);
    tree->Branch("jet_msv_energyTrkInJet", &v_jet_msv_energyTrkInJet);
    tree->Branch("jet_msv_nvsec", &v_jet_msv_nvsec);
    tree->Branch("jet_msv_normdist", &v_jet_msv_normdist);
    tree->Branch("jet_msv_vtx_cov0", &v_jet_msv_vtx_cov0);
    tree->Branch("jet_msv_vtx_cov1", &v_jet_msv_vtx_cov1);
    tree->Branch("jet_msv_vtx_cov2", &v_jet_msv_vtx_cov2);
    tree->Branch("jet_msv_vtx_cov3", &v_jet_msv_vtx_cov3);
    tree->Branch("jet_msv_vtx_cov4", &v_jet_msv_vtx_cov4);
    tree->Branch("jet_msv_vtx_cov5", &v_jet_msv_vtx_cov5);
    tree->Branch("jet_msv_vtx_mass", &v_jet_msv_vtx_mass);
    tree->Branch("jet_msv_vtx_efrc", &v_jet_msv_vtx_efrc);
    tree->Branch("jet_msv_vtx_ntrk", &v_jet_msv_vtx_ntrk);
    tree->Branch("jet_msv_vtx_pt", &v_jet_msv_vtx_pt);
    tree->Branch("jet_msv_vtx_eta", &v_jet_msv_vtx_eta);
    tree->Branch("jet_msv_vtx_phi", &v_jet_msv_vtx_phi);
    tree->Branch("jet_msv_vtx_dls", &v_jet_msv_vtx_dls);
    tree->Branch("jet_msv_vtx_x", &v_jet_msv_vtx_x);
    tree->Branch("jet_msv_vtx_y", &v_jet_msv_vtx_y);
    tree->Branch("jet_msv_vtx_z", &v_jet_msv_vtx_z);
    tree->Branch("jet_msv_vtx_chi", &v_jet_msv_vtx_chi);
    tree->Branch("jet_msv_vtx_ndf", &v_jet_msv_vtx_ndf);
  }

  if (!m_essentialInfo) tree->Branch("jet_ExKtbb_Hbb_DoubleMV2c20", &v_jet_ExKtbb_Hbb_DoubleMV2c20);
  if (!m_essentialInfo) tree->Branch("jet_ExKtbb_Hbb_SingleMV2c20", &v_jet_ExKtbb_Hbb_SingleMV2c20);
  if (!m_essentialInfo) tree->Branch("jet_ExKtbb_Hbb_MV2Only", &v_jet_ExKtbb_Hbb_MV2Only);
  if (!m_essentialInfo) tree->Branch("jet_ExKtbb_Hbb_MV2andJFDRSig", &v_jet_ExKtbb_Hbb_MV2andJFDRSig);
  if (!m_essentialInfo) tree->Branch("jet_ExKtbb_Hbb_MV2andTopos", &v_jet_ExKtbb_Hbb_MV2andTopos);

  if( m_saveJetLFCalibType) tree->Branch("jet_LFCalibType", &v_jet_LFCalibType);
  tree->Branch("jet_btag_ntrk", &v_jet_btag_ntrk);

  if (!m_reduceInfo) {
    m_track_branches.set_tree(*tree, "jet_trk_");
    tree->Branch("jet_trk_dr", &v_jet_trk_dr);
    tree->Branch("jet_trk_assoc_msv", &v_jet_trk_assoc_msv);    // mod nikola
    tree->Branch("jet_trk_algo", &v_jet_trk_algo);
    // tree->Branch("jet_trk_orig", &v_jet_trk_orig); moved to BHadronBranches
    if (!m_CPTrackingLooseLabel.empty()) {
      tree->Branch("jet_trk_is_tracking_cp_loose",
                   &v_jet_trk_is_tracking_cp_loose);
    }

    tree->Branch("jet_trk_vtx_X", &v_jet_trk_vtx_X);
    tree->Branch("jet_trk_vtx_Y", &v_jet_trk_vtx_Y);
    tree->Branch("jet_trk_vtx_Z", &v_jet_trk_vtx_Z);

    tree->Branch("jet_trk_d0_truth", &v_jet_trk_d0_truth);
    tree->Branch("jet_trk_z0_truth", &v_jet_trk_z0_truth);
    tree->Branch("jet_trk_ip3d_grade", &v_jet_trk_IP3D_grade);
    tree->Branch("jet_trk_ip3d_d0", &v_jet_trk_IP3D_d0);
    tree->Branch("jet_trk_ip3d_d02D", &v_jet_trk_IP3D_d02D);
    tree->Branch("jet_trk_ip3d_z0", &v_jet_trk_IP3D_z0);
    tree->Branch("jet_trk_ip3d_d0sig", &v_jet_trk_IP3D_d0sig);
    tree->Branch("jet_trk_ip3d_z0sig", &v_jet_trk_IP3D_z0sig);

    tree->Branch("jet_trk_ip3d_llr", &v_jet_trk_IP3D_llr);

    tree->Branch("jet_trk_pdg_id", &v_jet_trk_pdg_id);
    tree->Branch("jet_trk_barcode", &v_jet_trk_barcode);
    tree->Branch("jet_trk_parent_pdgid", &v_jet_trk_parent_pdgid);
  }

  if (!m_essentialInfo) tree->Branch("jet_sv1_ntrk",&v_jet_sv1_ntrk);
  if (!m_essentialInfo) tree->Branch("jet_ip3d_ntrk",&v_jet_ip3d_ntrk);
  if (!m_essentialInfo) tree->Branch("jet_jf_ntrk",&v_jet_jf_ntrk);

  // MVb variables
  if (!m_essentialInfo) tree->Branch("jet_width", &v_jet_width);
  if (!m_essentialInfo) tree->Branch("jet_n_trk_sigd0cut", &v_jet_n_trk_sigd0cut);
  if (!m_essentialInfo) tree->Branch("jet_trk3_d0sig", &v_jet_trk3_d0sig);
  if (!m_essentialInfo) tree->Branch("jet_trk3_z0sig", &v_jet_trk3_z0sig);
  if (!m_essentialInfo) tree->Branch("jet_sv_scaled_efc", &v_jet_sv_scaled_efc);
  if (!m_essentialInfo) tree->Branch("jet_jf_scaled_efc", &v_jet_jf_scaled_efc);

  // additions by andrea
  if (m_SMT) {
    tree->Branch("jet_mu_smt", &v_jet_mu_smt);
    tree->Branch("jet_mu_assJet_pt", &v_jet_mu_assJet_pt);
    tree->Branch("jet_mu_truthflav", &v_jet_mu_truthflav);
    tree->Branch("jet_mu_dR", &v_jet_mu_dR);
    tree->Branch("jet_mu_pTrel", &v_jet_mu_pTrel);
    tree->Branch("jet_mu_qOverPratio", &v_jet_mu_qOverPratio);
    tree->Branch("jet_mu_mombalsignif", &v_jet_mu_mombalsignif);
    tree->Branch("jet_mu_scatneighsignif", &v_jet_mu_scatneighsignif);
    tree->Branch("jet_mu_VtxTyp", &v_jet_mu_VtxTyp);
    tree->Branch("jet_mu_pt", &v_jet_mu_pt);
    tree->Branch("jet_mu_eta", &v_jet_mu_eta);
    tree->Branch("jet_mu_phi", &v_jet_mu_phi);
    tree->Branch("jet_mu_d0", &v_jet_mu_d0);
    tree->Branch("jet_mu_z0", &v_jet_mu_z0);
    tree->Branch("jet_mu_parent_pdgid", &v_jet_mu_parent_pdgid);
    tree->Branch("jet_mu_ID_qOverP_var", &v_jet_mu_ID_qOverP_var);
    tree->Branch("jet_mu_muonType", &v_jet_mu_muonType);
    // additions by nikola
    tree->Branch("jet_mu_fatjet_nMu", &v_jet_mu_fatjet_nMu);
    tree->Branch("jet_mu_fatjet_pTmax_pT", &v_jet_mu_fatjet_pTmax_pT);
    tree->Branch("jet_mu_fatjet_pTmax_pTrel", &v_jet_mu_fatjet_pTmax_pTrel);
    tree->Branch("jet_mu_fatjet_pTmax_pTrelFrac", &v_jet_mu_fatjet_pTmax_pTrelFrac);
  }

  clearvectors();
  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////
StatusCode btagIBLAnalysisAlg::finalize() {
  ATH_MSG_INFO ("Finalizing " << name() << "...");

  // Write tree into file
  tree->Write();

  // Clean up
  CHECK( m_jetCleaningTool.release() );
  CHECK( m_jetCalibrationTool.release() );

  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////

namespace {
  typedef std::vector<const xAOD::TrackParticle*> Tracks;
  typedef std::vector<const xAOD::IParticle*> Particles;
}

StatusCode btagIBLAnalysisAlg::execute() {
  typedef ElementLink<xAOD::TrackParticleContainer> TrackLink;
  typedef std::vector<TrackLink> TrackLinks;

  ATH_MSG_DEBUG ("Executing " << name() << "...");

  std::string triggerLogic = "HLT_j[0-9]+|L1_MBTS_1_1|L1_RD0_FILLED";

  if (v_L1triggerNames.size() == 0 && m_tdt != 0) {
    ///ATH_MSG_INFO ("Setting up trigger informations " << name() << "...");
    auto chainGroup = m_tdt->getChainGroup(m_triggerLogic);
    v_L1trigger = new bool[chainGroup->getListOfTriggers().size()];
    int count = -1;
    for (auto & trig : chainGroup->getListOfTriggers()) {
      count++;
      // std::cout << trig << std::endl;
      v_L1trigger[count] = 0;
      v_L1triggerNames.push_back(trig);
      tree->Branch(trig.c_str(), &(v_L1trigger[count]));
    }
  }
  clearvectors();

  //-------------------------
  // Event information
  //---------------------------
  const xAOD::EventInfo* eventInfo = 0;
  CHECK( evtStore()->retrieve(eventInfo) );

  // check if data or MC
  bool isData = !eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION);

  if (m_GRLname != "") {
    bool eventPasses = m_GRLSelectionTool->passRunLB(*eventInfo);
    if (!eventPasses) {
      ATH_MSG_DEBUG( " .... Event not passing GRL!!!");
      return StatusCode::SUCCESS;
    }
  }

  runnumber = eventInfo->runNumber();
  eventnumber = eventInfo->eventNumber();


  mcchannel = ( isData ? 0 : eventInfo->mcChannelNumber() );
  mcweight = ( isData ? 1 : eventInfo->mcEventWeight() );
  mu = eventInfo->averageInteractionsPerCrossing();
  Act_mu      = eventInfo->actualInteractionsPerCrossing();
  lbn = ( !isData ? 0 : eventInfo->lumiBlock() );

  larError = eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error;
  tileError = eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error;
  coreFlag = eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18);

  float tmpMu = m_PUtool->getCorrectedMu( *eventInfo );
  // std::cout << " origMu: " << mu << " newValue: " <<  tmpMu << std::endl;
  if (isData) mu = tmpMu;
  //  }

  // primary vertex
  const xAOD::VertexContainer *vertices = 0;
  CHECK( evtStore()->retrieve(vertices, "PrimaryVertices") );
  int* npv_p = 0;
  StatusCode vx_count_status = evtStore()->retrieve(npv_p, VX_COUNT_KEY);
  if (vx_count_status.isFailure()) {
    ATH_MSG_FATAL("could not find " + VX_COUNT_KEY + " in file");
    return StatusCode::FAILURE;
  }
  npv = *npv_p;
  if (npv < 1) {
    ATH_MSG_WARNING( ".... rejecting the event due to missing PV!!!!");
    return StatusCode::SUCCESS;
  }

  // get the vertex index (stored in BTaggingVertexAugmenter)
  int* indexPV_ptr = 0;
  CHECK(evtStore()->retrieve(indexPV_ptr, "BTaggingVertexIndex"));
  size_t indexPV = *indexPV_ptr;
  const xAOD::Vertex *myVertex = vertices->at(indexPV); // the (reco?) primary vertex
  PV_x = myVertex->x(); // VALERIO !!!!!!!!
  PV_y = myVertex->y(); // VALERIO !!!!!!!!
  PV_z = myVertex->z();

  //---------------------------------------------------------------------------------------------
  // Trigger part
  //---------------------------------------------------------------------------------------------
  if (m_tdt != 0) {
    bool passTrigger = false;
    for (unsigned int it = 0; it < v_L1triggerNames.size(); it++) {
      v_L1trigger[it] = m_tdt->isPassed(v_L1triggerNames[it]);
      if (v_L1trigger[it]) passTrigger = true;
    }

    if (isData) {
      if (!passTrigger) {
        return StatusCode::SUCCESS;
      }
    }
  }

  if (!isData) {
    const xAOD::TruthVertexContainer *truthVertices = 0;
    CHECK( evtStore()->retrieve(truthVertices, "TruthVertices") );
    for ( const auto* truthV : *truthVertices ) {
      // record the truth primary vertex position
      truth_PV_x = truthV->x();
      truth_PV_y = truthV->y();
      truth_PV_z = truthV->z();
      break;
    }
  }

  //---------------------------
  // Truth stuff
  //---------------------------
  const xAOD::TruthEventContainer *xTruthEventContainer = NULL;
  std::string truthevt = "TruthEvent";
  if (m_rel20) truthevt = "TruthEvents";
  if (!isData) CHECK( evtStore()->retrieve(xTruthEventContainer, truthevt) );

  std::vector<const xAOD::TruthParticle* > m_partonB;
  std::vector<const xAOD::TruthParticle* > m_partonC;
  std::vector<const xAOD::TruthParticle* > m_partonT;
  std::vector<TLorentzVector> truth_electrons;
  std::vector<TLorentzVector> truth_muons;

  if (!isData) {
    // select truth electrons for electron-jet overlap removal
    for ( const auto* truth : *xTruthEventContainer ) {
      const xAOD::TruthVertex *newTruthVertex = truth->signalProcessVertex();
      if (newTruthVertex != 0) {
        // std::cout << " mah: " << newTruthVertex << std::endl;
        // rewrite the truth primary vertex position
	truth_PV_x = newTruthVertex->x();
	truth_PV_y = newTruthVertex->y();
	truth_PV_z = newTruthVertex->z();
      }

      // loop over truth particles in the truth event container
      for(unsigned int i = 0; i < truth->nTruthParticles(); i++) {
	const xAOD::TruthParticle *particle = truth->truthParticle(i);

  if(!particle){ continue; } // June6,2017, in rel21 derivations some truth particles have been slimmed, and there is a null pointer here.

  // VALERIO !!!!!!!!
	if (particle->pt() > 3e3) {
	  if (fabs(particle->pdgId()) == 15) m_partonT.push_back(particle);
	  if (particle->isCharmHadron()) m_partonC.push_back(particle);
	  if (particle->isBottomHadron()) m_partonB.push_back(particle);
	}
	if (particle->pt() < 10e3)     continue; //was 15
	if (particle->status() != 1)   continue;
	if (particle->barcode() > 2e5) continue;

	if ( fabs(particle->pdgId()) != 11 && fabs(particle->pdgId()) != 13) continue;

	// see if this electron is coming from a W boson decay
	bool isfromW = isFromWZ( particle );

	if(!isfromW) continue;
	TLorentzVector telec;
	telec.SetPtEtaPhiM(particle->pt(), particle->eta(), particle->phi(), particle->m());

	if (fabs(particle->pdgId()) == 11) truth_electrons.push_back(telec);
	if (fabs(particle->pdgId()) == 13) truth_muons.push_back(telec);
      }
    }
  }

  //---------------------------
  // Jets
  //---------------------------
  const xAOD::JetContainer *jets = 0;
  CHECK( evtStore()->retrieve(jets, m_jetCollectionName) );
  const xAOD::JetContainer *truthjets = 0;
  if (!isData) {
    CHECK( evtStore()->retrieve( truthjets, "AntiKt4TruthJets") );
    truth_LeadJet_pt=-1;
    if (truthjets->size()!=0) truth_LeadJet_pt=(truthjets->at(0))->pt();
  }
  else truthjets = new xAOD::JetContainer();
  njets = 0;
  nbjets = 0;
  nbjets_q = 0;
  nbjets_HadI = 0;
  nbjets_HadF = 0;

  // keep in mind, jets will not have tracks if their eta is greater than 2.5
  float etaCut = 3.0;
  if (isData) etaCut = 2.5;

  // Loop over jets, apply calibration and only keep the ones with pT > XX
  std::vector<const xAOD::Jet*> selJets; selJets.reserve(10);

  bool badCleaning = false;
  for (const auto* jet : *jets) {
    xAOD::Jet *newjet = 0;

    if (m_calibrateJets) m_jetCalibrationTool->calibratedCopy(*jet, newjet);
    else newjet = new xAOD::Jet(*jet);

    if (newjet->pt() < m_jetPtCut) {
      delete newjet;
      continue;
    }
    if ( fabs(newjet->eta()) > etaCut ) {
      delete newjet;
      continue;
    }
    v_jet_pt_orig->push_back(jet->pt());
    v_jet_eta_orig->push_back(jet->eta());
    v_jet_phi_orig->push_back(jet->phi());
    v_jet_E_orig ->push_back(jet->e());
    selJets.push_back(newjet);
  }
  if (badCleaning) {
    for (unsigned int j = 0; j < selJets.size(); j++) {
      delete selJets.at(j);
    }
    return StatusCode::SUCCESS;
  }
  std::vector<const xAOD::Jet*> selJetsTMP = selJets;
  std::sort(selJets.begin(), selJets.end(), xaodJetPtSorting);
  std::vector<float> tmpPt = std::vector<float>(*v_jet_pt_orig);
  std::vector<float> tmpEta = std::vector<float>(*v_jet_eta_orig);
  std::vector<float> tmpPhi = std::vector<float>(*v_jet_phi_orig);
  std::vector<float> tmpE = std::vector<float>(*v_jet_E_orig);
  v_jet_pt_orig->clear();
  v_jet_eta_orig->clear();
  v_jet_phi_orig->clear();
  v_jet_E_orig->clear();
  for (unsigned int j = 0; j < selJets.size(); j++) {
    const xAOD::Jet *jet = selJets.at(j);
    int oIndex = -1;
    for (unsigned int j2 = 0; j2 < selJets.size(); j2++) {
      if (selJetsTMP.at(j2)->pt() == jet->pt()) {
	oIndex = j2;
	break;
      }
    }
    v_jet_pt_orig->push_back(tmpPt[oIndex]);
    v_jet_eta_orig->push_back(tmpEta[oIndex]);
    v_jet_phi_orig->push_back(tmpPhi[oIndex]);
    v_jet_E_orig->push_back(tmpE[oIndex]);
  }

  njets = selJets.size();
  ATH_MSG_DEBUG( "Total number of jets is: "<< njets );
  uint8_t getInt(0);   // for accessing summary information

  /////////////////////////////////////////////////////////////////////////
  // MAIN JET LOOP
  // Now run over the selected jets and do whatever else needs doing
  for (unsigned int j = 0; j < selJets.size(); j++) {
    const xAOD::Jet *jet = selJets.at(j);

    // protection against not properly filled objects ... IT SHOULD NEVER HAPPEN THOUGH
    try {
      const std::vector<ElementLink<xAOD::VertexContainer > > TMPvertices = jet->btagging()->auxdata<std::vector<ElementLink<xAOD::VertexContainer > > >("SV1_vertices");
    } catch(...) {
      std::vector< ElementLink< xAOD::TrackParticleContainer > > TMPassocTracks = jet->btagging()->auxdata<std::vector<ElementLink<xAOD::TrackParticleContainer> > >("BTagTrackToJetAssociator");
      ATH_MSG_WARNING(" THIS JET: " << jet->pt() << " eta: " << jet->eta() << " IS NOT FILLED PROPERLY!!!!!!!!!!!!! and contains: " << TMPassocTracks.size() << " tracks  .... PLEASE CHECK");
      njets-=1;
      continue;
    }

    //////////////////////////////////////////////////////////////////
    // flagging jets that overlap with electron
    bool iseljetoverlap = false;
    for(unsigned int i = 0; i < truth_electrons.size(); i++) {
      float dr = deltaR(jet->eta(), truth_electrons.at(i).Eta(), jet->phi(), truth_electrons.at(i).Phi());
      if (dr < 0.3) iseljetoverlap = true;
    }
    v_jet_aliveAfterOR->push_back( !iseljetoverlap );

    iseljetoverlap = false;
    for(unsigned int i= 0; i < truth_muons.size(); i++){
      float dr =deltaR(jet->eta(), truth_muons.at(i).Eta(),jet->phi(), truth_muons.at(i).Phi());
      if(dr < 0.3) iseljetoverlap = true;
    }
    v_jet_aliveAfterORmu->push_back( !iseljetoverlap );

    //////////////////////////////////////////////////////////////////
    // simple kinematic quantities + PUinfo + JVF/JVT
    v_jet_pt->push_back(jet->pt());
    v_jet_eta->push_back(jet->eta());
    v_jet_phi->push_back(jet->phi());
    v_jet_E->push_back(jet->e());
    v_jet_m->push_back(jet->m());
    v_jet_nConst->push_back( jet->numConstituents() );
    float dRiso = 10;
    for (unsigned int j2 = 0; j2 < selJets.size(); j2++) {
      if (j2 == j) continue;
      const xAOD::Jet *jet2 = selJets.at(j2);
      float dr = deltaR(jet->eta(), jet2->eta(), jet->phi(), jet2->phi());
      if (dr < dRiso) dRiso = dr;
    }
    v_jet_dRiso->push_back(dRiso);

    // new way to save b-hadron quantities ...
    m_bhadron_branches.fill(*jet, m_jetCollectionName);

    // KShort info
    if(m_kShortInfo){
      if (!isData){
        const xAOD::TruthEventContainer *xTruthEventContainer = NULL;
	std::string truthevt = "TruthEvent";
        if (m_rel20) truthevt = "TruthEvents";
        CHECK( evtStore()->retrieve(xTruthEventContainer, truthevt) );
        m_kshort_branches.fill(*jet, *xTruthEventContainer, PV_x, PV_y, PV_z);
      }
    }

    if(m_kShortRecoInfo){
      const xAOD::VertexContainer* v0s = NULL;
      CHECK( evtStore()->retrieve(v0s, "V0Candidates") );
      m_kshortreco_branches.fill(*jet, *v0s, PV_x, PV_y, PV_z);
    }

    // addition from Dan: fill clusters
    if (m_dumpCaloInfo) {
      m_cluster_branches.fill(jet->getConstituents());
      m_substructure_moment_branches.fill(*jet);
    }
    // and fill arbitrary branches
    if (m_arb_branches) m_arb_branches->fill(*jet);

    // matching reco jets to truth jets and recording the truth jet pT
    // picking the highest pT truth jet (with pT > 7GeV) that satisfies dR < 0.3
    // N.B. this assumes that truth jets are pT ordered
    int matchedPt = 0;
    float dRmatch = 100;
    for (const auto* tjet : *truthjets) {
      if (tjet->pt() < 10e3) continue;
      float dr = deltaR(jet->eta(), tjet->eta(), jet->phi(), tjet->phi());

      if (dr < 0.3) {
	dRmatch = dr;
	matchedPt = tjet->pt();
	break;
      }
    }
    if (dRmatch < 0.3) {
      v_jet_truthMatch->push_back(1);
      v_jet_truthPt->push_back(matchedPt);
    }
    else {
      v_jet_truthMatch->push_back(0);
      v_jet_truthPt->push_back(0);
    }

    // flagging as PU jet: no HS jets within 0.6
    bool truthFree = true;
    for (const auto* tjet : *truthjets) {
      if (tjet->pt() < 4e3) continue;
      float dr = deltaR(jet->eta(), tjet->eta(), jet->phi(), tjet->phi());

      if (dr < 0.6) {
	truthFree = false;
	break;
      }
    }
    v_jet_isPU->push_back(truthFree);

    // clean jets
    if (m_cleanJets) {
      const xAOD::Jet *jet_to_clean = jet;
      // additions by nikola
      if (m_clean_parent_jet) {
        jet_to_clean = GetParentJet(jet, "Parent");
      }
      v_jet_isBadMedium->push_back(!m_jetCleaningTool->keep(*jet_to_clean));
    }
    else v_jet_isBadMedium->push_back(0);

    // JVF
    float jvfV = 0;
    std::vector<float> testjvf;
    bool jetHasJVF = jet->getAttribute<std::vector<float> >(xAOD::JetAttribute::JVF, testjvf);
    if (jetHasJVF && testjvf.size() > indexPV) jvfV = testjvf.at(indexPV);

    v_jet_JVF->push_back(jvfV);
    // JVT
    float jvtV = 0;
    try {
      jvtV = jet->auxdata<float>("Jvt");
      float tmpJVT = jvtV;
      jvtV = m_jvt->updateJvt(*jet);
      if (tmpJVT != jvtV)  ATH_MSG_DEBUG(" initial: " << tmpJVT << " |  final: " << jvtV );
    } catch (...) {
      ATH_MSG_WARNING(" something went wrong with the JVT recalculation .... please investigate");
    };
    v_jet_JVT->push_back(jvtV);

    ///////////////////////////////////////////////////////////////////
    // all the possible labeling
    int thisJetTruthLabel = -1;
    if (m_rel20) thisJetTruthLabel = jetFlavourLabel(jet, xAOD::ConeFinalParton);
    else jet->getAttribute("TruthLabelID", thisJetTruthLabel);

    v_jet_truthflav->push_back(thisJetTruthLabel);
    if (thisJetTruthLabel == 5) nbjets++;

    int tmpLabel = GAFinalPartonFlavourLabel(jet);
    v_jet_GhostL_q->push_back(tmpLabel);
    if (tmpLabel == 5) nbjets_q++;

    tmpLabel = GAInitialHadronFlavourLabel(jet);
    v_jet_GhostL_HadI->push_back(tmpLabel);
    if (tmpLabel == 5) nbjets_HadI++;

    tmpLabel = GAFinalHadronFlavourLabel(jet);
    v_jet_GhostL_HadF->push_back(tmpLabel);
    if (tmpLabel == 5) nbjets_HadF++;

    tmpLabel = -1;
    try {
      jet->getAttribute("HadronConeExclTruthLabelID", tmpLabel);
    } catch(...) {};
    v_jet_LabDr_HadF->push_back(tmpLabel);

    int tmpLabelDoubleHadron = -1;
    try {
      jet->getAttribute("HadronConeExclExtendedTruthLabelID", tmpLabelDoubleHadron);
    } catch(...) {};
    v_jet_DoubleHadLabel->push_back(tmpLabelDoubleHadron);


    // requested by P.Berta
    float mindRtoB = 10;
    float mindRtoC = 10;
    float mindRtoT = 10;
    for (unsigned int ip = 0; ip < m_partonB.size(); ip++) {
      float dr = deltaR(jet->eta(), (m_partonB.at(ip))->eta(), jet->phi(), (m_partonB.at(ip))->phi());
      if (dr < mindRtoB) mindRtoB = dr;
    }
    for (unsigned int ip = 0; ip < m_partonC.size(); ip++) {
      float dr = deltaR(jet->eta(), (m_partonC.at(ip))->eta(), jet->phi(), (m_partonC.at(ip))->phi());
      if (dr < mindRtoC) mindRtoC = dr;
    }
    for (unsigned int ip = 0; ip < m_partonT.size(); ip++) {
      float dr = deltaR(jet->eta(), (m_partonT.at(ip))->eta(), jet->phi(), (m_partonT.at(ip))->phi());
      if (dr < mindRtoT) mindRtoT = dr;
    }
    v_jet_dRminToB->push_back(mindRtoB);
    v_jet_dRminToC->push_back(mindRtoC);
    v_jet_dRminToT->push_back(mindRtoT);

    //Flavour labelling used for LF calib (K, Lambda, HadrInteraction, Conversions)
    if(m_saveJetLFCalibType) 
      v_jet_LFCalibType->push_back( getLFCalibType(jet, myVertex) );
    
    // continue;  // VALERIO !!!!!!!!

    /////////////////////////////////////////////////////////////////
    // Get b-tag object and fill b-tagging information
    const xAOD::BTagging *bjet = jet->btagging();

    //std::cout << " got Btagging: " << std::endl;

    // IP2D
    std::vector< ElementLink< xAOD::TrackParticleContainer > > IP2DTracks;
    IP2DTracks = bjet->auxdata<std::vector<ElementLink< xAOD::TrackParticleContainer> > >("IP2D_TrackParticleLinks");
    if (IP2DTracks.size() > 0) {
      v_jet_ip2d_pb->push_back(bjet->IP2D_pb());
      v_jet_ip2d_pc->push_back(bjet->IP2D_pc());
      v_jet_ip2d_pu->push_back(bjet->IP2D_pu());
      v_jet_ip2d_llr->push_back(bjet->IP2D_loglikelihoodratio());
    }
    else {
      v_jet_ip2d_pb->push_back(-99);
      v_jet_ip2d_pc->push_back(-99);
      v_jet_ip2d_pu->push_back(-99);
      v_jet_ip2d_llr->push_back(-99);
    }

    // IP3D
    std::vector< ElementLink< xAOD::TrackParticleContainer > > IP3DTracks;
    IP3DTracks = bjet->auxdata<std::vector<ElementLink< xAOD::TrackParticleContainer> > >("IP3D_TrackParticleLinks");
    if (IP3DTracks.size()) {
      v_jet_ip3d_pb->push_back(bjet->IP3D_pb());
      v_jet_ip3d_pc->push_back(bjet->IP3D_pc());
      v_jet_ip3d_pu->push_back(bjet->IP3D_pu());
      v_jet_ip3d_llr->push_back(bjet->IP3D_loglikelihoodratio());
    }
    else {
      v_jet_ip3d_pb->push_back(-99);
      v_jet_ip3d_pc->push_back(-99);
      v_jet_ip3d_pu->push_back(-99);
      v_jet_ip3d_llr->push_back(-99);
    }

    float tmp_trkSum_SPt = 0;
    float tmp_trkSum_VPt = 0;
    float tmp_trkSum_VEta = 0;
    float tmp_trkSum_VPhi = 0;
    int tmp_trkSum_nTrk = 0;

    /////////////////////////////////////////////////////////////////////
    // getting manual distribution on the number of tracks per jet with the selection
    TrackLinks assocTracks = bjet->auxdata<TrackLinks>(m_track_associator);

    // temporary track loop - sums up the 4vectors of all valid b-tag tracks and outputs
    TLorentzVector pseudoTrackJet(0, 0, 0, 0);
    for (unsigned int iT = 0; iT < assocTracks.size(); iT++) {
      if (!assocTracks.at(iT).isValid()) continue;
      const xAOD::TrackParticle *tmpTrk = *(assocTracks.at(iT));


      if (m_InDetTrackSelectorTool->accept(*tmpTrk, myVertex) && m_TightTrackVertexAssociationTool->isCompatible(*tmpTrk, *myVertex) ) {
	TLorentzVector tmpTrack(0, 0, 0, 0);
	tmpTrack.SetPtEtaPhiM(tmpTrk->pt(), tmpTrk->eta(), tmpTrk->phi(), 0);
	pseudoTrackJet += tmpTrack;
	tmp_trkSum_SPt += tmpTrk->pt();
	tmp_trkSum_nTrk++;
      }
      else {
	//ATH_MSG_INFO(" .... track failed");
      }
    }
    tmp_trkSum_VPt = pseudoTrackJet.Pt();
    tmp_trkSum_VEta = pseudoTrackJet.Eta();
    tmp_trkSum_VPhi = pseudoTrackJet.Phi();
    v_jet_sumtrkS_pt->push_back(tmp_trkSum_SPt);
    v_jet_sumtrkV_pt->push_back(tmp_trkSum_VPt);
    v_jet_sumtrkV_eta->push_back(tmp_trkSum_VEta);
    v_jet_sumtrkV_phi->push_back(tmp_trkSum_VPhi);
    v_jet_sumtrk_ntrk->push_back(tmp_trkSum_nTrk);

    //std::cout << "  after track loop " << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // SV0 // VD: check the existence of the vertex and only then fill the variables
    // this mimics what's done in MV2
    for (auto* svx: m_svx_branches) {
      svx->fill(*bjet);
    }

    // SV1 // VD: check the existence of the vertex and only then fill the variables
    const std::vector<ElementLink<xAOD::VertexContainer > > SV1vertices = bjet->auxdata<std::vector<ElementLink<xAOD::VertexContainer > > >("SV1_vertices");
    if (SV1vertices.size() != 0) {
      v_jet_sv1_pb->push_back(bjet->SV1_pb());
      v_jet_sv1_pc->push_back(bjet->SV1_pc());
      v_jet_sv1_pu->push_back(bjet->SV1_pu());
      v_jet_sv1_llr->push_back(bjet->SV1_loglikelihoodratio());
    }
    else {
      v_jet_sv1_pb->push_back(-99);
      v_jet_sv1_pc->push_back(-99);
      v_jet_sv1_pu->push_back(-99);
      v_jet_sv1_llr->push_back(-99);
    }

    //JetFitter
    m_jetfitter_branches.fill(*jet);


    // DL1
    try {
      v_jet_dl1_pb->push_back(bjet->auxdata<double>("DL1_pb"));
      v_jet_dl1_pc->push_back(bjet->auxdata<double>("DL1_pc"));
      v_jet_dl1_pu->push_back(bjet->auxdata<double>("DL1_pu"));
    } catch(...) {}
    try {
      v_jet_dl1mu_pb->push_back(bjet->auxdata<double>("DL1mu_pb"));
      v_jet_dl1mu_pc->push_back(bjet->auxdata<double>("DL1mu_pc"));
      v_jet_dl1mu_pu->push_back(bjet->auxdata<double>("DL1mu_pu"));
    } catch(...) {}
      try {
      v_jet_dl1rnn_pb->push_back(bjet->auxdata<double>("DL1rnn_pb"));
      v_jet_dl1rnn_pc->push_back(bjet->auxdata<double>("DL1rnn_pc"));
      v_jet_dl1rnn_pu->push_back(bjet->auxdata<double>("DL1rnn_pu"));
    } catch(...) {}

    // Other
    v_jet_sv1ip3d->push_back(bjet->SV1plusIP3D_discriminant());
    try{
      v_jet_mv1    ->push_back(bjet->MV1_discriminant());
      v_jet_mv1c   ->push_back(bjet->auxdata<double>("MV1c_discriminant"));
    } catch(...){ }


      try { v_jet_mv2c00->push_back(bjet->auxdata<double>("MV2c00_discriminant")); } catch(...) { }
      try { v_jet_mv2c10->push_back(bjet->auxdata<double>("MV2c10_discriminant")); } catch(...) { }
      try { v_jet_mv2c10mu->push_back(bjet->auxdata<double>("MV2c10mu_discriminant")); } catch(...) { }
      try { v_jet_mv2c10rnn->push_back(bjet->auxdata<double>("MV2c10rnn_discriminant")); } catch(...) { }
      try { v_jet_mv2c20->push_back(bjet->auxdata<double>("MV2c20_discriminant")); } catch(...) { }
      try { v_jet_mv2c100->push_back(bjet->auxdata<double>("MV2c100_discriminant")); } catch(...) { }
      try { v_jet_mv2cl100->push_back(bjet->auxdata<double>("MV2cl100_discriminant")); } catch(...) { }
      try { v_jet_mv2c10flip->push_back(bjet->auxdata<double>("MV2c10Flip_discriminant")); } catch(...) { }
      try { v_jet_mv2c20flip->push_back(bjet->auxdata<double>("MV2c20Flip_discriminant")); } catch(...) { }


    try {
      v_jet_mv2m_pu->push_back(bjet->auxdata<double>("MV2m_pu"));
      v_jet_mv2m_pc->push_back(bjet->auxdata<double>("MV2m_pc"));
      v_jet_mv2m_pb->push_back(bjet->auxdata<double>("MV2m_pb"));
    } catch(...) { }

    try {
      v_jet_mvb    ->push_back(bjet->auxdata<double>("MV2c10b_discriminant"));
    } catch(...) { }


    if (m_doMSV){
      // MSV
      // need initial values if no msv vertex is find, fix in MSVvariablesFactory
      double msv1;
      double msv2;
      int msv_n2t;
      float msv_ergtkjet;
      int msv_nvsec;
      float msv_nrdist;

      bjet->variable<double>("MultiSVbb1", "discriminant", msv1);
      v_jet_multisvbb1->push_back(msv1);

      bjet->variable<double>("MultiSVbb2", "discriminant", msv2);
      v_jet_multisvbb2->push_back(msv2);

      bjet->variable<int>("MSV", "N2Tpair", msv_n2t);
      v_jet_msv_N2Tpair->push_back(msv_n2t);

      bjet->variable<float>("MSV", "energyTrkInJet", msv_ergtkjet);
      v_jet_msv_energyTrkInJet->push_back(msv_ergtkjet);

      bjet->variable<int>("MSV", "nvsec", msv_nvsec);
      v_jet_msv_nvsec->push_back(msv_nvsec);

      bjet->variable<float>("MSV", "normdist", msv_nrdist);
      v_jet_msv_normdist->push_back(msv_nrdist);

      std::vector< ElementLink< xAOD::VertexContainer > > msvVertices;
      bjet->variable<std::vector<ElementLink<xAOD::VertexContainer> > >("MSV", "vertices", msvVertices);

      //tmp vectors
      std::vector<float> j_msv_cov0 = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_cov1 = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_cov2 = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_cov3 = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_cov4 = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_cov5 = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_mass = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_efrc = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_ntrk = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_pt = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_eta = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_phi = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_dls = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_xp = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_yp = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_zp = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_chi = std::vector<float>(msv_nvsec, 0);
      std::vector<float> j_msv_ndf = std::vector<float>(msv_nvsec, 0);

      // loop over vertices
      int ivtx = 0;
      if (msvVertices.size() > 0) {
        const std::vector<ElementLink<xAOD::VertexContainer> >::const_iterator verticesEnd = msvVertices.end();
        for (std::vector<ElementLink<xAOD::VertexContainer> >::const_iterator vtxIter = msvVertices.begin(); vtxIter != verticesEnd; ++vtxIter){
          if (msvVertices.size() >= 10) continue;
          float mass = xAOD::SecVtxHelper::VertexMass(**vtxIter);
          float efrc = xAOD::SecVtxHelper::EnergyFraction(**vtxIter);
          int ntrk = xAOD::SecVtxHelper::VtxNtrk(**vtxIter);
          float pt = xAOD::SecVtxHelper::Vtxpt(**vtxIter);
          float eta = xAOD::SecVtxHelper::Vtxeta(**vtxIter);
          float phi = xAOD::SecVtxHelper::Vtxphi(**vtxIter);
          float dls = xAOD::SecVtxHelper::VtxnormDist(**vtxIter);
          float xp = (**vtxIter)->x();
          float yp = (**vtxIter)->y();
          float zp = (**vtxIter)->z();
          float chi = (**vtxIter)->chiSquared();
          float ndf = (**vtxIter)->numberDoF();
          std::vector<float> covariantMatrix = (**vtxIter)->covariance();

          if (ivtx < 10) {
            j_msv_cov0[ivtx] = covariantMatrix.at(0);
            j_msv_cov1[ivtx] = covariantMatrix.at(1);
            j_msv_cov2[ivtx] = covariantMatrix.at(2);
            j_msv_cov3[ivtx] = covariantMatrix.at(3);
            j_msv_cov4[ivtx] = covariantMatrix.at(4);
            j_msv_cov5[ivtx] = covariantMatrix.at(5);
            j_msv_mass[ivtx] = mass;
            j_msv_efrc[ivtx] = efrc;
            j_msv_ntrk[ivtx] = ntrk;
            j_msv_pt[ivtx] = pt;
            j_msv_eta[ivtx] = eta;
            j_msv_phi[ivtx] = phi;
            j_msv_dls[ivtx] = dls;
            j_msv_xp[ivtx] = xp;
            j_msv_yp[ivtx] = yp;
            j_msv_zp[ivtx] = zp;
            j_msv_chi[ivtx] = chi;
            j_msv_ndf[ivtx] = ndf;
            ivtx++;
          }
        }
      }

      // fill info per vertex
      v_jet_msv_vtx_cov0->push_back(j_msv_cov0);
      v_jet_msv_vtx_cov1->push_back(j_msv_cov1);
      v_jet_msv_vtx_cov2->push_back(j_msv_cov2);
      v_jet_msv_vtx_cov3->push_back(j_msv_cov3);
      v_jet_msv_vtx_cov4->push_back(j_msv_cov4);
      v_jet_msv_vtx_cov5->push_back(j_msv_cov5);
      v_jet_msv_vtx_mass->push_back(j_msv_mass);
      v_jet_msv_vtx_efrc->push_back(j_msv_efrc);
      v_jet_msv_vtx_ntrk->push_back(j_msv_ntrk);
      v_jet_msv_vtx_pt->push_back(j_msv_pt);
      v_jet_msv_vtx_eta->push_back(j_msv_eta);
      v_jet_msv_vtx_phi->push_back(j_msv_phi);
      v_jet_msv_vtx_dls->push_back(j_msv_dls);
      v_jet_msv_vtx_x->push_back(j_msv_xp);
      v_jet_msv_vtx_y->push_back(j_msv_yp);
      v_jet_msv_vtx_z->push_back(j_msv_eta);
      v_jet_msv_vtx_chi->push_back(j_msv_chi);
      v_jet_msv_vtx_ndf->push_back(j_msv_ndf);

    } // end m_doMSV


    ///////////some jetfitter stuff needed for MVb variables /////////////
    float jfefc       = -99;
    int   jfntrkAtVx  = -1;
    int   jfnvtx1t    = -1;

    std::vector<ElementLink<xAOD::BTagVertexContainer> > jfvertices;
    try {
      jfvertices =  bjet->auxdata<std::vector<ElementLink<xAOD::BTagVertexContainer> > >("JetFitter_JFvertices");
    } catch (...) {};

    int tmpNvtx = 0;
    int tmpNvtx1t = 0;

    bjet->taggerInfo(tmpNvtx, xAOD::JetFitter_nVTX);
    bjet->taggerInfo(tmpNvtx1t, xAOD::JetFitter_nSingleTracks); // 1 track vertices

    if (tmpNvtx > 0 || tmpNvtx1t > 0) {

      bjet->taggerInfo(jfefc, xAOD::JetFitter_energyFraction);
      bjet->taggerInfo(jfntrkAtVx, xAOD::JetFitter_nTracksAtVtx);
      jfnvtx1t = tmpNvtx1t;
    }
    /////////////////////////////////////////////////////////////////
    // Generating MVb variables (as in MV2Tag.cxx)
    bool trksOK=IP3DTracks.size();

    std::vector<float> vectD0, vectD0Signi, vectZ0, vectZ0Signi;    vectD0.clear(), vectD0Signi.clear(), vectZ0.clear(), vectZ0Signi.clear();
    trksOK &= bjet->variable< std::vector<float> > ("IP3D", "valD0wrtPVofTracks", vectD0     );
    trksOK &= bjet->variable< std::vector<float> > ("IP3D", "sigD0wrtPVofTracks", vectD0Signi);
    trksOK &= bjet->variable< std::vector<float> > ("IP3D", "valZ0wrtPVofTracks", vectZ0     );
    trksOK &= bjet->variable< std::vector<float> > ("IP3D", "sigZ0wrtPVofTracks", vectZ0Signi);
    if (vectD0.size() and vectD0Signi.size() and vectZ0.size() and vectZ0Signi.size()) {
      trksOK &= IP3DTracks.size() == vectD0.size();
      trksOK &= IP3DTracks.size() == vectZ0.size();
      trksOK &= IP3DTracks.size() == vectD0Signi.size();
      trksOK &= IP3DTracks.size() == vectZ0Signi.size();
    }
    //std::cout<<"debug: "<<IP3DTracks.size()<<" "<<vectD0.size()<<" "<<vectZ0.size()<<" "<<vectD0Signi.size()<<" "<<vectZ0Signi.size()<<std::endl;

    int ntrks = IP3DTracks.size();
    float width = 0;
    int   n_trk_d0cut = 0;
    float trk3_d0sig = -100;
    float trk3_z0sig = -100;
    float sv_scaled_efc = -1;
    float jf_scaled_efc = -1;
    if (trksOK) {
      ATH_MSG_VERBOSE("#BTAG# MV2: calculating MVb inputs.");

      float sum_pt = 0., sum_pt_dr = 0.;

      std::vector<std::pair<float, float> > trk_d0_z0;
      trk_d0_z0.reserve(IP3DTracks.size());

      unsigned trkIndex=0;
      for(auto trkIter = IP3DTracks.begin(); trkIter != IP3DTracks.end(); ++trkIter) {
        const xAOD::TrackParticle* aTemp = **trkIter;
        TLorentzVector trk;
        trk.SetPtEtaPhiM(aTemp->pt(), aTemp->eta(), aTemp->phi(), 0.);

        // no need for a dedicated selection here, the tracks are already
        // selected by the IP3D algorithm
        const float d0sig = vectD0Signi.at(trkIndex);
        const float z0sig = vectZ0Signi.at(trkIndex);
        trkIndex++;

        if (std::fabs(d0sig) > 1.8)
          n_trk_d0cut++;

        // track width components
        sum_pt += trk.Pt();
        const float dRtoJet = trk.DeltaR(jet->p4());
        sum_pt_dr += dRtoJet * trk.Pt();

        // for 3rd higest d0/z0 significance
        trk_d0_z0.push_back(std::make_pair(d0sig, z0sig));
      } //end of trk loop

      // sort by highest signed d0 sig
      std::sort(trk_d0_z0.begin(), trk_d0_z0.end(), [](const std::pair<float, float>& a, const std::pair<float, float>& b) {
        return a.first > b.first;
      } );

      //Assign MVb variables
      if (sum_pt > 0) width = sum_pt_dr / sum_pt;
      if (trk_d0_z0.size() > 2) trk3_d0sig = trk_d0_z0[2].first;
      if (trk_d0_z0.size() > 2) trk3_z0sig = trk_d0_z0[2].second;
      int sv1ntrkv = bjet->auxdata<int>("SV1_NGTinSvx");
      float sv1efc = bjet->auxdata<float>("SV1_efracsvx");
      if (sv1ntrkv>0) sv_scaled_efc  =  sv1efc * (static_cast<float>(ntrks) / sv1ntrkv);
      if (jfntrkAtVx + jfnvtx1t>0) jf_scaled_efc  =  jfefc * (static_cast<float>(ntrks) / (jfntrkAtVx + jfnvtx1t));
    }

    v_jet_width->push_back(width);
    v_jet_n_trk_sigd0cut->push_back(n_trk_d0cut);
    v_jet_trk3_d0sig->push_back(trk3_d0sig);
    v_jet_trk3_z0sig->push_back(trk3_z0sig);
    v_jet_sv_scaled_efc->push_back(sv_scaled_efc);
    v_jet_jf_scaled_efc->push_back(jf_scaled_efc);

    // ExKtbbTag
    if (bjet->isAvailable<double>("ExKtbb_Hbb_DoubleMV2c20")) {

      v_jet_ExKtbb_Hbb_DoubleMV2c20->push_back(bjet->auxdata<double>("ExKtbb_Hbb_DoubleMV2c20"));

      if (bjet->isAvailable<double>("ExKtbb_Hbb_SingleMV2c20")) {
        v_jet_ExKtbb_Hbb_SingleMV2c20->push_back(bjet->auxdata<double>("ExKtbb_Hbb_SingleMV2c20"));
      }
      if (bjet->isAvailable<double>("ExKtbb_Hbb_MV2Only")) {
        v_jet_ExKtbb_Hbb_MV2Only->push_back(bjet->auxdata<double>("ExKtbb_Hbb_MV2Only"));
      }
      if (bjet->isAvailable<double>("ExKtbb_Hbb_MV2andJFDRSig")) {
        v_jet_ExKtbb_Hbb_MV2andJFDRSig->push_back(bjet->auxdata<double>("ExKtbb_Hbb_MV2andJFDRSig"));
      }
      if (bjet->isAvailable<double>("ExKtbb_Hbb_MV2andTopos")) {
        v_jet_ExKtbb_Hbb_MV2andTopos->push_back(bjet->auxdata<double>("ExKtbb_Hbb_MV2andTopos"));
      }
    }
    else {
      ATH_MSG_DEBUG("WARNING! No ExKtbbTag run on " << m_jetCollectionName.c_str());
    }

    for (const auto& collection_branches: m_subjet_branches) {
      const auto& name = collection_branches.first;
      // std::vector<const xAOD::Jet*> subjet;
      auto subjet = jet->getAssociatedObjects<xAOD::Jet>(name);
      collection_branches.second->fill(subjet);
    }

    // now the tracking part: prepare all the tmpVectors
    int j_btag_ntrk = 0;
    int j_sv1_ntrk = 0;
    int j_ip3d_ntrk = 0;
    int j_jf_ntrk = 0;

    std::vector<float> j_trk_dr; // mod nikola
    std::vector<int> j_trk_assoc_msv; // mod nikola
    std::vector<int> j_trk_algo;
    //std::vector<int> j_trk_orig; //moved to BHadronBranches
    std::vector<int> j_trk_pdgid;
    std::vector<int> j_trk_barcode;
    std::vector<int> j_trk_parent_pdgid;
    std::vector<int> j_trk_cploose;
    std::vector<float> j_trk_vtx_X;
    std::vector<float> j_trk_vtx_Y;
    std::vector<float> j_trk_vtx_Z;
    std::vector<float> j_trk_vtx_dx;
    std::vector<float>  j_trk_vtx_dy;
    std::vector<float> j_trk_d0_truth;
    std::vector<float> j_trk_z0_truth;
    std::vector<int> j_trk_ip3d_grade;
    std::vector<float> j_trk_ip3d_d0;
    std::vector<float> j_trk_ip3d_d02D;
    std::vector<float> j_trk_ip3d_z0;
    std::vector<float> j_trk_ip3d_d0sig;
    std::vector<float> j_trk_ip3d_z0sig;
    std::vector<float> j_trk_ip3d_llr;

    // if (m_reduceInfo) continue;

    bool is8TeV = true;
    if (bjet->isAvailable<std::vector<ElementLink<xAOD::BTagVertexContainer> > >("JetFitter_JFvertices")) is8TeV = false;

    TrackLinks SV0Tracks ;
    TrackLinks SV1Tracks ;
    TrackLinks JFTracks;

    if (!is8TeV) {
      if (bjet->isAvailable<TrackLinks>("SV0_TrackParticleLinks")) {
        SV0Tracks = bjet->SV0_TrackParticleLinks();
      }
      SV1Tracks = bjet->SV1_TrackParticleLinks();
    }

    double jet_mu_dRmin_smt=999;
    float jet_mu_dRmin_pt=999,jet_mu_dRmin_dR=999,jet_mu_dRmin_truthflav=999,jet_mu_dRmin_eta=999,jet_mu_dRmin_phi=999,jet_mu_dRmin_assJet_pt=999,jet_mu_dRmin_qOverPratio=999,jet_mu_dRmin_mombalsignif=999,jet_mu_dRmin_scatneighsignif=999,jet_mu_dRmin_pTrel=999,jet_mu_dRmin_VtxTyp=999,jet_mu_dRmin_d0=999,jet_mu_dRmin_z0=999,jet_mu_dRmin_parent_pdgid=999,jet_mu_dRmin_ID_qOverP_var=999,jet_mu_dRmin_muonType=999;
    float jet_mu_fatjet_nMu = 0, jet_mu_fatjet_pTmax_pT = 999, jet_mu_fatjet_pTmax_pTrel = 999, jet_mu_fatjet_pTmax_pTrelFrac = 999;
    if (m_SMT) {
      // additions by nikola
      try {
	std::vector<ElementLink<xAOD::MuonContainer> > assocMuons;
        assocMuons = bjet->auxdata<std::vector<ElementLink<xAOD::MuonContainer> > >("Muons");
        if (assocMuons.size() != 0) {
          for (unsigned int iT = 0; iT < assocMuons.size(); iT++) {
            if (!assocMuons.at(iT).isValid()) continue;
            const xAOD::Muon *tmpMuon = *(assocMuons.at(iT));
            float dr = deltaR(tmpMuon->eta(), jet->eta(), tmpMuon->phi(), jet->phi());
            // const ElementLink< xAOD::TrackParticleContainer >& pMuIDTrack = tmpMuon->inDetTrackParticleLink();
            const ElementLink< xAOD::TrackParticleContainer >& pMuMSTrack = tmpMuon->muonSpectrometerTrackParticleLink();
            // const xAOD::Vertex *pVtx = (*pMuIDTrack)->vertex();
            // const std::vector<float>&cov = (*pMuIDTrack)->definingParametersCovMatrixVec();
            float momBalSignif0 = 999.;
            tmpMuon->parameter(momBalSignif0, xAOD::Muon::momentumBalanceSignificance);
            if (momBalSignif0 == 0) continue;
            if ((*pMuMSTrack)->qOverP() == 0) continue;
            if (dr >= 1.0) continue;
            jet_mu_fatjet_nMu += 1;
            if (tmpMuon->pt() > jet_mu_fatjet_pTmax_pT) {
              jet_mu_fatjet_pTmax_pT = tmpMuon->pt();
              TLorentzVector myjet, mymu;
              myjet.SetPtEtaPhiM(jet->pt(), jet->eta(), jet->phi(), 0);
              mymu.SetPtEtaPhiM(tmpMuon->pt(), tmpMuon->eta(), tmpMuon->phi(), 0);
              jet_mu_fatjet_pTmax_pTrel = myjet.Vect().Perp(mymu.Vect()) / 1000;
              jet_mu_fatjet_pTmax_pTrelFrac = jet_mu_fatjet_pTmax_pTrel / jet->pt();
            }
          }
        }
      } catch(...) {
        ATH_MSG_INFO("NO Muons found!");
      }

      // new from Valerio: if the variables are already available, do not calculate them
      if ( bjet->isAvailable<float>("SMT_mu_pt") ) {
	//std::cout << "SMT info already available, will get them from there ... " << std::endl;
	jet_mu_dRmin_smt            = bjet->auxdata<double>("SMT_discriminant");
	jet_mu_dRmin_dR             = bjet->auxdata<float>("SMT_dR");
	jet_mu_dRmin_pt             = bjet->auxdata<float>("SMT_mu_pt");
	jet_mu_dRmin_qOverPratio    = bjet->auxdata<float>("SMT_qOverPratio");
	jet_mu_dRmin_mombalsignif   = bjet->auxdata<float>("SMT_mombalsignif");
	jet_mu_dRmin_scatneighsignif= bjet->auxdata<float>("SMT_scatneighsignif");
	jet_mu_dRmin_pTrel          = bjet->auxdata<float>("SMT_pTrel");
	jet_mu_dRmin_d0             = bjet->auxdata<float>("SMT_mu_d0");
	jet_mu_dRmin_z0             = bjet->auxdata<float>("SMT_mu_z0");
	jet_mu_dRmin_ID_qOverP_var  = bjet->auxdata<float>("SMT_ID_qOverP");
	jet_mu_dRmin_assJet_pt=jet->pt()/1000; // ?? why is this variable needed?
	jet_mu_dRmin_truthflav=thisJetTruthLabel;

	ElementLink<xAOD::MuonContainer> tmpMuonLink= bjet->auxdata<ElementLink<xAOD::MuonContainer> >("SMT_mu_link");
	if ( tmpMuonLink.isValid() ) {
	  const xAOD::Muon* tmpMuon=(*tmpMuonLink);
	  //std::cout << " link is: " << tmpMuon << std::endl;
	  if ( tmpMuon!=0 ) {
	    jet_mu_dRmin_eta      =tmpMuon->eta();
	    jet_mu_dRmin_phi      =tmpMuon->phi();
	    jet_mu_dRmin_muonType =tmpMuon->muonType();
	    //std::cout << " .... after eta and friends" << std::endl;
	    const ElementLink< xAOD::TrackParticleContainer >& pMuIDTrack=tmpMuon->inDetTrackParticleLink();
	    //std::cout << "   the link is: " << pMuIDTrack << std::endl;
	    const xAOD::Vertex * pVtx=(*pMuIDTrack)->vertex();
	    if(pVtx!=NULL) {
	      jet_mu_dRmin_VtxTyp=pVtx->vertexType();
	    } else {jet_mu_dRmin_VtxTyp=999.;}
	    const xAOD::TruthParticle* matched_truth_muon=0;
	    if(tmpMuon->isAvailable<ElementLink<xAOD::TruthParticleContainer> >("truthParticleLink")) {
	      ElementLink<xAOD::TruthParticleContainer> link = tmpMuon->auxdata<ElementLink<xAOD::TruthParticleContainer> >("truthParticleLink");
	      if(link.isValid()) {
		matched_truth_muon = *link;
		int pdgid = parent_classify(matched_truth_muon);
		jet_mu_dRmin_parent_pdgid=pdgid;
	      } else {jet_mu_dRmin_parent_pdgid=999.;}
	    }
	  }
	}

      } else {
	try{
	  std::vector<ElementLink<xAOD::MuonContainer> > assocMuons;
	  assocMuons= bjet->auxdata<std::vector<ElementLink<xAOD::MuonContainer> > >("Muons");
	  if(assocMuons.size()!=0){
	    for (unsigned int iT=0; iT<assocMuons.size(); iT++) {
	      if (!assocMuons.at(iT).isValid()) continue;
	      const xAOD::Muon* tmpMuon= *(assocMuons.at(iT));
	      float dr = deltaR(tmpMuon->eta(),jet->eta(),tmpMuon->phi(),jet->phi());
	      if(dr>=0.4) continue;
	      const ElementLink< xAOD::TrackParticleContainer >& pMuIDTrack=tmpMuon->inDetTrackParticleLink();
	      const ElementLink< xAOD::TrackParticleContainer >& pMuMSTrack=tmpMuon->muonSpectrometerTrackParticleLink();
	      const xAOD::Vertex * pVtx=(*pMuIDTrack)->vertex();
	      const std::vector<float>&cov= (*pMuIDTrack)->definingParametersCovMatrixVec();
	      float momBalSignif0=999.;
	      tmpMuon->parameter(momBalSignif0, xAOD::Muon::momentumBalanceSignificance);
	      if(momBalSignif0==0) continue;
	      if((*pMuMSTrack)->qOverP()==0) continue;
	      if(dr<jet_mu_dRmin_dR){
		jet_mu_dRmin_dR=dr;
		jet_mu_dRmin_pt=tmpMuon->pt()/1000;
		jet_mu_dRmin_truthflav=thisJetTruthLabel;
		jet_mu_dRmin_eta=tmpMuon->eta();
		jet_mu_dRmin_phi=tmpMuon->phi();
		jet_mu_dRmin_assJet_pt=jet->pt()/1000;
		jet_mu_dRmin_qOverPratio=(*pMuIDTrack)->qOverP()/(*pMuMSTrack)->qOverP();
		float momBalSignif=999.;
		if(tmpMuon->parameter(momBalSignif, xAOD::Muon::momentumBalanceSignificance)) {
		  jet_mu_dRmin_mombalsignif=momBalSignif;
		} else {jet_mu_dRmin_mombalsignif=momBalSignif;}
		float scatNeighSignif=999.;
		if(tmpMuon->parameter(scatNeighSignif, xAOD::Muon::scatteringNeighbourSignificance)) {
		  jet_mu_dRmin_scatneighsignif=scatNeighSignif;
		} else {jet_mu_dRmin_scatneighsignif=scatNeighSignif;}
		TLorentzVector myjet, mymu;
		myjet.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),0);
		mymu.SetPtEtaPhiM(tmpMuon->pt(),tmpMuon->eta(),tmpMuon->phi(),0);
		jet_mu_dRmin_pTrel=myjet.Vect().Perp(mymu.Vect())/1000;
		if(pVtx!=NULL) {
		  jet_mu_dRmin_VtxTyp=pVtx->vertexType();
		} else {jet_mu_dRmin_VtxTyp=999.;}
		jet_mu_dRmin_d0=tmpMuon->primaryTrackParticle()->d0();
		jet_mu_dRmin_z0=tmpMuon->primaryTrackParticle()->z0();

		const xAOD::TruthParticle* matched_truth_muon=0;
		if(tmpMuon->isAvailable<ElementLink<xAOD::TruthParticleContainer> >("truthParticleLink")) {
		  ElementLink<xAOD::TruthParticleContainer> link = tmpMuon->auxdata<ElementLink<xAOD::TruthParticleContainer> >("truthParticleLink");
		  if(link.isValid()) {
		    matched_truth_muon = *link;
		    int pdgid = parent_classify(matched_truth_muon);
		    jet_mu_dRmin_parent_pdgid=pdgid;
		  } else {jet_mu_dRmin_parent_pdgid=999.;}
		}
		jet_mu_dRmin_ID_qOverP_var=cov[14];
		jet_mu_dRmin_muonType=tmpMuon->muonType();
	      }
	    }
	  }
	} catch(...) {
	  //std::cout << "NO Muons found!"<<std::endl;
	  //todo: write out some warning here but don't want to clog logfiles for now
	}
      }
     }

    if (m_SMT) {
      v_jet_mu_smt->push_back(jet_mu_dRmin_smt);
      v_jet_mu_assJet_pt->push_back(jet_mu_dRmin_assJet_pt);
      v_jet_mu_truthflav->push_back(jet_mu_dRmin_truthflav);
      v_jet_mu_pt->push_back(jet_mu_dRmin_pt);
      v_jet_mu_eta->push_back(jet_mu_dRmin_eta);
      v_jet_mu_phi->push_back(jet_mu_dRmin_phi);
      v_jet_mu_qOverPratio->push_back(jet_mu_dRmin_qOverPratio);
      v_jet_mu_mombalsignif->push_back(jet_mu_dRmin_mombalsignif);
      v_jet_mu_scatneighsignif->push_back(jet_mu_dRmin_scatneighsignif);
      v_jet_mu_dR->push_back(jet_mu_dRmin_dR);
      v_jet_mu_pTrel->push_back(jet_mu_dRmin_pTrel);
      v_jet_mu_VtxTyp->push_back(jet_mu_dRmin_VtxTyp);
      v_jet_mu_d0->push_back(jet_mu_dRmin_d0);
      v_jet_mu_z0->push_back(jet_mu_dRmin_z0);
      v_jet_mu_parent_pdgid->push_back(jet_mu_dRmin_parent_pdgid);
      v_jet_mu_ID_qOverP_var->push_back(jet_mu_dRmin_ID_qOverP_var);
      v_jet_mu_muonType->push_back(jet_mu_dRmin_muonType);
      v_jet_mu_fatjet_nMu->push_back(jet_mu_fatjet_nMu);
      v_jet_mu_fatjet_pTmax_pT->push_back(jet_mu_fatjet_pTmax_pT);
      v_jet_mu_fatjet_pTmax_pTrel->push_back(jet_mu_fatjet_pTmax_pTrel);
      v_jet_mu_fatjet_pTmax_pTrelFrac->push_back(jet_mu_fatjet_pTmax_pTrelFrac);
    }

    // if (m_reduceInfo) continue;


    for (unsigned int jfv = 0; jfv < jfvertices.size(); jfv++) {
      if (!jfvertices.at(jfv).isValid()) continue;
      const xAOD::BTagVertex *tmpVertex = *(jfvertices.at(jfv));
      const std::vector< ElementLink<xAOD::TrackParticleContainer> > tmpVect = tmpVertex->track_links(); // mod Remco
      JFTracks.insert(JFTracks.end(), tmpVect.begin(), tmpVect.end()); // mod Remco

    }

    j_btag_ntrk = 0; // assocTracks.size();
    j_sv1_ntrk = SV1Tracks.size();
    j_ip3d_ntrk = IP3DTracks.size();
    j_jf_ntrk = JFTracks.size();

    std::vector<int> tmpGrading = bjet->auxdata<std::vector<int> >("IP3D_gradeOfTracks");
    std::vector<float> tmpD0 = bjet->auxdata<std::vector<float> >("IP3D_valD0wrtPVofTracks");
    std::vector<float> tmpZ0 = bjet->auxdata<std::vector<float> >("IP3D_valZ0wrtPVofTracks");
    std::vector<float> tmpD0sig= bjet->auxdata<std::vector<float> >("IP3D_sigD0wrtPVofTracks");
    std::vector<float> tmpZ0sig= bjet->auxdata<std::vector<float> >("IP3D_sigZ0wrtPVofTracks");

    std::vector<float> tmpIP3DBwgt= bjet->auxdata<std::vector<float> >("IP3D_weightBofTracks");
    std::vector<float> tmpIP3DUwgt= bjet->auxdata<std::vector<float> >("IP3D_weightUofTracks");

    float ip3d_llr = MISSING_VALUE;

    j_ip3d_ntrk = tmpGrading.size();

    //////////////////////////////////////////////////////////////////////

    if (m_dumpGATracks) {
      fillGhostTracks(*jet, *myVertex);
    }

    if (m_reduceInfo) continue; // if using reduceInfo, don't output track information
    // jet direction:
    TLorentzVector jetDir;
    jetDir.SetPtEtaPhiE(v_jet_pt_orig->at(j), v_jet_eta_orig->at(j), v_jet_phi_orig->at(j), v_jet_E_orig->at(j));
    // Amg::Vector3D jetDirection(jet->px(),jet->py(),jet->pz());
    Amg::Vector3D jetDirection(jetDir.Px(),jetDir.Py(),jetDir.Pz());
    Amg::Vector3D unit = jetDirection.unit();


    // MAIN TRACK LOOP
    std::vector<const xAOD::TrackParticle*> selectedTracks;
    for (unsigned int iT = 0; iT < assocTracks.size(); iT++) {
      // std::cout << " .... trk link: " << iT << std::endl;
      if (!assocTracks.at(iT).isValid()) continue;
      const xAOD::TrackParticle *tmpTrk = *(assocTracks.at(iT));
      tmpTrk->summaryValue(getInt, xAOD::numberOfPixelHits);
      int nSi = getInt;
      tmpTrk->summaryValue(getInt, xAOD::numberOfSCTHits);
      nSi += getInt;
      if (nSi < 2) continue;
      selectedTracks.push_back(tmpTrk);
    }

    // fill cov branches
    Particles track_particles(selectedTracks.begin(), selectedTracks.end());
    m_track_cov_branches.fill(track_particles);
    m_track_branches.fill(track_particles);
    for (const auto* tmpTrk: selectedTracks) {
      j_btag_ntrk++;

      // addition by nikola
      TLorentzVector jet4vector, trk4vector;
      jet4vector.SetPtEtaPhiM(jet->pt(), jet->eta(), jet->phi(), jet->m());
      trk4vector.SetPtEtaPhiM(tmpTrk->pt(), tmpTrk->eta(), tmpTrk->phi(), 0);
      j_trk_dr.push_back(jet4vector.DeltaR(trk4vector));

      // if doing MSV, keep track of which MSV this track belongs to
      if (m_doMSV) {

           // build track to vertex maps (from track_to_vertex_associators)
          try{
              auto msv_vtx_map = trkvx::get_msv_map(*bjet);

              int trk_assoc_msv = msv_vtx_map.count(tmpTrk) ?
              msv_vtx_map.at(tmpTrk) : -1;
              j_trk_assoc_msv.push_back(trk_assoc_msv);
          } catch(...){ }
      }

      // algo
      unsigned int trackAlgo = 0;
      int index = -1;

      for (unsigned int iT = 0; iT < IP3DTracks.size(); iT++) {
        if (tmpTrk == *(IP3DTracks.at(iT))) {
          trackAlgo += 1 << IP3D;
          index = iT;
          break;
	}
      }

      if (index!=-1) {
	j_trk_ip3d_grade.push_back(tmpGrading.at(index));
	ip3d_llr=MISSING_VALUE;
	if (tmpIP3DUwgt.at(index)!=0) ip3d_llr = log(tmpIP3DBwgt.at(index)/tmpIP3DUwgt.at(index));
	j_trk_ip3d_llr.push_back(ip3d_llr);
      } else {
  j_trk_ip3d_grade.push_back(-10);
	j_trk_ip3d_llr.push_back(MISSING_VALUE);
      }
      if (particleInCollection(tmpTrk, IP2DTracks)) trackAlgo += 1 << IP2D;
      if (particleInCollection(tmpTrk, SV0Tracks)) trackAlgo +=1 << SV0;
      if (particleInCollection(tmpTrk, SV1Tracks)) trackAlgo +=1 << SV1;
      if (particleInCollection(tmpTrk, JFTracks)) trackAlgo +=1 << JF; // mod Remco
      j_trk_algo.push_back(trackAlgo);

      // origin
      // moved to BHadronBranches

      const xAOD::TruthParticle *truth = truthParticle(tmpTrk);

      if (truth) {

        j_trk_pdgid.push_back(truth->pdgId());
        j_trk_barcode.push_back(truth->barcode());
        j_trk_parent_pdgid.push_back( parent_classify(truth) );

        if (truth->prodVtx()) {
          j_trk_vtx_X.push_back(truth->prodVtx()->x());
          j_trk_vtx_Y.push_back(truth->prodVtx()->y());
          j_trk_vtx_Z.push_back(truth->prodVtx()->z());
	      }
        else {
          j_trk_vtx_X.push_back(-666);
          j_trk_vtx_Y.push_back(-666);
          j_trk_vtx_Z.push_back(-666);
	      }
      }
      else{

        j_trk_pdgid.push_back(MISSING_VALUE);
        j_trk_barcode.push_back(MISSING_VALUE);
        j_trk_parent_pdgid.push_back( MISSING_VALUE );
        j_trk_vtx_X.push_back(MISSING_VALUE);
        j_trk_vtx_Y.push_back(MISSING_VALUE);
        j_trk_vtx_Z.push_back(MISSING_VALUE);
      }


      if (!m_CPTrackingLooseLabel.empty()) {
        bool is_cp_loose = m_CPTrackingLooseLabel->accept(*tmpTrk, myVertex);
        j_trk_cploose.push_back(is_cp_loose ? 1 : 0);
      }
      // std::cout << " ..... after origin" << std::endl;
      try {
        // spatial coordinates: now from the tool:
        float d0 = m_track_accessors.d0(*tmpTrk);
        float d0Err = m_track_accessors.d0_sigma(*tmpTrk);
        float z0 = m_track_accessors.z0(*tmpTrk);
        float z0Err = m_track_accessors.z0_sigma(*tmpTrk);

        // sign of the impact parameter
        double signOfIP = get_3d_lifetime_sign(*tmpTrk, unit);
        double signOfIP2D = get_2d_lifetime_sign(*tmpTrk, unit);
        double signOfZIP = get_z_lifetime_sign(*tmpTrk, unit);

        // significances
        double sIP = signOfIP * fabs(d0);
        double significance = signOfIP * fabs(d0 / d0Err);
        double szIP = signOfZIP * fabs(z0);
        double z0Sig = signOfZIP * fabs(z0 / z0Err);

        j_trk_ip3d_d0.push_back(sIP);
        j_trk_ip3d_d02D.push_back(signOfIP2D*fabs(d0));
        j_trk_ip3d_z0.push_back(szIP);
        j_trk_ip3d_d0sig.push_back(significance);
        j_trk_ip3d_z0sig.push_back(z0Sig);
      } catch (SG::ExcBadAuxVar& exc) {
        std::string explanation(
          "missing information in a track: ");
        explanation.append(exc.what());
        explanation.append(
          ". You're probably not running the BTagTrackAugmenter."
          " Will fill with default values");
        ATH_MSG_WARNING(explanation);

        // fill with default values
        j_trk_ip3d_d0.push_back(MISSING_VALUE);
        j_trk_ip3d_d02D.push_back(MISSING_VALUE);
        j_trk_ip3d_z0.push_back(MISSING_VALUE);
        j_trk_ip3d_d0sig.push_back(MISSING_VALUE);
        j_trk_ip3d_z0sig.push_back(MISSING_VALUE);
      }

      // TRUTH track info ......
      if (!truth) {
        j_trk_d0_truth.push_back(MISSING_VALUE);
        j_trk_z0_truth.push_back(MISSING_VALUE);
      }
      else {
        float tmpd0T = MISSING_VALUE;
        float tmpz0T = MISSING_VALUE;
        try {
          tmpd0T = truth->auxdata< float >("d0");
          tmpz0T = truth->auxdata< float >("z0");
        } catch(...) {
          // todo: write out some warning here but don't want to clog logfiles for now
        }
        j_trk_d0_truth.push_back(tmpd0T);
        j_trk_z0_truth.push_back(tmpz0T);
      }
    } // end track loop

    v_jet_btag_ntrk->push_back(j_btag_ntrk);
    v_jet_sv1_ntrk->push_back(j_sv1_ntrk);
    v_jet_ip3d_ntrk->push_back(j_ip3d_ntrk);
    v_jet_jf_ntrk->push_back(j_jf_ntrk);
    v_jet_trk_dr->push_back(j_trk_dr);
    v_jet_trk_assoc_msv->push_back(j_trk_assoc_msv);
    v_jet_trk_algo->push_back(j_trk_algo);
    //v_jet_trk_orig->push_back(j_trk_orig); moved to BHadronBranches
    v_jet_trk_pdg_id->push_back(j_trk_pdgid);
    v_jet_trk_parent_pdgid->push_back(j_trk_parent_pdgid);
    v_jet_trk_barcode->push_back(j_trk_barcode);
    v_jet_trk_is_tracking_cp_loose->push_back(j_trk_cploose);
    v_jet_trk_vtx_X->push_back(j_trk_vtx_X);
    v_jet_trk_vtx_Y->push_back(j_trk_vtx_Y);
    v_jet_trk_vtx_Z->push_back(j_trk_vtx_Z);
    v_jet_trk_d0_truth->push_back(j_trk_d0_truth);
    v_jet_trk_z0_truth->push_back(j_trk_z0_truth);
    v_jet_trk_IP3D_grade->push_back(j_trk_ip3d_grade);
    v_jet_trk_IP3D_d0->push_back(j_trk_ip3d_d0);
    v_jet_trk_IP3D_d02D->push_back(j_trk_ip3d_d02D);
    v_jet_trk_IP3D_z0->push_back(j_trk_ip3d_z0);
    v_jet_trk_IP3D_d0sig->push_back(j_trk_ip3d_d0sig);
    v_jet_trk_IP3D_z0sig->push_back(j_trk_ip3d_z0sig);

    v_jet_trk_IP3D_llr->push_back(j_trk_ip3d_llr);

  } // end jet loop

  for (unsigned int j = 0; j < selJets.size(); j++) {
    delete selJets.at(j);
  }

  if (!isData) tree->Fill();
  else {
    if (njets == 0) return StatusCode::SUCCESS;
    if (v_jet_pt->at(0) < 30e3) return StatusCode::SUCCESS;
    tree->Fill();
  }

  // clear all the things that need clearing
  truth_electrons.clear();
  truth_muons.clear();
  selJets.clear();

  // addition from Dan: clear branch collections
  m_bhadron_branches.clear();
  m_kshort_branches.clear();
  m_kshortreco_branches.clear();
  m_jetfitter_branches.clear();
  m_cluster_branches.clear();
  m_substructure_moment_branches.clear();
  for (auto& coll_br: m_subjet_branches) {
    coll_br.second->clear();
  }
  for (auto& coll: m_svx_branches) {
    coll->clear();
  }
  m_track_branches.clear();
  m_track_cov_branches.clear();
  m_ga_track_branches.clear();
  m_ga_track_cov_branches.clear();

  if (m_arb_branches) m_arb_branches->clear();

  return StatusCode::SUCCESS;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// User defined functions
/////////////////////////
float btagIBLAnalysisAlg :: deltaR(float eta1, float eta2, float phi1, float phi2) {
  float DEta = fabs(eta1 - eta2);
  float DPhi = acos(cos(fabs(phi1 - phi2)));
  return sqrt(pow(DEta, 2) + pow(DPhi, 2));
}


const xAOD::Jet *btagIBLAnalysisAlg :: GetParentJet(const xAOD::Jet *Jet, std::string Keyname) {
  ElementLink<xAOD::JetContainer> el = Jet->auxdata<ElementLink<xAOD::JetContainer> >(Keyname);

  if(el.isValid()) {
    return *el;
  }
  else {
    ATH_MSG_WARNING("GetParentJet(): Unable to get parent link %s ! Null ptr is returned.");
    return 0;
  }
}


const xAOD::TruthParticle *btagIBLAnalysisAlg :: truthParticle(const xAOD::TrackParticle *trkPart) const {
  typedef ElementLink< xAOD::TruthParticleContainer > Link_t;
  static const char *NAME = "truthParticleLink";
  if( ! trkPart->isAvailable< Link_t >( NAME ) ) {
    return 0;
  }
  const Link_t& link = trkPart->auxdata< Link_t >( NAME );
  if (!link.isValid()) {
    return 0;
  }
  return *link;
}

bool GoesIntoC(const xAOD::TruthParticle *part) {
  if (!part) return false;
  if (!part->hasDecayVtx()) return false;
  const xAOD::TruthVertex *decayVtx = part->decayVtx();
  for (unsigned int ch = 0; ch < decayVtx->nOutgoingParticles(); ch++) {
    const xAOD::TruthParticle *tmpPart = decayVtx->outgoingParticle(ch);
    if (tmpPart->isCharmHadron()) return true;
  }
  return false;
}

// saving recursively only charged particle
void btagIBLAnalysisAlg :: GetParentTracks(const xAOD::TruthParticle *particle, std::vector<const xAOD::TruthParticle*> &tracksFromB, std::vector<const xAOD::TruthParticle*> &tracksFromC, bool isFromC, std::string indent) {
  if (!particle->hasDecayVtx()) return;
  const xAOD::TruthVertex *decayVtx = particle->decayVtx();
  indent += "  ";
  for (unsigned int ch = 0; ch < decayVtx->nOutgoingParticles(); ch++) {
    const xAOD::TruthParticle *tmpParticle = decayVtx->outgoingParticle(ch);
    if (tmpParticle->barcode() > 200e3) continue;
    if (!tmpParticle->isCharmHadron() && !(const_cast< xAOD::TruthParticle* >(tmpParticle))->isBottomHadron()) {
      if (tmpParticle->isCharged()) tracksFromB.push_back(tmpParticle);
      if (isFromC && tmpParticle->isCharged()) tracksFromC.push_back(tmpParticle);
    }

    if (isFromC) GetParentTracks(tmpParticle, tracksFromB, tracksFromC, true, indent);
    else GetParentTracks(tmpParticle, tracksFromB, tracksFromC, tmpParticle->isCharmHadron() && !GoesIntoC(tmpParticle), indent);
  }
}

void btagIBLAnalysisAlg :: clearvectors() {
  PV_x = MISSING_VALUE;
  PV_y = MISSING_VALUE;
  PV_z = MISSING_VALUE;
  truth_PV_x = MISSING_VALUE;
  truth_PV_y = MISSING_VALUE;
  truth_PV_z = MISSING_VALUE;
  truth_LeadJet_pt = 0;

  v_jet_pt->clear();
  v_jet_eta->clear();
  v_jet_phi->clear();
  v_jet_pt_orig->clear();
  v_jet_eta_orig->clear();
  v_jet_phi_orig->clear();
  v_jet_E_orig  ->clear();
  v_jet_sumtrkS_pt->clear();
  v_jet_sumtrkV_pt->clear();
  v_jet_sumtrkV_eta->clear();
  v_jet_sumtrkV_phi->clear();
  v_jet_sumtrk_ntrk->clear();
  v_jet_E->clear();
  v_jet_m->clear();
  v_jet_nConst->clear();
  v_jet_truthflav->clear();

  v_jet_GhostL_q->clear();
  v_jet_GhostL_HadI->clear();
  v_jet_GhostL_HadF->clear();
  v_jet_LabDr_HadF->clear();
  v_jet_DoubleHadLabel->clear();
  v_jet_aliveAfterOR->clear();
  v_jet_aliveAfterORmu->clear();
  v_jet_truthMatch->clear();
  v_jet_isBadMedium->clear();
  v_jet_isPU->clear();
  v_jet_truthPt->clear();
  v_jet_dRiso->clear();
  v_jet_JVT->clear();
  v_jet_JVF->clear();
  v_jet_dRminToB->clear();
  v_jet_dRminToC->clear();
  v_jet_dRminToT->clear();

  v_jet_ip2d_pb->clear();
  v_jet_ip2d_pc->clear();
  v_jet_ip2d_pu->clear();
  v_jet_ip2d_llr->clear();

  v_jet_ip3d_pb->clear();
  v_jet_ip3d_pc->clear();
  v_jet_ip3d_pu->clear();
  v_jet_ip3d_llr->clear();

  v_jet_sv1_pb->clear();
  v_jet_sv1_pc->clear();
  v_jet_sv1_pu->clear();
  v_jet_sv1_llr->clear();

  v_jet_dl1_pb->clear();
  v_jet_dl1_pc->clear();
  v_jet_dl1_pu->clear();

  v_jet_dl1mu_pb->clear();
  v_jet_dl1mu_pc->clear();
  v_jet_dl1mu_pu->clear();

  v_jet_dl1rnn_pb->clear();
  v_jet_dl1rnn_pc->clear();
  v_jet_dl1rnn_pu->clear();


  v_jet_sv1ip3d->clear();
  v_jet_mv1->clear();
  v_jet_mv1c->clear();
  v_jet_mv2c00->clear();
  v_jet_mv2c10->clear();
  v_jet_mv2c10mu->clear();
  v_jet_mv2c10rnn->clear();
  v_jet_mv2c20->clear();
  v_jet_mv2c100->clear();
  v_jet_mv2cl100->clear();
  v_jet_mv2m_pu->clear();
  v_jet_mv2m_pb->clear();
  v_jet_mv2m_pc->clear();
  v_jet_mvb->clear();
  v_jet_mv2c20flip->clear();
  v_jet_mv2c10flip->clear();

  v_jet_multisvbb1->clear();
  v_jet_multisvbb2->clear();
  v_jet_msv_N2Tpair->clear();
  v_jet_msv_energyTrkInJet->clear();
  v_jet_msv_nvsec->clear();
  v_jet_msv_normdist->clear();
  v_jet_msv_vtx_cov0->clear();
  v_jet_msv_vtx_cov1->clear();
  v_jet_msv_vtx_cov2->clear();
  v_jet_msv_vtx_cov3->clear();
  v_jet_msv_vtx_cov4->clear();
  v_jet_msv_vtx_cov5->clear();
  v_jet_msv_vtx_mass->clear();
  v_jet_msv_vtx_efrc->clear();
  v_jet_msv_vtx_ntrk->clear();
  v_jet_msv_vtx_pt->clear();
  v_jet_msv_vtx_eta->clear();
  v_jet_msv_vtx_phi->clear();
  v_jet_msv_vtx_dls->clear();
  v_jet_msv_vtx_x->clear();
  v_jet_msv_vtx_y->clear();
  v_jet_msv_vtx_z->clear();
  v_jet_msv_vtx_chi->clear();
  v_jet_msv_vtx_ndf->clear();

  v_jet_ExKtbb_Hbb_DoubleMV2c20->clear();
  v_jet_ExKtbb_Hbb_SingleMV2c20->clear();
  v_jet_ExKtbb_Hbb_MV2Only->clear();
  v_jet_ExKtbb_Hbb_MV2andJFDRSig->clear();
  v_jet_ExKtbb_Hbb_MV2andTopos->clear();

  v_jet_LFCalibType->clear();
  v_jet_btag_ntrk->clear();
  v_jet_trk_dr->clear();
  v_jet_trk_assoc_msv->clear();
  v_jet_trk_algo->clear();
  //v_jet_trk_orig->clear(); moved to BHadronBranches
  v_jet_trk_pdg_id->clear();
  v_jet_trk_parent_pdgid->clear();
  v_jet_trk_barcode->clear();
  v_jet_trk_is_tracking_cp_loose->clear();
  v_jet_trk_d0_truth->clear();
  v_jet_trk_z0_truth->clear();
  v_jet_trk_IP3D_grade->clear();
  v_jet_trk_IP3D_d0->clear();
  v_jet_trk_IP3D_d02D->clear();
  v_jet_trk_IP3D_z0->clear();
  v_jet_trk_IP3D_d0sig->clear();
  v_jet_trk_IP3D_z0sig->clear();

  v_jet_trk_vtx_X->clear();
  v_jet_trk_vtx_Y->clear();
  v_jet_trk_vtx_Z->clear();
  v_jet_trk_vtx_dx->clear();
  v_jet_trk_vtx_dy->clear();

  v_jet_trk_IP3D_llr->clear();

  v_jet_sv1_ntrk->clear();
  v_jet_ip3d_ntrk->clear();
  v_jet_jf_ntrk->clear();

  // MVb variables
  v_jet_width->clear();
  v_jet_n_trk_sigd0cut->clear();
  v_jet_trk3_d0sig->clear();
  v_jet_trk3_z0sig->clear();
  v_jet_sv_scaled_efc->clear();
  v_jet_jf_scaled_efc->clear();

  // additions by andrea
  v_jet_mu_smt->clear();
  v_jet_mu_pt->clear();
  v_jet_mu_eta->clear();
  v_jet_mu_phi->clear();
  v_jet_mu_qOverPratio->clear();
  v_jet_mu_dR->clear();
  v_jet_mu_d0->clear();
  v_jet_mu_z0->clear();
  v_jet_mu_VtxTyp->clear();
  v_jet_mu_mombalsignif->clear();
  v_jet_mu_scatneighsignif->clear();
  v_jet_mu_pTrel->clear();
  v_jet_mu_truthflav->clear();
  v_jet_mu_parent_pdgid->clear();
  v_jet_mu_ID_qOverP_var->clear();
  v_jet_mu_muonType->clear();
  v_jet_mu_assJet_pt->clear();
  // additions by nikola
  v_jet_mu_fatjet_nMu->clear();
  v_jet_mu_fatjet_pTmax_pT->clear();
  v_jet_mu_fatjet_pTmax_pTrel->clear();
  v_jet_mu_fatjet_pTmax_pTrelFrac->clear();
}

int btagIBLAnalysisAlg :: parent_classify(const xAOD::TruthParticle *theParticle) {
  const xAOD::TruthParticle *parent = 0; // the parent object
  Int_t particle_id = 999;
  Int_t parent_id = 999;

  if (theParticle == NULL) return parent_id;

  particle_id = theParticle->pdgId();
  parent = theParticle->parent(0);
  if (parent) parent_id = parent->pdgId();
  else return parent_id;

  while (fabs(parent_id) == fabs(particle_id) && fabs(parent_id) < 400 && fabs(parent_id) != 0) {
    parent = parent->parent(0);
    if (parent) parent_id = parent->pdgId();
    else break;
  }
  return parent_id;
}

void btagIBLAnalysisAlg::fillGhostTracks(const xAOD::Jet& jet,
                                         const xAOD::Vertex& vx) {

  typedef ElementLink<DataVector<xAOD::IParticle> > ParticleLink;
  typedef ElementLink<xAOD::TrackParticleContainer> TrackLink;
  typedef std::vector<ParticleLink> PartLinks;
  typedef std::vector<TrackLink> TrackLinks;
  typedef std::vector<const xAOD::IParticle*> PartVector;
  PartVector loose_tracks;
  for (const ParticleLink& trk: jet.auxdata<PartLinks>("GhostTrack")) {
    const auto* track = dynamic_cast<const xAOD::TrackParticle*>(*trk);
    if (!track) throw std::runtime_error("this isn't a track particle");
    if (m_CPTrackingLooseLabel.empty()) {
      throw std::logic_error("you need the tracking cp tool for loose");
    }
    bool is_cp_loose = m_CPTrackingLooseLabel->accept(*track, &vx);
    if (is_cp_loose) {
      loose_tracks.push_back(track);
    }
  }
  m_ga_track_branches.fill(loose_tracks);
  m_ga_track_cov_branches.fill(loose_tracks);
}


int btagIBLAnalysisAlg::getLFCalibType(const xAOD::Jet * inJet, const xAOD::Vertex* primVtx_ptr) { 
  
  // initialise flags
  int hasKShort     = 0;
  int hasLambda     = 0;
  int hasConversion = 0;
  int hasFake = 0;
  int hasHadMatInt  = 0;

  std::vector<float> j_trk_truthMatchProbability;
  std::vector<float> j_trk_truthDR;
  std::vector<float> j_trk_pt, j_trk_eta, j_trk_d0, j_trk_z0;

  // get BTagging object
  const xAOD::BTagging * tagInfo = inJet->btagging();
  if ( ! tagInfo ) {
    Error("btagIBLAnalysisAlg::getLFCalibType()", "Failed to retrieve BTagging object");
    return 0;
  }

  // // get tracks associated to BTagging object 
  // std::vector< ElementLink< xAOD::TrackParticleContainer > > assocTracks;
  // assocTracks = tagInfo->auxdata<std::vector<ElementLink<xAOD::TrackParticleContainer> > >("BTagTrackToJetAssociator");
  
  // prepare to get tracks   
  std::vector<const xAOD::TrackParticle *> trackParts;     

  // get tracks used in BTagging (SV1 + IP3D, can't get tracks from JetFitter..?)   
  for (size_t i = 0; i < tagInfo->nSV1_TrackParticles(); ++i) {
    const xAOD::TrackParticle * track = tagInfo->SV1_TrackParticle(i);
    if ( ! track ) continue;
    trackParts.push_back( track );
  }
  for (size_t i = 0; i < tagInfo->nIP3D_TrackParticles(); ++i) {
    const xAOD::TrackParticle * track = tagInfo->IP3D_TrackParticle(i);
    if ( ! track ) continue;
    bool skip = false;
    for (const xAOD::TrackParticle * tp : trackParts) {
      if (track == tp) {
        skip = true;   break;
      }
    }
    if ( ! skip ) trackParts.push_back( track );
  }

  //-------------------------------
  // Loop over tracks and get associated truth partices
  //-------------------------------
  for (const xAOD::TrackParticle * track : trackParts) {

    // check if track is valid
    if ( ! track ) continue;

    //Fill other quantities useful for the LF calib:
    j_trk_pt.push_back(track->pt());
    j_trk_eta.push_back(track->eta());
    j_trk_d0.push_back(track->d0());
    j_trk_z0.push_back(track->z0() + track->vz() - primVtx_ptr->z());

    float truthProb = track->auxdata< float >("truthMatchProbability");
    j_trk_truthMatchProbability.push_back(truthProb);
    if ( truthProb < 0.5) hasFake += 1; //consistent with https://cds.cern.ch/record/2266578?
    ////////////////////////

    // get associated truth particle
    typedef ElementLink< xAOD::TruthParticleContainer > Link_t;
    static const char* NAME = "truthParticleLink";
    if ( ! track->isAvailable< Link_t > ( NAME ) ) continue;
    const Link_t& link = track->auxdata< Link_t >( NAME );
    if ( ! link.isValid() ) continue;
    const xAOD::TruthParticle * truth = *link;

    // truth particle pdg id
    int pdgId = truth->pdgId();

    TLorentzVector v_tmp1, v_tmp2;
    v_tmp1.SetPtEtaPhiM( truth->pt(), truth->eta(), truth->phi(), truth->m() );
    v_tmp2.SetPtEtaPhiM( track->pt(), track->eta(), track->phi(), track->m() );
    j_trk_truthDR.push_back(v_tmp1.DeltaR(v_tmp2));

    //-------------------------------
    // check for conversion (gamma -> e+e-)
    //-------------------------------
    if ( ! hasConversion && abs(pdgId) == 11 ) {
      const xAOD::TruthParticle * part = truth;
      while ( part->nParents() == 1 ) {
        part = part->parent();
        int id = part->pdgId();
        if ( id == 22 ) {
          hasConversion += 1;
          break;
        }
        else if ( abs(id) != 11 ) break;
      }
    }

    //-------------------------------
    // check for Lambda
    //-------------------------------
    // lambda decays to proton + pi- 63.9% of the time (and to neutron + pi0 35.8% of the time)
    if ( ! hasLambda && (fabs(pdgId) == 2212 || fabs(pdgId) == 211) ) {
      for (size_t i = 0; i < truth->nParents(); ++i) {
        if (truth->parent(i) && fabs(truth->parent(i)->pdgId()) == 3122) {
          hasLambda += 1;
          break;
        }
      }
    }

    //-------------------------------
    // check for KShort: K_s-> pi+ pi- BR(69.2%) (and -> pi0 + pi0 BR(30.7%))
    // KShort is its own antiparticle, so there is only one PDGID
    //-------------------------------
    if ( ! hasKShort && abs(pdgId) == 211 ) {
      for (size_t i = 0; i < truth->nParents(); ++i) {
        if (truth->parent(i) && truth->parent(i)->pdgId() == 310) {
          hasKShort += 1;
          break;
        }
      }
    }

    //-------------------------------
    // hadronic material interactions
    //-------------------------------
    if ( ! hasHadMatInt ) {

      // geant4 barcode
      if ( truth->barcode() > 200000 ) {
      
        // check if from a KShort / Strange Baryon / Conversion
        int KShort = 0;
        int StrangeBaryon = 0;
        int Conversion = 0;
        int DecayInFlight = 0;
        for (size_t i = 0; i < truth->nParents(); ++i) {
          const xAOD::TruthParticle * part = truth->parent(i);
          int pdg = abs(part->pdgId());
          if (pdg == 310) KShort = 1;
          if (pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 3112 || pdg == 3322 || pdg == 3312 || pdg == 3334) StrangeBaryon = 1; // Lambda/Sigma{0,+,-}/Xi{0,+,-}/Omega{+,-}
          if (pdg == 22) Conversion = 1;
        }
	
        // check if need to look more
        if (Conversion + KShort + StrangeBaryon > 0) {
    
          // check if from decay in flight 
          const xAOD::TruthVertex * truth_vertex = truth->prodVtx();
          if ( ! truth_vertex ) continue;
          TLorentzVector v_out;
          for( size_t iout = 0; iout < truth_vertex->nOutgoingParticles(); iout++) {
            auto *p = truth_vertex->outgoingParticle( iout );
            if ( ! p ) continue;
            TLorentzVector v;
            v.SetPtEtaPhiM( p->pt(), p->eta(), p->phi(), p->m() );
            v_out += v;
          }    
          TLorentzVector v_in;
          for( size_t iin = 0; iin < truth_vertex->nIncomingParticles(); iin++) {
            auto *p = truth_vertex->incomingParticle( iin );
            if ( ! p ) continue;
            TLorentzVector v;
            v.SetPtEtaPhiM( p->pt(), p->eta(), p->phi(), p->m() );
            v_in += v;
          }
          double E_Balance = v_out.E() - v_in.E();
          if (E_Balance < 100.) DecayInFlight = 1;
	  
        }
	
        // set flag for hadronic material interaction
        if (Conversion + KShort + StrangeBaryon + DecayInFlight == 0) hasHadMatInt += 1;

      } 
    }    

  }
  
  // print result
  bool print = false;
  if ( print &&  
       (hasFake && hasKShort + hasLambda + hasConversion + hasHadMatInt) > 0 ) {
    std::cout << "==============  getLFCalibType()  ===============" << std::endl;
    std::cout << "Jet kinematics :" << std::endl;
    std::cout << "---> pt  = " << inJet->pt()  << std::endl;
    std::cout << "---> eta = " << inJet->eta() << std::endl;
    std::cout << "---> phi = " << inJet->phi() << std::endl;
    std::cout << "---> e   = " << inJet->e()   << std::endl;
    std::cout << "hasFake     = " << hasFake     << std::endl;
    std::cout << "hasKShort     = " << hasKShort     << std::endl;
    std::cout << "hasLambda     = " << hasLambda     << std::endl;
    std::cout << "hasConversion = " << hasConversion << std::endl;
    std::cout << "hasHadMatInt  = " << hasHadMatInt  << std::endl;
    std::cout << "================================================" << std::endl;
  }
  
  int outLFCalibType = 0x0;
  if(hasKShort) outLFCalibType |= LFCalibType::hasKShort;
  if(hasLambda) outLFCalibType |= LFCalibType::hasLambda;
  if(hasConversion) outLFCalibType |= LFCalibType::hasConversion;
  if(hasHadMatInt) outLFCalibType |= LFCalibType::hasHadMatInt;
  if(hasFake)  outLFCalibType |= LFCalibType::hasFake;

  //+++++++++++++++++++
  // v_jet_trk_truthMatchProbability->push_back(j_trk_truthMatchProbability);
  // v_jet_trk_truthDR->push_back(j_trk_truthDR);
  // v_jet_trk_pt->push_back(j_trk_pt);
  // v_jet_trk_eta->push_back(j_trk_eta);
  // v_jet_trk_d0->push_back(j_trk_d0);
  // v_jet_trk_z0->push_back(j_trk_z0);
  //+++++++++++++++++++

  return outLFCalibType;
  
}

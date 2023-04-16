///////////////////////////////////////
// btagAnalysisAlg.cxx
///////////////////////////////////////

#include "../btagAnalysis/track_to_vertex_associators.hh"
#include "../btagAnalysis/LifetimeSigning.hh"

#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ServiceHandle.h"

#include "../btagAnalysis/btagAnalysisAlg.h"
#include "../btagAnalysis/GAFlavourLabel.h"
#include "../btagAnalysis/ArbitraryJetBranches.hh"

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
#include "PileupReweighting/PileupReweightingTool.h"

#include "JetInterface/IJetUpdateJvt.h"

// some tracking mumbo jumbo
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <memory>
#include <utility>

using xAOD::IParticle;

// this is the key we use to keep track of how many primary vertices
// we found if it's missing we didn't run the BTagVertexAugmenter.
const std::string VX_COUNT_KEY = "BTaggingNumberOfPrimaryVertices";

bool xaodJetPtSorting(const xAOD::Jet *jet1, const xAOD::Jet *jet2) {
  return jet1->pt() > jet2->pt();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
btagAnalysisAlg::btagAnalysisAlg( const std::string& name, ISvcLocator *pSvcLocator ) :

  AthHistogramAlgorithm(name, pSvcLocator),
  m_stream("BTAGSTREAM"),
  m_jetCleaningTool("JetCleaningTool/JetCleaningTool", this),
  m_jetCalibrationTool(""),
  m_jvt(""),


  m_example_branches(),
  m_jet_properties_branches(),
  m_tagger_scores_branches(),
  m_impact_parameter_branches(),
  m_jetfitter_branches(),
  m_softmuon_branches(),
  m_kshort_branches(),
  m_kshortreco_branches(),
  m_bhadron_branches(),
  m_svx_collections(),
  m_track_branches(),
  m_track_cov_branches(),
  m_arb_branches(0),
  //subjet collection
  m_subjet_collections()
{

  //General Properties
  declareProperty( "Stream", m_stream );
  declareProperty( "JetCollectionName", m_jetCollectionName = "AntiKt4LCTopoJets" );
  declareProperty( "TrackAssociator",m_track_associator = "BTagTrackToJetAssociator");

  //Default handling
  declareProperty( "DefaultValDictionary", m_defaultDict );
  declareProperty( "ReplaceNanDefaults", m_replaceDefaults );
  //Tools
  declareProperty( "JetCleaningTool", m_jetCleaningTool );
  declareProperty( "JVTtool", m_jvt );

  //Flags
  declareProperty( "CalibrateJets", m_calibrateJets = true );
  declareProperty( "CleanJets", m_cleanJets = true );
  declareProperty( "CleanParentJet", m_clean_parent_jet = false );
  declareProperty( "doJVT", m_doJVT = true );

  declareProperty( "EventInfo", m_EventInfo = false );
  declareProperty( "retriveTruthJets", m_retriveTruthJets = false);
  declareProperty( "JetProperties", m_JetProperties = false );
  declareProperty( "TaggerScores", m_TaggerScoresInfo = false );
  declareProperty( "ImpactParameterInfo", m_ImpactParameterInfo = false );
  declareProperty( "SVInfo", m_SVInfo = false );
  declareProperty( "svxCollections", m_svx_collections);
  declareProperty( "JetFitterInfo", m_JetFitterInfo = false );
  declareProperty( "SoftMuoninfo", m_SoftMuoninfo = false );
  declareProperty( "bHadInfo", m_bHadInfo = false );
  declareProperty( "bHadExtraInfo", m_bHadExtraInfo = false);
  declareProperty( "subjetCollections", m_subjet_collections);
  declareProperty( "kshortInfo", m_kshortInfo = false );
  declareProperty( "TrackInfo", m_TrackInfo = false);
  declareProperty( "TrackCovariance", m_TrackCovariance = false );
  declareProperty( "branchDebug", m_showDebug = false );
  declareProperty( "AccessBtagObject", m_access_btag_object = true );
  declareProperty( "nRequiredSiHits", m_n_required_si_hits = 0 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
btagAnalysisAlg::~btagAnalysisAlg() {

  for (auto coll_br: m_subjet_branches) {
    delete coll_br.second;
    coll_br.second = 0;
  }

  if(m_SVInfo){
    for (auto coll: m_svx_branches) {
      delete coll;
      coll = 0;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
StatusCode btagAnalysisAlg::initialize() {
  ATH_MSG_INFO ("Initializing " << name() << "...");

  // Register output tree
  ServiceHandle<ITHistSvc> histSvc("THistSvc",name());

  CHECK( histSvc.retrieve() );

  ATH_MSG_INFO(m_jetCollectionName);

  m_tree = new TTree( ("bTag_" + m_jetCollectionName).c_str(),  ("bTag"  + m_jetCollectionName).c_str() );

  ATH_MSG_INFO ("registering tree in stream: " << m_stream);

  CHECK( histSvc->regTree("/" + m_stream + "/tree_" + m_jetCollectionName, m_tree) );

  CHECK( initializeTools() );


  //EventInfo
  if(m_EventInfo){
    m_tree->Branch("runnb",     &runnumber  );
    m_tree->Branch("eventnb",   &eventnumber);
    m_tree->Branch("mcchan",    &mcchannel  );
    m_tree->Branch("mcwg",      &mcweight   );
    m_tree->Branch("avgmu",     &mu         );
    m_tree->Branch("actmu",     &Act_mu     );
    m_tree->Branch("PVx",       &PV_x       );
    m_tree->Branch("PVy",       &PV_y       );
    m_tree->Branch("PVz",       &PV_z       );
    m_tree->Branch("truth_PVx", &truth_PV_x );
    m_tree->Branch("truth_PVy", &truth_PV_y );
    m_tree->Branch("truth_PVz", &truth_PV_z );
    m_tree->Branch("njets", &njets );
  }

  for (const auto& prefix_name_collection_name: m_subjet_collections) {
    const auto& prefix = prefix_name_collection_name.first;
    const auto& collection = prefix_name_collection_name.second;
    m_subjet_branches.emplace_back(collection, new SubjetBranches);
    m_subjet_branches.back().second->set_tree(*m_tree, prefix, m_showDebug);
  }

  if(m_JetProperties){ m_jet_properties_branches.set_tree(*m_tree, m_showDebug); }

  if(m_TaggerScoresInfo){ m_tagger_scores_branches.set_tree(*m_tree); }

  if(m_ImpactParameterInfo){ m_impact_parameter_branches.set_tree(*m_tree, m_defaultDict, m_replaceDefaults); }

  if(m_JetFitterInfo){ m_jetfitter_branches.set_tree(*m_tree, m_defaultDict, m_replaceDefaults); }

  if(m_SVInfo){
    for (const auto& prefix_name_edm_name: m_svx_collections) {
      const auto& ntuple_prefix = prefix_name_edm_name.first;
      const auto& edm_name = prefix_name_edm_name.second;
      m_svx_branches.emplace_back(new SVBranches(edm_name));
      m_svx_branches.back()->set_tree(*m_tree, ntuple_prefix, m_defaultDict, m_replaceDefaults);
    }
  }

  if(m_SoftMuoninfo){ m_softmuon_branches.set_tree(*m_tree); }

  if(m_bHadInfo){ m_bhadron_branches.set_tree(*m_tree, m_bHadExtraInfo, m_showDebug); }

  if( m_TrackInfo ){ m_track_branches.set_tree(*m_tree, "jet_trk_"); }

  if( m_TrackCovariance ){ m_track_cov_branches.set_tree(*m_tree, "jet_trk_"); }

  // A.X.: create some truth branches here
  mymc_pt = new std::vector<float>();
  mymc_eta = new std::vector<float>();
  mymc_phi = new std::vector<float>();
  mymc_pdgId = new std::vector<int>();
  mymc_decayVtx_x = new std::vector<float>();
  mymc_decayVtx_y = new std::vector<float>();
  mymc_decayVtx_z = new std::vector<float>();
  mymc_ix1 = new std::vector<int>();
  mymc_ix2 = new std::vector<int>();
  mymc1_pt = new std::vector<float>();
  mymc1_eta = new std::vector<float>();
  mymc1_phi = new std::vector<float>();
  mymc1_pdgId = new std::vector<int>();

  mymc_pt->clear();
  mymc_eta->clear();
  mymc_phi->clear();
  mymc_pdgId->clear();
  mymc_decayVtx_x->clear();
  mymc_decayVtx_y->clear();
  mymc_decayVtx_z->clear();
  mymc_ix1->clear();
  mymc_ix2->clear();
  mymc1_pt->clear();
  mymc1_eta->clear();
  mymc1_phi->clear();
  mymc1_pdgId->clear();

  m_tree->Branch("nmmc",&nmymc);
  m_tree->Branch("mymc_pt", &mymc_pt);
  m_tree->Branch("mymc_eta", &mymc_eta);
  m_tree->Branch("mymc_phi", &mymc_phi);
  m_tree->Branch("mymc_pdgId", &mymc_pdgId);
  m_tree->Branch("mymc_decayVtx_x", &mymc_decayVtx_x);
  m_tree->Branch("mymc_decayVtx_y", &mymc_decayVtx_y);
  m_tree->Branch("mymc_decayVtx_z", &mymc_decayVtx_z);
  m_tree->Branch("mymc_ix1", &mymc_ix1);
  m_tree->Branch("mymc_ix2", &mymc_ix2);
  m_tree->Branch("nmmc1",&nmymc1);
  m_tree->Branch("mymc1_pt", &mymc1_pt);
  m_tree->Branch("mymc1_eta", &mymc1_eta);
  m_tree->Branch("mymc1_phi", &mymc1_phi);
  m_tree->Branch("mymc1_pdgId", &mymc1_pdgId);

  return StatusCode::SUCCESS;
}


StatusCode btagAnalysisAlg::initializeTools(){

  CHECK( m_jetCleaningTool.retrieve() );
  if (m_doJVT) {
    CHECK( m_jvt.retrieve() );
  }

  m_jetCalibrationTool.setTypeAndName("JetCalibrationTool/BTagDumpAlg_" + m_jetCollectionName + "_JCalib");

  if (m_calibrateJets) CHECK( m_jetCalibrationTool.retrieve() );

  m_PUtool.setTypeAndName("CP::PileupReweightingTool/prw");

  CHECK( m_PUtool.retrieve() );

  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////
StatusCode btagAnalysisAlg::finalize() {
  ATH_MSG_INFO ("Finalizing " << name() << "...");

  // Write tree into file
  m_tree->Write();

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

StatusCode btagAnalysisAlg::execute() {

  typedef ElementLink<xAOD::TrackParticleContainer> TrackLink;
  typedef std::vector<TrackLink> TrackLinks;

  ATH_MSG_DEBUG ("Executing " << name() << "...");

  // A.X.: TODO: add https://gitlab.cern.ch/MultiBJets/MBJ_Analysis/-/blob/TruthLLP/src/MultibjetsAnalysis/Root/TreeMaker.cxx#L1354

  //-------------------------
  // Event information
  //---------------------------
    const xAOD::EventInfo* eventInfo = 0;

    ATH_MSG_DEBUG(" Retrieving Event Info " );

    CHECK( evtStore()->retrieve(eventInfo) );


    runnumber   = eventInfo->runNumber();
    eventnumber = eventInfo->eventNumber();
    mcchannel   = eventInfo->mcChannelNumber();
    mcweight    = eventInfo->mcEventWeight();
    mu          = eventInfo->averageInteractionsPerCrossing();
    Act_mu      = eventInfo->actualInteractionsPerCrossing();

    float tmpMu = m_PUtool->getCorrectedAverageInteractionsPerCrossing( *eventInfo );

    mu = tmpMu;


  // primary vertex
  const xAOD::VertexContainer *Primary_vertices = 0;

  ATH_MSG_DEBUG( " retrieve  PrimaryVertices " );
  CHECK( evtStore()->retrieve(Primary_vertices, "PrimaryVertices") );

  int* npv_p = 0;

  StatusCode vx_count_status = evtStore()->retrieve(npv_p, VX_COUNT_KEY);

  if (vx_count_status.isFailure()) {
    ATH_MSG_FATAL("could not find " + VX_COUNT_KEY + " in file");
    return StatusCode::FAILURE;
  }

 int npv = *npv_p;

  if (npv < 1) {
    ATH_MSG_WARNING( ".... rejecting the event due to missing PV!!!!");
    return StatusCode::SUCCESS;
  }

  ATH_MSG_DEBUG( " get the vertex index (stored in BTaggingVertexAugmenter) " );
  int* indexPV_ptr = 0;
  CHECK(evtStore()->retrieve(indexPV_ptr, "BTaggingVertexIndex"));

  size_t indexPV = *indexPV_ptr;

  const xAOD::Vertex *myVertex = Primary_vertices->at(indexPV);

  try{
  PV_x = myVertex->x();
  PV_y = myVertex->y();
  PV_z = myVertex->z();
  } catch (...) {
    ATH_MSG_DEBUG( " missing primary vertex! " );
    PV_x = -999;
    PV_y = -999;
    PV_z = -999;
  }

  // const xAOD::TruthVertexContainer *truthVertices = 0;
  // CHECK( evtStore()->retrieve(truthVertices, "TruthVertices") );
  // for ( const auto* truthV : *truthVertices ) {
  //     // record the truth primary vertex position
  //     truth_PV_x = truthV->x();
  //     truth_PV_y = truthV->y();
  //     truth_PV_z = truthV->z();
  //     break;
  // }

  ATH_MSG_DEBUG( " retrieve TruthEvents " );

  // //loop over truth particles, find prompt electrons and muons, and B,C and T partons
  const xAOD::TruthEventContainer *xTruthEventContainer = NULL;
  CHECK( evtStore()->retrieve(xTruthEventContainer, "TruthEvents") );

  std::vector<const xAOD::TruthParticle* > partonB;
  std::vector<const xAOD::TruthParticle* > partonC;
  std::vector<const xAOD::TruthParticle* > partonT;
  std::vector<TLorentzVector> truth_electrons;
  std::vector<TLorentzVector> truth_muons;

  // A.X.: LLP
  std::vector<Float_t> m_LLPLifeTime;
  std::vector<Float_t> m_LLPdecay_x;
  std::vector<Float_t> m_LLPdecay_y;
  std::vector<Float_t> m_LLPdecay_z;

  std::vector<const xAOD::TruthParticle*> LLPs;
  std::map<int, std::vector<const xAOD::TruthParticle*>> LLPDecayProd;
  int nLLP = 0;
  // end A.X.: LLP

  ATH_MSG_DEBUG( " loop over truth particles start " );
  // rewrite the truth primary vertex position
  for ( const auto* truth : *xTruthEventContainer ) {

      const xAOD::TruthVertex *trueSignalVertex = truth->signalProcessVertex();
      if (trueSignalVertex != 0) {

          truth_PV_x = trueSignalVertex->x();
          truth_PV_y = trueSignalVertex->y();
          truth_PV_z = trueSignalVertex->z();
    }

  ATH_MSG_DEBUG( " loop over truth particles in the truth event container " );
  // A.X.
  //std::cout << "A.X.: truth_PV_x = " << truth_PV_x << ", truth_PV_y = " << truth_PV_y << ", truth_PV_z = " << truth_PV_z << std::endl;
  nmymc = nmymc1 = 0;

  for(unsigned int i = 0; i < truth->nTruthParticles(); i++) {
  const xAOD::TruthParticle *particle = truth->truthParticle(i);

  if(!particle){ continue; } // June6,2017, in rel21 derivations some truth particles have been slimmed, and there is a null pointer here.

  // A.X.: detect and save useful MC info: SAMPLE DEPENDENT!!!
  int pdgid = particle->pdgId();
  if (pdgid == 1000022) // gluino 1000021, neutralino 1000022
  //if (pdgid == 35) // SUSY Higgs (H0)
  //if (pdgid == 36) // SUSY Higgs (A0)
  {
    //std::cout << "particle " << i << ": " << particle->pdgId() << " " << particle->status() << std::endl;
    const xAOD::TruthVertex* decayVtx = particle->decayVtx();
    if (decayVtx) {
      //std::cout << "decays at " << decayVtx->x() << " " << decayVtx->y() << " " << decayVtx->z() << " into";
      const unsigned int nOutgoing = ((xAOD::TruthVertex*)decayVtx) -> nOutgoingParticles();
      bool ok = false;
      for ( unsigned int iOut = 0; iOut < nOutgoing; ++iOut) {
	const xAOD::TruthParticle* child = (xAOD::TruthParticle*)decayVtx->outgoingParticle(iOut);
	if(!child) continue;
	//std::cout << " " << child->pdgId();
	if (iOut==0 && child->pdgId()!=particle->pdgId()) {
	  // keep this decay
	  ++nmymc;
	  ok = true;
	  mymc_pt->push_back(particle->pt());
	  mymc_eta->push_back(particle->eta());
	  mymc_phi->push_back(particle->phi());
	  mymc_pdgId->push_back(particle->pdgId());
	  mymc_decayVtx_x->push_back(decayVtx->x());
	  mymc_decayVtx_y->push_back(decayVtx->y());
	  mymc_decayVtx_z->push_back(decayVtx->z());
	  mymc_ix1->push_back(nmymc1);
	}
	if (ok) {
	  // save decay products in a separate list
	  ++nmymc1;
	  mymc1_pt->push_back(child->pt());
	  mymc1_eta->push_back(child->eta());
	  mymc1_phi->push_back(child->phi());
	  mymc1_pdgId->push_back(child->pdgId());
	}
      }
      if (ok) mymc_ix2->push_back(nmymc1);
      //std::cout << std::endl;

      // A.X.: LLP

      const xAOD::TruthParticle* p = particle;

      bool decayFilled = false;
      float lifetimeLab = -1.;
      LLPs.push_back( p );
      
      // Loop over its decay products and save them in a map
      // Sometimes neutralino decays into itself, 
      // so follow decay chain back a few steps to find real decay products
      // This method is not the smartest, but it seems to get the job done.
      for(size_t c = 0; c < p->nChildren(); c++) {
          const xAOD::TruthParticle* child = p->decayVtx()->outgoingParticle(c);
          if (child->pdgId()==1000022){
            for(size_t d = 0; d < child->nChildren(); d++) {
              const xAOD::TruthParticle* child2 = child->decayVtx()->outgoingParticle(d);
              if (child2->pdgId()==1000022){
                for(size_t e = 0; e < child2->nChildren(); e++) {
                  const xAOD::TruthParticle* child3 = child2->decayVtx()->outgoingParticle(e);
                  LLPDecayProd[nLLP].push_back( child3 );
                  if (!decayFilled){
                    m_LLPdecay_x.push_back( child3->prodVtx()->x() );
                    m_LLPdecay_y.push_back( child3->prodVtx()->y() );
                    m_LLPdecay_z.push_back( child3->prodVtx()->z() );
                    lifetimeLab = (child3->prodVtx()->v4().Vect()-trueSignalVertex->v4().Vect()).Mag()/( child2->p4().Beta()*child2->p4().Gamma()*TMath::C()  );
                    m_LLPLifeTime.push_back( lifetimeLab );
                    decayFilled = true;
                  }
                }
              }
              else {
                LLPDecayProd[nLLP].push_back( child2 );
                if (!decayFilled){
                  m_LLPdecay_x.push_back( child2->prodVtx()->x() );
                  m_LLPdecay_y.push_back( child2->prodVtx()->y() );
                  m_LLPdecay_z.push_back( child2->prodVtx()->z() );
                  lifetimeLab = (child2->prodVtx()->v4().Vect()-trueSignalVertex->v4().Vect()).Mag()/( child->p4().Beta()*child->p4().Gamma()*TMath::C()  );
                  m_LLPLifeTime.push_back( lifetimeLab );
                  decayFilled = true;
                }
              }
            }
          }
          else{
            LLPDecayProd[nLLP].push_back( child );
            if (!decayFilled){
              m_LLPdecay_x.push_back( child->prodVtx()->x() );
              m_LLPdecay_y.push_back( child->prodVtx()->y() );
              m_LLPdecay_z.push_back( child->prodVtx()->z() );
              lifetimeLab = (child->prodVtx()->v4().Vect()-trueSignalVertex->v4().Vect()).Mag()/( p->p4().Beta()*p->p4().Gamma()*TMath::C()  );
              m_LLPLifeTime.push_back( lifetimeLab );
              decayFilled = true;
            }
          }
        }
        nLLP++;
      }
      // end A.X.: LLP
  }
  // -- end A.X.

  if (particle->pt() > 3e3) {
    if (fabs(particle->pdgId()) == 15) partonT.push_back(particle);
    if (particle->isCharmHadron()    ) partonC.push_back(particle);
    if (particle->isBottomHadron()   ) partonB.push_back(particle);
  }

  if (particle->pt() < 10e3)     continue;
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

  ATH_MSG_DEBUG( " retrieve "+m_jetCollectionName );

  //---------------------------
  // Jets
  //---------------------------
  const xAOD::JetContainer *jets = 0;
  CHECK( evtStore()->retrieve(jets, m_jetCollectionName) );
  njets = jets->size();
  /// create vector of calibrated jet 4-vectors, for jet isolation calculation
  std::vector<TLorentzVector> calibJets;
  for (const auto* jet : *jets) {
    std::unique_ptr<xAOD::Jet> calib_jet(nullptr);

    if (m_calibrateJets){
      xAOD::Jet* jet_ptr(nullptr);
      m_jetCalibrationTool->calibratedCopy(*jet, jet_ptr);
      calib_jet.reset(jet_ptr);
    } else {
      calib_jet.reset(new xAOD::Jet(*jet));
    }

    TLorentzVector jetvector;
    jetvector.SetPtEtaPhiM(calib_jet->pt(), calib_jet->eta(), calib_jet->phi(), calib_jet->m());
    calibJets.push_back(jetvector);
  }

  const xAOD::JetContainer* truthjets(nullptr);
  if(m_retriveTruthJets){
  ATH_MSG_DEBUG( " retrieve AntiKt4TruthJets " );
  CHECK( evtStore()->retrieve( truthjets, "AntiKt4TruthJets") );
  }

  // A.X.: LLP

  // Loop now over LLP decay products and match them to closest truth jet in event that is within dR < 0.4
  std::map<int, std::vector<const xAOD::Jet*>> matchedTruthJets;

  for (int iLLP = 0; iLLP < nLLP; iLLP++){ // Loop over LLPs

      matchedTruthJets.insert(std::pair<int, std::vector<const xAOD::Jet*>>(iLLP, {}));
      for (const auto &p : LLPDecayProd.at(iLLP)){

        int pID = p->pdgId();
        double pEta = p->eta();
        double pPhi = p->phi();

        float closestJetdR = 10000.;
        int closestJet = 0;
        int jet_index = 0; 
        for(const auto &jet : *truthjets){
          Int_t jLabel = -1;
          jet->getAttribute("HadronConeExclTruthLabelID", jLabel);
          double jEta = jet->eta();
          double jPhi = jet->phi();
          double dRpj = sqrt( (pEta-jEta)*(pEta-jEta) + (pPhi-jPhi)*(pPhi-jPhi) );

          if ( (dRpj < 0.4) && (jLabel<7)){
            // Only match true HF jets with HF, true LF jets with LF
            if (jLabel==5 && !(std::fabs(pID)==5||std::fabs(pID)==6)) continue;
            if (jLabel<5 && !(std::fabs(pID)<5)) continue;
            if (dRpj<closestJetdR){ closestJet = jet_index; closestJetdR = dRpj;}
          }
          jet_index += 1;

        } // Truth jets
        if (closestJetdR<0.4){
          matchedTruthJets.at(iLLP).push_back(truthjets->at(closestJet));
        }
    }// LLPs decay products
  } // LLPs

  // end A.X.: LLP

  /////////////////////////////////////////////////////////////////////////
  // MAIN JET LOOP
  int jet_index = -1;

  ATH_MSG_DEBUG( " Main jet loop " );
  for (const auto* jet : *jets) {
    jet_index++;

    std::unique_ptr<xAOD::Jet> calib_jet(nullptr);

    if (m_calibrateJets){
      xAOD::Jet* jet_ptr(nullptr);
      m_jetCalibrationTool->calibratedCopy(*jet, jet_ptr);
      calib_jet.reset(jet_ptr);
    } else {
      calib_jet.reset(new xAOD::Jet(*jet));
    }

    int isBad = 0;
    // clean jets
    if (m_cleanJets) {
      const xAOD::Jet *jet_to_clean = calib_jet.get();
      isBad = (!m_jetCleaningTool->keep(*jet_to_clean));
                           // additions by nikola
      if (m_clean_parent_jet) {
          jet_to_clean = GetParentJet(jet, "Parent");
      }
    }

    // JVT
    float jvtV = NAN;

    if(m_doJVT){
        jvtV = calib_jet->auxdata<float>("Jvt");
        float tmpJVT = jvtV;
        jvtV = m_jvt->updateJvt(*calib_jet);
        if (tmpJVT != jvtV)  ATH_MSG_DEBUG(" initial: " << tmpJVT << " |  final: " << jvtV );
    }



    // Jet isolation
    float dRiso = 10;
    for (unsigned int j = 0; j < calibJets.size(); j++) {

      if (static_cast<int>(j) == jet_index) continue;
      TLorentzVector jet2 = calibJets.at(j);
      float dr = calib_jet->p4().DeltaR(jet2);
      if (dr < dRiso) dRiso = dr;

    }

    // A.X.: LLP

    // Try to match the jet with one of the truth jets that matched to a truth particle coming from the LLP
    float closestTruthJet = 100.;
    int LLPmatch = -1;
    Int_t truth_LLP_jetLabel = -1;
    float truth_LLP_decay_vertex[3];
    for (int iLLP = 0; iLLP < nLLP; iLLP++){

      for(auto truth_jet : matchedTruthJets.at(iLLP)){
        const auto &pLLP = LLPs.at(iLLP);
        bool matchRecoTruthJets = false;

        const xAOD::TruthVertex *trueSignalVertex = pLLP->prodVtx();
	if (trueSignalVertex==0) continue; // A.X.: fixed
        TVector3 truthVertex( trueSignalVertex->x(), trueSignalVertex->y(), trueSignalVertex->z() );
        TVector3 LLPDecayVertex( m_LLPdecay_x.at(iLLP), m_LLPdecay_y.at(iLLP), m_LLPdecay_z.at(iLLP) );

        double r_decay = ( LLPDecayVertex - truthVertex ).Mag();
        double r_calo = 1800.0; // Units in mm

        double truth_jet_theta = 2*atan(exp(-truth_jet->eta()));
        double pLLP_theta = 2*atan(exp(-pLLP->eta()));

        double dR_reco_truth = 1000.;
        if (r_decay*sin(pLLP_theta)>r_calo){
          // If LLP decays beyond start of calorimeter, don't apply correction
          dR_reco_truth = sqrt( pow( (jet->eta()-truth_jet->eta()), 2) + pow( (jet->phi()-truth_jet->phi()), 2 ) );
        }
        else { 
          double a = pow(sin(truth_jet_theta),2.);
          double b = 2.0 * r_decay * sin(pLLP_theta) * sin(truth_jet_theta) * ( cos(pLLP->phi())*cos(truth_jet->phi()) + sin(pLLP->phi())*sin(truth_jet->phi()) );
          double c = pow(r_decay*sin(pLLP_theta), 2) - r_calo*r_calo;

          double d_traveled = (-1.* b + sqrt(b*b - 4. * a * c))/(2. * a);

          double x_cal = r_decay * sin(pLLP_theta) * cos(pLLP->phi()) + d_traveled * sin(truth_jet_theta) * cos(truth_jet->phi());
          double y_cal = r_decay * sin(pLLP_theta) * sin(pLLP->phi()) + d_traveled * sin(truth_jet_theta) * sin(truth_jet->phi());
          double z_cal = r_decay * cos(pLLP_theta) + d_traveled * cos(truth_jet_theta);

          double new_truthJet_theta = acos(z_cal/(sqrt(x_cal*x_cal + y_cal*y_cal + z_cal*z_cal)));
          double new_truthJet_phi = atan2(y_cal, x_cal);
          double new_truthJet_eta = -log(std::fabs(tan(new_truthJet_theta/2.)));

          dR_reco_truth = sqrt( pow( (jet->eta()-new_truthJet_eta), 2) + pow( (jet->phi()-new_truthJet_phi), 2 ) );
        }

        if (dR_reco_truth < 0.3 && dR_reco_truth < closestTruthJet){
           matchRecoTruthJets = true;
           LLPmatch = iLLP;
           closestTruthJet = dR_reco_truth;
           truth_jet->getAttribute("HadronConeExclTruthLabelID", truth_LLP_jetLabel);
	   truth_LLP_decay_vertex[0] = m_LLPdecay_x.at(iLLP);
	   truth_LLP_decay_vertex[1] = m_LLPdecay_y.at(iLLP);
	   truth_LLP_decay_vertex[2] = m_LLPdecay_z.at(iLLP);
        }
      } // truth jets

    } // LLP 
    /*m_isLLPDecayProd.push_back(LLPmatch);
    if (closestTruthJet<0.4) m_dRtoLLPDecayProd.push_back(closestTruthJet);
    else m_dRtoLLPDecayProd.push_back(-1.);
    m_truthLLPJetLabel.push_back(truth_LLP_jetLabel);*/

    // end A.X.: LLP

    //Fill Tree Branches
    const xAOD::BTagging *bjet(nullptr);

    if(m_access_btag_object){
      bjet = calib_jet->btagging();
    }

    if(m_JetProperties){
      ATH_MSG_DEBUG( " filling JetProperties " );
      m_jet_properties_branches.fill(*jet, *calib_jet, jvtV, dRiso, isBad,
	LLPmatch, closestTruthJet, truth_LLP_jetLabel, truth_LLP_decay_vertex,
	truth_electrons, truth_muons,  truthjets, partonB, partonC, partonT);
    }

    if(m_TaggerScoresInfo){ ATH_MSG_DEBUG( " filling TaggerScores " );
      m_tagger_scores_branches.fill(*calib_jet); }

    if(m_ImpactParameterInfo){  ATH_MSG_DEBUG( " filling impact_parameter_branches " );
      m_impact_parameter_branches.fill(*calib_jet); }

    if(m_JetFitterInfo){ ATH_MSG_DEBUG( " filling JetFitterInfo " );
      m_jetfitter_branches.fill(*calib_jet); }

    if(m_SVInfo){
      for (auto* svx: m_svx_branches) {
        ATH_MSG_DEBUG( " filling SV info " );
        svx->fill(*bjet);
      }
    }

    if(m_SoftMuoninfo) { ATH_MSG_DEBUG( " filling Soft Muon info " ); m_softmuon_branches.fill(*calib_jet); }

    if(m_bHadInfo){  ATH_MSG_DEBUG( " filling b/c hadron info " ); m_bhadron_branches.fill(*calib_jet, m_jetCollectionName); }


    if(m_TrackInfo){

      typedef std::vector<const xAOD::IParticle*> Particles;

      typedef ElementLink<xAOD::TrackParticleContainer> TrackLink;
      typedef std::vector<TrackLink> TrackLinks;

      TrackLinks assocTracks = bjet->auxdata<TrackLinks>(m_track_associator);

      ATH_MSG_DEBUG( " starting track selection loop " );

      std::vector<const xAOD::TrackParticle*> selectedTracks;

      uint8_t getInt(0);

      for (unsigned int iT = 0; iT < assocTracks.size(); iT++) {
        // std::cout << " .... trk link: " << iT << std::endl;
        if (!assocTracks.at(iT).isValid()) continue;
        const xAOD::TrackParticle *tmpTrk = *(assocTracks.at(iT));
        tmpTrk->summaryValue(getInt, xAOD::numberOfPixelHits);
        int nSi = getInt;
        tmpTrk->summaryValue(getInt, xAOD::numberOfSCTHits);
        nSi += getInt;
        if (nSi < m_n_required_si_hits) continue;
        selectedTracks.push_back(tmpTrk);
      }

       ATH_MSG_DEBUG( " filling track branches " );

       Particles track_particles(selectedTracks.begin(), selectedTracks.end());

       m_track_branches.fill(track_particles, *bjet, *jet);

       if(m_TrackCovariance){
        m_track_cov_branches.fill(track_particles);
       }

       ATH_MSG_DEBUG( "deleting selected tracks" );

    }

    // fill in the subjets branches when running on the fatjet
    for (const auto& collection_branches: m_subjet_branches) {
      ATH_MSG_DEBUG( " filling subjet info " );
      const auto& name = collection_branches.first;
      auto subjet = calib_jet->getAssociatedObjects<xAOD::Jet>(name);
      std::sort(subjet.begin(), subjet.end(),
                [](auto* j1, auto* j2) { return j1->pt() > j2->pt();});
      ATH_MSG_DEBUG( " filling subjet info "<< name <<
                     " subjets size "<<subjet.size() );
      collection_branches.second->fill(subjet);
    }

  } // end jet loop



  m_tree->Fill();

  if(m_JetProperties){ m_jet_properties_branches.clear(); }
  if(m_TaggerScoresInfo){ m_tagger_scores_branches.clear(); }
  if(m_ImpactParameterInfo){ m_impact_parameter_branches.clear(); }
  if(m_JetFitterInfo){ m_jetfitter_branches.clear(); }
  if(m_SVInfo){
    for (auto& coll: m_svx_branches) {
      coll->clear();
    }
  }
  for (auto coll_br: m_subjet_branches) {
    coll_br.second->clear();
  }
  if(m_SoftMuoninfo){ m_softmuon_branches.clear();}
  if(m_bHadInfo){ m_bhadron_branches.clear(); }
  if(m_TrackInfo){ m_track_branches.clear(); }
  if(m_TrackCovariance){ m_track_cov_branches.clear(); }

  // A.X.: clear my vectors after the event is filled
  mymc_pt->clear();
  mymc_eta->clear();
  mymc_phi->clear();
  mymc_pdgId->clear();
  mymc_decayVtx_x->clear();
  mymc_decayVtx_y->clear();
  mymc_decayVtx_z->clear();
  mymc_ix1->clear();
  mymc_ix2->clear();
  mymc1_pt->clear();
  mymc1_eta->clear();
  mymc1_phi->clear();
  mymc1_pdgId->clear();

  return StatusCode::SUCCESS;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// User defined functions
/////////////////////////

const xAOD::Jet *btagAnalysisAlg::GetParentJet(const xAOD::Jet *Jet,
                                               std::string Keyname) {
  typedef ElementLink<xAOD::JetContainer> Jets;
  Jets el = Jet->auxdata<Jets>(Keyname);
  if(el.isValid()) {
    return *el;
  }
  else {
    ATH_MSG_WARNING(
      "GetParentJet(): Unable to get parent link " <<
      Keyname << " ! Null ptr is returned.");
    return 0;
  }
}




bool btagAnalysisAlg::isFromWZ( const xAOD::TruthParticle* particle ) {

  if ( particle==0 ) return false;

  if ( fabs(particle->pdgId())!= 11 && fabs(particle->pdgId())!= 13) return false;

  const xAOD::TruthVertex* prodvtx = particle->prodVtx();
  if ( prodvtx==0 ){ ATH_MSG_DEBUG(" isFromWZ : particle without prodvtx "); return false; }

  ATH_MSG_DEBUG(" isFromWZ : particle has prodVtx ");

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

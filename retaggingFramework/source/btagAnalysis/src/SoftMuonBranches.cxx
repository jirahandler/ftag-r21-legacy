#include "../btagAnalysis/SoftMuonBranches.hh"
#include "../btagAnalysis/SoftMuonBranchBuffer.hh"

#include "xAODJet/Jet.h"
#include "AthContainers/exceptions.h"
#include "TTree.h"



//!-----------------------------------------------------------------------------------------------------------------------------!//
SoftMuonBranches::SoftMuonBranches():
  m_branches(new SoftMuonBranchBuffer)
{
  // instantiate all the vectors here ...
  m_branches->v_jet_mu_smt = new std::vector<double>();
  m_branches->v_jet_mu_pt = new std::vector<float>();
  m_branches->v_jet_mu_eta = new std::vector<float>();
  m_branches->v_jet_mu_phi = new std::vector<float>();
  m_branches->v_jet_mu_qOverPratio = new std::vector<float>();
  m_branches->v_jet_mu_dR = new std::vector<float>();
  m_branches->v_jet_mu_d0 = new std::vector<float>();
  m_branches->v_jet_mu_z0 = new std::vector<float>();
  m_branches->v_jet_mu_VtxTyp = new std::vector<float>();
  m_branches->v_jet_mu_mombalsignif = new std::vector<float>();
  m_branches->v_jet_mu_scatneighsignif = new std::vector<float>();
  m_branches->v_jet_mu_pTrel = new std::vector<float>();
  m_branches->v_jet_mu_parent_pdgid = new std::vector<float>();
  m_branches->v_jet_mu_ID_qOverP_var = new std::vector<float>();
  m_branches->v_jet_mu_muonType = new std::vector<float>();



}

//!-----------------------------------------------------------------------------------------------------------------------------!//
SoftMuonBranches::~SoftMuonBranches() {
  // delete all the vectors here ...
  delete m_branches->v_jet_mu_smt;
  delete m_branches->v_jet_mu_pt;
  delete m_branches->v_jet_mu_eta;
  delete m_branches->v_jet_mu_phi;
  delete m_branches->v_jet_mu_qOverPratio;
  delete m_branches->v_jet_mu_dR;
  delete m_branches->v_jet_mu_d0;
  delete m_branches->v_jet_mu_z0;
  delete m_branches->v_jet_mu_VtxTyp;
  delete m_branches->v_jet_mu_mombalsignif;
  delete m_branches->v_jet_mu_scatneighsignif;
  delete m_branches->v_jet_mu_pTrel;
  delete m_branches->v_jet_mu_parent_pdgid;
  delete m_branches->v_jet_mu_ID_qOverP_var;
  delete m_branches->v_jet_mu_muonType;


  delete m_branches;
}

void SoftMuonBranches::set_tree(TTree& output_tree) const {

    output_tree.Branch("jet_mu_smt", &m_branches->v_jet_mu_smt);
    output_tree.Branch("jet_mu_dR", &m_branches->v_jet_mu_dR);
    output_tree.Branch("jet_mu_pTrel", &m_branches->v_jet_mu_pTrel);
    output_tree.Branch("jet_mu_qOverPratio", &m_branches->v_jet_mu_qOverPratio);
    output_tree.Branch("jet_mu_mombalsignif", &m_branches->v_jet_mu_mombalsignif);
    output_tree.Branch("jet_mu_scatneighsignif", &m_branches->v_jet_mu_scatneighsignif);
    output_tree.Branch("jet_mu_VtxTyp", &m_branches->v_jet_mu_VtxTyp);
    output_tree.Branch("jet_mu_pt", &m_branches->v_jet_mu_pt);
    output_tree.Branch("jet_mu_eta", &m_branches->v_jet_mu_eta);
    output_tree.Branch("jet_mu_phi", &m_branches->v_jet_mu_phi);
    output_tree.Branch("jet_mu_d0", &m_branches->v_jet_mu_d0);
    output_tree.Branch("jet_mu_z0", &m_branches->v_jet_mu_z0);
    output_tree.Branch("jet_mu_parent_pdgid", &m_branches->v_jet_mu_parent_pdgid);
    output_tree.Branch("jet_mu_ID_qOverP_var", &m_branches->v_jet_mu_ID_qOverP_var);
    output_tree.Branch("jet_mu_muonType", &m_branches->v_jet_mu_muonType);



}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void SoftMuonBranches::fill(const xAOD::Jet& jet) {

  const xAOD::BTagging *bjet = jet.btagging();

   double jet_mu_dRmin_smt      = NAN ;
   float jet_mu_dRmin_pt        = NAN ,
   jet_mu_dRmin_dR              = NAN ,
   jet_mu_dRmin_eta             = NAN ,
   jet_mu_dRmin_phi             = NAN ,
   jet_mu_dRmin_qOverPratio     = NAN ,
   jet_mu_dRmin_mombalsignif    = NAN ,
   jet_mu_dRmin_scatneighsignif = NAN ,
   jet_mu_dRmin_pTrel           = NAN ,
   jet_mu_dRmin_VtxTyp          = NAN ,
   jet_mu_dRmin_d0              = NAN ,
   jet_mu_dRmin_z0              = NAN ,
   jet_mu_dRmin_parent_pdgid    = NAN ,
   jet_mu_dRmin_ID_qOverP_var   = NAN ,
   jet_mu_dRmin_muonType        = NAN ;


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


  ElementLink<xAOD::MuonContainer> tmpMuonLink= bjet->auxdata<ElementLink<xAOD::MuonContainer> >("SMT_mu_link");
  if ( tmpMuonLink.isValid() ) {
    const xAOD::Muon* tmpMuon=(*tmpMuonLink);

    if ( tmpMuon!=0 ) {
      jet_mu_dRmin_eta      =tmpMuon->eta();
      jet_mu_dRmin_phi      =tmpMuon->phi();
      jet_mu_dRmin_muonType =tmpMuon->muonType();

      const ElementLink< xAOD::TrackParticleContainer >& pMuIDTrack=tmpMuon->inDetTrackParticleLink();

      const xAOD::Vertex * pVtx=(*pMuIDTrack)->vertex();

      if(pVtx!=NULL) {
        jet_mu_dRmin_VtxTyp=pVtx->vertexType();
      } else { jet_mu_dRmin_VtxTyp=999.;}

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


      m_branches->v_jet_mu_smt->push_back(jet_mu_dRmin_smt);
      m_branches->v_jet_mu_pt->push_back(jet_mu_dRmin_pt);
      m_branches->v_jet_mu_eta->push_back(jet_mu_dRmin_eta);
      m_branches->v_jet_mu_phi->push_back(jet_mu_dRmin_phi);
      m_branches->v_jet_mu_qOverPratio->push_back(jet_mu_dRmin_qOverPratio);
      m_branches->v_jet_mu_mombalsignif->push_back(jet_mu_dRmin_mombalsignif);
      m_branches->v_jet_mu_scatneighsignif->push_back(jet_mu_dRmin_scatneighsignif);
      m_branches->v_jet_mu_dR->push_back(jet_mu_dRmin_dR);
      m_branches->v_jet_mu_pTrel->push_back(jet_mu_dRmin_pTrel);
      m_branches->v_jet_mu_VtxTyp->push_back(jet_mu_dRmin_VtxTyp);
      m_branches->v_jet_mu_d0->push_back(jet_mu_dRmin_d0);
      m_branches->v_jet_mu_z0->push_back(jet_mu_dRmin_z0);
      m_branches->v_jet_mu_parent_pdgid->push_back(jet_mu_dRmin_parent_pdgid);
      m_branches->v_jet_mu_ID_qOverP_var->push_back(jet_mu_dRmin_ID_qOverP_var);
      m_branches->v_jet_mu_muonType->push_back(jet_mu_dRmin_muonType);


}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void SoftMuonBranches::clear() {
  // clear vectors

  m_branches->v_jet_mu_smt->clear();
  m_branches->v_jet_mu_pt->clear();
  m_branches->v_jet_mu_eta->clear();
  m_branches->v_jet_mu_phi->clear();
  m_branches->v_jet_mu_qOverPratio->clear();
  m_branches->v_jet_mu_dR->clear();
  m_branches->v_jet_mu_d0->clear();
  m_branches->v_jet_mu_z0->clear();
  m_branches->v_jet_mu_VtxTyp->clear();
  m_branches->v_jet_mu_mombalsignif->clear();
  m_branches->v_jet_mu_scatneighsignif->clear();
  m_branches->v_jet_mu_pTrel->clear();
  m_branches->v_jet_mu_parent_pdgid->clear();
  m_branches->v_jet_mu_ID_qOverP_var->clear();
  m_branches->v_jet_mu_muonType->clear();



}



int SoftMuonBranches :: parent_classify(const xAOD::TruthParticle *theParticle) {
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




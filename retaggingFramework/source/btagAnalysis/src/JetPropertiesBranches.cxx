#include "../btagAnalysis/JetPropertiesBranches.hh"
#include "../btagAnalysis/JetPropertiesBranchBuffer.hh"

#include "AthContainers/exceptions.h"
#include "TTree.h"


//!-----------------------------------------------------------------------------------------------------------------------------!//
JetPropertiesBranches::JetPropertiesBranches():
  m_branches(new JetPropertiesBranchBuffer)
{
  // instantiate all the vectors here ...
  m_branches->v_jet_pt             = new std::vector<float>();
  m_branches->v_jet_eta            = new std::vector<float>();
  m_branches->v_jet_phi            = new std::vector<float>();
  m_branches->v_jet_E              = new std::vector<float>();
  m_branches->v_jet_pt_orig        = new std::vector<float>();
  m_branches->v_jet_eta_orig       = new std::vector<float>();
  m_branches->v_jet_phi_orig       = new std::vector<float>();
  m_branches->v_jet_E_orig         = new std::vector<float>();
  m_branches->v_jet_LabDr_HadF     = new std::vector<int>();
  m_branches->v_jet_DoubleHadLabel = new std::vector<int>();
  m_branches->v_jet_JVT            = new std::vector<float>();

  m_branches->v_jet_m              = new std::vector<float>();
  m_branches->v_jet_nConst         = new std::vector<float>();
  m_branches->v_jet_dRiso          = new std::vector<float>();
  m_branches->v_jet_truthMatch     = new std::vector<int>();
  m_branches->v_jet_isPU           = new std::vector<int>();
  m_branches->v_jet_aliveAfterOR   = new std::vector<int>();
  m_branches->v_jet_aliveAfterORmu = new std::vector<int>();
  m_branches->v_jet_isBadMedium    = new std::vector<int>();
  m_branches->v_jet_truthPt        = new std::vector<float>();

  m_branches->v_jet_dRminToB  = new std::vector<float>();
  m_branches->v_jet_dRminToC  = new std::vector<float>();
  m_branches->v_jet_dRminToT  = new std::vector<float>();

  m_branches->v_jet_isLLPDecayProd   = new std::vector<int>();
  m_branches->v_jet_dRtoLLPDecayProd = new std::vector<float>();
  m_branches->v_jet_truthLLPJetLabel = new std::vector<int>();
  m_branches->v_jet_truthLLP_Decay_x = new std::vector<float>();
  m_branches->v_jet_truthLLP_Decay_y = new std::vector<float>();
  m_branches->v_jet_truthLLP_Decay_z = new std::vector<float>();

}

//!-----------------------------------------------------------------------------------------------------------------------------!//
JetPropertiesBranches::~JetPropertiesBranches() {
  // delete all the vectors here ...
  delete m_branches->v_jet_pt            ;
  delete m_branches->v_jet_eta           ;
  delete m_branches->v_jet_phi           ;
  delete m_branches->v_jet_E             ;
  delete m_branches->v_jet_pt_orig       ;
  delete m_branches->v_jet_eta_orig      ;
  delete m_branches->v_jet_phi_orig      ;
  delete m_branches->v_jet_E_orig        ;
  delete m_branches->v_jet_LabDr_HadF    ;
  delete m_branches->v_jet_DoubleHadLabel;
  delete m_branches->v_jet_JVT           ;
  delete m_branches->v_jet_m             ;
  delete m_branches->v_jet_nConst        ;
  delete m_branches->v_jet_dRiso         ;
  delete m_branches->v_jet_truthMatch    ;
  delete m_branches->v_jet_isPU          ;
  delete m_branches->v_jet_aliveAfterOR  ;
  delete m_branches->v_jet_aliveAfterORmu;
  delete m_branches->v_jet_isBadMedium   ;
  delete m_branches->v_jet_truthPt       ;
  delete m_branches->v_jet_dRminToB      ;
  delete m_branches->v_jet_dRminToC      ;
  delete m_branches->v_jet_dRminToT      ;

  delete m_branches->v_jet_isLLPDecayProd;
  delete m_branches->v_jet_dRtoLLPDecayProd;
  delete m_branches->v_jet_truthLLPJetLabel;
  delete m_branches->v_jet_truthLLP_Decay_x;
  delete m_branches->v_jet_truthLLP_Decay_y;
  delete m_branches->v_jet_truthLLP_Decay_z;
  delete m_branches;
}

void JetPropertiesBranches::set_tree(TTree& output_tree, bool show_debug){

  debug = show_debug;

  output_tree.Branch( "jet_pt" , &m_branches->v_jet_pt    );
  output_tree.Branch( "jet_eta" , &m_branches->v_jet_eta  );
  output_tree.Branch( "jet_phi" , &m_branches->v_jet_phi  );
  output_tree.Branch( "jet_E" , &m_branches->v_jet_E );
  output_tree.Branch( "jet_pt_orig" , &m_branches->v_jet_pt_orig );
  output_tree.Branch( "jet_eta_orig" , &m_branches->v_jet_eta_orig);
  output_tree.Branch( "jet_phi_orig" , &m_branches->v_jet_phi_orig);
  output_tree.Branch( "jet_E_orig" , &m_branches->v_jet_E_orig  );
  output_tree.Branch( "jet_LabDr_HadF" , &m_branches->v_jet_LabDr_HadF    );
  output_tree.Branch( "jet_DoubleHadLabel" , &m_branches->v_jet_DoubleHadLabel);
  output_tree.Branch( "jet_JVT" , &m_branches->v_jet_JVT     );
  output_tree.Branch( "jet_m" , &m_branches->v_jet_m );
  output_tree.Branch( "jet_nConst" , &m_branches->v_jet_nConst );
  output_tree.Branch( "jet_dRiso" , &m_branches->v_jet_dRiso );
  output_tree.Branch( "jet_truthMatch" , &m_branches->v_jet_truthMatch );
  output_tree.Branch( "jet_isPU" , &m_branches->v_jet_isPU );
  output_tree.Branch( "jet_aliveAfterOR" , &m_branches->v_jet_aliveAfterOR );
  output_tree.Branch( "jet_aliveAfterORmu" , &m_branches->v_jet_aliveAfterORmu );
  output_tree.Branch( "jet_isBadMedium" , &m_branches->v_jet_isBadMedium );
  output_tree.Branch( "jet_truthPt" , &m_branches->v_jet_truthPt );
  output_tree.Branch( "jet_dRminToB" , &m_branches->v_jet_dRminToB );
  output_tree.Branch( "jet_dRminToC" , &m_branches->v_jet_dRminToC );
  output_tree.Branch( "jet_dRminToT" , &m_branches->v_jet_dRminToT );

  output_tree.Branch( "jet_isLLPDecayProd" , &m_branches->v_jet_isLLPDecayProd );
  output_tree.Branch( "jet_dRtoLLPDecayProd" , &m_branches->v_jet_dRtoLLPDecayProd );
  output_tree.Branch( "jet_truthLLPJetLabel" , &m_branches->v_jet_truthLLPJetLabel );
  output_tree.Branch( "jet_truthLLP_Decay_x" , &m_branches->v_jet_truthLLP_Decay_x );
  output_tree.Branch( "jet_truthLLP_Decay_y" , &m_branches->v_jet_truthLLP_Decay_y );
  output_tree.Branch( "jet_truthLLP_Decay_z" , &m_branches->v_jet_truthLLP_Decay_z );

}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void JetPropertiesBranches::fill(const xAOD::Jet& orig_jet, const xAOD::Jet& calib_jet, float JVT, float dRiso, int isBad,
  int LLPmatch, float closestTruthJet, int truth_LLP_jetLabel, float* truth_LLP_decay_vertex,
  std::vector<TLorentzVector> truth_electrons, std::vector<TLorentzVector> truth_muons,
    const xAOD::JetContainer *truthjets,
    std::vector<const xAOD::TruthParticle* > m_partonB,
    std::vector<const xAOD::TruthParticle* > m_partonC,
    std::vector<const xAOD::TruthParticle* > m_partonT) {

  if(debug){ std::cout << " jet kinematics " << std::endl; }
  m_branches->v_jet_pt->push_back(calib_jet.pt());
  m_branches->v_jet_eta->push_back(calib_jet.eta());
  m_branches->v_jet_phi->push_back(calib_jet.phi());
  m_branches->v_jet_E->push_back(calib_jet.e());
  m_branches->v_jet_pt_orig->push_back(orig_jet.pt());
  m_branches->v_jet_eta_orig->push_back(orig_jet.eta());
  m_branches->v_jet_phi_orig->push_back(orig_jet.phi());
  m_branches->v_jet_E_orig->push_back(orig_jet.e());

  m_branches->v_jet_m->push_back( calib_jet.m() );
  m_branches->v_jet_nConst->push_back( calib_jet.numConstituents() );
  m_branches->v_jet_dRiso->push_back( dRiso );

  m_branches->v_jet_JVT->push_back(JVT);
  m_branches->v_jet_isBadMedium->push_back(isBad);

  if(debug){ std::cout << " jet truth labels " << std::endl; }
  int tmpLabel = -1;
    try {
      calib_jet.getAttribute("HadronConeExclTruthLabelID", tmpLabel);
    } catch(...) {};

  int tmpLabelDoubleHadron = -1;
    try {
      calib_jet.getAttribute("HadronConeExclExtendedTruthLabelID", tmpLabelDoubleHadron);
    } catch(...) {};

  m_branches->v_jet_LabDr_HadF->push_back(tmpLabel);
  m_branches->v_jet_DoubleHadLabel->push_back(tmpLabelDoubleHadron);



    // matching reco jets to truth jets and recording the truth jet pT
    // picking the highest pT truth jet (with pT > 7GeV) that satisfies dR < 0.3
    // N.B. this assumes that truth jets are pT ordered
    int matchedPt = 0;
    float dRmatch = 100;
    bool truthFree = true;
    if(debug){ std::cout << " truth jets  " << (truthjets ? "valid" : "empty" )<< std::endl; }
    if(truthjets){
    for (const auto* tjet : *truthjets) {

      if (tjet->pt() < 10e3) continue;
      float dr = calib_jet.p4().DeltaR(tjet->p4());

        if (dr < 0.6) {
          truthFree = false;
        }

      if (dr < 0.3 && dr < dRmatch) {
              dRmatch = dr;
              matchedPt = tjet->pt();
      }
    }
    }

    if (dRmatch < 0.3) {
        m_branches->v_jet_truthMatch->push_back(1);
        m_branches->v_jet_truthPt->push_back(matchedPt);
    }
    else {
        m_branches->v_jet_truthMatch->push_back(0);
        m_branches->v_jet_truthPt->push_back(0);
    }

     m_branches->v_jet_isPU->push_back(truthFree);


    //////////////////////////////////////////////////////////////////
    // flagging jets that overlap with electron
    if(debug){ std::cout << " filling overlap flags " << std::endl; }

    bool iseljetoverlap = false;
    for(unsigned int i = 0; i < truth_electrons.size(); i++) {
      float dr = calib_jet.p4().DeltaR(truth_electrons.at(i));
      if (dr < 0.3){ iseljetoverlap = true;
        break; }
    }
    m_branches->v_jet_aliveAfterOR->push_back( !iseljetoverlap );

    iseljetoverlap = false;
    for(unsigned int i= 0; i < truth_muons.size(); i++){
      float dr = calib_jet.p4().DeltaR(truth_muons.at(i));
      if(dr < 0.3){ iseljetoverlap = true;
        break;}
    }
    m_branches->v_jet_aliveAfterORmu->push_back( !iseljetoverlap );



    float mindRtoB = 10;
    float mindRtoC = 10;
    float mindRtoT = 10;

    if(debug){ std::cout << " jet purity flags " << std::endl; }

    for (unsigned int ip = 0; ip < m_partonB.size(); ip++) {
      float dr = calib_jet.p4().DeltaR(m_partonB.at(ip)->p4());
      if (dr < mindRtoB) mindRtoB = dr;
    }
    for (unsigned int ip = 0; ip < m_partonC.size(); ip++) {

      float dr = calib_jet.p4().DeltaR(m_partonC.at(ip)->p4());
      if (dr < mindRtoC) mindRtoC = dr;
    }
    for (unsigned int ip = 0; ip < m_partonT.size(); ip++) {
      float dr = calib_jet.p4().DeltaR(m_partonT.at(ip)->p4());
      if (dr < mindRtoT) mindRtoT = dr;
    }
    m_branches->v_jet_dRminToB->push_back(mindRtoB);
    m_branches->v_jet_dRminToC->push_back(mindRtoC);
    m_branches->v_jet_dRminToT->push_back(mindRtoT);

    m_branches->v_jet_isLLPDecayProd->push_back(LLPmatch);
    m_branches->v_jet_dRtoLLPDecayProd->push_back(closestTruthJet);
    m_branches->v_jet_truthLLPJetLabel->push_back(truth_LLP_jetLabel);
    m_branches->v_jet_truthLLP_Decay_x->push_back(truth_LLP_decay_vertex[0]);
    m_branches->v_jet_truthLLP_Decay_y->push_back(truth_LLP_decay_vertex[1]);
    m_branches->v_jet_truthLLP_Decay_z->push_back(truth_LLP_decay_vertex[2]);

}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void JetPropertiesBranches::clear() {
  // clear vectors
  m_branches->v_jet_pt->clear();
  m_branches->v_jet_eta->clear();
  m_branches->v_jet_phi->clear();
  m_branches->v_jet_E->clear();
  m_branches->v_jet_pt_orig->clear();
  m_branches->v_jet_eta_orig->clear();
  m_branches->v_jet_phi_orig->clear();
  m_branches->v_jet_E_orig->clear();
  m_branches->v_jet_LabDr_HadF->clear();
  m_branches->v_jet_DoubleHadLabel->clear();
  m_branches->v_jet_JVT->clear();
  m_branches->v_jet_m->clear();
  m_branches->v_jet_nConst->clear();
  m_branches->v_jet_dRiso->clear();
  m_branches->v_jet_truthMatch->clear();
  m_branches->v_jet_isPU->clear();
  m_branches->v_jet_aliveAfterOR->clear();
  m_branches->v_jet_aliveAfterORmu->clear();
  m_branches->v_jet_isBadMedium->clear();
  m_branches->v_jet_truthPt->clear();
  m_branches->v_jet_dRminToB->clear();
  m_branches->v_jet_dRminToC->clear();
  m_branches->v_jet_dRminToT->clear();

  m_branches->v_jet_isLLPDecayProd->clear();
  m_branches->v_jet_dRtoLLPDecayProd->clear();
  m_branches->v_jet_truthLLPJetLabel->clear();
  m_branches->v_jet_truthLLP_Decay_x->clear();
  m_branches->v_jet_truthLLP_Decay_y->clear();
  m_branches->v_jet_truthLLP_Decay_z->clear();
}

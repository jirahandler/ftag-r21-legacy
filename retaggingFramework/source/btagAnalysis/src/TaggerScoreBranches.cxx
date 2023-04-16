#include "../btagAnalysis/TaggerScoreBranches.hh"
#include "../btagAnalysis/TaggerScoreBranchBuffer.hh"

#include "xAODJet/Jet.h"
#include "AthContainers/exceptions.h"
#include "TTree.h"



//!-----------------------------------------------------------------------------------------------------------------------------!//
TaggerScoreBranches::TaggerScoreBranches():
  m_branches(new TaggerScoreBranchBuffer)
{
  // instantiate all the vectors here ...


  m_branches->v_jet_dl1_pb    = new std::vector<double>();
  m_branches->v_jet_dl1_pc    = new std::vector<double>();
  m_branches->v_jet_dl1_pu    = new std::vector<double>();
  m_branches->v_jet_dl1rmu_pb  = new std::vector<double>();
  m_branches->v_jet_dl1rmu_pc  = new std::vector<double>();
  m_branches->v_jet_dl1rmu_pu  = new std::vector<double>();
  m_branches->v_jet_dl1r_pb = new std::vector<double>();
  m_branches->v_jet_dl1r_pc = new std::vector<double>();
  m_branches->v_jet_dl1r_pu = new std::vector<double>();
  m_branches->v_jet_mv2c10    = new std::vector<double>();
  m_branches->v_jet_mv2c10mu  = new std::vector<double>();
  m_branches->v_jet_mv2c10rnn = new std::vector<double>();
  m_branches->v_jet_mv2c100   = new std::vector<double>();
  m_branches->v_jet_mv2cl100  = new std::vector<double>();
}

//!-----------------------------------------------------------------------------------------------------------------------------!//
TaggerScoreBranches::~TaggerScoreBranches() {
  // delete all the vectors here ...
  delete m_branches->v_jet_dl1_pb;
  delete m_branches->v_jet_dl1_pc;
  delete m_branches->v_jet_dl1_pu;
  delete m_branches->v_jet_dl1rmu_pb;
  delete m_branches->v_jet_dl1rmu_pc;
  delete m_branches->v_jet_dl1rmu_pu;
  delete m_branches->v_jet_dl1r_pb;
  delete m_branches->v_jet_dl1r_pc;
  delete m_branches->v_jet_dl1r_pu;
  delete m_branches->v_jet_mv2c10;
  delete m_branches->v_jet_mv2c10mu;
  delete m_branches->v_jet_mv2c10rnn;
  delete m_branches->v_jet_mv2c100;
  delete m_branches->v_jet_mv2cl100;

  delete m_branches;
}

void TaggerScoreBranches::set_tree(TTree& output_tree) const {

  output_tree.Branch( "jet_dl1_pb", &m_branches->v_jet_dl1_pb );
  output_tree.Branch( "jet_dl1_pc", &m_branches->v_jet_dl1_pc );
  output_tree.Branch( "jet_dl1_pu", &m_branches->v_jet_dl1_pu );
  output_tree.Branch( "jet_dl1rmu_pb", &m_branches->v_jet_dl1rmu_pb );
  output_tree.Branch( "jet_dl1rmu_pc", &m_branches->v_jet_dl1rmu_pc );
  output_tree.Branch( "jet_dl1rmu_pu", &m_branches->v_jet_dl1rmu_pu );
  output_tree.Branch( "jet_dl1r_pb", &m_branches->v_jet_dl1r_pb );
  output_tree.Branch( "jet_dl1r_pc", &m_branches->v_jet_dl1r_pc );
  output_tree.Branch( "jet_dl1r_pu", &m_branches->v_jet_dl1r_pu );
  output_tree.Branch( "jet_mv2c10", &m_branches->v_jet_mv2c10 );
  output_tree.Branch( "jet_mv2c10mu", &m_branches->v_jet_mv2c10mu );
  output_tree.Branch( "jet_mv2c10rnn", &m_branches->v_jet_mv2c10rnn );
  output_tree.Branch( "jet_mv2c100", &m_branches->v_jet_mv2c100 );
  output_tree.Branch( "jet_mv2cl100", &m_branches->v_jet_mv2cl100 );

}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void TaggerScoreBranches::fill(const xAOD::Jet& jet) {

    const xAOD::BTagging *bjet = jet.btagging();

    // DL1
    try {
      m_branches->v_jet_dl1_pb->push_back(bjet->auxdata<double>("DL1_pb"));
      m_branches->v_jet_dl1_pc->push_back(bjet->auxdata<double>("DL1_pc"));
      m_branches->v_jet_dl1_pu->push_back(bjet->auxdata<double>("DL1_pu"));
    } catch(...) {
      m_branches->v_jet_dl1_pb->push_back(-99);
      m_branches->v_jet_dl1_pc->push_back(-99);
      m_branches->v_jet_dl1_pu->push_back(-99);
    }

    try {
      m_branches->v_jet_dl1rmu_pb->push_back(bjet->auxdata<double>("DL1rmu_pb"));
      m_branches->v_jet_dl1rmu_pc->push_back(bjet->auxdata<double>("DL1rmu_pc"));
      m_branches->v_jet_dl1rmu_pu->push_back(bjet->auxdata<double>("DL1rmu_pu"));
    } catch(...) {
      m_branches->v_jet_dl1rmu_pb->push_back(-99);
      m_branches->v_jet_dl1rmu_pc->push_back(-99);
      m_branches->v_jet_dl1rmu_pu->push_back(-99);
    }

    try {
      m_branches->v_jet_dl1r_pb->push_back(bjet->auxdata<double>("DL1r_pb"));
      m_branches->v_jet_dl1r_pc->push_back(bjet->auxdata<double>("DL1r_pc"));
      m_branches->v_jet_dl1r_pu->push_back(bjet->auxdata<double>("DL1r_pu"));
    } catch(...) {
      m_branches->v_jet_dl1r_pb->push_back(-99);
      m_branches->v_jet_dl1r_pc->push_back(-99);
      m_branches->v_jet_dl1r_pu->push_back(-99);
    }

    try {
      m_branches->v_jet_mv2c10->push_back(bjet->auxdata<double>("MV2c10_discriminant"));
    } catch(...) {
      m_branches->v_jet_mv2c10->push_back(-99);
    }
    try {
      m_branches->v_jet_mv2c10mu->push_back(bjet->auxdata<double>("MV2c10mu_discriminant"));
    } catch(...) {
      m_branches->v_jet_mv2c10mu->push_back(-99);
    }
    try {
      m_branches->v_jet_mv2c10rnn->push_back(bjet->auxdata<double>("MV2c10rnn_discriminant"));
    } catch(...) {
      m_branches->v_jet_mv2c10rnn->push_back(-99);
    }

    try {
      m_branches->v_jet_mv2c100->push_back(bjet->auxdata<double>("MV2c100_discriminant"));
    } catch(...) {
      m_branches->v_jet_mv2c100->push_back(-99);
     }
    try {
      m_branches->v_jet_mv2cl100->push_back(bjet->auxdata<double>("MV2cl100_discriminant"));
    } catch(...) {
      m_branches->v_jet_mv2cl100->push_back(-99);
    }


}

//!-----------------------------------------------------------------------------------------------------------------------------!//
void TaggerScoreBranches::clear() {
  // clear vectors

  m_branches->v_jet_dl1_pb->clear();
  m_branches->v_jet_dl1_pc->clear();
  m_branches->v_jet_dl1_pu->clear();
  m_branches->v_jet_dl1rmu_pb->clear();
  m_branches->v_jet_dl1rmu_pc->clear();
  m_branches->v_jet_dl1rmu_pu->clear();
  m_branches->v_jet_dl1r_pb->clear();
  m_branches->v_jet_dl1r_pc->clear();
  m_branches->v_jet_dl1r_pu->clear();
  m_branches->v_jet_mv2c10->clear();
  m_branches->v_jet_mv2c10mu->clear();
  m_branches->v_jet_mv2c10rnn->clear();
  m_branches->v_jet_mv2c100->clear();
  m_branches->v_jet_mv2cl100->clear();


}






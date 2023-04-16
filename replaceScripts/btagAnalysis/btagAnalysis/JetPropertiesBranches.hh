#ifndef JETPROPERTIES_BRANCHES_HH
#define JETPROPERTIES_BRANCHES_HH

#include "TLorentzVector.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/Jet.h"
#include "xAODTruth/TruthParticle.h"

class TTree;

namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}

// branch buffers are stored as an external class to cut down on
// (re)compile time // VD I strongly disagree with it ;-)
struct JetPropertiesBranchBuffer;

class JetPropertiesBranches
{
public:
  // might want to add a prefix to the constructor for the tree branches
  JetPropertiesBranches();
  ~JetPropertiesBranches();

  // disable copying and assignment
  JetPropertiesBranches& operator=(JetPropertiesBranches) = delete;
  JetPropertiesBranches(const JetPropertiesBranches&) = delete;

  void set_tree(TTree& output_tree, bool show_debug);
  //void fill(const xAOD::Jet& orig_jet, const xAOD::Jet& calib_jet, float JVT, float dRiso, int isBad , std::vector<TLorentzVector> truth_electrons);

  void fill(const xAOD::Jet& orig_jet, const xAOD::Jet& calib_jet, float JVT, float dRiso, int isBad,
    int LLPmatch, float closestTruthJet, int truth_LLP_jetLabel, float* truth_LLP_decay_vertex,
    std::vector<TLorentzVector> truth_electrons,  std::vector<TLorentzVector> truth_muons,
    const xAOD::JetContainer *truthjets,
    std::vector<const xAOD::TruthParticle* > m_partonB,
    std::vector<const xAOD::TruthParticle* > m_partonC,
    std::vector<const xAOD::TruthParticle* > m_partonT);
  void clear();

private:

  bool debug;

  JetPropertiesBranchBuffer* m_branches;
};

#endif // JETPROPERTIES_BRANCHES_HH

#ifndef SOFTMUON_BRANCHES_HH
#define SOFTMUON_BRANCHES_HH

#include "xAODMuon/MuonContainer.h"
#include "xAODTruth/TruthEventContainer.h"


class TTree;

namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}

// branch buffers are stored as an external class to cut down on
// (re)compile time // VD I strongly disagree with it ;-)
struct SoftMuonBranchBuffer;

class SoftMuonBranches
{
public:
  // might want to add a prefix to the constructor for the tree branches
  SoftMuonBranches();
  ~SoftMuonBranches();

  // disable copying and assignment
  SoftMuonBranches& operator=(SoftMuonBranches) = delete;
  SoftMuonBranches(const SoftMuonBranches&) = delete;

  void set_tree(TTree& output_tree) const;
  void fill(const xAOD::Jet& jet);
  void clear();

private:

  SoftMuonBranchBuffer* m_branches;


  int parent_classify(const xAOD::TruthParticle *theParticle);

};

#endif // SOFTMUON_BRANCHES_HH

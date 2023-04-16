#ifndef TAGGERSCORE_BRANCHES_HH
#define TAGGERSCORE_BRANCHES_HH


class TTree;

namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}

// branch buffers are stored as an external class to cut down on
// (re)compile time // VD I strongly disagree with it ;-)
struct TaggerScoreBranchBuffer;

class TaggerScoreBranches
{
public:
  // might want to add a prefix to the constructor for the tree branches
  TaggerScoreBranches();
  ~TaggerScoreBranches();

  // disable copying and assignment
  TaggerScoreBranches& operator=(TaggerScoreBranches) = delete;
  TaggerScoreBranches(const TaggerScoreBranches&) = delete;

  void set_tree(TTree& output_tree) const;
  void fill(const xAOD::Jet& jet);
  void clear();

private:

  TaggerScoreBranchBuffer* m_branches;


};

#endif // TAGGERSCORE_BRANCHES_HH

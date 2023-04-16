#ifndef EXAMPLE_BRANCHES_HH
#define EXAMPLE_BRANCHES_HH


class TTree;

namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}

// branch buffers are stored as an external class to cut down on
// (re)compile time // VD I strongly disagree with it ;-)
struct ExampleBranchBuffer;

class ExampleBranches
{
public:
  // might want to add a prefix to the constructor for the tree branches
  ExampleBranches();
  ~ExampleBranches();

  // disable copying and assignment
  ExampleBranches& operator=(ExampleBranches) = delete;
  ExampleBranches(const ExampleBranches&) = delete;

  void set_tree(TTree& output_tree) const;
  void fill(const xAOD::Jet& jet);
  void clear();

private:

  ExampleBranchBuffer* m_branches;

  float ExampleFunction(float number);

};

#endif // EXAMPLE_BRANCHES_HH

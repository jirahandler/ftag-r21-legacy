#ifndef IMPACTPARAMETER_BRANCHES_HH
#define IMPACTPARAMETER_BRANCHES_HH

#include <map>

class TTree;

namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}

// branch buffers are stored as an external class to cut down on
// (re)compile time // VD I strongly disagree with it ;-)
struct ImpactParameterBranchBuffer;

class ImpactParameterBranches
{
public:
  // might want to add a prefix to the constructor for the tree branches
  ImpactParameterBranches();
  ~ImpactParameterBranches();

  // disable copying and assignment
  ImpactParameterBranches& operator=(ImpactParameterBranches) = delete;
  ImpactParameterBranches(const ImpactParameterBranches&) = delete;

  void set_tree(TTree& output_tree, std::map<std::string, double > defaultDict, bool replaceDefaults);
  void fill(const xAOD::Jet& jet);
  void clear();

private:

  ImpactParameterBranchBuffer* m_branches;
  std::map<std::string, double > m_defaultDict;
  bool m_replaceDefaults;

};

#endif // ImpactParameter_BRANCHES_HH

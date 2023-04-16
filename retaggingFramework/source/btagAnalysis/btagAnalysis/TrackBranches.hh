#ifndef TRACK_BRANCHES_HH
#define TRACK_BRANCHES_HH

#include "TLorentzVector.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/Jet.h"
#include "xAODTruth/TruthParticle.h"

namespace xAOD {
  class IParticle;
  class Jet_v1;
  typedef Jet_v1 Jet;
}


class TTree;

#include <vector>
#include <string>

// branch buffers are stored as an external class to cut down on
// (re)compile time
struct TrackBranchBuffer;
struct BTagTrackAccessors;

enum TAGGERALGO{ IP2D=0,
     IP3D,
     SV0,
     SV1,
     JF };

class TrackBranches
{
public:
  typedef std::vector<const xAOD::IParticle*> PartVector;

  TrackBranches();
  ~TrackBranches();

  // disable copying and assignment
  TrackBranches& operator=(TrackBranches) = delete;
  TrackBranches(const TrackBranches&) = delete;

  void set_tree(TTree& output_tree, const std::string& prefix) const;
  void fill(const PartVector& constituents, const xAOD::BTagging &btag, const xAOD::Jet& orig_jet);
  void clear();

private:
  TrackBranchBuffer* m_branches;
  BTagTrackAccessors* m_acc;


  bool particleInCollection( const xAOD::TrackParticle *trkPart, std::vector< ElementLink< xAOD::TrackParticleContainer > > trkColl );
  const xAOD::TruthParticle* truthParticle(const xAOD::TrackParticle *trkPart);
  int parent_classify(const xAOD::TruthParticle *theParticle);

  // short-circuit the branch filling if no tree is set
  mutable bool m_active;
};

#endif // TRACK_BRANCHES_HH

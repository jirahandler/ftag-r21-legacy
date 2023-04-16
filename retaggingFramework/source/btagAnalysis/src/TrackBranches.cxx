#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "AthContainers/exceptions.h"
#include "xAODBTagging/BTagging.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthEventContainer.h"


#include "../btagAnalysis/TrackBranches.hh"
#include "../btagAnalysis/TrackBranchBuffer.hh"
#include "../btagAnalysis/BTagTrackAccessors.hh"
#include "../btagAnalysis/LifetimeSigning.hh"

#include "TTree.h"

#include <string>
#include <stdexcept>
#include <memory>

TrackBranches::TrackBranches():
  m_branches(new TrackBranchBuffer),
  m_acc(new BTagTrackAccessors),
  m_active(false)
{
  typedef std::vector<std::vector<int> > vvi_t;
  typedef std::vector<std::vector<float> > vvf_t;

  m_branches->ntrk = new std::vector<int>;

  m_branches->pt = new vvf_t;
  m_branches->eta = new vvf_t;
  m_branches->theta = new vvf_t;
  m_branches->phi = new vvf_t;
  m_branches->qoverp = new vvf_t;
  m_branches->charge = new vvf_t;


}

TrackBranches::~TrackBranches()
{

  delete m_branches->ntrk;

  delete m_branches->pt;
  delete m_branches->eta;
  delete m_branches->theta;
  delete m_branches->phi;
  delete m_branches->qoverp;
  delete m_branches->charge;
  delete m_acc;
}

namespace {
  // branch name function (lowercase the variable name)
  std::string brnm(const std::string& pfx, std::string in) {
    // std::transform(in.begin(), in.end(), in.begin(), ::tolower);
    return pfx + in;
  }
}

void TrackBranches::set_tree(TTree& output_tree,
                             const std::string& prefix) const {
  if (m_active) throw std::logic_error("tried to set tree twice");
  m_active = true;

#define ADD_SIMPLE(nm) \
  output_tree.Branch(brnm(prefix, #nm).c_str(), &m_branches->nm)
  // basic kinematics
  ADD_SIMPLE(pt);
  ADD_SIMPLE(eta);
  ADD_SIMPLE(theta);
  ADD_SIMPLE(phi);
  ADD_SIMPLE(qoverp);
  ADD_SIMPLE(charge);

#undef ADD_SIMPLE


  output_tree.Branch("jet_trk_ntrk", &m_branches->ntrk);
}

namespace {
  // util function to get track values
  int get(const xAOD::TrackParticle& part, xAOD::SummaryType info) {
    uint8_t val;
    bool ok = part.summaryValue(val, info);
    if (!ok) throw std::logic_error("problem getting track summary value");
    return val;
  }
}

void TrackBranches::fill(const TrackBranches::PartVector& tracks, const xAOD::BTagging& btag, const xAOD::Jet& orig_jet) {
  if (!m_active) return;

  TLorentzVector jetDir;
  jetDir.SetPtEtaPhiE(orig_jet.pt(), orig_jet.eta(), orig_jet.phi(), orig_jet.e());
  Amg::Vector3D jetDirection(jetDir.Px(),jetDir.Py(),jetDir.Pz());
  Amg::Vector3D unit = jetDirection.unit();




  const int MISSING_VALUE = -999;


  typedef ElementLink<xAOD::TrackParticleContainer> TrackLink;
  typedef std::vector<TrackLink> TrackLinks;


  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> theta;
  std::vector<float> phi;
  std::vector<float> qoverp;
  std::vector<float> charge;

  for (const auto* part: tracks) {
    const auto* tmpTrk = dynamic_cast<const xAOD::TrackParticle*>(part);
    if (!tmpTrk) throw std::logic_error("This isn't a track particle");
    const auto& track = *tmpTrk;
    //pt.push_back(tmpTrk->pt());
    //eta.push_back(tmpTrk->eta());
    //theta.push_back(tmpTrk->theta());
    //phi.push_back(tmpTrk->phi());
    //qoverp.push_back(tmpTrk->qOverP());
    //charge.push_back(tmpTrk->charge());


#define GET(NAME) NAME.push_back(get(track, xAOD::NAME))

#undef GET

    const xAOD::TruthParticle* truth = truthParticle(tmpTrk);


  } // end of track loop

  // push back into member vectors
#define PUSH(var) m_branches->var->push_back(std::move(var))
  //PUSH(pt);
  //PUSH(eta);
  //PUSH(theta);
  //PUSH(phi);
  //PUSH(qoverp);
  //PUSH(charge);

#undef PUSH

  m_branches->ntrk->push_back(tracks.size());
}

void TrackBranches::clear() {
#define CLEAR(name) m_branches->name->clear()
  CLEAR(pt);
  CLEAR(eta);
  CLEAR(theta);
  CLEAR(phi);
  CLEAR(qoverp);
  CLEAR(charge);

#undef CLEAR

  m_branches->ntrk->clear();
}


bool TrackBranches::particleInCollection( const xAOD::TrackParticle *trkPart, std::vector< ElementLink< xAOD::TrackParticleContainer > > trkColl ) {
  for (unsigned int iT = 0; iT < trkColl.size(); iT++) {
    if (trkPart == *(trkColl.at(iT))) return true;
  }
  return false;
}



const xAOD::TruthParticle* TrackBranches::truthParticle(const xAOD::TrackParticle *trkPart){
  typedef ElementLink< xAOD::TruthParticleContainer > Link_t;
  static const char *NAME = "truthParticleLink";
  if( ! trkPart->isAvailable< Link_t >( NAME ) ) {
    return 0;
  }
  const Link_t& link = trkPart->auxdata< Link_t >( NAME );
  if (!link.isValid()) {
    return 0;
  }
  return *link;
}

int TrackBranches :: parent_classify(const xAOD::TruthParticle *theParticle) {
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


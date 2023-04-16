#include "TrackSystematicsAlg.h"
#include "InDetTrackSystematicsTools/IInDetTrackTruthFilterTool.h"
#include "InDetTrackSystematicsTools/IJetTrackFilterTool.h"
#include "InDetTrackSystematicsTools/IInDetTrackBiasingTool.h"
#include "InDetTrackSystematicsTools/IInDetTrackSmearingTool.h"
#include "InDetTrackSystematicsTools/InDetTrackSystematics.h"

#include "CxxUtils/make_unique.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/ShallowAuxInfo.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

using CxxUtils::make_unique;

namespace InDet {

  TrackSystematicsAlg::TrackSystematicsAlg(const std::string& name, ISvcLocator *pSvcLocator)
    : AthAlgorithm(name, pSvcLocator)
    , m_trackCollectionName("InDetTrackParticles")
    , m_outputTrackCollectionName("InDetTrackParticlesModified")
    , m_trackFilterTool("InDet::InDetTrackTruthFilterTool")
    , m_trackSmearingTool("InDet::InDetTrackSmearingTool")
    , m_trackBiasingTool("InDet::InDetTrackBiasingTool")
    , m_jetTrackFilterTool("InDet::JetTrackFilterTool")
  {
    declareProperty( "trackCollectionName", m_trackCollectionName );
    declareProperty( "outputTrackCollectionName", m_outputTrackCollectionName );
    declareProperty( "systematicVariations", m_systNames );
    declareProperty( "trackFilterTool", m_trackFilterTool );
    declareProperty( "trackSmearingTool", m_trackSmearingTool );
    declareProperty( "trackBiasingTool", m_trackBiasingTool );
    declareProperty( "jetTrackFilterTool", m_jetTrackFilterTool );
    declareProperty( "jetCollectionName", m_jetCollectionName );
  }

  TrackSystematicsAlg::~TrackSystematicsAlg() = default;

  StatusCode TrackSystematicsAlg::initialize() {

    ATH_CHECK( m_trackFilterTool.retrieve() );
    ATH_CHECK( m_jetTrackFilterTool.retrieve() );
    ATH_CHECK( m_trackSmearingTool.retrieve() );
    ATH_CHECK( m_trackBiasingTool.retrieve() );

    CP::SystematicSet activeSysts;
    for (const auto& name : m_systNames) {
      using pair_t = std::pair< InDet::TrackSystematic, CP::SystematicVariation >;
      auto it_systName = std::find_if(InDet::TrackSystematicMap.cbegin(),
				      InDet::TrackSystematicMap.cend(),
				      [&](const pair_t& pair){return pair.second.name() == name;});
      if (it_systName != InDet::TrackSystematicMap.cend()) {
	ATH_MSG_INFO( "Adding " << it_systName->second.name() );
	activeSysts.insert(it_systName->second);
      } else {
	ATH_MSG_ERROR( "Unrecognized systematic variation: " << name );
	return StatusCode::FAILURE;
      }
    }

    ATH_MSG_DEBUG( "activeSysts in " << name() << ":" );
    for (const auto& syst : activeSysts) {
      ATH_MSG_DEBUG( "  " << syst.name() );
    }

    auto codeSmear = m_trackSmearingTool->applySystematicVariation(activeSysts);
    if (codeSmear != CP::SystematicCode::Ok) {
      ATH_MSG_ERROR( "Failed to apply systematic variation to " << m_trackSmearingTool );
      return StatusCode::FAILURE;
    }
    auto codeFilter = m_trackFilterTool->applySystematicVariation(activeSysts);
    if (codeFilter != CP::SystematicCode::Ok) {
      ATH_MSG_ERROR( "Failed to apply systematic variation to " << m_trackFilterTool );
      return StatusCode::FAILURE;
    }
    const CP::SystematicVariation& tideSyst = InDet::TrackSystematicMap[InDet::TRK_EFF_LOOSE_TIDE];
    m_doJets = activeSysts.find(tideSyst) != activeSysts.end();
    if (m_doJets) ATH_MSG_INFO( "Jet collection set to " << m_jetCollectionName );
    auto codeJetFilter = m_jetTrackFilterTool->applySystematicVariation(activeSysts);
    if (codeJetFilter != CP::SystematicCode::Ok) {
      ATH_MSG_ERROR( "Failed to apply systematic variation to " << m_jetTrackFilterTool );
      return StatusCode::FAILURE;
    }
    auto codeBias = m_trackBiasingTool->applySystematicVariation(activeSysts);
    if (codeBias != CP::SystematicCode::Ok) {
      ATH_MSG_ERROR( "Failed to apply systematic variation to " << m_trackBiasingTool );
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  StatusCode TrackSystematicsAlg::execute() {

    const xAOD::TrackParticleContainer *tracks = nullptr;
    ATH_CHECK( evtStore()->retrieve( tracks, m_trackCollectionName) );

    const xAOD::JetContainer* jets = nullptr;
    if (m_doJets) {
      ATH_CHECK( evtStore()->retrieve( jets, m_jetCollectionName ) );
    }
    
    auto newTracks = make_unique<xAOD::TrackParticleContainer>();
    auto newTracksAux = make_unique<xAOD::AuxContainerBase>();
    newTracks->setStore(newTracksAux.get());

    using Link_t = ElementLink< xAOD::TrackParticleContainer >;
    using Deco_t = SG::AuxElement::Decorator< Link_t >;
    static Deco_t setOriginLink("originalTrackLink");
    static Deco_t setModifiedLink("correctedTrackLink");

    for ( const xAOD::TrackParticle* track : *tracks ) {

      if( !m_trackFilterTool->accept(track) ) continue;
      if( m_doJets && !m_jetTrackFilterTool->accept(track, jets) ) continue;

      xAOD::TrackParticle *newTrack = new xAOD::TrackParticle();
      newTracks->push_back(newTrack);
      *newTrack = *track; // copy over the information

      // apply bias correction to the track
      if (m_trackBiasingTool->applyCorrection(*newTrack) == CP::CorrectionCode::Error) {
	ATH_MSG_ERROR( "Could not apply InDetTrackBiasingTool." );
	return StatusCode::FAILURE;
      }

      // apply a smearing to the track
      if (m_trackSmearingTool->applyCorrection(*newTrack) == CP::CorrectionCode::Error) {
	ATH_MSG_ERROR( "Could not apply InDetTrackSmearingTool." );
	return StatusCode::FAILURE;
      }

      Link_t originLink(*tracks, track->index());
      setOriginLink(*newTrack) = originLink;
      
      Link_t modifiedLink(*newTracks, newTrack->index());
      setModifiedLink(*track) = modifiedLink;

    }

    ATH_MSG_VERBOSE( tracks->size() << " / " << newTracks->size() << " \toriginal / modified tracks" );

    ATH_CHECK( evtStore()->record(std::move(newTracks), m_outputTrackCollectionName) );
    ATH_CHECK( evtStore()->record(std::move(newTracksAux), m_outputTrackCollectionName+"Aux.") ); 

    return StatusCode::SUCCESS;
  }

  StatusCode TrackSystematicsAlg::finalize(){
    CHECK( m_trackFilterTool->release() );
    //    CHECK( m_jetTrackFilterTool->release() );
    CHECK( m_trackBiasingTool->release() );
    CHECK( m_trackSmearingTool->release() );
    return StatusCode::SUCCESS;
  }

}

#pragma once

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

namespace InDet {

  class IInDetTrackTruthFilterTool;
  class IInDetTrackSmearingTool;
  class IInDetTrackBiasingTool;
  class IJetTrackFilterTool;

/**
   @class TrackSystematicsAlg
   Algorithms that reads track container from StoreGate,
   create a copy of the tracks, apply tracking 
   variations and saves back the altered tracks.
   @author Remi Zaidan (remi.zaidan@cern.ch)
   @author Felix Clark (michael.ryan.clark@cern.ch)
*/
class TrackSystematicsAlg :  public AthAlgorithm { 
public:
  
  /** Constructors and destructors */
  TrackSystematicsAlg(const std::string& name, ISvcLocator *pSvcLocator);
  virtual ~TrackSystematicsAlg();
  
  /** Main routines specific to an ATHENA algorithm */
  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();
  
private:
  
  std::string m_trackCollectionName;
  std::string m_outputTrackCollectionName;

  // currently can only apply a single systematic variation at a time.
  // all the systematics in the set will be used at once.
  // should implement multiple variations
  std::vector< std::string > m_systNames;

  // the name of the jet collection to evaluate the TIDE systematic with
  std::string m_jetCollectionName = "AntiKt4EMTopoJets";
  bool m_doJets = false; // whether to retrieve a jet container. this is set automatically in initialize()

  ToolHandle< IInDetTrackTruthFilterTool > m_trackFilterTool;
  ToolHandle< IInDetTrackSmearingTool > m_trackSmearingTool; 
  ToolHandle< IInDetTrackBiasingTool > m_trackBiasingTool;
  ToolHandle< IJetTrackFilterTool > m_jetTrackFilterTool;

}; // class TrackSystematicsAlg

} // namespace InDet


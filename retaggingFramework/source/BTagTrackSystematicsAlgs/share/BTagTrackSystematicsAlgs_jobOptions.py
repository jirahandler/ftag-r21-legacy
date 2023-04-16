####
## example job options for using the track systematics algorithm

# the systematic variation we will apply.
# see https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TrackingID/InDetTrackSystematicsTools/trunk/InDetTrackSystematicsTools/InDetTrackSystematics.h for a full list
systName = "TRK_BIAS_D0_WM"
# systName = "TRK_EFF_LOOSE_TIDE"

tracksInputKey = "InDetTrackParticles"
# name the output track collection with the systematic variation you want to apply
tracksOutputKey = tracksInputKey + "_" + systName
slimTracks = True # set this to True to cut down on the size of the output containers significantly


import AthenaPoolCnvSvc.ReadAthenaPool
algSeq = CfgMgr.AthSequencer("AthAlgSeq")
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags

asgTestDir = "/afs/cern.ch/user/a/asgbase/patspace/xAODs/r6630/"
filename = asgTestDir+"mc15_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s2576_s2132_r6630_tid05358802_00/AOD.05358802._004308.pool.root.1"
# filename = "/afs/cern.ch/work/r/rjansky/public/4Felix/group.perf-idtracking.426045.HerwigppEvtGen_UEEE5_CTEQ6L1_jetjet_JZ5W.recon.DAOD_IDTIDE.e4410_s2608_r6869.v2_EXT0/group.perf-idtracking.7307936.EXT0._000208.DOAD_TIDE.pool.root"
athenaCommonFlags.FilesInput = [filename]
svcMgr.EventSelector.InputCollections = athenaCommonFlags.FilesInput()

# in order to slim tracks, we need to use a generic AuxContainer instead of a TrackParticleContainer
# the following algorithm will do the conversion for us, but must be scheduled at the very beginning.
if (slimTracks):
    algSeq += CfgMgr.xAODMaker__AuxStoreWrapper( "MyAuxStoreWrapperAlg",
                                                 OutputLevel=INFO )
    algSeq.MyAuxStoreWrapperAlg.SGKeys = [ tracksInputKey+"Aux." ]

from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackTruthOriginTool
trackOriginTool = InDet__InDetTrackTruthOriginTool()
ToolSvc += trackOriginTool

from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackTruthFilterTool
trackFilterTool = InDet__InDetTrackTruthFilterTool()
trackFilterTool.Seed = 1234
trackFilterTool.trackOriginTool = trackOriginTool
ToolSvc += trackFilterTool

from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__JetTrackFilterTool
jetTrackFilterTool = InDet__JetTrackFilterTool(OutputLevel=DEBUG)
jetTrackFilterTool.Seed = 4321
ToolSvc += jetTrackFilterTool

from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackSmearingTool
trackSmearingTool = InDet__InDetTrackSmearingTool()
trackSmearingTool.Seed = 12345
ToolSvc += trackSmearingTool

from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackBiasingTool
trackBiasingTool = InDet__InDetTrackBiasingTool()
ToolSvc += trackBiasingTool


from BTagTrackSystematicsAlgs.BTagTrackSystematicsAlgsConf import BTag__TrackSystematicsAlg
alg = BTag__TrackSystematicsAlg("TrackSystematicsAlg_"+systName)
#alg.OutputLevel = DEBUG
alg.trackFilterTool = trackFilterTool
alg.jetTrackFilterTool = jetTrackFilterTool
alg.trackSmearingTool = trackSmearingTool
alg.trackBiasingTool = trackBiasingTool

# all of these variations will be applied at once.
# The multiple systmatics here are for testing, and probably aren't what you would apply in an analysis
# for an analysis you should only use a single systematic at a time.
alg.systematicVariations = [systName]
                            # "TRK_RES_D0_MEAS",
                            # "TRK_BIAS_Z0_WM",
                            # "TRK_EFF_LOOSE_GLOBAL", "TRK_EFF_LOOSE_IBL", "TRK_EFF_LOOSE_PP0"
# unimplementing the TIDE systematic until the final data file is ready
#                            "TRK_EFF_LOOSE_TIDE"

alg.trackCollectionName = tracksInputKey
alg.outputTrackCollectionName = tracksOutputKey
alg.jetCollectionName = "AntiKt4EMTopoJets"

algSeq += alg

# create an output stream

streamName = "xAODSkimStream"

from OutputStreamAthenaPool.MultipleStreamManager import MSMgr
xaodStream = MSMgr.NewPoolRootStream( streamName, "TrackingSyst_" + systName + ".pool.root" )

# we definitely want to keep the event info
xaodStream.AddItem( "xAOD::EventInfo#EventInfo" )
xaodStream.AddItem( "xAOD::EventAuxInfo#EventInfoAux." )

# the vertices are usually relevant for tracking analyses as well
xaodStream.AddItem( "xAOD::VertexContainer#PrimaryVertices" )
xaodStream.AddItem( "xAOD::VertexAuxContainer#PrimaryVerticesAux." )

xaodStream.AddItem( "xAOD::TrackParticleContainer#" + tracksInputKey ) # record original track particle container
xaodStream.AddItem( "xAOD::TrackParticleContainer#" + tracksOutputKey ) #record modified container

auxString = "Aux."
if (slimTracks):
    auxString += "d0.z0.phi.theta.qOverP" # this can be extended with as many variables as you wish

xaodStream.AddItem( "xAOD::AuxContainerBase#"+tracksInputKey+auxString )
xaodStream.AddItem( "xAOD::AuxContainerBase#"+tracksOutputKey+auxString )


# set up basic reader tool
import AthenaRootComps.ReadAthenaxAODHybrid

# utilities used below
from AthenaCommon.AlgSequence import AlgSequence

from btagAnalysis.configHelpers import get_short_name
from btagAnalysis.configHelpers import setupTools
from btagAnalysis.configHelpers import get_calibration_tool

# set up inputs
svcMgr.EventSelector.InputCollections = [
    '/afs/cern.ch/work/d/dguest/data/jue/DAOD_FTAG4.hhbbbb.pool.root',
   ]
svcMgr += CfgMgr.AthenaEventLoopMgr(EventPrintoutInterval=100)

# define some common tools that are used by this algorithm
setupTools(ToolSvc, CfgMgr)

# start the alg sequence here
algSeq = AlgSequence()
algSeq += CfgMgr.BTagVertexAugmenter()

# set up output stream
JetCollections = [
 #   'AntiKt4EMTopoJets',
# 'AntiKt10LCTopoTrimmedPtFrac5SmallR20ExKt2SubJets'
    'AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets'
    ]
associated_subjets = {
    'jet_vrtrkjet_': 'AntiKt10LCTopoTrimmedPtFrac5SmallR20ExKt2SubJets',
}
svcMgr += CfgMgr.THistSvc()
for jet in JetCollections:
    shortJetName=get_short_name(jet)
    svcMgr.THistSvc.Output += [
        shortJetName+" DATAFILE='flav_"+shortJetName+".root' OPT='RECREATE'"]


### Main Ntuple Dumper Algorithm
for JetCollection in JetCollections:
    shortJetName=get_short_name(JetCollection)
    alg = CfgMgr.btagAnalysisAlg(
        "BTagDumpAlg_"+JetCollection,
        OutputLevel= DEBUG,
        Stream=shortJetName,
#    JVTtool=ToolSvc.JVT,
    )

    alg.JetCollectionName = JetCollection
    if "TrackJets" in JetCollection or "Truth" in JetCollection:
        alg.CleanJets     = False
        alg.CalibrateJets = False
    alg.CleanJets     = False
    alg.CalibrateJets = False

    alg.JetCleaningTool.CutLevel= "LooseBad"
    alg.JetCleaningTool.DoUgly  = True
    alg.subjetCollections = associated_subjets

    ## what to include in ntuple ####
    alg.EventInfo = True
    alg.JetProperties = False
    #taggers (MV2, DL1)
    alg.TaggerScores = False
    alg.CleanParentJet = False
    ##IPxD+RNNIP
    alg.ImpactParameterInfo = False #True
    ##SV1
    alg.SVInfo = False# True
    alg.branchDebug = True
    alg.AccessBtagObject = True
    alg.svxCollections = {'jet_sv1_': 'SV1'}
    ##JetFitter
    alg.JetFitterInfo = False
    ###SoftMuonTagger
    alg.SoftMuoninfo = False#True
    alg.doJVT = False
    ## b and c hadron truth info
    alg.bHadInfo = False#True
    #include all b and c decay products, and trk_origin
    alg.bHadExtraInfo = False
    ## kshort
    alg.kshortInfo = False
    #example
    #alg.exampleBranchInfo = False

    algSeq += alg

    ToolSvc += get_calibration_tool(CfgMgr, JetCollection, False)


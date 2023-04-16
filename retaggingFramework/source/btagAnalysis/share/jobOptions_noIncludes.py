# set up basic reader tool
import AthenaRootComps.ReadAthenaxAODHybrid

# utilities used below
from AthenaCommon.AlgSequence import AlgSequence

from btagAnalysis.configHelpers import get_short_name
from btagAnalysis.configHelpers import setupTools
from btagAnalysis.configHelpers import get_calibration_tool

# set up inputs
svcMgr.EventSelector.InputCollections = [
    '/afs/cern.ch/work/d/dguest/data/jue/DAOD_FTAG4.ttbar.pool.root',
    # '/afs/cern.ch/work/d/dguest/public/rachael/DAOD_FTAG2.ttbar_mc16.pool.root'
]
svcMgr += CfgMgr.AthenaEventLoopMgr(EventPrintoutInterval=100)

# define some common tools that are used by this algorithm
setupTools(ToolSvc, CfgMgr)

# start the alg sequence here
algSeq = AlgSequence()
algSeq += CfgMgr.BTagVertexAugmenter()

# set up output stream
JetCollections = [
    'AntiKt4EMTopoJets',
    ]
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
        OutputLevel=INFO,
        Stream=shortJetName,
        JVTtool=ToolSvc.JVT,
    )

    alg.JetCollectionName = JetCollection
    if "TrackJets" in JetCollection or "Truth" in JetCollection:
        alg.CleanJets     = False
        alg.CalibrateJets = False

    alg.JetCleaningTool.CutLevel= "LooseBad"
    alg.JetCleaningTool.DoUgly  = True

    ## what to include in ntuple ####
    alg.EventInfo = True
    alg.JetProperties = True
    #taggers (MV2, DL1)
    alg.TaggerScores = True
    ##IPxD+RNNIP
    alg.ImpactParameterInfo = True
    ##SV1
    alg.SVInfo = True
    alg.svxCollections = {'jet_sv1_': 'SV1'}
    ##JetFitter
    alg.JetFitterInfo = True
    ###SoftMuonTagger
    alg.SoftMuoninfo = True
    ## b and c hadron truth info
    alg.bHadInfo = True
    #include all b and c decay products, and trk_origin
    alg.bHadExtraInfo = False
    ## kshort
    alg.kshortInfo = False
    #example
    #alg.exampleBranchInfo = False

    algSeq += alg

    ToolSvc += get_calibration_tool(CfgMgr, JetCollection, False)


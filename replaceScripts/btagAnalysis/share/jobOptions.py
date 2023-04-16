##########################################################################################################################################################
##########################################################################################################################################################

doRetag = False # no re-tagging

############### For Nominal SF
systName = ''
#systName = "TRK_RES_D0_MEAS"
#systName = "TRK_RES_Z0_MEAS"
#systName = "TRK_FAKE_RATE_LOOSE"


############### For sytematics of SF
#systName = "TRK_RES_D0_MEAS_UP"
#systName = "TRK_RES_D0_MEAS_DOWN"
#systName = "TRK_RES_Z0_MEAS_UP"
#systName = "TRK_RES_Z0_MEAS_DOWN"
#systName = "TRK_RES_D0_MEAS_DiJetMap"
#systName = "TRK_RES_Z0_MEAS_DiJetMap"
#systName = "TRK_RES_D0Z0_MEAS"
#systName = "TRK_RES_D0Z0Corl_MEAS"
#systName = 'TRK_FAKE_RATE_TIGHT'

############### Didn't have the following yet
#systName = 'TRK_TAIL_RATE'
#systName = 'TRK_STRANGEHAD_RATE'
#systName = 'TRK_CONV_HADINT_RATE'
#systName = 'TRK_TAIL_PAR_SYS'
#systName = 'TRK_TAIL_RATE_SYS'

from MV2defaults import default_values

#JetCollections = [ 'AntiKt4EMTopoJets' ]
JetCollections = [ 'AntiKt4EMPFlowJets' ]

#########################################################################################################################################################
#########################################################################################################################################################
### Define input xAOD and output ntuple file name
import glob
from AthenaCommon.AthenaCommonFlags import jobproperties as jp
jp.AthenaCommonFlags.EvtMax.set_Value_and_Lock( vars().get('EVTMAX', -1) )
#jp.AthenaCommonFlags.SkipEvents.set_Value_and_Lock(16400)
#jp.AthenaCommonFlags.EvtMax.set_Value_and_Lock(10)

runData = False

jp.AthenaCommonFlags.FilesInput = [

##For liverpool
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000001.pool.root.1'
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000002.pool.root.1',
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000003.pool.root.1',
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000004.pool.root.1',
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000005.pool.root.1'
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000006.pool.root.1',
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000007.pool.root.1',
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000008.pool.root.1',
#'/hepstore/awychan/PhD/QT/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000009.pool.root.1'

##for lxplus
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000001.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000002.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000003.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000004.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000005.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000006.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000007.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000008.pool.root.1',
#'/afs/cern.ch/work/w/wachan/public/FTAG1-ttbar-DAOD/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_s3126_r9364_p3703/DAOD_FTAG1.16200019._000009.pool.root.1'

#'/afs/cern.ch/work/k/khanov/public/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.deriv.DAOD_FTAG1.e6337_a875_r10201_p4133/DAOD_FTAG1.21024631._005916.pool.root.1'

#'/afs/cern.ch/work/k/khanov/public/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.merge.AOD.e6337_e5984_s3126_r11590_r10726/AOD.19032588._000097.pool.root.1'

'/tmp/khanov/mc16_13TeV.503822.MGPy8EG_A14NNPDF23LO_GG_rpvHF_1100_200_0p01ns.deriv.DAOD_TOPQ1.e8342_s3126_r10724_p3832/DAOD_TOPQ1.26850254._000001.pool.root.1',
'/tmp/khanov/mc16_13TeV.503822.MGPy8EG_A14NNPDF23LO_GG_rpvHF_1100_200_0p01ns.deriv.DAOD_TOPQ1.e8342_s3126_r10724_p3832/DAOD_TOPQ1.26850254._000002.pool.root.1',
'/tmp/khanov/mc16_13TeV.503822.MGPy8EG_A14NNPDF23LO_GG_rpvHF_1100_200_0p01ns.deriv.DAOD_TOPQ1.e8342_s3126_r10724_p3832/DAOD_TOPQ1.26850254._000003.pool.root.1',
'/tmp/khanov/mc16_13TeV.503822.MGPy8EG_A14NNPDF23LO_GG_rpvHF_1100_200_0p01ns.deriv.DAOD_TOPQ1.e8342_s3126_r10724_p3832/DAOD_TOPQ1.26850254._000004.pool.root.1'

]

# from PyUtils import AthFile
# af = AthFile.fopen( jp.AthenaCommonFlags.FilesInput()[0] )

evtPrintoutInterval = vars().get('EVTPRINT', 5000)
svcMgr += CfgMgr.AthenaEventLoopMgr( EventPrintoutInterval=evtPrintoutInterval )

svcMgr += CfgMgr.THistSvc()

for jet in JetCollections:

  shortJetName=jet.replace("AntiKt","Akt").replace("TopoJets","To").replace("TrackJets","Tr").replace("PFlowJets","Pf")
  svcMgr.THistSvc.Output += [ shortJetName+" DATAFILE='flav_"+shortJetName+systName+".root' OPT='RECREATE'"]

##########################################################################################################################################################
##########################################################################################################################################################

from RecExConfig.RecFlags import rec
rec.doESD.set_Value_and_Lock        (False)
rec.doWriteESD.set_Value_and_Lock   (False)
rec.doAOD.set_Value_and_Lock        (False)
rec.doWriteAOD.set_Value_and_Lock   (False)
rec.doWriteTAG.set_Value_and_Lock   (False)
rec.doDPD.set_Value_and_Lock        (False)
rec.doTruth.set_Value_and_Lock      (False)

rec.doApplyAODFix.set_Value_and_Lock(False)
include ("RecExCommon/RecExCommon_topOptions.py")

from AthenaCommon.AlgSequence import AlgSequence
algSeq = AlgSequence()


from AnaAlgorithm.DualUseConfig import createAlgorithm
algSeq += createAlgorithm( 'CP::SysListLoaderAlg', 'SysLoaderAlg' )
algSeq.SysLoaderAlg.sigmaRecommended = 1

from AsgAnalysisAlgorithms.OverlapAnalysisSequence import makeOverlapAnalysisSequence

overlapSequence = makeOverlapAnalysisSequence("mc",
                                                doMuPFJetOR = True, doMuons = False, doElectrons = False, doJets = True, doTaus = False, doPhotons = False, doFatJets = False)
overlapSequence.configure(
    inputName = {'muons': 'Muons', 'jets': 'AntiKt4EMPFlowJets'},
    outputName = {'muons': 'Muons_OR', 'jets': 'AntiKt4EMPFlowJets_OR'},
    affectingSystematics = {'jets': ''})

algSeq += overlapSequence


#------------------------------------------------------------------------------------------------------------------------
# Apply systematics to build new track collection
#------------------------------------------------------------------------------------------------------------------------
tracksInputKey = "InDetTrackParticles"
tracksOutputKey = tracksInputKey

if systName is not "":

  doRetag = True

  # name the output track collection with the systematic variation you want to apply
  tracksOutputKey = tracksInputKey + "_" + systName

  from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackTruthOriginTool
  trackOriginTool = InDet__InDetTrackTruthOriginTool()
  ToolSvc += trackOriginTool

  from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackTruthFilterTool
  trackFilterTool = InDet__InDetTrackTruthFilterTool()
  trackFilterTool.Seed = 1234
  trackFilterTool.trackOriginTool = trackOriginTool
  ToolSvc += trackFilterTool

  from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__JetTrackFilterTool
  jetTrackFilterTool = InDet__JetTrackFilterTool(OutputLevel=INFO)
  jetTrackFilterTool.Seed = 4321
  ToolSvc += jetTrackFilterTool

  from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackSmearingTool
  trackSmearingTool = InDet__InDetTrackSmearingTool()
  trackSmearingTool.Seed = 12345
#  trackSmearingTool.calibFileIP_lowpt = "trackIPAlign_lowpt.root"
#  trackSmearingTool.calibFileIP_highpt = "trackIPAlign_highpt.root"
  #trackSmearingTool.systSF = float(systSF) #+++++++THIS CAN BE UNCOMMENTED IF CUSTOMIZED TRACKING SYS IS USED

  ToolSvc += trackSmearingTool

  from InDetTrackSystematicsTools.InDetTrackSystematicsToolsConf import InDet__InDetTrackBiasingTool
  trackBiasingTool = InDet__InDetTrackBiasingTool()
  ToolSvc += trackBiasingTool

  from BTagTrackSystematicsAlgs.BTagTrackSystematicsAlgsConf import InDet__TrackSystematicsAlg
  alg = InDet__TrackSystematicsAlg("TrackSystematicsAlg_"+systName)
  alg.OutputLevel = INFO
  alg.trackFilterTool = trackFilterTool
  alg.jetTrackFilterTool = jetTrackFilterTool
  alg.trackSmearingTool = trackSmearingTool
  alg.trackBiasingTool = trackBiasingTool

  # all of these variations will be applied at once.
  # The multiple systmatics here are for testing, and probably aren't what you would apply in an analysis
  # for an analysis you should only use a single systematic at a time.
  alg.systematicVariations = [systName]
  # "TRK_RES_D0_MEAS", "TRK_BIAS_Z0_WM", "TRK_EFF_LOOSE_GLOBAL", "TRK_EFF_LOOSE_IBL", "TRK_EFF_LOOSE_PP0"
  # unimplementing the TIDE systematic until the final data file is ready: "TRK_EFF_LOOSE_TIDE"

  alg.trackCollectionName = tracksInputKey
  alg.outputTrackCollectionName = tracksOutputKey
  #alg.jetCollectionName = "AntiKt4EMTopoJets"
  alg.jetCollectionName = "AntiKt4EMPFlowJets"

  algSeq += alg


if doRetag:

  #------------------------------
  # ASC
  #------------------------------
  from BTagging.BTaggingFlags import BTaggingFlags
  print 'READING RetagFragment'
  include("RetagFragment.py")
  print 'FINISH READING RetagFragment'

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

##########################################################################################################################################################
### Tools

jvt = CfgMgr.JetVertexTaggerTool('JVT')
ToolSvc += jvt

ToolSvc += CfgMgr.CP__PileupReweightingTool("prw",
                                            OutputLevel = INFO,
                                            UsePeriodConfig= "MC16"
                                            )


from TrkVertexFitterUtils.TrkVertexFitterUtilsConf import Trk__TrackToVertexIPEstimator
ToolSvc+=Trk__TrackToVertexIPEstimator("trkIPEstimator")


##########################################################################################################################################################

algSeq += CfgMgr.BTagVertexAugmenter()

### Main Ntuple Dumper Algorithm
for JetCollection in JetCollections:


  shortJetName=JetCollection.replace("AntiKt","Akt").replace("TopoJets","To").replace("TrackJets","Tr").replace("PFlowJets","Pf")
  alg = CfgMgr.btagAnalysisAlg("BTagDumpAlg_"+JetCollection,
                                  OutputLevel=INFO, #DEBUG
                                  Stream=shortJetName,
                                  JVTtool=ToolSvc.JVT,
                                  )


  alg.JetCollectionName = JetCollection
  alg.doJVT = True #if this is false JVT is NAN, if true an attempt is made to update JVT after calibration


  alg.DefaultValDictionary = default_values
  alg.ReplaceNanDefaults = True

  ## A.X.: cleaning already done

  #if "TrackJets" in JetCollection or "Truth" in JetCollection:

  #  alg.CleanJets     = False
  #  alg.CalibrateJets = False
  #  alg.doJVT = False

  alg.CleanJets     = False
  alg.CalibrateJets = False
  alg.doJVT = False

  ##

  alg.JetCleaningTool.CutLevel= "LooseBad"
  alg.JetCleaningTool.DoUgly  = True



## what to include in ntuple ####
  #example
  #alg.exampleBranchInfo = False

  if runData:
    alg.isData = True
  alg.EventInfo = True
  alg.retriveTruthJets = True
  alg.JetProperties = True
  ##alg.JetPt = 15000
  #taggers (MV2, DL1)
  alg.TaggerScores = True
  ##IPxD+RNNIP
  alg.ImpactParameterInfo = False # was True 
  ##SV1
  alg.SVInfo = False # was True
  alg.svxCollections = {'jet_sv1_': 'SV1'}
  ##JetFitter
  alg.JetFitterInfo = False
  ###SoftMuonTagger
  alg.SoftMuoninfo = False
  ## b and c hadron truth info
  alg.bHadInfo = True
  alg.bHadExtraInfo = False #include all b and c decay products, and trk_origin
  ## kshort
  alg.kshortInfo = False

  #show debug info for branches
  alg.branchDebug = False

  #track information
  alg.TrackInfo = False
  alg.nRequiredSiHits = 2 #number of hits required to save a track

  #you can disable the track augmenter if youre not filling the track branches
  """
  algSeq += CfgMgr.BTagTrackAugmenter(
    "BTagTrackAugmenter_" + JetCollection,
    OutputLevel=INFO,
    JetCollectionName = JetCollection,
    TrackToVertexIPEstimator = ToolSvc.trkIPEstimator,
    SaveTrackVectors = True,
  )
  """

  alg.AccessBtagObject = True # for fatjets, turn this to False

  algSeq += alg

  from btagAnalysis.configHelpers import get_calibration_tool
  ToolSvc += get_calibration_tool(CfgMgr, JetCollection, False)

  # from btagAnalysis.configHelpers import get_calibration_tool_2016_calib
  # ToolSvc += get_calibration_tool_2016_calib(CfgMgr, JetCollection, False)


from PerfMonComps.PerfMonFlags import jobproperties as PerfMon_jp
PerfMon_jp.PerfMonFlags.doMonitoring = False
PerfMon_jp.PerfMonFlags.doFastMon = False

theApp.EvtMax = -1   #says how many events to run over. Set to -1 for all events
#theApp.EvtMax = 10

###########################################################################################################################################################################

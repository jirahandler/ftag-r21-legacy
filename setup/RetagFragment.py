##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
### THIS is the full retagging configuration
if doRetag:
  suffix = "retag"
  tracksKey = tracksOutputKey
  print 'RetagFragment:: TracksKey', tracksKey
  
  JetCollectionList = [ (JetCollection,
                         JetCollection.replace('ZTrack', 'Track').replace('PV0Track', 'Track'))
                        for JetCollection in JetCollections ]

  print 'JetCollectionList ',JetCollectionList
  #BTaggingFlags.CalibrationChannelAliases += [ "AntiKt4TopoEM->AntiKt4EMTopo" ]
  BTaggingFlags.CalibrationChannelAliases += [ "AntiKt4EMPFlow->AntiKt4EMPFlow" ]

  BTaggingFlags.Jets = [ name[1][:-4] for name in JetCollectionList]
  btag = "BTagging_"
  AuthorSubString = [ btag+name[1][:-4] for name in JetCollectionList]
  tmpSVname = "SecVtx"
  tmpJFVxname = "JFVtx"
  SA = 'standalone_'

  print 'AuthorSubString ', AuthorSubString

  from BTagging.BTaggingConfiguration import getConfiguration
  BTagConf = getConfiguration()
  BTagConf.doNotCheckForTaggerObstacles()
  NotInJetToolManager = [] # For jet collections

  from BTagging.BTaggingConf import Analysis__StandAloneJetBTaggerAlg as StandAloneJetBTaggerAlg

  for i, jet in enumerate(JetCollectionList):

    # set up algorithms to (re)match tracks and muons to jets
    from ParticleJetTools.ParticleJetToolsConf import JetAssocConstAlg
    from ParticleJetTools.ParticleJetToolsConf import JetParticleShrinkingConeAssociation, JetParticleFixedConeAssociation

    trackAssoc = JetParticleShrinkingConeAssociation(
      "DefaultBTaggingTrackAssoc",
      InputParticleCollectionName=tracksKey,
      OutputCollectionName="MatchedTracks" + '_' + suffix,
      coneSizeFitPar1=+0.239,
      coneSizeFitPar2=-1.220,
      coneSizeFitPar3=-1.64e-5
      )

    muonAssoc = JetParticleFixedConeAssociation(
      "DefaultBTaggingMuonAssoc",
      InputParticleCollectionName="Muons",
      OutputCollectionName="MatchedMuons" + '_' + suffix,
      coneSize=0.4,
      )

    assocalg = JetAssocConstAlg(
      "BTaggingRetagAssocAlg",
      JetCollections=JetCollections,
      Associators=[trackAssoc, muonAssoc]
      )

    algSeq += assocalg

    # set up the tool to link tracks to the b-tagging object
    from BTagging.BTaggingConf import Analysis__BTagTrackAssociation
    assoc = Analysis__BTagTrackAssociation(
      'thisBTagTrackAssociation_MANUAL_TOOL',
      AssociatedTrackLinks = "MatchedTracks" + '_' + suffix,
      AssociatedMuonLinks = "MatchedMuons" + '_' + suffix,
      TrackContainerName = tracksKey,
      MuonContainerName = "Muons",
      TrackAssociationName = "BTagTrackToJetAssociator",
      MuonAssociationName = "Muons",
      OutputLevel = VERBOSE
      )

    ToolSvc += assoc

    # set up the b-tagging tool itself
    btagger = BTagConf.setupJetBTaggerTool(ToolSvc = ToolSvc,
                                           JetCollection=jet[1][:-4],
                                           AddToToolSvc = True,
                                           Verbose = True,
                                           options = {"name"         : AuthorSubString[i].lower(),
                                                      "BTagName"     : AuthorSubString[i],
                                                      "BTagJFVtxName": tmpJFVxname,
                                                      "BTagSVName"   : tmpSVname,
                                                      "BTagTrackAssocTool" : assoc,
                                                      "OutputLevel" : INFO
                                                      }
                                           )

    SAbtagger = StandAloneJetBTaggerAlg(name=SA + AuthorSubString[i].lower(),
                                        JetBTaggerTool=btagger,
                                        JetCollectionName = jet[0],
                                        OutputLevel = INFO
                                        )

    # add the btagger algorithm to the sequencer; this is code that will be run upstream of the bTagAnalysis code!
    algSeq += SAbtagger
    print SAbtagger

  if len(NotInJetToolManager) > 0:
    AuthorSubString = list(set(AuthorSubString) - set(NotInJetToolManager))

  # Both standard and aux container must be listed explicitly. For release 19, the container version must be explicit.
  BaseName    = "xAOD::BTaggingContainer_v1#"
  BaseAuxName = "xAOD::BTaggingAuxContainer_v1#"
  #AOD list
  BTaggingFlags.btaggingAODList += [ BaseName + author for author in AuthorSubString]
  BTaggingFlags.btaggingAODList += [ BaseAuxName + author + 'Aux.' for author in AuthorSubString]
  BTaggingFlags.btaggingAODList += [ BaseName + author + 'AOD' for author in AuthorSubString]
  BTaggingFlags.btaggingAODList += [ BaseAuxName + author + 'AODAux.' for author in AuthorSubString]
  #ESD list
  BTaggingFlags.btaggingESDList += [ BaseName + author for author in AuthorSubString]
  BTaggingFlags.btaggingESDList += [ BaseAuxName + author + 'Aux.' for author in AuthorSubString]
  #AOD list SeCVert
  BaseNameSecVtx    = "xAOD::VertexContainer_v1#"
  BaseAuxNameSecVtx = "xAOD::VertexAuxContainer_v1#"
  BTaggingFlags.btaggingAODList += [ BaseNameSecVtx + author + tmpSVname for author in AuthorSubString]
  BTaggingFlags.btaggingAODList += [ BaseAuxNameSecVtx + author + tmpSVname + 'Aux.-vxTrackAtVertex' for author in AuthorSubString]
  #ESD list
  BTaggingFlags.btaggingESDList += [ BaseNameSecVtx + author + tmpSVname for author in AuthorSubString]
  BTaggingFlags.btaggingESDList += [ BaseAuxNameSecVtx + author + tmpSVname + 'Aux.-vxTrackAtVertex' for author in AuthorSubString]
  #AOD list JFSeCVert
  BaseNameJFSecVtx    = "xAOD::BTagVertexContainer_v1#"
  BaseAuxNameJFSecVtx = "xAOD::BTagVertexAuxContainer_v1#"
  BTaggingFlags.btaggingAODList += [ BaseNameJFSecVtx + author + tmpJFVxname for author in AuthorSubString]
  BTaggingFlags.btaggingAODList += [ BaseAuxNameJFSecVtx + author + tmpJFVxname + 'Aux.' for author in AuthorSubString]
  #ESD list
  BTaggingFlags.btaggingESDList += [ BaseNameJFSecVtx + author + tmpJFVxname for author in AuthorSubString]
  BTaggingFlags.btaggingESDList += [ BaseAuxNameJFSecVtx + author + tmpJFVxname + 'Aux.' for author in AuthorSubString]

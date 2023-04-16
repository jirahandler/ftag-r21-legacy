################################################################################
# Package: btagAnalysis
################################################################################

# Declare the package name:
atlas_subdir( btagAnalysis )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          GaudiKernel
                          PRIVATE
                          Control/AthenaBaseComps
                          DataQuality/GoodRunsLists
                          Event/xAOD/xAODBTagging
                          Event/xAOD/xAODEventInfo
                          Event/xAOD/xAODEventShape
                          Event/xAOD/xAODJet
                          Event/xAOD/xAODMuon
                          Event/xAOD/xAODTracking
                          Event/xAOD/xAODTrigger
                          Event/xAOD/xAODTruth
                          DetectorDescription/GeoPrimitives
                          InnerDetector/InDetRecTools/InDetTrackSelectionTool
                          InnerDetector/InDetRecTools/TrackVertexAssociationTool			  
                          PhysicsAnalysis/AnalysisCommon/ParticleJetTools
                          PhysicsAnalysis/AnalysisCommon/PileupReweighting
                          Reconstruction/Jet/JetCalibTools
                          Reconstruction/Jet/JetInterface
                          Reconstruction/RecoTools/ITrackToVertex
                          Reconstruction/Jet/JetSubStructureUtils
                          Trigger/TrigAnalysis/TrigDecisionTool
                          PhysicsAnalysis/DerivationFramework/DerivationFrameworkFlavourTag 
			  PhysicsAnalysis/TrackingID/InDetTrackSystematicsTools			  
)

# External dependencies:
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

atlas_add_dictionary(StlDict src/Stl.h src/selection.xml)

# Component(s) in the package:
atlas_add_component( btagAnalysis
                     btagAnalysis/*
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                     LINK_LIBRARIES ${ROOT_LIBRARIES} GaudiKernel AthenaBaseComps GoodRunsListsLib xAODBTagging xAODEventInfo xAODEventShape xAODJet xAODMuon xAODTracking xAODTrigger xAODTruth TrackVertexAssociationToolLib ParticleJetToolsLib JetSubStructureUtils PileupReweightingLib JetInterface ITrackToVertex TrkVertexFitterInterfaces TrigDecisionToolLib DerivationFrameworkFlavourTagLib InDetTrackSystematicsToolsLib)

# Install files from the package:
atlas_install_python_modules( python/*.py )


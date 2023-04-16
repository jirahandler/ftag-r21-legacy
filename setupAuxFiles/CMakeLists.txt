################################################################################
# Package: BTagTrackSystematicsAlgs
################################################################################

atlas_subdir( BTagTrackSystematicsAlgs )

atlas_depends_on_subdirs(
   PRIVATE
   Control/AthenaBaseComps
   GaudiKernel
   PhysicsAnalysis/TrackingID/InDetTrackSystematicsTools
   Control/CxxUtils
   Event/xAOD/xAODCore
   Event/xAOD/xAODTracking
)

# create the component library
atlas_add_component( BTagTrackSystematicsAlgs
   src/*.h src/*.cxx src/components/*.cxx
   INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${TBB_INCLUDE_DIRS}
   LINK_LIBRARIES AthenaBaseComps InDetTrackSystematicsToolsLib CxxUtils xAODCore xAODTracking
)

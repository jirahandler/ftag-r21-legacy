setupATLAS
lsetup git

mkdir retaggingFramework
cd retaggingFramework
mkdir source build run

cd source

git clone https://:@gitlab.cern.ch:8443/atlas-flavor-tagging-tools/FlavourTagPerformanceFramework.git
cd FlavourTagPerformanceFramework
git checkout freshstart
#source scripts/setup.sh
cd ..
mv FlavourTagPerformanceFramework/btagAnalysis .

git atlas init-workdir https://:@gitlab.cern.ch:8443/atlas/athena.git
cd athena
git checkout release/21.2.62.0
git atlas addpkg InDetTrackSystematicsTools
cd ..
mv athena/PhysicsAnalysis/TrackingID/InDetTrackSystematicsTools .
rm -rf athena

git clone https://:@gitlab.cern.ch:8443/atlas-ftag-calibration/BTagTrackSystematicsAlgs.git
cd BTagTrackSystematicsAlgs
cp ../../../setupAuxFiles/CMakeLists.txt .
cd ../../..

#cp setup/btagIBLAnalysisAlg.* retaggingFramework/source/FlavourTagPerformanceFramework/src/
cp setup/InDetTrackSystematics.h retaggingFramework/source/InDetTrackSystematicsTools/InDetTrackSystematicsTools/
cp setup/InDetTrackSmearingTool.h retaggingFramework/source/InDetTrackSystematicsTools/InDetTrackSystematicsTools/
cp setup/InDetTrackSmearingTool.cxx retaggingFramework/source/InDetTrackSystematicsTools/Root/

rm -rf retaggingFramework/source/btagAnalysis/share/jobOptions.py
rm -rf retaggingFramework/source/btagAnalysis/src/btagAnalysisAlg.cxx
rm -rf retaggingFramework/source/btagAnalysis/btagAnalysis/btagAnalysisAlg.h

cp setup/jobOptions.py retaggingFramework/source/btagAnalysis/share/.
cp setup/RetagFragment.py retaggingFramework/source/btagAnalysis/share/.
cp setup/SFCalculation.sh retaggingFramework/run/.
cp setup/btagAnalysisAlg.cxx retaggingFramework/source/btagAnalysis/src/.
cp setup/btagAnalysisAlg.h retaggingFramework/source/btagAnalysis/btagAnalysis/.


cd retaggingFramework/build
asetup AthDerivation,21.2.62.0,here
mv CMakeLists.txt ../source
cmake ../source
#make clean
make
#source x86_64-slc6-gcc62-opt/setup.sh
source x86_64-centos7-gcc62-opt/setup.sh

cd ../run
ln -s ../source/btagAnalysis/share/jobOptions.py .
ln -s ../source/btagAnalysis/share/RetagFragment.py .
ln -s ../source/btagAnalysis/share/MV2defaults.py .

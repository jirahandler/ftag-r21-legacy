export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

asetup AthDerivation,21.2.87.0,here

mkdir -p run

cd run

for JO in ../btagAnalysis/share/*.py; do
    ln -sf $JO
done

cd ..

if [[ ! -d build ]] ; then
    ./scripts/build.sh
else
    echo 'already built, run `./scripts/build.sh` to rebuild'
fi

source build/x86_*/setup.sh 

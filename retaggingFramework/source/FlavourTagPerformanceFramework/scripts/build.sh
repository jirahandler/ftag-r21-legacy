
mkdir -p build

cd build

cmake ..

make -j 4

source ./x86_*/setup.sh 

cd ..
# What does this package do?
It uses legacy CERN ATLAS FTAG tools to get flavor tagging information from FTAG1/FTAG2 DAODs in the form of ntuples.

##Setup Instructions
==================

First make an new directory and do:
```
    setupATLAS -q
	lsetup git
    git clone git@github.com:jirahandler/ftag-r21-legacy.git 
```
Then to setup the package for the first time (note the changes in `memo` file) :
```
    cd ftag-r21-legacy
    source memo

    #This is legacy setup, doesn't work now
	#source setup.sh

```
Once the setup is done, you are now in the `run/

Every time (except the first time which you just clone this package) you need to setup the ackage in the `run/` directory by doing:
```
    source memo
```

##To submit samples on the grid:
==============================

To submit the ntuple to grid, see https://gitlab.cern.ch/atlas-flavor-tagging-tools/FlavourTagPerformanceFramework/blob/freshstart/README.md for more details.

You need to do the submit under `run/` directory by first doing:

```
voms-proxy-init --voms atlas
lsetup rucio
lsetup panda
```
Then from the `run` directory edit the `../source/FlavourTagPerformanceFramework/scripts/grid-submit.sh` with the datasets you want to run over. 

Then from the `run` directory itself excecute `grid-submit.sh` :

```
../source/FlavourTagPerformanceFramework/scripts/grid-submit.sh
```

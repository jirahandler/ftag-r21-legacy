Setup Instructions
==================

```sh
git clone ssh://git@gitlab.cern.ch:7999/atlas-flavor-tagging-tools/FlavourTagPerformanceFramework.git
cd FlavourTagPerformanceFramework
git checkout freshstart
source scripts/setup.sh
```

to run:

Once you've built this package, it should contain a "run" directory

```sh
cd run
athena jobOptions.py
```

To submit samples on the grid:
==============================

edit the `scripts/grid_submit/grid-submit.sh` with the datasets you want to run over, and various options (see comments inside the script for details).

commit all changes (even if the change was only to `grid-submit.sh`) and create a tag
```sh
git tag <name of tag> -m'describe what the sample is for'
git push origin <name of tag>
```

(make sure to add a message to your tag, if you don't add one it won't be picked up by the grid script)

setup panda, and then excecute `grid-submit.sh` from the `run` directory (where the job options are)

```sh
./../scripts/grid_submit/grid-submit.sh

```






#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME} #(or /star/u/${LOGNAME} or /sphenix/u/${LOGNAME} or /eic/u/${LOGNAME})
 
source /opt/sphenix/core/bin/sphenix_setup.sh -n
# source /opt/sphenix/core/bin/setup_local.sh $HOME/install
source /opt/sphenix/core/bin/setup_local.sh $HOME/install_dual

# print the environment - needed for debugging
printenv

# this is how you run your Fun4All_G4_sPHENIX.C macro in batch: 
nevnt=$1
oname=$2
pythseed=$3


echo ../bin/main ${nevnt} ${oname}  ${pythseed}
../bin/main ${nevnt} ${oname}  ${pythseed}

echo all done

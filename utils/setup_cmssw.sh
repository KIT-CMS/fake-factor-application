#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

THIS_PWD=$PWD
cd CMSSW_10_2_14/src
eval `scramv1 runtime -sh`
cd $THIS_PWD

#!/bin/bash

# Export CMSSW
export SCRAM_ARCH=slc6_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project CMSSW CMSSW_10_2_14
cd CMSSW_10_2_14/src

# Setup fake factor tools
git clone ssh://git@github.com/CMS-HTT/Jet2TauFakes.git HTTutilities/Jet2TauFakes
cd HTTutilities/Jet2TauFakes
git checkout v0.2.2
scram b

# Get fake factor database
mkdir data_2016
git clone -b 2016 ssh://git@gitlab.cern.ch:7999/cms-htt/Jet2TauFakesFiles.git data_2016
mkdir data_2017
git clone -b 2017 ssh://git@gitlab.cern.ch:7999/cms-htt/Jet2TauFakesFiles.git data_2017
mkdir data_2018
git clone -b 2018 ssh://git@gitlab.cern.ch:7999/cms-htt/Jet2TauFakesFiles.git data_2018

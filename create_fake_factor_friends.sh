#!/bin/bash

ERA=$1
CONFIGKEY=$2
CATEGORYMODE=$3
OUTPUTDIR=$4

source fake-factor-application/utils/setup_cmssw.sh
source utils/setup_python.sh
source utils/setup_samples.sh $ERA

if [[ $ERA == *"2016"* ]]
then
    FF_database_ET=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2016/SM2016_ML/tight/et/fakeFactors_tight.root
    FF_database_MT=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2016/SM2016_ML/tight/mt/fakeFactors_tight.root
    FF_database_TT=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2016/SM2016_ML/tight/tt/fakeFactors_tight.root
    FF_workspce=fake-factor-application/htt_ff_fractions_2016.xroot
    USE_WORKSPACE=-w #change to -w to use fracctions from workspace
elif [[ $ERA == *"2017"* ]]
then
    FF_database_ET=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2017/SM2017/tight/vloose/et/fakeFactors.root
    FF_database_MT=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2017/SM2017/tight/vloose/mt/fakeFactors.root
    FF_database_TT=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2017/SM2017/tight/vloose/tt/fakeFactors.root
    FF_workspce=fake-factor-application/htt_ff_fractions_2017.xroot
    USE_WORKSPACE=-w
elif [[ $ERA == *"2018"* ]]
then
    FF_database_ET=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2018/SM2018/tight/vloose/et/fakeFactors.root
    FF_database_MT=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2018/SM2018/tight/vloose/mt/fakeFactors.root
    FF_database_TT=CMSSW_10_2_14/src/HTTutilities/Jet2TauFakes/data_2018/SM2018/tight/vloose/tt/fakeFactors.root
    FF_workspce=fake-factor-application/htt_ff_fractions_2017.xroot
    USE_WORKSPACE=-w
fi

python fake-factor-application/create_fake_factor_friends.py --era $ERA \
        -i $ARTUS_OUTPUTS \
        -o $OUTPUTDIR \
        --et-friend-directories $ARTUS_FRIENDS_ET \
        --mt-friend-directories $ARTUS_FRIENDS_MT \
        --tt-friend-directories $ARTUS_FRIENDS_TT \
        -c $CONFIGKEY \
        --et-fake-factor-directory $FF_database_ET \
        --mt-fake-factor-directory $FF_database_MT \
        --tt-fake-factor-directory $FF_database_TT \
        --num-threads 32 \
        --category-mode $CATEGORYMODE \
        --workspace $FF_workspce $USE_WORKSPACE

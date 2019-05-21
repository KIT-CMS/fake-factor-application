#!/bin/bash

ERA=$1
CONFIGKEY=$2
CATEGORYMODE=$3
OUTPUTDIR=$4

#./fake-factor-application/produce_shapes.sh $ERA $CONFIGKEY
./fake-factor-application/create_fake_factor_friends.sh $ERA $CONFIGKEY $CATEGORYMODE $OUTPUTDIR

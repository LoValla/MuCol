#!/bin/bash

source ../setup.sh

python ./pigun_generator.py \
--output="piguns.slcio" \
--numberOfEvents=5

ddsim --inputFile="piguns.slcio" \
--outputFile="output_sim.slcio" \
--steeringFile="steer_sim.py"


#ln -s ../mucoll-benchmarks/reconstruction/marlin/PandoraSettings ./
ACTS_PATH=$(echo $MARLIN_DLL | tr ':' '\n' | grep actstracking | sed "s:/lib.*::")

export MARLIN_DLL=$MARLIN_DLL:\
/eos/user/l/lvalla/MuColl/TrackPerf/build/libTrackPerf.so

mkdir ./Outputs

Marlin ../mucoll-benchmarks/digirecolctuple.xml \
--global.LCIOInputFiles="Outputs/output_sim.slcio" \
--DD4hep.DD4hepXMLFile="$MUCOLL_GEO"


#Marlin ../mucoll-benchmarks/lctupleandanalysis.xml \
#--global.LCIOInputFiles="Outputs/output_reco.slcio" \
#--DD4hep.DD4hepXMLFile="$MUCOLL_GEO"

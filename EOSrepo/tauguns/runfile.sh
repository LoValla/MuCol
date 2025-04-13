#!/bin/bash

ln -s /cvmfs/muoncollider.cern.ch/release/2.8-patch2/setup.sh ./
source setup.sh

python EOSrepo/taugun_generator.py \
--output="tauguns.slcio" \
--numberOfEvents=500

ddsim --inputFile="tauguns.slcio" \
--outputFile="output_sim.slcio" \
--steeringFile="EOSrepo/steer_sim.py"

ln -s EOSrepo/mucoll-benchmarks/reconstruction/marlin/PandoraSettings ./

Marlin EOSrepo/mucoll-benchmarks/digirecolctuple.xml \
--global.LCIOInputFiles="output_sim.slcio" \
--DD4hep.DD4hepXMLFile="$MUCOLL_GEO"
#!/bin/bash 
ddsim --steeringFile steer_baseline.py
echo "EXECUTION OF DDSIM COMPLETED"

Marlin ../mucoll-benchmarks/digitisation/marlin/digi_steer_BIB.xml \
--global.LCIOInputFiles="Outputs/output_sim.slcio" \
--DD4hep.DD4hepXMLFile="$MUCOLL_GEO" 
echo "DIGITISATION COMPLETED"

ln -s ../mucoll-benchmarks/reconstruction/marlin/PandoraSettings ./
ACTS_PATH=$(echo $MARLIN_DLL | tr ':' '\n' | grep actstracking | sed "s:/lib.*::")
Marlin /eos/user/l/lvalla/MuColl/mucoll-benchmarks/reconstruction/marlin/reco_steer.xml \
--CKFTracking.MatFile="${ACTS_PATH}/share/ACTSTracking/data/material-maps.json" \
--CKFTracking.TGeoFile="${ACTS_PATH}/share/ACTSTracking/data/MuColl_v1.root" \
--DD4hep.DD4hepXMLFile="$MUCOLL_GEO"
echo "RECONSTRUCTION COMPLETED"

export MARLIN_DLL=$MARLIN_DLL:\
/eos/user/l/lvalla/MuColl/MyTauFinder/build/libMyTauFinder.so

Marlin --global.LCIOInputFiles="Outputs/output_reco.slcio" \
/eos/user/l/lvalla/MuColl/MyTauFinder/share/MyTauFinder.xml

Marlin --global.LCIOInputFiles="Outputs/output_taufinder.slcio" \
/eos/user/l/lvalla/MuColl/MyTauFinder/share/MyEvaluateTauFinderGun.xml

root -l -b ../TauAnalysis.cxx

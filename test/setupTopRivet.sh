#!/bin/bash
source $CMSSW_BASE/src/GeneratorInterface/RivetInterface/test/rivetSetup.sh

RIVETTOP=$CMSSW_BASE/src/GeneratorInterface/RivetTop
for I in `\ls $RIVETTOP/data/*.yoda`; do
    F=`basename $I`
    cd $CMSSW_BASE/src/GeneratorInterface/RivetInterface/data
    if [ ! -f $F ]; then 
      cp $RIVETTOP/data/$F .
    fi
    cd -
done

export RIVET_REF_PATH=$RIVET_REF_PATH:$CMSSW_BASE/src/GeneratorInterface/RivetTop/data
export RIVET_INFO_PATH=$RIVET_INFO_PATH:$CMSSW_BASE/src/GeneratorInterface/RivetTop/data
export RIVET_PLOT_PATH=$RIVET_PLOT_PATH:$CMSSW_BASE/src/GeneratorInterface/RivetTop/data

which yodamerge &> /dev/null || GETYODA=1
if [ $GETYODA -eq 1 ]; then
  eval `scram tool info yoda | grep YODA_BASE`
  export PATH=$PATH:$YODA_BASE/bin
fi

#!/bin/bash
source $CMSSW_BASE/src/GeneratorInterface/RivetInterface/test/rivetSetup.sh

RIVETTOP=$CMSSW_BASE/src/GeneratorInterface/RivetTop
for I in `\ls $RIVETTOP/data/*.yoda`; do
    F=`basename $I`
    cd $CMSSW_BASE/src/GeneratorInterface/RivetInterface/data
    if [ ! -f $F ]; then 
      echo ln -s $RIVETTOP/data/$F
      ln -s $RIVETTOP/data/$F
    fi
    cd -
done

export RIVET_REF_PATH=$RIVET_REF_PATH:$CMSSW_BASE/src/GeneratorInterface/RivetTop/data
export RIVET_INFO_PATH=$RIVET_INFO_PATH:$CMSSW_BASE/src/GeneratorInterface/RivetTop/data
export RIVET_PLOT_PATH=$RIVET_PLOT_PATH:$CMSSW_BASE/src/GeneratorInterface/RivetTop/data

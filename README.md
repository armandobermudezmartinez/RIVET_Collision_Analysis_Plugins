# RivetTop

This GitLab repository contains Rivet routines for Top PAG related analyses. Some may not be validated.

## Organisation

The `master` branch contains all the plugins and data. Configuration files are in the `test/` repository, while `.yoda` and `.plot` (and possibly `.info`) files are in the `data/` repository.

## Installation

<pre>
cmsrel CMSSW_7_5_0
cd CMSSW_7_5_0/src
cmsenv
git-cms-init
git-cms-addpkg GeneratorInterface/RivetInterface
git-cms-merge-topic jhgoh:RivetConsumesMigration75
git-cms-merge-topic jhgoh:RivetRefHistFromEnvVar
git-cms-merge-topic jhgoh:LHEweight
mkdir TopMonteCarlo
cd TopMonteCarlo
git clone ssh://git@gitlab.cern.ch:7999/CMS-TOP-Rivet/RivetTop
git clone ssh://git@gitlab.cern.ch:7999/CMS-TOP-Rivet/Configuration
cd RivetTop/test
source setupTopRivet.sh
cd $CMSSW_BASE/src
scram b -j8
</pre>

---------------------------------------

# The Les Houches plugin

## Description 

This plugin has been developped for the Les Houches ttbar working group. 

It fullfills the requirements from:

http://phystev.cnrs.fr/wiki/2015:groups:tools:ttjets

More details can be found here:

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TOPRivetForLesHouches

and the result there:

http://ebouvier.web.cern.ch/ebouvier/TOPRivetForLesHouches/plots/CMS_AN_PseudoTop/index.html



## How-to-setup (and run) 

    scramv1 project -n CMSSW_7_5_0-GitLab CMSSW_7_5_0
    cd CMSSW_7_5_0-GitLab/src/
    cmsenv
    git-cms-init
    git-cms-addpkg GeneratorInterface/RivetInterface
    git-cms-addpkg /Configuration/Generator
    wget -P Configuration/Generator/python/ http://test-efeafs.web.cern.ch/test-efeafs/Hadronizer_pythia8_cff.py 
    cd GeneratorInterface/
    git clone ssh://git@gitlab.cern.ch:7999/CMS-TOP-Rivet/RivetTop
    cd ..
    scram b -j6
    cd GeneratorInterface/RivetTop/test/
    cmsRun runRivietWithPythia8.py
    cmsRun runRivietWithHerwig.py
    rivet-mkhtml -c ../data/CMS_LesHouches2015.plot Pythia8.yoda:'Powheg+Pythia 8' Herwig.yoda:'Powheg+Herwig++'
    firefox plots/index.hmtl &
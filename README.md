# RivetTop

This GitLab repository contains Rivet routines for Top PAG related analyses. Some may not be validated.

## Organisation

The `master` branch contains all the plugins and data.
Configuration files are in the `test/` repository, while `.yoda` and `.plot` (and possibly `.info`) files are in the `data/` repository.

Common interfaces to be used in analysis modules are available.

  * PartonTop : Build parton level top quarks and its decay products.
  * PseudoTop : Build particle level top quarks (pseudo top) and its decay products.
  * CMSGenParticle (not validated) : Emulate CMS genParticlesForJets in Run-I

Some analyses done in CMS are also available.

  * TOP-12-028 (CMS\_2015\_I1370682.cc) : Differential cross section measurement, under validation
  * TOP-12-042 (CMS\_TOP\_12\_042.cc) : Jet multiplicity measurement, under validation

Other modules are also available for the MC studies

  * Les Houches working group (CMS\_LesHouches2015.cc) : Analysis for the Les Houches ttbar working group using the dilepton channel (see the section below)
  * MC\_TTBAR\_HADRON.cc : Lepton+jets channel study

## Installation

    cmsrel CMSSW_7_5_0
    cd CMSSW_7_5_0/src
    cmsenv

    git-cms-init
    git-cms-addpkg GeneratorInterface/RivetInterface
    git-cms-addpkg Configuration/Generator
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

---------------------------------------

## The Les Houches plugin
### Description
This plugin has been developped for the Les Houches ttbar working group.

It fullfills the requirements from:

http://phystev.cnrs.fr/wiki/2015:groups:tools:ttjets

More details can be found here:

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TOPRivetForLesHouches

and the result there:

http://ebouvier.web.cern.ch/ebouvier/TOPRivetForLesHouches/plots/CMS\_AN\_PseudoTop/index.html

### How-to-setup (and run)

    wget -P $CMSSW_BASE/src/Configuration/Generator/python http://test-efeafs.web.cern.ch/test-efeafs/Hadronizer_pythia8_cff.py
    cd $CMSSW_BASE/src/GeneratorInterface/RivetTop/test
    cmsRun runRivietWithPythia8.py
    cmsRun runRivietWithHerwig.py
    rivet-mkhtml -c ../data/CMS_LesHouches2015.plot Pythia8.yoda:'Powheg+Pythia 8' Herwig.yoda:'Powheg+Herwig++'
    firefox plots/index.html &
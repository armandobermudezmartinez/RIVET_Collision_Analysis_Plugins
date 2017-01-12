# CMS-TOP-Rivet

This GitLab repository contains Rivet routines for Top PAG related analyses and configuration files to run them. 
Routines are gather in the `CMS-TOP-Rivet/RivetTop` subpackage. Some may not be validated.
Configuration files to run them are in the `CMS-TOP-Rivet/Configuration` subpackage.
The `Installation` section gather instructions to get both subpackages.

## Installation

Many of the plugins in `CMS-TOP-Rivet/RivetTop` and configuration files in `CMS-TOP-Rivet/Configuration` have been developed in previous versions of CMSSW and Rivet. 
Please, be aware that bugs could have been introduced while porting to newer versions.

    cmsrel CMSSW_8_1_0
    cd CMSSW_8_1_0/src
    cmsenv

    git-cms-init
    git-cms-addpkg GeneratorInterface/RivetInterface
    git-cms-addpkg Configuration/Generator
    git-cms-merge-topic jhgoh:TOP-RIVET-80X

    mkdir TopMonteCarlo
    cd TopMonteCarlo
    git clone https://:@gitlab.cern.ch:8443/CMS-TOP-Rivet/RivetTop.git
    git clone https://:@gitlab.cern.ch:8443/CMS-TOP-Rivet/Configuration.git

    cd RivetTop/test
    source setupTopRivet.sh
    cd $CMSSW_BASE/src
    scram b -j8

## Organisation of `CMS-TOP-Rivet/RivetTop` 

The `master` branch contains all the plugins and data.
The `scr/` and `interface/` repositories contain `.cc` and `.h` files respectively for these plugins.
The `data/` repository contains `.yoda` (for data and MC) and `.plot` files.
The `test/` repository contains only small scripts that run on `.yoda` files (for corrections, normalization...).

Common interfaces to be used in analysis modules are available.

  * PartonTop : Build parton level top quarks and its decay products.
  * PseudoTop : Build particle level top quarks (pseudo top) and its decay products.
  * CMSGenParticle (not validated) : Emulate CMS genParticlesForJets in Run-I

Some analyses done in CMS are also available.

  * TOP-12-028 (CMS\_2015\_I1370682.cc) : Differential cross section measurement, under validation
  * TOP-12-042 (CMS\_TOP\_12\_042.cc) : Differential cross section using event variables, under validation

Analysis plugins that rely on parton level information are developed for internal use and correction function extraction.
  * TOP-12-028 (CMS\_2015\_I1370682\_internal.cc) : Differential cross section (synchronized with MadGraph in the paper)
  * TOP-12-041 (CMS\_TOP\_12\_041\_internal.cc) : Jet multiplicity measurement (under validation)

Other modules are also available for the MC studies

  * Les Houches working group (CMS\_LesHouches2015.cc) : Analysis for the Les Houches ttbar working group using the dilepton channel (see the section below)
  * MC\_TTBAR\_HADRON.cc : Lepton+jets channel study

## Organisation of `CMS-TOP-Rivet/Configuration`

The `analysis/` repository contains configuration file for MC generators.
The `plugins/` repository contains `.cc` and `.h` files for plugins.
The `python/` repository contains `.py` files to run plugins and Rivet routines.
The `crab/` repository contains configuration files to run Crab3

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

    cp $CMSSW_BASE/src/TopMonteCarlo/Configuration/analysis/Hadronizer_pythia8_cff.py $CMSSW_BASE/src/TopMonteCarlo/Configuration/python/Hadronizer_pythia8_cff.py
    cd $CMSSW_BASE/src/TopMonteCarlo/Configuration/python
    cmsRun runRivietWithPythia8.py
    cmsRun runRivietWithHerwig.py
    rivet-mkhtml -c $CMSSW_BASE/src/TopMonteCarlo/RivetTop/data/CMS_LesHouches2015.plot Pythia8.yoda:'Powheg+Pythia 8' Herwig.yoda:'Powheg+Herwig++'
    firefox plots/index.html &

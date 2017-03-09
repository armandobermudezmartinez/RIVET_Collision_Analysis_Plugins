# CMS Rivet repository

This GitLab repository contains Rivet routines for CMS analyses and sample configuration files to run them.

Please consult the README files in the subdirectories for information about the plugins contained there.

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

    git clone https://:@gitlab.cern.ch:8443/CMS-TOP-Rivet/RivetTop.git Rivet

    source Rivet/rivetSetup.sh
    scram b -j8

Hint: this repository will be renamed.

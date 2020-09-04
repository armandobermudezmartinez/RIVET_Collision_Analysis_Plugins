# CMS Rivet repository

This GitLab repository contains Rivet routines for CMS analyses and sample configuration files to run them.

Please consult the README files in the subdirectories for information about the plugins contained there.

Please follow the [contribution guide](CONTRIBUTING.md) for developing your plugins.

## Installation

First, create a personal fork of this repository: https://gitlab.cern.ch/cms-gen/Rivet/forks/new

    cmsrel CMSSW_10_6_0
    cd CMSSW_10_6_0/src
    cmsenv

    git-cms-init
    git-cms-addpkg GeneratorInterface/RivetInterface
    git-cms-addpkg Configuration/Generator

    git clone ssh://git@gitlab.cern.ch:7999/${USER}/Rivet.git
    cd Rivet
    git remote add cms-gen ssh://git@gitlab.cern.ch:7999/cms-gen/Rivet.git
    git fetch cms-gen master-rivet2
    git checkout master-rivet2
    git pull cms-gen master-rivet2

    source Rivet/rivetSetup.sh
    scram b -j8


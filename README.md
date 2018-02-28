# CMS Rivet repository

This GitLab repository contains Rivet routines for CMS analyses and sample configuration files to run them.

Please consult the README files in the subdirectories for information about the plugins contained there.

Please follow the [contribution guide](CONTRIBUTING.md) for developing your plugins.

## Installation

    cmsrel CMSSW_10_0_0
    cd CMSSW_10_0_0/src
    cmsenv

    git-cms-init
    git-cms-addpkg GeneratorInterface/RivetInterface
    git-cms-addpkg Configuration/Generator

    git clone https://:@gitlab.cern.ch:8443/${USER}/Rivet.git

    source Rivet/rivetSetup.sh
    scram b -j8


# CMS Rivet repository

This GitLab repository contains Rivet routines for CMS analyses and sample configuration files to run them.

Please consult the README files in the subdirectories for information about the plugins contained there.

Please follow the [contribution guide](CONTRIBUTING.md) for developing your plugins.

## Installation

    cmsrel CMSSW_8_1_0
    cd CMSSW_8_1_0/src
    cmsenv

    git-cms-init
    git-cms-addpkg GeneratorInterface/RivetInterface
    git-cms-addpkg Configuration/Generator
    git-cms-merge-topic jhgoh:TOP-RIVET-80X

    git clone https://:@gitlab.cern.ch:8443/cms-gen/Rivet.git

    source Rivet/rivetSetup.sh
    scram b -j8

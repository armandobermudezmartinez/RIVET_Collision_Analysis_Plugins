== Description ==

This branch contains a plugin developped for the Les Houches ttbar working group. It fullfills the requirements from:
http://phystev.cnrs.fr/wiki/2015:groups:tools:ttjets
More details can be found here:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/TOPRivetForLesHouches
and the result there:
http://ebouvier.web.cern.ch/ebouvier/TOPRivetForLesHouches/plots/CMS_AN_PseudoTop/index.html


== How-to-setup (and run) ==

    scramv1 project -n CMSSW_7_5_0-GitLab CMSSW_7_5_0
    cd CMSSW_7_5_0-GitLab/src/
    cmsenv
    git-cms-init
    git-cms-addpkg GeneratorInterface/RivetInterface
    git-cms-addpkg /Configuration/Generator
    wget -P Configuration/Generator/python/ http://test-efeafs.web.cern.ch/test-efeafs/Hadronizer_pythia8_cff.py 
    cd GeneratorInterface/
    git clone ssh://git@gitlab.cern.ch:7999/CMS-TOP-Rivet/RivetTop
    cd RivetTop/
    git checkout -b LesHouches2015 origin/LesHouches2015
    cd ../../
    scram b -j6
    cd GeneratorInterface/RivetTop/test/
    cmsRun runRivietWithPythia8.py
    cmsRun runRivietWithHerwig.py
    rivet-mkhtml -c ../data/CMS_LesHouches2015.plot Pythia8.yoda:'Powheg+Pythia 8' Herwig.yoda:'Powheg+Herwig++'
    firefox plots/index.hmtl &



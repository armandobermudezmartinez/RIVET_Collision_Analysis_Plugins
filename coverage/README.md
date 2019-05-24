# Rivet CMS analysis coverage

Current coverage: https://cms-rivet.web.cern.ch/cms-rivet/rivet-coverage-cms.html

## Instructions for updating

API token needed for including merge requests, get it here: https://gitlab.cern.ch/profile/personal_access_tokens

    get-marcxml-inspire-cms
    get-rivethd-marcxml inspire-cms-*.marc.xml
    wget --no-check-certificate https://bitbucket.org/heprivet/rivet/raw/release-2-7-x/doc/rivet-coverage-cms.rank
    mk-coverage-html-cms inspire-cms-*.json -r rivet-coverage-cms.rank -R --token <you personal token>
    cp rivet-coverage-cms.html /eos/project/c/cmsweb/www/generators/Rivet/
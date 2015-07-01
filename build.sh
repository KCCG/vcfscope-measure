#!/bin/bash
set -x -e -o pipefail

function dxlogin {
  source $HOME/dx-toolkit/environment
  dx login --noprojects --token=${bamboo_DX_TOKEN}
}

function mkproject {
  proj_id=`dx new project -s ${PROJ_NAME} --brief`
  dx api ${proj_id} invite '{"invitee":"org-garvan_kccg_developers", "level":"VIEW", "suppressEmailNotification":true}'
}

function rmproject {
  dx rmproject -y $PROJ_NAME
}

# Deploy, but do not publish an app, for testing
function build_smoketest {
  dx build . --version bamboo_smoketest --app --yes
}

# Run a basic test of a deployed app.
# The tests themselves are run by the script (-iruntests=true),
# and are housed in resources/home/dnanexus/test-calcs.R
function run_smoketest {
  vcfgz=project-BVJz7k0098GX43GV9ZPFXVY2:file-Bf9b5qj098GppvyPX1P2K8Kf
  region=project-BVJz7k0098GX43GV9ZPFXVY2:file-Bf9b5zQ098Ggq5zZ3kf846ZP
  jobid=$(dx run app-kccg-validation-reporter/bamboo_smoketest -ivcfgz=$vcfgz -iregion=$region -iextended=true -iruntests=true --yes --brief)

  # returns non-0 if there's a problem.
  dx watch -q --no-job-info  -f '{msg}' ${jobid}
  dx describe ${jobid}

  echo -n "Total job price: "
  dx describe --json ${jobid} | jq .totalPrice
}

# Deploy and publish an app. This is irreversible (but the release can be deprecated).
function publish_app {
  if [[ -n "${VERSION}" ]]; then
    dx build . --version ${VERSION} --no-version-autonumbering --no-update --app --publish --yes
  else
    echo >&2 "Error, $$VERSION not set. Looks like you are trying to deploy from a build agent"
    exit 1
  fi
}


if [[ $1 = "smoketest" ]]; then 
  dxlogin
  mkproject
  build_smoketest
  run_smoketest
elif [[ $1 = "publish" ]]; then
  dxlogin
  publish_app
elif [[ $1 = "clean" ]]; then
  dxlogin
  rmproject
fi

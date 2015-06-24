#!/bin/bash
set -x -e -o pipefail

function dxlogin {
  source $HOME/dx-toolkit/environment
  dx login --noprojects --token=${bamboo_DX_TOKEN}
}

# Deploy, but do not publish an app, for testing
function build_smoketest_app {
  dx build . --version smoketest --app --yes
}

# Deploy and publish an app. This is irreversible (but the release can be deprecated).
function publish_app {
  if [[ -n "${VERSION}" ]]; then
    dx build . --version ${VERSION} --app --publish --yes
  else
    echo >&2 "Error, $$VERSION not set. Looks like you are trying to deploy from a build agent"
    exit 1
  fi
}

# Run a basic test of a deployed app.
function run_smoketest {
  vcfgz=project-Bb9KVk8029vp1qzXz4yx4xB3:file-BbBg1g00z1f8vYJq606Yjv7K
  region=project-BZ4JvjQ0K74XK3bP71gykXKQ:file-Bf4JFV00K74xkqBP2Qzbq1p4
  jobid=$(dx run app-kccg-validation-reporter/smoketest -ivcfgz=$vcfgz -iregion=$region -iextended=true --yes --brief)

  #
  # consistent to all applets
  #
  # returns non-0 if there's a problem.
  dx watch -q --no-job-info  -f '{msg}' ${jobid}
  dx describe ${jobid}

  echo -n "Total job price: "
  dx describe --json ${jobid} | jq .totalPrice
}


if [[ $1 = "smoketest" ]]; then 
  dxlogin
  build_smoketest_app
  run_smoketest
elif [[ $1 = "publish" ]]; then
  dxlogin
  publish_app
fi

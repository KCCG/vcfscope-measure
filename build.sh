#!/bin/bash
# test harness to buid an applet, and run a small test
# usage $0 <run|deploy|rm>

# exit on any error
set -x -e -o pipefail

function dxlogin {
#  which dx
  source $HOME/dx-toolkit/environment
  dx login --noprojects --token=${bamboo_DX_TOKEN}
}

function create_project {
  dx new project -s ${PROJ_NAME}
}

# build an applet.
function build_applet {
  dx build .
  dx mkdir /assets
  dx cp kccg:/assets/kccg_validation_reporter_resources_bundle-1.0.tar.gz /assets/
  dx cp kccg:/assets/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz /assets/
}

# deploy and publish an applet. - ie install it to kccg
function publish_applet {
  dx select kccg
  dx build -f
}

# deploy, but do not publish an app. This can be edited prior to publish()ing
function build_app {
  dx build . --version 9.9.9 --app --yes
}

# deploy and publish an app. This is irreversible (but the release can be deprecated).
function publish_app {
  if [[ -n "${VERSION}" ]]; then
    dx build . --version ${VERSION} --app --publish --yes
  else
    echo >&2 "Error, $$VERSION not set. Looks like you are trying to deploy from a build agent"
    exit 1
  fi
}

# Note, there are more tests within test/smoke_tests.py
function smoketest {
  vcfgz=project-Bb9KVk8029vp1qzXz4yx4xB3:file-BbBg1g00z1f8vYJq606Yjv7K
  region=project-BZ4JvjQ0K74XK3bP71gykXKQ:file-Bf4JFV00K74xkqBP2Qzbq1p4
  jobid=$(dx run /kccg-validation-reporter-dx -ivcfgz=$vcfgz -iregion=$region --yes --brief)

  #
  # consistent to all applets
  #
  # returns non-0 if there's a problem.
  dx watch -q --no-job-info  -f '{msg}' ${jobid}
  dx describe ${jobid}

  echo -n "Total job price: "
  dx describe --json ${jobid} | jq .totalPrice
}


function rmproject {
   dx rmproject -y $PROJ_NAME
}

if [[ $1 = "run" ]]; then 
  dxlogin
  create_project
  build_applet
  smoketest
elif [[ $1 = "publish" ]]; then
  dxlogin
  publish_applet
elif [[ $1 = "rm" ]]; then
  dxlogin
  rmproject
fi

########################################
# Creating BAMBOO Build & Deploy Plans #
########################################
# SOP for Creating a Build Plan
# Create : Clone an existing Deploy plan: choose kccg-bwa-mem's deploy plan
# Edit stash repository
# Double check the stages
# run.

# SOP for Creating a Deployment Project
# Creating a Build Plan
# Environment: DNANexus
#
# 1) Clean Working Dir
#
# 2) Artifact Download
#    Artifact Name: app
#
# 3) Script
#    Task description: make executable
#    Script body: chmod +x resources/usr/bin/*
#
# 4) Command
#    Task description: publish app
#    Argument: bash build.sh publish
#    Environment variables: VERSION=${bamboo.deploy.version}

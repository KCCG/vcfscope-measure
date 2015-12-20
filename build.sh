#!/bin/bash
# test harness to buid an applet, and run a small test
# usage $0 <run|deploy|rm>

# exit on any error
set -x -e -o pipefail

function dxlogin {
  source $HOME/dx-toolkit/environment
  dx login --noprojects --token=${bamboo_DX_TOKEN}
}

function createproject {
    proj_id=`dx new project -s ${PROJ_NAME} --brief`
    dx api ${proj_id} invite '{"invitee":"org-garvan_kccg_developers", "level":"VIEW", "suppressEmailNotification":true}'
}

# deploy, but do not publish an app. This can be edited prior to publish()ing
function build_app {
  dx build . --version 9.9.7 --app --yes
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

function rmproject {
  dx rmproject -y $PROJ_NAME
}

if [[ $1 = "run" ]]; then 
    dxlogin
    createproject
    build_app
elif [[ $1 = "publish" ]]; then
    dxlogin
    publish_app
elif [[ $1 = "rm" ]]; then
    dxlogin
    rmproject
fi

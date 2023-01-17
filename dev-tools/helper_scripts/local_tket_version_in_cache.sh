#!/usr/bin/env bash

# This script returns yes if the tket requirement of the given conanfile is present in cache or remote
#   First argument: Path to conanfile

set -euo pipefail

CONANFILE_PATH=$1
REQUIRED_TKET_VERSION="tket/$(cat ${CONANFILE_PATH} | grep version | egrep -o "\d+\.\d+\.\d+")@tket/stable"
conan get ${REQUIRED_TKET_VERSION} &> /dev/null && echo "true" && exit
echo "false"

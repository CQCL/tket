#!/usr/bin/env bash

# This script updates the tket version locally in the corresponding conan configuration
# files for tket, tests, proptests and pytket
#   First argument: Path to top-level repository directory
#   Second argument: Semantic version to set or 'next'

set -euo pipefail

TKET_SRC_DIR=$1
NEW_TKET_VERSION=$2

CURRENT_TKET_VERSION=$(cat ${TKET_SRC_DIR}/recipes/tket/conanfile.py  | grep version | grep -E -o "[0-9]+\.[0-9]+\.[0-9]+")

set_local_tket_version() {
  echo "Setting local tket version to: ${NEW_TKET_VERSION}"

  sed -i'' -e "s/${CURRENT_TKET_VERSION}/${NEW_TKET_VERSION}/" "${TKET_SRC_DIR}/recipes/tket/conanfile.py"
  sed -i'' -e "s/tket\/.*\@/tket\/${NEW_TKET_VERSION}\@/" "${TKET_SRC_DIR}/recipes/tket-tests/conanfile.py"
  sed -i'' -e "s/tket\/.*\@/tket\/${NEW_TKET_VERSION}\@/" "${TKET_SRC_DIR}/recipes/tket-proptests/conanfile.py"
  sed -i'' -e "s/tket\/.*\@/tket\/${NEW_TKET_VERSION}\@/" "${TKET_SRC_DIR}/pytket/conanfile.txt"

}

if [[ "${NEW_TKET_VERSION}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] # if given version is a three part semver
then
  set_local_tket_version
  exit
elif [[ "${NEW_TKET_VERSION}" =~ ^[n|N][e|E][x|X][t|T]$ ]] # if given version is 'next' (capitalization doesn't matter)
then
  NEXT_TKET_PATCH=$(awk -F. '/[0-9]+\./{$NF++;print}' OFS=. <<< ${CURRENT_TKET_VERSION})
  NEW_TKET_VERSION=$NEXT_TKET_PATCH
  set_local_tket_version
  exit
fi

echo "Error: Given version ${NEW_TKET_VERSION} is not a three-part semver or 'next'" >> /dev/stderr
exit 1

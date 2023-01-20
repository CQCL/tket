# This makefile bundles commands for streamlined development
# Use `make help` for an overview of commands
SHELL=/usr/bin/env bash

build_type=Debug
cmake_build_dir=cmake-build-debug
conan_profile_name=tket
ifeq ($(build_type), Release)
  cmake_build_dir=cmake-build-release
  conan_profile_name=tket-release
endif
tket_remote_repository=https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
libs_dir=../libs
tket_recipe_dir=../recipes/tket
tket_tests_recipe_dir=../recipes/tket-tests
tket_proptests_recipe_dir=../recipes/tket-proptests
cmake_source_dir=.

n_cpus=$(shell python -c 'import multiprocessing as m; print(m.cpu_count())')

##@ General

.PHONY: help
help: ## Display this help
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_0-9-]+:.*?##/ { printf "  \033[36m%-25s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Setup Conan Profile

.PHONY: conan-profile
conan-profile: ## Create and configure conan profile
	-conan profile new $(conan_profile_name) --detect
	-conan remote add tket-libs $(tket_remote_repository)
	conan config set general.revisions_enabled=1
	conan profile update options.tket:shared=True $(conan_profile_name)
	conan profile update options.tklog:shared=True $(conan_profile_name)
	conan profile update settings.build_type=$(build_type) $(conan_profile_name)

.PHONY: linux-profile-settings
linux-profile-settings: ## Set suggested settings for conan profile on linux
	conan profile update settings.compiler.libcxx=libstdc++11 $(conan_profile_name) # see warning about gcc ABI compatibility at profile creation

##@ Install Dependencies and Setup Cmake

 # Currently the b2 package (a dependency of boost) in the remote tket-libs does not include compiler info
 # which can lead to errors when building boost from source. Conan will try to use
 # the remote binary even if it is incompatible with the compiler, leading to linking issues.
 # This fixes that problem by forcing a b2 build from source
.PHONY: b2
b2: # install and build b2 (boost build requirements) from source
	conan install --profile=$(conan_profile_name) "b2/[>4.0.0]@" --build

.PHONY: local-libs
local-libs: # create conan packages for local libs from this repository
	conan create --profile=$(conan_profile_name) $(libs_dir)/tklog tket/stable --build=missing
	conan create --profile=$(conan_profile_name) $(libs_dir)/tkassert tket/stable --build=missing
	conan create --profile=$(conan_profile_name) $(libs_dir)/tkrng tket/stable --build=missing
	conan create --profile=$(conan_profile_name) $(libs_dir)/tktokenswap tket/stable --build=missing
	conan create --profile=$(conan_profile_name) $(libs_dir)/tkwsm tket/stable --build=missing

.PHONY: package-local-tket-if-necessary # if the required tket version package is not in cache, package it (skip test phase with --test-folder=None)
required_tket_version_available=$(shell ./helper_scripts/local_tket_version_in_cache.sh $(tket_recipe_dir)/conanfile.py)
package-local-tket-if-necessary:
ifeq ($(required_tket_version_available), false)
	-conan create --profile=$(conan_profile_name) $(tket_recipe_dir) tket/stable --test-folder=None
endif

.PHONY: dev-env-from-source
dev-env-from-source: b2 local-libs ## Install tket, tket-tests, and tket-proptests into a local build directory and configure cmake (build from source if not in conan cache or remotes)
	conan install --profile=$(conan_profile_name) $(tket_recipe_dir)  --install-folder=$(cmake_build_dir)/tket/src --build=missing
	make package-local-tket-if-necessary
	conan install --profile=$(conan_profile_name) $(tket_tests_recipe_dir) --install-folder=$(cmake_build_dir)/tket/tests --build=missing
	conan install --profile=$(conan_profile_name) $(tket_proptests_recipe_dir) --install-folder=$(cmake_build_dir)/tket/proptests --build=missing
	cmake -DCMAKE_BUILD_TYPE=$(build_type) -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S $(cmake_source_dir) -B $(cmake_build_dir)

.PHONY: dev-env
dev-env: ## Install tket, tket-tests, and tket-proptests into a local build directory and configure cmake (dependencies from cache or remote)
	conan install --profile=$(conan_profile_name) $(tket_recipe_dir)  --install-folder=$(cmake_build_dir)/tket/src
	make package-local-tket-if-necessary
	conan install --profile=$(conan_profile_name) $(tket_tests_recipe_dir) --install-folder=$(cmake_build_dir)/tket/tests
	conan install --profile=$(conan_profile_name) $(tket_proptests_recipe_dir) --install-folder=$(cmake_build_dir)/tket/proptests
	cmake -DCMAKE_BUILD_TYPE=$(build_type) -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S $(cmake_source_dir) -B $(cmake_build_dir)

##@ Build and Test

.PHONY: build-tests
build-tests: ## build test_tket binary
	cmake --build $(cmake_build_dir) --target test_tket -j $(n_cpus)

.PHONY: build-proptests
build-proptests: ## build proptests binary
	cmake --build $(cmake_build_dir) --target proptest -j $(n_cpus)

.PHONY: build
build: build-tests build-proptests ## build all

.PHONY: clean
clean: ## clean all
	cmake --build $(cmake_build_dir) --target clean -j $(n_cpus)

test_args="~[latex]"
.PHONY: test
test: build-tests ## build and run tket tests, override arguments to test binary using test_args variable, e.g., `make test_args='-r compact "[long]"' test`
	-$(cmake_build_dir)/tket/tests/bin/test_tket $(test_args)

.PHONY: test-file
File=""
file_test_filter=$(patsubst %.cpp,%, $(notdir $(File)))
test-file: build-tests ## build and run tket tests from a specific test file (usage: `make File=<test_file> test-file`)
	-$(cmake_build_dir)/tket/tests/bin/test_tket -# -r compact "[#$(file_test_filter)]"

.PHONY: proptests
proptests: build-proptests ## build and run tket proptests
	$(cmake_build_dir)/tket/proptests/bin/proptest

##@ Helpers

.PHONY: tket-version-update
new_tket_version=next
tket-version-update: ## update locally configured tket version (in all conan files) to next patch version or a semantic version specified using `new_tket_version=<semver>`
	./helper_scripts/update_tket_version.sh .. $(new_tket_version)


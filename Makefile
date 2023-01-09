# This makefile bundles commands for streamlined development
# Use `make help` for an overview of commands

conan_profile_name=tket
tket_remote_repository=https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
build_type=Debug
ifeq ($(build_type), Release)
  conan_build_dir=cmake-build-release
else
  conan_build_dir=cmake-build-debug
endif

n_cpus=$(shell python -c 'import multiprocessing as m; print(m.cpu_count())')

##@ General

.PHONY: help
help: ## Display this help
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_0-9-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Setup Conan Profile

.PHONY: conan-profile
conan-profile: ## Create and configure conan profile
	-conan profile new $(conan_profile_name) --detect
	-conan remote add tket-libs $(tket_remote_repository)
	conan config set general.revisions_enabled=1
	conan profile update options.tket:shared=True $(conan_profile_name)
	conan profile update options.tklog:shared=True $(conan_profile_name)
	conan profile update settings.build_type=$(build_type) $(conan_profile_name)

##@ Install Dependencies

.PHONY: pybind
pybind: ## install custom pybind package
	conan profile update settings.build_type=$(build_type) $(conan_profile_name)
	conan remove -f "pybind11/*"
	conan create --profile=$(conan_profile_name) recipes/pybind11

 # Currently the b2 package (a dependency of boost) in the remote tket-libs does not include compiler info
 # which can lead to errors when building boost from source. Conan will try to use
 # the remote binary even if it is incompatible with the compiler, leading to linking issues.
 # This fixes that problem by forcing a b2 build from source
.PHONY: b2
b2: # install and build b2 (for boost builds) from source
	conan profile update settings.build_type=$(build_type) $(conan_profile_name)
	conan install --profile=$(conan_profile_name) b2/4.9.2@ --build

.PHONY: local-libs
local-libs: # create conan packages for local libs from this repository
	conan profile update settings.build_type=$(build_type) $(conan_profile_name)
	conan create --profile=$(conan_profile_name) libs/tklog tket/stable --build=missing
	conan create --profile=$(conan_profile_name) libs/tkassert tket/stable --build=missing
	conan create --profile=$(conan_profile_name) libs/tkrng tket/stable --build=missing
	conan create --profile=$(conan_profile_name) libs/tktokenswap tket/stable --build=missing
	conan create --profile=$(conan_profile_name) libs/tkwsm tket/stable --build=missing

.PHONY: dev-env-from-scratch
dev-env-from-scratch: b2 local-libs ## Install tket, tket-tests, and tket-proptests into a local build directory and configure cmake (build everything from source)
	conan profile update settings.build_type=$(build_type) $(conan_profile_name)
	conan install --profile=$(conan_profile_name) recipes/tket  --install-folder=$(conan_build_dir)/tket/src --build=missing
	conan install --profile=$(conan_profile_name) recipes/tket-tests --install-folder=$(conan_build_dir)/tket/tests --build=missing
	conan install --profile=$(conan_profile_name) recipes/tket-proptests --install-folder=$(conan_build_dir)/tket/proptests --build=missing
	cmake -DCMAKE_BUILD_TYPE=$(build_type) -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S . -B $(conan_build_dir)

.PHONY: dev-env
dev-env: ## Install tket, tket-tests, and tket-proptests into a local build directory and configure cmake (dependencies from cache or remote)
	conan profile update settings.build_type=$(build_type) $(conan_profile_name)
	conan install --profile=$(conan_profile_name) recipes/tket  --install-folder=$(conan_build_dir)/tket/src
	conan create --profile=$(conan_profile_name) recipes/tket  tket/stable
	conan install --profile=$(conan_profile_name) recipes/tket-tests --install-folder=$(conan_build_dir)/tket/tests
	conan install --profile=$(conan_profile_name) recipes/tket-proptests --install-folder=$(conan_build_dir)/tket/proptests
	cmake -DCMAKE_BUILD_TYPE=$(build_type) -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S . -B $(conan_build_dir)

##@ Build

.PHONY: build
build: ## build with cmake
	cmake --build $(conan_build_dir) --target test_tket -j $(n_cpus)
	cmake --build $(conan_build_dir) --target proptest -j $(n_cpus)

##@ Test

test_args="~[latex]"
.PHONY: test
test: build ## run tket tests, override arguments to test binary using test_args variable, e.g., `make test_args='-r compact "[long]"' test`
	-$(conan_build_dir)/tket/tests/bin/test_tket $(test_args)

.PHONY: test-file
File=""
file_test_filter=$(patsubst %.cpp,%, $(notdir $(File)))
test-file: ## run tket tests from a specific test file (usage: `make File=<test_file> test-file`)
	-$(conan_build_dir)/tket/tests/bin/test_tket -# -r compact "[#$(file_test_filter)]"

.PHONY: proptests
proptests: ## run tket proptests
	$(conan_build_dir)/tket/proptests/bin/proptest

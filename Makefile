conan_profile=tket
conan_build_type=Release
conan_build_dir=cmake-build-debug

.PHONY: setup
 setup:
	-conan profile new $(conan_profile) --detect
	-conan remote add tket-libs https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs

.PHONY: configure
configure:
	conan config set general.revisions_enabled=1
	conan profile update options.tket:shared=True tket
	conan profile update options.tklog:shared=True tket
	conan profile update settings.build_type=$(conan_build_type) tket

.PHONY: pybind
pybind:
	conan remove -f "pybind11/*"
	conan create --profile=$(conan_profile) recipes/pybind11

 # currently the b2 package in the remote tket-libs does not include compiler info
 # which can lead to errors when building boost from source, since conan will try to use
 # the remote binary even if it is incompatible with the compiler
 # this fixes that
.PHONY: build-b2-from-source
build-b2-from-source:
	conan install --profile=$(conan_profile) b2/4.9.2@ --build

.PHONY: build-local-libs
build-local-libs:
	conan create --profile=$(conan_profile) libs/tklog tket/stable --build=missing
	conan create --profile=$(conan_profile) libs/tkassert tket/stable --build=missing
	conan create --profile=$(conan_profile) libs/tkrng tket/stable --build=missing
	conan create --profile=$(conan_profile) libs/tktokenswap tket/stable --build=missing
	conan create --profile=$(conan_profile) libs/tkwsm tket/stable --build=missing

.PHONY: install-tket
install-tket:
	conan install --profile=$(conan_profile) recipes/tket  --install-folder=$(conan_build_dir)/tket/src -build=missing
	conan install --profile=$(conan_profile) recipes/tket-tests --install-folder=$(conan_build_dir)/tket/tests --build=missing

.PHONY: configure-tket
configure-tket:
	cmake -DCMAKE_BUILD_TYPE=$(conan_build_type) -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S . -B $(conan_build_dir)

.PHONY: build-tket
build-tket:
	cmake --build $(conan_build_dir) --target test_tket -j 6

.PHONY: test
test:
	$(conan_build_dir)/tket/tests/bin/test_tket "~[latex]"

.PHONY: all
all: setup configure pybind build-b2-from-source build-local-libs install-tket configure-tket test
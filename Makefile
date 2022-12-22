
.PHONY: setup
 setup:
	-conan profile new tket --detect
	conan remote clean
	-conan remote add tket-libs https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs

.PHONY: configure
configure:
	conan config set general.revisions_enabled=1
	conan profile update options.tket:shared=True tket
	conan profile update options.tklog:shared=True tket
	conan profile update settings.build_type=Debug tket

.PHONY: pybind
pybind:
	conan remove -f "pybind11/*"
	conan create --profile=tket recipes/pybind11

 # currently the b2 package in the remote tket-libs does not include compiler info
 # which leads to errors when building boost from source, since conan will try to use
 # the remote binary even when it is incompatible with the compiler
 # this fixes that
.PHONY: build-b2-from-source
build-b2-from-source:
	conan install --profile=tket b2/4.9.2@ --build

.PHONY: build-others-from-source
build-others-from-source:
	conan install --profile=tket gmp/6.2.1@ --build
	conan install --profile=tket m4/1.4.19@ --build


.PHONY: build-local-libs
build-local-libs:
	conan create --profile=tket libs/tklog tket/stable --build=missing
	conan create --profile=tket libs/tkassert tket/stable --build=missing
	conan create --profile=tket libs/tkrng tket/stable --build=missing
	conan create --profile=tket libs/tktokenswap tket/stable --build=missing
	conan create --profile=tket libs/tkwsm tket/stable --build=missing

.PHONY: build-tket
build-tket:
	conan create --profile=tket recipes/tket tket/stable --build=missing

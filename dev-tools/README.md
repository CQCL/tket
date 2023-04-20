# Development Tools

*NB* The notes and scripts in this directory are obsolete since the migration
from conan 1 to conan 2. They should be either removed or updated.

This directory contains a [Makefile](Makefile) and helper scripts to bundle some common development tasks.
It currently supports building and running the projects tket-tests, tket-proptests. Building and testing
`pytket` is not (yet) supported.

## Requirements
- `cmake`
- `conan`
- `ccache`
- `ninja`
- `python`

## Usage
- Run `make` or `make help` within this directory to see a list of targets and their functionality
- To run the commands from the `tket` root directory, use `make -C dev-tools <target>`

## Build and Test Setup

Dependencies of `tket`, `tket-tests`, and `tket-proptests` are installed using `conan`. The `tket`, `tket-tests`,
and `tket-proptests` projects themselves are built using `cmake` directly through a top-level `CMakeLists.txt` file in
this directory. A cmake build directory (default: `cmake-build-debug`) is created within this directory for building
purposes.

### IDE Integration
If you are using an IDE that supports `cmake`, point it to this top-level `CMakeLists.txt` and build directory
in order to get functioning method lookup, linting, building, and testing. For some IDE's, more configuration may be necessary.

### Important configuration variables
Variables can be overriden using `make <variable_name>=<desired_value> <target>`

- `build_type`: `Debug` or `Release`, default: `Release` (Overriding with other values not intended) 
- `cmake_build_dir`: defaults to `cmake-build-debug` or `cmake-build-release` depending on `build_type`
- `conan_profile_name`: defaults to `tket` if `build_type=Debug`; `tket-release` if `build_type=Release`
- `tket_remote_repository`: conan remote for tket packages, default `https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs`
- `n_cpus`: number of threads to use for builds, defaults to value returned by `python -c 'import multiprocessing as m; print(m.cpu_count())'`

### Intended development process

1. Setup recommended conan profile: `make conan-profile` (or `make build_type=Release conan-profile` for release builds)
2. For linux systems: `make linux-profile-settings` is recommended
3. Install required conan packages and configure `cmake` build:
   1. `make dev-env-from-source` -> Build missing packages from source (e.g., for unsupported compilers, build settings) 
   2. `make dev-env` -> Dependencies from conan cache or remote (for supported conan configurations)
4. Build tket-tests and tket-proptest binaries: `make build`
5. Build and run tests or proptests: `make test` or `make proptests`

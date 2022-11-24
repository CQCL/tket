# Copyright 2019-2022 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from conans import ConanFile, CMake, tools
from conans.errors import ConanInvalidConfiguration
import os
import shutil


class TketConan(ConanFile):
    name = "tket"
    version = "1.0.31"
    license = "CQC Proprietary"
    homepage = "https://github.com/CQCL/tket"
    url = "https://github.com/conan-io/conan-center-index"
    description = "Quantum SDK"
    topics = ("quantum", "computation", "compiler")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "profile_coverage": [True, False],
    }
    default_options = {"shared": False, "profile_coverage": False}
    generators = "cmake"
    exports_sources = ["../../tket/src/*", "!*/build/*"]
    requires = (
        # libraries installed from remote:
        # https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
        "boost/1.80.0",
        "symengine/0.9.0",
        "eigen/3.4.0",
        "nlohmann_json/3.11.2",
        "tklog/0.1.2@tket/stable",
        "tkassert/0.1.1@tket/stable",
        "tkrng/0.1.2@tket/stable",
        "tktokenswap/0.1.2@tket/stable",
        "tkwsm/0.2.1@tket/stable",
    )

    # List of components in a topological sort according to dependencies:
    comps = [
        "Utils",
        "OpType",
        "Clifford",
        "Ops",
        "Graphs",
        "Gate",
        "PauliGraph",
        "Circuit",
        "Architecture",
        "Simulation",
        "Diagonalisation",
        "Characterisation",
        "ZX",
        "Converters",
        "Placement",
        "ArchAwareSynth",
        "Mapping",
        "MeasurementSetup",
        "Transformations",
        "Predicates",
    ]

    _cmake = None

    def _configure_cmake(self):
        if self._cmake is None:
            self._cmake = CMake(self)
            self._cmake.definitions["PROFILE_COVERAGE"] = self.options.profile_coverage
            self._cmake.configure()
        return self._cmake

    def validate(self):
        if self.options.profile_coverage and self.settings.compiler != "gcc":
            raise ConanInvalidConfiguration(
                "`profile_coverage` option only available with gcc"
            )

    def configure(self):
        # Disable features that are still under the LGPL.
        self.options["eigen"].MPL2_only = True

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy("LICENSE", dst="licenses", src="../..")
        for comp in self.comps:
            self.copy(f"{comp}/include/*.hpp", dst=f"include/{comp}", keep_path=False)
        self.copy("*.dll", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.dylib", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)

    def package_info(self):
        # These must be ordered "top to bottom", so that when building statically
        # unresolved symbols in libraries earlier in the list get resolved by later
        # libraries.
        self.cpp_info.libs = [f"tket-{comp}" for comp in reversed(self.comps)]

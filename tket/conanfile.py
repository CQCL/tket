# Copyright 2019-2023 Cambridge Quantum Computing
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

from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.errors import ConanInvalidConfiguration


class TketConan(ConanFile):
    name = "tket"
    version = "1.2.8"
    package_type = "library"
    license = "Apache 2"
    homepage = "https://github.com/CQCL/tket"
    url = "https://github.com/conan-io/conan-center-index"
    description = "Quantum SDK"
    topics = ("quantum", "computation", "compiler")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "profile_coverage": [True, False],
    }
    default_options = {"shared": False, "fPIC": True, "profile_coverage": False}
    exports_sources = "CMakeLists.txt", "src/*", "include/*"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        # Disable features that are still under the LGPL.
        self.options["eigen"].MPL2_only = True
        # Only need header libraries from boost.
        self.options["boost"].header_only = True

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.variables["PROFILE_COVERAGE"] = self.options.profile_coverage
        tc.generate()

    def validate(self):
        if self.options.profile_coverage and self.settings.compiler != "gcc":
            raise ConanInvalidConfiguration(
                "`profile_coverage` option only available with gcc"
            )

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["tket"]

    def requirements(self):
        # libraries installed from remote:
        # https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
        self.requires("boost/1.81.0", transitive_headers=True)
        self.requires("symengine/0.10.1", transitive_headers=True)
        self.requires("eigen/3.4.0", transitive_headers=True)
        self.requires("nlohmann_json/3.11.2", transitive_headers=True)
        self.requires("tklog/0.3.3@tket/stable")
        self.requires("tkassert/0.3.3@tket/stable", transitive_headers=True)
        self.requires("tkrng/0.3.3@tket/stable")
        self.requires("tktokenswap/0.3.3@tket/stable")
        self.requires("tkwsm/0.3.3@tket/stable")

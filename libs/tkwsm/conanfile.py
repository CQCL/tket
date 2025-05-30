# Copyright Quantinuum
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


class TkwsmConan(ConanFile):
    name = "tkwsm"
    version = "0.3.12"
    package_type = "library"
    license = "Apache 2"
    url = "https://github.com/CQCL/tket"
    description = "Weighted-subgraph-monomorphism algorithms library"
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "profile_coverage": [True, False],
    }
    default_options = {"shared": False, "fPIC": True, "profile_coverage": False}
    # a hash of the conanfile.py + the files listed here builds the conan revision
    exports_sources = "CMakeLists.txt", "cmake/*", "src/*", "include/*"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
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
        self.cpp_info.libs = ["tkwsm"]

    def requirements(self):
        self.requires("tkassert/0.3.4@tket/stable")
        self.requires("tkrng/0.3.3@tket/stable")
        self.requires(
            "boost/tci-1.88.0@tket/stable",
            transitive_headers=True,
            transitive_libs=False,
        )

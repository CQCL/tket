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

from conans import ConanFile, CMake
from conans.errors import ConanInvalidConfiguration


class TkassertConan(ConanFile):
    name = "tkassert"
    version = "0.1.1"
    license = "Apache 2"
    url = "https://github.com/CQCL/tket"
    description = "Assertions"
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "profile_coverage": [True, False],
    }
    default_options = {"shared": False, "fPIC": True, "profile_coverage": False}
    generators = "cmake"
    exports_sources = "src/*"
    requires = ["tklog/0.1.2@tket/stable"]

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def validate(self):
        if self.options.profile_coverage and self.settings.compiler != "gcc":
            raise ConanInvalidConfiguration(
                "`profile_coverage` option only available with gcc"
            )

    def build(self):
        cmake = CMake(self)
        cmake.definitions["PROFILE_COVERAGE"] = self.options.profile_coverage
        cmake.configure(source_folder="src")
        cmake.build()

    def package(self):
        self.copy("include/*.hpp", dst="include/tkassert", src="src", keep_path=False)
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.dll", dst="lib", keep_path=False)
        self.copy("*.dylib*", dst="lib", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["tkassert"]

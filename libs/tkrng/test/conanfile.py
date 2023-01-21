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

from conans import ConanFile, CMake
import platform


class TestTkrngConan(ConanFile):
    name = "test-tkrng"
    version = "0.1.2"
    license = "Apache 2"
    url = "https://github.com/CQCL/tket"
    description = "Unit tests for tkrng"
    settings = "os", "compiler", "build_type", "arch"
    options = {"with_coverage": [True, False]}
    default_options = {"with_coverage": False}
    generators = "cmake"
    exports_sources = "*"
    requires = ["tkrng/0.1.2", "catch2/3.3.0"]

    _cmake = None

    def _configure_cmake(self):
        if self._cmake is None:
            self._cmake = CMake(self)
            self._cmake.definitions["WITH_COVERAGE"] = self.options.with_coverage
            self._cmake.configure()
        return self._cmake

    def validate(self):
        if self.options.with_coverage and self.settings.compiler != "gcc":
            raise ConanInvalidConfiguration(
                "`with_coverage` option only available with gcc"
            )

    def configure(self):
        if self.options.with_coverage:
            self.options["tkrng"].profile_coverage = True

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        executable_filename = "bin/test-tkrng"
        if platform.system() == "Windows":
            executable_filename = executable_filename + ".exe"
        self.copy(executable_filename)

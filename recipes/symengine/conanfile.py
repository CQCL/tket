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
import os

required_conan_version = ">=1.33.0"


class SymengineConan(ConanFile):
    name = "symengine"
    version = "0.8.1.1"
    description = "A fast symbolic manipulation library, written in C++"
    license = "MIT"
    topics = ("symbolic", "algebra")
    homepage = "https://symengine.org/"
    exports_sources = ["CMakeLists.txt"]
    generators = "cmake"
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "integer_class": ["boostmp", "gmp"],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "integer_class": "boostmp",
    }
    short_paths = True

    _cmake = None

    def requirements(self):
        if self.options.integer_class == "boostmp":
            self.requires("boost/1.78.0")
        else:
            self.requires("gmp/6.2.1")

    def source(self):
        git = tools.Git()
        git.clone(
            "https://github.com/symengine/symengine",
            branch="2eb109a0e554b62683662cc5559fccf2ea0c0348",
            shallow=True,
        )

    def _configure_cmake(self):
        if self._cmake is None:
            self._cmake = CMake(self)
            self._cmake.definitions["BUILD_TESTS"] = False
            self._cmake.definitions["BUILD_BENCHMARKS"] = False
            self._cmake.definitions["INTEGER_CLASS"] = self.options.integer_class
            self._cmake.definitions["MSVC_USE_MT"] = False
            self._cmake.configure()
        return self._cmake

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        if self.options.shared:
            del self.options.fPIC

    def build(self):
        tools.replace_in_file(
            os.path.join(self.source_folder, "CMakeLists.txt"),
            "project(symengine)",
            """project(symengine)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()""",
        )
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy("LICENSE", dst="licenses", src=self.source_folder)
        cmake = self._configure_cmake()
        cmake.install()
        cmake.patch_config_paths()
        # [CMAKE-MODULES-CONFIG-FILES (KB-H016)]
        tools.remove_files_by_mask(self.package_folder, "*.cmake")
        # [DEFAULT PACKAGE LAYOUT (KB-H013)]
        tools.rmdir(os.path.join(self.package_folder, "CMake"))

    def package_info(self):
        self.cpp_info.libs = ["symengine"]
        if any("teuchos" in v for v in tools.collect_libs(self)):
            self.cpp_info.libs.append("teuchos")
        self.cpp_info.names["cmake_find_package"] = "symengine"
        # FIXME: symengine exports a non-namespaced `symengine` target.
        self.cpp_info.names["cmake_find_package_multi"] = "symengine"

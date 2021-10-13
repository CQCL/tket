# Copyright 2019-2021 Cambridge Quantum Computing
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
    version = "1.0.1"
    license = "CQC Proprietary"
    author = "Alec Edgington <alec.edgington@cambridgequantum.com>"
    url = "https://github.com/CQCL-DEV/tket"
    description = "Quantum SDK"
    topics = ("quantum", "computation", "compiler")
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True], "profile_coverage": [True, False]}
    default_options = {"shared": True, "profile_coverage": False}
    generators = "cmake"
    # Putting "patches" in both "exports_sources" and "exports" means that this works
    # in either the CI workflow (`conan create`) or the development workflow
    # (`conan build`). Maybe there is a better way?
    exports_sources = ["../../tket/src/*", "!*/build/*", "patches/*"]
    exports = ["patches/*"]
    requires = (
        "boost/1.77.0",
        "symengine/0.8.1",
        "eigen/3.4.0",
        "spdlog/1.9.2",
        "nlohmann_json/3.10.3",
    )

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
        # Use header-only version of spdlog/fmt to avoid linker warnings.
        self.options["fmt"].header_only = True
        self.options["spdlog"].header_only = True

    def build(self):
        # Build with boost patches
        boost_include_path = self.deps_cpp_info["boost"].include_paths[0]
        curdir = os.path.dirname(os.path.realpath(__file__))
        patches = {
            # TKET-1407
            # If and when the boost package is fixed we will remove this.
            os.path.join(
                boost_include_path, "boost", "graph", "detail", "adjacency_list.hpp"
            ): os.path.join(curdir, "patches", "adjacency_list.diff"),
            # TKET-1376
            # This will be submitted as a PR to boost. If it is accepted we will remove
            # the patch here. If not, we should consider switching to another library
            # for subgraph matching, or writing our own code.
            os.path.join(
                boost_include_path, "boost", "graph", "vf2_sub_graph_iso.hpp"
            ): os.path.join(curdir, "patches", "vf2_sub_graph_iso.diff"),
        }
        for filepath, patch_file in patches.items():
            print("Patching " + filepath)
            shutil.copyfile(filepath, filepath + ".original")
            tools.patch(base_path=boost_include_path, patch_file=patch_file)
        cmake = self._configure_cmake()
        try:
            print("Building")
            cmake.build()
        finally:
            for filepath in patches.keys():
                print("Restoring " + filepath)
                shutil.move(filepath + ".original", filepath)

    def package(self):
        self.copy("*.hpp", dst="include")
        self.copy("*.dll", dst="lib", keep_path=False)
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.dylib", dst="lib", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["tket"]

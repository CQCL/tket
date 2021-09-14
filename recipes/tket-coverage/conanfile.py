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
import os
import shutil


class TketCoverageConan(ConanFile):
    name = "tket-coverage"
    version = "0.6.3"
    license = "CQC Proprietary"
    author = "Melf Johannsen <melf.johannsen@cambridgequantum.com>"
    url = "https://github.com/CQCL-DEV/tket"
    description = "Unit tests for tket"
    topics = ("quantum", "computation", "compiler")
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    exports_sources = ["../../bubble/*", "../tket/patches/*"]
    exports = ["../tket/patches/*"]
    symlinks = True
    requires = (
        "boost/1.77.0",
        "symengine/0.8.1",
        "eigen/3.4.0",
        "spdlog/1.9.2",
        "catch2/2.13.7",
        "nlohmann_json/3.10.2",
    )

    def build(self):
        self.run("ln -s ../src coverage/src")
        self.run("ln -s ../tests coverage/tests")
        cmake = CMake(self)
        cmake.configure(source_folder="coverage")

        # Build with boost patches
        boost_include_path = self.deps_cpp_info["boost"].include_paths[0]
        curdir = os.path.dirname(os.path.realpath(__file__))
        patches = {
            # TKET-1407
            # If and when the boost package is fixed we will remove this.
            os.path.join(
                boost_include_path, "boost", "graph", "detail", "adjacency_list.hpp"
            ): os.path.join(curdir, "adjacency_list.diff"),
            # TKET-1376
            # This will be submitted as a PR to boost. If it is accepted we will remove
            # the patch here. If not, we should consider switching to another library
            # for subgraph matching, or writing our own code.
            os.path.join(
                boost_include_path, "boost", "graph", "vf2_sub_graph_iso.hpp"
            ): os.path.join(curdir, "vf2_sub_graph_iso.diff"),
        }
        for filepath, patch_file in patches.items():
            print("Patching " + filepath)
            shutil.copyfile(filepath, filepath + ".original")
            tools.patch(base_path=boost_include_path, patch_file=patch_file)
        try:
            print("Building")
            self.run("cmake . -DCMAKE_BUILD_TYPE=Debug -Wno-dev")
            self.run("cmake --build . -- coverage")
            self.run("rm coverage/src")
            self.run("rm coverage/tests")
        finally:
            for filepath in patches.keys():
                print("Restoring " + filepath)
                shutil.move(filepath + ".original", filepath)

    def package(self):
        self.copy("bin/test_tket")
        self.copy("random_angles.txt", dst="bin", keep_path=False)

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
import platform


class TketProptestsConan(ConanFile):
    name = "tket-proptests"
    version = "0.6.3"
    license = "CQC Proprietary"
    author = "Alec Edgington <alec.edgington@cambridgequantum.com>"
    url = "https://github.com/CQCL/tket"
    description = "Property tests for tket"
    topics = ("quantum", "computation", "compiler")
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    exports_sources = "../../tket/proptests/*"
    requires = (
        "tket/1.0.1",
        "rapidcheck/cci.20220514",
    )

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        executable_filename = "bin/proptest"
        if platform.system() == "Windows":
            executable_filename = executable_filename + ".exe"
        self.copy(executable_filename)

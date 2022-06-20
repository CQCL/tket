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


class TktokenswapConan(ConanFile):
    name = "tktokenswap"
    version = "0.1.0"
    license = "Apache 2"
    url = "https://github.com/CQCL/tket"
    description = "Token swapping algorithms library"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}
    generators = "cmake"
    exports_sources = "src/*"
    requires = [
        "tklog/0.1.1@tket/stable",
        "tkassert/0.1.0@tket/stable",
        "tkrng/0.1.1@tket/stable",
        "boost/1.79.0",
    ]

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def build(self):
        cmake = CMake(self)
        cmake.configure(source_folder="src")
        cmake.build()

    def package(self):
        self.copy(
            "include/*.hpp", dst="include/tktokenswap", src="src", keep_path=False
        )
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.dll", dst="lib", keep_path=False)
        self.copy("*.dylib*", dst="lib", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["tktokenswap"]

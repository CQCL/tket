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
from conan.tools.cmake import CMake, CMakeToolchain
from conan.tools.layout import basic_layout
from conan.tools.files import get, copy, replace_in_file, rm, rmdir
import os


class PyBind11Conan(ConanFile):
    name = "pybind11"
    version = "2.11.1"
    description = "Seamless operability between C++11 and Python"
    topics = "pybind11", "python", "binding"
    homepage = "https://github.com/pybind/pybind11"
    license = "BSD-3-Clause"
    settings = "os", "arch", "compiler", "build_type"
    no_copy_source = True

    def layout(self):
        basic_layout(self, src_folder="src")

    def source(self):
        get(
            self,
            f"https://github.com/pybind/pybind11/archive/refs/tags/v{self.version}.tar.gz",
            destination=self.source_folder,
            strip_root=True,
        )

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["PYBIND11_INSTALL"] = True
        tc.variables["PYBIND11_TEST"] = False
        tc.variables["PYBIND11_CMAKECONFIG_INSTALL_DIR"] = "lib/cmake/pybind11"
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )
        cmake = CMake(self)
        cmake.install()
        for filename in [
            "pybind11Targets.cmake",
            "pybind11Config.cmake",
            "pybind11ConfigVersion.cmake",
        ]:
            rm(
                self,
                filename,
                os.path.join(self.package_folder, "lib", "cmake", "pybind11"),
            )

        rmdir(self, os.path.join(self.package_folder, "share"))

        replace_in_file(
            self,
            os.path.join(
                self.package_folder, "lib", "cmake", "pybind11", "pybind11Common.cmake"
            ),
            "if(TARGET pybind11::lto)",
            "if(FALSE)",
        )
        replace_in_file(
            self,
            os.path.join(
                self.package_folder, "lib", "cmake", "pybind11", "pybind11Common.cmake"
            ),
            "add_library(",
            "# add_library(",
        )
        replace_in_file(
            self,
            os.path.join(
                self.package_folder, "lib", "cmake", "pybind11", "pybind11Common.cmake"
            ),
            """add_library(pybind11::embed IMPORTED INTERFACE ${optional_global})
set_property(
  TARGET pybind11::embed
  APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES pybind11::pybind11)""",
            """# add_library(pybind11::embed IMPORTED INTERFACE ${optional_global})
# set_property(
#   TARGET pybind11::embed
#   APPEND
#   PROPERTY INTERFACE_LINK_LIBRARIES pybind11::pybind11)""",
        )
        replace_in_file(
            self,
            os.path.join(
                self.package_folder, "lib", "cmake", "pybind11", "pybind11Tools.cmake"
            ),
            """  target_link_libraries(pybind11::embed INTERFACE pybind11::pybind11
                                                  pybind11::_ClassicPythonLibraries)""",
            """  # target_link_libraries(pybind11::embed INTERFACE pybind11::pybind11
                 #                                  pybind11::_ClassicPythonLibraries)""",
        )

    def package_id(self):
        self.info.clear()

    def package_info(self):
        cmake_base_path = os.path.join("lib", "cmake", "pybind11")
        self.cpp_info.set_property("cmake_target_name", "pybind11_all_do_not_use")
        self.cpp_info.components["headers"].includedirs = ["include"]
        self.cpp_info.components["pybind11_"].set_property(
            "cmake_target_name", "pybind11::pybind11"
        )
        self.cpp_info.components["pybind11_"].set_property(
            "cmake_module_file_name", "pybind11"
        )
        self.cpp_info.components["pybind11_"].names["cmake_find_package"] = "pybind11"
        self.cpp_info.components["pybind11_"].builddirs = [cmake_base_path]
        self.cpp_info.components["pybind11_"].requires = ["headers"]
        cmake_file = os.path.join(cmake_base_path, "pybind11Common.cmake")
        self.cpp_info.set_property("cmake_build_modules", [cmake_file])
        for generator in ["cmake_find_package", "cmake_find_package_multi"]:
            self.cpp_info.components["pybind11_"].build_modules[generator].append(
                cmake_file
            )
        # self.cpp_info.components["embed"].requires = ["pybind11_"]
        self.cpp_info.components["module"].requires = ["pybind11_"]
        self.cpp_info.components["python_link_helper"].requires = ["pybind11_"]
        self.cpp_info.components["windows_extras"].requires = ["pybind11_"]
        self.cpp_info.components["lto"].requires = ["pybind11_"]
        self.cpp_info.components["thin_lto"].requires = ["pybind11_"]
        self.cpp_info.components["opt_size"].requires = ["pybind11_"]
        self.cpp_info.components["python2_no_register"].requires = ["pybind11_"]

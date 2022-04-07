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

from conans import CMake, ConanFile, tools
from conans.errors import ConanInvalidConfiguration
import os

required_conan_version = ">=1.43.0"


class SpdlogConan(ConanFile):
    name = "spdlog"
    version = "1.10.0"
    description = "Fast C++ logging library"
    url = "https://github.com/conan-io/conan-center-index"
    homepage = "https://github.com/gabime/spdlog"
    topics = ("spdlog", "logging", "header-only")
    license = "MIT"

    settings = "os", "arch", "compiler", "build_type"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "header_only": [True, False],
        "wchar_support": [True, False],
        "wchar_filenames": [True, False],
        "no_exceptions": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "header_only": False,
        "wchar_support": False,
        "wchar_filenames": False,
        "no_exceptions": False,
    }

    exports_sources = "CMakeLists.txt"
    generators = "cmake", "cmake_find_package", "cmake_find_package_multi"
    _cmake = None

    @property
    def _source_subfolder(self):
        return "source_subfolder"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        if self.options.shared:
            del self.options.fPIC
        if self.options.header_only:
            del self.options.shared
            del self.options.fPIC

    def requirements(self):
        pass

    def validate(self):
        if self.settings.os != "Windows" and (
            self.options.wchar_support or self.options.wchar_filenames
        ):
            raise ConanInvalidConfiguration("wchar is only supported under windows")
        if self.options.get_safe("shared", False):
            if self.settings.os == "Windows" and tools.Version(self.version) < "1.6.0":
                raise ConanInvalidConfiguration(
                    "spdlog shared lib is not yet supported under windows"
                )
            if (
                self.settings.compiler == "Visual Studio"
                and "MT" in self.settings.compiler.runtime
            ):
                raise ConanInvalidConfiguration(
                    "Visual Studio build for shared library with MT runtime is not supported"
                )

    def package_id(self):
        if self.options.header_only:
            self.info.header_only()

    def source(self):
        tools.get(
            **self.conan_data["sources"][self.version],
            destination=self._source_subfolder,
            strip_root=True
        )

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake
        self._cmake = CMake(self)
        self._cmake.definitions["SPDLOG_BUILD_EXAMPLE"] = False
        self._cmake.definitions["SPDLOG_BUILD_EXAMPLE_HO"] = False
        self._cmake.definitions["SPDLOG_BUILD_TESTS"] = False
        self._cmake.definitions["SPDLOG_BUILD_TESTS_HO"] = False
        self._cmake.definitions["SPDLOG_BUILD_BENCH"] = False
        self._cmake.definitions["SPDLOG_BUILD_SHARED"] = (
            not self.options.header_only and self.options.shared
        )
        self._cmake.definitions["SPDLOG_WCHAR_SUPPORT"] = self.options.wchar_support
        self._cmake.definitions["SPDLOG_WCHAR_FILENAMES"] = self.options.wchar_filenames
        self._cmake.definitions["SPDLOG_INSTALL"] = True
        self._cmake.definitions["SPDLOG_NO_EXCEPTIONS"] = self.options.no_exceptions
        if self.settings.os in ("iOS", "tvOS", "watchOS"):
            self._cmake.definitions["SPDLOG_NO_TLS"] = True
        if tools.cross_building(self):
            # Workaround to find CMake config files in some cross-build scenario
            # TODO: to remove if something like https://github.com/conan-io/conan/issues/9427#issuecomment-993685376
            # is implemented
            self._cmake.definitions["CMAKE_FIND_ROOT_PATH_MODE_PACKAGE"] = "BOTH"
        self._cmake.configure()
        return self._cmake

    def _disable_werror(self):
        tools.replace_in_file(
            os.path.join(self._source_subfolder, "cmake", "utils.cmake"), "/WX", ""
        )

    def build(self):
        self._disable_werror()
        if not self.options.header_only:
            cmake = self._configure_cmake()
            cmake.build()

    def package(self):
        self.copy("LICENSE", dst="licenses", src=self._source_subfolder)
        if self.options.header_only:
            self.copy(
                pattern="*.h",
                dst="include",
                src=os.path.join(self._source_subfolder, "include"),
            )
        else:
            cmake = self._configure_cmake()
            cmake.install()
            tools.rmdir(os.path.join(self.package_folder, "lib", "cmake"))
            tools.rmdir(os.path.join(self.package_folder, "lib", "pkgconfig"))
            tools.rmdir(os.path.join(self.package_folder, "lib", "spdlog", "cmake"))

    def package_info(self):
        target = "spdlog_header_only" if self.options.header_only else "spdlog"
        self.cpp_info.set_property("cmake_file_name", "spdlog")
        self.cpp_info.set_property("cmake_target_name", "spdlog::{}".format(target))
        self.cpp_info.set_property("pkg_config_name", "spdlog")

        self.cpp_info.names["cmake_find_package"] = "spdlog"
        self.cpp_info.names["cmake_find_package_multi"] = "spdlog"
        self.cpp_info.components["libspdlog"].names["cmake_find_package"] = target
        self.cpp_info.components["libspdlog"].names["cmake_find_package_multi"] = target

        if not self.options.header_only:
            self.cpp_info.components["libspdlog"].libs = tools.collect_libs(self)
            self.cpp_info.components["libspdlog"].defines.append("SPDLOG_COMPILED_LIB")
        if self.options.wchar_support:
            self.cpp_info.components["libspdlog"].defines.append(
                "SPDLOG_WCHAR_TO_UTF8_SUPPORT"
            )
        if self.options.wchar_filenames:
            self.cpp_info.components["libspdlog"].defines.append(
                "SPDLOG_WCHAR_FILENAMES"
            )
        if self.options.no_exceptions:
            self.cpp_info.components["libspdlog"].defines.append("SPDLOG_NO_EXCEPTIONS")
        if self.settings.os in ["Linux", "FreeBSD"]:
            self.cpp_info.components["libspdlog"].system_libs = ["pthread"]
        if self.options.header_only and self.settings.os in ("iOS", "tvOS", "watchOS"):
            self.cpp_info.components["libspdlog"].defines.append("SPDLOG_NO_TLS")

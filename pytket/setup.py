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

import os
import platform
import re
import subprocess
import sys
import json
import shutil
from multiprocessing import cpu_count
from distutils.version import LooseVersion
from concurrent.futures import ThreadPoolExecutor as Pool
from shutil import which
import setuptools  # type: ignore
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext  # type: ignore
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


def libfile(name):
    sysname = platform.system()
    if sysname == "Linux":
        return "lib" + name + ".so"
    elif sysname == "Darwin":
        return "lib" + name + ".dylib"
    elif sysname == "Windows":
        return name + ".dll"
    else:
        return None


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        self.cfg = "Release"
        self.check_extensions_list(self.extensions)
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(self.extensions[0].name))
        )
        extsource = self.extensions[0].sourcedir
        self.cmake_config(extdir, extsource)

        if platform.system() == "Windows":
            for ext in self.extensions:
                self.build_extension(ext)
        else:
            num_jobs = int(os.environ.get("MAKE_N_THREADS", default=cpu_count()))
            with Pool(num_jobs) as pool:
                _future = list(pool.map(self.build_extension, self.extensions))

        if platform.system() in ["Darwin", "Windows"]:
            # Hack to put the tket library alongside the extension libraries
            conan_tket_profile = os.getenv("CONAN_TKET_PROFILE", default="tket")
            jsondump = "conaninfo.json"
            subprocess.run(
                [
                    "conan",
                    "info",
                    "--profile",
                    conan_tket_profile,
                    "--path",
                    "--json",
                    jsondump,
                    ".",
                ],
                cwd=extsource,
            )
            with open(jsondump) as f:
                conaninfo = dict([(comp["reference"], comp) for comp in json.load(f)])
            os.remove(jsondump)
            reqs = conaninfo["conanfile.txt"]["requires"]
            tket_reqs = [req for req in reqs if req.startswith("tket/")]
            assert len(tket_reqs) == 1
            tket_req = tket_reqs[0]
            directory = conaninfo[tket_req]["package_folder"]
            tket_libs = [
                "tket-Utils",
                "tket-ZX",
                "tket-OpType",
                "tket-Clifford",
                "tket-Ops",
                "tket-Graphs",
                "tket-Gate",
                "tket-PauliGraph",
                "tket-Circuit",
                "tket-Architecture",
                "tket-Simulation",
                "tket-Diagonalisation",
                "tket-Characterisation",
                "tket-Converters",
                "tket-TokenSwapping",
                "tket-Placement",
                "tket-Mapping",
                "tket-MeasurementSetup",
                "tket-Transformations",
                "tket-ArchAwareSynth",
                "tket-Predicates",
            ]
            for tket_lib in tket_libs:
                shutil.copy(os.path.join(directory, "lib", libfile(tket_lib)), extdir)

    def cmake_config(self, extdir, extsource):

        env = os.environ.copy()

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        tket_lib_dir = os.environ.get("TKET_LIB_FILES", default=None)
        if tket_lib_dir:
            cmake_args.append("-DTKET_LIB_FILES=" + tket_lib_dir)

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(
                    self.cfg.upper(), extdir
                )
            ]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + self.cfg]
            if which("ninja") is not None:
                cmake_args += ["-G", "Ninja"]

        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        conan_cmd = os.getenv("CONAN_CMD", default="conan")
        conan_tket_profile = os.getenv("CONAN_TKET_PROFILE", default="tket")
        install_cmd = [
            conan_cmd,
            "install",
            "--profile=" + conan_tket_profile,
            "--build=missing",
            "-o",
            "tket:shared=True",
            extsource,
        ]
        if platform.system() == "Darwin" and platform.processor() == "arm":
            install_cmd.extend(
                [
                    "-o",
                    "boost:without_fiber=True",
                    "-o",
                    "boost:without_json=True",
                    "-o",
                    "boost:without_nowide=True",
                ]
            )
        subprocess.check_call(install_cmd, cwd=self.build_temp, env=env)

        subprocess.check_call(
            ["cmake", extsource] + cmake_args, cwd=self.build_temp, env=env
        )

    def build_extension(self, ext):
        build_args = ["--config", self.cfg]
        if platform.system() == "Windows":
            build_args += ["--", "/m"]
        subprocess.check_call(
            ["cmake", "--build", ".", "--target", ext.name.split(".")[2]] + build_args,
            cwd=self.build_temp,
        )


binders = [
    "logging",
    "utils_serialization",
    "circuit",
    "passes",
    "predicates",
    "partition",
    "pauli",
    "mapping",
    "transform",
    "tailoring",
    "tableau",
    "zx",
    "placement",
    "architecture",
]

setup_dir = os.path.abspath(os.path.dirname(__file__))
plat_name = os.getenv("WHEEL_PLAT_NAME")


class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        if plat_name is not None:
            print(f"Overriding plat_name to {plat_name}")
            self.plat_name = plat_name
            self.plat_name_supplied = True


setup(
    name="pytket",
    author="TKET development team",
    author_email="tket-support@cambridgequantum.com",
    python_requires=">=3.8",
    url="https://cqcl.github.io/tket/pytket/api/",
    description="Python module for interfacing with the CQC tket library of quantum "
    "software",
    long_description=open("package.md", "r").read(),
    long_description_content_type="text/markdown",
    license="Apache 2",
    packages=setuptools.find_packages(),
    install_requires=[
        "sympy ~=1.6",
        "numpy >=1.21.4, <2.0",
        "lark-parser ~=0.7",
        "scipy >=1.7.2, <2.0",
        "networkx ~= 2.4",
        "graphviz ~= 0.14",
        "jinja2 ~= 3.0",
        "types-pkg_resources",
        "typing-extensions ~= 4.2",
    ],
    ext_modules=[
        CMakeExtension("pytket._tket.{}".format(binder)) for binder in binders
    ],
    cmdclass={"build_ext": CMakeBuild, "bdist_wheel": bdist_wheel},
    classifiers=[
        "Environment :: Console",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    include_package_data=True,
    package_data={"pytket": ["py.typed"]},
    zip_safe=False,
    use_scm_version={
        "root": os.path.dirname(setup_dir),
        "write_to": os.path.join(setup_dir, "pytket", "_version.py"),
        "write_to_template": "__version__ = '{version}'",
    },
)

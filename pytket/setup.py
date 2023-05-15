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

import os
import platform
import re
import subprocess
import sys
import json
import shutil
from distutils.version import LooseVersion
import setuptools  # type: ignore
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext  # type: ignore
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class ConanBuild(build_ext):
    def run(self):
        self.check_extensions_list(self.extensions)
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(self.extensions[0].name))
        )
        extsource = self.extensions[0].sourcedir

        jsonstr = subprocess.check_output(
            [
                "conan",
                "create",
                ".",
                "--build=missing",
                "-o",
                "boost/*:header_only=True",
                "-o",
                "tket/*:shared=True",
                "-o",
                "tklog/*:shared=True",
                "--format",
                "json",
            ],
            cwd=extsource,
        )

        # Collect the paths to the libraries to package together
        conaninfo = json.loads(jsonstr)
        nodes = conaninfo["graph"]["nodes"]
        if os.path.exists(extdir):
            shutil.rmtree(extdir)
        os.makedirs(extdir)
        for comp in ["tklog", "tket", "pytket"]:
            compnodes = [node for node in nodes if node["ref"].startswith(comp + "/")]
            assert len(compnodes) == 1
            compnode = compnodes[0]
            lib_folder = os.path.join(compnode["package_folder"], "lib")
            for lib in os.listdir(lib_folder):
                libpath = os.path.join(lib_folder, lib)
                # Don't copy the `cmake` directory.
                if not os.path.isdir(libpath):
                    shutil.copy(libpath, extdir)


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
    python_requires=">=3.9",
    project_urls={
        "Documentation": "https://cqcl.github.io/tket/pytket/api/index.html",
        "Source": "https://github.com/CQCL/tket",
        "Tracker": "https://github.com/CQCL/tket/issues",
    },
    description="Python module for interfacing with the CQC tket library of quantum "
    "software",
    long_description=open("package.md", "r").read(),
    long_description_content_type="text/markdown",
    license="Apache 2",
    packages=setuptools.find_packages() + ["pytket.qasm.includes"],
    install_requires=[
        "sympy ~=1.6",
        "numpy >=1.21.4, <2.0",
        "lark-parser ~=0.7",
        "scipy >=1.7.2, <2.0",
        "networkx ~= 3.0",
        "graphviz ~= 0.14",
        "jinja2 ~= 3.0",
        "types-pkg_resources",
        "typing-extensions ~= 4.2",
        "qwasm ~= 1.0",
    ],
    extras_require={
        "ZX": ["quimb ~= 1.5"],
    },
    ext_modules=[
        CMakeExtension("pytket._tket.{}".format(binder)) for binder in binders
    ],
    cmdclass={"build_ext": ConanBuild, "bdist_wheel": bdist_wheel},
    classifiers=[
        "Environment :: Console",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
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

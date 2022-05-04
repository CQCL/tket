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

from os.path import exists
import base64
import hashlib

from typing import List
from pytket import Circuit


class WasmFileHandler:
    """Add a wasm file to your workflow, stores a copy of the file and
    checks the function signatures of the file. Offers function to add
    a wasm op to a circuit"""

    def __init__(self, filepath: str, checksignatures: bool = False):
        """construct a wasm file handler"""
        self._filepath = filepath
        self._checksignatures = checksignatures

        if not exists(self._filepath):
            raise ValueError("wasm file not found at given path")

        with open(self._filepath, "rb") as wasm_file:
            self._wasm_file_encoded = base64.b64encode(wasm_file.read())

        self._wasmuid = hashlib.md5(self._wasm_file_encoded).hexdigest()

    def __str__(self) -> str:
        """str representation of the wasm file"""
        return self._wasmuid

    def add_wasmop_to_circuit(
        self,
        circ: Circuit,
        funcname: str,
        i32list_i: List[int],
        i32list_o: List[int],
        args: List[int],
    ) -> None:
        """function to add a wasm op to the circuit with the given name of the
        function and the classical bits assigend to the input of the function"""

        # add check of function signatures here

        circ.add_wasm(funcname, self._wasmuid, i32list_i, i32list_o, args)

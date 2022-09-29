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
from wasm import (  # type: ignore
    decode_module,
    SEC_TYPE,
    SEC_FUNCTION,
    SEC_EXPORT,
    LANG_TYPE_I32,
    LANG_TYPE_I64,
    LANG_TYPE_F32,
    LANG_TYPE_F64,
    LANG_TYPE_EMPTY,
)


class WasmFileHandler:
    """Add a wasm file to your workflow, stores a copy of the file and
    checks the function signatures of the file. Offers function to add
    a wasm op to a circuit"""

    type_lookup = {
        LANG_TYPE_I32: "i32",
        LANG_TYPE_I64: "i64",
        LANG_TYPE_F32: "f32",
        LANG_TYPE_F64: "f64",
        LANG_TYPE_EMPTY: None,
    }

    def __init__(self, filepath: str):
        """construct a wasm file handler"""
        self._filepath = filepath

        if not exists(self._filepath):
            raise ValueError("wasm file not found at given path")

        with open(self._filepath, "rb") as wasm_file:
            self._wasm_file_encoded = base64.b64encode(wasm_file.read())

        self._wasmuid = hashlib.md5(self._wasm_file_encoded).hexdigest()

        self._function_signatures: list = []
        self._function_names: list = []

        # stores the names of the functions mapped
        #  to the number of parameters and the number of return values
        self._functions = dict()

        # contains the list of functions that are not allowed
        # to use in pytket (because of types that are not i32)
        self._unvalid_functions = []

        self._wasm_file = base64.decodebytes(self._wasm_file_encoded)

        mod_iter = iter(decode_module(self._wasm_file))
        _, _ = next(mod_iter)

        for _, cur_sec_data in mod_iter:
            # read in list of function signatures
            if cur_sec_data.id == SEC_TYPE:
                for idx, entry in enumerate(cur_sec_data.payload.entries):
                    self._function_signatures.append({})
                    self._function_signatures[idx]["parameter_types"] = [
                        self.type_lookup[pt] for pt in entry.param_types
                    ]
                    if entry.return_count > 1:
                        if (
                            isinstance(entry.return_type, list)
                            and len(entry.return_type) == entry.return_count
                        ):
                            self._function_signatures[idx]["return_types"] = [
                                self.type_lookup[rt] for rt in entry.return_type
                            ]
                        elif isinstance(entry.return_type, int):
                            self._function_signatures[idx]["return_types"] = [
                                self.type_lookup[entry.return_type]
                            ] * entry.return_count
                        else:
                            raise ValueError(
                                f"Only parameter and return values of i32 types are"
                                + f"allowed, found different type: {entry.return_type}"
                            )
                    elif entry.return_count == 1:
                        self._function_signatures[idx]["return_types"] = [
                            self.type_lookup[entry.return_type]
                        ]
                    else:
                        self._function_signatures[idx]["return_types"] = []

            # read in list of function names
            elif cur_sec_data.id == SEC_EXPORT:
                for entry in cur_sec_data.payload.entries:
                    if entry.kind == 0:
                        self._function_names.append(entry.field_str.tobytes().decode())

            # read in map of function signatures to function names
            elif cur_sec_data.id == SEC_FUNCTION:
                self._function_types = cur_sec_data.payload.types

        for i, x in enumerate(self._function_names):

            # check for only i32 type in parameters and return values
            valid = True
            for t in self._function_signatures[self._function_types[i]][
                "parameter_types"
            ]:
                if t != "i32":
                    valid = False
            for t in self._function_signatures[self._function_types[i]]["return_types"]:
                if t != "i32":
                    valid = False

            if valid:
                self._functions[x] = (
                    len(
                        self._function_signatures[self._function_types[i]][
                            "parameter_types"
                        ]
                    ),
                    len(
                        self._function_signatures[self._function_types[i]][
                            "return_types"
                        ]
                    ),
                )

            if not valid:
                self._unvalid_functions.append(x)

    def __str__(self) -> str:
        """str representation of the wasm file"""
        return self._wasmuid

    def __repr__(self) -> str:
        """str representation of the containment of the wasm file"""
        result = f"Function in wasm file with the uid {self._wasmuid}:\n"
        for x in self._functions:
            result += f"function '{x}' with {self._functions[x][0]} parameter(s)"
            result += f" and {self._functions[x][1]} return value(s)\n"

        for x in self._unvalid_functions:
            result += f"unusable function: {x} \n"

        return result

    def check_function(
        self, function_name: str, number_of_parameters: int, number_of_returns: int
    ) -> bool:
        """checks a given function name and signature if it is included
        :param function_name: name of the function that is checked
        :type function_name: str
        :param number_of_parameters: number of i32 parameters of the function
        :type number_of_parameters: int
        :param number_of_returns: number of i32 return values of the function
        :type number_of_returns: int
        """
        if (
            (function_name in self._functions)
            and (self._functions[function_name][0] == number_of_parameters)
            and (self._functions[function_name][1] == number_of_returns)
        ):
            return True
        return False

# Copyright 2019-2024 Cambridge Quantum Computing
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

from qwasm import (  # type: ignore
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

    def __init__(self, filepath: str, check_file: bool = True, int_size: int = 32):
        """
        Construct a wasm file handler

        :param filepath: Path to the wasm file
        :type filepath: str
        :param check_file: If ``True`` checks file for compatibility with wasm
          standards. If ``False`` checks are skipped.
        :type check_file: bool
        :param int_size: length of the integer that is used in the wasm file
        :type int_size: int
        """
        self._int_size = int_size
        if int_size == 32:
            self._int_type = self.type_lookup[LANG_TYPE_I32]
        elif int_size == 64:
            self._int_type = self.type_lookup[LANG_TYPE_I64]
        else:
            raise ValueError(
                "given integer length not valid, only 32 and 64 are allowed"
            )

        self._filepath = filepath

        function_signatures: list = []
        function_names: list = []
        _func_lookup = {}

        if not exists(self._filepath):
            raise ValueError("wasm file not found at given path")

        with open(self._filepath, "rb") as file:
            self._wasm_file: bytes = file.read()

        self._wasm_file_encoded = base64.b64encode(self._wasm_file)

        self._wasmfileuid = hashlib.sha256(self._wasm_file_encoded).hexdigest()

        self._check_file = check_file

        # stores the names of the functions mapped
        #  to the number of parameters and the number of return values
        self._functions = dict()

        # contains the list of functions that are not allowed
        # to use in pytket (because of types that are not i32)
        self._unsupported_function = []

        if self._check_file:
            mod_iter = iter(decode_module(self._wasm_file))
            _, _ = next(mod_iter)

            for _, cur_sec_data in mod_iter:
                # read in list of function signatures
                if cur_sec_data.id == SEC_TYPE:
                    for idx, entry in enumerate(cur_sec_data.payload.entries):
                        function_signatures.append({})
                        function_signatures[idx]["parameter_types"] = [
                            self.type_lookup[pt] for pt in entry.param_types
                        ]
                        if entry.return_count > 1:
                            if (
                                isinstance(entry.return_type, list)
                                and len(entry.return_type) == entry.return_count
                            ):
                                function_signatures[idx]["return_types"] = [
                                    self.type_lookup[rt] for rt in entry.return_type
                                ]
                            elif isinstance(entry.return_type, int):
                                function_signatures[idx]["return_types"] = [
                                    self.type_lookup[entry.return_type]
                                ] * entry.return_count
                            else:
                                raise ValueError(
                                    f"Only parameter and return values of "
                                    + f"i{self._int_size} types are"
                                    + f" allowed, found type: {entry.return_type}"
                                )
                        elif entry.return_count == 1:
                            function_signatures[idx]["return_types"] = [
                                self.type_lookup[entry.return_type]
                            ]
                        else:
                            function_signatures[idx]["return_types"] = []

                # read in list of function names
                elif cur_sec_data.id == SEC_EXPORT:
                    f_idx = 0
                    for _, entry in enumerate(cur_sec_data.payload.entries):
                        if entry.kind == 0:
                            f_name = entry.field_str.tobytes().decode()
                            function_names.append(f_name)
                            _func_lookup[f_name] = (f_idx, entry.index)
                            f_idx += 1

                # read in map of function signatures to function names
                elif cur_sec_data.id == SEC_FUNCTION:
                    self._function_types = cur_sec_data.payload.types

            for x in function_names:
                # check for only integer type in parameters and return values
                supported_function = True
                idx = _func_lookup[x][1]

                if idx >= len(self._function_types):
                    raise ValueError("invalid wasm file")

                for t in function_signatures[self._function_types[idx]][
                    "parameter_types"
                ]:
                    if t != self._int_type:
                        supported_function = False
                for t in function_signatures[self._function_types[idx]]["return_types"]:
                    if t != self._int_type:
                        supported_function = False

                if (
                    len(function_signatures[self._function_types[idx]]["return_types"])
                    > 1
                ):
                    supported_function = False

                if supported_function:
                    self._functions[x] = (
                        len(
                            function_signatures[self._function_types[idx]][
                                "parameter_types"
                            ]
                        ),
                        len(
                            function_signatures[self._function_types[idx]][
                                "return_types"
                            ]
                        ),
                    )

                if not supported_function:
                    self._unsupported_function.append(x)

            if "init" not in self._functions:
                raise ValueError("wasm file needs to contain a function named 'init'")

            if self._functions["init"][0] != 0:
                raise ValueError("init function should not have any parameter")

            if self._functions["init"][1] != 0:
                raise ValueError("init function should not have any results")

    def __str__(self) -> str:
        """str representation of the wasm file"""
        return self._wasmfileuid

    def __repr__(self) -> str:
        """str representation of the contents of the wasm file"""
        if self._check_file:
            result = f"Functions in wasm file with the uid {self._wasmfileuid}:\n"
            for x in self._functions:
                result += f"function '{x}' with "
                result += f"{self._functions[x][0]} i{self._int_size} parameter(s)"
                result += (
                    f" and {self._functions[x][1]} i{self._int_size} return value(s)\n"
                )

            for x in self._unsupported_function:
                result += (
                    f"unsupported function with invalid "
                    f"parameter or result type: '{x}' \n"
                )

            return result
        else:
            raise ValueError(
                """the content of the wasm file can only be printed if the
wasm file was checked"""
            )

    def bytecode(self) -> bytes:
        """The WASM content as bytecode"""
        return self._wasm_file

    def check_function(
        self, function_name: str, number_of_parameters: int, number_of_returns: int
    ) -> bool:
        """
        Checks a given function name and signature if it is included

        :param function_name: name of the function that is checked
        :type function_name: str
        :param number_of_parameters: number of integer parameters of the function
        :type number_of_parameters: int
        :param number_of_returns: number of integer return values of the function
        :type number_of_returns: int
        :return: true if the signature and the name of the function is correct"""

        return not self._check_file or (
            (function_name in self._functions)
            and (self._functions[function_name][0] == number_of_parameters)
            and (self._functions[function_name][1] == number_of_returns)
        )

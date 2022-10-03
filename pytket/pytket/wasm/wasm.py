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

from wasmer import Store, Module, Instance, Function  # type: ignore


class WasmFileHandler:
    """Add a wasm file to your workflow, stores a copy of the file and
    checks the function signatures of the file. Offers function to add
    a wasm op to a circuit"""

    def __init__(self, filepath: str):
        """construct a wasm file handler"""
        self._filepath = filepath

        if not exists(self._filepath):
            raise ValueError("wasm file not found at given path")

        with open(self._filepath, "rb") as file:
            wasm_file: bytes = file.read()

        self._wasm_file_encoded = base64.b64encode(wasm_file)

        self._wasmuid = hashlib.md5(self._wasm_file_encoded).hexdigest()

        # check if the file is valid to run
        if not Module.validate(Store(), wasm_file):
            raise ValueError("wasm file not valid")

        wasm_module = Module(Store(), wasm_file)

        instance = Instance(wasm_module)

        # stores the names of the functions mapped
        #  to the number of parameters and the number of return values
        self._functions = dict()

        # contains the list of functions that are not allowed
        # to use in pytket (because of types that are not i32)
        self._unsupported_function = []

        self._wasm_file = base64.decodebytes(self._wasm_file_encoded)

        for wasm_obj in instance.exports:
            if len(wasm_obj) > 1 and isinstance(wasm_obj[1], Function):
                supported_function = True
                wasm_function = wasm_obj[1]

                # the direct evaluation of the types converts to python
                #  ints which remove the information if we are working
                #  with i32 or i64, so we are unfortunately required
                #  with the str version of the signature in the wasm file
                if (
                    (len(str(wasm_function.type).split("FunctionType(params: [")) == 2)
                    and (len(str(wasm_function.type).split("], results: [")) == 2)
                    and (len(str(wasm_function.type).split("])")) == 2)
                ):
                    wasm_parameter = (
                        str(wasm_function.type)
                        .split("FunctionType(params: [")[1]
                        .split("], results: [")[0]
                        .split(", ")
                    )

                    if wasm_parameter == [""]:
                        wasm_parameter = []

                    # special handling for no parameters
                    for t in wasm_parameter:
                        if t != "I32":
                            supported_function = False

                    wasm_results = (
                        str(wasm_function.type)
                        .split("], results: [")[1]
                        .split("])")[0]
                        .split(", ")
                    )

                    # special handling for void return
                    if wasm_results == [""]:
                        wasm_results = []

                    for t in wasm_results:
                        if t != "I32":
                            supported_function = False

                    if supported_function:
                        self._functions[wasm_obj[0]] = (
                            len(wasm_parameter),
                            len(wasm_results),
                        )
                else:
                    supported_function = False

                if not supported_function:
                    self._unsupported_function.append(wasm_obj[0])

    def __str__(self) -> str:
        """str representation of the wasm file"""
        return self._wasmuid

    def __repr__(self) -> str:
        """str representation of the containment of the wasm file"""
        result = f"Functions in wasm file with the uid {self._wasmuid}:\n"
        for x in self._functions:
            result += f"function '{x}' with {self._functions[x][0]} i32 parameter(s)"
            result += f" and {self._functions[x][1]} i32 return value(s)\n"

        for x in self._unsupported_function:
            result += (
                f"unsupported function with unvalid parameter or result type: '{x}' \n"
            )

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
        :return: true if the signature and the name of the function is correct"""
        return (
            (function_name in self._functions)
            and (self._functions[function_name][0] == number_of_parameters)
            and (self._functions[function_name][1] == number_of_returns)
        )

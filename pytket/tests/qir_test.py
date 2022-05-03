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

from pathlib import Path
import platform
import pytest
from pytket import Circuit

if platform.machine() == "x86_64":
    from pytket.qir.qir import (
        circuit_to_qir,
        circuit_from_qir,
        ExtendedModule,
        QIRParser,
        QUANTINUUM_GATES,
    )

    def test_raise_quantinuum_gateset_keyerror() -> None:
        c = Circuit(2)
        c.CY(0, 1)
        with pytest.raises(KeyError):
            circuit_to_qir(c, "RaiseError.ll", QUANTINUUM_GATES)

    def test_extended_module_for_quantinuum_gateset(
        ext_module_quantinuum_gateset: ExtendedModule,
    ) -> None:
        em = ext_module_quantinuum_gateset
        em_ir_str = em.module.ir()
        call_h = f"call void @__quantinuum__qis__h__body(%Qubit* null)"
        call_x = f"call void @__quantinuum__qis__x__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
        call_y = f"call void @__quantinuum__qis__y__body(%Qubit* null)"
        call_z = f"call void @__quantinuum__qis__z__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
        call_rx = f"call void @__quantinuum__qis__rx__body(double 0.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_ry = (
            f"call void @__quantinuum__qis__ry__body(double 1.000000e+00, %Qubit* null)"
        )
        call_rz = f"call void @__quantinuum__qis__rz__body(double 2.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_phx = f"call void @__quantinuum__qis__phx__body(double 1.000000e+00, double 2.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_cnot = f"call void @__quantinuum__qis__cnot__body(%Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_zzmax = f"call void @__quantinuum__qis__zzmax__body(%Qubit* inttoptr (i64 1 to %Qubit*), %Qubit* null)"
        call_zzph = f"call void @__quantinuum__qis__zzph__body(double 1.000000e+00, %Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_mz = (
            f"call void @__quantinuum__qis__mz__body(%Qubit* null, %Result* %zero)"
        )
        assert call_h in em_ir_str
        assert call_x in em_ir_str
        assert call_y in em_ir_str
        assert call_z in em_ir_str
        assert call_rx in em_ir_str
        assert call_ry in em_ir_str
        assert call_rz in em_ir_str
        assert call_phx in em_ir_str
        assert call_cnot in em_ir_str
        assert call_zzmax in em_ir_str
        assert call_zzph in em_ir_str
        assert call_mz in em_ir_str

    def test_qir_from_pytket_circuit_and_quantinuum_gateset(
        circuit_quantinuum_gateset, file_name: str
    ) -> None:
        with open(file_name, "r") as input:
            data = input.read()
        call_h = f"call void @__quantinuum__qis__h__body(%Qubit* null)"
        call_x = f"call void @__quantinuum__qis__x__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
        call_y = f"call void @__quantinuum__qis__y__body(%Qubit* null)"
        call_z = f"call void @__quantinuum__qis__z__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
        call_rx = f"call void @__quantinuum__qis__rx__body(double 0.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_ry = (
            f"call void @__quantinuum__qis__ry__body(double 1.000000e+00, %Qubit* null)"
        )
        call_rz = f"call void @__quantinuum__qis__rz__body(double 2.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_phx = f"call void @__quantinuum__qis__phx__body(double 1.500000e+00, double 5.000000e-01, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_cnot = f"call void @__quantinuum__qis__cnot__body(%Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_zzmax = f"call void @__quantinuum__qis__zzmax__body(%Qubit* inttoptr (i64 1 to %Qubit*), %Qubit* null)"
        call_zzph = f"call void @__quantinuum__qis__zzph__body(double 1.000000e+00, %Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_mz = f"call void @__quantinuum__qis__mz__body(%Qubit* inttoptr (i64 1 to %Qubit*), %Result* %zero)"

        assert call_h in data
        assert call_x in data
        assert call_y in data
        assert call_z in data
        assert call_rx in data
        assert call_ry in data
        assert call_rz in data
        assert call_phx in data
        assert call_cnot in data
        assert call_zzmax in data
        assert call_zzph in data
        assert call_mz in data

    def test_raise_pyqir_gateset_keyerror() -> None:
        c = Circuit(2)
        c.CY(0, 1)
        with pytest.raises(KeyError):
            circuit_to_qir(c, "RaiseError.ll")

    @pytest.mark.skip
    def test_qir_from_pytket_circuit_and_pyqir_gateset(
        circuit_pyqir_gateset, file_name: str
    ):
        with open(file_name, "r") as input:
            data = input.read()
        print(data)
        call_h = f"call void @__quantum__qis__h__body(%Qubit* null)"
        call_x = (
            f"call void @__quantum__qis__x__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
        )
        call_y = f"call void @__quantum__qis__y__body(%Qubit* null)"
        call_z = (
            f"call void @__quantum__qis__z__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
        )
        call_s = f"call void @__quantum__qis__s__body(%Qubit* null)"
        call_s_adj = (
            f"call void @__quantum__qis__s__adj(%Qubit* inttoptr (i64 1 to %Qubit*))"
        )
        call_t = f"call void @__quantum__qis__t__body(%Qubit* null)"
        call_t_adj = (
            f"call void @__quantum__qis__t__adj(%Qubit* inttoptr (i64 1 to %Qubit*))"
        )
        call_reset = f"call void @__quantum__qis__reset__body(%Qubit* null)"
        call_cnot = f"call void @__quantum__qis__cnot__body(%Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_cz = f"call void @__quantum__qis__cz__body(%Qubit* inttoptr (i64 1 to %Qubit*), %Qubit* null)"
        call_rx = f"call void @__quantum__qis__rx__body(double 0.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_ry = (
            f"call void @__quantum__qis__ry__body(double 1.000000e+00, %Qubit* null)"
        )
        call_rz = f"call void @__quantum__qis__rz__body(double 2.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
        call_m = f"call %Result* @__quantum__qis__m__body(%Qubit* inttoptr (i64 1 to %Qubit*))"

        assert call_h in data
        assert call_x in data
        assert call_y in data
        assert call_z in data
        assert call_s in data
        assert call_s_adj in data
        assert call_t in data
        assert call_t_adj in data
        assert call_reset in data
        assert call_cnot in data
        assert call_cz in data
        assert call_m in data
        assert call_rx in data
        assert call_ry in data
        assert call_rz in data

    @pytest.mark.skip
    def test_bitwise_ops_for_quantinuum_gateset(
        circuit_bitwise_ops: Circuit, bitwise_file: str
    ) -> None:
        with open(bitwise_file, "r") as input:
            data = input.read()
        print(data)
        call_and = f"call void @__quantinuum__cis__and__body(%Result* %zero, %Result* %zero1, %Result* %zero2)"
        call_or = f"call void @__quantinuum__cis__or__body(%Result* %zero3, %Result* %zero4, %Result* %zero5)"
        call_xor = f"call void @__quantinuum__cis__xor__body(%Result* %zero6, %Result* %zero7, %Result* %zero8)"

        assert call_and in data
        assert call_or in data
        assert call_xor in data

    def test_qir_parser(qir_bc_file_path: Path) -> None:
        qir_parser = QIRParser(str(qir_bc_file_path))
        func_name = (
            "Microsoft__Quantum__Samples__SimpleGrover__SearchForMarkedInput__Interop"
        )
        func = qir_parser.module.get_func_by_name(func_name)
        assert func.name == func_name
        assert qir_parser.get_required_qubits() == 3
        assert qir_parser.get_required_results() == 2

    def test_circuit_from_qir_from_bitcode(
        qir_bc_file_path: Path, qir_circuit: Circuit
    ) -> None:
        c = circuit_from_qir(qir_bc_file_path)
        assert c == qir_circuit

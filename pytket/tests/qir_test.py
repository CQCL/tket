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

from pytket import Circuit

from pytket.qir.qir import (
    circuit_to_qir_str,
    ExtendedModule,
    QUANTINUUM_GATES
)


def test_extended_module_for_quantinuum_gateset() -> None:
    em = ExtendedModule(
        name="Simple module for Quantinuum gateset.",
        num_qubits=2,
        num_results=1,
        gateset=QUANTINUUM_GATES
    )
    em.module.builder.call(em.h, [em.module.qubits[0]])
    em.module.builder.call(em.x, [em.module.qubits[1]])
    em.module.builder.call(em.y, [em.module.qubits[0]])
    em.module.builder.call(em.z, [em.module.qubits[1]])
    em.module.builder.call(
        em.rx, [
            0.0,
            em.module.qubits[1]
        ]
    )
    em.module.builder.call(
        em.ry, [
            1.0,
            em.module.qubits[0]
        ]
    )
    em.module.builder.call(
        em.rz, [
            2.0,
            em.module.qubits[1]
        ]
    )
    em.module.builder.call(
        em.phx, [
            1.0,
            2.0,
            em.module.qubits[1]
        ]
    )
    em.module.builder.call(
        em.cnot, [
            em.module.qubits[0],
            em.module.qubits[1]
        ]
    )
    em.module.builder.call(
        em.zzmax, [
            em.module.qubits[1],
            em.module.qubits[0]
        ]
    )
    em.module.builder.call(
        em.zzph, [
            1.0,
            em.module.qubits[0],
            em.module.qubits[1]
        ]
    )
    em.module.builder.call(
        em.mz, [
            em.module.qubits[0],
            em.module.results[0]
        ]
    )
    em_ir_str = em.module.ir()
    print(em_ir_str)
    call_h = f"call void @__quantinuum__qis__h__body(%Qubit* null)"
    call_x = f"call void @__quantinuum__qis__x__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
    call_y = f"call void @__quantinuum__qis__y__body(%Qubit* null)"
    call_z = f"call void @__quantinuum__qis__z__body(%Qubit* inttoptr (i64 1 to %Qubit*))"
    call_rx = f"call void @__quantinuum__qis__rx__body(double 0.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
    call_ry = f"call void @__quantinuum__qis__ry__body(double 1.000000e+00, %Qubit* null)"
    call_rz = f"call void @__quantinuum__qis__rz__body(double 2.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
    call_phx = f"call void @__quantinuum__qis__phx__body(double 1.000000e+00, double 2.000000e+00, %Qubit* inttoptr (i64 1 to %Qubit*))"
    call_cnot = f"call void @__quantinuum__qis__cnot__body(%Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
    call_zzmax = f"call void @__quantinuum__qis__zzmax__body(%Qubit* inttoptr (i64 1 to %Qubit*), %Qubit* null)"
    call_zzph = f"call void @__quantinuum__qis__zzph__body(double 1.000000e+00, %Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
    call_mz = f"call void @__quantinuum__qis__mz__body(%Qubit* null, %Result* %zero)"
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
    



def test_pyqir_builtin_gates_from_pytket_circuit() -> None:
    c = Circuit(2)
    c.CX(0, 1)
    mod_name = "Simple CX circuit"
    c_qir_str = circuit_to_qir_str(c, mod_name)
    print(c_qir_str)

def test_qir_str_measure() -> None:
    c = Circuit(1, 1)
    c.Measure(0, 0)
    module_name = "Simple Measure circuit"
    c_qir_str = circuit_to_qir_str(c, module_name)
    call = f"call %Result* @__quantum__qis__m__body(%Qubit* null)"
    assert call in c_qir_str


def test_qir_str_rz() -> None:
    c = Circuit(1)
    c.Rz(0.0, 0)
    module_name = "Simple Rz circuit"
    c_qir_str = circuit_to_qir_str(c, module_name)
    call = f"call void @__quantum__qis__rz__body(double 0.000000e+00, %Qubit* null)"
    assert call in c_qir_str


def test_extended_module() -> None:
    em = ExtendedModule('Extended Module', 2, 1)
    em.module.builder.call(em.zzmax, [em.module.qubits[0],em.module.qubits[1]])
    em_ir_str = em.module.ir()
    call = f"call void @__quantum__qis__zzmax__body(%Qubit* null, %Qubit* inttoptr (i64 1 to %Qubit*))"
    assert call in em_ir_str

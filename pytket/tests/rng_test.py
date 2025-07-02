# Copyright Quantinuum
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

import pytest

from pytket.circuit import Circuit, OpType
from pytket.qasm import circuit_to_qasm_str


def test_rng_seed() -> None:
    circ = Circuit()
    creg = circ.add_c_register("c", 64)
    circ.add_c_setbits([True, True], [creg[3], creg[11]])
    circ.set_rng_seed(creg)
    assert circ.n_gates_of_type(OpType.RNGSeed) == 1
    dreg = circ.add_c_register("d", 63)
    with pytest.raises(ValueError):
        circ.set_rng_seed(dreg)


def test_rng_bound() -> None:
    circ = Circuit()
    creg = circ.add_c_register("c", 32)
    circ.add_c_setbits([True, True], [creg[3], creg[11]])
    circ.set_rng_bound(creg)
    assert circ.n_gates_of_type(OpType.RNGBound) == 1
    dreg = circ.add_c_register("d", 64)
    with pytest.raises(ValueError):
        circ.set_rng_bound(dreg)


def test_rng_index() -> None:
    circ = Circuit()
    creg = circ.add_c_register("c", 32)
    circ.add_c_setbits([True, True], [creg[3], creg[11]])
    circ.set_rng_index(creg)
    assert circ.n_gates_of_type(OpType.RNGIndex) == 1
    dreg = circ.add_c_register("d", 64)
    with pytest.raises(ValueError):
        circ.set_rng_index(dreg)


def test_rng_num() -> None:
    circ = Circuit()
    creg = circ.add_c_register("c", 32)
    circ.get_rng_num(creg)
    assert circ.n_gates_of_type(OpType.RNGNum) == 1
    dreg = circ.add_c_register("d", 64)
    with pytest.raises(ValueError):
        circ.get_rng_num(dreg)


def test_rng() -> None:
    circ = Circuit(1, 1)
    seed = circ.add_c_register("seed", 64)
    bound = circ.add_c_register("bound", 32)
    index = circ.add_c_register("index", 32)
    num = circ.add_c_register("num", 32)
    circ.add_c_setbits([True, True, True], [seed[1], seed[3], seed[5]])
    circ.set_rng_seed(seed)
    circ.add_c_setbits([True], [bound[1]])
    circ.set_rng_bound(bound)
    circ.set_rng_index(index)
    circ.get_rng_num(num)
    circ.H(0, condition_bits=[num[0]], condition_value=1)
    circ.measure_all()
    assert circ.n_gates == 8
    qasm = circuit_to_qasm_str(circ, header="hqslib1", maxwidth=64)
    assert (
        qasm
        == """OPENQASM 2.0;
include "hqslib1.inc";

qreg q[1];
creg bound[32];
creg c[1];
creg index[32];
creg num[32];
creg seed[64];
bound[1] = 1;
seed[1] = 1;
seed[3] = 1;
seed[5] = 1;
RNGseed(seed);
RNGbound(bound);
RNGindex(index);
num = RNGnum();
if(num[0]==1) h q[0];
measure q[0] -> c[0];
"""
    )


def test_shot_num() -> None:
    circ = Circuit()
    creg = circ.add_c_register("c", 32)
    circ.get_job_shot_num(creg)
    assert circ.n_gates_of_type(OpType.JobShotNum) == 1
    dreg = circ.add_c_register("d", 16)
    with pytest.raises(ValueError):
        circ.get_job_shot_num(dreg)
    qasm = circuit_to_qasm_str(circ, header="hqslib1")
    assert (
        qasm
        == """OPENQASM 2.0;
include "hqslib1.inc";

creg c[32];
creg d[16];
c = JOB_shotnum;
"""
    )

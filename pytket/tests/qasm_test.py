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

from io import StringIO
import re
from pathlib import Path
from typing import List

import pytest

from pytket._tket.unit_id import _TEMP_BIT_NAME, _TEMP_BIT_REG_BASE
from pytket.circuit import (
    Circuit,
    OpType,
    fresh_symbol,
    Qubit,
    Bit,
    reg_eq,
    reg_neq,
    reg_lt,
    reg_gt,
    reg_leq,
    reg_geq,
    if_not_bit,
    BitRegister,
    CustomGate,
)
from pytket.circuit.decompose_classical import DecomposeClassicalError
from pytket.qasm import (
    circuit_from_qasm,
    circuit_to_qasm,
    circuit_from_qasm_str,
    circuit_to_qasm_str,
    circuit_from_qasm_wasm,
)
from pytket.qasm.qasm import QASMParseError, QASMUnsupportedError
from pytket.qasm.includes.load_includes import (
    _get_declpath,
    _get_files,
    _get_defpath,
    _write_defs,
    _write_decls,
    _load_gdict,
)
from pytket.transform import Transform
from pytket.passes import DecomposeClassicalExp, DecomposeBoxes

curr_file_path = Path(__file__).resolve().parent


def register_to_list(br: BitRegister) -> list[Bit]:
    return [br[i] for i in range(br.size)]


def test_qasm_correct() -> None:
    fname = str(curr_file_path / "qasm_test_files/test1.qasm")
    c = circuit_from_qasm(fname)
    assert c.n_qubits == 4
    assert c.depth() == 8
    coms = c.get_commands()
    assert len(coms) == 13
    correct_str = "[XXPhase(0.0375) q[0], q[1];, Rz(1.5) q[3];, ZZPhase(0.0375) q[0], q[1];, Rx(0.0375) q[3];, Rz(0.5) q[3];, CX q[0], q[3];, CZ q[0], q[1];, Rz(1.5) q[3];, Rx(1.9625) q[3];, CCX q[3], q[1], q[2];, Barrier q[0], q[3], q[2];, CU1(0.8) q[0], q[1];, U3(1, 0.5, 0.3) q[2];]"
    assert str(coms) == correct_str
    # TKET-871
    fname2 = str(curr_file_path / "qasm_test_files/test9.qasm")
    c2 = circuit_from_qasm(fname2)
    assert c2.n_qubits == 4
    assert c2.depth() == 2
    coms2 = c2.get_commands()
    assert len(coms2) == 3
    correct_str2 = "[CRz(0.3) q[0], q[1];, CRy(0.5) q[3], q[0];, CRx(0.5) q[2], q[1];]"
    assert str(coms2) == correct_str2
    # TKET-957
    fname3 = str(curr_file_path / "qasm_test_files/test10.qasm")
    c3 = circuit_from_qasm(fname3)
    assert c3.n_qubits == 2
    assert c3.depth() == 3
    coms3 = c3.get_commands()
    assert len(coms3) == 4
    correct_str3 = "[SX q[0];, X q[1];, SXdg q[1];, CSX q[0], q[1];]"
    assert str(coms3) == correct_str3


def test_qasm_qubit() -> None:
    with pytest.raises(RuntimeError) as errorinfo:
        fname = str(curr_file_path / "qasm_test_files/test2.qasm")
        circuit_from_qasm(fname)
    assert "Circuit does not contain unit with id: q[4]" in str(errorinfo.value)


def test_qasm_whitespace() -> None:
    fname = str(curr_file_path / "qasm_test_files/test3.qasm")
    c = circuit_from_qasm(fname)
    assert c.n_qubits == 10
    coms = c.get_commands()
    correct_str = "[Rz(0.5) q[3];, Rz(1.5) q[4];, Rx(0.085) q[7];, CX q[0], q[3];, CZ q[0], q[5];, Rz(1.5) q[3];, Rx(2.25) q[3];]"
    assert str(coms) == correct_str


def test_qasm_gate() -> None:
    with pytest.raises(QASMParseError) as errorinfo:
        fname = str(curr_file_path / "qasm_test_files/test4.qasm")
        circuit_from_qasm(fname)
    errormsg = str(errorinfo.value)
    assert "Cannot parse gate of type: gatedoesntexist" in errormsg
    assert "test4.qasm" in errormsg
    assert "Line:6" in errormsg


# @pytest.mark.skip()
def test_qasm_measure() -> None:
    fname = str(curr_file_path / "qasm_test_files/test5.qasm")
    circ = circuit_from_qasm(fname)
    assert str(circ.get_commands()) == "[X q[10];, Measure q[10] --> c[10];]"
    outf = str(curr_file_path / "qasm_test_files/testout4.qasm")
    circuit_to_qasm(circ, outf)
    circ = circuit_from_qasm(outf)
    assert str(circ.get_commands()) == "[X q[10];, Measure q[10] --> c[10];]"


def test_qasm_roundtrip() -> None:
    fname = str(curr_file_path / "qasm_test_files/test1.qasm")
    c = circuit_from_qasm(fname)
    qasm_out = str(curr_file_path / "qasm_test_files/testout.qasm")
    circuit_to_qasm(c, qasm_out)
    c2 = circuit_from_qasm(qasm_out)
    assert c == c2
    # TKET-871
    fname2 = str(curr_file_path / "qasm_test_files/test9.qasm")
    c3 = circuit_from_qasm(fname2)
    circuit_to_qasm(c3, qasm_out)
    c4 = circuit_from_qasm(qasm_out)
    assert c3 == c4


def test_qasm_str_roundtrip() -> None:
    with open(curr_file_path / "qasm_test_files/test1.qasm", "r") as f:
        c = circuit_from_qasm_str(f.read())
        qasm_str = circuit_to_qasm_str(c)
        c2 = circuit_from_qasm_str(qasm_str)
    assert c == c2


def test_qasm_str_roundtrip_oqc() -> None:
    with open(curr_file_path / "qasm_test_files/test15.qasm", "r") as f:
        c = circuit_from_qasm_str(f.read())
        qasm_str = circuit_to_qasm_str(c, "oqclib1")
        c2 = circuit_from_qasm_str(qasm_str)
    assert c == c2


def test_readout() -> None:
    fname = str(curr_file_path / "qasm_test_files/testout2.qasm")
    circ = Circuit(3)
    circ.X(0)
    circ.CX(1, 2)
    circ.H(2)
    circuit_to_qasm(circ, fname)
    circ2 = circuit_from_qasm(fname)
    assert circ2.depth() == 2
    assert circ2.n_gates == 3
    assert str(circ2.get_commands()) == "[X q[0];, CX q[1], q[2];, H q[2];]"


def test_symbolic_write() -> None:
    fname = str(curr_file_path / "qasm_test_files/testout3.qasm")
    circ = Circuit(2)
    a = fresh_symbol("a")
    circ.Rz(a, 0)
    circ.Rz(0.5, 0)
    circ.CX(0, 1)
    circ.Rx(0.5, 0)
    Transform.OptimisePostRouting().apply(circ)
    circuit_to_qasm(circ, fname)
    circ2 = circuit_from_qasm(fname)
    coms2 = circ2.get_commands()
    new_params = coms2[0].op.params
    assert new_params[2] == a + 1.0


def test_custom_gate() -> None:
    fname = str(curr_file_path / "qasm_test_files/test6.qasm")
    c = circuit_from_qasm(fname)
    assert str(c.get_commands()) == "[mygate(alpha,0.2) q[0], q[1];]"

    with open(curr_file_path / "qasm_test_files/test6_output.qasm") as f:
        # test custom gates are unrolled
        assert circuit_to_qasm_str(c) == f.read()
    Transform.DecomposeBoxes().apply(c)
    assert str(c.get_commands()) == "[Rz(alpha) q[0];, CX q[1], q[0];, Rx(0.2) q[1];]"


def test_register_commands() -> None:
    fname = str(curr_file_path / "qasm_test_files/test7.qasm")
    c = circuit_from_qasm(fname)
    assert (
        str(c.get_commands())
        == "[Barrier q[0], q[3], p[0], p[1];, U1(0.3) p[0];, U1(0.3) p[1];, CX p[0], r[0];, CX p[1], r[1];, Measure r[0] --> c[0];, Measure r[1] --> c[1];]"
    )


def test_conditional_gates() -> None:
    circ = Circuit(2, 2)
    circ.X(0)
    circ.Measure(0, 0)
    circ.Measure(1, 1)
    circ.Z(0, condition_bits=[0, 1], condition_value=2)
    circ.add_conditional_barrier([0, 1], [], [0, 1], 1)
    circ.Measure(0, 0, condition_bits=[0, 1], condition_value=1)
    qasm_out = str(curr_file_path / "qasm_test_files/testout5.qasm")
    circuit_to_qasm(circ, qasm_out)
    c2 = circuit_from_qasm(qasm_out)
    assert circ == c2


def test_named_conditional_barrier() -> None:
    circ = Circuit(2, 2)
    circ.add_bit(Bit("test", 3))
    circ.Z(0, condition_bits=[0, 1], condition_value=2)
    circ.add_conditional_barrier(
        [Qubit("q", 0), Bit("test", 3)],
        [Bit("c", 0), Bit("c", 1)],
        0,
        data="cond_barrier",
    )
    qs_str: str = circuit_to_qasm_str(circ)
    c_from_qs: Circuit = circuit_from_qasm_str(qs_str)
    assert qs_str == circuit_to_qasm_str(c_from_qs)


def test_hqs_conditional() -> None:
    c = Circuit(1)
    a = c.add_c_register("a", 8)
    b = c.add_c_register("b", 10)
    d = c.add_c_register("c", 10)

    c.add_c_setbits([True], [a[0]])
    c.add_c_setbits([False, True] + [False] * 6, register_to_list(a))
    c.add_c_setbits([True, True] + [False] * 8, register_to_list(b))

    c.X(0, condition=reg_eq(a ^ b, 1))
    c.X(0, condition=(a[0] ^ b[0]))
    c.X(0, condition=reg_eq(a & b, 1))
    c.X(0, condition=reg_eq(a | b, 1))

    c.X(0, condition=a[0])
    c.X(0, condition=reg_neq(a, 1))
    c.X(0, condition=if_not_bit(a[0]))
    c.X(0, condition=reg_gt(a, 1))
    c.X(0, condition=reg_lt(a, 1))
    c.X(0, condition=reg_geq(a, 1))
    c.X(0, condition=reg_leq(a, 1))

    DecomposeClassicalExp().apply(c)

    assert circuit_to_qasm_str(c, header="hqslib1")
    with pytest.raises(Exception) as errorinfo:
        circuit_to_qasm_str(c)
        assert "Complex classical gates not supported with qelib1" in str(
            errorinfo.value
        )

    copy1 = c.copy()
    copy1.X(0, condition_bits=[0, 1], condition_value=1)
    with pytest.raises(Exception) as errorinfo:
        circuit_to_qasm_str(copy1, header="hqslib1")
        assert "HQS OpenQASM conditions must act on one bit" in str(errorinfo.value)

    copy2 = c.copy()
    copy2.add_c_range_predicate(2, 5, register_to_list(d), b[0])
    with pytest.raises(Exception) as errorinfo:
        circuit_to_qasm_str(copy1, header="hqslib1")
        assert "Range can only be bounded on one side" in str(errorinfo.value)


def test_barrier() -> None:
    c = Circuit(3, 3)
    c.H(0)
    c.H(2)
    c.add_barrier([0], [0], "comment")
    result = """OPENQASM 2.0;\ninclude "hqslib1_dev.inc";\n\nqreg q[3];
creg c[3];\nh q[0];\nh q[2];\ncomment q[0],c[0];\n"""
    assert result == circuit_to_qasm_str(c, header="hqslib1_dev")


def test_barrier_2() -> None:
    c = Circuit(3, 3)
    c.H(0)
    c.H(2)
    c.add_barrier([0], [0], "different_comment )@#-(")

    result = """OPENQASM 2.0;\ninclude "hqslib1_dev.inc";\n\nqreg q[3];
creg c[3];\nh q[0];\nh q[2];\ndifferent_comment )@#-( q[0],c[0];\n"""
    assert result == circuit_to_qasm_str(c, header="hqslib1_dev")


def test_hqs_conditional_params() -> None:
    # https://github.com/CQCL/tket/issues/17
    c = Circuit(1, 1)
    c.add_gate(OpType.PhasedX, [1, 0], [0], condition_bits=[0])
    s = circuit_to_qasm_str(c, header="hqslib1")
    assert "U1q(1.0*pi,0.0*pi)" in s


def test_input_error_modes() -> None:
    with pytest.raises(QASMParseError):
        circuit_from_qasm_str('OPENQASM 2.0;\ninclude "qelib2.inc";')
    with pytest.raises(Exception) as errorinfo:
        circuit_from_qasm_str(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\n\ngate anrz(p) a {\n    rz(p) a;'
        )
    with pytest.raises(QASMParseError) as errorinfo:
        circuit_from_qasm_str(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nifrz(1.5*pi) q[3];'
        )
    assert "Cannot parse gate of type: ifrz" in str(errorinfo.value)
    with pytest.raises(Exception) as errorinfo:
        circuit_from_qasm_str(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\ncreg c[4];\nrogue q[3] -> c[2];'
        )
    assert "Unexpected token" in str(errorinfo.value)
    with pytest.raises(Exception) as errorinfo:
        circuit_from_qasm_str(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nrz(@£hog) q[2];'
        )
    assert "No terminal matches" in str(errorinfo.value)
    with pytest.raises(Exception) as errorinfo:
        circuit_from_qasm_str(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nnotanrz(0.3) q[2];'
        )
    assert "Cannot parse gate of type: notanrz" in str(errorinfo.value)
    with pytest.raises(Exception) as errorinfo:
        circuit_from_qasm_str(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nnotanx q[2];'
        )
    assert "Cannot parse gate of type: notanx" in str(errorinfo.value)


def test_output_error_modes() -> None:
    with pytest.raises(Exception) as errorinfo:
        c = Circuit()
        c.add_qubit(Qubit("q"))
        circuit_to_qasm_str(c)

    assert "OPENQASM registers must use a single index" in str(errorinfo.value)
    with pytest.raises(Exception) as errorinfo:
        c = Circuit(2, 2)
        c.CX(0, 1, condition_bits=[0], condition_value=0)
        circuit_to_qasm_str(c)
    assert "OpenQASM conditions must be an entire classical register" in str(
        errorinfo.value
    )
    with pytest.raises(Exception) as errorinfo:
        c = Circuit()
        qreg = c.add_q_register("q", 2)
        creg = c.add_c_register("c", 2)
        dreg = c.add_c_register("d", 2)
        c.CX(qreg[0], qreg[1], condition_bits=[creg[0], dreg[0]], condition_value=0)
        circuit_to_qasm_str(c)
    assert "OpenQASM conditions must be a single classical register" in str(
        errorinfo.value
    )


def test_header_stops_gate_definition() -> None:
    c = Circuit(2)
    c.add_gate(OpType.ZZMax, [0, 1])
    # adds a custom gate, "zzmax"
    qasm_str_qelib1 = circuit_to_qasm_str(c, header="qelib1")
    # adds "ZZ"
    qasm_str_hqslib1 = circuit_to_qasm_str(c, header="hqslib1")
    assert "gate zzmax zzmaxq0,zzmaxq1 {" in qasm_str_qelib1
    assert "}" in qasm_str_qelib1
    assert "zzmax q[0],q[1];" in qasm_str_qelib1
    assert "ZZ q[0],q[1];" in qasm_str_hqslib1
    assert circuit_from_qasm_str(qasm_str_qelib1) == circuit_from_qasm_str(
        qasm_str_hqslib1
    )


def test_builtin_gates() -> None:
    # QASM file using the built-in gates "CX" and "U" (output by staq)
    fname = str(curr_file_path / "qasm_test_files/test8.qasm")
    c = circuit_from_qasm(fname)
    assert "U3(0, 0, 1.75) q[2];" in str(c.get_commands())


def test_new_qelib1_aliases() -> None:
    # Check that aliases added to qelib1 parse into an equivalent instruction.
    fname = str(curr_file_path / "qasm_test_files/test16.qasm")
    c = circuit_from_qasm(fname)
    commands_str = str(c.get_commands())
    assert "U1(0) q[0]" in commands_str
    assert "U3(0, 0, 0) q[0]" in commands_str


def test_h1_rzz() -> None:
    c = (
        Circuit(2)
        .ZZPhase(0.3, 0, 1)
        .ZZPhase(2.4, 0, 1)
        .ZZPhase(1.4, 0, 1)
        .ZZPhase(1.0, 0, 1)
        .ZZPhase(-2.3, 0, 1)
        .ZZPhase(-1.4, 0, 1)
        .ZZPhase(-1.0, 0, 1)
    )
    assert "rzz" in circuit_to_qasm_str(c, header="qelib1")
    hqs_qasm_str = circuit_to_qasm_str(c, header="hqslib1")
    assert "RZZ" in hqs_qasm_str

    with open(str(curr_file_path / "qasm_test_files/zzphase.qasm"), "r") as f:
        fread_str = str(f.read())
        assert str(hqs_qasm_str) == fread_str


def test_extended_qasm() -> None:
    fname = str(curr_file_path / "qasm_test_files/test17.qasm")
    out_fname = str(curr_file_path / "qasm_test_files/test17_output.qasm")
    c = circuit_from_qasm_wasm(fname, "testfile.wasm")

    out_qasm = circuit_to_qasm_str(c, "hqslib1")
    with open(out_fname) as f:
        assert out_qasm == f.read()

    c2 = circuit_from_qasm_wasm(out_fname, "testfile.wasm")

    assert circuit_to_qasm_str(c2, "hqslib1")

    with pytest.raises(DecomposeClassicalError) as e:
        DecomposeClassicalExp().apply(c)


def test_decomposable_extended() -> None:
    fname = str(curr_file_path / "qasm_test_files/test18.qasm")
    out_fname = str(curr_file_path / "qasm_test_files/test18_output.qasm")

    c = circuit_from_qasm_wasm(fname, "testfile.wasm")
    DecomposeClassicalExp().apply(c)

    out_qasm = circuit_to_qasm_str(c, "hqslib1")
    with open(out_fname) as f:
        assert out_qasm == f.read()


def test_opaque() -> None:
    c = circuit_from_qasm_str(
        'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\nopaque myopaq() q1, q2;\n myopaq() q[0], q[1];'
    )
    with pytest.raises(QASMUnsupportedError) as e:
        circuit_to_qasm_str(c)
    assert "CustomGate myopaq has empty definition" in str(e.value)


def test_alternate_encoding() -> None:
    encoded_files = {
        "utf-16": str(curr_file_path / "qasm_test_files/utf16.qasm"),
        "utf-32": str(curr_file_path / "qasm_test_files/utf32.qasm"),
    }
    for enc, fil in encoded_files.items():
        with pytest.raises(Exception) as e:
            circuit_from_qasm(fil)
        c = circuit_from_qasm(fil, encoding=enc)
        assert c.n_gates == 6


def test_opaque_gates() -> None:
    # test opaque handled as barrier
    input_qasm = """OPENQASM 2.0;\ninclude "hqslib1_dev.inc";\n\nqreg q[4];
order2() q[0],q[2];\nrz(1.5*pi) q[3];\nsleep(1) q[0];\nrx(0.0375*pi) q[3];
cu1(0.8*pi) q[0],q[1];\nsleep(1,2) q[3];\n"""
    c = circuit_from_qasm_str(input_qasm)

    result_circ_qasm = circuit_to_qasm_str(c, "hqslib1_dev")
    assert input_qasm == result_circ_qasm


@pytest.mark.parametrize("fi", _get_files())
def test_include_loading(fi: Path) -> None:
    # If this test fails, you may have updated
    #  an .inc file and need to re-run load_includes.py
    id_match = r"'id': '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}'"
    gdict = _load_gdict(fi)
    for write_method, get_method in [
        (_write_defs, _get_defpath),
        (_write_decls, _get_declpath),
    ]:
        defs_stream = StringIO()
        write_method(defs_stream, gdict)
        defs_string = defs_stream.getvalue()
        with open(get_method(fi)) as f_def:
            fil_content = f_def.read()
        # get rid of uuids in the definitions
        defs_string = re.sub(id_match, "'id': ''", defs_string)
        fil_content = re.sub(id_match, "'id': ''", fil_content)
        # bytes comparison much faster than str comparison
        assert bytes(defs_string, "utf-8") == bytes(fil_content, "utf-8")


def test_custom_gate_with_barrier() -> None:
    c = circuit_from_qasm_str(
        """
    OPENQASM 2.0;
    include "qelib1.inc";

    gate my_gate q
    {
        barrier q;
    }

    qreg r[1];
    my_gate r[0];
    """
    )
    cmds = c.get_commands()
    assert len(cmds) == 1
    op = cmds[0].op
    assert op.type == OpType.CustomGate
    assert isinstance(op, CustomGate)
    opcirc = op.get_circuit()
    opcmds = opcirc.get_commands()
    assert len(opcmds) == 1
    assert opcmds[0].op.type == OpType.Barrier


def test_non_lib_gates() -> None:
    c = Circuit(3)
    c.add_gate(OpType.TK2, [0.2, 0.5, 0.7], [0, 1])
    c.add_gate(OpType.CV, [1, 2])
    c.add_gate(OpType.CVdg, [0, 1])
    c.add_gate(OpType.BRIDGE, [0, 1, 2])
    c.add_gate(OpType.ISWAP, [0.5], [0, 1])
    c.add_gate(OpType.PhasedISWAP, [0.2, 0.3], [0, 1])
    c.add_gate(OpType.YYPhase, [0.3], [0, 1])
    c.add_gate(OpType.XXPhase3, [0.4], [2, 0, 1])
    c.add_gate(OpType.ZZMax, [0, 1])
    c.add_gate(OpType.ESWAP, [0.7], [0, 1])
    c.add_gate(OpType.FSim, [0.3, 0.4], [1, 2])
    c.add_gate(OpType.ISWAPMax, [1, 2])
    c.add_gate(OpType.ECR, [2, 0])
    # add copies
    c.add_gate(OpType.TK2, [0.3, 0.6, 0.8], [2, 1])
    c.add_gate(OpType.CV, [1, 0])
    c.add_gate(OpType.CVdg, [1, 0])
    c.add_gate(OpType.BRIDGE, [0, 2, 1])
    c.add_gate(OpType.ISWAP, [0.7], [1, 2])
    c.add_gate(OpType.PhasedISWAP, [0.1, 0.2], [2, 1])
    c.add_gate(OpType.YYPhase, [0.4], [1, 2])
    c.add_gate(OpType.XXPhase3, [0.5], [2, 1, 0])
    c.add_gate(OpType.ZZMax, [1, 2])
    c.add_gate(OpType.ESWAP, [0.8], [2, 1])
    c.add_gate(OpType.FSim, [0.9, 0.2], [0, 2])
    c.add_gate(OpType.ISWAPMax, [0, 2])
    c.add_gate(OpType.ECR, [0, 1])

    qs = circuit_to_qasm_str(c)
    c2 = circuit_from_qasm_str(qs)
    qs2 = circuit_to_qasm_str(c2)
    assert qs == qs2


def test_scratch_bits_filtering() -> None:
    # test removing unused scratch register
    c = Circuit(1)
    a = c.add_c_register("a", 2)
    b = c.add_c_register("b", 3)
    c.add_c_copybits([b[1]], [a[1]], condition=reg_neq(b, 2))
    assert c.get_c_register(_TEMP_BIT_NAME)
    qstr = circuit_to_qasm_str(c, "hqslib1")
    assert _TEMP_BIT_NAME not in qstr
    qasm_out = str(curr_file_path / "qasm_test_files/testout6.qasm")
    circuit_to_qasm(c, qasm_out, "hqslib1")
    with open(qasm_out, "r") as f:
        assert _TEMP_BIT_NAME not in f.read()

    # test keeping used
    c = Circuit(1)
    a = c.add_c_register("a", 4)
    b = c.add_c_register("b", 3)
    c.X(0, condition=(a[0] ^ b[0]))
    assert c.get_c_register(_TEMP_BIT_NAME)
    qstr = circuit_to_qasm_str(c, "hqslib1")
    assert _TEMP_BIT_NAME in qstr
    circuit_to_qasm(c, qasm_out, "hqslib1")
    with open(qasm_out, "r") as f:
        assert _TEMP_BIT_NAME in f.read()

    # test multiple scratch registers
    c = circuit_from_qasm_str(
        f"""
    OPENQASM 2.0;
    include "hqslib1.inc";
    qreg q[1];
    creg a[1];
    creg b[1];
    creg d[2];
    creg {_TEMP_BIT_NAME}[32];
    creg {_TEMP_BIT_NAME}_1[32];
    {_TEMP_BIT_NAME}[0] = (a[0] ^ b[0]);
    if({_TEMP_BIT_NAME}[0]==1) x q[0];
    """
    )
    assert c.get_c_register(_TEMP_BIT_NAME)
    assert c.get_c_register(f"{_TEMP_BIT_NAME}_1")
    qstr = circuit_to_qasm_str(c, "hqslib1")
    assert _TEMP_BIT_NAME in qstr
    assert f"{_TEMP_BIT_NAME}_1" not in qstr
    circuit_to_qasm(c, qasm_out, "hqslib1")
    with open(qasm_out, "r") as f:
        fstr = f.read()
        assert _TEMP_BIT_NAME in fstr
        assert f"{_TEMP_BIT_NAME}_1" not in fstr

    # test leaving _TEMP_BIT_REG_BASE untouched
    qasm_str = f"""OPENQASM 2.0;
    include "hqslib1.inc";
    qreg q[1];
    creg a[8];
    creg b[10];
    creg d[10];
    creg {_TEMP_BIT_REG_BASE}_0[32];
    {_TEMP_BIT_REG_BASE}_0 = (a ^ b);
    """
    c = circuit_from_qasm_str(qasm_str)
    assert c.get_c_register(f"{_TEMP_BIT_REG_BASE}_0")
    qstr = circuit_to_qasm_str(c, "hqslib1")
    assert f"creg {_TEMP_BIT_REG_BASE}_0[32]" in qstr
    circuit_to_qasm(c, qasm_out, "hqslib1")
    with open(qasm_out, "r") as f:
        fstr = f.read()
        assert f"creg {_TEMP_BIT_REG_BASE}_0[32]" in fstr


def test_qasm_phase() -> None:
    c = Circuit(1, 1)
    c.H(0)
    c.Phase(0.5)
    c.Phase(0.5, condition_bits=[Bit(0)], condition_value=1)
    c.add_barrier([0])
    c.H(0)
    c.add_phase(0.5)
    c0 = Circuit(1, 1)
    c0.H(0)
    c0.add_barrier([0])
    c0.H(0)
    for header in ["qelib1", "hqslib1"]:
        qstr = circuit_to_qasm_str(c, header)
        c1 = circuit_from_qasm_str(qstr)
        assert c1 == c0


def test_CopyBits() -> None:
    input_qasm = """OPENQASM 2.0;\ninclude "hqslib1.inc";\n\ncreg c0[1];
creg c1[3];\nc0[0] = c1[1];\n"""
    c = circuit_from_qasm_str(input_qasm)
    result_circ_qasm = circuit_to_qasm_str(c, "hqslib1")
    assert input_qasm == result_circ_qasm

    c = Circuit()
    creg0 = c.add_c_register("c0", 1)
    creg1 = c.add_c_register("c1", 2)
    creg2 = c.add_c_register("c2", 2)
    # should output bit-wise assignment
    c.add_c_copybits([creg1[1]], [creg0[0]])
    # should output register-wise assignment
    c.add_c_copybits([creg2[0], creg2[1]], [creg1[0], creg1[1]])
    # should output bit-wise assignment
    c.add_c_copybits([creg2[0], creg2[1]], [creg1[1], creg1[0]])
    result_circ_qasm = circuit_to_qasm_str(c, "hqslib1")

    correct_qasm = """OPENQASM 2.0;\ninclude "hqslib1.inc";\n\ncreg c0[1];
creg c1[2];\ncreg c2[2];\nc0[0] = c1[1];\nc1 = c2;\nc1[1] = c2[0];\nc1[0] = c2[1];\n"""
    assert result_circ_qasm == correct_qasm


def test_RZZ_read_from() -> None:
    c = circuit_from_qasm_str(
        """
    OPENQASM 2.0;
    include "hqslib1.inc";

    qreg q[2];
    RZZ(0.5*pi) q[0],q[1];
    """
    )
    DecomposeBoxes().apply(c)
    assert "RZZ(0.5*pi) q[0],q[1];" in circuit_to_qasm_str(c, header="hqslib1")
    assert "rzz(0.5*pi) q[0],q[1];" in circuit_to_qasm_str(c)


def test_conditional_expressions() -> None:
    def cond_circ(bits: List[int]) -> Circuit:
        c = Circuit(4, 4)
        c.X(0)
        c.X(1)
        c.X(2)
        c.Measure(0, 0)
        c.Measure(1, 1)
        c.Measure(2, 2)
        c.X(3, condition_bits=bits, condition_value=3)
        c.Measure(3, 3)
        return c

    c1 = cond_circ([1])
    with pytest.raises(QASMUnsupportedError):
        circuit_to_qasm_str(c1)
    assert "if(c[1]==3)" in circuit_to_qasm_str(c1, header="hqslib1")
    c12 = cond_circ([1, 2])
    with pytest.raises(QASMUnsupportedError):
        circuit_to_qasm_str(c12)
    with pytest.raises(QASMUnsupportedError):
        circuit_to_qasm_str(c12, header="hqslib1")
    c0123 = cond_circ([0, 1, 2, 3])
    assert "if(c==3)" in circuit_to_qasm_str(c0123)
    assert "if(c==3)" in circuit_to_qasm_str(c0123, header="hqslib1")
    c0132 = cond_circ([0, 1, 3, 2])
    with pytest.raises(QASMUnsupportedError):
        circuit_to_qasm_str(c0132)
    with pytest.raises(QASMUnsupportedError):
        circuit_to_qasm_str(c0132, header="hqslib1")


def test_tk2_definition() -> None:
    c = Circuit(2).TK2(0.2, 0.2, 0.2, 0, 1)
    qs = circuit_to_qasm_str(c)
    assert "gate tk2 (param0, param1, param2) tk2q0,tk2q1 {" in qs
    assert "u3(3.5*pi,(param0/pi + 3.0)*pi,0.5*pi) tk2q0;" in qs
    assert "u3(1.0*pi,0.0*pi,(param1/pi + 1.0)*pi) tk2q1;" in qs
    assert "u3(0.0*pi,1.5*pi,(param2/pi + 0.5)*pi) tk2q1;" in qs
    assert "tk2(0.2*pi,0.2*pi,0.2*pi) q[0],q[1];" in qs


def test_rxxyyzz_conversion() -> None:
    c = circuit_from_qasm_str(
        """
    OPENQASM 2.0;
    include "hqslib1.inc";

    qreg q[2];
    Rxxyyzz(0.3*pi, 0.4*pi, 0.5*pi) q[0],q[1];
    """
    )
    c1 = Circuit(2).TK2(0.3, 0.4, 0.5, 0, 1)
    assert c == c1

    qelib_qs = circuit_to_qasm_str(c, header="qelib1")
    assert "gate tk2 (param0, param1, param2) tk2q0,tk2q1 {" in qelib_qs
    assert "tk2(0.3*pi,0.4*pi,0.5*pi) q[0],q[1];" in qelib_qs

    hqslib_qs = circuit_to_qasm_str(c, header="hqslib1")
    correct_qasm = """OPENQASM 2.0;\ninclude "hqslib1.inc";\n\nqreg q[2];\nRxxyyzz(0.3*pi,0.4*pi,0.5*pi) q[0],q[1];\n"""
    assert hqslib_qs == correct_qasm


def test_classical_assignment_order() -> None:
    # https://github.com/CQCL/tket/issues/1013
    circ = Circuit(1)
    reg_meas = circ.add_c_register("c0", 1)
    reg_cond = circ.add_c_register("c1", 1)
    circ.add_c_setreg(0, reg_cond)
    circ.X(0, condition=reg_eq(reg_cond, 1))
    circ.add_c_setreg(1, reg_cond)
    circ.X(0, condition=reg_eq(reg_cond, 1))
    circ.Measure(Qubit(0), reg_meas[0])
    qasm = circuit_to_qasm_str(circ, header="hqslib1").split("\n")
    right_order = [
        "c1 = 0;",
        "if(c1==1) tk_SCRATCH_BIT[0] = 1;",
        "c1 = 1;",
        "if(tk_SCRATCH_BIT[0]==1) x q[0];",
        "if(c1==1) x q[0];",
    ]
    posns = [qasm.index(line) for line in right_order]
    for i in range(len(right_order) - 1):
        assert posns[i] < posns[i + 1]


def test_classical_assignment_order_1() -> None:
    circ = Circuit(1)
    reg_meas = circ.add_c_register("c0", 1)
    reg_cond = circ.add_c_register("c1", 1)
    reg_aux = circ.add_c_register("c2", 1)
    circ.add_c_setreg(0, reg_cond)
    circ.add_c_setreg(1, reg_aux)
    circ.add_c_xor(reg_cond[0], reg_aux[0], reg_cond[0])
    circ.X(0, condition=reg_eq(reg_cond, 1))
    circ.add_c_xor(reg_cond[0], reg_aux[0], reg_cond[0])
    circ.Measure(Qubit(0), reg_meas[0])
    qasm = circuit_to_qasm_str(circ, header="hqslib1")
    correct_qasm = """OPENQASM 2.0;
include "hqslib1.inc";

qreg q[1];
creg c0[1];
creg c1[1];
creg c2[1];
creg tk_SCRATCH_BIT[1];
c1 = 0;
c2 = 1;
c1[0] = c2[0] ^ c1[0];
if(c1==1) tk_SCRATCH_BIT[0] = 1;
if(c1!=1) tk_SCRATCH_BIT[0] = 0;
c1[0] = c2[0] ^ c1[0];
if(tk_SCRATCH_BIT[0]==1) x q[0];
measure q[0] -> c0[0];
"""
    assert qasm == correct_qasm


def test_max_reg_width() -> None:
    circ_in = Circuit(1, 33)
    circ_in.H(0).Measure(0, 32)
    with pytest.raises(QASMUnsupportedError):
        circuit_to_qasm_str(circ_in)
    qasm_out = circuit_to_qasm_str(circ_in, maxwidth=64)
    assert "measure q[0] -> c[32];" in qasm_out
    qasm_in = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[1];
creg c[33];
h q[0];
measure q[0] -> c[32];"""
    with pytest.raises(QASMUnsupportedError):
        circuit_from_qasm_str(qasm_in)
    circ_out = circuit_from_qasm_str(qasm_in, maxwidth=64)
    assert len(circ_out.bits) == 33


def test_classical_expbox_arg_order() -> None:
    qasm = """
    OPENQASM 2.0;
    include "hqslib1.inc";
    
    qreg q[1];
    
    creg a[4];
    creg b[4];
    creg c[4];
    creg d[4];
    
    c = a ^ b | d;
    """

    circ = circuit_from_qasm_str(qasm)
    args = circ.get_commands()[0].args
    expected_symbol_order = ["a", "b", "d", "c"]
    expected_index_order = [0, 1, 2, 3]
    assert len(args) == 4 * 4
    arg_index = 0
    for symbol in expected_symbol_order:
        for index in expected_index_order:
            assert args[arg_index].reg_name == symbol
            assert args[arg_index].index[0] == index
            arg_index += 1


def test_register_name_check() -> None:
    # register names must have the expression [a-z][a-zA-Z0-9_]*
    qasm = """
    OPENQASM 2.0;
    include "hqslib1.inc";
    
    qreg Q[1];
    """
    with pytest.raises(QASMParseError) as e:
        circ = circuit_from_qasm_str(qasm)
    err_msg = "Invalid register definition 'Q[1]'"
    assert err_msg in str(e.value)

    c = Circuit()
    qb = Qubit("Q", 0)
    c.add_qubit(qb)
    c.H(qb)
    with pytest.raises(QASMUnsupportedError) as e2:
        qasm = circuit_to_qasm_str(c)
    err_msg = "Invalid register name 'Q'"
    assert err_msg in str(e2.value)


def test_conditional_custom() -> None:
    # https://github.com/CQCL/tket/issues/1299
    qasm0 = """
    OPENQASM 2.0;
    include "qelib1.inc";

    gate cx_o0 q0,q1 { x q0; cx q0,q1; x q0; }

    qreg q[2];
    creg c[1];

    if(c==1) cx_o0 q[1],q[0];
    measure q[0] -> c[0];
    """

    circ = circuit_from_qasm_str(qasm0)
    qasm1 = circuit_to_qasm_str(circ)
    assert (
        qasm1
        == """OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];
creg c[1];
if(c==1) x q[1];
if(c==1) cx q[1],q[0];
if(c==1) x q[1];
measure q[0] -> c[0];
"""
    )


if __name__ == "__main__":
    test_qasm_correct()
    test_qasm_qubit()
    test_qasm_whitespace()
    test_qasm_gate()
    test_qasm_measure()
    test_qasm_roundtrip()
    test_qasm_str_roundtrip()
    test_qasm_str_roundtrip_oqc()
    test_readout()
    test_symbolic_write()
    test_custom_gate()
    test_custom_gate_with_barrier()
    test_input_error_modes()
    test_output_error_modes()
    test_builtin_gates()
    test_new_qelib1_aliases()
    test_h1_rzz()
    test_opaque()
    test_opaque_gates()
    test_non_lib_gates()
    test_scratch_bits_filtering()
    test_extended_qasm()
    test_register_commands()
    test_conditional_gates()
    test_named_conditional_barrier()
    test_hqs_conditional()
    test_hqs_conditional_params()
    test_barrier()
    test_barrier_2()
    test_decomposable_extended()
    test_alternate_encoding()
    test_header_stops_gate_definition()
    test_tk2_definition()
    test_rxxyyzz_conversion()
    test_classical_expbox_arg_order()
    test_register_name_check()

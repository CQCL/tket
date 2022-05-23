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

from io import StringIO
import re
from pathlib import Path

import pytest  # type: ignore

from pytket.circuit import (  # type: ignore
    Circuit,
    OpType,
    fresh_symbol,
    Qubit,
    reg_eq,
    reg_neq,
    reg_lt,
    reg_gt,
    reg_leq,
    reg_geq,
    if_not_bit,
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
from pytket.transform import Transform  # type: ignore
from pytket.passes import DecomposeClassicalExp  # type: ignore

curr_file_path = Path(__file__).resolve().parent


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
    circ.Measure(0, 0, condition_bits=[0, 1], condition_value=1)
    qasm_out = str(curr_file_path / "qasm_test_files/testout5.qasm")
    circuit_to_qasm(circ, qasm_out)
    c2 = circuit_from_qasm(qasm_out)
    assert circ == c2


def test_hqs_conditional() -> None:
    c = Circuit(1)
    a = c.add_c_register("a", 8)
    b = c.add_c_register("b", 10)
    d = c.add_c_register("c", 10)

    c.add_c_setbits([True], [a[0]])
    c.add_c_setbits([False, True] + [False] * 6, list(a))
    c.add_c_setbits([True, True] + [False] * 8, list(b))

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
        assert "Complex classical gates only supported with hqslib1" in str(
            errorinfo.value
        )

    copy1 = c.copy()
    copy1.X(0, condition_bits=[0, 1], condition_value=1)
    with pytest.raises(Exception) as errorinfo:
        circuit_to_qasm_str(copy1, header="hqslib1")
        assert "HQS OpenQASM conditions must act on one bit" in str(errorinfo.value)

    copy2 = c.copy()
    copy2.add_c_range_predicate(2, 5, list(d), b[0])
    with pytest.raises(Exception) as errorinfo:
        circuit_to_qasm_str(copy1, header="hqslib1")
        assert "Range can only be bounded on one side" in str(errorinfo.value)


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
        c = Circuit(2)
        c.CV(0, 1)
        circuit_to_qasm_str(c)
    assert "Cannot print command of type: CV" in str(errorinfo.value)
    with pytest.raises(Exception) as errorinfo:
        c = Circuit(2)
        c.add_gate(OpType.ZZMax, [0, 1])
        circuit_to_qasm_str(c, header="qelib1")
    assert "Gate of type ZZ is not defined in header qelib1.inc" in str(errorinfo.value)
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
    c = Circuit(2)
    c.add_gate(OpType.ZZPhase, [0.1], [0, 1])
    assert "rzz" in circuit_to_qasm_str(c, header="qelib1")
    assert "RZZ" in circuit_to_qasm_str(c, header="hqslib1")


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
    assert "CustomGate myopaq has empty definition"


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


if __name__ == "__main__":
    test_qasm_correct()
    test_qasm_qubit()
    test_qasm_whitespace()
    test_qasm_gate()
    test_qasm_measure()
    test_qasm_roundtrip()
    test_qasm_str_roundtrip_oqc()
    test_readout()
    test_symbolic_write()
    test_custom_gate()
    test_input_error_modes()
    test_output_error_modes()
    test_builtin_gates()
    test_new_qelib1_aliases()

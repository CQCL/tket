OPENQASM 2.0;
include "qelib1.inc";
//Test controlled rotation gates
qreg q[4];
crz(0.3 * pi) q[0],q[1];
crx(0.5 * pi) q[2],q[1];
cry(0.5 * pi) q[3],q[0];

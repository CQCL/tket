OPENQASM 2.0;
include "qelib1.inc";
//test SX, SXdg, CSX gates
qreg q[2];
sx q[0];
x q[1];
sxdg q[1];
csx q[0],q[1];

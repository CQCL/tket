OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];
rz((alpha)*pi) q[0];
cx q[1],q[0];
rx(0.2*pi) q[1];

OPENQASM 2.0;
include "qelib1.inc";

qreg q[10];
creg c[4];
rz(1.5*pi) q[4];
rx(0.085*pi) q[7];
rz(0.5*pi) q[3];
cx q[0], q[3];
rz(1.5*pi) q[3];
rx(2.25*pi) q[3];
cz q[0] ,q[5];
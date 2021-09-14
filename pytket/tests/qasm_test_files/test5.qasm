OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];

x q[10];
measure q[10] -> c[10];
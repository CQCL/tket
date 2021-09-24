OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];
qreg p[2];
qreg r[2];
creg c[2];
barrier q[0],q[3],p;
u1(0.3*pi) p;
cx p, r;
measure r -> c;
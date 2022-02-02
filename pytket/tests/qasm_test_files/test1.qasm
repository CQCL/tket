OPENQASM 2.0;
include "qelib1.inc";
//some comments
qreg q[4];
rz(1.5*pi) q[3];
rx(0.0375*pi) q[3];
rxx(0.0375*pi) q[0],q[1];
rz(0.5*pi) q[3];
rzz(0.0375*pi) q[0],q[1];
cx q[0],q[3];
rz(1.5*pi) q[3];
rx(1.9625*pi) q[3];
cz q[0] ,q[1]; //hey look ma its a cz
ccx q[3],q[1],q[2];
barrier q[0],q[3],q[2];
u3(3.141596, 0.5* pi ,0.3*pi) q[2];
cu1(0.8*pi) q[0],q[1];

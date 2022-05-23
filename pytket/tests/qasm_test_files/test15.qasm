OPENQASM 2.0;
include "oqclib1.inc";
//some comments
qreg q[4];
rz(1.5*pi) q[3];
sx q[3];
rz(0.5*pi) q[3];
ecr q[0],q[3];
rz(1.5*pi) q[3];
sx q[3];
ecr q[0] ,q[1];
barrier q[0],q[3],q[2];
sx q[2];
ecr q[0],q[1];
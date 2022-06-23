OPENQASM 2.0;
include "hqslib1.inc";

qreg q[2];
RZZ(0.3*pi) q[0],q[1];
RZZ(0.4*pi) q[0],q[1];
RZZ(-0.6*pi) q[0],q[1];
RZZ(1.0*pi) q[0],q[1];
RZZ(-0.2999999999999998*pi) q[0],q[1];
RZZ(0.6*pi) q[0],q[1];
RZZ(1.0*pi) q[0],q[1];

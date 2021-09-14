OPENQASM 2.0;
include "qelib1.inc";

gate anrz(p) a {
    rz(p) a;
}

gate mygate(theta, phi) a, b {
    anrz(theta) a;
    cx b, a;
    rx(phi) b;
}

qreg q[2];
mygate(alpha*pi,0.2*pi) q[0], q[1];
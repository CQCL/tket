// file: oqclib1.inc

// --- Primitives ---

// 1-qubit Z-axis rotation
// NOTE: For HQS compiler / simulation, only one of these can be active at a
// time
opaque rz(lam) q;
// gate Rz(lam) q
//{
//    U(0,0,lam) q;
// }

// 1-qubit rotation gate
// NOTE: For HQS compiler / simulation, only one of these can be active at a
// time
opaque sx() q;
// gate U1q(theta, phi) q
//{
//    U(theta, phi-pi/2, pi/2-phi) q;
// }

// Unitary 2-qubit gate
// NOTE: For HQS compiler / simulation, only one of these can be active at a
// time
opaque ecr() q1, q2;
// gate ZZ() q1,q2
//{
//	U1q(pi/2, pi/2) q2;
//	CX q1, q2;
//	Rz(-pi/2) q1;
//	U1q(pi/2, 0) q2;
//	U1q(pi/2, -pi/2) q2;
// }

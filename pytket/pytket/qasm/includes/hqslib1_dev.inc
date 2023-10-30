// file: hqslib1_dev.inc

// --- Primitives ---

// 1-qubit Z-axis rotation
// NOTE: For HQS compiler / simulation, only one of these can be active at a time
opaque Rz(lam) q;
//gate Rz(lam) q
//{
//   U(0,0,lam) q;
//}

// 1-qubit rotation gate
// NOTE: For HQS compiler / simulation, only one of these can be active at a time
opaque U1q(theta, phi) q;
//gate U1q(theta, phi) q
//{
//   U(theta, phi-pi/2, pi/2-phi) q;
//}

// Unitary 2-qubit gate
// NOTE: For HQS compiler / simulation, only one of these can be active at a time
opaque ZZ() q1,q2;
//gate ZZ() q1,q2
//{
//	U1q(pi/2, pi/2) q2;
//	CX q1, q2;
//	Rz(-pi/2) q1;
//	U1q(pi/2, 0) q2;
//	U1q(pi/2, -pi/2) q2;
//}

opaque RZZ(theta) q1,q2;

opaque Rxxyyzz(alpha, beta, gamma) q1,q2;

opaque Rxxyyzz_zphase(alpha, beta, gamma, z0, z1) q1,q2;

// --- QE Hardware primitives ---

// 3-parameter 2-pulse single qubit gate
gate U(theta,phi,lam) q {
     Rz(lam) q;
     U1q(theta, pi/2) q;
     Rz(phi) q;
}

// --- Standard Gates ---

// Pauli gate: bit-flip
gate x() a
{
   U1q(pi, 0) a;
}
// Pauli gate: bit and phase flip
gate y() a
{
   U1q(pi, pi/2) a;
}
// Pauli gate: phase flip
gate z() a
{
   Rz(pi) a;
}
// Clifford gate: CNOT
gate CX() c,t
{
   U1q(-pi/2, pi/2) t;
   ZZ() c, t;
   Rz(-pi/2) c;
   U1q(pi/2, pi) t;
   Rz(-pi/2) t;
}
gate cx() c,t
{
   CX c,t;
}
// Clifford gate: Hadamard
gate h() a
{
   U1q(pi/2, -pi/2) a;
   Rz(pi) a;
}
// Clifford gate: sqrt(Z) phase gate
gate s() a
{
   Rz(pi/2) a;
}
// Clifford gate: conjugate of sqrt(Z)
gate sdg() a
{
   Rz(-pi/2) a;
}
// C3 gate: sqrt(S) phase gate
gate t() a
{
   Rz(pi/4) a;
}
// C3 gate: conjugate of sqrt(S)
gate tdg() a
{
   Rz(-pi/4) a;
}

// --- Standard rotations ---

// Rotation around X-axis
gate rx(theta) a
{
   U1q(theta, 0) a;
}
// Rotation around Y-axis
gate ry(theta) a
{
   U1q(theta, pi/2) a;
}
// Rotation around Z-axis
gate rz(phi) a
{
   Rz(phi) a;
}

// --- QE Standard User-Defined Gates  ---

// controlled-Phase
gate cz() a,b
{
   h b;
   cx a,b;
   h b;
}
// controlled-Y
gate cy() a,b
{
   sdg b;
   cx a,b;
   s b;
}
// controlled-H
gate ch() a,b
{
   h b; sdg b;
   cx a,b;
   h b; t b;
   cx a,b;
   t b; h b; s b; x b; s a;
}
// C3 gate: Toffoli
gate ccx() a,b,c
{
   h c;
   cx b,c; tdg c;
   cx a,c; t c;
   cx b,c; tdg c;
   cx a,c; t b; t c; h c;
   cx a,b; t a; tdg b;
   cx a,b;
}
// controlled rz rotation
gate crz(lam) a,b
{
   Rz(lam/2) b;
   cx a, b;
   Rz(-lam/2) b;
   cx a, b;
}
// controlled phase rotation
gate cu1(lam) a,b
{
   Rz(lam/2) a;
   cx a, b;
   Rz(-lam/2) b;
   cx a, b;
   Rz(lam/2) b;
}
// controlled-U
gate cu3(theta, phi, lam) c, t
{
   Rz((lam-phi)/2) t;
   cx c, t;
   Rz(-(phi+lam)/2) t;
   U1q(-theta/2, pi/2) t;
   cx c, t;
   U1q(theta/2, pi/2) t;
   Rz(phi) t;
}
// controlled-sqrt(X)
gate csx a,b { h b; cu1(pi/2) a,b; h b; }
// controlled-S
gate cs a,b {cu1(pi/2) a,b;}
// controlled-conjugate S
gate csdg a,b {cu1(-pi/2) a,b;}
// controlled-U gate
gate cu(theta,phi,lam,gamma) c, t
{ p(gamma) c;
  p((lam+phi)/2) c;
  p((lam-phi)/2) t;
  cx c,t;
  u(-theta/2,0,-(phi+lam)/2) t;
  cx c,t;
  u(theta/2,phi,0) t;
}
// generic single qubit gate
gate u(theta,phi,lam) q { U(theta,phi,lam) q; }
// 3-parameter 2-pulse single qubit gate
gate u3(theta,phi,lam) q {
     U(theta, phi, lam) q;
}
// 2-parameter 1-pulse single qubit gate
gate u2(phi,lam) q {
     u3(pi/2, phi, lam) q;
}
// 1-parameter 0-pulse single qubit gate
gate u1(lam) q { Rz(lam) q; }
// idle gate (identity)
gate id a { u3(0,0,0) a; }
// phase gate
gate p(lam) q { U(0,0,lam) q; }
// sqrt(X)
gate sx a { sdg a; h a; sdg a; }
// inverse sqrt(X)
gate sxdg a { s a; h a; s a; }
// swap
gate swap a,b { cx a,b; cx b,a; cx a,b; }
// cswap (Fredkin)
gate cswap a,b,c
{
  cx c,b;
  ccx a,b,c;
  cx c,b;
}
// controlled rx rotation
gate crx(lam) a,b
{
  u1(pi/2) b;
  cx a,b;
  u3(-lam/2,0,0) b;
  cx a,b;
  u3(lam/2,-pi/2,0) b;
}
// controlled ry rotation
gate cry(lam) a,b
{
  ry(lam/2) b;
  cx a,b;
  ry(-lam/2) b;
  cx a,b;
}
gate cp(lam) a,b
{
  p(lam/2) a;
  cx a,b;
  p(-lam/2) b;
  cx a,b;
  p(lam/2) b;
}
// two-qubit XX rotation
gate rxx(theta) a,b
{
  u3(pi/2, theta, 0) a;
  h b;
  cx a,b;
  u1(-theta) b;
  cx a,b;
  h b;
  u2(-pi, pi-theta) a;
}
gate rzz(theta) a,b
{
   RZZ(theta) a,b;
}
// relative-phase CCX
gate rccx a,b,c
{
  u2(0,pi) c;
  u1(pi/4) c;
  cx b, c;
  u1(-pi/4) c;
  cx a, c;
  u1(pi/4) c;
  cx b, c;
  u1(-pi/4) c;
  u2(0,pi) c;
}
// relative-phase 3-controlled X gate
gate rc3x a,b,c,d
{
  u2(0,pi) d;
  u1(pi/4) d;
  cx c,d;
  u1(-pi/4) d;
  u2(0,pi) d;
  cx a,d;
  u1(pi/4) d;
  cx b,d;
  u1(-pi/4) d;
  cx a,d;
  u1(pi/4) d;
  cx b,d;
  u1(-pi/4) d;
  u2(0,pi) d;
  u1(pi/4) d;
  cx c,d;
  u1(-pi/4) d;
  u2(0,pi) d;
}
// 3-controlled X gate
gate c3x a,b,c,d
{
    h d;
    p(pi/8) a;
    p(pi/8) b;
    p(pi/8) c;
    p(pi/8) d;
    cx a, b;
    p(-pi/8) b;
    cx a, b;
    cx b, c;
    p(-pi/8) c;
    cx a, c;
    p(pi/8) c;
    cx b, c;
    p(-pi/8) c;
    cx a, c;
    cx c, d;
    p(-pi/8) d;
    cx b, d;
    p(pi/8) d;
    cx c, d;
    p(-pi/8) d;
    cx a, d;
    p(pi/8) d;
    cx c, d;
    p(-pi/8) d;
    cx b, d;
    p(pi/8) d;
    cx c, d;
    p(-pi/8) d;
    cx a, d;
    h d;
}
// 3-controlled sqrt(X) gate, this equals the C3X gate where the CU1 rotations are -pi/8 not -pi/4
gate c3sqrtx a,b,c,d
{
    h d; cu1(-pi/8) a,d; h d;
    cx a,b;
    h d; cu1(pi/8) b,d; h d;
    cx a,b;
    h d; cu1(-pi/8) b,d; h d;
    cx b,c;
    h d; cu1(pi/8) c,d; h d;
    cx a,c;
    h d; cu1(-pi/8) c,d; h d;
    cx b,c;
    h d; cu1(pi/8) c,d; h d;
    cx a,c;
    h d; cu1(-pi/8) c,d; h d;
}
// 4-controlled X gate
gate c4x a,b,c,d,e
{
    h e; cu1(-pi/2) d,e; h e;
    c3x a,b,c,d;
    h e; cu1(pi/2) d,e; h e;
    c3x a,b,c,d;
    c3sqrtx a,b,c,e;
}
// --- Primitives ---

opaque sleep(t) q;

opaque order2() q0,q1;
opaque order3() q0,q1,q2;
opaque order4() q0,q1,q2,q3;
opaque order5() q0,q1,q2,q3,q4;
opaque order6() q0,q1,q2,q3,q4,q5;
opaque order7() q0,q1,q2,q3,q4,q5,q6;
opaque order8() q0,q1,q2,q3,q4,q5,q6,q7;
opaque order9() q0,q1,q2,q3,q4,q5,q6,q7,q8;
opaque order10() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9;
opaque order11() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10;
opaque order12() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11;
opaque order13() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12;
opaque order14() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13;
opaque order15() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14;
opaque order16() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15;
opaque order17() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16;
opaque order18() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17;
opaque order19() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18;
opaque order20() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,q19;

opaque group2() q0,q1;
opaque group3() q0,q1,q2;
opaque group4() q0,q1,q2,q3;
opaque group5() q0,q1,q2,q3,q4;
opaque group6() q0,q1,q2,q3,q4,q5;
opaque group7() q0,q1,q2,q3,q4,q5,q6;
opaque group8() q0,q1,q2,q3,q4,q5,q6,q7;
opaque group9() q0,q1,q2,q3,q4,q5,q6,q7,q8;
opaque group10() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9;
opaque group11() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10;
opaque group12() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11;
opaque group13() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12;
opaque group14() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13;
opaque group15() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14;
opaque group16() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15;
opaque group17() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16;
opaque group18() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17;
opaque group19() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18;
opaque group20() q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,q19;

//include "hqslib1.inc";

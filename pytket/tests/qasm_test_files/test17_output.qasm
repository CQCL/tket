OPENQASM 2.0;
include "hqslib1.inc";

qreg q[1];
creg a[2];
creg b[3];
creg c[4];
creg d[1];
creg tk_SCRATCH_BIT[2];
c = 2;
d = (a << 1);
c = a;
if(b!=2) c[1] = ((b[1] & a[1]) | a[0]);
c = (b & a);
b = (a + b);
b[1] = (b[0] ^ (~ b[2]));
c = (a - (b ** c));
d = (c >> 2);
c[0] = 1;
d = (a[0] ^ 1);
b = ((a * c) / b);
CCE1(c);
a = CCE2(a, b);
if(c>=2) h q[0];
if(d[0]==1) rx(1.0*pi) q[0];

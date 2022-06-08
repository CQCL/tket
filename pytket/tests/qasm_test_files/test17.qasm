OPENQASM 2.0;
include "hqslib1.inc";

extern [32] = CCE2([32], [32]);


qreg q[1];
creg c[4];
creg a[2];
creg b[3];
creg d[1];

c = 2;
c = a;
if (b != 2) c[1] = b[1] & a[1] | a[0];
c = b & a;
b = a + b;
b[1] = b[0] + b[2];
c = a - (b**c);
d = a << 1;
d = c >> 2;
b = a * c / b;
d[0] = a[0] ^ 1;
if(c>=2) h q[0];
if(d == 1) rx((0.5+0.5)*pi) q[0];
a = CCE2(a, b);
CCE1(c);

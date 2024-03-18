# Copyright 2019-2024 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


grammar = r"""
prog: oqasm? (incl|creg|qreg|gdef|opaq|reset|meas|barr|extern|ifc|cop|mixedcall|ccall)*

oqasm:  "OPENQASM" version ";"
?version: DECIMAL

incl: "include" "\"" INAM "\"" ";"
INAM: CNAME ("." CNAME)*

creg: "creg" _iarg ";"
qreg:  "qreg" _iarg ";"
mixedcall: _arg pars? margs ";"
ccall: _arg pars ";"
gatecall: id pars? args? ";"
gdef:  "gate" id pars? args "{" gatecall* "}"
opaq:  "opaque" id pars? args ";"
reset: "reset" _marg ";"
meas:  "measure" _marg "->" _marg ";"
barr:  "barrier" margs ";"
extern: "extern" bwargs "=" id "(" bwargs ")" ";"
ifc:   "if" (cond) (mixedcall|ccall|meas|reset|cop|barr)
cop:   assign
assign: margs "=" _exp ";"

!cond: "(" (_arg|_iarg) ("<"|">"|"<="|">="|"!="|"==") pint ")"

_exp: b_or | _xor_exp
b_or: _exp "|" _xor_exp

_xor_exp: xor | _b_and_exp
xor: _xor_exp "^" _b_and_exp

_b_and_exp: b_and | _shift_exp
b_and: _b_and_exp "&" _shift_exp

_shift_exp: lshift | rshift | _add_exp
lshift: _shift_exp "<<" _add_exp
rshift: _shift_exp ">>" _add_exp

_add_exp: add | sub | _prod_exp
add: _add_exp "+" _prod_exp
sub: _add_exp "-" _prod_exp

_prod_exp: mul | div | _neg_exp
mul: _prod_exp "*" _neg_exp
div: _prod_exp "/" _neg_exp

_neg_exp: neg | _b_not_exp
neg: "-" _neg_exp

_b_not_exp: b_not | _pow_exp
b_not: "~" _b_not_exp

_pow_exp: ipow | _atom_exp
ipow: _pow_exp "**" _atom_exp

_atom_exp: "(" _exp ")" 
    | cce_call
    | _atom

cce_call: ARG exp_args
_atom: _arg | _iarg | INT

_par_exp: par_pow | _par_add_exp

par_pow: _par_exp "**" _par_add_exp

_par_add_exp: par_add | par_sub | _par_prod_exp
par_add: _par_add_exp "+" _par_prod_exp
par_sub: _par_add_exp "-" _par_prod_exp

_par_prod_exp: par_mul | par_div | _par_neg_exp
par_mul: _par_prod_exp "*" _par_neg_exp
par_div: _par_prod_exp "/" _par_neg_exp

_par_neg_exp: par_neg | _par_atom_exp
par_neg: "-" _par_neg_exp

_par_atom_exp: "(" _par_exp ")"
    | sin | cos | tan | ln | sqrt
    | _fatom
sin:  "sin"  "(" _par_exp ")"
cos:  "cos"  "(" _par_exp ")"
tan:  "tan"  "(" _par_exp ")"
ln:   "ln"   "(" _par_exp ")"
sqrt: "sqrt" "(" _par_exp ")"

_fatom: snum | ARG

pars:  "(" _par_exp? ("," _par_exp)* ")"
exp_args: "(" _exp? ("," _exp)* ")"
args:  _arg ("," _arg)*
iargs: _iarg ("," _iarg)*
margs: _marg ("," _marg)*
bwargs: _bwarg ("," _bwarg)*

?id: CNAME
_arg: ARG
_iarg: IARG
_marg: (ARG|IARG)
_bwarg: "[" INT "]"
?snum: SIGNED_NUMBER
?pint: INT
ARG: CNAME
IARG: ARG ("[" INT "]")

COMMENT: "//" /[^\n]/*
%import common (WS, DECIMAL, INT, SIGNED_NUMBER, CNAME)
%ignore WS
%ignore COMMENT
"""

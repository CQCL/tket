OPENQASM 2.0;
include "hqslib1.inc";
// Generated using: PECOS version 0.6.0.dev3
creg m_bell[2];
creg m_out[1];
qreg sin_d[7];
qreg sin_a[3];
creg sin_c[32];
creg sin_syn_meas[32];
creg sin_last_raw_syn_x[32];
creg sin_last_raw_syn_z[32];
creg sin_scratch[32];
creg sin_flag_x[3];
creg sin_flags_z[3];
creg sin_flags[3];
creg sin_raw_meas[7];
creg sin_syn_x[3];
creg sin_syn_z[3];
creg sin_syndromes[3];
creg sin_verify_prep[32];
qreg smid_d[7];
qreg smid_a[3];
creg smid_c[32];
creg smid_syn_meas[32];
creg smid_last_raw_syn_x[32];
creg smid_last_raw_syn_z[32];
creg smid_scratch[32];
creg smid_flag_x[3];
creg smid_flags_z[3];
creg smid_flags[3];
creg smid_raw_meas[7];
creg smid_syn_x[3];
creg smid_syn_z[3];
creg smid_syndromes[3];
creg smid_verify_prep[32];
qreg sout_d[7];
qreg sout_a[3];
creg sout_c[32];
creg sout_syn_meas[32];
creg sout_last_raw_syn_x[32];
creg sout_last_raw_syn_z[32];
creg sout_scratch[32];
creg sout_flag_x[3];
creg sout_flags_z[3];
creg sout_flags[3];
creg sout_raw_meas[7];
creg sout_syn_x[3];
creg sout_syn_z[3];
creg sout_syndromes[3];
creg sout_verify_prep[32];

barrier smid_d[0], smid_d[1], smid_d[2], smid_d[3], smid_d[4], smid_d[5], smid_d[6], smid_a[0];

reset smid_d;
reset smid_a[0];
barrier smid_d, smid_a[0];
h smid_d[0];
h smid_d[4];
h smid_d[6];

cx smid_d[4], smid_d[5];
cx smid_d[0], smid_d[1];
cx smid_d[6], smid_d[3];
cx smid_d[4], smid_d[2];
cx smid_d[6], smid_d[5];
cx smid_d[0], smid_d[3];
cx smid_d[4], smid_d[1];
cx smid_d[3], smid_d[2];

barrier smid_a[0],smid_d[1],smid_d[3],smid_d[5];
//verification step
cx smid_d[5],smid_a[0];
cx smid_d[1],smid_a[0];
cx smid_d[3],smid_a[0];
measure smid_a[0] -> smid_c[0];


if(smid_c[0] == 1) barrier smid_d[0], smid_d[1], smid_d[2], smid_d[3], smid_d[4], smid_d[5], smid_d[6], smid_a[0];

if(smid_c[0] == 1) reset smid_d;
if(smid_c[0] == 1) reset smid_a[0];
if(smid_c[0] == 1) barrier smid_d, smid_a[0];
if(smid_c[0] == 1) h smid_d[0];
if(smid_c[0] == 1) h smid_d[4];
if(smid_c[0] == 1) h smid_d[6];

if(smid_c[0] == 1) cx smid_d[4], smid_d[5];
if(smid_c[0] == 1) cx smid_d[0], smid_d[1];
if(smid_c[0] == 1) cx smid_d[6], smid_d[3];
if(smid_c[0] == 1) cx smid_d[4], smid_d[2];
if(smid_c[0] == 1) cx smid_d[6], smid_d[5];
if(smid_c[0] == 1) cx smid_d[0], smid_d[3];
if(smid_c[0] == 1) cx smid_d[4], smid_d[1];
if(smid_c[0] == 1) cx smid_d[3], smid_d[2];

if(smid_c[0] == 1) barrier smid_a[0],smid_d[1],smid_d[3],smid_d[5];
//verification step
if(smid_c[0] == 1) cx smid_d[5],smid_a[0];
if(smid_c[0] == 1) cx smid_d[1],smid_a[0];
if(smid_c[0] == 1) cx smid_d[3],smid_a[0];
if(smid_c[0] == 1) measure smid_a[0] -> smid_c[0];


if(smid_c[0] == 1) barrier smid_d[0], smid_d[1], smid_d[2], smid_d[3], smid_d[4], smid_d[5], smid_d[6], smid_a[0];

if(smid_c[0] == 1) reset smid_d;
if(smid_c[0] == 1) reset smid_a[0];
if(smid_c[0] == 1) barrier smid_d, smid_a[0];
if(smid_c[0] == 1) h smid_d[0];
if(smid_c[0] == 1) h smid_d[4];
if(smid_c[0] == 1) h smid_d[6];

if(smid_c[0] == 1) cx smid_d[4], smid_d[5];
if(smid_c[0] == 1) cx smid_d[0], smid_d[1];
if(smid_c[0] == 1) cx smid_d[6], smid_d[3];
if(smid_c[0] == 1) cx smid_d[4], smid_d[2];
if(smid_c[0] == 1) cx smid_d[6], smid_d[5];
if(smid_c[0] == 1) cx smid_d[0], smid_d[3];
if(smid_c[0] == 1) cx smid_d[4], smid_d[1];
if(smid_c[0] == 1) cx smid_d[3], smid_d[2];

if(smid_c[0] == 1) barrier smid_a[0],smid_d[1],smid_d[3],smid_d[5];
//verification step
if(smid_c[0] == 1) cx smid_d[5],smid_a[0];
if(smid_c[0] == 1) cx smid_d[1],smid_a[0];
if(smid_c[0] == 1) cx smid_d[3],smid_a[0];
if(smid_c[0] == 1) measure smid_a[0] -> smid_c[0];



barrier sout_d[0], sout_d[1], sout_d[2], sout_d[3], sout_d[4], sout_d[5], sout_d[6], sout_a[0];

reset sout_d;
reset sout_a[0];
barrier sout_d, sout_a[0];
h sout_d[0];
h sout_d[4];
h sout_d[6];

cx sout_d[4], sout_d[5];
cx sout_d[0], sout_d[1];
cx sout_d[6], sout_d[3];
cx sout_d[4], sout_d[2];
cx sout_d[6], sout_d[5];
cx sout_d[0], sout_d[3];
cx sout_d[4], sout_d[1];
cx sout_d[3], sout_d[2];

barrier sout_a[0],sout_d[1],sout_d[3],sout_d[5];
//verification step
cx sout_d[5],sout_a[0];
cx sout_d[1],sout_a[0];
cx sout_d[3],sout_a[0];
measure sout_a[0] -> sout_c[0];


if(sout_c[0] == 1) barrier sout_d[0], sout_d[1], sout_d[2], sout_d[3], sout_d[4], sout_d[5], sout_d[6], sout_a[0];

if(sout_c[0] == 1) reset sout_d;
if(sout_c[0] == 1) reset sout_a[0];
if(sout_c[0] == 1) barrier sout_d, sout_a[0];
if(sout_c[0] == 1) h sout_d[0];
if(sout_c[0] == 1) h sout_d[4];
if(sout_c[0] == 1) h sout_d[6];

if(sout_c[0] == 1) cx sout_d[4], sout_d[5];
if(sout_c[0] == 1) cx sout_d[0], sout_d[1];
if(sout_c[0] == 1) cx sout_d[6], sout_d[3];
if(sout_c[0] == 1) cx sout_d[4], sout_d[2];
if(sout_c[0] == 1) cx sout_d[6], sout_d[5];
if(sout_c[0] == 1) cx sout_d[0], sout_d[3];
if(sout_c[0] == 1) cx sout_d[4], sout_d[1];
if(sout_c[0] == 1) cx sout_d[3], sout_d[2];

if(sout_c[0] == 1) barrier sout_a[0],sout_d[1],sout_d[3],sout_d[5];
//verification step
if(sout_c[0] == 1) cx sout_d[5],sout_a[0];
if(sout_c[0] == 1) cx sout_d[1],sout_a[0];
if(sout_c[0] == 1) cx sout_d[3],sout_a[0];
if(sout_c[0] == 1) measure sout_a[0] -> sout_c[0];


if(sout_c[0] == 1) barrier sout_d[0], sout_d[1], sout_d[2], sout_d[3], sout_d[4], sout_d[5], sout_d[6], sout_a[0];

if(sout_c[0] == 1) reset sout_d;
if(sout_c[0] == 1) reset sout_a[0];
if(sout_c[0] == 1) barrier sout_d, sout_a[0];
if(sout_c[0] == 1) h sout_d[0];
if(sout_c[0] == 1) h sout_d[4];
if(sout_c[0] == 1) h sout_d[6];

if(sout_c[0] == 1) cx sout_d[4], sout_d[5];
if(sout_c[0] == 1) cx sout_d[0], sout_d[1];
if(sout_c[0] == 1) cx sout_d[6], sout_d[3];
if(sout_c[0] == 1) cx sout_d[4], sout_d[2];
if(sout_c[0] == 1) cx sout_d[6], sout_d[5];
if(sout_c[0] == 1) cx sout_d[0], sout_d[3];
if(sout_c[0] == 1) cx sout_d[4], sout_d[1];
if(sout_c[0] == 1) cx sout_d[3], sout_d[2];

if(sout_c[0] == 1) barrier sout_a[0],sout_d[1],sout_d[3],sout_d[5];
//verification step
if(sout_c[0] == 1) cx sout_d[5],sout_a[0];
if(sout_c[0] == 1) cx sout_d[1],sout_a[0];
if(sout_c[0] == 1) cx sout_d[3],sout_a[0];
if(sout_c[0] == 1) measure sout_a[0] -> sout_c[0];


barrier smid_d, sout_d;
// Logical H
h smid_d;
// Transversal Logical CX
barrier smid_d, sout_d;
cx smid_d[0], sout_d[0];
cx smid_d[1], sout_d[1];
cx smid_d[2], sout_d[2];
cx smid_d[3], sout_d[3];
cx smid_d[4], sout_d[4];
cx smid_d[5], sout_d[5];
cx smid_d[6], sout_d[6];
barrier smid_d, sout_d;

smid_flag_x = 0;
smid_flags_z = 0;

// X check 1, Z check 2, Z check 3
// ===============================

reset smid_a[0];
reset smid_a[1];
reset smid_a[2];

h smid_a[0];
h smid_a[1];
h smid_a[2];

cx smid_a[0],smid_d[3];  // 5 -> 4
cz smid_a[1],smid_d[5];  // 6 -> 6
cz smid_a[2],smid_d[2];  // 7 -> 3

barrier smid_a[0],smid_a[1];
cz smid_a[0],smid_a[1];
barrier smid_a[0],smid_a[1];

cx smid_a[0],smid_d[0];  // 1 -> 1
cz smid_a[1],smid_d[4];  // 2 -> 5
cz smid_a[2],smid_d[3];  // 5 -> 4

cx smid_a[0],smid_d[1];  // 3 -> 2
cz smid_a[1],smid_d[2];  // 7 -> 3
cz smid_a[2],smid_d[6];  // 4 -> 7

barrier smid_a[0],smid_a[2];
cz smid_a[0],smid_a[2];
barrier smid_a[0],smid_a[2];

cx smid_a[0],smid_d[2];  // 7 -> 3
cz smid_a[1],smid_d[1];  // 3 -> 2
cz smid_a[2],smid_d[5];  // 6 -> 6

h smid_a[0];
h smid_a[1];
h smid_a[2];

measure smid_a[0] -> smid_flag_x[0];
measure smid_a[1] -> smid_flags_z[1];
measure smid_a[2] -> smid_flags_z[2];

smid_flag_x[0] = smid_flag_x[0] ^ smid_last_raw_syn_x[0];
smid_flags_z[1] = smid_flags_z[1] ^ smid_last_raw_syn_z[1];
smid_flags_z[2] = smid_flags_z[2] ^ smid_last_raw_syn_z[2];

smid_flags = smid_flag_x | smid_flags_z;


// Z check 1, X check 2, X check 3
// ===============================

if(smid_flags == 0) reset smid_a[0];
if(smid_flags == 0) reset smid_a[1];
if(smid_flags == 0) reset smid_a[2];

if(smid_flags == 0) h smid_a[0];
if(smid_flags == 0) h smid_a[1];
if(smid_flags == 0) h smid_a[2];


if(smid_flags == 0) barrier smid_a[0],smid_d[3];
if(smid_flags == 0) cz smid_a[0],smid_d[3];
if(smid_flags == 0) barrier smid_a[0],smid_d[3];

if(smid_flags == 0) barrier smid_a[1],smid_d[5];
if(smid_flags == 0) cx smid_a[1],smid_d[5];
if(smid_flags == 0) barrier smid_a[1],smid_d[5];

if(smid_flags == 0) barrier smid_a[2],smid_d[2];
if(smid_flags == 0) cx smid_a[2],smid_d[2];
if(smid_flags == 0) barrier smid_a[2],smid_d[2];



if(smid_flags == 0) barrier smid_a[0], smid_d[0], smid_d[1], smid_d[2], smid_d[3], smid_d[4], smid_d[5], smid_d[6], smid_a[1], smid_a[2];
if(smid_flags == 0) cz smid_a[1],smid_a[0];
if(smid_flags == 0) barrier smid_a[0], smid_d[0], smid_d[1], smid_d[2], smid_d[3], smid_d[4], smid_d[5], smid_d[6], smid_a[1], smid_a[2];


if(smid_flags == 0) barrier smid_a[0],smid_d[0];
if(smid_flags == 0) cz smid_a[0],smid_d[0];
if(smid_flags == 0) barrier smid_a[0],smid_d[0];

if(smid_flags == 0) barrier smid_a[1],smid_d[4];
if(smid_flags == 0) cx smid_a[1],smid_d[4];
if(smid_flags == 0) barrier smid_a[1],smid_d[4];

if(smid_flags == 0) barrier smid_a[2],smid_d[3];
if(smid_flags == 0) cx smid_a[2],smid_d[3];
if(smid_flags == 0) barrier smid_a[2],smid_d[3];



if(smid_flags == 0) barrier smid_a[0],smid_d[1];
if(smid_flags == 0) cz smid_a[0],smid_d[1];
if(smid_flags == 0) barrier smid_a[0],smid_d[1];

if(smid_flags == 0) barrier smid_a[1],smid_d[2];
if(smid_flags == 0) cx smid_a[1],smid_d[2];
if(smid_flags == 0) barrier smid_a[1],smid_d[2];

if(smid_flags == 0) barrier smid_a[2],smid_d[6];
if(smid_flags == 0) cx smid_a[2],smid_d[6];
if(smid_flags == 0) barrier smid_a[2],smid_d[6];


if(smid_flags == 0) barrier smid_a[0], smid_d[0], smid_d[1], smid_d[2], smid_d[3], smid_d[4], smid_d[5], smid_d[6], smid_a[1], smid_a[2];
if(smid_flags == 0) cz smid_a[2],smid_a[0];
if(smid_flags == 0) barrier smid_a[0], smid_d[0], smid_d[1], smid_d[2], smid_d[3], smid_d[4], smid_d[5], smid_d[6], smid_a[1], smid_a[2];



if(smid_flags == 0) barrier smid_a[0],smid_d[2];
if(smid_flags == 0) cz smid_a[0],smid_d[2];
if(smid_flags == 0) barrier smid_a[0],smid_d[2];

if(smid_flags == 0) barrier smid_a[1],smid_d[1];
if(smid_flags == 0) cx smid_a[1],smid_d[1];
if(smid_flags == 0) barrier smid_a[1],smid_d[1];

if(smid_flags == 0) barrier smid_a[2],smid_d[5];
if(smid_flags == 0) cx smid_a[2],smid_d[5];
if(smid_flags == 0) barrier smid_a[2],smid_d[5];


if(smid_flags == 0) h smid_a[0];
if(smid_flags == 0) h smid_a[1];
if(smid_flags == 0) h smid_a[2];



if(smid_flags == 0) measure smid_a[0] -> smid_flags_z[0];
if(smid_flags == 0) measure smid_a[1] -> smid_flag_x[1];
if(smid_flags == 0) measure smid_a[2] -> smid_flag_x[2];

// XOR flags/syndromes
if(smid_flags == 0) smid_flags_z[0] = smid_flags_z[0] ^ smid_last_raw_syn_z[0];
if(smid_flags == 0) smid_flag_x[1] = smid_flag_x[1] ^ smid_last_raw_syn_x[1];
if(smid_flags == 0) smid_flag_x[2] = smid_flag_x[2] ^ smid_last_raw_syn_x[2];

if(smid_flags == 0) smid_flags = smid_flag_x | smid_flags_z;


// Run the 6 non-flagged checks (if non-trivial flags)
// ===================================================
// // X check 1, Z check 2, Z check 3

if(smid_flags != 0) smid_syn_x = 0;
if(smid_flags != 0) smid_syn_z = 0;

if(smid_flags != 0) reset smid_a[0];
if(smid_flags != 0) reset smid_a[1];
if(smid_flags != 0) reset smid_a[2];

if(smid_flags != 0) h smid_a[0];
if(smid_flags != 0) h smid_a[1];
if(smid_flags != 0) h smid_a[2];

if(smid_flags != 0) cx smid_a[0],smid_d[2];
if(smid_flags != 0) cz smid_a[1],smid_d[5];
if(smid_flags != 0) cz smid_a[2],smid_d[6];

if(smid_flags != 0) cx smid_a[0],smid_d[1];
if(smid_flags != 0) cz smid_a[1],smid_d[2];
if(smid_flags != 0) cz smid_a[2],smid_d[5];

if(smid_flags != 0) cx smid_a[0],smid_d[3];
if(smid_flags != 0) cz smid_a[1],smid_d[1];
if(smid_flags != 0) cz smid_a[2],smid_d[2];

if(smid_flags != 0) cx smid_a[0],smid_d[0];
if(smid_flags != 0) cz smid_a[1],smid_d[4];
if(smid_flags != 0) cz smid_a[2],smid_d[3];

if(smid_flags != 0) h smid_a[0];
if(smid_flags != 0) h smid_a[1];
if(smid_flags != 0) h smid_a[2];

if(smid_flags != 0) measure smid_a[0] -> smid_syn_x[0];
if(smid_flags != 0) measure smid_a[1] -> smid_syn_z[1];
if(smid_flags != 0) measure smid_a[2] -> smid_syn_z[2];

// // Z check 1, X check 2, X check 3

if(smid_flags != 0) reset smid_a[0];
if(smid_flags != 0) reset smid_a[1];
if(smid_flags != 0) reset smid_a[2];

if(smid_flags != 0) h smid_a[0];
if(smid_flags != 0) h smid_a[1];
if(smid_flags != 0) h smid_a[2];

if(smid_flags != 0) cz smid_a[0],smid_d[2];
if(smid_flags != 0) cx smid_a[1],smid_d[5];
if(smid_flags != 0) cx smid_a[2],smid_d[6];

if(smid_flags != 0) cz smid_a[0],smid_d[1];
if(smid_flags != 0) cx smid_a[1],smid_d[2];
if(smid_flags != 0) cx smid_a[2],smid_d[5];

if(smid_flags != 0) cz smid_a[0],smid_d[3];
if(smid_flags != 0) cx smid_a[1],smid_d[1];
if(smid_flags != 0) cx smid_a[2],smid_d[2];

if(smid_flags != 0) cz smid_a[0],smid_d[0];
if(smid_flags != 0) cx smid_a[1],smid_d[4];
if(smid_flags != 0) cx smid_a[2],smid_d[3];

if(smid_flags != 0) h smid_a[0];
if(smid_flags != 0) h smid_a[1];
if(smid_flags != 0) h smid_a[2];

if(smid_flags != 0) measure smid_a[0] -> smid_syn_z[0];
if(smid_flags != 0) measure smid_a[1] -> smid_syn_x[1];
if(smid_flags != 0) measure smid_a[2] -> smid_syn_x[2];


// =========================
// BEGIN Run X decoder
// =========================

if(smid_flags!=0) smid_syndromes = smid_syn_x ^ smid_last_raw_syn_x;
if(smid_flags==0) smid_syndromes = 0;

// apply corrections
if(smid_syndromes == 2) smid_c[4] = smid_c[4] ^ 1;
if(smid_syndromes == 4) smid_c[4] = smid_c[4] ^ 1;
if(smid_syndromes == 6) smid_c[4] = smid_c[4] ^ 1;

// alter correction based on flags
// ===============================

// 1&2 (1 -> 2)
// ------------
smid_scratch = 0;
if(smid_flag_x == 1) smid_scratch[0] = 1;
if(smid_syndromes == 2) smid_scratch[1] = 1;

smid_scratch[2] = smid_scratch[0] & smid_scratch[1];
if(smid_scratch[2] == 1) smid_c[4] = smid_c[4] ^ 1;

// 1&4 (1 -> 3)
// ------------
smid_scratch = 0;
if(smid_flag_x == 1) smid_scratch[0] = 1;
if(smid_syndromes == 4) smid_scratch[1] = 1;

smid_scratch[2] = smid_scratch[0] & smid_scratch[1];
if(smid_scratch[2] == 1) smid_c[4] = smid_c[4] ^ 1;


// 6&4 (2,3 -> 3)
// ------------
smid_scratch = 0;
if(smid_flag_x == 6) smid_scratch[0] = 1;
if(smid_syndromes == 4) smid_scratch[1] = 1;

smid_scratch[2] = smid_scratch[0] & smid_scratch[1];
if(smid_scratch[2] == 1) smid_c[4] = smid_c[4] ^ 1;

if(smid_flags!=0) smid_last_raw_syn_x = smid_syn_x;

// =========================
// END Run X decoder
// =========================



// ACTIVE ERROR CORRECTION FOR X SYNDROMES

smid_scratch  = 0;

if(smid_syndromes[0] == 1) smid_scratch = smid_scratch  ^ 1;  // only part that differs for X vs Z syns
if(smid_syndromes[1] == 1) smid_scratch  = smid_scratch  ^ 12;
if(smid_syndromes[2] == 1) smid_scratch  = smid_scratch  ^ 48;

if(smid_c[4]==1) smid_scratch  = smid_scratch  ^ 112;  // logical operator

if(smid_scratch[0] == 1) z smid_d[0];
// if(smid_scratch[1] == 1) z smid_d[1];  // not possible for X stabilizers
if(smid_scratch[2] == 1) z smid_d[2];
if(smid_scratch[3] == 1) z smid_d[3];
if(smid_scratch[4] == 1) z smid_d[4];
if(smid_scratch[5] == 1) z smid_d[5];
if(smid_scratch[6] == 1) z smid_d[6];

smid_c[4] = 0;
// smid_syndromes = 0;
smid_last_raw_syn_x = 0;
// smid_syn_x = 0;
// smid_flag_x = 0;
// smid_flags = 0;



// =========================
// BEGIN Run Z decoder
// =========================

if(smid_flags!=0) smid_syndromes = smid_syn_z ^ smid_last_raw_syn_z;
if(smid_flags==0) smid_syndromes = 0;

// apply corrections
if(smid_syndromes == 2) smid_c[3] = smid_c[3] ^ 1;
if(smid_syndromes == 4) smid_c[3] = smid_c[3] ^ 1;
if(smid_syndromes == 6) smid_c[3] = smid_c[3] ^ 1;

// alter correction based on flags
// ===============================

// 1&2 (1 -> 2)
// ------------
smid_scratch = 0;
if(smid_flags_z == 1) smid_scratch[0] = 1;
if(smid_syndromes == 2) smid_scratch[1] = 1;

smid_scratch[2] = smid_scratch[0] & smid_scratch[1];
if(smid_scratch[2] == 1) smid_c[3] = smid_c[3] ^ 1;

// 1&4 (1 -> 3)
// ------------
smid_scratch = 0;
if(smid_flags_z == 1) smid_scratch[0] = 1;
if(smid_syndromes == 4) smid_scratch[1] = 1;

smid_scratch[2] = smid_scratch[0] & smid_scratch[1];
if(smid_scratch[2] == 1) smid_c[3] = smid_c[3] ^ 1;


// 6&4 (2,3 -> 3)
// ------------
smid_scratch = 0;
if(smid_flags_z == 6) smid_scratch[0] = 1;
if(smid_syndromes == 4) smid_scratch[1] = 1;

smid_scratch[2] = smid_scratch[0] & smid_scratch[1];
if(smid_scratch[2] == 1) smid_c[3] = smid_c[3] ^ 1;

if(smid_flags!=0) smid_last_raw_syn_z = smid_syn_z;

// =========================
// END Run Z decoder
// =========================



// ACTIVE ERROR CORRECTION FOR Z SYNDROMES

smid_scratch  = 0;

if(smid_syndromes[0] == 1) smid_scratch = smid_scratch  ^ 14;  // only part that differs for X vs Z syns
if(smid_syndromes[1] == 1) smid_scratch  = smid_scratch  ^ 12;
if(smid_syndromes[2] == 1) smid_scratch  = smid_scratch  ^ 48;

if(smid_c[3]==1) smid_scratch  = smid_scratch  ^ 112;  // logical operator

// if(smid_scratch[0] == 1) z smid_d[0]; // not possible for X stabilizers
if(smid_scratch[1] == 1) x smid_d[1];
if(smid_scratch[2] == 1) x smid_d[2];
if(smid_scratch[3] == 1) x smid_d[3];
if(smid_scratch[4] == 1) x smid_d[4];
if(smid_scratch[5] == 1) x smid_d[5];
if(smid_scratch[6] == 1) x smid_d[6];

smid_c[3] = 0;
// smid_syndromes = 0;
smid_last_raw_syn_z = 0;
// smid_syn_z = 0;
// smid_flags_z = 0;
// smid_flags = 0;


sout_flag_x = 0;
sout_flags_z = 0;

// X check 1, Z check 2, Z check 3
// ===============================

reset sout_a[0];
reset sout_a[1];
reset sout_a[2];

h sout_a[0];
h sout_a[1];
h sout_a[2];

cx sout_a[0],sout_d[3];  // 5 -> 4
cz sout_a[1],sout_d[5];  // 6 -> 6
cz sout_a[2],sout_d[2];  // 7 -> 3

barrier sout_a[0],sout_a[1];
cz sout_a[0],sout_a[1];
barrier sout_a[0],sout_a[1];

cx sout_a[0],sout_d[0];  // 1 -> 1
cz sout_a[1],sout_d[4];  // 2 -> 5
cz sout_a[2],sout_d[3];  // 5 -> 4

cx sout_a[0],sout_d[1];  // 3 -> 2
cz sout_a[1],sout_d[2];  // 7 -> 3
cz sout_a[2],sout_d[6];  // 4 -> 7

barrier sout_a[0],sout_a[2];
cz sout_a[0],sout_a[2];
barrier sout_a[0],sout_a[2];

cx sout_a[0],sout_d[2];  // 7 -> 3
cz sout_a[1],sout_d[1];  // 3 -> 2
cz sout_a[2],sout_d[5];  // 6 -> 6

h sout_a[0];
h sout_a[1];
h sout_a[2];

measure sout_a[0] -> sout_flag_x[0];
measure sout_a[1] -> sout_flags_z[1];
measure sout_a[2] -> sout_flags_z[2];

sout_flag_x[0] = sout_flag_x[0] ^ sout_last_raw_syn_x[0];
sout_flags_z[1] = sout_flags_z[1] ^ sout_last_raw_syn_z[1];
sout_flags_z[2] = sout_flags_z[2] ^ sout_last_raw_syn_z[2];

sout_flags = sout_flag_x | sout_flags_z;


// Z check 1, X check 2, X check 3
// ===============================

if(sout_flags == 0) reset sout_a[0];
if(sout_flags == 0) reset sout_a[1];
if(sout_flags == 0) reset sout_a[2];

if(sout_flags == 0) h sout_a[0];
if(sout_flags == 0) h sout_a[1];
if(sout_flags == 0) h sout_a[2];


if(sout_flags == 0) barrier sout_a[0],sout_d[3];
if(sout_flags == 0) cz sout_a[0],sout_d[3];
if(sout_flags == 0) barrier sout_a[0],sout_d[3];

if(sout_flags == 0) barrier sout_a[1],sout_d[5];
if(sout_flags == 0) cx sout_a[1],sout_d[5];
if(sout_flags == 0) barrier sout_a[1],sout_d[5];

if(sout_flags == 0) barrier sout_a[2],sout_d[2];
if(sout_flags == 0) cx sout_a[2],sout_d[2];
if(sout_flags == 0) barrier sout_a[2],sout_d[2];



if(sout_flags == 0) barrier sout_a[0], sout_d[0], sout_d[1], sout_d[2], sout_d[3], sout_d[4], sout_d[5], sout_d[6], sout_a[1], sout_a[2];
if(sout_flags == 0) cz sout_a[1],sout_a[0];
if(sout_flags == 0) barrier sout_a[0], sout_d[0], sout_d[1], sout_d[2], sout_d[3], sout_d[4], sout_d[5], sout_d[6], sout_a[1], sout_a[2];


if(sout_flags == 0) barrier sout_a[0],sout_d[0];
if(sout_flags == 0) cz sout_a[0],sout_d[0];
if(sout_flags == 0) barrier sout_a[0],sout_d[0];

if(sout_flags == 0) barrier sout_a[1],sout_d[4];
if(sout_flags == 0) cx sout_a[1],sout_d[4];
if(sout_flags == 0) barrier sout_a[1],sout_d[4];

if(sout_flags == 0) barrier sout_a[2],sout_d[3];
if(sout_flags == 0) cx sout_a[2],sout_d[3];
if(sout_flags == 0) barrier sout_a[2],sout_d[3];



if(sout_flags == 0) barrier sout_a[0],sout_d[1];
if(sout_flags == 0) cz sout_a[0],sout_d[1];
if(sout_flags == 0) barrier sout_a[0],sout_d[1];

if(sout_flags == 0) barrier sout_a[1],sout_d[2];
if(sout_flags == 0) cx sout_a[1],sout_d[2];
if(sout_flags == 0) barrier sout_a[1],sout_d[2];

if(sout_flags == 0) barrier sout_a[2],sout_d[6];
if(sout_flags == 0) cx sout_a[2],sout_d[6];
if(sout_flags == 0) barrier sout_a[2],sout_d[6];


if(sout_flags == 0) barrier sout_a[0], sout_d[0], sout_d[1], sout_d[2], sout_d[3], sout_d[4], sout_d[5], sout_d[6], sout_a[1], sout_a[2];
if(sout_flags == 0) cz sout_a[2],sout_a[0];
if(sout_flags == 0) barrier sout_a[0], sout_d[0], sout_d[1], sout_d[2], sout_d[3], sout_d[4], sout_d[5], sout_d[6], sout_a[1], sout_a[2];



if(sout_flags == 0) barrier sout_a[0],sout_d[2];
if(sout_flags == 0) cz sout_a[0],sout_d[2];
if(sout_flags == 0) barrier sout_a[0],sout_d[2];

if(sout_flags == 0) barrier sout_a[1],sout_d[1];
if(sout_flags == 0) cx sout_a[1],sout_d[1];
if(sout_flags == 0) barrier sout_a[1],sout_d[1];

if(sout_flags == 0) barrier sout_a[2],sout_d[5];
if(sout_flags == 0) cx sout_a[2],sout_d[5];
if(sout_flags == 0) barrier sout_a[2],sout_d[5];


if(sout_flags == 0) h sout_a[0];
if(sout_flags == 0) h sout_a[1];
if(sout_flags == 0) h sout_a[2];



if(sout_flags == 0) measure sout_a[0] -> sout_flags_z[0];
if(sout_flags == 0) measure sout_a[1] -> sout_flag_x[1];
if(sout_flags == 0) measure sout_a[2] -> sout_flag_x[2];

// XOR flags/syndromes
if(sout_flags == 0) sout_flags_z[0] = sout_flags_z[0] ^ sout_last_raw_syn_z[0];
if(sout_flags == 0) sout_flag_x[1] = sout_flag_x[1] ^ sout_last_raw_syn_x[1];
if(sout_flags == 0) sout_flag_x[2] = sout_flag_x[2] ^ sout_last_raw_syn_x[2];

if(sout_flags == 0) sout_flags = sout_flag_x | sout_flags_z;


// Run the 6 non-flagged checks (if non-trivial flags)
// ===================================================
// // X check 1, Z check 2, Z check 3

if(sout_flags != 0) sout_syn_x = 0;
if(sout_flags != 0) sout_syn_z = 0;

if(sout_flags != 0) reset sout_a[0];
if(sout_flags != 0) reset sout_a[1];
if(sout_flags != 0) reset sout_a[2];

if(sout_flags != 0) h sout_a[0];
if(sout_flags != 0) h sout_a[1];
if(sout_flags != 0) h sout_a[2];

if(sout_flags != 0) cx sout_a[0],sout_d[2];
if(sout_flags != 0) cz sout_a[1],sout_d[5];
if(sout_flags != 0) cz sout_a[2],sout_d[6];

if(sout_flags != 0) cx sout_a[0],sout_d[1];
if(sout_flags != 0) cz sout_a[1],sout_d[2];
if(sout_flags != 0) cz sout_a[2],sout_d[5];

if(sout_flags != 0) cx sout_a[0],sout_d[3];
if(sout_flags != 0) cz sout_a[1],sout_d[1];
if(sout_flags != 0) cz sout_a[2],sout_d[2];

if(sout_flags != 0) cx sout_a[0],sout_d[0];
if(sout_flags != 0) cz sout_a[1],sout_d[4];
if(sout_flags != 0) cz sout_a[2],sout_d[3];

if(sout_flags != 0) h sout_a[0];
if(sout_flags != 0) h sout_a[1];
if(sout_flags != 0) h sout_a[2];

if(sout_flags != 0) measure sout_a[0] -> sout_syn_x[0];
if(sout_flags != 0) measure sout_a[1] -> sout_syn_z[1];
if(sout_flags != 0) measure sout_a[2] -> sout_syn_z[2];

// // Z check 1, X check 2, X check 3

if(sout_flags != 0) reset sout_a[0];
if(sout_flags != 0) reset sout_a[1];
if(sout_flags != 0) reset sout_a[2];

if(sout_flags != 0) h sout_a[0];
if(sout_flags != 0) h sout_a[1];
if(sout_flags != 0) h sout_a[2];

if(sout_flags != 0) cz sout_a[0],sout_d[2];
if(sout_flags != 0) cx sout_a[1],sout_d[5];
if(sout_flags != 0) cx sout_a[2],sout_d[6];

if(sout_flags != 0) cz sout_a[0],sout_d[1];
if(sout_flags != 0) cx sout_a[1],sout_d[2];
if(sout_flags != 0) cx sout_a[2],sout_d[5];

if(sout_flags != 0) cz sout_a[0],sout_d[3];
if(sout_flags != 0) cx sout_a[1],sout_d[1];
if(sout_flags != 0) cx sout_a[2],sout_d[2];

if(sout_flags != 0) cz sout_a[0],sout_d[0];
if(sout_flags != 0) cx sout_a[1],sout_d[4];
if(sout_flags != 0) cx sout_a[2],sout_d[3];

if(sout_flags != 0) h sout_a[0];
if(sout_flags != 0) h sout_a[1];
if(sout_flags != 0) h sout_a[2];

if(sout_flags != 0) measure sout_a[0] -> sout_syn_z[0];
if(sout_flags != 0) measure sout_a[1] -> sout_syn_x[1];
if(sout_flags != 0) measure sout_a[2] -> sout_syn_x[2];


// =========================
// BEGIN Run X decoder
// =========================

if(sout_flags!=0) sout_syndromes = sout_syn_x ^ sout_last_raw_syn_x;
if(sout_flags==0) sout_syndromes = 0;

// apply corrections
if(sout_syndromes == 2) sout_c[4] = sout_c[4] ^ 1;
if(sout_syndromes == 4) sout_c[4] = sout_c[4] ^ 1;
if(sout_syndromes == 6) sout_c[4] = sout_c[4] ^ 1;

// alter correction based on flags
// ===============================

// 1&2 (1 -> 2)
// ------------
sout_scratch = 0;
if(sout_flag_x == 1) sout_scratch[0] = 1;
if(sout_syndromes == 2) sout_scratch[1] = 1;

sout_scratch[2] = sout_scratch[0] & sout_scratch[1];
if(sout_scratch[2] == 1) sout_c[4] = sout_c[4] ^ 1;

// 1&4 (1 -> 3)
// ------------
sout_scratch = 0;
if(sout_flag_x == 1) sout_scratch[0] = 1;
if(sout_syndromes == 4) sout_scratch[1] = 1;

sout_scratch[2] = sout_scratch[0] & sout_scratch[1];
if(sout_scratch[2] == 1) sout_c[4] = sout_c[4] ^ 1;


// 6&4 (2,3 -> 3)
// ------------
sout_scratch = 0;
if(sout_flag_x == 6) sout_scratch[0] = 1;
if(sout_syndromes == 4) sout_scratch[1] = 1;

sout_scratch[2] = sout_scratch[0] & sout_scratch[1];
if(sout_scratch[2] == 1) sout_c[4] = sout_c[4] ^ 1;

if(sout_flags!=0) sout_last_raw_syn_x = sout_syn_x;

// =========================
// END Run X decoder
// =========================



// ACTIVE ERROR CORRECTION FOR X SYNDROMES

sout_scratch  = 0;

if(sout_syndromes[0] == 1) sout_scratch = sout_scratch  ^ 1;  // only part that differs for X vs Z syns
if(sout_syndromes[1] == 1) sout_scratch  = sout_scratch  ^ 12;
if(sout_syndromes[2] == 1) sout_scratch  = sout_scratch  ^ 48;

if(sout_c[4]==1) sout_scratch  = sout_scratch  ^ 112;  // logical operator

if(sout_scratch[0] == 1) z sout_d[0];
// if(sout_scratch[1] == 1) z sout_d[1];  // not possible for X stabilizers
if(sout_scratch[2] == 1) z sout_d[2];
if(sout_scratch[3] == 1) z sout_d[3];
if(sout_scratch[4] == 1) z sout_d[4];
if(sout_scratch[5] == 1) z sout_d[5];
if(sout_scratch[6] == 1) z sout_d[6];

sout_c[4] = 0;
// sout_syndromes = 0;
sout_last_raw_syn_x = 0;
// sout_syn_x = 0;
// sout_flag_x = 0;
// sout_flags = 0;



// =========================
// BEGIN Run Z decoder
// =========================

if(sout_flags!=0) sout_syndromes = sout_syn_z ^ sout_last_raw_syn_z;
if(sout_flags==0) sout_syndromes = 0;

// apply corrections
if(sout_syndromes == 2) sout_c[3] = sout_c[3] ^ 1;
if(sout_syndromes == 4) sout_c[3] = sout_c[3] ^ 1;
if(sout_syndromes == 6) sout_c[3] = sout_c[3] ^ 1;

// alter correction based on flags
// ===============================

// 1&2 (1 -> 2)
// ------------
sout_scratch = 0;
if(sout_flags_z == 1) sout_scratch[0] = 1;
if(sout_syndromes == 2) sout_scratch[1] = 1;

sout_scratch[2] = sout_scratch[0] & sout_scratch[1];
if(sout_scratch[2] == 1) sout_c[3] = sout_c[3] ^ 1;

// 1&4 (1 -> 3)
// ------------
sout_scratch = 0;
if(sout_flags_z == 1) sout_scratch[0] = 1;
if(sout_syndromes == 4) sout_scratch[1] = 1;

sout_scratch[2] = sout_scratch[0] & sout_scratch[1];
if(sout_scratch[2] == 1) sout_c[3] = sout_c[3] ^ 1;


// 6&4 (2,3 -> 3)
// ------------
sout_scratch = 0;
if(sout_flags_z == 6) sout_scratch[0] = 1;
if(sout_syndromes == 4) sout_scratch[1] = 1;

sout_scratch[2] = sout_scratch[0] & sout_scratch[1];
if(sout_scratch[2] == 1) sout_c[3] = sout_c[3] ^ 1;

if(sout_flags!=0) sout_last_raw_syn_z = sout_syn_z;

// =========================
// END Run Z decoder
// =========================



// ACTIVE ERROR CORRECTION FOR Z SYNDROMES

sout_scratch  = 0;

if(sout_syndromes[0] == 1) sout_scratch = sout_scratch  ^ 14;  // only part that differs for X vs Z syns
if(sout_syndromes[1] == 1) sout_scratch  = sout_scratch  ^ 12;
if(sout_syndromes[2] == 1) sout_scratch  = sout_scratch  ^ 48;

if(sout_c[3]==1) sout_scratch  = sout_scratch  ^ 112;  // logical operator

// if(sout_scratch[0] == 1) z sout_d[0]; // not possible for X stabilizers
if(sout_scratch[1] == 1) x sout_d[1];
if(sout_scratch[2] == 1) x sout_d[2];
if(sout_scratch[3] == 1) x sout_d[3];
if(sout_scratch[4] == 1) x sout_d[4];
if(sout_scratch[5] == 1) x sout_d[5];
if(sout_scratch[6] == 1) x sout_d[6];

sout_c[3] = 0;
// sout_syndromes = 0;
sout_last_raw_syn_z = 0;
// sout_syn_z = 0;
// sout_flags_z = 0;
// sout_flags = 0;


barrier sin_d[0], sin_d[1], sin_d[2], sin_d[3], sin_d[4], sin_d[5], sin_d[6], sin_a[0];

reset sin_d;
reset sin_a[0];
barrier sin_d, sin_a[0];
h sin_d[0];
h sin_d[4];
h sin_d[6];

cx sin_d[4], sin_d[5];
cx sin_d[0], sin_d[1];
cx sin_d[6], sin_d[3];
cx sin_d[4], sin_d[2];
cx sin_d[6], sin_d[5];
cx sin_d[0], sin_d[3];
cx sin_d[4], sin_d[1];
cx sin_d[3], sin_d[2];

barrier sin_a[0],sin_d[1],sin_d[3],sin_d[5];
//verification step
cx sin_d[5],sin_a[0];
cx sin_d[1],sin_a[0];
cx sin_d[3],sin_a[0];
measure sin_a[0] -> sin_c[0];


if(sin_c[0] == 1) barrier sin_d[0], sin_d[1], sin_d[2], sin_d[3], sin_d[4], sin_d[5], sin_d[6], sin_a[0];

if(sin_c[0] == 1) reset sin_d;
if(sin_c[0] == 1) reset sin_a[0];
if(sin_c[0] == 1) barrier sin_d, sin_a[0];
if(sin_c[0] == 1) h sin_d[0];
if(sin_c[0] == 1) h sin_d[4];
if(sin_c[0] == 1) h sin_d[6];

if(sin_c[0] == 1) cx sin_d[4], sin_d[5];
if(sin_c[0] == 1) cx sin_d[0], sin_d[1];
if(sin_c[0] == 1) cx sin_d[6], sin_d[3];
if(sin_c[0] == 1) cx sin_d[4], sin_d[2];
if(sin_c[0] == 1) cx sin_d[6], sin_d[5];
if(sin_c[0] == 1) cx sin_d[0], sin_d[3];
if(sin_c[0] == 1) cx sin_d[4], sin_d[1];
if(sin_c[0] == 1) cx sin_d[3], sin_d[2];

if(sin_c[0] == 1) barrier sin_a[0],sin_d[1],sin_d[3],sin_d[5];
//verification step
if(sin_c[0] == 1) cx sin_d[5],sin_a[0];
if(sin_c[0] == 1) cx sin_d[1],sin_a[0];
if(sin_c[0] == 1) cx sin_d[3],sin_a[0];
if(sin_c[0] == 1) measure sin_a[0] -> sin_c[0];


if(sin_c[0] == 1) barrier sin_d[0], sin_d[1], sin_d[2], sin_d[3], sin_d[4], sin_d[5], sin_d[6], sin_a[0];

if(sin_c[0] == 1) reset sin_d;
if(sin_c[0] == 1) reset sin_a[0];
if(sin_c[0] == 1) barrier sin_d, sin_a[0];
if(sin_c[0] == 1) h sin_d[0];
if(sin_c[0] == 1) h sin_d[4];
if(sin_c[0] == 1) h sin_d[6];

if(sin_c[0] == 1) cx sin_d[4], sin_d[5];
if(sin_c[0] == 1) cx sin_d[0], sin_d[1];
if(sin_c[0] == 1) cx sin_d[6], sin_d[3];
if(sin_c[0] == 1) cx sin_d[4], sin_d[2];
if(sin_c[0] == 1) cx sin_d[6], sin_d[5];
if(sin_c[0] == 1) cx sin_d[0], sin_d[3];
if(sin_c[0] == 1) cx sin_d[4], sin_d[1];
if(sin_c[0] == 1) cx sin_d[3], sin_d[2];

if(sin_c[0] == 1) barrier sin_a[0],sin_d[1],sin_d[3],sin_d[5];
//verification step
if(sin_c[0] == 1) cx sin_d[5],sin_a[0];
if(sin_c[0] == 1) cx sin_d[1],sin_a[0];
if(sin_c[0] == 1) cx sin_d[3],sin_a[0];
if(sin_c[0] == 1) measure sin_a[0] -> sin_c[0];



sin_flag_x = 0;
sin_flags_z = 0;

// X check 1, Z check 2, Z check 3
// ===============================

reset sin_a[0];
reset sin_a[1];
reset sin_a[2];

h sin_a[0];
h sin_a[1];
h sin_a[2];

cx sin_a[0],sin_d[3];  // 5 -> 4
cz sin_a[1],sin_d[5];  // 6 -> 6
cz sin_a[2],sin_d[2];  // 7 -> 3

barrier sin_a[0],sin_a[1];
cz sin_a[0],sin_a[1];
barrier sin_a[0],sin_a[1];

cx sin_a[0],sin_d[0];  // 1 -> 1
cz sin_a[1],sin_d[4];  // 2 -> 5
cz sin_a[2],sin_d[3];  // 5 -> 4

cx sin_a[0],sin_d[1];  // 3 -> 2
cz sin_a[1],sin_d[2];  // 7 -> 3
cz sin_a[2],sin_d[6];  // 4 -> 7

barrier sin_a[0],sin_a[2];
cz sin_a[0],sin_a[2];
barrier sin_a[0],sin_a[2];

cx sin_a[0],sin_d[2];  // 7 -> 3
cz sin_a[1],sin_d[1];  // 3 -> 2
cz sin_a[2],sin_d[5];  // 6 -> 6

h sin_a[0];
h sin_a[1];
h sin_a[2];

measure sin_a[0] -> sin_flag_x[0];
measure sin_a[1] -> sin_flags_z[1];
measure sin_a[2] -> sin_flags_z[2];

sin_flag_x[0] = sin_flag_x[0] ^ sin_last_raw_syn_x[0];
sin_flags_z[1] = sin_flags_z[1] ^ sin_last_raw_syn_z[1];
sin_flags_z[2] = sin_flags_z[2] ^ sin_last_raw_syn_z[2];

sin_flags = sin_flag_x | sin_flags_z;


// Z check 1, X check 2, X check 3
// ===============================

if(sin_flags == 0) reset sin_a[0];
if(sin_flags == 0) reset sin_a[1];
if(sin_flags == 0) reset sin_a[2];

if(sin_flags == 0) h sin_a[0];
if(sin_flags == 0) h sin_a[1];
if(sin_flags == 0) h sin_a[2];


if(sin_flags == 0) barrier sin_a[0],sin_d[3];
if(sin_flags == 0) cz sin_a[0],sin_d[3];
if(sin_flags == 0) barrier sin_a[0],sin_d[3];

if(sin_flags == 0) barrier sin_a[1],sin_d[5];
if(sin_flags == 0) cx sin_a[1],sin_d[5];
if(sin_flags == 0) barrier sin_a[1],sin_d[5];

if(sin_flags == 0) barrier sin_a[2],sin_d[2];
if(sin_flags == 0) cx sin_a[2],sin_d[2];
if(sin_flags == 0) barrier sin_a[2],sin_d[2];



if(sin_flags == 0) barrier sin_a[0], sin_d[0], sin_d[1], sin_d[2], sin_d[3], sin_d[4], sin_d[5], sin_d[6], sin_a[1], sin_a[2];
if(sin_flags == 0) cz sin_a[1],sin_a[0];
if(sin_flags == 0) barrier sin_a[0], sin_d[0], sin_d[1], sin_d[2], sin_d[3], sin_d[4], sin_d[5], sin_d[6], sin_a[1], sin_a[2];


if(sin_flags == 0) barrier sin_a[0],sin_d[0];
if(sin_flags == 0) cz sin_a[0],sin_d[0];
if(sin_flags == 0) barrier sin_a[0],sin_d[0];

if(sin_flags == 0) barrier sin_a[1],sin_d[4];
if(sin_flags == 0) cx sin_a[1],sin_d[4];
if(sin_flags == 0) barrier sin_a[1],sin_d[4];

if(sin_flags == 0) barrier sin_a[2],sin_d[3];
if(sin_flags == 0) cx sin_a[2],sin_d[3];
if(sin_flags == 0) barrier sin_a[2],sin_d[3];



if(sin_flags == 0) barrier sin_a[0],sin_d[1];
if(sin_flags == 0) cz sin_a[0],sin_d[1];
if(sin_flags == 0) barrier sin_a[0],sin_d[1];

if(sin_flags == 0) barrier sin_a[1],sin_d[2];
if(sin_flags == 0) cx sin_a[1],sin_d[2];
if(sin_flags == 0) barrier sin_a[1],sin_d[2];

if(sin_flags == 0) barrier sin_a[2],sin_d[6];
if(sin_flags == 0) cx sin_a[2],sin_d[6];
if(sin_flags == 0) barrier sin_a[2],sin_d[6];


if(sin_flags == 0) barrier sin_a[0], sin_d[0], sin_d[1], sin_d[2], sin_d[3], sin_d[4], sin_d[5], sin_d[6], sin_a[1], sin_a[2];
if(sin_flags == 0) cz sin_a[2],sin_a[0];
if(sin_flags == 0) barrier sin_a[0], sin_d[0], sin_d[1], sin_d[2], sin_d[3], sin_d[4], sin_d[5], sin_d[6], sin_a[1], sin_a[2];



if(sin_flags == 0) barrier sin_a[0],sin_d[2];
if(sin_flags == 0) cz sin_a[0],sin_d[2];
if(sin_flags == 0) barrier sin_a[0],sin_d[2];

if(sin_flags == 0) barrier sin_a[1],sin_d[1];
if(sin_flags == 0) cx sin_a[1],sin_d[1];
if(sin_flags == 0) barrier sin_a[1],sin_d[1];

if(sin_flags == 0) barrier sin_a[2],sin_d[5];
if(sin_flags == 0) cx sin_a[2],sin_d[5];
if(sin_flags == 0) barrier sin_a[2],sin_d[5];


if(sin_flags == 0) h sin_a[0];
if(sin_flags == 0) h sin_a[1];
if(sin_flags == 0) h sin_a[2];



if(sin_flags == 0) measure sin_a[0] -> sin_flags_z[0];
if(sin_flags == 0) measure sin_a[1] -> sin_flag_x[1];
if(sin_flags == 0) measure sin_a[2] -> sin_flag_x[2];

// XOR flags/syndromes
if(sin_flags == 0) sin_flags_z[0] = sin_flags_z[0] ^ sin_last_raw_syn_z[0];
if(sin_flags == 0) sin_flag_x[1] = sin_flag_x[1] ^ sin_last_raw_syn_x[1];
if(sin_flags == 0) sin_flag_x[2] = sin_flag_x[2] ^ sin_last_raw_syn_x[2];

if(sin_flags == 0) sin_flags = sin_flag_x | sin_flags_z;


// Run the 6 non-flagged checks (if non-trivial flags)
// ===================================================
// // X check 1, Z check 2, Z check 3

if(sin_flags != 0) sin_syn_x = 0;
if(sin_flags != 0) sin_syn_z = 0;

if(sin_flags != 0) reset sin_a[0];
if(sin_flags != 0) reset sin_a[1];
if(sin_flags != 0) reset sin_a[2];

if(sin_flags != 0) h sin_a[0];
if(sin_flags != 0) h sin_a[1];
if(sin_flags != 0) h sin_a[2];

if(sin_flags != 0) cx sin_a[0],sin_d[2];
if(sin_flags != 0) cz sin_a[1],sin_d[5];
if(sin_flags != 0) cz sin_a[2],sin_d[6];

if(sin_flags != 0) cx sin_a[0],sin_d[1];
if(sin_flags != 0) cz sin_a[1],sin_d[2];
if(sin_flags != 0) cz sin_a[2],sin_d[5];

if(sin_flags != 0) cx sin_a[0],sin_d[3];
if(sin_flags != 0) cz sin_a[1],sin_d[1];
if(sin_flags != 0) cz sin_a[2],sin_d[2];

if(sin_flags != 0) cx sin_a[0],sin_d[0];
if(sin_flags != 0) cz sin_a[1],sin_d[4];
if(sin_flags != 0) cz sin_a[2],sin_d[3];

if(sin_flags != 0) h sin_a[0];
if(sin_flags != 0) h sin_a[1];
if(sin_flags != 0) h sin_a[2];

if(sin_flags != 0) measure sin_a[0] -> sin_syn_x[0];
if(sin_flags != 0) measure sin_a[1] -> sin_syn_z[1];
if(sin_flags != 0) measure sin_a[2] -> sin_syn_z[2];

// // Z check 1, X check 2, X check 3

if(sin_flags != 0) reset sin_a[0];
if(sin_flags != 0) reset sin_a[1];
if(sin_flags != 0) reset sin_a[2];

if(sin_flags != 0) h sin_a[0];
if(sin_flags != 0) h sin_a[1];
if(sin_flags != 0) h sin_a[2];

if(sin_flags != 0) cz sin_a[0],sin_d[2];
if(sin_flags != 0) cx sin_a[1],sin_d[5];
if(sin_flags != 0) cx sin_a[2],sin_d[6];

if(sin_flags != 0) cz sin_a[0],sin_d[1];
if(sin_flags != 0) cx sin_a[1],sin_d[2];
if(sin_flags != 0) cx sin_a[2],sin_d[5];

if(sin_flags != 0) cz sin_a[0],sin_d[3];
if(sin_flags != 0) cx sin_a[1],sin_d[1];
if(sin_flags != 0) cx sin_a[2],sin_d[2];

if(sin_flags != 0) cz sin_a[0],sin_d[0];
if(sin_flags != 0) cx sin_a[1],sin_d[4];
if(sin_flags != 0) cx sin_a[2],sin_d[3];

if(sin_flags != 0) h sin_a[0];
if(sin_flags != 0) h sin_a[1];
if(sin_flags != 0) h sin_a[2];

if(sin_flags != 0) measure sin_a[0] -> sin_syn_z[0];
if(sin_flags != 0) measure sin_a[1] -> sin_syn_x[1];
if(sin_flags != 0) measure sin_a[2] -> sin_syn_x[2];


// =========================
// BEGIN Run X decoder
// =========================

if(sin_flags!=0) sin_syndromes = sin_syn_x ^ sin_last_raw_syn_x;
if(sin_flags==0) sin_syndromes = 0;

// apply corrections
if(sin_syndromes == 2) sin_c[4] = sin_c[4] ^ 1;
if(sin_syndromes == 4) sin_c[4] = sin_c[4] ^ 1;
if(sin_syndromes == 6) sin_c[4] = sin_c[4] ^ 1;

// alter correction based on flags
// ===============================

// 1&2 (1 -> 2)
// ------------
sin_scratch = 0;
if(sin_flag_x == 1) sin_scratch[0] = 1;
if(sin_syndromes == 2) sin_scratch[1] = 1;

sin_scratch[2] = sin_scratch[0] & sin_scratch[1];
if(sin_scratch[2] == 1) sin_c[4] = sin_c[4] ^ 1;

// 1&4 (1 -> 3)
// ------------
sin_scratch = 0;
if(sin_flag_x == 1) sin_scratch[0] = 1;
if(sin_syndromes == 4) sin_scratch[1] = 1;

sin_scratch[2] = sin_scratch[0] & sin_scratch[1];
if(sin_scratch[2] == 1) sin_c[4] = sin_c[4] ^ 1;


// 6&4 (2,3 -> 3)
// ------------
sin_scratch = 0;
if(sin_flag_x == 6) sin_scratch[0] = 1;
if(sin_syndromes == 4) sin_scratch[1] = 1;

sin_scratch[2] = sin_scratch[0] & sin_scratch[1];
if(sin_scratch[2] == 1) sin_c[4] = sin_c[4] ^ 1;

if(sin_flags!=0) sin_last_raw_syn_x = sin_syn_x;

// =========================
// END Run X decoder
// =========================



// ACTIVE ERROR CORRECTION FOR X SYNDROMES

sin_scratch  = 0;

if(sin_syndromes[0] == 1) sin_scratch = sin_scratch  ^ 1;  // only part that differs for X vs Z syns
if(sin_syndromes[1] == 1) sin_scratch  = sin_scratch  ^ 12;
if(sin_syndromes[2] == 1) sin_scratch  = sin_scratch  ^ 48;

if(sin_c[4]==1) sin_scratch  = sin_scratch  ^ 112;  // logical operator

if(sin_scratch[0] == 1) z sin_d[0];
// if(sin_scratch[1] == 1) z sin_d[1];  // not possible for X stabilizers
if(sin_scratch[2] == 1) z sin_d[2];
if(sin_scratch[3] == 1) z sin_d[3];
if(sin_scratch[4] == 1) z sin_d[4];
if(sin_scratch[5] == 1) z sin_d[5];
if(sin_scratch[6] == 1) z sin_d[6];

sin_c[4] = 0;
// sin_syndromes = 0;
sin_last_raw_syn_x = 0;
// sin_syn_x = 0;
// sin_flag_x = 0;
// sin_flags = 0;



// =========================
// BEGIN Run Z decoder
// =========================

if(sin_flags!=0) sin_syndromes = sin_syn_z ^ sin_last_raw_syn_z;
if(sin_flags==0) sin_syndromes = 0;

// apply corrections
if(sin_syndromes == 2) sin_c[3] = sin_c[3] ^ 1;
if(sin_syndromes == 4) sin_c[3] = sin_c[3] ^ 1;
if(sin_syndromes == 6) sin_c[3] = sin_c[3] ^ 1;

// alter correction based on flags
// ===============================

// 1&2 (1 -> 2)
// ------------
sin_scratch = 0;
if(sin_flags_z == 1) sin_scratch[0] = 1;
if(sin_syndromes == 2) sin_scratch[1] = 1;

sin_scratch[2] = sin_scratch[0] & sin_scratch[1];
if(sin_scratch[2] == 1) sin_c[3] = sin_c[3] ^ 1;

// 1&4 (1 -> 3)
// ------------
sin_scratch = 0;
if(sin_flags_z == 1) sin_scratch[0] = 1;
if(sin_syndromes == 4) sin_scratch[1] = 1;

sin_scratch[2] = sin_scratch[0] & sin_scratch[1];
if(sin_scratch[2] == 1) sin_c[3] = sin_c[3] ^ 1;


// 6&4 (2,3 -> 3)
// ------------
sin_scratch = 0;
if(sin_flags_z == 6) sin_scratch[0] = 1;
if(sin_syndromes == 4) sin_scratch[1] = 1;

sin_scratch[2] = sin_scratch[0] & sin_scratch[1];
if(sin_scratch[2] == 1) sin_c[3] = sin_c[3] ^ 1;

if(sin_flags!=0) sin_last_raw_syn_z = sin_syn_z;

// =========================
// END Run Z decoder
// =========================



// ACTIVE ERROR CORRECTION FOR Z SYNDROMES

sin_scratch  = 0;

if(sin_syndromes[0] == 1) sin_scratch = sin_scratch  ^ 14;  // only part that differs for X vs Z syns
if(sin_syndromes[1] == 1) sin_scratch  = sin_scratch  ^ 12;
if(sin_syndromes[2] == 1) sin_scratch  = sin_scratch  ^ 48;

if(sin_c[3]==1) sin_scratch  = sin_scratch  ^ 112;  // logical operator

// if(sin_scratch[0] == 1) z sin_d[0]; // not possible for X stabilizers
if(sin_scratch[1] == 1) x sin_d[1];
if(sin_scratch[2] == 1) x sin_d[2];
if(sin_scratch[3] == 1) x sin_d[3];
if(sin_scratch[4] == 1) x sin_d[4];
if(sin_scratch[5] == 1) x sin_d[5];
if(sin_scratch[6] == 1) x sin_d[6];

sin_c[3] = 0;
// sin_syndromes = 0;
sin_last_raw_syn_z = 0;
// sin_syn_z = 0;
// sin_flags_z = 0;
// sin_flags = 0;

// Transversal Logical CX
barrier sin_d, smid_d;
cx sin_d[0], smid_d[0];
cx sin_d[1], smid_d[1];
cx sin_d[2], smid_d[2];
cx sin_d[3], smid_d[3];
cx sin_d[4], smid_d[4];
cx sin_d[5], smid_d[5];
cx sin_d[6], smid_d[6];
barrier sin_d, smid_d;
// Logical H
h sin_d;
// Destructive logical Z measurement

barrier sin_d;

measure sin_d[0] -> sin_raw_meas[0];
measure sin_d[1] -> sin_raw_meas[1];
measure sin_d[2] -> sin_raw_meas[2];
measure sin_d[3] -> sin_raw_meas[3];
measure sin_d[4] -> sin_raw_meas[4];
measure sin_d[5] -> sin_raw_meas[5];
measure sin_d[6] -> sin_raw_meas[6];

// determine raw logical output
// ============================
sin_c[1] = sin_raw_meas[4] ^ sin_raw_meas[5] ^ sin_raw_meas[6];



// =================== //
// PROCESS MEASUREMENT //
// =================== //

// Determine correction to get logical output
// ==========================================
sin_syn_meas[0] = sin_raw_meas[0] ^ sin_raw_meas[1] ^ sin_raw_meas[2] ^ sin_raw_meas[3];
sin_syn_meas[1] = sin_raw_meas[1] ^ sin_raw_meas[2] ^ sin_raw_meas[4] ^ sin_raw_meas[5];
sin_syn_meas[2] = sin_raw_meas[2] ^ sin_raw_meas[3] ^ sin_raw_meas[5] ^ sin_raw_meas[6];

// XOR syndromes
sin_syn_meas = sin_syn_meas ^ sin_last_raw_syn_z;

// Correct logical output based on measured out syndromes
sin_c[2] = sin_c[1];
if(sin_syn_meas == 2) sin_c[2] = sin_c[2] ^ 1;
if(sin_syn_meas == 4) sin_c[2] = sin_c[2] ^ 1;
if(sin_syn_meas == 6) sin_c[2] = sin_c[2] ^ 1;

// Apply Pauli frame update (flip the logical output)
// Update for logical Z out
sin_c[2] = sin_c[2] ^ sin_c[3];
m_bell[0] = sin_c[2];
// Destructive logical Z measurement

barrier smid_d;

measure smid_d[0] -> smid_raw_meas[0];
measure smid_d[1] -> smid_raw_meas[1];
measure smid_d[2] -> smid_raw_meas[2];
measure smid_d[3] -> smid_raw_meas[3];
measure smid_d[4] -> smid_raw_meas[4];
measure smid_d[5] -> smid_raw_meas[5];
measure smid_d[6] -> smid_raw_meas[6];

// determine raw logical output
// ============================
smid_c[1] = smid_raw_meas[4] ^ smid_raw_meas[5] ^ smid_raw_meas[6];



// =================== //
// PROCESS MEASUREMENT //
// =================== //

// Determine correction to get logical output
// ==========================================
smid_syn_meas[0] = smid_raw_meas[0] ^ smid_raw_meas[1] ^ smid_raw_meas[2] ^ smid_raw_meas[3];
smid_syn_meas[1] = smid_raw_meas[1] ^ smid_raw_meas[2] ^ smid_raw_meas[4] ^ smid_raw_meas[5];
smid_syn_meas[2] = smid_raw_meas[2] ^ smid_raw_meas[3] ^ smid_raw_meas[5] ^ smid_raw_meas[6];

// XOR syndromes
smid_syn_meas = smid_syn_meas ^ smid_last_raw_syn_z;

// Correct logical output based on measured out syndromes
smid_c[2] = smid_c[1];
if(smid_syn_meas == 2) smid_c[2] = smid_c[2] ^ 1;
if(smid_syn_meas == 4) smid_c[2] = smid_c[2] ^ 1;
if(smid_syn_meas == 6) smid_c[2] = smid_c[2] ^ 1;

// Apply Pauli frame update (flip the logical output)
// Update for logical Z out
smid_c[2] = smid_c[2] ^ smid_c[3];
m_bell[1] = smid_c[2];
// Logical X
if(m_bell[1] == 0) x sout_d[4];
if(m_bell[1] == 0) x sout_d[5];
if(m_bell[1] == 0) x sout_d[6];
// Logical Z
if(m_bell[0] == 0) z sout_d[4];
if(m_bell[0] == 0) z sout_d[5];
if(m_bell[0] == 0) z sout_d[6];
// Destructive logical Z measurement

barrier sout_d;

measure sout_d[0] -> sout_raw_meas[0];
measure sout_d[1] -> sout_raw_meas[1];
measure sout_d[2] -> sout_raw_meas[2];
measure sout_d[3] -> sout_raw_meas[3];
measure sout_d[4] -> sout_raw_meas[4];
measure sout_d[5] -> sout_raw_meas[5];
measure sout_d[6] -> sout_raw_meas[6];

// determine raw logical output
// ============================
sout_c[1] = sout_raw_meas[4] ^ sout_raw_meas[5] ^ sout_raw_meas[6];



// =================== //
// PROCESS MEASUREMENT //
// =================== //

// Determine correction to get logical output
// ==========================================
sout_syn_meas[0] = sout_raw_meas[0] ^ sout_raw_meas[1] ^ sout_raw_meas[2] ^ sout_raw_meas[3];
sout_syn_meas[1] = sout_raw_meas[1] ^ sout_raw_meas[2] ^ sout_raw_meas[4] ^ sout_raw_meas[5];
sout_syn_meas[2] = sout_raw_meas[2] ^ sout_raw_meas[3] ^ sout_raw_meas[5] ^ sout_raw_meas[6];

// XOR syndromes
sout_syn_meas = sout_syn_meas ^ sout_last_raw_syn_z;

// Correct logical output based on measured out syndromes
sout_c[2] = sout_c[1];
if(sout_syn_meas == 2) sout_c[2] = sout_c[2] ^ 1;
if(sout_syn_meas == 4) sout_c[2] = sout_c[2] ^ 1;
if(sout_syn_meas == 6) sout_c[2] = sout_c[2] ^ 1;

// Apply Pauli frame update (flip the logical output)
// Update for logical Z out
sout_c[2] = sout_c[2] ^ sout_c[3];
m_out[0] = sout_c[2];

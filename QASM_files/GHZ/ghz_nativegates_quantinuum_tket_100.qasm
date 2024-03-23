// Benchmark was created by MQT Bench on 2023-06-29
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: v1.0.0
// TKET version: 1.16.0
// Used Gate Set: ['rz', 'ry', 'rx', 'rzz', 'measure']

OPENQASM 2.0;
include "qelib1.inc";

qreg q[100];
creg meas[100];
rz(1.0*pi) q[0];
rz(1.0*pi) q[1];
rz(1.0*pi) q[2];
rz(1.0*pi) q[3];
rz(1.0*pi) q[4];
rz(1.0*pi) q[5];
rz(1.0*pi) q[6];
rz(1.0*pi) q[7];
rz(1.0*pi) q[8];
rz(1.0*pi) q[9];
rz(1.0*pi) q[10];
rz(1.0*pi) q[11];
rz(1.0*pi) q[12];
rz(1.0*pi) q[13];
rz(1.0*pi) q[14];
rz(1.0*pi) q[15];
rz(1.0*pi) q[16];
rz(1.0*pi) q[17];
rz(1.0*pi) q[18];
rz(1.0*pi) q[19];
rz(1.0*pi) q[20];
rz(1.0*pi) q[21];
rz(1.0*pi) q[22];
rz(1.0*pi) q[23];
rz(1.0*pi) q[24];
rz(1.0*pi) q[25];
rz(1.0*pi) q[26];
rz(1.0*pi) q[27];
rz(1.0*pi) q[28];
rz(1.0*pi) q[29];
rz(1.0*pi) q[30];
rz(1.0*pi) q[31];
rz(1.0*pi) q[32];
rz(1.0*pi) q[33];
rz(1.0*pi) q[34];
rz(1.0*pi) q[35];
rz(1.0*pi) q[36];
rz(1.0*pi) q[37];
rz(1.0*pi) q[38];
rz(1.0*pi) q[39];
rz(1.0*pi) q[40];
rz(1.0*pi) q[41];
rz(1.0*pi) q[42];
rz(1.0*pi) q[43];
rz(1.0*pi) q[44];
rz(1.0*pi) q[45];
rz(1.0*pi) q[46];
rz(1.0*pi) q[47];
rz(1.0*pi) q[48];
rz(1.0*pi) q[49];
rz(1.0*pi) q[50];
rz(1.0*pi) q[51];
rz(1.0*pi) q[52];
rz(1.0*pi) q[53];
rz(1.0*pi) q[54];
rz(1.0*pi) q[55];
rz(1.0*pi) q[56];
rz(1.0*pi) q[57];
rz(1.0*pi) q[58];
rz(1.0*pi) q[59];
rz(1.0*pi) q[60];
rz(1.0*pi) q[61];
rz(1.0*pi) q[62];
rz(1.0*pi) q[63];
rz(1.0*pi) q[64];
rz(1.0*pi) q[65];
rz(1.0*pi) q[66];
rz(1.0*pi) q[67];
rz(1.0*pi) q[68];
rz(1.0*pi) q[69];
rz(1.0*pi) q[70];
rz(1.0*pi) q[71];
rz(1.0*pi) q[72];
rz(1.0*pi) q[73];
rz(1.0*pi) q[74];
rz(1.0*pi) q[75];
rz(1.0*pi) q[76];
rz(1.0*pi) q[77];
rz(1.0*pi) q[78];
rz(1.0*pi) q[79];
rz(1.0*pi) q[80];
rz(1.0*pi) q[81];
rz(1.0*pi) q[82];
rz(1.0*pi) q[83];
rz(1.0*pi) q[84];
rz(1.0*pi) q[85];
rz(1.0*pi) q[86];
rz(1.0*pi) q[87];
rz(1.0*pi) q[88];
rz(1.0*pi) q[89];
rz(1.0*pi) q[90];
rz(1.0*pi) q[91];
rz(1.0*pi) q[92];
rz(1.0*pi) q[93];
rz(1.0*pi) q[94];
rz(1.0*pi) q[95];
rz(1.0*pi) q[96];
rz(1.0*pi) q[97];
rz(1.0*pi) q[98];
rx(0.5*pi) q[99];
rx(3.5*pi) q[0];
rx(3.5*pi) q[1];
rx(3.5*pi) q[2];
rx(3.5*pi) q[3];
rx(3.5*pi) q[4];
rx(3.5*pi) q[5];
rx(3.5*pi) q[6];
rx(3.5*pi) q[7];
rx(3.5*pi) q[8];
rx(3.5*pi) q[9];
rx(3.5*pi) q[10];
rx(3.5*pi) q[11];
rx(3.5*pi) q[12];
rx(3.5*pi) q[13];
rx(3.5*pi) q[14];
rx(3.5*pi) q[15];
rx(3.5*pi) q[16];
rx(3.5*pi) q[17];
rx(3.5*pi) q[18];
rx(3.5*pi) q[19];
rx(3.5*pi) q[20];
rx(3.5*pi) q[21];
rx(3.5*pi) q[22];
rx(3.5*pi) q[23];
rx(3.5*pi) q[24];
rx(3.5*pi) q[25];
rx(3.5*pi) q[26];
rx(3.5*pi) q[27];
rx(3.5*pi) q[28];
rx(3.5*pi) q[29];
rx(3.5*pi) q[30];
rx(3.5*pi) q[31];
rx(3.5*pi) q[32];
rx(3.5*pi) q[33];
rx(3.5*pi) q[34];
rx(3.5*pi) q[35];
rx(3.5*pi) q[36];
rx(3.5*pi) q[37];
rx(3.5*pi) q[38];
rx(3.5*pi) q[39];
rx(3.5*pi) q[40];
rx(3.5*pi) q[41];
rx(3.5*pi) q[42];
rx(3.5*pi) q[43];
rx(3.5*pi) q[44];
rx(3.5*pi) q[45];
rx(3.5*pi) q[46];
rx(3.5*pi) q[47];
rx(3.5*pi) q[48];
rx(3.5*pi) q[49];
rx(3.5*pi) q[50];
rx(3.5*pi) q[51];
rx(3.5*pi) q[52];
rx(3.5*pi) q[53];
rx(3.5*pi) q[54];
rx(3.5*pi) q[55];
rx(3.5*pi) q[56];
rx(3.5*pi) q[57];
rx(3.5*pi) q[58];
rx(3.5*pi) q[59];
rx(3.5*pi) q[60];
rx(3.5*pi) q[61];
rx(3.5*pi) q[62];
rx(3.5*pi) q[63];
rx(3.5*pi) q[64];
rx(3.5*pi) q[65];
rx(3.5*pi) q[66];
rx(3.5*pi) q[67];
rx(3.5*pi) q[68];
rx(3.5*pi) q[69];
rx(3.5*pi) q[70];
rx(3.5*pi) q[71];
rx(3.5*pi) q[72];
rx(3.5*pi) q[73];
rx(3.5*pi) q[74];
rx(3.5*pi) q[75];
rx(3.5*pi) q[76];
rx(3.5*pi) q[77];
rx(3.5*pi) q[78];
rx(3.5*pi) q[79];
rx(3.5*pi) q[80];
rx(3.5*pi) q[81];
rx(3.5*pi) q[82];
rx(3.5*pi) q[83];
rx(3.5*pi) q[84];
rx(3.5*pi) q[85];
rx(3.5*pi) q[86];
rx(3.5*pi) q[87];
rx(3.5*pi) q[88];
rx(3.5*pi) q[89];
rx(3.5*pi) q[90];
rx(3.5*pi) q[91];
rx(3.5*pi) q[92];
rx(3.5*pi) q[93];
rx(3.5*pi) q[94];
rx(3.5*pi) q[95];
rx(3.5*pi) q[96];
rx(3.5*pi) q[97];
rx(3.5*pi) q[98];
rz(0.5*pi) q[99];
rz(0.5*pi) q[0];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rz(0.5*pi) q[3];
rz(0.5*pi) q[4];
rz(0.5*pi) q[5];
rz(0.5*pi) q[6];
rz(0.5*pi) q[7];
rz(0.5*pi) q[8];
rz(0.5*pi) q[9];
rz(0.5*pi) q[10];
rz(0.5*pi) q[11];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rz(0.5*pi) q[14];
rz(0.5*pi) q[15];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rz(0.5*pi) q[18];
rz(0.5*pi) q[19];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rz(0.5*pi) q[25];
rz(0.5*pi) q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[28];
rz(0.5*pi) q[29];
rz(0.5*pi) q[30];
rz(0.5*pi) q[31];
rz(0.5*pi) q[32];
rz(0.5*pi) q[33];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rz(0.5*pi) q[36];
rz(0.5*pi) q[37];
rz(0.5*pi) q[38];
rz(0.5*pi) q[39];
rz(0.5*pi) q[40];
rz(0.5*pi) q[41];
rz(0.5*pi) q[42];
rz(0.5*pi) q[43];
rz(0.5*pi) q[44];
rz(0.5*pi) q[45];
rz(0.5*pi) q[46];
rz(0.5*pi) q[47];
rz(0.5*pi) q[48];
rz(0.5*pi) q[49];
rz(0.5*pi) q[50];
rz(0.5*pi) q[51];
rz(0.5*pi) q[52];
rz(0.5*pi) q[53];
rz(0.5*pi) q[54];
rz(0.5*pi) q[55];
rz(0.5*pi) q[56];
rz(0.5*pi) q[57];
rz(0.5*pi) q[58];
rz(0.5*pi) q[59];
rz(0.5*pi) q[60];
rz(0.5*pi) q[61];
rz(0.5*pi) q[62];
rz(0.5*pi) q[63];
rz(0.5*pi) q[64];
rz(0.5*pi) q[65];
rz(0.5*pi) q[66];
rz(0.5*pi) q[67];
rz(0.5*pi) q[68];
rz(0.5*pi) q[69];
rz(0.5*pi) q[70];
rz(0.5*pi) q[71];
rz(0.5*pi) q[72];
rz(0.5*pi) q[73];
rz(0.5*pi) q[74];
rz(0.5*pi) q[75];
rz(0.5*pi) q[76];
rz(0.5*pi) q[77];
rz(0.5*pi) q[78];
rz(0.5*pi) q[79];
rz(0.5*pi) q[80];
rz(0.5*pi) q[81];
rz(0.5*pi) q[82];
rz(0.5*pi) q[83];
rz(0.5*pi) q[84];
rz(0.5*pi) q[85];
rz(0.5*pi) q[86];
rz(0.5*pi) q[87];
rz(0.5*pi) q[88];
rz(0.5*pi) q[89];
rz(0.5*pi) q[90];
rz(0.5*pi) q[91];
rz(0.5*pi) q[92];
rz(0.5*pi) q[93];
rz(0.5*pi) q[94];
rz(0.5*pi) q[95];
rz(0.5*pi) q[96];
rz(0.5*pi) q[97];
rz(0.5*pi) q[98];
rx(0.5*pi) q[99];
rx(0.5*pi) q[0];
rx(0.5*pi) q[1];
rx(0.5*pi) q[2];
rx(0.5*pi) q[3];
rx(0.5*pi) q[4];
rx(0.5*pi) q[5];
rx(0.5*pi) q[6];
rx(0.5*pi) q[7];
rx(0.5*pi) q[8];
rx(0.5*pi) q[9];
rx(0.5*pi) q[10];
rx(0.5*pi) q[11];
rx(0.5*pi) q[12];
rx(0.5*pi) q[13];
rx(0.5*pi) q[14];
rx(0.5*pi) q[15];
rx(0.5*pi) q[16];
rx(0.5*pi) q[17];
rx(0.5*pi) q[18];
rx(0.5*pi) q[19];
rx(0.5*pi) q[20];
rx(0.5*pi) q[21];
rx(0.5*pi) q[22];
rx(0.5*pi) q[23];
rx(0.5*pi) q[24];
rx(0.5*pi) q[25];
rx(0.5*pi) q[26];
rx(0.5*pi) q[27];
rx(0.5*pi) q[28];
rx(0.5*pi) q[29];
rx(0.5*pi) q[30];
rx(0.5*pi) q[31];
rx(0.5*pi) q[32];
rx(0.5*pi) q[33];
rx(0.5*pi) q[34];
rx(0.5*pi) q[35];
rx(0.5*pi) q[36];
rx(0.5*pi) q[37];
rx(0.5*pi) q[38];
rx(0.5*pi) q[39];
rx(0.5*pi) q[40];
rx(0.5*pi) q[41];
rx(0.5*pi) q[42];
rx(0.5*pi) q[43];
rx(0.5*pi) q[44];
rx(0.5*pi) q[45];
rx(0.5*pi) q[46];
rx(0.5*pi) q[47];
rx(0.5*pi) q[48];
rx(0.5*pi) q[49];
rx(0.5*pi) q[50];
rx(0.5*pi) q[51];
rx(0.5*pi) q[52];
rx(0.5*pi) q[53];
rx(0.5*pi) q[54];
rx(0.5*pi) q[55];
rx(0.5*pi) q[56];
rx(0.5*pi) q[57];
rx(0.5*pi) q[58];
rx(0.5*pi) q[59];
rx(0.5*pi) q[60];
rx(0.5*pi) q[61];
rx(0.5*pi) q[62];
rx(0.5*pi) q[63];
rx(0.5*pi) q[64];
rx(0.5*pi) q[65];
rx(0.5*pi) q[66];
rx(0.5*pi) q[67];
rx(0.5*pi) q[68];
rx(0.5*pi) q[69];
rx(0.5*pi) q[70];
rx(0.5*pi) q[71];
rx(0.5*pi) q[72];
rx(0.5*pi) q[73];
rx(0.5*pi) q[74];
rx(0.5*pi) q[75];
rx(0.5*pi) q[76];
rx(0.5*pi) q[77];
rx(0.5*pi) q[78];
rx(0.5*pi) q[79];
rx(0.5*pi) q[80];
rx(0.5*pi) q[81];
rx(0.5*pi) q[82];
rx(0.5*pi) q[83];
rx(0.5*pi) q[84];
rx(0.5*pi) q[85];
rx(0.5*pi) q[86];
rx(0.5*pi) q[87];
rx(0.5*pi) q[88];
rx(0.5*pi) q[89];
rx(0.5*pi) q[90];
rx(0.5*pi) q[91];
rx(0.5*pi) q[92];
rx(0.5*pi) q[93];
rx(0.5*pi) q[94];
rx(0.5*pi) q[95];
rx(0.5*pi) q[96];
rx(0.5*pi) q[97];
rx(0.5*pi) q[98];
rz(0.5*pi) q[99];
rz(0.5*pi) q[0];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rz(0.5*pi) q[3];
rz(0.5*pi) q[4];
rz(0.5*pi) q[5];
rz(0.5*pi) q[6];
rz(0.5*pi) q[7];
rz(0.5*pi) q[8];
rz(0.5*pi) q[9];
rz(0.5*pi) q[10];
rz(0.5*pi) q[11];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rz(0.5*pi) q[14];
rz(0.5*pi) q[15];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rz(0.5*pi) q[18];
rz(0.5*pi) q[19];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rz(0.5*pi) q[25];
rz(0.5*pi) q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[28];
rz(0.5*pi) q[29];
rz(0.5*pi) q[30];
rz(0.5*pi) q[31];
rz(0.5*pi) q[32];
rz(0.5*pi) q[33];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rz(0.5*pi) q[36];
rz(0.5*pi) q[37];
rz(0.5*pi) q[38];
rz(0.5*pi) q[39];
rz(0.5*pi) q[40];
rz(0.5*pi) q[41];
rz(0.5*pi) q[42];
rz(0.5*pi) q[43];
rz(0.5*pi) q[44];
rz(0.5*pi) q[45];
rz(0.5*pi) q[46];
rz(0.5*pi) q[47];
rz(0.5*pi) q[48];
rz(0.5*pi) q[49];
rz(0.5*pi) q[50];
rz(0.5*pi) q[51];
rz(0.5*pi) q[52];
rz(0.5*pi) q[53];
rz(0.5*pi) q[54];
rz(0.5*pi) q[55];
rz(0.5*pi) q[56];
rz(0.5*pi) q[57];
rz(0.5*pi) q[58];
rz(0.5*pi) q[59];
rz(0.5*pi) q[60];
rz(0.5*pi) q[61];
rz(0.5*pi) q[62];
rz(0.5*pi) q[63];
rz(0.5*pi) q[64];
rz(0.5*pi) q[65];
rz(0.5*pi) q[66];
rz(0.5*pi) q[67];
rz(0.5*pi) q[68];
rz(0.5*pi) q[69];
rz(0.5*pi) q[70];
rz(0.5*pi) q[71];
rz(0.5*pi) q[72];
rz(0.5*pi) q[73];
rz(0.5*pi) q[74];
rz(0.5*pi) q[75];
rz(0.5*pi) q[76];
rz(0.5*pi) q[77];
rz(0.5*pi) q[78];
rz(0.5*pi) q[79];
rz(0.5*pi) q[80];
rz(0.5*pi) q[81];
rz(0.5*pi) q[82];
rz(0.5*pi) q[83];
rz(0.5*pi) q[84];
rz(0.5*pi) q[85];
rz(0.5*pi) q[86];
rz(0.5*pi) q[87];
rz(0.5*pi) q[88];
rz(0.5*pi) q[89];
rz(0.5*pi) q[90];
rz(0.5*pi) q[91];
rz(0.5*pi) q[92];
rz(0.5*pi) q[93];
rz(0.5*pi) q[94];
rz(0.5*pi) q[95];
rz(0.5*pi) q[96];
rz(0.5*pi) q[97];
rz(0.5*pi) q[98];
rzz(0.5*pi) q[99],q[98];
rz(0.5*pi) q[98];
rz(0.5*pi) q[99];
rx(0.5*pi) q[98];
rx(0.5*pi) q[99];
rz(0.5*pi) q[98];
rz(0.5*pi) q[99];
rx(0.5*pi) q[98];
rz(0.5*pi) q[99];
rz(0.5*pi) q[98];
rx(0.5*pi) q[99];
rz(0.5*pi) q[98];
rz(0.5*pi) q[99];
rx(0.5*pi) q[98];
rz(0.5*pi) q[98];
rzz(0.5*pi) q[98],q[97];
rz(0.5*pi) q[97];
rz(0.5*pi) q[98];
rx(0.5*pi) q[97];
rx(0.5*pi) q[98];
rz(0.5*pi) q[97];
rz(0.5*pi) q[98];
rx(0.5*pi) q[97];
rz(0.5*pi) q[98];
rz(0.5*pi) q[97];
rx(0.5*pi) q[98];
rz(0.5*pi) q[97];
rz(0.5*pi) q[98];
rx(0.5*pi) q[97];
rz(0.5*pi) q[97];
rzz(0.5*pi) q[97],q[96];
rz(0.5*pi) q[96];
rz(0.5*pi) q[97];
rx(0.5*pi) q[96];
rx(0.5*pi) q[97];
rz(0.5*pi) q[96];
rz(0.5*pi) q[97];
rx(0.5*pi) q[96];
rz(0.5*pi) q[97];
rz(0.5*pi) q[96];
rx(0.5*pi) q[97];
rz(0.5*pi) q[96];
rz(0.5*pi) q[97];
rx(0.5*pi) q[96];
rz(0.5*pi) q[96];
rzz(0.5*pi) q[96],q[95];
rz(0.5*pi) q[95];
rz(0.5*pi) q[96];
rx(0.5*pi) q[95];
rx(0.5*pi) q[96];
rz(0.5*pi) q[95];
rz(0.5*pi) q[96];
rx(0.5*pi) q[95];
rz(0.5*pi) q[96];
rz(0.5*pi) q[95];
rx(0.5*pi) q[96];
rz(0.5*pi) q[95];
rz(0.5*pi) q[96];
rx(0.5*pi) q[95];
rz(0.5*pi) q[95];
rzz(0.5*pi) q[95],q[94];
rz(0.5*pi) q[94];
rz(0.5*pi) q[95];
rx(0.5*pi) q[94];
rx(0.5*pi) q[95];
rz(0.5*pi) q[94];
rz(0.5*pi) q[95];
rx(0.5*pi) q[94];
rz(0.5*pi) q[95];
rz(0.5*pi) q[94];
rx(0.5*pi) q[95];
rz(0.5*pi) q[94];
rz(0.5*pi) q[95];
rx(0.5*pi) q[94];
rz(0.5*pi) q[94];
rzz(0.5*pi) q[94],q[93];
rz(0.5*pi) q[93];
rz(0.5*pi) q[94];
rx(0.5*pi) q[93];
rx(0.5*pi) q[94];
rz(0.5*pi) q[93];
rz(0.5*pi) q[94];
rx(0.5*pi) q[93];
rz(0.5*pi) q[94];
rz(0.5*pi) q[93];
rx(0.5*pi) q[94];
rz(0.5*pi) q[93];
rz(0.5*pi) q[94];
rx(0.5*pi) q[93];
rz(0.5*pi) q[93];
rzz(0.5*pi) q[93],q[92];
rz(0.5*pi) q[92];
rz(0.5*pi) q[93];
rx(0.5*pi) q[92];
rx(0.5*pi) q[93];
rz(0.5*pi) q[92];
rz(0.5*pi) q[93];
rx(0.5*pi) q[92];
rz(0.5*pi) q[93];
rz(0.5*pi) q[92];
rx(0.5*pi) q[93];
rz(0.5*pi) q[92];
rz(0.5*pi) q[93];
rx(0.5*pi) q[92];
rz(0.5*pi) q[92];
rzz(0.5*pi) q[92],q[91];
rz(0.5*pi) q[91];
rz(0.5*pi) q[92];
rx(0.5*pi) q[91];
rx(0.5*pi) q[92];
rz(0.5*pi) q[91];
rz(0.5*pi) q[92];
rx(0.5*pi) q[91];
rz(0.5*pi) q[92];
rz(0.5*pi) q[91];
rx(0.5*pi) q[92];
rz(0.5*pi) q[91];
rz(0.5*pi) q[92];
rx(0.5*pi) q[91];
rz(0.5*pi) q[91];
rzz(0.5*pi) q[91],q[90];
rz(0.5*pi) q[90];
rz(0.5*pi) q[91];
rx(0.5*pi) q[90];
rx(0.5*pi) q[91];
rz(0.5*pi) q[90];
rz(0.5*pi) q[91];
rx(0.5*pi) q[90];
rz(0.5*pi) q[91];
rz(0.5*pi) q[90];
rx(0.5*pi) q[91];
rz(0.5*pi) q[90];
rz(0.5*pi) q[91];
rx(0.5*pi) q[90];
rz(0.5*pi) q[90];
rzz(0.5*pi) q[90],q[89];
rz(0.5*pi) q[89];
rz(0.5*pi) q[90];
rx(0.5*pi) q[89];
rx(0.5*pi) q[90];
rz(0.5*pi) q[89];
rz(0.5*pi) q[90];
rx(0.5*pi) q[89];
rz(0.5*pi) q[90];
rz(0.5*pi) q[89];
rx(0.5*pi) q[90];
rz(0.5*pi) q[89];
rz(0.5*pi) q[90];
rx(0.5*pi) q[89];
rz(0.5*pi) q[89];
rzz(0.5*pi) q[89],q[88];
rz(0.5*pi) q[88];
rz(0.5*pi) q[89];
rx(0.5*pi) q[88];
rx(0.5*pi) q[89];
rz(0.5*pi) q[88];
rz(0.5*pi) q[89];
rx(0.5*pi) q[88];
rz(0.5*pi) q[89];
rz(0.5*pi) q[88];
rx(0.5*pi) q[89];
rz(0.5*pi) q[88];
rz(0.5*pi) q[89];
rx(0.5*pi) q[88];
rz(0.5*pi) q[88];
rzz(0.5*pi) q[88],q[87];
rz(0.5*pi) q[87];
rz(0.5*pi) q[88];
rx(0.5*pi) q[87];
rx(0.5*pi) q[88];
rz(0.5*pi) q[87];
rz(0.5*pi) q[88];
rx(0.5*pi) q[87];
rz(0.5*pi) q[88];
rz(0.5*pi) q[87];
rx(0.5*pi) q[88];
rz(0.5*pi) q[87];
rz(0.5*pi) q[88];
rx(0.5*pi) q[87];
rz(0.5*pi) q[87];
rzz(0.5*pi) q[87],q[86];
rz(0.5*pi) q[86];
rz(0.5*pi) q[87];
rx(0.5*pi) q[86];
rx(0.5*pi) q[87];
rz(0.5*pi) q[86];
rz(0.5*pi) q[87];
rx(0.5*pi) q[86];
rz(0.5*pi) q[87];
rz(0.5*pi) q[86];
rx(0.5*pi) q[87];
rz(0.5*pi) q[86];
rz(0.5*pi) q[87];
rx(0.5*pi) q[86];
rz(0.5*pi) q[86];
rzz(0.5*pi) q[86],q[85];
rz(0.5*pi) q[85];
rz(0.5*pi) q[86];
rx(0.5*pi) q[85];
rx(0.5*pi) q[86];
rz(0.5*pi) q[85];
rz(0.5*pi) q[86];
rx(0.5*pi) q[85];
rz(0.5*pi) q[86];
rz(0.5*pi) q[85];
rx(0.5*pi) q[86];
rz(0.5*pi) q[85];
rz(0.5*pi) q[86];
rx(0.5*pi) q[85];
rz(0.5*pi) q[85];
rzz(0.5*pi) q[85],q[84];
rz(0.5*pi) q[84];
rz(0.5*pi) q[85];
rx(0.5*pi) q[84];
rx(0.5*pi) q[85];
rz(0.5*pi) q[84];
rz(0.5*pi) q[85];
rx(0.5*pi) q[84];
rz(0.5*pi) q[85];
rz(0.5*pi) q[84];
rx(0.5*pi) q[85];
rz(0.5*pi) q[84];
rz(0.5*pi) q[85];
rx(0.5*pi) q[84];
rz(0.5*pi) q[84];
rzz(0.5*pi) q[84],q[83];
rz(0.5*pi) q[83];
rz(0.5*pi) q[84];
rx(0.5*pi) q[83];
rx(0.5*pi) q[84];
rz(0.5*pi) q[83];
rz(0.5*pi) q[84];
rx(0.5*pi) q[83];
rz(0.5*pi) q[84];
rz(0.5*pi) q[83];
rx(0.5*pi) q[84];
rz(0.5*pi) q[83];
rz(0.5*pi) q[84];
rx(0.5*pi) q[83];
rz(0.5*pi) q[83];
rzz(0.5*pi) q[83],q[82];
rz(0.5*pi) q[82];
rz(0.5*pi) q[83];
rx(0.5*pi) q[82];
rx(0.5*pi) q[83];
rz(0.5*pi) q[82];
rz(0.5*pi) q[83];
rx(0.5*pi) q[82];
rz(0.5*pi) q[83];
rz(0.5*pi) q[82];
rx(0.5*pi) q[83];
rz(0.5*pi) q[82];
rz(0.5*pi) q[83];
rx(0.5*pi) q[82];
rz(0.5*pi) q[82];
rzz(0.5*pi) q[82],q[81];
rz(0.5*pi) q[81];
rz(0.5*pi) q[82];
rx(0.5*pi) q[81];
rx(0.5*pi) q[82];
rz(0.5*pi) q[81];
rz(0.5*pi) q[82];
rx(0.5*pi) q[81];
rz(0.5*pi) q[82];
rz(0.5*pi) q[81];
rx(0.5*pi) q[82];
rz(0.5*pi) q[81];
rz(0.5*pi) q[82];
rx(0.5*pi) q[81];
rz(0.5*pi) q[81];
rzz(0.5*pi) q[81],q[80];
rz(0.5*pi) q[80];
rz(0.5*pi) q[81];
rx(0.5*pi) q[80];
rx(0.5*pi) q[81];
rz(0.5*pi) q[80];
rz(0.5*pi) q[81];
rx(0.5*pi) q[80];
rz(0.5*pi) q[81];
rz(0.5*pi) q[80];
rx(0.5*pi) q[81];
rz(0.5*pi) q[80];
rz(0.5*pi) q[81];
rx(0.5*pi) q[80];
rz(0.5*pi) q[80];
rzz(0.5*pi) q[80],q[79];
rz(0.5*pi) q[79];
rz(0.5*pi) q[80];
rx(0.5*pi) q[79];
rx(0.5*pi) q[80];
rz(0.5*pi) q[79];
rz(0.5*pi) q[80];
rx(0.5*pi) q[79];
rz(0.5*pi) q[80];
rz(0.5*pi) q[79];
rx(0.5*pi) q[80];
rz(0.5*pi) q[79];
rz(0.5*pi) q[80];
rx(0.5*pi) q[79];
rz(0.5*pi) q[79];
rzz(0.5*pi) q[79],q[78];
rz(0.5*pi) q[78];
rz(0.5*pi) q[79];
rx(0.5*pi) q[78];
rx(0.5*pi) q[79];
rz(0.5*pi) q[78];
rz(0.5*pi) q[79];
rx(0.5*pi) q[78];
rz(0.5*pi) q[79];
rz(0.5*pi) q[78];
rx(0.5*pi) q[79];
rz(0.5*pi) q[78];
rz(0.5*pi) q[79];
rx(0.5*pi) q[78];
rz(0.5*pi) q[78];
rzz(0.5*pi) q[78],q[77];
rz(0.5*pi) q[77];
rz(0.5*pi) q[78];
rx(0.5*pi) q[77];
rx(0.5*pi) q[78];
rz(0.5*pi) q[77];
rz(0.5*pi) q[78];
rx(0.5*pi) q[77];
rz(0.5*pi) q[78];
rz(0.5*pi) q[77];
rx(0.5*pi) q[78];
rz(0.5*pi) q[77];
rz(0.5*pi) q[78];
rx(0.5*pi) q[77];
rz(0.5*pi) q[77];
rzz(0.5*pi) q[77],q[76];
rz(0.5*pi) q[76];
rz(0.5*pi) q[77];
rx(0.5*pi) q[76];
rx(0.5*pi) q[77];
rz(0.5*pi) q[76];
rz(0.5*pi) q[77];
rx(0.5*pi) q[76];
rz(0.5*pi) q[77];
rz(0.5*pi) q[76];
rx(0.5*pi) q[77];
rz(0.5*pi) q[76];
rz(0.5*pi) q[77];
rx(0.5*pi) q[76];
rz(0.5*pi) q[76];
rzz(0.5*pi) q[76],q[75];
rz(0.5*pi) q[75];
rz(0.5*pi) q[76];
rx(0.5*pi) q[75];
rx(0.5*pi) q[76];
rz(0.5*pi) q[75];
rz(0.5*pi) q[76];
rx(0.5*pi) q[75];
rz(0.5*pi) q[76];
rz(0.5*pi) q[75];
rx(0.5*pi) q[76];
rz(0.5*pi) q[75];
rz(0.5*pi) q[76];
rx(0.5*pi) q[75];
rz(0.5*pi) q[75];
rzz(0.5*pi) q[75],q[74];
rz(0.5*pi) q[74];
rz(0.5*pi) q[75];
rx(0.5*pi) q[74];
rx(0.5*pi) q[75];
rz(0.5*pi) q[74];
rz(0.5*pi) q[75];
rx(0.5*pi) q[74];
rz(0.5*pi) q[75];
rz(0.5*pi) q[74];
rx(0.5*pi) q[75];
rz(0.5*pi) q[74];
rz(0.5*pi) q[75];
rx(0.5*pi) q[74];
rz(0.5*pi) q[74];
rzz(0.5*pi) q[74],q[73];
rz(0.5*pi) q[73];
rz(0.5*pi) q[74];
rx(0.5*pi) q[73];
rx(0.5*pi) q[74];
rz(0.5*pi) q[73];
rz(0.5*pi) q[74];
rx(0.5*pi) q[73];
rz(0.5*pi) q[74];
rz(0.5*pi) q[73];
rx(0.5*pi) q[74];
rz(0.5*pi) q[73];
rz(0.5*pi) q[74];
rx(0.5*pi) q[73];
rz(0.5*pi) q[73];
rzz(0.5*pi) q[73],q[72];
rz(0.5*pi) q[72];
rz(0.5*pi) q[73];
rx(0.5*pi) q[72];
rx(0.5*pi) q[73];
rz(0.5*pi) q[72];
rz(0.5*pi) q[73];
rx(0.5*pi) q[72];
rz(0.5*pi) q[73];
rz(0.5*pi) q[72];
rx(0.5*pi) q[73];
rz(0.5*pi) q[72];
rz(0.5*pi) q[73];
rx(0.5*pi) q[72];
rz(0.5*pi) q[72];
rzz(0.5*pi) q[72],q[71];
rz(0.5*pi) q[71];
rz(0.5*pi) q[72];
rx(0.5*pi) q[71];
rx(0.5*pi) q[72];
rz(0.5*pi) q[71];
rz(0.5*pi) q[72];
rx(0.5*pi) q[71];
rz(0.5*pi) q[72];
rz(0.5*pi) q[71];
rx(0.5*pi) q[72];
rz(0.5*pi) q[71];
rz(0.5*pi) q[72];
rx(0.5*pi) q[71];
rz(0.5*pi) q[71];
rzz(0.5*pi) q[71],q[70];
rz(0.5*pi) q[70];
rz(0.5*pi) q[71];
rx(0.5*pi) q[70];
rx(0.5*pi) q[71];
rz(0.5*pi) q[70];
rz(0.5*pi) q[71];
rx(0.5*pi) q[70];
rz(0.5*pi) q[71];
rz(0.5*pi) q[70];
rx(0.5*pi) q[71];
rz(0.5*pi) q[70];
rz(0.5*pi) q[71];
rx(0.5*pi) q[70];
rz(0.5*pi) q[70];
rzz(0.5*pi) q[70],q[69];
rz(0.5*pi) q[69];
rz(0.5*pi) q[70];
rx(0.5*pi) q[69];
rx(0.5*pi) q[70];
rz(0.5*pi) q[69];
rz(0.5*pi) q[70];
rx(0.5*pi) q[69];
rz(0.5*pi) q[70];
rz(0.5*pi) q[69];
rx(0.5*pi) q[70];
rz(0.5*pi) q[69];
rz(0.5*pi) q[70];
rx(0.5*pi) q[69];
rz(0.5*pi) q[69];
rzz(0.5*pi) q[69],q[68];
rz(0.5*pi) q[68];
rz(0.5*pi) q[69];
rx(0.5*pi) q[68];
rx(0.5*pi) q[69];
rz(0.5*pi) q[68];
rz(0.5*pi) q[69];
rx(0.5*pi) q[68];
rz(0.5*pi) q[69];
rz(0.5*pi) q[68];
rx(0.5*pi) q[69];
rz(0.5*pi) q[68];
rz(0.5*pi) q[69];
rx(0.5*pi) q[68];
rz(0.5*pi) q[68];
rzz(0.5*pi) q[68],q[67];
rz(0.5*pi) q[67];
rz(0.5*pi) q[68];
rx(0.5*pi) q[67];
rx(0.5*pi) q[68];
rz(0.5*pi) q[67];
rz(0.5*pi) q[68];
rx(0.5*pi) q[67];
rz(0.5*pi) q[68];
rz(0.5*pi) q[67];
rx(0.5*pi) q[68];
rz(0.5*pi) q[67];
rz(0.5*pi) q[68];
rx(0.5*pi) q[67];
rz(0.5*pi) q[67];
rzz(0.5*pi) q[67],q[66];
rz(0.5*pi) q[66];
rz(0.5*pi) q[67];
rx(0.5*pi) q[66];
rx(0.5*pi) q[67];
rz(0.5*pi) q[66];
rz(0.5*pi) q[67];
rx(0.5*pi) q[66];
rz(0.5*pi) q[67];
rz(0.5*pi) q[66];
rx(0.5*pi) q[67];
rz(0.5*pi) q[66];
rz(0.5*pi) q[67];
rx(0.5*pi) q[66];
rz(0.5*pi) q[66];
rzz(0.5*pi) q[66],q[65];
rz(0.5*pi) q[65];
rz(0.5*pi) q[66];
rx(0.5*pi) q[65];
rx(0.5*pi) q[66];
rz(0.5*pi) q[65];
rz(0.5*pi) q[66];
rx(0.5*pi) q[65];
rz(0.5*pi) q[66];
rz(0.5*pi) q[65];
rx(0.5*pi) q[66];
rz(0.5*pi) q[65];
rz(0.5*pi) q[66];
rx(0.5*pi) q[65];
rz(0.5*pi) q[65];
rzz(0.5*pi) q[65],q[64];
rz(0.5*pi) q[64];
rz(0.5*pi) q[65];
rx(0.5*pi) q[64];
rx(0.5*pi) q[65];
rz(0.5*pi) q[64];
rz(0.5*pi) q[65];
rx(0.5*pi) q[64];
rz(0.5*pi) q[65];
rz(0.5*pi) q[64];
rx(0.5*pi) q[65];
rz(0.5*pi) q[64];
rz(0.5*pi) q[65];
rx(0.5*pi) q[64];
rz(0.5*pi) q[64];
rzz(0.5*pi) q[64],q[63];
rz(0.5*pi) q[63];
rz(0.5*pi) q[64];
rx(0.5*pi) q[63];
rx(0.5*pi) q[64];
rz(0.5*pi) q[63];
rz(0.5*pi) q[64];
rx(0.5*pi) q[63];
rz(0.5*pi) q[64];
rz(0.5*pi) q[63];
rx(0.5*pi) q[64];
rz(0.5*pi) q[63];
rz(0.5*pi) q[64];
rx(0.5*pi) q[63];
rz(0.5*pi) q[63];
rzz(0.5*pi) q[63],q[62];
rz(0.5*pi) q[62];
rz(0.5*pi) q[63];
rx(0.5*pi) q[62];
rx(0.5*pi) q[63];
rz(0.5*pi) q[62];
rz(0.5*pi) q[63];
rx(0.5*pi) q[62];
rz(0.5*pi) q[63];
rz(0.5*pi) q[62];
rx(0.5*pi) q[63];
rz(0.5*pi) q[62];
rz(0.5*pi) q[63];
rx(0.5*pi) q[62];
rz(0.5*pi) q[62];
rzz(0.5*pi) q[62],q[61];
rz(0.5*pi) q[61];
rz(0.5*pi) q[62];
rx(0.5*pi) q[61];
rx(0.5*pi) q[62];
rz(0.5*pi) q[61];
rz(0.5*pi) q[62];
rx(0.5*pi) q[61];
rz(0.5*pi) q[62];
rz(0.5*pi) q[61];
rx(0.5*pi) q[62];
rz(0.5*pi) q[61];
rz(0.5*pi) q[62];
rx(0.5*pi) q[61];
rz(0.5*pi) q[61];
rzz(0.5*pi) q[61],q[60];
rz(0.5*pi) q[60];
rz(0.5*pi) q[61];
rx(0.5*pi) q[60];
rx(0.5*pi) q[61];
rz(0.5*pi) q[60];
rz(0.5*pi) q[61];
rx(0.5*pi) q[60];
rz(0.5*pi) q[61];
rz(0.5*pi) q[60];
rx(0.5*pi) q[61];
rz(0.5*pi) q[60];
rz(0.5*pi) q[61];
rx(0.5*pi) q[60];
rz(0.5*pi) q[60];
rzz(0.5*pi) q[60],q[59];
rz(0.5*pi) q[59];
rz(0.5*pi) q[60];
rx(0.5*pi) q[59];
rx(0.5*pi) q[60];
rz(0.5*pi) q[59];
rz(0.5*pi) q[60];
rx(0.5*pi) q[59];
rz(0.5*pi) q[60];
rz(0.5*pi) q[59];
rx(0.5*pi) q[60];
rz(0.5*pi) q[59];
rz(0.5*pi) q[60];
rx(0.5*pi) q[59];
rz(0.5*pi) q[59];
rzz(0.5*pi) q[59],q[58];
rz(0.5*pi) q[58];
rz(0.5*pi) q[59];
rx(0.5*pi) q[58];
rx(0.5*pi) q[59];
rz(0.5*pi) q[58];
rz(0.5*pi) q[59];
rx(0.5*pi) q[58];
rz(0.5*pi) q[59];
rz(0.5*pi) q[58];
rx(0.5*pi) q[59];
rz(0.5*pi) q[58];
rz(0.5*pi) q[59];
rx(0.5*pi) q[58];
rz(0.5*pi) q[58];
rzz(0.5*pi) q[58],q[57];
rz(0.5*pi) q[57];
rz(0.5*pi) q[58];
rx(0.5*pi) q[57];
rx(0.5*pi) q[58];
rz(0.5*pi) q[57];
rz(0.5*pi) q[58];
rx(0.5*pi) q[57];
rz(0.5*pi) q[58];
rz(0.5*pi) q[57];
rx(0.5*pi) q[58];
rz(0.5*pi) q[57];
rz(0.5*pi) q[58];
rx(0.5*pi) q[57];
rz(0.5*pi) q[57];
rzz(0.5*pi) q[57],q[56];
rz(0.5*pi) q[56];
rz(0.5*pi) q[57];
rx(0.5*pi) q[56];
rx(0.5*pi) q[57];
rz(0.5*pi) q[56];
rz(0.5*pi) q[57];
rx(0.5*pi) q[56];
rz(0.5*pi) q[57];
rz(0.5*pi) q[56];
rx(0.5*pi) q[57];
rz(0.5*pi) q[56];
rz(0.5*pi) q[57];
rx(0.5*pi) q[56];
rz(0.5*pi) q[56];
rzz(0.5*pi) q[56],q[55];
rz(0.5*pi) q[55];
rz(0.5*pi) q[56];
rx(0.5*pi) q[55];
rx(0.5*pi) q[56];
rz(0.5*pi) q[55];
rz(0.5*pi) q[56];
rx(0.5*pi) q[55];
rz(0.5*pi) q[56];
rz(0.5*pi) q[55];
rx(0.5*pi) q[56];
rz(0.5*pi) q[55];
rz(0.5*pi) q[56];
rx(0.5*pi) q[55];
rz(0.5*pi) q[55];
rzz(0.5*pi) q[55],q[54];
rz(0.5*pi) q[54];
rz(0.5*pi) q[55];
rx(0.5*pi) q[54];
rx(0.5*pi) q[55];
rz(0.5*pi) q[54];
rz(0.5*pi) q[55];
rx(0.5*pi) q[54];
rz(0.5*pi) q[55];
rz(0.5*pi) q[54];
rx(0.5*pi) q[55];
rz(0.5*pi) q[54];
rz(0.5*pi) q[55];
rx(0.5*pi) q[54];
rz(0.5*pi) q[54];
rzz(0.5*pi) q[54],q[53];
rz(0.5*pi) q[53];
rz(0.5*pi) q[54];
rx(0.5*pi) q[53];
rx(0.5*pi) q[54];
rz(0.5*pi) q[53];
rz(0.5*pi) q[54];
rx(0.5*pi) q[53];
rz(0.5*pi) q[54];
rz(0.5*pi) q[53];
rx(0.5*pi) q[54];
rz(0.5*pi) q[53];
rz(0.5*pi) q[54];
rx(0.5*pi) q[53];
rz(0.5*pi) q[53];
rzz(0.5*pi) q[53],q[52];
rz(0.5*pi) q[52];
rz(0.5*pi) q[53];
rx(0.5*pi) q[52];
rx(0.5*pi) q[53];
rz(0.5*pi) q[52];
rz(0.5*pi) q[53];
rx(0.5*pi) q[52];
rz(0.5*pi) q[53];
rz(0.5*pi) q[52];
rx(0.5*pi) q[53];
rz(0.5*pi) q[52];
rz(0.5*pi) q[53];
rx(0.5*pi) q[52];
rz(0.5*pi) q[52];
rzz(0.5*pi) q[52],q[51];
rz(0.5*pi) q[51];
rz(0.5*pi) q[52];
rx(0.5*pi) q[51];
rx(0.5*pi) q[52];
rz(0.5*pi) q[51];
rz(0.5*pi) q[52];
rx(0.5*pi) q[51];
rz(0.5*pi) q[52];
rz(0.5*pi) q[51];
rx(0.5*pi) q[52];
rz(0.5*pi) q[51];
rz(0.5*pi) q[52];
rx(0.5*pi) q[51];
rz(0.5*pi) q[51];
rzz(0.5*pi) q[51],q[50];
rz(0.5*pi) q[50];
rz(0.5*pi) q[51];
rx(0.5*pi) q[50];
rx(0.5*pi) q[51];
rz(0.5*pi) q[50];
rz(0.5*pi) q[51];
rx(0.5*pi) q[50];
rz(0.5*pi) q[51];
rz(0.5*pi) q[50];
rx(0.5*pi) q[51];
rz(0.5*pi) q[50];
rz(0.5*pi) q[51];
rx(0.5*pi) q[50];
rz(0.5*pi) q[50];
rzz(0.5*pi) q[50],q[49];
rz(0.5*pi) q[49];
rz(0.5*pi) q[50];
rx(0.5*pi) q[49];
rx(0.5*pi) q[50];
rz(0.5*pi) q[49];
rz(0.5*pi) q[50];
rx(0.5*pi) q[49];
rz(0.5*pi) q[50];
rz(0.5*pi) q[49];
rx(0.5*pi) q[50];
rz(0.5*pi) q[49];
rz(0.5*pi) q[50];
rx(0.5*pi) q[49];
rz(0.5*pi) q[49];
rzz(0.5*pi) q[49],q[48];
rz(0.5*pi) q[48];
rz(0.5*pi) q[49];
rx(0.5*pi) q[48];
rx(0.5*pi) q[49];
rz(0.5*pi) q[48];
rz(0.5*pi) q[49];
rx(0.5*pi) q[48];
rz(0.5*pi) q[49];
rz(0.5*pi) q[48];
rx(0.5*pi) q[49];
rz(0.5*pi) q[48];
rz(0.5*pi) q[49];
rx(0.5*pi) q[48];
rz(0.5*pi) q[48];
rzz(0.5*pi) q[48],q[47];
rz(0.5*pi) q[47];
rz(0.5*pi) q[48];
rx(0.5*pi) q[47];
rx(0.5*pi) q[48];
rz(0.5*pi) q[47];
rz(0.5*pi) q[48];
rx(0.5*pi) q[47];
rz(0.5*pi) q[48];
rz(0.5*pi) q[47];
rx(0.5*pi) q[48];
rz(0.5*pi) q[47];
rz(0.5*pi) q[48];
rx(0.5*pi) q[47];
rz(0.5*pi) q[47];
rzz(0.5*pi) q[47],q[46];
rz(0.5*pi) q[46];
rz(0.5*pi) q[47];
rx(0.5*pi) q[46];
rx(0.5*pi) q[47];
rz(0.5*pi) q[46];
rz(0.5*pi) q[47];
rx(0.5*pi) q[46];
rz(0.5*pi) q[47];
rz(0.5*pi) q[46];
rx(0.5*pi) q[47];
rz(0.5*pi) q[46];
rz(0.5*pi) q[47];
rx(0.5*pi) q[46];
rz(0.5*pi) q[46];
rzz(0.5*pi) q[46],q[45];
rz(0.5*pi) q[45];
rz(0.5*pi) q[46];
rx(0.5*pi) q[45];
rx(0.5*pi) q[46];
rz(0.5*pi) q[45];
rz(0.5*pi) q[46];
rx(0.5*pi) q[45];
rz(0.5*pi) q[46];
rz(0.5*pi) q[45];
rx(0.5*pi) q[46];
rz(0.5*pi) q[45];
rz(0.5*pi) q[46];
rx(0.5*pi) q[45];
rz(0.5*pi) q[45];
rzz(0.5*pi) q[45],q[44];
rz(0.5*pi) q[44];
rz(0.5*pi) q[45];
rx(0.5*pi) q[44];
rx(0.5*pi) q[45];
rz(0.5*pi) q[44];
rz(0.5*pi) q[45];
rx(0.5*pi) q[44];
rz(0.5*pi) q[45];
rz(0.5*pi) q[44];
rx(0.5*pi) q[45];
rz(0.5*pi) q[44];
rz(0.5*pi) q[45];
rx(0.5*pi) q[44];
rz(0.5*pi) q[44];
rzz(0.5*pi) q[44],q[43];
rz(0.5*pi) q[43];
rz(0.5*pi) q[44];
rx(0.5*pi) q[43];
rx(0.5*pi) q[44];
rz(0.5*pi) q[43];
rz(0.5*pi) q[44];
rx(0.5*pi) q[43];
rz(0.5*pi) q[44];
rz(0.5*pi) q[43];
rx(0.5*pi) q[44];
rz(0.5*pi) q[43];
rz(0.5*pi) q[44];
rx(0.5*pi) q[43];
rz(0.5*pi) q[43];
rzz(0.5*pi) q[43],q[42];
rz(0.5*pi) q[42];
rz(0.5*pi) q[43];
rx(0.5*pi) q[42];
rx(0.5*pi) q[43];
rz(0.5*pi) q[42];
rz(0.5*pi) q[43];
rx(0.5*pi) q[42];
rz(0.5*pi) q[43];
rz(0.5*pi) q[42];
rx(0.5*pi) q[43];
rz(0.5*pi) q[42];
rz(0.5*pi) q[43];
rx(0.5*pi) q[42];
rz(0.5*pi) q[42];
rzz(0.5*pi) q[42],q[41];
rz(0.5*pi) q[41];
rz(0.5*pi) q[42];
rx(0.5*pi) q[41];
rx(0.5*pi) q[42];
rz(0.5*pi) q[41];
rz(0.5*pi) q[42];
rx(0.5*pi) q[41];
rz(0.5*pi) q[42];
rz(0.5*pi) q[41];
rx(0.5*pi) q[42];
rz(0.5*pi) q[41];
rz(0.5*pi) q[42];
rx(0.5*pi) q[41];
rz(0.5*pi) q[41];
rzz(0.5*pi) q[41],q[40];
rz(0.5*pi) q[40];
rz(0.5*pi) q[41];
rx(0.5*pi) q[40];
rx(0.5*pi) q[41];
rz(0.5*pi) q[40];
rz(0.5*pi) q[41];
rx(0.5*pi) q[40];
rz(0.5*pi) q[41];
rz(0.5*pi) q[40];
rx(0.5*pi) q[41];
rz(0.5*pi) q[40];
rz(0.5*pi) q[41];
rx(0.5*pi) q[40];
rz(0.5*pi) q[40];
rzz(0.5*pi) q[40],q[39];
rz(0.5*pi) q[39];
rz(0.5*pi) q[40];
rx(0.5*pi) q[39];
rx(0.5*pi) q[40];
rz(0.5*pi) q[39];
rz(0.5*pi) q[40];
rx(0.5*pi) q[39];
rz(0.5*pi) q[40];
rz(0.5*pi) q[39];
rx(0.5*pi) q[40];
rz(0.5*pi) q[39];
rz(0.5*pi) q[40];
rx(0.5*pi) q[39];
rz(0.5*pi) q[39];
rzz(0.5*pi) q[39],q[38];
rz(0.5*pi) q[38];
rz(0.5*pi) q[39];
rx(0.5*pi) q[38];
rx(0.5*pi) q[39];
rz(0.5*pi) q[38];
rz(0.5*pi) q[39];
rx(0.5*pi) q[38];
rz(0.5*pi) q[39];
rz(0.5*pi) q[38];
rx(0.5*pi) q[39];
rz(0.5*pi) q[38];
rz(0.5*pi) q[39];
rx(0.5*pi) q[38];
rz(0.5*pi) q[38];
rzz(0.5*pi) q[38],q[37];
rz(0.5*pi) q[37];
rz(0.5*pi) q[38];
rx(0.5*pi) q[37];
rx(0.5*pi) q[38];
rz(0.5*pi) q[37];
rz(0.5*pi) q[38];
rx(0.5*pi) q[37];
rz(0.5*pi) q[38];
rz(0.5*pi) q[37];
rx(0.5*pi) q[38];
rz(0.5*pi) q[37];
rz(0.5*pi) q[38];
rx(0.5*pi) q[37];
rz(0.5*pi) q[37];
rzz(0.5*pi) q[37],q[36];
rz(0.5*pi) q[36];
rz(0.5*pi) q[37];
rx(0.5*pi) q[36];
rx(0.5*pi) q[37];
rz(0.5*pi) q[36];
rz(0.5*pi) q[37];
rx(0.5*pi) q[36];
rz(0.5*pi) q[37];
rz(0.5*pi) q[36];
rx(0.5*pi) q[37];
rz(0.5*pi) q[36];
rz(0.5*pi) q[37];
rx(0.5*pi) q[36];
rz(0.5*pi) q[36];
rzz(0.5*pi) q[36],q[35];
rz(0.5*pi) q[35];
rz(0.5*pi) q[36];
rx(0.5*pi) q[35];
rx(0.5*pi) q[36];
rz(0.5*pi) q[35];
rz(0.5*pi) q[36];
rx(0.5*pi) q[35];
rz(0.5*pi) q[36];
rz(0.5*pi) q[35];
rx(0.5*pi) q[36];
rz(0.5*pi) q[35];
rz(0.5*pi) q[36];
rx(0.5*pi) q[35];
rz(0.5*pi) q[35];
rzz(0.5*pi) q[35],q[34];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rx(0.5*pi) q[34];
rx(0.5*pi) q[35];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rx(0.5*pi) q[34];
rz(0.5*pi) q[35];
rz(0.5*pi) q[34];
rx(0.5*pi) q[35];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rx(0.5*pi) q[34];
rz(0.5*pi) q[34];
rzz(0.5*pi) q[34],q[33];
rz(0.5*pi) q[33];
rz(0.5*pi) q[34];
rx(0.5*pi) q[33];
rx(0.5*pi) q[34];
rz(0.5*pi) q[33];
rz(0.5*pi) q[34];
rx(0.5*pi) q[33];
rz(0.5*pi) q[34];
rz(0.5*pi) q[33];
rx(0.5*pi) q[34];
rz(0.5*pi) q[33];
rz(0.5*pi) q[34];
rx(0.5*pi) q[33];
rz(0.5*pi) q[33];
rzz(0.5*pi) q[33],q[32];
rz(0.5*pi) q[32];
rz(0.5*pi) q[33];
rx(0.5*pi) q[32];
rx(0.5*pi) q[33];
rz(0.5*pi) q[32];
rz(0.5*pi) q[33];
rx(0.5*pi) q[32];
rz(0.5*pi) q[33];
rz(0.5*pi) q[32];
rx(0.5*pi) q[33];
rz(0.5*pi) q[32];
rz(0.5*pi) q[33];
rx(0.5*pi) q[32];
rz(0.5*pi) q[32];
rzz(0.5*pi) q[32],q[31];
rz(0.5*pi) q[31];
rz(0.5*pi) q[32];
rx(0.5*pi) q[31];
rx(0.5*pi) q[32];
rz(0.5*pi) q[31];
rz(0.5*pi) q[32];
rx(0.5*pi) q[31];
rz(0.5*pi) q[32];
rz(0.5*pi) q[31];
rx(0.5*pi) q[32];
rz(0.5*pi) q[31];
rz(0.5*pi) q[32];
rx(0.5*pi) q[31];
rz(0.5*pi) q[31];
rzz(0.5*pi) q[31],q[30];
rz(0.5*pi) q[30];
rz(0.5*pi) q[31];
rx(0.5*pi) q[30];
rx(0.5*pi) q[31];
rz(0.5*pi) q[30];
rz(0.5*pi) q[31];
rx(0.5*pi) q[30];
rz(0.5*pi) q[31];
rz(0.5*pi) q[30];
rx(0.5*pi) q[31];
rz(0.5*pi) q[30];
rz(0.5*pi) q[31];
rx(0.5*pi) q[30];
rz(0.5*pi) q[30];
rzz(0.5*pi) q[30],q[29];
rz(0.5*pi) q[29];
rz(0.5*pi) q[30];
rx(0.5*pi) q[29];
rx(0.5*pi) q[30];
rz(0.5*pi) q[29];
rz(0.5*pi) q[30];
rx(0.5*pi) q[29];
rz(0.5*pi) q[30];
rz(0.5*pi) q[29];
rx(0.5*pi) q[30];
rz(0.5*pi) q[29];
rz(0.5*pi) q[30];
rx(0.5*pi) q[29];
rz(0.5*pi) q[29];
rzz(0.5*pi) q[29],q[28];
rz(0.5*pi) q[28];
rz(0.5*pi) q[29];
rx(0.5*pi) q[28];
rx(0.5*pi) q[29];
rz(0.5*pi) q[28];
rz(0.5*pi) q[29];
rx(0.5*pi) q[28];
rz(0.5*pi) q[29];
rz(0.5*pi) q[28];
rx(0.5*pi) q[29];
rz(0.5*pi) q[28];
rz(0.5*pi) q[29];
rx(0.5*pi) q[28];
rz(0.5*pi) q[28];
rzz(0.5*pi) q[28],q[27];
rz(0.5*pi) q[27];
rz(0.5*pi) q[28];
rx(0.5*pi) q[27];
rx(0.5*pi) q[28];
rz(0.5*pi) q[27];
rz(0.5*pi) q[28];
rx(0.5*pi) q[27];
rz(0.5*pi) q[28];
rz(0.5*pi) q[27];
rx(0.5*pi) q[28];
rz(0.5*pi) q[27];
rz(0.5*pi) q[28];
rx(0.5*pi) q[27];
rz(0.5*pi) q[27];
rzz(0.5*pi) q[27],q[26];
rz(0.5*pi) q[26];
rz(0.5*pi) q[27];
rx(0.5*pi) q[26];
rx(0.5*pi) q[27];
rz(0.5*pi) q[26];
rz(0.5*pi) q[27];
rx(0.5*pi) q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[26];
rx(0.5*pi) q[27];
rz(0.5*pi) q[26];
rz(0.5*pi) q[27];
rx(0.5*pi) q[26];
rz(0.5*pi) q[26];
rzz(0.5*pi) q[26],q[25];
rz(0.5*pi) q[25];
rz(0.5*pi) q[26];
rx(0.5*pi) q[25];
rx(0.5*pi) q[26];
rz(0.5*pi) q[25];
rz(0.5*pi) q[26];
rx(0.5*pi) q[25];
rz(0.5*pi) q[26];
rz(0.5*pi) q[25];
rx(0.5*pi) q[26];
rz(0.5*pi) q[25];
rz(0.5*pi) q[26];
rx(0.5*pi) q[25];
rz(0.5*pi) q[25];
rzz(0.5*pi) q[25],q[24];
rz(0.5*pi) q[24];
rz(0.5*pi) q[25];
rx(0.5*pi) q[24];
rx(0.5*pi) q[25];
rz(0.5*pi) q[24];
rz(0.5*pi) q[25];
rx(0.5*pi) q[24];
rz(0.5*pi) q[25];
rz(0.5*pi) q[24];
rx(0.5*pi) q[25];
rz(0.5*pi) q[24];
rz(0.5*pi) q[25];
rx(0.5*pi) q[24];
rz(0.5*pi) q[24];
rzz(0.5*pi) q[24],q[23];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rx(0.5*pi) q[23];
rx(0.5*pi) q[24];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rx(0.5*pi) q[23];
rz(0.5*pi) q[24];
rz(0.5*pi) q[23];
rx(0.5*pi) q[24];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rx(0.5*pi) q[23];
rz(0.5*pi) q[23];
rzz(0.5*pi) q[23],q[22];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rx(0.5*pi) q[22];
rx(0.5*pi) q[23];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rx(0.5*pi) q[22];
rz(0.5*pi) q[23];
rz(0.5*pi) q[22];
rx(0.5*pi) q[23];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rx(0.5*pi) q[22];
rz(0.5*pi) q[22];
rzz(0.5*pi) q[22],q[21];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rx(0.5*pi) q[21];
rx(0.5*pi) q[22];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rx(0.5*pi) q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[21];
rx(0.5*pi) q[22];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rx(0.5*pi) q[21];
rz(0.5*pi) q[21];
rzz(0.5*pi) q[21],q[20];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rx(0.5*pi) q[20];
rx(0.5*pi) q[21];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rx(0.5*pi) q[20];
rz(0.5*pi) q[21];
rz(0.5*pi) q[20];
rx(0.5*pi) q[21];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rx(0.5*pi) q[20];
rz(0.5*pi) q[20];
rzz(0.5*pi) q[20],q[19];
rz(0.5*pi) q[19];
rz(0.5*pi) q[20];
rx(0.5*pi) q[19];
rx(0.5*pi) q[20];
rz(0.5*pi) q[19];
rz(0.5*pi) q[20];
rx(0.5*pi) q[19];
rz(0.5*pi) q[20];
rz(0.5*pi) q[19];
rx(0.5*pi) q[20];
rz(0.5*pi) q[19];
rz(0.5*pi) q[20];
rx(0.5*pi) q[19];
rz(0.5*pi) q[19];
rzz(0.5*pi) q[19],q[18];
rz(0.5*pi) q[18];
rz(0.5*pi) q[19];
rx(0.5*pi) q[18];
rx(0.5*pi) q[19];
rz(0.5*pi) q[18];
rz(0.5*pi) q[19];
rx(0.5*pi) q[18];
rz(0.5*pi) q[19];
rz(0.5*pi) q[18];
rx(0.5*pi) q[19];
rz(0.5*pi) q[18];
rz(0.5*pi) q[19];
rx(0.5*pi) q[18];
rz(0.5*pi) q[18];
rzz(0.5*pi) q[18],q[17];
rz(0.5*pi) q[17];
rz(0.5*pi) q[18];
rx(0.5*pi) q[17];
rx(0.5*pi) q[18];
rz(0.5*pi) q[17];
rz(0.5*pi) q[18];
rx(0.5*pi) q[17];
rz(0.5*pi) q[18];
rz(0.5*pi) q[17];
rx(0.5*pi) q[18];
rz(0.5*pi) q[17];
rz(0.5*pi) q[18];
rx(0.5*pi) q[17];
rz(0.5*pi) q[17];
rzz(0.5*pi) q[17],q[16];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rx(0.5*pi) q[16];
rx(0.5*pi) q[17];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rx(0.5*pi) q[16];
rz(0.5*pi) q[17];
rz(0.5*pi) q[16];
rx(0.5*pi) q[17];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rx(0.5*pi) q[16];
rz(0.5*pi) q[16];
rzz(0.5*pi) q[16],q[15];
rz(0.5*pi) q[15];
rz(0.5*pi) q[16];
rx(0.5*pi) q[15];
rx(0.5*pi) q[16];
rz(0.5*pi) q[15];
rz(0.5*pi) q[16];
rx(0.5*pi) q[15];
rz(0.5*pi) q[16];
rz(0.5*pi) q[15];
rx(0.5*pi) q[16];
rz(0.5*pi) q[15];
rz(0.5*pi) q[16];
rx(0.5*pi) q[15];
rz(0.5*pi) q[15];
rzz(0.5*pi) q[15],q[14];
rz(0.5*pi) q[14];
rz(0.5*pi) q[15];
rx(0.5*pi) q[14];
rx(0.5*pi) q[15];
rz(0.5*pi) q[14];
rz(0.5*pi) q[15];
rx(0.5*pi) q[14];
rz(0.5*pi) q[15];
rz(0.5*pi) q[14];
rx(0.5*pi) q[15];
rz(0.5*pi) q[14];
rz(0.5*pi) q[15];
rx(0.5*pi) q[14];
rz(0.5*pi) q[14];
rzz(0.5*pi) q[14],q[13];
rz(0.5*pi) q[13];
rz(0.5*pi) q[14];
rx(0.5*pi) q[13];
rx(0.5*pi) q[14];
rz(0.5*pi) q[13];
rz(0.5*pi) q[14];
rx(0.5*pi) q[13];
rz(0.5*pi) q[14];
rz(0.5*pi) q[13];
rx(0.5*pi) q[14];
rz(0.5*pi) q[13];
rz(0.5*pi) q[14];
rx(0.5*pi) q[13];
rz(0.5*pi) q[13];
rzz(0.5*pi) q[13],q[12];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(0.5*pi) q[12];
rx(0.5*pi) q[13];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(0.5*pi) q[12];
rz(0.5*pi) q[13];
rz(0.5*pi) q[12];
rx(0.5*pi) q[13];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(0.5*pi) q[12];
rz(0.5*pi) q[12];
rzz(0.5*pi) q[12],q[11];
rz(0.5*pi) q[11];
rz(0.5*pi) q[12];
rx(0.5*pi) q[11];
rx(0.5*pi) q[12];
rz(0.5*pi) q[11];
rz(0.5*pi) q[12];
rx(0.5*pi) q[11];
rz(0.5*pi) q[12];
rz(0.5*pi) q[11];
rx(0.5*pi) q[12];
rz(0.5*pi) q[11];
rz(0.5*pi) q[12];
rx(0.5*pi) q[11];
rz(0.5*pi) q[11];
rzz(0.5*pi) q[11],q[10];
rz(0.5*pi) q[10];
rz(0.5*pi) q[11];
rx(0.5*pi) q[10];
rx(0.5*pi) q[11];
rz(0.5*pi) q[10];
rz(0.5*pi) q[11];
rx(0.5*pi) q[10];
rz(0.5*pi) q[11];
rz(0.5*pi) q[10];
rx(0.5*pi) q[11];
rz(0.5*pi) q[10];
rz(0.5*pi) q[11];
rx(0.5*pi) q[10];
rz(0.5*pi) q[10];
rzz(0.5*pi) q[10],q[9];
rz(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[9];
rx(0.5*pi) q[10];
rz(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[9];
rz(0.5*pi) q[10];
rz(0.5*pi) q[9];
rx(0.5*pi) q[10];
rz(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[9];
rz(0.5*pi) q[9];
rzz(0.5*pi) q[9],q[8];
rz(0.5*pi) q[8];
rz(0.5*pi) q[9];
rx(0.5*pi) q[8];
rx(0.5*pi) q[9];
rz(0.5*pi) q[8];
rz(0.5*pi) q[9];
rx(0.5*pi) q[8];
rz(0.5*pi) q[9];
rz(0.5*pi) q[8];
rx(0.5*pi) q[9];
rz(0.5*pi) q[8];
rz(0.5*pi) q[9];
rx(0.5*pi) q[8];
rz(0.5*pi) q[8];
rzz(0.5*pi) q[8],q[7];
rz(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[7];
rx(0.5*pi) q[8];
rz(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[7];
rz(0.5*pi) q[8];
rz(0.5*pi) q[7];
rx(0.5*pi) q[8];
rz(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[7];
rz(0.5*pi) q[7];
rzz(0.5*pi) q[7],q[6];
rz(0.5*pi) q[6];
rz(0.5*pi) q[7];
rx(0.5*pi) q[6];
rx(0.5*pi) q[7];
rz(0.5*pi) q[6];
rz(0.5*pi) q[7];
rx(0.5*pi) q[6];
rz(0.5*pi) q[7];
rz(0.5*pi) q[6];
rx(0.5*pi) q[7];
rz(0.5*pi) q[6];
rz(0.5*pi) q[7];
rx(0.5*pi) q[6];
rz(0.5*pi) q[6];
rzz(0.5*pi) q[6],q[5];
rz(0.5*pi) q[5];
rz(0.5*pi) q[6];
rx(0.5*pi) q[5];
rx(0.5*pi) q[6];
rz(0.5*pi) q[5];
rz(0.5*pi) q[6];
rx(0.5*pi) q[5];
rz(0.5*pi) q[6];
rz(0.5*pi) q[5];
rx(0.5*pi) q[6];
rz(0.5*pi) q[5];
rz(0.5*pi) q[6];
rx(0.5*pi) q[5];
rz(0.5*pi) q[5];
rzz(0.5*pi) q[5],q[4];
rz(0.5*pi) q[4];
rz(0.5*pi) q[5];
rx(0.5*pi) q[4];
rx(0.5*pi) q[5];
rz(0.5*pi) q[4];
rz(0.5*pi) q[5];
rx(0.5*pi) q[4];
rz(0.5*pi) q[5];
rz(0.5*pi) q[4];
rx(0.5*pi) q[5];
rz(0.5*pi) q[4];
rz(0.5*pi) q[5];
rx(0.5*pi) q[4];
rz(0.5*pi) q[4];
rzz(0.5*pi) q[4],q[3];
rz(0.5*pi) q[3];
rz(0.5*pi) q[4];
rx(0.5*pi) q[3];
rx(0.5*pi) q[4];
rz(0.5*pi) q[3];
rz(0.5*pi) q[4];
rx(0.5*pi) q[3];
rz(0.5*pi) q[4];
rz(0.5*pi) q[3];
rx(0.5*pi) q[4];
rz(0.5*pi) q[3];
rz(0.5*pi) q[4];
rx(0.5*pi) q[3];
rz(0.5*pi) q[3];
rzz(0.5*pi) q[3],q[2];
rz(0.5*pi) q[2];
rz(0.5*pi) q[3];
rx(0.5*pi) q[2];
rx(0.5*pi) q[3];
rz(0.5*pi) q[2];
rz(0.5*pi) q[3];
rx(0.5*pi) q[2];
rz(0.5*pi) q[3];
rz(0.5*pi) q[2];
rx(0.5*pi) q[3];
rz(0.5*pi) q[2];
rz(0.5*pi) q[3];
rx(0.5*pi) q[2];
rz(0.5*pi) q[2];
rzz(0.5*pi) q[2],q[1];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[1];
rz(0.5*pi) q[2];
rz(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[1];
rz(0.5*pi) q[1];
rzz(0.5*pi) q[1],q[0];
rz(0.5*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[0];
rz(0.5*pi) q[1];
rz(1.0*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[1];
rz(0.5*pi) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15],q[16],q[17],q[18],q[19],q[20],q[21],q[22],q[23],q[24],q[25],q[26],q[27],q[28],q[29],q[30],q[31],q[32],q[33],q[34],q[35],q[36],q[37],q[38],q[39],q[40],q[41],q[42],q[43],q[44],q[45],q[46],q[47],q[48],q[49],q[50],q[51],q[52],q[53],q[54],q[55],q[56],q[57],q[58],q[59],q[60],q[61],q[62],q[63],q[64],q[65],q[66],q[67],q[68],q[69],q[70],q[71],q[72],q[73],q[74],q[75],q[76],q[77],q[78],q[79],q[80],q[81],q[82],q[83],q[84],q[85],q[86],q[87],q[88],q[89],q[90],q[91],q[92],q[93],q[94],q[95],q[96],q[97],q[98],q[99];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
measure q[4] -> meas[4];
measure q[5] -> meas[5];
measure q[6] -> meas[6];
measure q[7] -> meas[7];
measure q[8] -> meas[8];
measure q[9] -> meas[9];
measure q[10] -> meas[10];
measure q[11] -> meas[11];
measure q[12] -> meas[12];
measure q[13] -> meas[13];
measure q[14] -> meas[14];
measure q[15] -> meas[15];
measure q[16] -> meas[16];
measure q[17] -> meas[17];
measure q[18] -> meas[18];
measure q[19] -> meas[19];
measure q[20] -> meas[20];
measure q[21] -> meas[21];
measure q[22] -> meas[22];
measure q[23] -> meas[23];
measure q[24] -> meas[24];
measure q[25] -> meas[25];
measure q[26] -> meas[26];
measure q[27] -> meas[27];
measure q[28] -> meas[28];
measure q[29] -> meas[29];
measure q[30] -> meas[30];
measure q[31] -> meas[31];
measure q[32] -> meas[32];
measure q[33] -> meas[33];
measure q[34] -> meas[34];
measure q[35] -> meas[35];
measure q[36] -> meas[36];
measure q[37] -> meas[37];
measure q[38] -> meas[38];
measure q[39] -> meas[39];
measure q[40] -> meas[40];
measure q[41] -> meas[41];
measure q[42] -> meas[42];
measure q[43] -> meas[43];
measure q[44] -> meas[44];
measure q[45] -> meas[45];
measure q[46] -> meas[46];
measure q[47] -> meas[47];
measure q[48] -> meas[48];
measure q[49] -> meas[49];
measure q[50] -> meas[50];
measure q[51] -> meas[51];
measure q[52] -> meas[52];
measure q[53] -> meas[53];
measure q[54] -> meas[54];
measure q[55] -> meas[55];
measure q[56] -> meas[56];
measure q[57] -> meas[57];
measure q[58] -> meas[58];
measure q[59] -> meas[59];
measure q[60] -> meas[60];
measure q[61] -> meas[61];
measure q[62] -> meas[62];
measure q[63] -> meas[63];
measure q[64] -> meas[64];
measure q[65] -> meas[65];
measure q[66] -> meas[66];
measure q[67] -> meas[67];
measure q[68] -> meas[68];
measure q[69] -> meas[69];
measure q[70] -> meas[70];
measure q[71] -> meas[71];
measure q[72] -> meas[72];
measure q[73] -> meas[73];
measure q[74] -> meas[74];
measure q[75] -> meas[75];
measure q[76] -> meas[76];
measure q[77] -> meas[77];
measure q[78] -> meas[78];
measure q[79] -> meas[79];
measure q[80] -> meas[80];
measure q[81] -> meas[81];
measure q[82] -> meas[82];
measure q[83] -> meas[83];
measure q[84] -> meas[84];
measure q[85] -> meas[85];
measure q[86] -> meas[86];
measure q[87] -> meas[87];
measure q[88] -> meas[88];
measure q[89] -> meas[89];
measure q[90] -> meas[90];
measure q[91] -> meas[91];
measure q[92] -> meas[92];
measure q[93] -> meas[93];
measure q[94] -> meas[94];
measure q[95] -> meas[95];
measure q[96] -> meas[96];
measure q[97] -> meas[97];
measure q[98] -> meas[98];
measure q[99] -> meas[99];

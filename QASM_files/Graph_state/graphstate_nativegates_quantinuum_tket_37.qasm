// Benchmark was created by MQT Bench on 2024-03-17
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 1.1.0
// TKET version: 1.25.0
// Used Gate Set: ['rzz', 'rz', 'ry', 'rx', 'measure', 'barrier']

OPENQASM 2.0;
include "qelib1.inc";

qreg q[37];
creg meas[37];
rx(1.0*pi) q[0];
rz(3.0*pi) q[1];
rx(1.0*pi) q[2];
rz(3.0*pi) q[3];
rx(1.0*pi) q[4];
rz(3.0*pi) q[5];
rx(1.0*pi) q[6];
rz(3.0*pi) q[7];
rx(1.0*pi) q[8];
rz(3.0*pi) q[9];
rx(1.0*pi) q[10];
rz(3.0*pi) q[11];
rz(3.0*pi) q[12];
rx(1.0*pi) q[13];
rz(3.0*pi) q[14];
rx(1.0*pi) q[15];
rz(3.0*pi) q[16];
rz(3.0*pi) q[17];
rz(3.0*pi) q[18];
rx(1.0*pi) q[19];
rz(3.0*pi) q[20];
rz(3.0*pi) q[21];
rz(3.0*pi) q[22];
rz(3.0*pi) q[23];
rz(3.0*pi) q[24];
rx(1.0*pi) q[25];
rz(3.0*pi) q[26];
rz(3.0*pi) q[27];
rz(3.0*pi) q[28];
rz(3.0*pi) q[29];
rz(3.0*pi) q[30];
rz(3.0*pi) q[31];
rz(3.0*pi) q[32];
rz(3.0*pi) q[33];
rz(3.0*pi) q[34];
rz(3.0*pi) q[35];
rz(3.0*pi) q[36];
rz(0.5*pi) q[0];
rx(1.0*pi) q[1];
rz(0.5*pi) q[2];
rx(1.0*pi) q[3];
rz(0.5*pi) q[4];
rx(1.0*pi) q[5];
rz(0.5*pi) q[6];
rx(1.0*pi) q[7];
rz(0.5*pi) q[8];
rx(1.0*pi) q[9];
rz(0.5*pi) q[10];
rx(1.0*pi) q[11];
rx(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(1.0*pi) q[14];
rz(0.5*pi) q[15];
rx(0.5*pi) q[16];
rx(1.0*pi) q[17];
rx(1.0*pi) q[18];
rz(0.5*pi) q[19];
rx(1.0*pi) q[20];
rx(0.5*pi) q[21];
rx(1.0*pi) q[22];
rx(0.5*pi) q[23];
rx(1.0*pi) q[24];
rz(0.5*pi) q[25];
rx(0.5*pi) q[26];
rx(1.0*pi) q[27];
rx(1.0*pi) q[28];
rx(1.0*pi) q[29];
rx(1.0*pi) q[30];
rx(0.5*pi) q[31];
rx(1.0*pi) q[32];
rx(0.5*pi) q[33];
rx(0.5*pi) q[34];
rx(0.5*pi) q[35];
rx(0.5*pi) q[36];
rx(0.5*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[3];
rx(0.5*pi) q[4];
rz(0.5*pi) q[5];
rx(0.5*pi) q[6];
rz(0.5*pi) q[7];
rx(0.5*pi) q[8];
rz(0.5*pi) q[9];
rx(0.5*pi) q[10];
rz(0.5*pi) q[11];
rz(0.5*pi) q[12];
rx(0.5*pi) q[13];
rz(0.5*pi) q[14];
rx(0.5*pi) q[15];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rz(0.5*pi) q[18];
rx(0.5*pi) q[19];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rx(0.5*pi) q[25];
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
rz(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[3];
rz(0.5*pi) q[4];
rx(0.5*pi) q[5];
rz(0.5*pi) q[6];
rx(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[11];
rx(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(0.5*pi) q[14];
rz(0.5*pi) q[15];
rx(0.5*pi) q[16];
rx(0.5*pi) q[17];
rx(0.5*pi) q[18];
rz(0.5*pi) q[19];
rx(0.5*pi) q[20];
rx(0.5*pi) q[21];
rx(0.5*pi) q[22];
rx(0.5*pi) q[23];
rx(0.5*pi) q[24];
rz(0.5*pi) q[25];
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
rz(0.5*pi) q[1];
rz(0.5*pi) q[3];
rz(0.5*pi) q[5];
rz(0.5*pi) q[7];
rz(0.5*pi) q[9];
rz(0.5*pi) q[11];
rz(0.5*pi) q[12];
rz(0.5*pi) q[14];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rz(0.5*pi) q[18];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
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
rzz(0.5*pi) q[0],q[1];
rzz(0.5*pi) q[2],q[3];
rzz(0.5*pi) q[4],q[5];
rzz(0.5*pi) q[6],q[7];
rzz(0.5*pi) q[8],q[9];
rzz(0.5*pi) q[10],q[11];
rzz(0.5*pi) q[13],q[14];
rzz(0.5*pi) q[19],q[20];
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
rx(0.5*pi) q[13];
rx(0.5*pi) q[14];
rx(0.5*pi) q[19];
rx(0.5*pi) q[20];
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
rz(0.5*pi) q[13];
rz(0.5*pi) q[14];
rz(0.5*pi) q[19];
rz(0.5*pi) q[20];
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
rx(0.5*pi) q[13];
rx(0.5*pi) q[14];
rx(0.5*pi) q[19];
rx(0.5*pi) q[20];
rz(0.5*pi) q[0];
rz(1.0*pi) q[1];
rz(0.5*pi) q[2];
rz(1.0*pi) q[3];
rz(0.5*pi) q[4];
rz(1.0*pi) q[5];
rz(0.5*pi) q[6];
rz(1.0*pi) q[7];
rz(0.5*pi) q[8];
rz(1.0*pi) q[9];
rz(0.5*pi) q[10];
rz(1.0*pi) q[11];
rz(0.5*pi) q[13];
rz(1.0*pi) q[14];
rz(0.5*pi) q[19];
rz(1.0*pi) q[20];
rx(0.5*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[3];
rx(0.5*pi) q[4];
rz(0.5*pi) q[5];
rx(0.5*pi) q[6];
rz(0.5*pi) q[7];
rx(0.5*pi) q[8];
rz(0.5*pi) q[9];
rx(0.5*pi) q[10];
rz(0.5*pi) q[11];
rx(0.5*pi) q[13];
rz(0.5*pi) q[14];
rx(0.5*pi) q[19];
rz(0.5*pi) q[20];
rz(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[3];
rz(0.5*pi) q[4];
rx(0.5*pi) q[5];
rz(0.5*pi) q[6];
rx(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[11];
rz(0.5*pi) q[13];
rx(0.5*pi) q[14];
rz(0.5*pi) q[19];
rx(0.5*pi) q[20];
rzz(0.5*pi) q[0],q[33];
rz(0.5*pi) q[1];
rzz(0.5*pi) q[2],q[29];
rz(0.5*pi) q[3];
rzz(0.5*pi) q[4],q[12];
rz(0.5*pi) q[5];
rzz(0.5*pi) q[6],q[26];
rz(0.5*pi) q[7];
rzz(0.5*pi) q[8],q[17];
rz(0.5*pi) q[9];
rzz(0.5*pi) q[10],q[27];
rz(0.5*pi) q[11];
rzz(0.5*pi) q[13],q[36];
rz(0.5*pi) q[14];
rz(0.5*pi) q[20];
rx(0.5*pi) q[0];
rzz(0.5*pi) q[1],q[16];
rx(0.5*pi) q[2];
rzz(0.5*pi) q[3],q[21];
rx(0.5*pi) q[4];
rx(0.5*pi) q[6];
rzz(0.5*pi) q[7],q[34];
rx(0.5*pi) q[8];
rzz(0.5*pi) q[9],q[35];
rx(0.5*pi) q[10];
rzz(0.5*pi) q[11],q[18];
rx(0.5*pi) q[12];
rx(0.5*pi) q[13];
rzz(0.5*pi) q[14],q[24];
rx(0.5*pi) q[17];
rzz(0.5*pi) q[20],q[22];
rx(0.5*pi) q[26];
rx(0.5*pi) q[27];
rx(0.5*pi) q[29];
rx(0.5*pi) q[33];
rx(0.5*pi) q[36];
rz(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[3];
rz(0.5*pi) q[4];
rz(0.5*pi) q[6];
rx(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[11];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(0.5*pi) q[14];
rx(0.5*pi) q[16];
rz(0.5*pi) q[17];
rx(0.5*pi) q[18];
rx(0.5*pi) q[20];
rx(0.5*pi) q[21];
rx(0.5*pi) q[22];
rx(0.5*pi) q[24];
rz(0.5*pi) q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[29];
rz(0.5*pi) q[33];
rx(0.5*pi) q[34];
rx(0.5*pi) q[35];
rz(0.5*pi) q[36];
rx(0.5*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[3];
rx(0.5*pi) q[4];
rx(0.5*pi) q[6];
rz(0.5*pi) q[7];
rx(0.5*pi) q[8];
rz(0.5*pi) q[9];
rx(0.5*pi) q[10];
rz(0.5*pi) q[11];
rx(0.5*pi) q[12];
rx(0.5*pi) q[13];
rz(0.5*pi) q[14];
rz(0.5*pi) q[16];
rx(0.5*pi) q[17];
rz(0.5*pi) q[18];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[24];
rx(0.5*pi) q[26];
rx(0.5*pi) q[27];
rx(0.5*pi) q[29];
rx(0.5*pi) q[33];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rx(0.5*pi) q[36];
rz(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[3];
rz(0.5*pi) q[4];
rz(0.5*pi) q[6];
rx(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[11];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(0.5*pi) q[14];
rx(0.5*pi) q[16];
rz(1.0*pi) q[17];
rx(0.5*pi) q[18];
rx(0.5*pi) q[20];
rx(0.5*pi) q[21];
rx(0.5*pi) q[22];
rx(0.5*pi) q[24];
rz(0.5*pi) q[26];
rz(1.0*pi) q[27];
rz(1.0*pi) q[29];
rz(0.5*pi) q[33];
rx(0.5*pi) q[34];
rx(0.5*pi) q[35];
rz(0.5*pi) q[36];
rx(0.5*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[3];
rx(0.5*pi) q[4];
rx(0.5*pi) q[6];
rz(0.5*pi) q[7];
rx(0.5*pi) q[8];
rz(0.5*pi) q[9];
rx(0.5*pi) q[10];
rz(0.5*pi) q[11];
rx(0.5*pi) q[12];
rx(0.5*pi) q[13];
rz(0.5*pi) q[14];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rz(1.0*pi) q[18];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rz(1.0*pi) q[22];
rz(1.0*pi) q[24];
rx(0.5*pi) q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[29];
rx(0.5*pi) q[33];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rx(0.5*pi) q[36];
rz(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[3];
rz(0.5*pi) q[4];
rz(0.5*pi) q[6];
rx(0.5*pi) q[7];
rz(0.5*pi) q[8];
rx(0.5*pi) q[9];
rz(0.5*pi) q[10];
rx(0.5*pi) q[11];
rz(0.5*pi) q[12];
rz(0.5*pi) q[13];
rx(0.5*pi) q[14];
rx(0.5*pi) q[16];
rx(0.5*pi) q[17];
rz(0.5*pi) q[18];
rx(0.5*pi) q[20];
rx(0.5*pi) q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[24];
rz(0.5*pi) q[26];
rx(0.5*pi) q[27];
rx(0.5*pi) q[29];
rz(0.5*pi) q[33];
rx(0.5*pi) q[34];
rx(0.5*pi) q[35];
rz(0.5*pi) q[36];
rz(0.5*pi) q[1];
rz(0.5*pi) q[3];
rzz(0.5*pi) q[5],q[12];
rz(0.5*pi) q[7];
rz(0.5*pi) q[9];
rz(0.5*pi) q[11];
rz(0.5*pi) q[14];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rx(0.5*pi) q[18];
rz(0.5*pi) q[20];
rz(0.5*pi) q[21];
rx(0.5*pi) q[22];
rx(0.5*pi) q[24];
rzz(0.5*pi) q[25],q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[29];
rz(0.5*pi) q[34];
rz(0.5*pi) q[35];
rx(0.5*pi) q[5];
rx(0.5*pi) q[12];
rzz(0.5*pi) q[15],q[16];
rzz(0.5*pi) q[17],q[35];
rz(0.5*pi) q[18];
rzz(0.5*pi) q[19],q[21];
rz(0.5*pi) q[22];
rz(0.5*pi) q[24];
rx(0.5*pi) q[25];
rx(0.5*pi) q[26];
rzz(0.5*pi) q[27],q[32];
rzz(0.5*pi) q[29],q[30];
rz(0.5*pi) q[5];
rz(0.5*pi) q[12];
rx(0.5*pi) q[15];
rx(0.5*pi) q[16];
rx(0.5*pi) q[17];
rzz(0.5*pi) q[18],q[23];
rx(0.5*pi) q[19];
rx(0.5*pi) q[21];
rzz(0.5*pi) q[24],q[36];
rz(0.5*pi) q[25];
rz(0.5*pi) q[26];
rx(0.5*pi) q[27];
rx(0.5*pi) q[29];
rx(0.5*pi) q[30];
rx(0.5*pi) q[32];
rx(0.5*pi) q[35];
rx(0.5*pi) q[5];
rx(0.5*pi) q[12];
rz(0.5*pi) q[15];
rz(0.5*pi) q[16];
rz(0.5*pi) q[17];
rx(0.5*pi) q[18];
rz(0.5*pi) q[19];
rz(0.5*pi) q[21];
rx(0.5*pi) q[23];
rx(0.5*pi) q[24];
rx(0.5*pi) q[25];
rx(0.5*pi) q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[29];
rz(0.5*pi) q[30];
rz(0.5*pi) q[32];
rz(0.5*pi) q[35];
rx(0.5*pi) q[36];
rz(0.5*pi) q[5];
rz(3.5*pi) q[12];
rx(0.5*pi) q[15];
rx(0.5*pi) q[16];
rx(0.5*pi) q[17];
rz(0.5*pi) q[18];
rx(0.5*pi) q[19];
rx(0.5*pi) q[21];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rz(0.5*pi) q[25];
rz(3.5*pi) q[26];
rx(0.5*pi) q[27];
rx(0.5*pi) q[29];
rx(0.5*pi) q[30];
rx(0.5*pi) q[32];
rx(0.5*pi) q[35];
rz(0.5*pi) q[36];
rx(0.5*pi) q[5];
rx(0.5*pi) q[12];
rz(0.5*pi) q[15];
rz(3.5*pi) q[16];
rz(0.5*pi) q[17];
rx(0.5*pi) q[18];
rz(0.5*pi) q[19];
rz(3.5*pi) q[21];
rx(0.5*pi) q[23];
rx(0.5*pi) q[24];
rx(0.5*pi) q[25];
rx(0.5*pi) q[26];
rz(0.5*pi) q[27];
rz(0.5*pi) q[29];
rz(1.0*pi) q[30];
rz(1.0*pi) q[32];
rz(3.5*pi) q[35];
rx(0.5*pi) q[36];
rz(0.5*pi) q[5];
rx(0.5*pi) q[15];
rx(0.5*pi) q[16];
rx(0.5*pi) q[17];
rz(0.5*pi) q[18];
rx(0.5*pi) q[19];
rx(0.5*pi) q[21];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rz(0.5*pi) q[25];
rx(0.5*pi) q[27];
rx(0.5*pi) q[29];
rz(0.5*pi) q[30];
rz(0.5*pi) q[32];
rx(0.5*pi) q[35];
rz(3.5*pi) q[36];
rz(0.5*pi) q[15];
rz(0.5*pi) q[17];
rx(0.5*pi) q[18];
rz(0.5*pi) q[19];
rx(0.5*pi) q[23];
rx(0.5*pi) q[24];
rzz(0.5*pi) q[25],q[34];
rz(0.5*pi) q[27];
rz(0.5*pi) q[29];
rx(0.5*pi) q[30];
rx(0.5*pi) q[32];
rx(0.5*pi) q[36];
rzz(0.5*pi) q[15],q[28];
rz(0.5*pi) q[18];
rz(0.5*pi) q[23];
rz(0.5*pi) q[24];
rx(0.5*pi) q[25];
rz(0.5*pi) q[30];
rz(0.5*pi) q[32];
rx(0.5*pi) q[34];
rx(0.5*pi) q[15];
rzz(0.5*pi) q[22],q[23];
rz(0.5*pi) q[25];
rx(0.5*pi) q[28];
rzz(0.5*pi) q[32],q[33];
rz(0.5*pi) q[34];
rz(0.5*pi) q[15];
rx(0.5*pi) q[22];
rx(0.5*pi) q[23];
rx(0.5*pi) q[25];
rz(0.5*pi) q[28];
rx(0.5*pi) q[32];
rx(0.5*pi) q[33];
rx(0.5*pi) q[34];
rx(0.5*pi) q[15];
rz(0.5*pi) q[22];
rz(0.5*pi) q[23];
rz(0.5*pi) q[25];
rx(0.5*pi) q[28];
rz(0.5*pi) q[32];
rz(0.5*pi) q[33];
rz(3.5*pi) q[34];
rz(0.5*pi) q[15];
rx(0.5*pi) q[22];
rx(0.5*pi) q[23];
rx(0.5*pi) q[25];
rz(1.0*pi) q[28];
rx(0.5*pi) q[32];
rx(0.5*pi) q[33];
rx(0.5*pi) q[34];
rx(0.5*pi) q[15];
rz(0.5*pi) q[22];
rz(3.5*pi) q[23];
rz(0.5*pi) q[25];
rz(0.5*pi) q[28];
rz(0.5*pi) q[32];
rz(3.5*pi) q[33];
rz(0.5*pi) q[15];
rx(0.5*pi) q[22];
rx(0.5*pi) q[23];
rx(0.5*pi) q[28];
rx(0.5*pi) q[32];
rx(0.5*pi) q[33];
rz(0.5*pi) q[22];
rz(0.5*pi) q[28];
rz(0.5*pi) q[32];
rzz(0.5*pi) q[28],q[31];
rx(0.5*pi) q[28];
rx(0.5*pi) q[31];
rz(0.5*pi) q[28];
rz(0.5*pi) q[31];
rx(0.5*pi) q[28];
rx(0.5*pi) q[31];
rz(0.5*pi) q[28];
rz(0.5*pi) q[31];
rx(0.5*pi) q[28];
rx(0.5*pi) q[31];
rz(0.5*pi) q[28];
rz(0.5*pi) q[31];
rzz(0.5*pi) q[30],q[31];
rx(0.5*pi) q[30];
rx(0.5*pi) q[31];
rz(0.5*pi) q[30];
rz(0.5*pi) q[31];
rx(0.5*pi) q[30];
rx(0.5*pi) q[31];
rz(0.5*pi) q[30];
rz(3.5*pi) q[31];
rx(0.5*pi) q[30];
rx(0.5*pi) q[31];
rz(0.5*pi) q[30];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15],q[16],q[17],q[18],q[19],q[20],q[21],q[22],q[23],q[24],q[25],q[26],q[27],q[28],q[29],q[30],q[31],q[32],q[33],q[34],q[35],q[36];
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

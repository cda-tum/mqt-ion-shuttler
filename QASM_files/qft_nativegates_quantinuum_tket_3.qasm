// Benchmark was created by MQT Bench on 2023-06-29
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: v1.0.0
// TKET version: 1.16.0
// Used Gate Set: ['rz', 'ry', 'rx', 'rzz', 'measure']

OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];
creg c[3];
creg meas[3];
rz(0.8750000000000009*pi) q[0];
rx(1.5*pi) q[1];
rz(1.0*pi) q[2];
rx(2.5*pi) q[0];
rz(1.5*pi) q[1];
rx(1.25*pi) q[2];
rz(0.5*pi) q[0];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rz(0.5*pi) q[0];
rx(0.5*pi) q[1];
rx(0.5*pi) q[2];
rx(0.5*pi) q[0];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rz(0.5*pi) q[0];
rzz(0.25*pi) q[2],q[1];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[1];
rz(0.5*pi) q[2];
rz(3.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(1.75*pi) q[1];
rz(0.5*pi) q[2];
rzz(0.125*pi) q[2],q[0];
rz(0.5*pi) q[1];
rz(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[2];
rx(0.5*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[2];
rz(0.5*pi) q[0];
rz(0.5*pi) q[2];
rz(0.5*pi) q[0];
rz(3.5*pi) q[2];
rx(0.5*pi) q[0];
rx(0.5*pi) q[2];
rz(0.5*pi) q[0];
rz(2.1249999999999996*pi) q[2];
rzz(0.25*pi) q[1],q[0];
rz(0.5*pi) q[0];
rz(0.5*pi) q[1];
rx(0.5*pi) q[0];
rx(0.5*pi) q[1];
rz(0.5*pi) q[0];
rz(0.5*pi) q[1];
rz(3.5*pi) q[1];
rx(0.5*pi) q[1];
rz(3.75*pi) q[1];

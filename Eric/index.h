/*
	Created by Guanzheng Xu
*/

#ifndef INDEX_H
#define INDEX_H

const int PQ       = 1;   // load
const int PV       = 2;   // generator
const int REF      = 3;   // reference generator
const int NONE     = 4;

//Bus Data Format
const int BUS_I       = 0;    // bus number (positive integer)
const int BUS_TYPE    = 1;    // bus type (1-PQ bus, 2-PV bus, 3-reference bus, 4-isolated bus)
const int PD          = 2;    // Pd, real power demand (MW)
const int QD          = 3;    // Qd, reactive power demand (MVAr)
const int GS          = 4;    // Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
const int BS          = 5;    // Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
const int BUS_AREA    = 6;    // area number (positive integer)
const int VM          = 7;    // Vm, voltage magnitude (p.u.)
const int VA          = 8;    // Va, voltage angle (degrees)
const int BASE_KV     = 9;    // baseKV, base voltage (kV)
const int ZONE        = 10;   // zone, loss zone (positive integer)
const int VMAX        = 11;   // maxVm, maximum voltage magnitude (p.u.)  (not in PTI format)
const int VMIN        = 12;   // minVm, minimum voltage magnitude (p.u.)  (not in PTI format)

//Generator Data Format
const int GEN_BUS     = 0;    // bus number
const int PG          = 1;    // Pg, real power output (MW)
const int QG          = 2;    // Qg, reactive power output (MVAr)
const int QMAX        = 3;    // Qmax, maximum reactive power output (MVAr)
const int QMIN        = 4;    // Qmin, minimum reactive power output (MVAr)
const int VG          = 5;    // Vg, voltage magnitude setpoint (p.u.)
const int MBASE       = 6;    // mBase, total MVA base of this machine, defaults to baseMVA
const int GEN_STATUS  = 7;    // status, > 0 - machine in service, <= 0 - machine out of service
const int PMAX        = 8;    // Pmax, maximum real power output (MW)
const int PMIN        = 9;    // Pmin, minimum real power output (MW)
const int PC1         = 10;   // Pc1, lower real power output of PQ capability curve (MW)
const int PC2         = 11;   // Pc2, upper real power output of PQ capability curve (MW)
const int QC1MIN      = 12;   // Qc1min, minimum reactive power output at Pc1 (MVAr)
const int QC1MAX      = 13;   // Qc1max, maximum reactive power output at Pc1 (MVAr)
const int QC2MIN      = 14;   // Qc2min, minimum reactive power output at Pc2 (MVAr)
const int QC2MAX      = 15;   // Qc2max, maximum reactive power output at Pc2 (MVAr)
const int RAMP_AGC    = 16;   // ramp rate for load following/AGC (MW/min)
const int RAMP_10     = 17;   // ramp rate for 10 minute reserves (MW)
const int RAMP_30     = 18;   // ramp rate for 30 minute reserves (MW)
const int RAMP_Q      = 19;   // ramp rate for reactive power (2 sec timescale) (MVAr/min)
const int APF         = 20;   // APF, area participation factor

//Branch Data Format
const int F_BUS       = 0;    // f, from bus number
const int T_BUS       = 1;    // t, to bus number
const int BR_R        = 2;    // r, resistance (p.u.)
const int BR_X        = 3;    // x, reactance (p.u.)
const int BR_B        = 4;    // b, total line charging susceptance (p.u.)
const int RATE_A      = 5;    // rateA, MVA rating A (long term rating)
const int RATE_B      = 6;    // rateB, MVA rating B (short term rating)
const int RATE_C      = 7;    // rateC, MVA rating C (emergency rating)
const int TAP         = 8;    // ratio, transformer off nominal turns ratio
const int SHIFT       = 9;    // angle, transformer phase shift angle (degrees)
const int BR_STATUS   = 10;   // initial branch status, 1-in service, 0-out of service
const int ANGMIN      = 11;   // minimum angle difference, angle(Vf) - angle(Vt) (degrees)
const int ANGMAX      = 12;   // maximum angle difference, angle(Vf) - angle(Vt) (degrees)

#endif
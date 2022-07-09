#include <iostream>
#include <cstring>
#include <math.h>
#include <vector>
#include <particles.h>

using namespace std;
double const MMM = 1;
double Cm = 0.09;
int NPART = 500000;
int id, LX, LY, LZ, LXf, LYf, LZf, iwall, dum, dumCh, dumx, iter = 0, cell, iterStep, LostParticles = 0;
double Xp1, Yp1, Zp1, Up1, Vp1, Wp1, PDIAM, dFP, Ckg, VolCell;
double Xf1, Yf1, Zf1, Uf1, Vf1, Wf1, Xfm1, Yfm1, Zfm1, Ufm1, Vfm1, Wfm1;
double UF, VF, WF, UFM = 0., VFM = 0., WFM = 0., UFPM, VFPM, WFPM, UP1, VP1, WP1;
double Kturb, Eturb, Le, VC;
double DSEED, GGG, x;
string IDIR, IDIR2, IDIRfluidPar;
double xTline = 0., yTline = 0., zTline = 0.;
vector<vector<vector<int>>> ii;
vector<vector<vector<double>>> Uflow, Vflow, Wflow, ED, TE, CONC1m, Concentration;
vector<vector<vector<int>>> iwallVec;
vector<double> XPpast, YPpast, ZPpast,
    UPpast, VPpast, WPpast, pdiami;
vector<double> XPFpast, YPFpast, ZPFpast,
    UPFpast, VPFpast, WPFpast;
vector<string> IDpast, IDFpast;
vector<double> Xgrid, Ygrid, Zgrid, DX, DY, DZ, Xload, Yload, Zload, SDU, SNS, SEW, XloadNI, YloadNJ, ZloadNK, CONV;
vector<double> Uflow1d, Uflow1dd, Vflow1d, Vflow1dd, Wflow1d, Wflow1dd, UflowLoad, VflowLoad, WflowLoad, TEflow1dd, EDflow1dd, EDflowLoad, TEflowLoad;
vector<int> Iwallload, iwallV1, CELLpast, CELLFpast, LXpast, LYpast, LZpast, LXFpast, LYFpast, LZFpast;
vector<double> PINAX;
string str, str2;
vector<particles> particle, FluidParticle;
//// currierInterp Function Declaration
int LX1, LY1, LZ1, LX2, LY2, LZ2, LY3, LY4, LY5, LZ3, LZ4, LZ5, LX3, LX4, LX5;
int pC, pN, pS, pE, pW, pD, pU, p5;
double QN, QS, QC, QE, QW, QD, QU, Q1, Q2, Q3, Q4, Q5, Q;
double POL1, POL2, POL3, POL4, POL5, POL;
double Xmax, Xmin, Ymax, Ymin, Zmax, Zmin, xp, yp, zp;

string idNum;
ofstream Results, tecRes2, iwallPri, tecResn, Traj1, Traj2, Traj3, Traj4, Traj5, Traj6, Traj7, Traj8, myNODES, Stk_num, testCase, CONC_, IwallTest2, SDWtest;
ifstream CCC;

#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "particles.h"
#include "Grid.h"

using namespace std;

struct GRID
{
    int NI, NJ, NK;
    vector<double> X, Y, Z;
};

class iodata
{
private:
    GRID grid;
    vector<particles> particle;
    vector<Nodes> myNodes;
    int NI, NJ, NK;
    int N, L, ISTART0, IEND0, JSTART0, JEND0, KSTART0, KEND0, LLL;
    double den, visc;
    double KLM, PDEN0, PANU, PAE, SUDEN, SUNU, SUE, SEW, DGAM, YELAST, DMA, DNA, DPA, NPDIA,
        ZEMIS1, ZEMIS2, IEMIS, JEMIS, TURINTU, TURINTV, TURINTW, THCKP, THCKG, RGAS, SUDEN2, SUNU2, SUE2,
        DGAM2, A0, GRAV, DisCoef;
    int MINTSPSV, MAXTSPSV, MAXTSPENTRA;
    vector<double> prop, X, XC, DXEP, DXPW,
        WFE, WFW, XPEF, XEEF, XWWF, XWPF;
    vector<double> Y, YC, DYPS, DYNP, SNS,
        WFN, WFS, YPNF, YNNF, YSSF, YSPF;
    vector<double> Z, ZC, DZPU, DZDP, SDU,
        WFD, WFU, ZPDF, ZDDF, ZUUF, ZUPF;
    vector<double> pdiami, PMFRA;
    vector<vector<vector<double>>> U, V, W, TE, ED, TEMP, UPAST, VPAST, WPAST, TEPAST, EDPAST, TEMPAST, VIS, P, DEN,
        DENPAST, CN, CE, CW, CS, CD, CU, DPU, DPV, DPW, THCK, CCP, CMU;
    vector<double> Xgrid, Ygrid, Zgrid, U1d, V1d, W1d, TE1d, ED1d;

public:
    iodata();
    iodata(int NI, int NJ, int NK)
    {
        this->NI = NI;
        this->NJ = NJ;
        this->NK = NK;
    }
    void resizeVec();
    void readVel(string fileName);
    void readFluidPro(string fileName);
    void readPartPro(string fileName);
    int getZEMIS1() { return ZEMIS1; }
    int getZEMIS2() { return ZEMIS2; }
    int getIEMIS() { return IEMIS; }
    int getJEMIS() { return JEMIS; }

    //////////////////////////////////////////////////////
    vector<vector<vector<double>>> getU() { return U; }
    vector<vector<vector<double>>> getV() { return V; }
    vector<vector<vector<double>>> getW() { return W; }
    vector<vector<vector<double>>> getTE() { return TE; }
    vector<vector<vector<double>>> getED() { return ED; }

    vector<double> getU1d() { return U1d; }
    vector<double> getV1d() { return V1d; }
    vector<double> getW1d() { return W1d; }
    vector<double> getTE1d() { return TE1d; }
    vector<double> getED1d() { return ED1d; }

    ////////////////////////////////////////////////////////

    double getVisc() { return visc; }
    double getDen() { return den; }
    vector<double> getPdiami() { return pdiami; }
    double getNPDIA() { return NPDIA; }
    int getMAXTSPENTRA() { return MAXTSPENTRA; }
    int getMAXTSPSV() { return MAXTSPSV; }
    double getPDEN0() { return PDEN0; }
    double getGRAV() { return GRAV; }
    double getDisCoef() { return DisCoef; }
    int getLLL() { return LLL; }
    double getTURINTU() { return TURINTU; }
    double getTURINTV() { return TURINTV; }
    double getTURINTW() { return TURINTW; }
    ~iodata();
};
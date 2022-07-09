#include <vector>
#include <cstring>
#include <math.h>
#include "iodata.h"
#include "Grid.h"
#include "declare.h"
#include "GNodes.h"

double step = 0.01;
void ode(double UF, double VF, double WF, vector<struct particles> &particle, int count,
         double den, double VISCS, double GRAV, double PDEN);
void odeFluidP(double UF, double VF, double WF, vector<struct particles> &FluidParticle, int count,
               double den, double VISCS, double GRAV, double PDEN);
void currierInterp(vector<struct Nodes> &myNodes, vector<vector<vector<int>>> &ii, vector<std::vector<std::vector<double>>> &Uflow, vector<std::vector<std::vector<double>>> &Vflow,
                   vector<std::vector<std::vector<double>>> &Wflow, double Xp1, double Yp1, double Zp1, int &LX, int &LY, int &LZ, double &UFM, double &VFM, double &WFM);
void tline(int id, vector<struct Nodes> &myNodes, vector<vector<vector<int>>> &ii, vector<struct particles> &particle,
           double Xm1, double Ym1, double Zm1, string &IDIR, double &xb, double &yb, double &zb);
void ParticleSP(int id, vector<struct Gnodes> &myNodes, vector<vector<vector<int>>> &ii, vector<struct particles> &particle, string &IDIR, vector<Nodes> &iiNodes);
void FPlost(int id, vector<struct particles> &Fluidparticle, vector<struct particles> &particle, vector<struct Gnodes> &myNodes,
            vector<std::vector<std::vector<double>>> &Uflow, vector<std::vector<std::vector<double>>> &Vflow, vector<std::vector<std::vector<double>>> &Wflow, int NI, int NJ, int NK);
void iwallCheck(int id, vector<struct Gnodes> &myNodes, vector<struct particles> &particle, string IDIR, double &xTline, double &yTline, double &zTline, vector<int> &LXpast, vector<int> &LYpast, vector<int> &LZpast, vector<int> &CELLpast);

template <class T>
T min(T &a, T &b)
{
    return (a < b ? a : b);
}

bool isNotNumber(const string &str)
{
    for (char const &c : str)
    {
        if (std::isdigit(c) == 0)
            return true;
    }
    return false;
}


//// Sovle the euler scheme for Particles
void ode(double UF, double VF, double WF, vector<struct particles> &particle, int count,
         double den, double VISCS, double GRAV, double PDEN)
{
    double DR = 0.1, UPM1, VPM1, WPM1;
    double B, REYND, CD, A, GR;
    B = sqrt((pow((UF - particle[count].up), 2) + pow((VF - particle[count].vp), 2) + pow((WF - particle[count].wp), 2)));
    REYND = den * B * particle[count].pdiam / VISCS;
    CD = .44;
    if (REYND < 1000.)
        CD = (1. + .15 * pow(REYND, .687)) / (REYND / 24.);
    CD = CD * DR;
    A = 3. / 4. * VISCS / (particle[count].pdiam * PDEN) * CD * REYND;
    GR = GRAV - GRAV * (den / PDEN);
    UPM1 = particle[count].up;
    VPM1 = particle[count].vp;
    WPM1 = particle[count].wp;
    particle[count].up = UPM1 + A * (UF - UPM1) * step;
    particle[count].vp = VPM1 + (A * (VF - VPM1) + GR) * step;
    particle[count].wp = WPM1 + A * (WF - WPM1) * step;
    particle[count].xp = particle[count].xp + .5 * (particle[count].up + UPM1) * step;
    particle[count].yp = particle[count].yp + .5 * (particle[count].vp + VPM1) * step;
    particle[count].zp = particle[count].zp + .5 * (particle[count].wp + WPM1) * step;
}


//// Sovle the euler scheme for Fluid Particle
void odeFluidP(double UF, double VF, double WF, vector<struct particles> &FluidParticle, int count,
               double den, double VISCS, double GRAV, double PDEN)
{
    double DR = 0.1, UPM1, VPM1, WPM1;
    double B, REYND, CD, A, GR;
    B = sqrt((pow((UF - FluidParticle[count].up), 2) + pow((VF - FluidParticle[count].vp), 2) + pow((WF - FluidParticle[count].wp), 2)));
    if (B == 0)
        B = 0.0001;
    REYND = den * B * FluidParticle[count].pdiam / VISCS;
    CD = .44;
    if (REYND < 1000.)
        CD = (1. + .15 * pow(REYND, .687)) / (REYND / 24.);
    CD = CD * DR;
    A = 3. / 4. * VISCS / (FluidParticle[count].pdiam * PDEN) * CD * REYND;
    GR = GRAV - (den / PDEN);
    UPM1 = FluidParticle[count].up;
    VPM1 = FluidParticle[count].vp;
    WPM1 = FluidParticle[count].wp;
    FluidParticle[count].up = UPM1 + A * (UF - UPM1) * step;
    FluidParticle[count].vp = VPM1 + (A * (VF - VPM1) + GR) * step;
    FluidParticle[count].wp = WPM1 + A * (WF - WPM1) * step;
    FluidParticle[count].xp = FluidParticle[count].xp + .5 * (FluidParticle[count].up + UPM1) * step;
    FluidParticle[count].yp = FluidParticle[count].yp + .5 * (FluidParticle[count].vp + VPM1) * step;
    FluidParticle[count].zp = FluidParticle[count].zp + .5 * (FluidParticle[count].wp + WPM1) * step;
    if (FluidParticle[count].xp != FluidParticle[count].xp)
        cout << " " << endl; //////   Check for NaN error (NaN values have the property that comparisons involving them are always false)
}

//-----CARRIER PHASE PROPERTIES AT PARTICLE POSITION-------------------
//-----INTERPOLATION COEFFICIENTS
void currierInterp(vector<struct Gnodes> &myNodes, vector<vector<vector<int>>> &ii, vector<std::vector<std::vector<double>>> &Uflow, vector<std::vector<std::vector<double>>> &Vflow,
                   vector<std::vector<std::vector<double>>> &Wflow, double Xp1, double Yp1, double Zp1, int &LX, int &LY, int &LZ, double &UFM, double &VFM, double &WFM)
{
    pC = ii[LX][LY][LZ];
    pN = ii[LX][LY + 1][LZ];
    pS = ii[LX][LY - 1][LZ];
    pE = ii[LX][LY][LZ + 1];
    pW = ii[LX][LY][LZ - 1];
    pD = ii[LX - 1][LY][LZ];
    pU = ii[LX + 1][LY][LZ];

    QN = sqrt(pow(Xp1 - myNodes[pN].x, 2.) + pow(Yp1 - myNodes[pN].y, 2.) + pow(Zp1 - myNodes[pN].z, 2.));
    QS = sqrt(pow(Xp1 - myNodes[pS].x, 2.) + pow(Yp1 - myNodes[pS].y, 2.) + pow(Zp1 - myNodes[pS].z, 2.));
    QC = sqrt(pow(Xp1 - myNodes[pC].x, 2.) + pow(Yp1 - myNodes[pC].y, 2.) + pow(Zp1 - myNodes[pC].z, 2.));
    QE = sqrt(pow(Xp1 - myNodes[pE].x, 2.) + pow(Yp1 - myNodes[pE].y, 2.) + pow(Zp1 - myNodes[pE].z, 2.));
    QW = sqrt(pow(Xp1 - myNodes[pW].x, 2.) + pow(Yp1 - myNodes[pW].y, 2.) + pow(Zp1 - myNodes[pW].z, 2.));
    QD = sqrt(pow(Xp1 - myNodes[pD].x, 2.) + pow(Yp1 - myNodes[pD].y, 2.) + pow(Zp1 - myNodes[pD].z, 2.));
    QU = sqrt(pow(Xp1 - myNodes[pU].x, 2.) + pow(Yp1 - myNodes[pU].y, 2.) + pow(Zp1 - myNodes[pU].z, 2.));
    // cout << myNodes[pU].x << " " << myNodes[pU].y << " " << myNodes[pU].z << " "
    Q1 = QC;
    LX1 = LX;
    LY1 = LY;
    LZ1 = LZ;
    LX2 = LX;
    LZ2 = LZ;
    if (QN < QS)
    {
        Q2 = QN;
        LY2 = LY + 1;
    }
    else
    {
        Q2 = QS;
        LY2 = LY - 1;
    }
    //////////////////////
    LY3 = LY;
    LX3 = LX;
    if (QE < QW)
    {
        Q3 = QE;
        LZ3 = LZ + 1;
    }
    else
    {
        Q3 = QW;
        LZ3 = LZ - 1;
    }
    ////////////////////
    LZ4 = LZ;
    LY4 = LY;
    if (QU < QD)
    {
        Q4 = QU;
        LX4 = LX + 1;
    }
    else
    {
        Q4 = QU;
        LX4 = LX - 1;
    }
    /////////////////////
    LX5 = LX3;
    LY5 = LY2;
    LZ5 = LZ4;
    p5 = ii[LX3][LY2][LZ4];
    Q5 = sqrt(pow(Xp1 - myNodes[p5].x, 2.) + pow(Yp1 - myNodes[p5].y, 2.) + pow(Zp1 - myNodes[p5].z, 2.));
    Q = Q1 + Q2 + Q3 + Q4 + Q5;
    if (Q1 == 0)
        Q1 = 0.00001;
    if (Q2 == 0)
        Q2 = 0.00001;
    if (Q3 == 0)
        Q3 = 0.00001;
    if (Q4 == 0)
        Q4 = 0.00001;
    if (Q5 == 0)
        Q5 = 0.00001;

    POL1 = (Q - Q1) * Q2 * Q3 * Q4 * Q5;
    POL2 = (Q - Q2) * Q1 * Q3 * Q4 * Q5;
    POL3 = (Q - Q3) * Q1 * Q2 * Q4 * Q5;
    POL4 = (Q - Q4) * Q1 * Q2 * Q3 * Q5;
    POL5 = (Q - Q5) * Q1 * Q2 * Q3 * Q4;
    POL = POL1 + POL2 + POL3 + POL4 + POL5;

    UFM = (POL1 * Uflow[LX1][LY1][LZ1] + POL2 * Uflow[LX2][LY2][LZ2] + POL3 * Uflow[LX3][LY3][LZ3] +
           POL4 * Uflow[LX4][LY4][LZ4] + POL5 * Uflow[LX5][LY5][LZ5]) /
          POL;
    VFM = (POL1 * Vflow[LX1][LY1][LZ1] + POL2 * Vflow[LX2][LY2][LZ2] + POL3 * Vflow[LX3][LY3][LZ3] +
           POL4 * Vflow[LX4][LY4][LZ4] + POL5 * Vflow[LX5][LY5][LZ5]) /
          POL;
    WFM = (POL1 * Wflow[LX1][LY1][LZ1] + POL2 * Wflow[LX2][LY2][LZ2] + POL3 * Wflow[LX3][LY3][LZ3] +
           POL4 * Wflow[LX4][LY4][LZ4] + POL5 * Wflow[LX5][LY5][LZ5]) /
          POL;

    if (UFM != UFM)
        cout << " " << endl; //////   Check for NaN error (NaN values have the property that comparisons involving them are always false)
}


/////////   Check if the Fluid Paticle is out of the box
void FPlost(int id, vector<struct particles> &FluidParticle, vector<struct particles> &particle, vector<struct Gnodes> &myNodes,
            vector<std::vector<std::vector<double>>> &Uflow, vector<std::vector<std::vector<double>>> &Vflow, vector<std::vector<std::vector<double>>> &Wflow, int NI, int NJ, int NK)
{
    if (FluidParticle[id].lx >= NI - 2 || FluidParticle[id].ly >= NJ - 2 || FluidParticle[id].lz >= NK - 2 || FluidParticle[id].lx <= 1 || FluidParticle[id].ly < 2 || FluidParticle[id].lz <= 1) // if f.particle is out of box, change values ----  ////////////////////////////
    {
        int dum = ii[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].lx = particle[id].lx;
        FluidParticle[id].ly = particle[id].ly;
        FluidParticle[id].lz = particle[id].lz;
        FluidParticle[id].xp = particle[id].xp;
        FluidParticle[id].yp = particle[id].yp;
        FluidParticle[id].zp = particle[id].zp;
        FluidParticle[id].up = Uflow[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].vp = Vflow[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].wp = Wflow[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].idNum = particle[id].idNum;
    }
    if (FluidParticle[id].up == 0 || FluidParticle[id].vp == 0 || FluidParticle[id].wp == 0) // if f.particle is out of box, change values ----  ////////////////////////////
    {
        dum = ii[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].lx = particle[id].lx;
        FluidParticle[id].ly = particle[id].ly;
        FluidParticle[id].lz = particle[id].lz;
        FluidParticle[id].xp = particle[id].xp;
        FluidParticle[id].yp = particle[id].yp;
        FluidParticle[id].zp = particle[id].zp;
        FluidParticle[id].up = Uflow[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].vp = Vflow[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].wp = Wflow[particle[id].lx][particle[id].ly][particle[id].lz];
        FluidParticle[id].idNum = particle[id].idNum;
        cout << Uflow[particle[id].lx][particle[id].ly][particle[id].lz] << " " << Vflow[particle[id].lx][particle[id].ly][particle[id].lz] << " " << FluidParticle[id].up << " " << FluidParticle[id].vp << endl;
    }
}
bool Checkcell(vector<struct Gnodes> &myNodes, vector<struct particles> &particle, int cell)
{
    int dum = cell;
    double Xmax = myNodes[dum].x + myNodes[dum].sew / 2, Xmin = myNodes[dum].x - myNodes[dum].sew / 2,
           Ymax = myNodes[dum].y + myNodes[dum].sns / 2, Ymin = myNodes[dum].y - myNodes[dum].sns / 2,
           Zmax = myNodes[dum].z + myNodes[dum].sdu / 2, Zmin = myNodes[dum].z - myNodes[dum].sdu / 2;
    double xp = particle[id].xp, yp = particle[id].yp, zp = particle[id].zp;
    if (xp < Xmax && xp > Xmin && yp < Ymax && yp > Ymin && zp < Zmax && zp > Zmin)
        return true;
    else
        return false;
}


//// Check if the particle hit the walls and change the values of the particle
void iwallCheck(int id, vector<struct Gnodes> &myNodes, vector<struct particles> &particle, string IDIR, double &xTline, double &yTline, double &zTline, vector<int> &LXpast, vector<int> &LYpast, vector<int> &LZpast, vector<int> &CELLpast)
{
    int dum = ii[particle[id].lx][particle[id].ly][particle[id].lz];
    int iwall = myNodes[dum].iwall;
    if (iwall == 1 || iwall == 5) 
    {
        particle[id].xp = xTline;
        particle[id].yp = yTline + 1.E-6;
        particle[id].zp = zTline;
        if (IDIR == "E" || IDIR == "W")
            particle[id].up = -particle[id].up;
        if (IDIR == "N" || IDIR == "S")
            particle[id].vp = -particle[id].vp;
        if (IDIR == "U" || IDIR == "D")
            particle[id].wp = -particle[id].wp;
        particle[id].lx = LXpast[id];
        particle[id].ly = LYpast[id];
        particle[id].lz = LZpast[id];
        particle[id].cell = CELLpast[id];
    }
}
//// Find IDIR variable
void ParticleSP(int id, vector<struct Gnodes> &myNodes, vector<vector<vector<int>>> &ii, vector<struct particles> &particle, string &IDIR, vector<Nodes> &iiNodes)
{
    int dum = ii[particle[id].lx][particle[id].ly][particle[id].lz];
    if (Checkcell(myNodes, particle, dum))
        IDIR = "O";
    if (Checkcell(myNodes, particle, iiNodes[dum].iiw))
        particle[id].lx--, IDIR = "W", particle[id].cell = iiNodes[dum].iiw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iie))
        particle[id].lx++, IDIR = "E", particle[id].cell = iiNodes[dum].iie;
    if (Checkcell(myNodes, particle, iiNodes[dum].iin))
        particle[id].ly++, IDIR = "N", particle[id].cell = iiNodes[dum].iin;
    if (Checkcell(myNodes, particle, iiNodes[dum].iis))
        particle[id].ly--, IDIR = "S", particle[id].cell = iiNodes[dum].iis;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiu))
        particle[id].lz++, IDIR = "U", particle[id].cell = iiNodes[dum].iiu;
    if (Checkcell(myNodes, particle, iiNodes[dum].iid))
        particle[id].lz--, IDIR = "D", particle[id].cell = iiNodes[dum].iid;
    if (Checkcell(myNodes, particle, iiNodes[dum].iidw))
        particle[id].lx--, particle[id].lz--, IDIR = "DW", particle[id].cell = iiNodes[dum].iidw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iidn))
        particle[id].ly++, particle[id].lz--, IDIR = "DN", particle[id].cell = iiNodes[dum].iidn;
    if (Checkcell(myNodes, particle, iiNodes[dum].iius))
        particle[id].ly--, particle[id].lz++, IDIR = "US", particle[id].cell = iiNodes[dum].iius;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiun))
        particle[id].ly++, particle[id].lz++, IDIR = "UN", particle[id].cell = iiNodes[dum].iiun;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiuw))
        particle[id].lx--, particle[id].lz++, IDIR = "UW", particle[id].cell = iiNodes[dum].iiuw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiue))
        particle[id].lx++, particle[id].lz++, IDIR = "UE", particle[id].cell = iiNodes[dum].iiue;
    if (Checkcell(myNodes, particle, iiNodes[dum].iisw))
        particle[id].lx--, particle[id].ly--, IDIR = "SW", particle[id].cell = iiNodes[dum].iisw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iinw))
        particle[id].lx--, particle[id].ly++, IDIR = "NW", particle[id].cell = iiNodes[dum].iinw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iise))
        particle[id].lx++, particle[id].ly--, IDIR = "SE", particle[id].cell = iiNodes[dum].iise;
    if (Checkcell(myNodes, particle, iiNodes[dum].iine))
        particle[id].lx++, particle[id].ly++, IDIR = "NE", particle[id].cell = iiNodes[dum].iine;
    if (Checkcell(myNodes, particle, iiNodes[dum].iidse))
        particle[id].lx++, particle[id].ly--, particle[id].lz--, IDIR = "DSE", particle[id].cell = iiNodes[dum].iidse;
    if (Checkcell(myNodes, particle, iiNodes[dum].iidsw))
        particle[id].lx--, particle[id].ly--, particle[id].lz--, IDIR = "DSW", particle[id].cell = iiNodes[dum].iidsw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iidnw))
        particle[id].lx--, particle[id].ly++, particle[id].lz--, IDIR = "DNW", particle[id].cell = iiNodes[dum].iidnw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iidne))
        particle[id].lx++, particle[id].ly++, particle[id].lz--, IDIR = "DNE", particle[id].cell = iiNodes[dum].iidne;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiusw))
        particle[id].lx--, particle[id].ly--, particle[id].lz++, IDIR = "USW", particle[id].cell = iiNodes[dum].iiusw;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiuse))
        particle[id].lx++, particle[id].ly--, particle[id].lz++, IDIR = "USE", particle[id].cell = iiNodes[dum].iiuse;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiuwn))
        particle[id].lx--, particle[id].ly++, particle[id].lz++, IDIR = "UWN", particle[id].cell = iiNodes[dum].iiuwn;
    if (Checkcell(myNodes, particle, iiNodes[dum].iiune))
        particle[id].lx++, particle[id].ly++, particle[id].lz++, IDIR = "UNE", particle[id].cell = iiNodes[dum].iiune;
}

// Find crosspoint of particle motion and wall
void tline(int id, vector<struct Gnodes> &myNodes, vector<vector<vector<int>>> &ii, vector<struct particles> &particle,
           double Xm1, double Ym1, double Zm1, string &IDIR, double &xb, double &yb, double &zb)

{
    // xb, yb, zb; /// Οι θέσεις τις τομής ευθείας (σωματιδιου) - επιφάνεια χωρίου
    int LX = particle[id].lx;
    int LY = particle[id].ly;
    int LZ = particle[id].lz;
    int dumx = ii[LX][LY][LZ];
    double Xmax = myNodes[dumx].x + myNodes[dumx].sdu / 2, Xmin = myNodes[dumx].x - myNodes[dumx].sdu / 2,
           Ymax = myNodes[dumx].y + myNodes[dumx].sns / 2, Ymin = myNodes[dumx].y - myNodes[dumx].sns / 2,
           Zmax = myNodes[dumx].z + myNodes[dumx].sew / 2, Zmin = myNodes[dumx].z - myNodes[dumx].sew / 2;
    double xp = particle[id].xp, yp = particle[id].yp, zp = particle[id].zp;

    double a = xp - Xm1, b = yp - Ym1, c = zp - Zm1;
    double t, XN, YN, ZN, dN, XS, YS, ZS, dS, XU, YU, ZU, dU, XD, YD, ZD, dD, XE, YE, ZE, dE, XW, YW, ZW, dW;
    /// Distance from North plane
    XN = Xm1 + a / b * (Ymin - Ym1);
    YN = Ym1 + (Ymin - Ym1);
    ZN = Zm1 + c / b * (Ymin - Ym1);
    dN = sqrt((XN - Xm1) * (XN - Xm1) + (YN - Ym1) * (YN - Ym1) + (ZN - Zm1) * (ZN - Zm1));
    /////////////////////////////
    /// Distance from South plane
    XS = Xm1 + a / b * (Ymax - Ym1);
    YS = Ym1 + (Ymax - Ym1);
    ZS = Zm1 + c / b * (Ymax - Ym1);
    dS = sqrt((XS - Xm1) * (XS - Xm1) + (YS - Ym1) * (YS - Ym1) + (ZS - Zm1) * (ZS - Zm1));
    /////////////////////////////
    /// Distance from Upper plane
    XU = Xm1 + (Xmin - Xm1);
    YU = Ym1 + b / a * (Xmin - Xm1);
    ZU = Zm1 + c / a * (Xmin - Xm1);
    dU = sqrt((XU - Xm1) * (XU - Xm1) + (YU - Ym1) * (YU - Ym1) + (ZU - Zm1) * (ZU - Zm1));
    /////////////////////////////
    /// Distance from Down plane
    XD = Xm1 + (Xmax - Xm1);
    YD = Ym1 + b / a * (Xmax - Xm1);
    ZD = Zm1 + c / a * (Xmax - Xm1);
    dD = sqrt((XD - Xm1) * (XD - Xm1) + (YD - Ym1) * (YD - Ym1) + (ZD - Zm1) * (ZD - Zm1));
    /////////////////////////////
    /// Distance from East plane
    XE = Xm1 + a / c * (Zmin - Zm1);
    YE = Ym1 + b / c * (Zmin - Zm1);
    ZE = Zm1 + (Zmin - Zm1);
    dE = sqrt((XE - Xm1) * (XE - Xm1) + (YE - Ym1) * (YE - Ym1) + (ZE - Zm1) * (ZE - Zm1));
    /////////////////////////////
    /// Distance from West plane
    XW = Xm1 + a / c * (Zmax - Zm1);
    YW = Ym1 + b / c * (Zmax - Zm1);
    ZW = Zm1 + (Zmax - Zm1);
    dW = sqrt((XW - Xm1) * (XW - Xm1) + (YW - Ym1) * (YW - Ym1) + (ZW - Zm1) * (ZW - Zm1));
    /////////////////////////////
    if (IDIR == "U")
    {
        xb = XU;
        yb = YU;
        zb = ZU;
    }
    if (IDIR == "D")
    {
        xb = XD;
        yb = YD;
        zb = ZD;
    }
    if (IDIR == "S")
    {
        xb = XS;
        yb = YS;
        zb = ZS;
    }
    if (IDIR == "N")
    {
        xb = XN;
        yb = YN;
        zb = ZN;
    }
    if (IDIR == "E")
    {
        xb = XE;
        yb = YE;
        zb = ZE;
    }
    if (IDIR == "W")
    {
        xb = XW;
        yb = YW;
        zb = ZW;
    }
    if (IDIR == "NE")
    {
        if (dN < dE)
        {
            xb = XN;
            yb = YN;
            zb = ZN;
            IDIR = "N";
        }
        else
        {
            xb = XE;
            yb = YE;
            zb = ZE;
            IDIR = "E";
        }
    }
    if (IDIR == "SE")
    {
        if (dS < dE)
        {
            xb = XS;
            yb = YS;
            zb = ZS;
            IDIR = "S";
        }
        else
        {
            xb = XE;
            yb = YE;
            zb = ZE;
            IDIR = "E";
        }
    }
    if (IDIR == "NW")
    {
        if (dN < dW)
        {
            xb = XN;
            yb = YN;
            zb = ZN;
            IDIR = "N";
        }
        else
        {
            xb = XW;
            yb = YW;
            zb = ZW;
            IDIR = "W";
        }
    }
    if (IDIR == "SW")
    {
        if (dS < dW)
        {
            xb = XS;
            yb = YS;
            zb = ZN;
            IDIR = "S";
        }
        else
        {
            xb = XW;
            yb = YW;
            zb = ZW;
            IDIR = "W";
        }
    }

    if (IDIR == "UN")
    {
        if (dU < dN)
        {
            xb = XU;
            yb = YU;
            zb = ZU;
            IDIR = "U";
        }
        else
        {
            xb = XN;
            yb = YN;
            zb = ZN;
            IDIR = "N";
        }
    }
    if (IDIR == "DN")
    {
        if (dD < dN)
        {
            xb = XD;
            yb = YD;
            zb = ZD;
            IDIR = "D";
        }
        else
        {
            xb = XN;
            yb = YN;
            zb = ZN;
            IDIR = "N";
        }
    }
    if (IDIR == "US")
    {
        if (dU < dS)
        {
            xb = XU;
            yb = YU;
            zb = ZU;
            IDIR = "U";
        }
        else
        {
            xb = XS;
            yb = YS;
            zb = ZS;
            IDIR = "S";
        }
    }
    if (IDIR == "DS")
    {
        if (dD < dS)
        {
            xb = XD;
            yb = YD;
            zb = ZD;
            IDIR = "D";
        }
        else
        {
            xb = XS;
            yb = YS;
            zb = ZS;
            IDIR = "S";
        }
    }
    if (IDIR == "UE")
    {
        if (dU < dE)
        {
            xb = XU;
            yb = YU;
            zb = ZU;
            IDIR = "U";
        }
        else
        {
            xb = XE;
            yb = YE;
            zb = ZE;
            IDIR = "E";
        }
    }
    if (IDIR == "DE")
    {
        if (dD < dE)
        {
            xb = XD;
            yb = YD;
            zb = ZD;
            IDIR = "D";
        }
        else
        {
            xb = XE;
            yb = YE;
            zb = ZE;
            IDIR = "E";
        }
    }
    if (IDIR == "DW")
    {
        if (dD < dW)
        {
            xb = XD;
            yb = YD;
            zb = ZD;
            IDIR = "D";
        }
        else
        {
            xb = XW;
            yb = YW;
            zb = ZW;
            IDIR = "W";
        }
    }
    if (IDIR == "UW")
    {
        if (dU < dW)
        {
            xb = XU;
            yb = YU;
            zb = ZU;
            IDIR = "U";
        }
        else
        {
            xb = XW;
            yb = YW;
            zb = ZW;
            IDIR = "W";
        }
    }
}

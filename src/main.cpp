#include <iostream>
#include <cstring>
#include <math.h>
#include <vector>
#include "particles.h"
#include "iodata.h"
#include "Grid.h"
#include "GNodes.h"
#include "functions.h"
#include "RW_BinaryFormats.h"
#include "GAUSS.h"
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
using namespace std;

int main()
{
    Traj1.open("Traj1.dat"), CONC_.open("CONC.dat"), Stk_num.open("Stokes_num.dat"), CCC.open("CCC");
    while (CCC >> x)
    {
        PINAX.push_back(x);
    }
    srand(time(NULL));
    /* generate secret number between 0 and 1: */
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    /// initialize vectors
    LXpast.resize(NPART), LYpast.resize(NPART), LZpast.resize(NPART), XPpast.resize(NPART), YPpast.resize(NPART),
        ZPpast.resize(NPART), UPpast.resize(NPART), VPpast.resize(NPART), WPpast.resize(NPART), IDpast.resize(NPART), CELLpast.resize(NPART);
    LXFpast.resize(NPART), LYFpast.resize(NPART), LZFpast.resize(NPART), XPFpast.resize(NPART), YPFpast.resize(NPART),
        ZPFpast.resize(NPART), UPFpast.resize(NPART), VPFpast.resize(NPART), WPFpast.resize(NPART), IDFpast.resize(NPART), CELLFpast.resize(NPART);

    //////////////////////////////////////////////////////////////////////////////////
    // Grid initialization  //////////////////////////////////////////////////////////
    StructuredMesh GridInfo;
    GridInfo.readGrid("inputData/3DGRID_00");
    GridInfo.DefineMesh();
    GridInfo.initializeVariables();
    GridInfo.StringIndexing();
    vector<Gnodes> myNodes;
    SDU = GridInfo.getSDU();
    SNS = GridInfo.getSNS();
    SEW = GridInfo.getSEW();

    vector<Nodes> iiNodes;
    iiNodes = GridInfo.getiiNodes();
    int NI = GridInfo.getNI();
    int NJ = GridInfo.getNJ();
    int NK = GridInfo.getNK();

    iwallVec = GridInfo.getiwall();
    ii = GridInfo.getii();
    Xgrid = GridInfo.getXgrid();
    Ygrid = GridInfo.getYgrid();
    Zgrid = GridInfo.getZgrid();
    iwallV1 = GridInfo.getIwall1D(); //from ASCII filw
    int NIJKM1 = (NI - 1) * (NJ - 1) * (NK - 1);
    int NIJKM2 = (NI - 2) * (NJ - 2) * (NK - 2);
    // Vectors Grid_values;
    // initialize_doubles(Grid_values.A, NIJKM1, Xgrid);
    // initialize_doubles(Grid_values.B, NIJKM1, Ygrid);
    // initialize_doubles(Grid_values.C, NIJKM1, Zgrid);

    // // write to file
    // Grid_values.saveGrid("dataGrid.bin");

    // read back in to memory
    Vectors loaded_Grid_values;
    loaded_Grid_values.loadGrid("dataGrid.bin");
    Xload = loaded_Grid_values.getXgrid();
    Yload = loaded_Grid_values.getYgrid();
    Zload = loaded_Grid_values.getZgrid();
    Iwallload = loaded_Grid_values.getIwall1d();
    //////////////////////////////////////////////////////////////////////////////////
    ////// Create stuct "myNodes" ////////////////////////////////////////////////////
    Gnodes Gnode;
    int counter = 0;
    for (int k = 1; k <= NK - 1; k++)
    {
        Gnode.k = k;
        for (int j = 1; j <= NJ - 1; j++)
        {
            Gnode.j = j;
            for (int i = 1; i <= NI - 1; i++)
            {
                Gnode.x = Xload[counter];
                Gnode.y = Yload[counter];
                Gnode.z = Zload[counter];
                Gnode.iwall = iwallV1[counter];
                Gnode.sdu = SDU[k];
                Gnode.sew = SEW[i];
                Gnode.sns = SNS[j];
                Gnode.i = i;
                myNodes.push_back(Gnode);
                counter++;
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    // End - Grid initialization  ////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    //Fluid Flow initialization //////////////////////////////////////////////////////
    iodata ProInfo(NI, NJ, NK);
    ProInfo.resizeVec();
    ProInfo.readFluidPro("inputData/FLUIDINPUT"); //
    ProInfo.readPartPro("inputData/PRTCLINPUT");  //

    double VISCS = ProInfo.getVisc();
    double den = ProInfo.getDen();

    pdiami = ProInfo.getPdiami();
    double PDEN = ProInfo.getPDEN0();
    double GRAV = ProInfo.getGRAV();
    int NPDIA = ProInfo.getNPDIA();

    int MAXTSPENTRA = ProInfo.getMAXTSPENTRA();
    int JEMIS = ProInfo.getJEMIS();
    int IEMIS = ProInfo.getIEMIS();
    int ZEMIS1 = ProInfo.getZEMIS1();
    int ZEMIS2 = ProInfo.getZEMIS2();
    int NTIMST = 1;
    double UPFS = ProInfo.getTURINTU();
    double VPFS = ProInfo.getTURINTV();
    double WPFS = ProInfo.getTURINTW();
    int MAXTSPSV = ProInfo.getMAXTSPSV();
    double DisCoef = ProInfo.getDisCoef();
    int LLL = ProInfo.getLLL();

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////// Use these valaus for saving ASCII data in BINARY format /////////////////////////////////
    // Uflow1dd = ProInfo.getU1d();
    // Vflow1dd = ProInfo.getV1d();
    // Wflow1dd = ProInfo.getW1d();
    // TEflow1dd = ProInfo.getTE1d();
    // EDflow1dd = ProInfo.getED1d();
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Vectors Flow_values;
    // initialize_doubles(Flow_values.U, NIJKM2, Uflow1dd);
    // initialize_doubles(Flow_values.V, NIJKM2, Vflow1dd);
    // initialize_doubles(Flow_values.W, NIJKM2, Wflow1dd);
    // initialize_doubles(Flow_values.T, NIJKM2, TEflow1dd);
    // initialize_doubles(Flow_values.E, NIJKM2, EDflow1dd);
    // // // // write to file
    // Flow_values.saveFlow("dataFlow.bin");

    // // read back in to memory
    Vectors loaded_Flow_values, loaded_Flow_KE;
    loaded_Flow_values.loadFlow("dataFlow.bin");
    // loaded_Flow_KE.loadFlow("dataFlow.bin");

    UflowLoad = loaded_Flow_values.getUflow();
    VflowLoad = loaded_Flow_values.getVflow();
    WflowLoad = loaded_Flow_values.getWflow();
    TEflowLoad = loaded_Flow_values.getTEflow();
    EDflowLoad = loaded_Flow_values.getEDflow();
    Uflow.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    Vflow.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    Wflow.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    TE.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    // ED.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    Concentration.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    int dd1 = 0;
    for (int i = 1; i < NI - 1; i++)
    {
        for (int j = 1; j < NJ - 1; j++)
        {
            for (int k = 1; k < NK - 1; k++)
            {
                Uflow[i][j][k] = UflowLoad[dd1];
                Vflow[i][j][k] = VflowLoad[dd1];
                Wflow[i][j][k] = WflowLoad[dd1];
                TE[i][j][k] = TEflowLoad[dd1];
                // ED[i][j][k] = EDflowLoad[dd1];
                dd1++;
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // End Fluid Flow initialization ///////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    /////  Stokes Num   ////////////////////////////////////////////
    double VM = sqrt(Uflow[5][5][5] * Uflow[5][5][5] + Vflow[5][5][5] * Vflow[5][5][5] + Wflow[5][5][5] * Wflow[5][5][5]);
    double fb_d = 2.5;
    int duum = 0;
    double FluidParDiameter = 1.E-6;
    double FluidParDensity = 1.225;
    double tp = PDEN * (pdiami[0] * pdiami[0]) / (18 * VISCS); // response time of the particle
    double tmf = VM / fb_d;                                    //time scale of the fluid based on main flow
    double stk_mf = tp / tmf;
    Stk_num << "Stokes Number = " << stk_mf << endl;
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    double dmdtEXP = 0.0066783;                      // Mass Flow Rate of MUST
    double ParNum = MMM * LLL * MAXTSPENTRA * NPDIA; // Maximum num of particles / trajectory
    double dmdtP = dmdtEXP / (MMM * ParNum);         // Mass Flow Rate of each particle / trajectory
    bool iloop = true;
    int dum;
    int imaxtsp;
    while (iloop)
    {
        id = 0;
        imaxtsp = min(MAXTSPENTRA, NTIMST);
        for (int istep = 0; istep < imaxtsp; istep++)
        {
            for (int m = 0; m < MMM; m++)
            {
                for (int nd = 0; nd < NPDIA; nd++)
                {
                    for (int l = 0; l < LLL; l++)
                    {
                        if (istep == NTIMST - 1)
                        {
                            LX = IEMIS;
                            LY = JEMIS;
                            LZ = ZEMIS1;
                            LXf = IEMIS;
                            LYf = JEMIS;
                            LZf = ZEMIS1;
                            cell = ii[LX][LY][LZ];
                            dum = ii[LX][LY][LZ];
                            idNum = std::to_string(id);
                            DSEED = distribution(generator);
                            GGG = GAUSS(DSEED, PINAX);
                            Kturb = TE[LX][LY][LZ];
                            // Eturb = ED[LX][LY][LZ];
                            Eturb = EDflowLoad[dum];
                            UPFS = sqrt(2 * Kturb / 3);
                            UP1 = Uflow[LX][LY][LZ] + DisCoef * UPFS * GGG; //+ UPFS * GGG;
                            VP1 = Vflow[LX][LY][LZ] + DisCoef * UPFS * GGG;
                            WP1 = Wflow[LX][LY][LZ] + DisCoef * UPFS * GGG;
                            particle.push_back(particles(LX, LY, LZ, myNodes[dum].x, myNodes[dum].y, myNodes[dum].z,
                                                         UP1, VP1, WP1, pdiami[nd], idNum, cell));
                            particle[id].setParDensity(PDEN);
                            particle[id].setParMass();
                            FluidParticle.push_back(particles(LX, LY, LZ, myNodes[dum].x, myNodes[dum].y, myNodes[dum].z,
                                                              UP1, VP1, WP1, FluidParDiameter, idNum, cell));
                            particle[id].sumtime = 1;
                        }
                        if (istep < NTIMST - 1)
                        {
                            cell = ii[LX][LY][LZ];
                            LX = LXpast[id];
                            LY = LYpast[id];
                            LZ = LZpast[id];
                            particle[id].cell = CELLpast[id];
                            particle[id].lx = LX;
                            particle[id].ly = LY;
                            particle[id].lz = LZ;
                            particle[id].xp = XPpast[id];
                            particle[id].yp = YPpast[id];
                            particle[id].zp = ZPpast[id];
                            particle[id].up = UPpast[id];
                            particle[id].vp = VPpast[id];
                            particle[id].wp = WPpast[id];
                            particle[id].idNum = IDpast[id];
                            //////////////////////////////
                            //////////////////////////////
                            LXf = LXFpast[id];
                            LYf = LYFpast[id];
                            LZf = LZFpast[id];
                            FluidParticle[id].cell = CELLFpast[id];
                            FluidParticle[id].lx = LXf;
                            FluidParticle[id].ly = LYf;
                            FluidParticle[id].lz = LZf;
                            FluidParticle[id].xp = XPFpast[id];
                            FluidParticle[id].yp = YPFpast[id];
                            FluidParticle[id].zp = ZPFpast[id];
                            FluidParticle[id].up = UPFpast[id];
                            FluidParticle[id].vp = VPFpast[id];
                            FluidParticle[id].wp = WPFpast[id];
                            FluidParticle[id].idNum = IDFpast[id];
                        }
                        str = to_string(id) + ".Out_of_box";
                        if (particle[id].idNum == str)
                        {
                            goto particle_out;
                        }
                        str2 = to_string(id) + ".Stuck_in_Filed";
                        if (particle[id].idNum == str2)
                        {
                            goto particle_out;
                        }

                        currierInterp(myNodes, ii, Uflow, Vflow, Wflow, XPpast[id], YPpast[id], ZPpast[id], LX, LY, LZ, UFM, VFM, WFM);          // Lagrange Interp
                        currierInterp(myNodes, ii, Uflow, Vflow, Wflow, XPFpast[id], YPFpast[id], ZPFpast[id], LXf, LYf, LZf, UFPM, VFPM, WFPM); // Lagrange Interp for Fluid Particle
                        ////----------------------------------------------------------------------
                        ode(UFM, VFM, WFM, particle, id, den, VISCS, GRAV, PDEN);                          // Euler Scheme
                        odeFluidP(UFPM, VFPM, WFPM, FluidParticle, id, den, VISCS, GRAV, FluidParDensity); //Euler Scheme for the fluid particle motion
                        ////----------------------------------------------------------------------
                        ParticleSP(id, myNodes, ii, particle, IDIR, iiNodes);              // Find IDIR variable
                        ParticleSP(id, myNodes, ii, FluidParticle, IDIRfluidPar, iiNodes); //Find IDIR variable for Fluid Par
                        ////----------------------------------------------------------------------
                        if (IDIR != "O")
                            tline(id, myNodes, ii, particle, XPpast[id], YPpast[id], ZPpast[id], IDIR, xTline, yTline, zTline); // crosspoints xTline, yTline, zTline from trajector line and solid wall.
                        dum = ii[particle[id].lx][particle[id].ly][particle[id].lz];
                        ////----------------------------------------------------------------------
                        iwallCheck(id, myNodes, particle, IDIR, xTline, yTline, zTline, LXpast, LYpast, LZpast, CELLpast);
                        Le = pow(Cm, 0.5) * pow(Kturb, 1.5) / Eturb; // constant eddy length scale
                        dFP = sqrt((particle[id].xp - FluidParticle[id].xp) * (particle[id].xp - FluidParticle[id].xp) + (particle[id].yp - FluidParticle[id].yp) * (particle[id].yp - FluidParticle[id].yp) + (particle[id].zp - FluidParticle[id].zp) * (particle[id].zp - FluidParticle[id].zp));
                        if (dFP > Le) // Check if the leghth of particle and Fludi particle is bigger than the eddy length scale
                        {
                            dum = ii[particle[id].lx][particle[id].ly][particle[id].lz];
                            DSEED = distribution(generator);
                            GGG = GAUSS(DSEED, PINAX);
                            Kturb = TE[particle[id].lx][particle[id].ly][particle[id].lz];
                            Eturb = EDflowLoad[dum];
                            UPFS = sqrt(2 * Kturb / 3);
                            particle[id].up = particle[id].up + DisCoef * UPFS * GGG; //;
                            particle[id].vp = particle[id].vp + DisCoef * UPFS * GGG;
                            particle[id].wp = particle[id].wp + DisCoef * UPFS * GGG;
                            FluidParticle[id].lx = particle[id].lx;
                            FluidParticle[id].ly = particle[id].ly;
                            FluidParticle[id].lz = particle[id].lz;
                            FluidParticle[id].xp = particle[id].xp;
                            FluidParticle[id].yp = particle[id].yp;
                            FluidParticle[id].zp = particle[id].zp;
                            FluidParticle[id].up = Uflow[particle[id].lx][particle[id].ly][particle[id].lz];
                            FluidParticle[id].vp = Vflow[particle[id].lx][particle[id].ly][particle[id].lz];
                            FluidParticle[id].wp = Wflow[particle[id].lx][particle[id].ly][particle[id].lz];
                            FluidParticle[id].cell = particle[id].cell;
                        }
                        ///////////////////////////////////////////////////////////////////////////////
                        // Concentration Calc   ////////////////////////////////////////////////////////////////////
                        if (particle[id].cell != CELLpast[id])
                        {
                            VolCell = myNodes[dumCh].sdu * myNodes[dumCh].sns * myNodes[dumCh].sew; /// Volume of Cell
                            Ckg = dmdtP * step / VolCell;
                            iterStep = 1;
                            Concentration[LX][LY][LZ] += Ckg;
                        }
                        else
                        {
                            Ckg = dmdtP * step / VolCell;
                            iterStep++;
                            Concentration[LX][LY][LZ] += Ckg;
                        }
                        if (isinf(Concentration[LX][LY][LZ]))
                        {
                            cout << " Error :Concentration is Infinity Num" << endl;
                        }

                        ///////////////////////////////////////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////
                        /////////// Check if the particle is out of the box Or Stuck in the Field//
                        if (LX >= NI - 3 || LY >= NJ - 3 || LZ >= NK - 3 || LX <= 1 || LY < 2 || LZ <= 1)
                        {
                            particle[id].idNum += ".Out_of_box";
                        }
                        ///////////////////////////////////////////////////////////////////////////
                        FPlost(id, FluidParticle, particle, myNodes, Uflow, Vflow, Wflow, NI, NJ, NK);
                        ///////////////////////////////////////////////////////////////////////////
                        if (particle[id].cell != CELLpast[id])
                            particle[id].sumtime = 1;
                        else
                            particle[id].sumtime++;
                        if (particle[id].sumtime > 200)
                            particle[id].idNum += ".Stuck_in_Filed";
                        ///////////////////////////////////////////////////////////////////////////
                        if (isNotNumber(particle[id].idNum))
                            LostParticles++; // Check the num of lost particles
                        ///////////////////////////////////////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////
                        // write_data:
                        //////////////////// Write_Data for the traces of particles////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////////////////////////////////
                        // if (nd == 0 && iter % 8 == 0)
                        Traj1 << particle[id].xp << " " << particle[id].yp << " " << particle[id].zp << endl;
                        ////////////////////////////////////////////////////////////
                        //--SAVE VALUES OF Kth PARTICLE
                        LXpast[id] = particle[id].lx;
                        LYpast[id] = particle[id].ly;
                        LZpast[id] = particle[id].lz;
                        CELLpast[id] = particle[id].cell;
                        XPpast[id] = particle[id].xp;
                        YPpast[id] = particle[id].yp;
                        ZPpast[id] = particle[id].zp;
                        UPpast[id] = particle[id].up;
                        VPpast[id] = particle[id].vp;
                        WPpast[id] = particle[id].wp;
                        IDpast[id] = particle[id].idNum;
                        //--SAVE VALUES OF Kth FLUID PARTICLE
                        LXFpast[id] = FluidParticle[id].lx;
                        LYFpast[id] = FluidParticle[id].ly;
                        LZFpast[id] = FluidParticle[id].lz;
                        CELLFpast[id] = FluidParticle[id].cell;
                        XPFpast[id] = FluidParticle[id].xp;
                        YPFpast[id] = FluidParticle[id].yp;
                        ZPFpast[id] = FluidParticle[id].zp;
                        UPFpast[id] = FluidParticle[id].up;
                        VPFpast[id] = FluidParticle[id].vp;
                        WPFpast[id] = FluidParticle[id].wp;
                        IDFpast[id] = FluidParticle[id].idNum;
                        ///////////////////////////////////////////////////////////////////////
                        ///////////////////////////////////////////////////////////////////////
                        if ((iter % (10 * MAXTSPSV)) == 0)
                        {
                            cout << "  Number of time Step: " << NTIMST << endl;
                            cout << " Particle : " << particle[id].idNum << "     "
                                 << " Lost Particles "
                                 << " " << LostParticles << "     "
                                 << " X - Y - Z :" << particle[id].xp << " " << particle[id].yp << " " << particle[id].zp << endl;
                            cout << " U - V - W :"
                                 << " "
                                 << " " << particle[id].up << " " << particle[id].vp << " " << particle[id].wp << endl;
                            cout << " " << endl;
                        }

                        //////////////////// write_concentration///////////////////////////////
                        ///////////////////////////////////////////////////////////////////////

                        if ((NTIMST >= 10000) && (iter % (20 * MAXTSPSV)) == 0)
                        {
                            std::ofstream CONC_("CONC.dat", std::ofstream::trunc);
                            // int ddum1 = 0;
                            cout << " Write Concentration : "
                                 << " " << endl;
                            for (int k = 0; k < NK - 1; k++)
                            {
                                for (int j = 0; j < NJ - 1; j++)
                                {
                                    for (int i = 0; i < NI - 1; i++)
                                    {
                                        // CONC_ << CONC1m[i][j][k] << endl;
                                        CONC_ << Concentration[i][j][k] << endl;

                                        // ddum1++;
                                    }
                                }
                            }
                            cout << " End of write_concentration " << endl;
                        }

                        CONC_.close();

                    particle_out:;
                        id++;

                        if (FluidParticle[id].xp != FluidParticle[id].xp)
                            cout << " " << endl; //////   Check for NaN error (NaN values have the property that comparisons involving them are always false)

                        iter++;
                    };
                };
            };
        };
        NTIMST++;
        if (NTIMST > 250000000)
            iloop = false;
    };

    cout << iloop << " " << NTIMST << endl;
    cin.get();
}
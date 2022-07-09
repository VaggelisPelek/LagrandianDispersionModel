#include "Grid.h"

StructuredMesh::StructuredMesh()
{
}

void StructuredMesh::DefineMesh()
{
    int dum = 0;
    double tmp1, tmp2, tmp3;
    double dxe, dxw, dys, dyp, dzu, dzd, sew, sns, sdu;
    int iwal;
    iwall.resize(NI, vector<vector<int>>(NJ, vector<int>(NK)));
    ifstream myfile("inputData/XYZcord.dat", ios::in);
    ifstream myfileSEW("inputData/SEW.dat", ios::in);
    ifstream myfileSNS("inputData/SNS.dat", ios::in);
    ifstream myfileSDU("inputData/SDU.dat", ios::in);
    ifstream myfileIWALL("inputData/IWALL.dat", ios::in);

    double number;

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    // X-cordinate, Y-cordinate, Z-cordinate
    // for (int k = 1; k <= NK - 1; k++)
    // {
    //     for (int j = 1; j <= NJ - 1; j++)
    //     {
    //         for (int i = 1; i <= NI - 1; i++)
    //         {
    //             myfile >> tmp1;
    //             myfile >> tmp2;
    //             myfile >> tmp3;
    //             XGrid.push_back(tmp1);
    //             YGrid.push_back(tmp2);
    //             ZGrid.push_back(tmp3);
    //         }
    //     }
    // }
    // IWALL
    double xx, yy, zz;
    for (int k = 1; k <= NK - 1; k++)
    {
        for (int j = 1; j <= NJ - 1; j++)
        {
            for (int i = 1; i <= NI - 1; i++)
            {

                myfileIWALL >> iwal;
                iwall1d.push_back(iwal);
                dum++;
            }
        }
    }

    // DXEP - DXPW
    for (int i = 1; i < NI; i++)
    {
        myfileSEW >> sew;
        SEW.push_back(sew);
    }
    //// DYPS - DYNP
    for (int j = 0; j < NJ; j++)
    {
        myfileSNS >> sns;
        SNS.push_back(sns);
    }
    //// DZPU - DZDP
    for (int k = 0; k < NK; k++)
    {
        myfileSDU >> sdu;
        SDU.push_back(sdu);
    }
    ////////////////////////////////
    //////////////////////////////////
}

void StructuredMesh ::readGrid(const string fileName) //Read from 3DGRID_00
{
    ifstream Tridgrid;
    int size = 0;

    // Tridgrid.open((name.str()).c_str(), ios::in | ios::binary);
    Tridgrid.open(fileName, std::ios::in | std::ios::binary);
    Tridgrid.seekg(0, ios::end);
    size = (int)Tridgrid.tellg();
    Tridgrid.seekg(0, ios::beg);
    const int Header_Length = sizeof(int);  // UNFORMATTED form of FORTRAN has a header and trailer in every record of the binary file
    const int Trailer_Length = sizeof(int); //
    char header[Header_Length];
    char trailer[Trailer_Length];

    Tridgrid.read(header, sizeof(header)); // read header
    Tridgrid.read((char *)&(NI), sizeof(int));
    Tridgrid.read((char *)&(NJ), sizeof(int));
    Tridgrid.read((char *)&(NK), sizeof(int));
    Tridgrid.read(trailer, sizeof(trailer));
    Tridgrid.read(header, sizeof(header));

    Tridgrid.close();
}

void StructuredMesh::initializeVariables()
{
    ii.resize(NI, vector<vector<int>>(NJ, vector<int>(NK)));
}

void StructuredMesh::StringIndexing()
{
    ofstream testii, iivec;
    testii.open("iitest"), iivec.open("ii");
    int ismax = 0;
    iiNodes;
    Nodes node;
    node.iie = 0, node.iiw = 0, node.iis = 0, node.iin = 0,
    node.iid = 0, node.iiu = 0, node.iids = 0, node.iidn = 0,
    node.iide = 0, node.iidw = 0, node.iius = 0, node.iiun = 0,
    node.iiuw = 0, node.iiue = 0, node.iisw = 0, node.iinw = 0,
    node.iise = 0, node.iine = 0, node.iidsw = 0, node.iidse = 0,
    node.iiuse = 0, node.iiusw = 0, node.iidnw = 0, node.iidne = 0,
    node.iiuwn = 0, node.iiune = 0;
    for (int iz = 1; iz < NK; iz++)
    {
        for (int iy = 1; iy < NJ; iy++)
        {
            for (int ix = 1; ix < NI; ix++)
            {
                ii[ix][iy][iz] = ismax;
                iiNodes.push_back(node);
                ismax++;
            }
        }
    }

    // /////////////  Creat ii[][][]  matrix
    int dum = 0;
    int L = 0;
    for (int iz = 2; iz < NK - 1; iz++)
    {
        for (int iy = 2; iy < NJ - 1; iy++)
        {
            for (int ix = 2; ix < NI - 1; ix++)
            {
                L = ii[ix][iy][iz];
                if (ix < NI)
                {
                    iiNodes[L].iie = ii[ix + 1][iy][iz];
                    if (iy < NJ)
                    {
                        iiNodes[L].iine = ii[ix + 1][iy + 1][iz];
                    }
                    else
                    {
                        iiNodes[L].iine = 0;
                    }
                    if (iy > 1)
                    {
                        iiNodes[L].iise = ii[ix + 1][iy - 1][iz];
                    }
                    else
                    {
                        iiNodes[L].iise = 0;
                    }
                }
                else
                {
                    iiNodes[L].iie = 0;
                }

                if (ix > 1)
                {
                    iiNodes[L].iiw = ii[ix - 1][iy][iz];
                    if (iy < NJ)
                    {
                        iiNodes[L].iinw = ii[ix - 1][iy + 1][iz];
                    }
                    else
                    {
                        iiNodes[L].iinw = 0;
                    }
                    if (iy > 1)
                    {
                        iiNodes[L].iisw = ii[ix - 1][iy - 1][iz];
                    }
                    else
                    {
                        iiNodes[L].iisw = 0;
                    }
                }
                else
                {
                    iiNodes[L].iiw = 0;
                }

                if (iy < NJ)
                {
                    iiNodes[L].iin = ii[ix][iy + 1][iz];
                }
                else
                {
                    iiNodes[L].iin = 0;
                }
                if (iy > 1)
                {
                    iiNodes[L].iis = ii[ix][iy - 1][iz];
                }
                else
                {
                    iiNodes[L].iis = 0;
                }

                ///// U - plane ////////////////////////////////////////////////////////////////////
                if (iz < NK)
                {
                    iiNodes[L].iiu = ii[ix][iy][iz + 1];
                }
                else
                {
                    iiNodes[L].iiu = 0;
                }

                if (iz < NK && ix < NI)
                {
                    iiNodes[L].iiue = ii[ix + 1][iy][iz + 1];
                }
                else
                {
                    iiNodes[L].iiue = 0;
                }

                if (iz < NK && ix > 1)
                {
                    iiNodes[L].iiuw = ii[ix - 1][iy][iz + 1];
                }
                else
                {
                    iiNodes[L].iiuw = 0;
                }

                if (iz < NK && iy < NJ)
                {
                    iiNodes[L].iiun = ii[ix][iy + 1][iz + 1];
                }
                else
                {
                    iiNodes[L].iiun = 0;
                }

                if (iz < NK && iy < NJ && ix < NI)
                {
                    iiNodes[L].iiune = ii[ix + 1][iy + 1][iz + 1];
                }
                else
                {
                    iiNodes[L].iiune = 0;
                }

                if (iz < NK && iy < NJ && ix > 1)
                {
                    iiNodes[L].iiuwn = ii[ix - 1][iy + 1][iz + 1];
                }
                else
                {
                    iiNodes[L].iiuwn = 0;
                }

                if (iz < NK && iy > 1)
                {
                    iiNodes[L].iius = ii[ix][iy - 1][iz + 1];
                }
                else
                {
                    iiNodes[L].iius = 0;
                }

                if (iz < NK && iy > 1 && ix < NI)
                {
                    iiNodes[L].iiuse = ii[ix + 1][iy - 1][iz + 1];
                }
                else
                {
                    iiNodes[L].iiuse = 0;
                }

                if (iz < NK && iy > 1 && ix > 1)
                {
                    iiNodes[L].iiusw = ii[ix - 1][iy - 1][iz + 1];
                }
                else
                {
                    iiNodes[L].iiusw = 0;
                }

                ///// D - plane ///////////////////////////////////////////////////////////
                if (iz > 1)
                {
                    iiNodes[L].iid = ii[ix][iy][iz - 1];
                }
                else
                {
                    iiNodes[L].iid = 0;
                }

                if (iz > 1 && ix < NI)
                {
                    iiNodes[L].iide = ii[ix + 1][iy][iz - 1];
                }
                else
                {
                    iiNodes[L].iide = 0;
                }

                if (iz > 1 && ix > 1)
                {
                    iiNodes[L].iidw = ii[ix - 1][iy][iz - 1];
                }
                else
                {
                    iiNodes[L].iidw = 0;
                }

                if (iz > 1 && iy < NJ)
                {
                    iiNodes[L].iidn = ii[ix][iy + 1][iz - 1];
                }
                else
                {
                    iiNodes[L].iidn = 0;
                }

                if (iz > 1 && iy < NJ && ix < NI)
                {
                    iiNodes[L].iidne = ii[ix + 1][iy + 1][iz - 1];
                }
                else
                {
                    iiNodes[L].iidne = 0;
                }

                if (iz > 1 && iy < NJ && ix > 1)
                {
                    iiNodes[L].iidnw = ii[ix - 1][iy + 1][iz - 1];
                }
                else
                {
                    iiNodes[L].iidnw = 0;
                }

                if (iz > 1 && iy > 1)
                {
                    iiNodes[L].iids = ii[ix][iy - 1][iz - 1];
                }
                else
                {
                    iiNodes[L].iids = 0;
                }

                if (iz > 1 && iy > 1 && ix < NI)
                {
                    iiNodes[L].iidse = ii[ix + 1][iy - 1][iz - 1];
                }
                else
                {
                    iiNodes[L].iidse = 0;
                }

                if (iz > 1 && iy > 1 && ix > 1)
                {
                    iiNodes[L].iidsw = ii[ix - 1][iy - 1][iz - 1];
                }
                else
                {
                    iiNodes[L].iidsw = 0;
                }
                dum++;
            }
        }
    }
}

StructuredMesh::~StructuredMesh()
{
}

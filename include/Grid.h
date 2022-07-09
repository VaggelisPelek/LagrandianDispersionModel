#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

struct Nodes
{
    int iie, iiw, iin, iis, iinw, iine, iisw, iise, iiu, iiuw, iiue, iius, iiun,
        iiuwn, iiune, iiusw, iiuse, iid, iidw, iide, iids, iidn, iidnw, iidne, iidsw, iidse;
};

class StructuredMesh
{
private:
    int NON;
    int NI, NJ, NK;
    int N, L, ISTART0, IEND0, JSTART0, JEND0, KSTART0, KEND0;

    //////////////////////////////
    std::vector<Nodes> iiNodes;
    vector<double> XGrid, YGrid, ZGrid, DXEP, DXPW, DYPS, DYNP, DZPU, DZDP, SDU, SNS, SEW;
    vector<vector<vector<int>>> ii;
    vector<vector<vector<int>>> iwall3d, iwall;
    std::vector<int> iwall1d;

public:
    StructuredMesh();
    //read file with nodes
    void readGrid(const string s);
    // Define Your Mesh Specs
    void DefineMesh(); //const string sstring fileName
    //initialize lists and tables
    void initializeVariables();
    // string indexing
    void StringIndexing();

    //Get input Details
    int getNON()
    {
        NON = NI * NJ * NK;
        return NON;
    }
    int getNI() { return NI; }
    int getNJ() { return NJ; }
    int getNK() { return NK; }
    std::vector<Nodes> getiiNodes() { return iiNodes; }
    std::vector<vector<vector<int>>> getii() { return ii; }
    std::vector<vector<vector<int>>> getiwall() { return iwall; }

    std::vector<double> getSDU() { return SDU; }
    std::vector<double> getSNS() { return SNS; }
    std::vector<double> getSEW() { return SEW; }

    std::vector<double> getXgrid() { return XGrid; }
    std::vector<double> getYgrid() { return YGrid; }
    std::vector<double> getZgrid() { return ZGrid; }
    std::vector<int> getIwall1D() { return iwall1d; }
    ~StructuredMesh();
};
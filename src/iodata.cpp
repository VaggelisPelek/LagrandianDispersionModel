#include "iodata.h"

iodata::iodata()
{
}

void iodata::readVel(string fileName)
{
    double u, v, w, te, ed;
    //// ASCII ////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    ifstream intVel;
    intVel.open(fileName.c_str(), std::ios::in);
    ifstream intKE("inputData/fluidProp.dat", ios::in);

    int dum = 0;
    for (int i = 1; i < NI - 1; i++)
    {
        for (int j = 1; j < NJ - 1; j++)
        {
            for (int k = 1; k < NK - 1; k++)
            {
                intVel >> u;
                U[i][j][k] = u;
                U1d[dum] = u;

                dum++;
            }
        }
    }

    int duum = 0;
    for (int i = 1; i < NI - 1; i++)
    {
        for (int j = 1; j < NJ - 1; j++)
        {
            for (int k = 1; k < NK - 1; k++)
            {
                intVel >> v;
                V[i][j][k] = v;
                V1d[duum] = v;
                duum++;
            }
        }
    }
    int duuum = 0;
    for (int i = 1; i < NI - 1; i++)
    {
        for (int j = 1; j < NJ - 1; j++)
        {
            for (int k = 1; k < NK - 1; k++)
            {
                intVel >> w;
                W[i][j][k] = w;
                W1d[duuum] = w;

                duuum++;
            }
        }
    }

    int dumED = 0;
    for (int i = 1; i < NI - 1; i++)
    {
        for (int j = 1; j < NJ - 1; j++)
        {
            for (int k = 1; k < NK - 1; k++)
            {
                intKE >> ed;
                // ED[i][j][k] = ed;
                ED1d[dumED] = ed;
                dumED++;
            }
        }
    }
    int dumTE = 0;
    for (int i = 1; i < NI - 1; i++)
    {
        for (int j = 1; j < NJ - 1; j++)
        {
            for (int k = 1; k < NK - 1; k++)
            {
                intKE >> te;
                TE1d[dumTE] = te;
                dumTE++;
            }
        }
    }
}

// Read from file FLUIDINPUT
void iodata::readFluidPro(string fileName)
{
    ifstream intflow;
    intflow.open(fileName.c_str(), std::ios::in); //// ASCII
    int i = 0;
    vector<double> content;
    content.resize(2);
    if (!fileName.c_str())
        cout << "File cannot be opened";
    if (!intflow)
    {
        cout << "Unable to open FLUIDINPUT";
    }

    std::string line;
    while (getline(intflow, line))
    {
        std::string::size_type n = line.find("//");
        if (n != std::string::npos)
            line.erase(n);

        std::istringstream iss(line);

        double num;
        while (iss >> num)
            content[i] = num;
        i++;
    }
    den = content[0];
    visc = content[1];
    intflow.close();
}

// Read from file PRTCLINPUT
void iodata::readPartPro(string fileName)
{
    ifstream intPart;
    intPart.open(fileName.c_str(), std::ios::in); //// ASCII
    vector<double> content2;
    // content2.resize(47);
    double cont[49];
    int i = 0;
    if (!intPart)
    {
        cout << "Unable to open PRTCLINPUT";
    }

    std::string line;
    while (getline(intPart, line))
    {
        std::string::size_type n = line.find("//");
        if (n != std::string::npos)
            line.erase(n);

        std::istringstream iss(line);

        double num;
        while (iss >> num)
            content2.push_back(num);
        i++;
    }

    intPart.close();

    KLM = content2[0];
    PDEN0 = content2[1];
    PANU = content2[2];
    PAE = content2[3];
    SUDEN = content2[4];
    SUNU = content2[5];
    SUE = content2[6];
    DGAM = content2[7];
    YELAST = content2[8];
    DMA = content2[9];
    DNA = content2[10];
    DPA = content2[11];
    NPDIA = content2[12];
    int dum = 13;
    for (int i = 1; i < 8 + 1; i++) //npdia = 8 in PRTCLINPUT
    {
        pdiami.push_back(content2[dum]);
        dum++;
    }
    int duum = 21;
    for (int i = 0; i < 8; i++) //npdia = 8 in PRTCLINPUT
    {
        PMFRA.push_back(content2[duum]);
        duum++;
    }
    ZEMIS1 = content2[29];
    ZEMIS2 = content2[30];
    IEMIS = content2[31];
    JEMIS = content2[32];
    TURINTU = content2[33];
    TURINTV = content2[34];
    TURINTW = content2[35];
    THCKP = content2[36];
    THCKG = content2[37];
    RGAS = content2[38];
    SUDEN2 = content2[39];
    SUNU2 = content2[40];
    SUE2 = content2[41];
    DGAM2 = content2[42];
    A0 = content2[43];
    LLL = content2[44];
    GRAV = content2[45];
    DisCoef = content2[46];
    MINTSPSV = content2[47];
    MAXTSPSV = content2[48];
    MAXTSPENTRA = content2[49];
}

void iodata::resizeVec()
{

    X.resize(NI * NJ * NK);
    Y.resize(NI * NJ * NK);
    Z.resize(NI * NJ * NK);
    int nijk = (NI - 2) * (NJ - 2) * (NK - 2);
    U.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    V.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    W.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    U1d.resize(nijk);
    V1d.resize(nijk);
    W1d.resize(nijk);
    TE1d.resize(nijk);
    ED1d.resize(nijk);
    TE.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    ED.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    TEMP.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    UPAST.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    VPAST.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    WPAST.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    VIS.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    P.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
    DEN.resize(NI, vector<vector<double>>(NJ, vector<double>(NK)));
}

iodata::~iodata()
{
}
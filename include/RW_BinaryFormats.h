/// Conversion ASCII .dat files to Binary format .bin
////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

struct Vectors
{
    std::vector<double> A, B, C, U, V, W, E, T;
    std::vector<int> D;
    bool saveGrid(const char *filename);
    bool saveFlow(const char *filename);
    bool loadGrid(const char *filename);
    bool loadGrid_test(const char *filename);
    bool loadFlow(const char *filename);
    bool operator==(const Vectors &rhs);
    std::vector<double> getXgrid() { return A; }
    std::vector<double> getYgrid() { return B; }
    std::vector<double> getZgrid() { return C; }
    std::vector<int> getIwall1d() { return D; }
    std::vector<double> getUflow() { return U; }
    std::vector<double> getVflow() { return V; }
    std::vector<double> getWflow() { return W; }
    std::vector<double> getEDflow() { return E; }
    std::vector<double> getTEflow() { return T; }
};

void initialize_doubles(std::vector<double> &v, int size, vector<double> &Grid)
{
    v.resize(size);
    for (int n = 0; n < size; ++n)
    {
        v[n] = Grid[n];
        if (n % 1000 == 0)
            cout << Grid[n] << endl;
    }
}
void initialize_ints(std::vector<int> &v, int size, vector<int> &Grid)
{
    v.resize(size);
    for (int n = 0; n < size; ++n)
    {
        v[n] = Grid[n];
        if (n % 1000 == 0)
            cout << Grid[n] << endl;
    }
}

bool Vectors::saveGrid(const char *filename)
{
    std::ofstream out(filename, std::ios::binary);
    int a = A.size(), b = B.size(),
        c = C.size(), d = D.size();
    out.write(reinterpret_cast<const char *>(&a), sizeof(a));
    out.write(reinterpret_cast<const char *>(&b), sizeof(b));
    out.write(reinterpret_cast<const char *>(&c), sizeof(c));

    out.write(reinterpret_cast<const char *>(&A[0]), sizeof(double) * A.size());
    out.write(reinterpret_cast<const char *>(&B[0]), sizeof(double) * B.size());
    out.write(reinterpret_cast<const char *>(&C[0]), sizeof(double) * C.size());
    return true;
}

bool Vectors::saveFlow(const char *filename)
{
    std::ofstream out(filename, std::ios::binary);
    int a = U.size(), b = V.size(),
        c = W.size(), d = E.size(), e = T.size();
    out.write(reinterpret_cast<const char *>(&a), sizeof(a));
    out.write(reinterpret_cast<const char *>(&b), sizeof(b));
    out.write(reinterpret_cast<const char *>(&c), sizeof(c));
    out.write(reinterpret_cast<const char *>(&d), sizeof(d));
    out.write(reinterpret_cast<const char *>(&e), sizeof(e));

    out.write(reinterpret_cast<const char *>(&U[0]), sizeof(double) * U.size());
    out.write(reinterpret_cast<const char *>(&V[0]), sizeof(double) * V.size());
    out.write(reinterpret_cast<const char *>(&W[0]), sizeof(double) * W.size());
    out.write(reinterpret_cast<const char *>(&E[0]), sizeof(double) * E.size());
    out.write(reinterpret_cast<const char *>(&T[0]), sizeof(double) * T.size());
    return true;
}

bool Vectors::loadGrid(const char *filename)
{
    std::ifstream in(filename, std::ios::binary);
    int a, b, c, d; //, c, d, e, f;
    in.read(reinterpret_cast<char *>(&a), sizeof(a));
    in.read(reinterpret_cast<char *>(&b), sizeof(b));
    in.read(reinterpret_cast<char *>(&c), sizeof(c));
    A.resize(a);
    B.resize(b);
    C.resize(c);

    in.read(reinterpret_cast<char *>(&A[0]), sizeof(double) * A.size());
    in.read(reinterpret_cast<char *>(&B[0]), sizeof(double) * B.size());
    in.read(reinterpret_cast<char *>(&C[0]), sizeof(double) * C.size());
    return true;
}

bool Vectors::loadFlow(const char *filename)
{
    std::ifstream in(filename, std::ios::binary);
    int a, b, c, d, e; //, f;
    in.read(reinterpret_cast<char *>(&a), sizeof(a));
    in.read(reinterpret_cast<char *>(&b), sizeof(b));
    in.read(reinterpret_cast<char *>(&c), sizeof(c));
    in.read(reinterpret_cast<char *>(&d), sizeof(d));
    in.read(reinterpret_cast<char *>(&e), sizeof(e));
    U.resize(a);
    V.resize(b);
    W.resize(c);
    E.resize(d);
    T.resize(e);

    in.read(reinterpret_cast<char *>(&U[0]), sizeof(double) * U.size());
    in.read(reinterpret_cast<char *>(&V[0]), sizeof(double) * V.size());
    in.read(reinterpret_cast<char *>(&W[0]), sizeof(double) * W.size());
    in.read(reinterpret_cast<char *>(&E[0]), sizeof(double) * E.size());
    in.read(reinterpret_cast<char *>(&T[0]), sizeof(double) * T.size());
    return true;
}
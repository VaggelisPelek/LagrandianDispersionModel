#pragma once
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
const double pi = atan(1) * 4;

class particles
{
private:
public:
    int lx, ly, lz, cell;
    double xp, yp, zp, up, vp, wp, visc, den, pdiam, mass, sumtime;
    string idNum;
    particles() {}
    particles(int lx, int ly, int lz, double x, double y, double z, double u, double v, double w, double pdiam, string id, int cell)
    {
        this->lx = lx;
        this->ly = ly;
        this->lz = lz;
        this->xp = x;
        this->yp = y;
        this->zp = z;
        this->up = u;
        this->vp = v;
        this->wp = w;
        this->pdiam = pdiam;
        this->idNum = id;
        this->cell = cell;
    }

    void setParViscosity(double viscosity) { visc = viscosity; }
    void setParDensity(double Density) { den = Density; }
    void setParMass() { mass = den * pi / 3 * pdiam * pdiam; }

    double getXf() { return xp; }
    double getYf() { return yp; }
    double getZf() { return zp; }
    double getUf() { return up; }
    double getVf() { return vp; }
    double getWf() { return wp; }
};
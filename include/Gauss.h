#pragma once
#include <iostream>
#include <cstring>
#include <math.h>
#include <vector>
#include <random>
using namespace std;
constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;

long double RNDG(long double DSDS)
{
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);
    double RAND = distr(eng);
    const int nrolls = 10000;  // number of experiments
    const int nstars = 95;     // maximum number of stars to distribute
    const int nintervals = 10; // number of intervals
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return number;
}

double GAUSS(double DSEED, vector<double> PINAX)
{
    double GGG;
    double PPP = DSEED; //RAND;
    if (PPP < 0.001)
        PPP = 0.001;
    if (PPP > 0.999)
        PPP = 0.999;
    int KKK4 = int(PPP * 1000.);
    double XIN = PINAX[KKK4];
    double XOUT = PINAX[KKK4 + 1];
    double PIN = 1. / 1000. * KKK4;
    double POUT = 1. / 1000. * (KKK4 + 1);
    double X = (PPP * (XOUT - XIN) - PIN * XOUT + POUT * XIN) / (POUT - PIN);
    return GGG = X;
}
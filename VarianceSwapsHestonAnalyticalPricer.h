#ifndef VARIANCESWAPSHESTONANALYTICALPRICER_H
#define VARIANCESWAPSHESTONANALYTICALPRICER_H

#include "VarianceSwapsPricer.h"
#include <complex>

class VarianceSwapsHestonAnalyticalPricer : public VarianceSwapsHestonPricer
{
private:

public:
    VarianceSwapsHestonAnalyticalPricer(const HestonModel& hestonModel);
    ~VarianceSwapsHestonAnalyticalPricer();
    complex<double> aTerm (double omega) const ;
    complex<double> bTerm (double omega) const ;
    complex<double> gTerm (double omega) const ;
    complex<double> FunctionC (double tau, double omega) const;
    static complex<double> FunctionD (double tau, double omega) const;
    complex<double> FunctionDPrime (double tau, double omega) const;
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif

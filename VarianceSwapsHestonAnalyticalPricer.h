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
    std::complex<double> aTerm (double omega) const ;
    std::complex<double> bTerm (double omega) const ;
    std::complex<double> gTerm (double omega) const ;
    std::complex<double> FunctionC (double tau, double omega) const;
    std::complex<double> FunctionD (double tau, double omega) const;
    std::complex<double> FunctionDPrime (double tau, double omega) const;
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif

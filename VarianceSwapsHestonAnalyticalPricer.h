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
    double qtilde () const;
    std::complex<double> functionC (double tau, double omega) const;
    std::complex<double> functionD (double tau, double omega) const;
    std::complex<double> functionDPrime (double tau, double omega) const;
    std::complex<double> functionCPrime (double tau, double omega) const;
    std::complex<double> functionDSecond (double tau, double omega) const;
    std::complex<double> functionCSecond (double tau, double omega) const;
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif

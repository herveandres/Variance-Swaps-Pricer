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
    std::complex<double> functionC (double tau, double omega) const;
    std::complex<double> functionD (double tau, double omega) const;
    std::complex<double> functionDPrime (double tau, double omega) const;
    std::complex<double> functionCPrime (double tau, double omega) const;
    std::complex<double> functionDSecond (double tau, double omega) const;
    std::complex<double> functionCSecond (double tau, double omega) const;
    double qtilde () const;
    double cTerm (size_t i, const VarianceSwap& varianceSwap) const;
    double wTerm (size_t i, const VarianceSwap& varianceSwap) const;
    double u_1Term (const VarianceSwap& varianceSwap) const;
    double u_iTerm (size_t i, const VarianceSwap& varianceSwap) const;
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif

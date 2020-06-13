#ifndef VARIANCESWAPSHESTONANALYTICALPRICER_H
#define VARIANCESWAPSHESTONANALYTICALPRICER_H

#include "VarianceSwapsPricer.h"
#include <complex>

class VarianceSwapsHestonAnalyticalPricer : public VarianceSwapsHestonPricer
{
private:
    HestonModel* hestonModel_;
    std::complex<double> aTerm (double omega) const ;
    std::complex<double> bTerm (double omega) const ;
    std::complex<double> gTerm (double omega) const ;
    std::complex<double> functionC (double tau, double omega) const;
    std::complex<double> functionD (double tau, double omega) const;
    std::complex<double> functionDPrime (double tau, double omega) const;
    std::complex<double> functionCPrime (double tau, double omega) const;
    std::complex<double> functionDSecond (double tau, double omega) const;
    std::complex<double> functionCSecond (double tau, double omega) const;
    double qtildeTerm () const;
    double cTerm (double t) const;
    double wTerm (double t) const;
    double u1Term (double t1, double t2) const;
    double uiTerm (double t1, double t2) const;
public:
    VarianceSwapsHestonAnalyticalPricer(const HestonModel& hestonModel);
    ~VarianceSwapsHestonAnalyticalPricer();
    double price(const VarianceSwap& varianceSwap) const override;
    double continousPrice(const VarianceSwap& varianceSwap);
};

#endif

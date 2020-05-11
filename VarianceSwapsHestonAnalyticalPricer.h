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
    complex<double> aTerm (const double omega) const ;
    complex<double> bTerm (const double omega) const ;
    complex<double> gTerm (const double omega) const ;
    complex<double> FunctionC (const double tau, const double omega) const;
    complex<double> FunctionD (const double tau, const double omega) const;
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif

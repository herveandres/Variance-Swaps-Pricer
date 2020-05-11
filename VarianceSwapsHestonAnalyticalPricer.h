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
    complex<double> a_term (const double omega) const ;
    complex<double> b_term (const double omega) const ;
    complex<double> g_term (const double omega) const ;
    complex<double> Function_C (const double tau, const double omega) const;
    complex<double> Function_D (const double tau, const double omega) const;
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif

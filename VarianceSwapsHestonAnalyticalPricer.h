#ifndef VARIANCESWAPSHESTONANALYTICALPRICER_H
#define VARIANCESWAPSHESTONANALYTICALPRICER_H

#include "VarianceSwapsPricer.h"

class VarianceSwapsHestonAnalyticalPricer : public VarianceSwapsHestonPricer
{
private:

public:
    VarianceSwapsHestonAnalyticalPricer(const HestonModel& hestonModel);
    ~VarianceSwapsHestonAnalyticalPricer();
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif
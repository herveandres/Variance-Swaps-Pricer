#ifndef VARIANCESWAPSHESTONMONTECARLOPRICER_H
#define VARIANCESWAPSHESTONMONTECARLOPRICER_H

#include "VarianceSwapsPricer.h"
#include "HestonLogSpotPathSimulator.h"

class VarianceSwapsHestonMonteCarloPricer : public VarianceSwapsHestonPricer
{
private:
    HestonLogSpotPathSimulator* hestonPathSimulator_;
public:
    VarianceSwapsHestonMonteCarloPricer(const HestonModel& hestonModel,
                                        const HestonLogSpotPathSimulator& hestonPathSimulator);
    ~VarianceSwapsHestonMonteCarloPricer();
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif
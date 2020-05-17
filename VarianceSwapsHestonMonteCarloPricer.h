#ifndef VARIANCESWAPSHESTONMONTECARLOPRICER_H
#define VARIANCESWAPSHESTONMONTECARLOPRICER_H

#include "VarianceSwapsPricer.h"
#include "HestonLogSpotPathSimulator.h"

class VarianceSwapsHestonMonteCarloPricer : public VarianceSwapsHestonPricer
{
private:
    HestonLogSpotPathSimulator* hestonPathSimulator_;
    size_t nbSimulations_;
    double pathPrice(std::vector<double> path, double maturity) const;
public:
    VarianceSwapsHestonMonteCarloPricer(const HestonModel& hestonModel,
                                        const HestonLogSpotPathSimulator& hestonPathSimulator,
                                        std::size_t nbSimulations);
    ~VarianceSwapsHestonMonteCarloPricer();
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif
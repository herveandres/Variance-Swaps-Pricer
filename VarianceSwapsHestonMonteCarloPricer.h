#ifndef VARIANCESWAPSHESTONMONTECARLOPRICER_H
#define VARIANCESWAPSHESTONMONTECARLOPRICER_H

#include "VarianceSwapsPricer.h"
#include "HestonLogSpotPathSimulator.h"

class VarianceSwapsHestonMonteCarloPricer : public VarianceSwapsHestonPricer
{
private:
    HestonLogSpotPathSimulator* hestonPathSimulator_;
    size_t nbSimulations_;
    //Method computing the price of a variance swap for a given path of the underlying
    double pathPrice(std::vector<double> path, double maturity) const;
public:
    VarianceSwapsHestonMonteCarloPricer(const HestonLogSpotPathSimulator& hestonPathSimulator,
                                        std::size_t nbSimulations);
    ~VarianceSwapsHestonMonteCarloPricer();
    VarianceSwapsHestonMonteCarloPricer(const VarianceSwapsHestonMonteCarloPricer& mcPricer);
    VarianceSwapsHestonMonteCarloPricer& operator=(
                        const VarianceSwapsHestonMonteCarloPricer& mcPricer);
    //Method returning the Monte Carlo price of the variance swap given as argument
    double price(const VarianceSwap& varianceSwap) const override;
};

#endif
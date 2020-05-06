#include "VarianceSwapsHestonMonteCarloPricer.h"

VarianceSwapsHestonMonteCarloPricer::VarianceSwapsHestonMonteCarloPricer
                                (const HestonModel& hestonModel,
                                const HestonLogSpotPathSimulator& hestonPathSimulator):
            VarianceSwapsHestonPricer(hestonModel), hestonPathSimulator_(hestonPathSimulator.clone())
{

}

VarianceSwapsHestonMonteCarloPricer::~VarianceSwapsHestonMonteCarloPricer()
{
    delete hestonPathSimulator_;
}

double VarianceSwapsHestonMonteCarloPricer::price(const VarianceSwap& varianceSwap) const{
    //A compléter
    return 0;
}
#include "VarianceSwapsHestonMonteCarloPricer.h"

VarianceSwapsHestonMonteCarloPricer::VarianceSwapsHestonMonteCarloPricer(
                                            const HestonModel& hestonModel):
                                    VarianceSwapsHestonPricer(hestonModel)
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
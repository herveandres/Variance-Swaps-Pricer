#include "VarianceSwapsHestonAnalyticalPricer.h"

VarianceSwapsHestonAnalyticalPricer::VarianceSwapsHestonAnalyticalPricer(
                                            const HestonModel& hestonModel):
                                            VarianceSwapsHestonPricer(hestonModel)
{
    
}

VarianceSwapsHestonAnalyticalPricer::~VarianceSwapsHestonAnalyticalPricer()
{

}

double VarianceSwapsHestonAnalyticalPricer::price(const VarianceSwap& varianceSwap) const{
    //A compléter
    return 0;
}
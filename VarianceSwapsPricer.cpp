#include "VarianceSwapsPricer.h"



VarianceSwapsPricer::~VarianceSwapsPricer()
{

}

VarianceSwapsHestonPricer::VarianceSwapsHestonPricer(const HestonModel& hestonModel):
                                        hestonModel_(new HestonModel(hestonModel))
{
    
}

VarianceSwapsHestonPricer::~VarianceSwapsHestonPricer()
{
    delete hestonModel_;
}
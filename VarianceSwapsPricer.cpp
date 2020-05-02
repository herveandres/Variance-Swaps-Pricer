#include "VarianceSwapsPricer.h"


VarianceSwapsPricer::VarianceSwapsPricer(const Model& model): model_(model.clone())
{

}

VarianceSwapsPricer::~VarianceSwapsPricer()
{
    delete model_;
}

VarianceSwapsHestonPricer::VarianceSwapsHestonPricer(const HestonModel& hestonModel):
                                       VarianceSwapsPricer(hestonModel)
{
    
}

VarianceSwapsHestonPricer::~VarianceSwapsHestonPricer()
{

}
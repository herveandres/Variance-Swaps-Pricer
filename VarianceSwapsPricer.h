#ifndef VARIANCESWAPSPRICER_H
#define VARIANCESWAPSPRICER_H

#include "VarianceSwap.h"
#include "Model.h"

//Abstract class
class VarianceSwapsPricer
{
protected:

public:
    virtual ~VarianceSwapsPricer();
    virtual double price(const VarianceSwap& varianceSwap) const = 0;
};

//Abstract class
class VarianceSwapsHestonPricer : public VarianceSwapsPricer
{
public:
    virtual ~VarianceSwapsHestonPricer();
    virtual double price(const VarianceSwap& varianceSwap) const = 0;
};

#endif

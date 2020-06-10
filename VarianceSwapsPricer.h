#ifndef VARIANCESWAPSPRICER_H
#define VARIANCESWAPSPRICER_H

#include "VarianceSwap.h"
#include "Model.h"

class VarianceSwapsPricer
{
protected:

public:
    virtual ~VarianceSwapsPricer();
    virtual double price(const VarianceSwap& varianceSwap) const = 0;
};

class VarianceSwapsHestonPricer : public VarianceSwapsPricer
{
public:
    virtual ~VarianceSwapsHestonPricer();
    virtual double price(const VarianceSwap& varianceSwap) const = 0;
};

#endif

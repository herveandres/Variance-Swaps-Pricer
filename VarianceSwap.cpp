#include "VarianceSwap.h"
#include "MathFunctions.h"

VarianceSwap::VarianceSwap(double maturity, std::size_t nbOfObservations) 
{
    dates_ = MathFunctions::buildLinearSpace(0,maturity,nbOfObservations);
}


std::vector<double> VarianceSwap::getDates() const
{
    return dates_;
}
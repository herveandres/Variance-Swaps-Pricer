#include "VarianceSwap.h"

VarianceSwap::VarianceSwap(double maturity, size_t nbOfObservations) 
    : maturity_(maturity), nbOfObservations_(nbOfObservations)
{
    double delta_t = maturity/nbOfObservations;
    for(size_t i = 0; i <= nbOfObservations; i++)
    {
        dates_.push_back(i*delta_t);
    }
}
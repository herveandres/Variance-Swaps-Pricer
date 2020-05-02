#include "HestonVariancePathSimulator.h"

HestonVariancePathSimulator::HestonVariancePathSimulator(const HestonModel& hestonModel):
                            PathSimulator(hestonModel)
{

}

HestonVariancePathSimulator::~HestonVariancePathSimulator()
{
    
}

TruncatedGaussianScheme::TruncatedGaussianScheme(const HestonModel& hestonModel):
                        HestonVariancePathSimulator(hestonModel)
{

}

QuadraticExponentialScheme::QuadraticExponentialScheme(const HestonModel& hestonModel):
                        HestonVariancePathSimulator(hestonModel)
{

}

TruncatedGaussianScheme* TruncatedGaussianScheme::clone() const
{
    return new TruncatedGaussianScheme(*this);
}

QuadraticExponentialScheme* QuadraticExponentialScheme::clone() const
{
    return new QuadraticExponentialScheme(*this);
}
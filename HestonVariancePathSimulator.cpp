#include "HestonVariancePathSimulator.h"

HestonVariancePathSimulator::HestonVariancePathSimulator(const HestonModel& hestonModel):
                            PathSimulator(hestonModel)
{

}

HestonVariancePathSimulator::~HestonVariancePathSimulator()
{
    
}

TruncatedGaussian::TruncatedGaussian(const HestonModel& hestonModel):
                        HestonVariancePathSimulator(hestonModel)
{

}

QuadraticExponential::QuadraticExponential(const HestonModel& hestonModel):
                        HestonVariancePathSimulator(hestonModel)
{

}

TruncatedGaussian* TruncatedGaussian::clone() const
{
    return new TruncatedGaussian(*this);
}

QuadraticExponential* QuadraticExponential::clone() const
{
    return new QuadraticExponential(*this);
}
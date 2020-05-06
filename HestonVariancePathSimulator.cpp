#include "HestonVariancePathSimulator.h"

HestonVariancePathSimulator::HestonVariancePathSimulator(
                                        const std::vector<double>& timePoints,
                                        const HestonModel& hestonModel):
                            PathSimulator(hestonModel.getInitialVolatility(),
                                          timePoints),
                            hestonModel_(new HestonModel(hestonModel))
{

}

HestonVariancePathSimulator::~HestonVariancePathSimulator()
{
    delete hestonModel_;
}

TruncatedGaussianScheme::TruncatedGaussianScheme(const std::vector<double>& timePoints,
                                                const HestonModel& hestonModel):
                        HestonVariancePathSimulator(timePoints,hestonModel)
{

}

TruncatedGaussianScheme::TruncatedGaussianScheme(const TruncatedGaussianScheme& truncatedGaussianScheme):
                    HestonVariancePathSimulator(truncatedGaussianScheme.timePoints_,
                                            *truncatedGaussianScheme.hestonModel_)
{

}

TruncatedGaussianScheme* TruncatedGaussianScheme::clone() const
{
    return new TruncatedGaussianScheme(*this);
}

double TruncatedGaussianScheme::nextStep(std::size_t currentIndex, double currentValue) const
{
    //A compléter
}

QuadraticExponentialScheme::QuadraticExponentialScheme(const std::vector<double>& timePoints,
                                                    const HestonModel& hestonModel):
                        HestonVariancePathSimulator(timePoints,hestonModel)
{

}

QuadraticExponentialScheme::QuadraticExponentialScheme(const QuadraticExponentialScheme&
                                                        quadraticExponentialScheme):
                        HestonVariancePathSimulator(quadraticExponentialScheme.timePoints_,
                                                *quadraticExponentialScheme.hestonModel_)
{

}

QuadraticExponentialScheme* QuadraticExponentialScheme::clone() const
{
    return new QuadraticExponentialScheme(*this);
}

double QuadraticExponentialScheme::nextStep(std::size_t currentIndex, double currentValue) const
{
    //A compléter
}

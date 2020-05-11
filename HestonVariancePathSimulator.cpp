#include <cmath>
#include "HestonVariancePathSimulator.h"
#include "MathFunctions.h"

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
                                                const HestonModel& hestonModel, double confidenceMultiplier):
                        HestonVariancePathSimulator(timePoints,hestonModel), 
                        confidenceMultiplier_(confidenceMultiplier)
{
    preComputations();
}

TruncatedGaussianScheme::TruncatedGaussianScheme(const TruncatedGaussianScheme& truncatedGaussianScheme):
                    HestonVariancePathSimulator(truncatedGaussianScheme.timePoints_,
                                            *truncatedGaussianScheme.hestonModel_),
                    confidenceMultiplier_(truncatedGaussianScheme.confidenceMultiplier_)
{
    preComputations();
}

TruncatedGaussianScheme* TruncatedGaussianScheme::clone() const
{
    return new TruncatedGaussianScheme(*this);
}

void TruncatedGaussianScheme::preComputations()
{
    double theta = hestonModel_->getMeanReversionLevel();
    double kappa = hestonModel_->getMeanReversionSpeed();
    double eps = hestonModel_->getVolOfVol();
    double delta;
    double expMinusKappaDelta;

    for(std::size_t i = 0; i < timePoints_.size()-1; i++)
    {
        delta = timePoints_[i+1] - timePoints_[i];
        expMinusKappaDelta = exp(-kappa*delta);
        k1_.push_back(expMinusKappaDelta);
        k2_.push_back(theta*(1-expMinusKappaDelta));
        k3_.push_back(eps*eps*expMinusKappaDelta*(1-expMinusKappaDelta)/kappa);
        k4_.push_back(theta*eps*eps*(1-expMinusKappaDelta)*(1-expMinusKappaDelta)/(2*kappa));
    }
}

double TruncatedGaussianScheme::nextStep(std::size_t currentIndex, double currentValue) const
{
    double m = k1_[currentIndex]*currentValue + k2_[currentIndex];
    double s2 = k3_[currentIndex]*currentValue + k4_[currentIndex];

}   

double TruncatedGaussianScheme::h(double r, double psi)
{
    double phi = MathFunctions::normalPDF(r);
    double Phi = MathFunctions::normalCDF(r);

    return r*phi+Phi*(1+r*r)-(1+psi)*(phi+r*Phi)*(phi+r*Phi);
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

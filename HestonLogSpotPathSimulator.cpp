#include <cmath>
#include "HestonLogSpotPathSimulator.h"

HestonLogSpotPathSimulator::HestonLogSpotPathSimulator(
                                const std::vector<double>& timePoints,
                                const HestonModel& hestonModel,
                                const HestonVariancePathSimulator& variancePathSimulator):
                                PathSimulator(log(hestonModel.getInitialAssetValue()),
                                              timePoints),
                                hestonModel_(new HestonModel(hestonModel)),
                                variancePathSimulator_(variancePathSimulator.clone())
{

}

HestonLogSpotPathSimulator::~HestonLogSpotPathSimulator()
{
    delete hestonModel_;
    delete variancePathSimulator_;
}

BroadieKayaScheme::BroadieKayaScheme(const std::vector<double>& timePoints,
                    const HestonModel& hestonModel,
                    const HestonVariancePathSimulator& variancePathSimulator):
        HestonLogSpotPathSimulator(timePoints, hestonModel, variancePathSimulator)
{

}

BroadieKayaScheme::BroadieKayaScheme(const BroadieKayaScheme& broadieKayaScheme):
        HestonLogSpotPathSimulator(broadieKayaScheme.timePoints_,
                                    *broadieKayaScheme.hestonModel_,
                                    *broadieKayaScheme.variancePathSimulator_)
{

}

BroadieKayaScheme* BroadieKayaScheme::clone() const{
    return new BroadieKayaScheme(*this);
}

double BroadieKayaScheme::nextStep(std::size_t currentIndex, double currentValue) const
{
    //A compléter
}


void DiscretizationScheme::preComputations()
{
    double rho = hestonModel_->getCorrelation();
    double theta = hestonModel_->getMeanReversionLevel();
    double kappa = hestonModel_->getMeanReversionSpeed();
    double eps = hestonModel_->getVolOfVol();
    double delta;
    for(std::size_t i = 0; i < timePoints_.size()-1; i++)
    {
        delta = timePoints_[i+1] - timePoints_[i];
        k0_.push_back(-rho*kappa*theta*delta/eps);
        k1_.push_back(gamma1_*delta*(kappa*rho/eps-0.5)-rho/eps);
        k2_.push_back(gamma2_*delta*(kappa*rho/eps-0.5)+rho/eps);
        k3_.push_back(gamma1_*delta*(1.-rho*rho));
        k4_.push_back(gamma2_*delta*(1.-rho*rho));
    }
}

// void DiscretizationScheme::setGamma(double gamma1, double gamma2){
//     gamma1_ = gamma1;
//     gamma2_ = gamma2;
// }

DiscretizationScheme::DiscretizationScheme(const DiscretizationScheme& broadieKayaScheme):
        HestonLogSpotPathSimulator(broadieKayaScheme.timePoints_,
                                    *broadieKayaScheme.hestonModel_,
                                    *broadieKayaScheme.variancePathSimulator_)
{

}

DiscretizationScheme* DiscretizationScheme::clone() const{
    return new DiscretizationScheme(*this);
}

double DiscretizationScheme::nextStep(std::size_t currentIndex, double currentValue) const
{
    //A compléter
}


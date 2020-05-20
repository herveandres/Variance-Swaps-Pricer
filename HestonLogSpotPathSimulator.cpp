#include <cmath>
#include "HestonLogSpotPathSimulator.h"
#include "MathFunctions.h"

HestonLogSpotPathSimulator::HestonLogSpotPathSimulator(
                                const std::vector<double>& timePoints,
                                const HestonModel& hestonModel,
                                const HestonVariancePathSimulator& variancePathSimulator):
                                PathSimulator(std::log(hestonModel.getInitialAssetValue()),
                                              timePoints),
                                hestonModel_(new HestonModel(hestonModel)),
                                variancePathSimulator_(variancePathSimulator.clone())
{ 
    // Rajouter éventuellement une exception qui permet de gérer le cas où
    //les timePoints de variancePathSimulator correspondent pas à ceux du logspot
}

HestonLogSpotPathSimulator::~HestonLogSpotPathSimulator()
{
    delete hestonModel_;
    delete variancePathSimulator_;
}

std::vector<double> HestonLogSpotPathSimulator::path() const
{
    std::vector<double> logSpotPath {initialValue_};
    std::vector<double>  variancePath = variancePathSimulator_->path();
	for (std::size_t index = 0; index < timePoints_.size() - 1; ++index)
		logSpotPath.push_back(nextStep(index, logSpotPath[index], variancePath));

	return logSpotPath;
}


BroadieKayaScheme::BroadieKayaScheme(const std::vector<double>& timePoints,
                    const HestonModel& hestonModel,
                    const HestonVariancePathSimulator& variancePathSimulator,
                    double gamma1,
                    double gamma2):
        HestonLogSpotPathSimulator(timePoints, hestonModel, variancePathSimulator),
        gamma1_(gamma1), gamma2_(gamma2)
{
    preComputations();
}

BroadieKayaScheme::BroadieKayaScheme(const BroadieKayaScheme& broadieKayaScheme):
        HestonLogSpotPathSimulator(broadieKayaScheme.timePoints_,
                                    *broadieKayaScheme.hestonModel_,
                                    *broadieKayaScheme.variancePathSimulator_),
        gamma1_(broadieKayaScheme.gamma1_), gamma2_(broadieKayaScheme.gamma2_),
        k0_(broadieKayaScheme.k0_), k1_(broadieKayaScheme.k1_),
        k2_(broadieKayaScheme.k2_), k3_(broadieKayaScheme.k3_),
        k4_(broadieKayaScheme.k4_)
{

}

void BroadieKayaScheme::preComputations()
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
BroadieKayaScheme* BroadieKayaScheme::clone() const{
    return new BroadieKayaScheme(*this);
}

double BroadieKayaScheme::nextStep(std::size_t currentIndex, double currentValue, const std::vector<double>& variancePath) const
{
    double U = MathFunctions::simulateUniformRandomVariable();
    double Z = MathFunctions::normalCDFInverse(U);
    //Ajouter drift
    return currentValue+k0_[currentIndex]+k1_[currentIndex]*variancePath[currentIndex]
            +k2_[currentIndex]*variancePath[currentIndex+1]
            +std::sqrt(k3_[currentIndex]*variancePath[currentIndex]
                        +k4_[currentIndex]*variancePath[currentIndex+1])*Z;
}


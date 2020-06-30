#include <cmath>
#include "HestonLogSpotPathSimulator.h"
#include "MathFunctions.h"

HestonLogSpotPathSimulator::HestonLogSpotPathSimulator(
                                const HestonVariancePathSimulator& variancePathSimulator):
        PathSimulator(std::log(variancePathSimulator.getHestonModel().getInitialAssetValue()),
                        variancePathSimulator.getTimePoints()),
        variancePathSimulator_(variancePathSimulator.clone())
{ 

}

HestonLogSpotPathSimulator::HestonLogSpotPathSimulator(
                        const HestonLogSpotPathSimulator& logSpotPathSimulator):
        PathSimulator(logSpotPathSimulator.initialValue_, logSpotPathSimulator.timePoints_),
        variancePathSimulator_(logSpotPathSimulator.variancePathSimulator_->clone())               
{

}

HestonLogSpotPathSimulator::~HestonLogSpotPathSimulator()
{
    delete variancePathSimulator_;
}

HestonLogSpotPathSimulator& HestonLogSpotPathSimulator::operator=(
                        const HestonLogSpotPathSimulator& logSpotPathSimulator)
{
    if (this == &logSpotPathSimulator)
		return *this;
	else
	{
		delete variancePathSimulator_;												
		variancePathSimulator_ = logSpotPathSimulator.variancePathSimulator_->clone();	

        initialValue_ = logSpotPathSimulator.initialValue_;
        timePoints_ = logSpotPathSimulator.timePoints_;
	}
	return *this;
}

std::vector<double> HestonLogSpotPathSimulator::path() const
{
    std::vector<double> logSpotPath {initialValue_};
    //We compute the variance path at this stage once for all and we pass it to nextStep
    //It's not a class attribute in order to avoid that path is non const
    std::vector<double>  variancePath = variancePathSimulator_->path();
	for (std::size_t index = 0; index < timePoints_.size() - 1; ++index)
		logSpotPath.push_back(nextStep(index, logSpotPath[index], variancePath));

	return logSpotPath;
}


BroadieKayaScheme::BroadieKayaScheme(
                    const HestonVariancePathSimulator& variancePathSimulator,
                    double gamma1,
                    double gamma2):
        HestonLogSpotPathSimulator(variancePathSimulator),
        gamma1_(gamma1), gamma2_(gamma2)
{
    preComputations();
}

BroadieKayaScheme::BroadieKayaScheme(const BroadieKayaScheme& broadieKayaScheme):
        HestonLogSpotPathSimulator(*broadieKayaScheme.variancePathSimulator_),
        gamma1_(broadieKayaScheme.gamma1_), gamma2_(broadieKayaScheme.gamma2_),
        k0_(broadieKayaScheme.k0_), k1_(broadieKayaScheme.k1_),
        k2_(broadieKayaScheme.k2_), k3_(broadieKayaScheme.k3_),
        k4_(broadieKayaScheme.k4_)
{

}

void BroadieKayaScheme::preComputations()
{
    HestonModel hestonModel = variancePathSimulator_->getHestonModel();
    double rho = hestonModel.getCorrelation();
    double theta = hestonModel.getMeanReversionLevel();
    double kappa = hestonModel.getMeanReversionSpeed();
    double eps = hestonModel.getVolOfVol();
    double delta;

    //NB : we allow the time grid to be non-equidistant so that the computed quantities are time dependent
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
    double Z = MathFunctions::simulateGaussianRandomVariable();
    //Ajouter drift
    return currentValue
            +variancePathSimulator_->getHestonModel().getDrift()*(timePoints_[currentIndex+1] - timePoints_[currentIndex])
            +k0_[currentIndex]+k1_[currentIndex]*variancePath[currentIndex]
            +k2_[currentIndex]*variancePath[currentIndex+1]
            +std::sqrt(k3_[currentIndex]*variancePath[currentIndex]
                        +k4_[currentIndex]*variancePath[currentIndex+1])*Z;
}


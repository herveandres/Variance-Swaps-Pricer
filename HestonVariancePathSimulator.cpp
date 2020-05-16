#include <cmath>
#include <algorithm>
#include <iostream>
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
                                                const HestonModel& hestonModel, 
                                                double confidenceMultiplier,
                                                std::size_t psiGridSize,
                                                double initialGuess):
                        HestonVariancePathSimulator(timePoints,hestonModel), 
                        confidenceMultiplier_(confidenceMultiplier),
                        psiGrid_(std::vector<double>(psiGridSize)),
                        initialGuess_(initialGuess)
{
    preComputations();
}

TruncatedGaussianScheme::TruncatedGaussianScheme(const TruncatedGaussianScheme& truncatedGaussianScheme):
                    HestonVariancePathSimulator(truncatedGaussianScheme.timePoints_,
                                            *truncatedGaussianScheme.hestonModel_),
                    confidenceMultiplier_(truncatedGaussianScheme.confidenceMultiplier_),
                    psiGrid_(truncatedGaussianScheme.psiGrid_),
                    initialGuess_(truncatedGaussianScheme.initialGuess_)
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

    //Pre-computation of k1, k2, k3, k4 s.t. m = k1 V + k2 and s² = k3 V + k4
    for(std::size_t i = 0; i < timePoints_.size()-1; i++)
    {
        delta = timePoints_[i+1] - timePoints_[i];
        expMinusKappaDelta = exp(-kappa*delta);
        k1_.push_back(expMinusKappaDelta);
        k2_.push_back(theta*(1-expMinusKappaDelta));
        k3_.push_back(eps*eps*expMinusKappaDelta*(1-expMinusKappaDelta)/kappa);
        k4_.push_back(theta*eps*eps*(1-expMinusKappaDelta)*(1-expMinusKappaDelta)/(2*kappa));
    }


    //Pre-computation of f_mu and f_sigma
    double min = 1.0/(confidenceMultiplier_*confidenceMultiplier_);
    double max = eps*eps/(2*kappa*theta);
    psiGrid_ = MathFunctions::buildLinearSpace(min,max,psiGrid_.size());
    double r, psi, phi, Phi;
    for(std::size_t i = 0; i < psiGrid_.size(); i++)
    {
        psi = psiGrid_[i];
        r = MathFunctions::newtonMethod(initialGuess_,
                                        [psi](double r){return h(r,psi);},
                                        [psi](double r){return hPrime(r,psi);});
        phi = MathFunctions::normalPDF(r);
        Phi = MathFunctions::normalCDF(r);
        fmu_.push_back(r/(phi+r*Phi));
        fsigma_.push_back(pow(psi,-0.5)/(phi+r*Phi));
    }
}

double TruncatedGaussianScheme::nextStep(std::size_t currentIndex, double currentValue) const
{
    double m = k1_[currentIndex]*currentValue + k2_[currentIndex];
    double s2 = k3_[currentIndex]*currentValue + k4_[currentIndex];
    double psi = s2/(m*m);
    double mu, sigma;
    if(psi < 1.0/confidenceMultiplier_*confidenceMultiplier_)
    {
        mu = m;
        sigma = sqrt(s2);
    }
    else
    {
        std::size_t idxPhi = MathFunctions::binarySearch(psiGrid_,psi);
        double psi0 = psiGrid_[idxPhi], psi1 = psiGrid_[idxPhi+1];
        //Linear interpolation of fmu and fsigma using the pre-computed values
        double fmu = (fmu_[idxPhi]*(psi1-psi)+fmu_[idxPhi+1]*(psi-psi0))/(psi1-psi0);
        double fsigma = (fsigma_[idxPhi]*(psi1-psi)+fsigma_[idxPhi+1]*(psi-psi0))/(psi1-psi0);

        mu = fmu*m;
        sigma = fsigma*sqrt(s2);
    }
    double z = MathFunctions::simulateGaussianRandomVariable();
    return std::max(mu+sigma*z,0.0);
}   

double TruncatedGaussianScheme::h(double r, double psi)
{
    double phi = MathFunctions::normalPDF(r);
    double Phi = MathFunctions::normalCDF(r);

    return r*phi+Phi*(1+r*r)-(1+psi)*(phi+r*Phi)*(phi+r*Phi);
}

double TruncatedGaussianScheme::hPrime(double r, double psi)
{
    double phi = MathFunctions::normalPDF(r);
    double Phi = MathFunctions::normalCDF(r);

    return 2*phi+2*r*Phi-2*(1+psi)*Phi*(phi+r*Phi);
}


QuadraticExponentialScheme::QuadraticExponentialScheme(const std::vector<double>& timePoints,
                                                    const HestonModel& hestonModel, double psiC):
                        HestonVariancePathSimulator(timePoints,hestonModel), psiC_(psiC)
{


}

void QuadraticExponentialScheme::preComputations()
{
    double theta = hestonModel_->getMeanReversionLevel();
    double kappa = hestonModel_->getMeanReversionSpeed();
    double eps = hestonModel_->getVolOfVol();
    double delta;
    double expMinusKappaDelta;

    //Pre-computation of k1, k2, k3, k4 s.t. m = k1 V + k2 and s² = k3 V + k4
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

double QuadraticExponentialScheme::nextStep(std::size_t currentIndex, double currentValue) const{
    double m = k1_[currentIndex]*currentValue + k2_[currentIndex];
    double s2 = k3_[currentIndex]*currentValue + k4_[currentIndex];
    double psi = s2/(m*m);
    double p = (psi-1.)/(psi+1.);
    double U = MathFunctions::simulateGaussianRandomVariable();
    if (psi<psiC_){
        double temp_value = 2/psi;
        double b = sqrt(temp_value - 1. + sqrt(temp_value*(temp_value-1.)));
        double a = m/(1+b*b);
        double Zv = MathFunctions::normalCDF(U);
        return a*(b+Zv)*(b+Zv);
    }
    else {
        if (U<p){
            return 0.;
        }
        else{
            double beta = (1-p)/m;
            return log((1-p)/(1-U))/beta;
        }
    }
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

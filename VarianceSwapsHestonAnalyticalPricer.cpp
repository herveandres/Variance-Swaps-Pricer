#include <complex>
#include "VarianceSwapsHestonAnalyticalPricer.h"
#include "MathFunctions.h"

VarianceSwapsHestonAnalyticalPricer::VarianceSwapsHestonAnalyticalPricer(
                                            const HestonModel& hestonModel):
                                            VarianceSwapsHestonPricer(hestonModel)
{
    
}

VarianceSwapsHestonAnalyticalPricer::~VarianceSwapsHestonAnalyticalPricer()
{

}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
aTerm(double omega) const {
    std::complex<double> j(0.,1.);
    return hestonModel_->getMeanReversionSpeed() - hestonModel_->getCorrelation() * hestonModel_->getVolOfVol() * omega * j  ;
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
bTerm(double omega) const {
    std::complex<double> j(0.,1.);
    return sqrt(aTerm(omega)+ hestonModel_->getVolOfVol() * hestonModel_->getVolOfVol() * (j*omega + omega*omega) )  ;
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
gTerm(double omega) const {
    return (aTerm(omega)-bTerm(omega))/(aTerm(omega)+bTerm(omega));
}

double VarianceSwapsHestonAnalyticalPricer::
qtilde() const {
    return (2*hestonModel_->getMeanReversionSpeed()/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()));
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionC(double tau, double omega) const {
    std::complex<double> j(0.,1.);
    return tau * hestonModel_->getDrift() * (j * omega - 1.)+ hestonModel_->getMeanReversionSpeed() *  hestonModel_->getMeanReversionLevel() / (hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()) * ((aTerm(omega)-bTerm(omega))*tau-2.*log((1.-gTerm(omega)*exp(-bTerm(omega)*tau))/(1.-gTerm(omega)))) ;
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionD(double tau, double omega) const {
    std::complex<double> D_init;
    D_init=(aTerm(omega)-bTerm(omega))/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol());
    return D_init*(aTerm(omega)+bTerm(omega))*(1.-exp(-bTerm(omega)*tau))/(1.-gTerm(omega)*exp(-bTerm(omega)*tau));
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionDPrime(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionD(tau,omega);
    };
    return MathFunctions::differencesFinies(f,omega,tau);
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionCPrime(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionC(tau,omega);
    };
    return MathFunctions::differencesFinies(f,omega,tau);
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionCSecond(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionCPrime(tau,omega);
    };
    return MathFunctions::differencesFinies(f,omega,tau);
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionDSecond(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionDPrime(tau,omega);
    };
    return MathFunctions::differencesFinies(f,omega,tau);
}

double VarianceSwapsHestonAnalyticalPricer::price(const VarianceSwap& varianceSwap) const{
    //A compléter
    return 0;
}

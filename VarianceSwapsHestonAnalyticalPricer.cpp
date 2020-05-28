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
//useful variables to compute function C and D
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

//functions C and D and their derivatives
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

//useful terms in respect of the Chi2 law
double VarianceSwapsHestonAnalyticalPricer::
qtilde() const {
    return (2*hestonModel_->getMeanReversionSpeed()/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()));
}

double VarianceSwapsHestonAnalyticalPricer::
cTerm(size_t i, const VarianceSwap& varianceSwap) const {
    std::vector<double> dates = varianceSwap.getDates();
    return 2 * hestonModel_->getMeanReversionSpeed() / (hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()*(1.-exp(hestonModel_->getMeanReversionSpeed()*dates[i-1])));
}

double VarianceSwapsHestonAnalyticalPricer::
wTerm(size_t i, const VarianceSwap &varianceSwap) const {
    std::vector<double> dates = varianceSwap.getDates();
    return cTerm(i, varianceSwap)*exp(-hestonModel_->getMeanReversionSpeed()*dates[i-1])*hestonModel_->getInitialVolatility();
}

//log²(Sti/Sti-1)
double VarianceSwapsHestonAnalyticalPricer::
u_iTerm(size_t i, const VarianceSwap &varianceSwap) const {
    std::vector<double> dates = varianceSwap.getDates();
    double first_term = -(( qtilde() + 2*wTerm(i,varianceSwap))+qtilde() + wTerm(i,varianceSwap)*wTerm(i,varianceSwap)) * functionDPrime(dates[i]-dates[i-1],0.).real() / (cTerm(i, varianceSwap)*cTerm(i, varianceSwap));
    double second_term = -(qtilde()+wTerm(i,varianceSwap))*(2*functionCPrime(dates[i]-dates[i-1],0.).real() * functionDPrime(dates[i]-dates[i-1],0.).real() * functionDSecond(dates[i]-dates[i-1],0.).real()) / cTerm(i,varianceSwap);
    double third_term= -(functionCPrime(dates[i]-dates[i-1],0.).real() * functionCPrime(dates[i]-dates[i-1],0.).real() + functionCSecond(dates[i]-dates[i-1],0.).real());
    return exp(hestonModel_->getDrift()*(dates[i]-dates[i-1])) * (first_term + second_term + third_term) ;
}

double VarianceSwapsHestonAnalyticalPricer::price(const VarianceSwap& varianceSwap) const{
    double price=0;
    std::vector<double> dates = varianceSwap.getDates();
    for (size_t i = 1; i < dates.size()+1; i++){
        price = price + u_iTerm(i,varianceSwap);
    }
    price = price * 10000. / (dates.size()*(dates[1]-dates[0]));
    return price;
}

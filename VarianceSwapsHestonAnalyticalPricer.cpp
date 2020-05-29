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
    //return (aTerm(omega)-bTerm(omega))/(aTerm(omega)+bTerm(omega)); //defintion during class
    return (aTerm(omega)+bTerm(omega))/(aTerm(omega)-bTerm(omega)); //definition in PDF
}

//functions C and D and their derivatives
std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionC(double tau, double omega) const {
    std::complex<double> j(0.,1.);
    //return tau * hestonModel_->getDrift() * (j * omega - 1.)+ hestonModel_->getMeanReversionSpeed() *  hestonModel_->getMeanReversionLevel() / (hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()) * ((aTerm(omega)-bTerm(omega))*tau-2.*log((1.-gTerm(omega)*exp(-bTerm(omega)*tau))/(1.-gTerm(omega)))) ; //definition in class
    return tau * hestonModel_->getDrift() * (j * omega - 1.)+ hestonModel_->getMeanReversionSpeed() *  hestonModel_->getMeanReversionLevel() / (hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()) * ((aTerm(omega)+bTerm(omega))*tau-2.*log((1.-gTerm(omega)*exp(bTerm(omega)*tau))/(1.-gTerm(omega)))) ; // definition in PDF
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionD(double tau, double omega) const {
    std::complex<double> D_init;
    //D_init=(aTerm(omega)-bTerm(omega))/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()); // definition in class
    D_init=(aTerm(omega)+bTerm(omega))/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()); //definition in PDF
    //return D_init*(aTerm(omega)+bTerm(omega))*(1.-exp(-bTerm(omega)*tau))/(1.-gTerm(omega)*exp(-bTerm(omega)*tau)); // in class
    return D_init*(aTerm(omega)+bTerm(omega))*(1.-exp(bTerm(omega)*tau))/(1.-gTerm(omega)*exp(bTerm(omega)*tau)); // in PDF
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
    return 2 * hestonModel_->getMeanReversionSpeed() / (hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()*(1.-exp(-hestonModel_->getMeanReversionSpeed()*dates[i-1])));
}

double VarianceSwapsHestonAnalyticalPricer::
wTerm(size_t i, const VarianceSwap &varianceSwap) const {
    std::vector<double> dates = varianceSwap.getDates();
    return cTerm(i, varianceSwap)*exp(-hestonModel_->getMeanReversionSpeed()*dates[i-1])*hestonModel_->getInitialVolatility();
}

//log²(Sti/Sti-1)
double VarianceSwapsHestonAnalyticalPricer::
u_1Term(const VarianceSwap &varianceSwap) const {
    std::vector<double> dates = varianceSwap.getDates();
    double first_term = (functionDPrime(dates[1],0.).real() * hestonModel_->getVolOfVol())*(functionDPrime(dates[1],0.).real() * hestonModel_->getVolOfVol());
    double second_term = hestonModel_->getVolOfVol()*(2*functionCPrime(dates[1],0.).real() * functionDPrime(dates[1],0.).real() - functionDSecond(dates[1],0.).real());
    double third_term = functionCPrime(dates[1],0.).real() * functionCPrime(dates[1],0.).real() - functionCSecond(dates[1],0.).real();
    return (first_term + second_term + third_term) ;
}

double VarianceSwapsHestonAnalyticalPricer::
u_iTerm(size_t i, const VarianceSwap &varianceSwap) const {
    std::vector<double> dates = varianceSwap.getDates();
    double first_term = (( qtilde() + 2*wTerm(i,varianceSwap))+(qtilde() + wTerm(i,varianceSwap))*(qtilde() + wTerm(i,varianceSwap))) * functionDPrime(dates[i]-dates[i-1],0.).real() * functionDPrime(dates[i]-dates[i-1],0.).real() / (cTerm(i, varianceSwap)*cTerm(i, varianceSwap));
    double second_term = (qtilde()+wTerm(i,varianceSwap))*(2*functionCPrime(dates[i]-dates[i-1],0.).real() * functionDPrime(dates[i]-dates[i-1],0.).real() - functionDSecond(dates[i]-dates[i-1],0.).real()) / cTerm(i,varianceSwap);
    double third_term= (functionCPrime(dates[i]-dates[i-1],0.).real() * functionCPrime(dates[i]-dates[i-1],0.).real() - functionCSecond(dates[i]-dates[i-1],0.).real());
    return (first_term + second_term + third_term) ;
}

double VarianceSwapsHestonAnalyticalPricer::price(const VarianceSwap& varianceSwap) const{
    double price=0;
    price = price + u_1Term(varianceSwap);
    std::vector<double> dates = varianceSwap.getDates();
    for (std::size_t i = 2; i < dates.size()+1; i++){
        price = price + u_iTerm(i,varianceSwap);
    }
    price = price * 10000. / (dates.size()*(dates[1]-dates[0]));
    return price;
}

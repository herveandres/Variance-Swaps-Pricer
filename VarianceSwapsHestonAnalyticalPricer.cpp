#include "VarianceSwapsHestonAnalyticalPricer.h"
#include <complex>

VarianceSwapsHestonAnalyticalPricer::VarianceSwapsHestonAnalyticalPricer(
                                            const HestonModel& hestonModel):
                                            VarianceSwapsHestonPricer(hestonModel)
{
    
}

VarianceSwapsHestonAnalyticalPricer::~VarianceSwapsHestonAnalyticalPricer()
{

}

complex<double> VarianceSwapsHestonAnalyticalPricer::
aTerm(const double omega) const {
    complex<double> j(0.,1.);
    return hestonModel_->getMeanReversionSpeed() - hestonModel_->getCorrelation() * hestonModel_->getVolOfVol() * omega * j  ;
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
bTerm(const double omega) const {
    complex<double> j(0.,1.);
    return sqrt(aTerm(omega)+ hestonModel_->getVolOfVol() * hestonModel_->getVolOfVol() * (j*omega + omega*omega) )  ;
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
gTerm(const double omega) const {
    return (aTerm(omega)-bTerm(omega))/(aTerm(omega)+bTerm(omega));
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
FunctionC(const double tau, const double omega) const {
    complex<double> j(0.,1.);
    return tau * hestonModel_->getDrift() * (j * omega - 1.)+ hestonModel_->getMeanReversionSpeed() *  hestonModel_->getMeanReversionLevel() / (hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()) * ((aTerm(omega)-bTerm(omega))*tau-2.*log((1.-gTerm(omega)*exp(-bTerm(omega)*tau))/(1.-gTerm(omega)))) ;
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
FunctionD(const double tau, const double omega) const {
    complex<double> D_init;
    D_init=(aTerm(omega)-bTerm(omega))/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol());
    return D_init*(aTerm(omega)+bTerm(omega))*(1.-exp(-bTerm(omega)*tau))/(1.-gTerm(omega)*exp(-bTerm(omega)*tau));
}

double VarianceSwapsHestonAnalyticalPricer::price(const VarianceSwap& varianceSwap) const{
    //A compléter
    return 0;
}

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
a_term(const double omega) const {
    complex<double> j(0.,1.);
    return hestonModel_->getMeanReversionSpeed() - hestonModel_->getCorrelation() * hestonModel_->getVolOfVol() * omega * j  ;
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
b_term(const double omega) const {
    complex<double> j(0.,1.);
    return sqrt(a_term(omega)+ hestonModel_->getVolOfVol() * hestonModel_->getVolOfVol() * (j*omega + omega*omega) )  ;
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
g_term(const double omega) const {
    return (a_term(omega)-b_term(omega))/(a_term(omega)+b_term(omega));
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
Function_C(const double tau, const double omega) const {
    complex<double> j(0.,1.);
    return tau * hestonModel_->getDrift() * (j * omega - 1.)+ hestonModel_->getMeanReversionSpeed() *  hestonModel_->getMeanReversionLevel() / (hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()) * ((a_term(omega)-b_term(omega))*tau-2.*log((1.-g_term(omega)*exp(-b_term(omega)*tau))/(1.-g_term(omega)))) ;
}

complex<double> VarianceSwapsHestonAnalyticalPricer::
Function_D(const double tau, const double omega) const {
    complex<double> D_init;
    D_init=(a_term(omega)-b_term(omega))/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol());
    return D_init*(a_term(omega)+b_term(omega))*(1.-exp(-b_term(omega)*tau))/(1.-g_term(omega)*exp(-b_term(omega)*tau));
}

double VarianceSwapsHestonAnalyticalPricer::price(const VarianceSwap& varianceSwap) const{
    //A compléter
    return 0;
}

#ifndef VARIANCESWAPSHESTONANALYTICALPRICER_H
#define VARIANCESWAPSHESTONANALYTICALPRICER_H

#include "VarianceSwapsPricer.h"
#include <complex>

class VarianceSwapsHestonAnalyticalPricer : public VarianceSwapsHestonPricer
{
private:
    HestonModel* hestonModel_;

    //useful variables to compute function C and D
    std::complex<double> aTerm (double omega) const ;
    std::complex<double> bTerm (double omega) const ;
    std::complex<double> gTerm (double omega) const ;

    //functions C and D and their derivatives
    std::complex<double> functionC (double tau, double omega) const;
    std::complex<double> functionD (double tau, double omega) const;
    std::complex<double> functionDPrime (double tau, double omega) const;
    std::complex<double> functionCPrime (double tau, double omega) const;
    std::complex<double> functionDSecond (double tau, double omega) const;
    std::complex<double> functionCSecond (double tau, double omega) const;

    //useful terms with respect to the Chi2 law
    double qtildeTerm () const;
    double cTerm (double t) const;
    double wTerm (double t) const;

    //E[log²(Sti/Sti-1)]
    double u1Term (double t1, double t2) const; //Case i = 1
    double uiTerm (double t1, double t2) const; //Case i > 1
public:
    VarianceSwapsHestonAnalyticalPricer(const HestonModel& hestonModel);

    // Copy constructor, Assignement operator and Destructor are needed because one of the member variable is a pointer
    ~VarianceSwapsHestonAnalyticalPricer();
    VarianceSwapsHestonAnalyticalPricer(const VarianceSwapsHestonAnalyticalPricer& analyticalPricer);
    VarianceSwapsHestonAnalyticalPricer& operator=(
                        const VarianceSwapsHestonAnalyticalPricer& analyticalPricer);

    //Method returning the analytical price of the variance swap given as argument
    double price(const VarianceSwap& varianceSwap) const override;

    //Method returning the analytical price in the continuous case
    double continousPrice(const VarianceSwap& varianceSwap);
};

#endif

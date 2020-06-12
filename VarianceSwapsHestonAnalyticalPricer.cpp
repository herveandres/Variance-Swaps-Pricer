#include <complex>
#include "VarianceSwapsHestonAnalyticalPricer.h"
#include "MathFunctions.h"

VarianceSwapsHestonAnalyticalPricer::VarianceSwapsHestonAnalyticalPricer(
                                            const HestonModel& hestonModel):
        hestonModel_(new HestonModel(hestonModel))
{
    
}

VarianceSwapsHestonAnalyticalPricer::~VarianceSwapsHestonAnalyticalPricer()
{
    delete hestonModel_;
}

std::complex<double>j(0.,1.);
//useful variables to compute function C and D
std::complex<double> VarianceSwapsHestonAnalyticalPricer::
aTerm(double omega) const {
    return hestonModel_->getMeanReversionSpeed() - hestonModel_->getCorrelation() * hestonModel_->getVolOfVol() * omega * j  ;
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
bTerm(double omega) const {
    double sigma = hestonModel_->getVolOfVol();
    std::complex<double> a = aTerm(omega);
    return std::sqrt(a*a+ sigma*sigma * (j*omega + omega*omega) )  ;
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
gTerm(double omega) const {
    std::complex<double> a = aTerm(omega),
                    b = bTerm(omega);
    // std::cout << a << " " << b << std::endl;
    // std::cout << a-b << std::endl;
    return (a-b)/(a+b); //defintion during class
    // return (a+b)/(a-b); //definition in PDF
}

//functions C and D and their derivatives
std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionC(double tau, double omega) const {
    double r = hestonModel_->getDrift(),
           kappa = hestonModel_->getMeanReversionSpeed(),
           theta = hestonModel_->getMeanReversionLevel(),
           sigma = hestonModel_->getVolOfVol();
    std::complex<double> a = aTerm(omega),
                        b = bTerm(omega),
                        g = gTerm(omega);
    // std::cout << a << " " << b << " " << g << std::endl;
    return tau * r * (j * omega - 1.)+ kappa * theta / (sigma*sigma) * ((a-b)*tau-2.*std::log((1.-g*std::exp(-b*tau))/(1.-g))) ; // definition in class
    // return tau * r * (j * omega - 1.)+ kappa * theta / (sigma*sigma) * ((a+b)*tau-2.*std::log((1.-g*std::exp(b*tau))/(1.-g))) ; // definition in PDF
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionD(double tau, double omega) const {
    double sigma = hestonModel_->getVolOfVol();
    std::complex<double> a = aTerm(omega), 
                         b = bTerm(omega),
                         g = gTerm(omega);
    std::complex<double> D_init = (a-b)/(sigma*sigma); //definition during class
    // std::complex<double> D_init = (a+b)/(sigma*sigma); //definition in PDF
    //D_init=(aTerm(omega)-bTerm(omega))/(hestonModel_->getVolOfVol()*hestonModel_->getVolOfVol()); // definition in class
    return D_init*(1.-std::exp(-b*tau))/(1.-g*std::exp(-b*tau)); // in class
    // return D_init*(1.-std::exp(b*tau))/(1.-g*std::exp(b*tau)); // in PDF
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionDPrime(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionD(tau,omega);
    };
    return MathFunctions::finiteDifference(f,omega,tau);
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionCPrime(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionC(tau,omega);
    };
    return MathFunctions::finiteDifference(f,omega,tau);
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionCSecond(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionCPrime(tau,omega);
    };
    return MathFunctions::finiteDifference(f,omega,tau);
}

std::complex<double> VarianceSwapsHestonAnalyticalPricer::
functionDSecond(double tau, double omega) const {
    std::function<std::complex<double>(double,double)> f = [=](double tau, double omega){
        return this->functionDPrime(tau,omega);
    };
    return MathFunctions::finiteDifference(f,omega,tau);
}

//useful terms in respect of the Chi2 law
double VarianceSwapsHestonAnalyticalPricer::
qtildeTerm() const {
    double kappa = hestonModel_->getMeanReversionSpeed(), 
           theta = hestonModel_->getMeanReversionLevel(),
           sigma = hestonModel_->getVolOfVol();
    return (2*kappa*theta/(sigma*sigma));
}

double VarianceSwapsHestonAnalyticalPricer::
cTerm(double t) const {
    double kappa = hestonModel_->getMeanReversionSpeed(),
           sigma = hestonModel_->getVolOfVol();
    return 2 * kappa / (sigma*sigma*(1.-std::exp(-kappa*t)));
}

double VarianceSwapsHestonAnalyticalPricer::
wTerm(double t) const {
    return cTerm(t)*std::exp(-hestonModel_->getMeanReversionSpeed()*t)*hestonModel_->getInitialVolatility();
}

//log²(Sti/Sti-1)
double VarianceSwapsHestonAnalyticalPricer::
u1Term(double t1, double t2) const {
    double delta = t2-t1,
           v0 = hestonModel_->getInitialVolatility();
    std::complex<double> C1 = functionCPrime(delta,0.),
                        C2 = functionCSecond(delta,0.),
                        D1 = functionDPrime(delta,0.),
                        D2 = functionDSecond(delta,0.);
    std::complex<double> firstTerm = D1*D1*v0*v0,
                        secondTerm = v0*(2.0*C1*D1-D2),
                        thirdTerm = C1*C1 - C2;
    std::complex<double> u1 = firstTerm+secondTerm+thirdTerm; 
    return u1.real();
}

double VarianceSwapsHestonAnalyticalPricer::
uiTerm(double t1, double t2) const {
    double delta = t2-t1,
           qtilde = qtildeTerm(),
           Wi = wTerm(t1),
           ci = cTerm(t1);
    std::complex<double> C1 = functionCPrime(delta,0.),
           C2 = functionCSecond(delta,0.),
           D1 = functionDPrime(delta,0.),
           D2 = functionDSecond(delta,0.);
    std::complex<double> firstTerm = (qtilde+2*Wi+(qtilde + Wi)*(qtilde + Wi)) 
                                        *D1*D1/(ci*ci),
                        secondTerm = (qtilde+Wi)*(2.0*C1*D1 - D2) / ci,
                        thirdTerm= C1*C1 - C2;
    std::complex<double> ui = firstTerm+secondTerm+thirdTerm;
    return ui.real();
}

double VarianceSwapsHestonAnalyticalPricer::price(const VarianceSwap& varianceSwap) const{
    double price=0;
    std::vector<double> dates = varianceSwap.getDates();

    price = price + u1Term(dates[0],dates[1]);
    for (std::size_t i = 2; i < dates.size(); i++)
    {   
        // std::cout << uiTerm(dates[i-1],dates[i]) << std::endl;
        price = price + uiTerm(dates[i-1],dates[i]);
    }
    price = price * 10000. / dates.back();
    return price;
}

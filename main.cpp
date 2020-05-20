#include <iostream>
#include "HestonLogSpotPathSimulator.h"
#include "HestonVariancePathSimulator.h"
#include "VarianceSwap.h"
#include "MathFunctions.h"
#include "VarianceSwapsHestonMonteCarloPricer.h"

int main()
{
    //Heston model parameters
    double drift = 0, kappa = 0.5, theta = 0.04, eps = 1.0, rho = -0.9,
            V0 = 0.04, X0 = 100;

    HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

    //Variance swap parameters
    double maturity = 10.0;
    std::size_t nbOfObservations = maturity*2+1;

    VarianceSwap varianceSwap(maturity,nbOfObservations);
    std::vector<double> dates = varianceSwap.getDates();

    //Simulation parameters
    size_t nbSimulations = 100000, nbTimePoints = 100;
    std::vector<double> timePoints, temp;
    for(std::size_t i = 0; i < dates.size()-1; i++)
    {   
        temp = MathFunctions::buildLinearSpace(dates[i],dates[i+1],nbTimePoints);
        timePoints.insert(timePoints.end(), temp.begin(), temp.end()-1);
        if(i == dates.size()-2)
            timePoints.push_back(temp.back());
    }
    std::cout << "Computation of the price using TG + BroadieKaya" << std::endl;
    TruncatedGaussianScheme truncatedGaussianScheme(timePoints,hestonModel);
    BroadieKayaScheme broadieKayaSchemeTG(timePoints,hestonModel,truncatedGaussianScheme);
    VarianceSwapsHestonMonteCarloPricer mcPricerBKTG(hestonModel,broadieKayaSchemeTG,nbSimulations);
    std::cout << mcPricerBKTG.price(varianceSwap) << std::endl << std::endl;

    std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;
    QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModel);
    BroadieKayaScheme broadieKayaSchemeQE(timePoints,hestonModel,quadraticExponentialScheme);
    VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(hestonModel,broadieKayaSchemeQE,nbSimulations);
    std::cout << mcPricerBKQE.price(varianceSwap) << std::endl;
    return 0;
}
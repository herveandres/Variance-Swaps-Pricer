#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include "HestonLogSpotPathSimulator.h"
#include "HestonVariancePathSimulator.h"
#include "VarianceSwap.h"
#include "MathFunctions.h"
#include "VarianceSwapsHestonMonteCarloPricer.h"
#include "VarianceSwapsHestonAnalyticalPricer.h"

void testNbOfObservations()
{
    //Heston model parameters
    double drift = 0, kappa = 0.5, theta = 0.04, eps = 1, rho = -0.9,
            V0 = 0.04, X0 = 100;

    HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

    //Variance swap parameters
    double maturity = 10.0;
    double nbOfObservations;

    std::ofstream file;
    file.open ("../Tests/test_convergence_nb_of_observations_analytical.csv");
    file << "Nombre d'observations;Prix Analytique \n";

    for (size_t i=2 ; i<31 ; i=i+2)
    {
        nbOfObservations= i*maturity+1;

        VarianceSwap varianceSwap(maturity,nbOfObservations);

        std::cout << "Analytical computation of the price for " << nbOfObservations << " observations" << std::endl;
        VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
        double analyticalPrice = anPricer.price(varianceSwap);
        std::cout << analyticalPrice << std::endl << std::endl;

        file << nbOfObservations << ";";
        file << analyticalPrice << "\n";

    }
    file.close();
}

void testNbOfSimulations()
{
    //Heston model parameters
    double drift = 0, kappa = 0.5, theta = 0.04, eps = 1, rho = -0.9,
            V0 = 0.04, X0 = 100;

    HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

    //Variance swap parameters
    double maturity = 10.0;
    double nbOfObservations = 2*maturity+1;

    VarianceSwap varianceSwap(maturity,nbOfObservations);

    std::cout << "Analytical computation of the price" << std::endl;
    VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
    double analyticalPrice = anPricer.price(varianceSwap);
    std::cout << analyticalPrice << std::endl << std::endl;

    size_t nbSimulationsMin=10000;
    size_t nbSimulationsMax=210000;
    size_t pasSimulations=40000;
    size_t nbTimePoints = 200;
    std::vector<double> dates = varianceSwap.getDates();

    std::vector<double> timePoints, temp;
    for(std::size_t j = 0; j < dates.size()-1; j++)
    {
        temp = MathFunctions::buildLinearSpace(dates[j],dates[j+1],nbTimePoints);
        timePoints.insert(timePoints.end(), temp.begin(), temp.end()-1);
        if(j == dates.size()-2)
            timePoints.push_back(temp.back());
    }

    TruncatedGaussianScheme truncatedGaussianScheme(timePoints,hestonModel);
    BroadieKayaScheme broadieKayaSchemeTG(timePoints,hestonModel,truncatedGaussianScheme);

    QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModel);
    BroadieKayaScheme broadieKayaSchemeQE(timePoints,hestonModel,quadraticExponentialScheme);

    std::ofstream file;
    file.open ("../Tests/test_convergence_nb_of_simulations.csv");
    file << "Nombre de simulations;Prix BKTG;Prix BKQE \n";

    for (size_t nbSimulations = nbSimulationsMin; nbSimulations<nbSimulationsMax+1 ; nbSimulations=nbSimulations+pasSimulations ){

        std::cout << "---------- Nombre de simulations : " << nbSimulations << " ----------" << std::endl << std::endl;
        std::cout << "Computation of the price using TG + BroadieKaya" << std::endl;

        VarianceSwapsHestonMonteCarloPricer mcPricerBKTG(hestonModel,broadieKayaSchemeTG,nbSimulations);
        double BKTGprice = mcPricerBKTG.price(varianceSwap);
        std::cout << BKTGprice << std::endl << std::endl;

        file << nbSimulations << ";";
        file << BKTGprice << ";";

        std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;

        VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(hestonModel,broadieKayaSchemeQE,nbSimulations);
        double BKQEprice = mcPricerBKQE.price(varianceSwap);
        std::cout << BKQEprice << std::endl << std::endl << std::endl;

        file << BKQEprice << "\n";
    }
    file.close();

}

void testThreeParametersSets()
{
    std::vector<std::map<std::string,double> > parametersSets;

    //Heston model parameters and variance 
    //Case I
    std::map<std::string,double> parameters;
    parameters["drift"] = 0; parameters["kappa"] = 0.5; parameters["theta"] = 0.04;
    parameters["eps"] = 1; parameters["rho"] = -0.9; parameters["V0"] = 0.04;
    parameters["X0"] = 100; parameters["maturity"] = 10.0;
    parametersSets.push_back(parameters);

    //Case II
    parameters["drift"] = 0; parameters["kappa"] = 0.3; parameters["theta"] = 0.04;
    parameters["eps"] = 0.9; parameters["rho"] = -0.5; parameters["V0"] = 0.04;
    parameters["X0"] = 100; parameters["maturity"] = 15.0;
    parametersSets.push_back(parameters);

    //Case III
    parameters["drift"] = 0; parameters["kappa"] = 1; parameters["theta"] = 0.09;
    parameters["eps"] = 1; parameters["rho"] = -0.3; parameters["V0"] = 0.09;
    parameters["X0"] = 100; parameters["maturity"] = 5.0;
    parametersSets.push_back(parameters);


    //Creation of Heston models and variance swaps.
    std::vector<HestonModel> hestonModels;
    std::vector<VarianceSwap> varianceSwaps;

    size_t nbOfObservations;
    for(size_t i = 0; i < parametersSets.size(); i++)
    {
        hestonModels.push_back(HestonModel(parametersSets[i]["drift"],
                                        parametersSets[i]["kappa"],
                                        parametersSets[i]["theta"],
                                        parametersSets[i]["eps"],
                                        parametersSets[i]["rho"],
                                        parametersSets[i]["V0"],
                                        parametersSets[i]["X0"]));
        nbOfObservations = 2*parametersSets[i]["maturity"]+1;
        varianceSwaps.push_back(VarianceSwap(parametersSets[i]["maturity"],
                                             nbOfObservations));
        
    }

    //Computations of prices using the three methods.
    std::vector<double> dates;
    size_t nbSimulations = 10000, nbTimePoints = 200;
    for(size_t i = 0; i < hestonModels.size(); i++)
    {
        std::vector<double> timePoints, temp;
        dates = varianceSwaps[i].getDates();
        for(std::size_t j = 0; j < dates.size()-1; j++)
        {   
            temp = MathFunctions::buildLinearSpace(dates[j],dates[j+1],nbTimePoints);
            timePoints.insert(timePoints.end(), temp.begin(), temp.end()-1);
            if(j == dates.size()-2)
                timePoints.push_back(temp.back());
        }

        std::cout << "------------- Case " << i+1 << " -------------" << std::endl;
        std::cout << "Analytical computation of the price" << std::endl;
        VarianceSwapsHestonAnalyticalPricer anPricer(hestonModels[i]);
        std::cout << anPricer.price(varianceSwaps[i]) << std::endl << std::endl;


        std::cout << "Computation of the price using TG + BroadieKaya" << std::endl;
        TruncatedGaussianScheme truncatedGaussianScheme(timePoints,hestonModels[i]);
        BroadieKayaScheme broadieKayaSchemeTG(timePoints,hestonModels[i],truncatedGaussianScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKTG(hestonModels[i],broadieKayaSchemeTG,nbSimulations);
        std::cout << mcPricerBKTG.price(varianceSwaps[i]) << std::endl << std::endl;

        std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;
        QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModels[i]);
        BroadieKayaScheme broadieKayaSchemeQE(timePoints,hestonModels[i],quadraticExponentialScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(hestonModels[i],broadieKayaSchemeQE,nbSimulations);
        std::cout << mcPricerBKQE.price(varianceSwaps[i]) << std::endl;
        std::cout << std::endl << std::endl;
    }
}

void testEvolutionMCPricesWithDiscretizationTimestep()
{
    //Heston model parameters
    double drift = 0, kappa = 0.5, theta = 0.04, eps = 1, rho = -0.9,
            V0 = 0.04, X0 = 100;

    HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

    //Variance swap parameters
    double maturity = 10.0;
    std::size_t nbOfObservations = maturity*2+1;

    VarianceSwap varianceSwap(maturity,nbOfObservations);


    std::ofstream file;
    file.open ("../Tests/test_convergence_timestep_wider.csv");
    //Pricing analytique
    std::cout << "Analytical computation of the price" << std::endl;
    VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
    double analyticalPrice = anPricer.price(varianceSwap);
    std::cout << analyticalPrice << std::endl << std::endl;

    //Pricing Monte-Carlo
    file << "Nombre de points;Prix Analytique;Prix BKTG;Prix BKQE \n";
    std::vector<double> dates = varianceSwap.getDates();
    size_t nbSimulations = 10000;
    std::vector<double> nbTimePoints{100,500,1000,1500,2000,2500,3000};

    for(std::size_t i = 0; i < nbTimePoints.size(); i++)
    {
        file << nbTimePoints[i] << ";"; 
        file << analyticalPrice << ";";
        std::vector<double> timePoints, temp;  
        for(std::size_t j = 0; j < dates.size()-1; j++)
        {   
            temp = MathFunctions::buildLinearSpace(dates[j],dates[j+1],nbTimePoints[i]);
            timePoints.insert(timePoints.end(), temp.begin(), temp.end()-1);
            if(j == dates.size()-2)
                timePoints.push_back(temp.back());
        }
        std::cout << "---------- Nombre de points : " << nbTimePoints[i] << " ----------" << std::endl << std::endl;
        std::cout << "Computation of the price using TG + BroadieKaya" << std::endl;
        TruncatedGaussianScheme truncatedGaussianScheme(timePoints,hestonModel);
        BroadieKayaScheme broadieKayaSchemeTG(timePoints,hestonModel,truncatedGaussianScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKTG(hestonModel,broadieKayaSchemeTG,nbSimulations);
        double BKTGprice = mcPricerBKTG.price(varianceSwap);
        std::cout << BKTGprice << std::endl << std::endl;
        file << BKTGprice << ";";

        std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;
        QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModel);
        BroadieKayaScheme broadieKayaSchemeQE(timePoints,hestonModel,quadraticExponentialScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(hestonModel,broadieKayaSchemeQE,nbSimulations);
        double BKQEprice = mcPricerBKQE.price(varianceSwap);
        std::cout << BKQEprice << std::endl << std::endl << std::endl;
        file << BKQEprice << "\n";
    }

    file.close();
}

int main()
{   
    // testThreeParametersSets();
    //testEvolutionMCPricesWithDiscretizationTimestep();
    //testNbOfObservations();
    testNbOfSimulations();
    return 0;
}

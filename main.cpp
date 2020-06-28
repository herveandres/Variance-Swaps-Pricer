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

void testKappaParameter(){
    //Heston model parameters
    double drift = 0, theta = 0.04, eps = 1, rho = -0.9,
            V0 = 0.04, X0 = 100;
    double kappaInit = 0.2;

    //Variance swap parameters
    double maturity = 1.0;
    double nbOfObservations = 2*maturity+1;

    VarianceSwap varianceSwap(maturity,nbOfObservations);

    std::ofstream file;
    file.open ("C:/Users/Arnaud/Documents/GitHub/Variance-Swaps-Pricer/Tests/test_kappa_influence.csv");
    file << "Kappa;Prix Analytique \n";

    for (size_t i=0 ; i<15 ; i=i+1){
        double kappa = kappaInit + i * 0.05;
        HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

        std::cout << "Analytical computation of the price for kappa = " << kappa << std::endl;
        VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
        double analyticalPrice = anPricer.price(varianceSwap);
        std::cout << analyticalPrice << std::endl << std::endl;

        file << kappa << ";";
        file << analyticalPrice << "\n";

    }
    file.close();
}

void testMaturityParameter(){
    //Heston model parameters
    double drift = 0, theta = 0.04, eps = 1, rho = -0.9,
            V0 = 0.04, X0 = 100, kappa = 0.5 ;

    HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

    //Variance swap parameters
    double maturityInit = 0.5;


    std::ofstream file;
    file.open ("C:/Users/Arnaud/Documents/GitHub/Variance-Swaps-Pricer/Tests/test_maturity_influence.csv");
    file << "Kappa;Prix Analytique \n";

    for (size_t i=0 ; i<20 ; i=i+1){
        double maturity = maturityInit + i * 0.5;
        double nbOfObservations = 2*maturity+1;

        VarianceSwap varianceSwap(maturity,nbOfObservations);

        std::cout << "Analytical computation of the price for maturity = " << maturity << std::endl;
        VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
        double analyticalPrice = anPricer.price(varianceSwap);
        std::cout << analyticalPrice << std::endl << std::endl;

        file << maturity << ";";
        file << analyticalPrice << "\n";

    }
    file.close();
}

void testNbOfObservations()
{
    //Heston model parameters
    double drift = 0, kappa = 0.5, theta = 0.04, eps = 1, rho = -0.9,
            V0 = 0.04, X0 = 100;

    HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

    //Variance swap parameters
    double maturity = 10.0;
    double nbOfObservations;
    double nbOfObservationsMax= 60;

    std::ofstream file;
    file.open ("../Tests/test_convergence_nb_of_observations_analytical.csv");
    file << "Nombre d'observations;Prix Analytique;Prix Analytique Continu; Différence \n";

    for (size_t i=2 ; i<nbOfObservationsMax+1 ; i=i+20)
    {
        nbOfObservations= i*maturity+1;

        VarianceSwap varianceSwap(maturity,nbOfObservations);

        std::cout << "Analytical computation of the price for " << nbOfObservations << " observations" << std::endl;
        VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
        double analyticalPrice = anPricer.price(varianceSwap);
        std::cout << analyticalPrice << std::endl << std::endl;

        std::cout << "Analytical computation of the continous price " << std::endl;
        double continuousAnalyticalPrice = anPricer.continousPrice(varianceSwap);
        std::cout << continuousAnalyticalPrice << std::endl << std::endl;

        std::cout << "Difference" << std::endl;
        std::cout << analyticalPrice-continuousAnalyticalPrice << std::endl << std::endl;

        file << nbOfObservations << ";";
        file << analyticalPrice << ";";
        file << continuousAnalyticalPrice << ";";
        file << analyticalPrice-continuousAnalyticalPrice << "\n";

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
    double maturity = 1.0;
    double nbOfObservations = 2*maturity+1;

    VarianceSwap varianceSwap(maturity,nbOfObservations);

    std::cout << "Analytical computation of the price" << std::endl;
    VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
    double analyticalPrice = anPricer.price(varianceSwap);
    std::cout << analyticalPrice << std::endl << std::endl;

    size_t nbSimulationsMin=100000;
    size_t nbSimulationsMax=1000000;
    size_t pasSimulations=100000;
    size_t nbTimePoints = 1000;
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
    BroadieKayaScheme broadieKayaSchemeTG(truncatedGaussianScheme);

    QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModel);
    BroadieKayaScheme broadieKayaSchemeQE(quadraticExponentialScheme);

    std::ofstream file;
    file.open ("../Tests/test_convergence_nb_of_simulations.csv");
    file << "Nombre de simulations;Prix analytique;Prix BKTG;Prix BKQE \n";

    for (size_t nbSimulations = nbSimulationsMin; nbSimulations<nbSimulationsMax+1 ; nbSimulations=nbSimulations+pasSimulations ){

        std::cout << "---------- Nombre de simulations : " << nbSimulations << " ----------" << std::endl << std::endl;
        std::cout << "Computation of the price using TG + BroadieKaya" << std::endl;

        VarianceSwapsHestonMonteCarloPricer mcPricerBKTG(broadieKayaSchemeTG,nbSimulations);
        double BKTGprice = mcPricerBKTG.price(varianceSwap);
        std::cout << BKTGprice << std::endl << std::endl;

        file << nbSimulations << ";";
        file << analyticalPrice << ";";
        file << BKTGprice << ";";

        std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;

        VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(broadieKayaSchemeQE,nbSimulations);
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
        BroadieKayaScheme broadieKayaSchemeTG(truncatedGaussianScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKTG(broadieKayaSchemeTG,nbSimulations);
        std::cout << mcPricerBKTG.price(varianceSwaps[i]) << std::endl << std::endl;

        std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;
        QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModels[i]);
        BroadieKayaScheme broadieKayaSchemeQE(quadraticExponentialScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(broadieKayaSchemeQE,nbSimulations);
        std::cout << mcPricerBKQE.price(varianceSwaps[i]) << std::endl;
        std::cout << std::endl << std::endl;
    }
}

void testDiscretizationTimestep()
{
    //Heston model parameters
    double drift = 0, kappa = 0.5, theta = 0.04, eps = 1, rho = -0.9,
            V0 = 0.04, X0 = 100;

    HestonModel hestonModel(drift,kappa,theta,eps,rho,V0,X0);

    //Variance swap parameters
    double maturity = 1.0;
    std::size_t nbOfObservations = maturity*2+1;

    VarianceSwap varianceSwap(maturity,nbOfObservations);


    std::ofstream file;
    file.open ("../Tests/test_convergence_timestep_temp.csv");
    //Pricing analytique
    std::cout << "Analytical computation of the price" << std::endl;
    VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
    double analyticalPrice = anPricer.price(varianceSwap);
    std::cout << analyticalPrice << std::endl << std::endl;

    //Pricing Monte-Carlo
    file << "Nombre de points;Prix Analytique;Prix BKTG;Prix BKQE \n";
    std::vector<double> dates = varianceSwap.getDates();
    size_t nbSimulations = 10000;
    std::vector<double> nbTimePoints{100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};

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
        BroadieKayaScheme broadieKayaSchemeTG(truncatedGaussianScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKTG(broadieKayaSchemeTG,nbSimulations);
        double BKTGprice = mcPricerBKTG.price(varianceSwap);
        std::cout << BKTGprice << std::endl << std::endl;
        file << BKTGprice << ";";

        std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;
        QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModel);
        BroadieKayaScheme broadieKayaSchemeQE(quadraticExponentialScheme);
        VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(broadieKayaSchemeQE,nbSimulations);
        double BKQEprice = mcPricerBKQE.price(varianceSwap);
        std::cout << BKQEprice << std::endl << std::endl << std::endl;
        file << BKQEprice << "\n";
    }

    file.close();
}


void test_QEMC(){
    //Heston model parameters
       double drift = 0, kappa = 0.5, theta = 0.04, eps = 0.5, rho = -0.9,
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
       std::cout << "Analytical computation of the price" << std::endl;
       VarianceSwapsHestonAnalyticalPricer anPricer(hestonModel);
       std::cout << anPricer.price(varianceSwap) << std::endl << std::endl;

    std::cout << "Computation of the price using QE MC + BroadieKaya" << std::endl;
    QuadraticExponentialMartingaleCorrectionScheme quadraticExponentialMCScheme(timePoints,hestonModel, true);
    BroadieKayaScheme broadieKayaSchemeQEMC(quadraticExponentialMCScheme);
    VarianceSwapsHestonMonteCarloPricer mcPricerBKQEMC(broadieKayaSchemeQEMC,nbSimulations);
    std::cout << mcPricerBKQEMC.price(varianceSwap) << std::endl;

    std::cout << "Computation of the price using QE + BroadieKaya" << std::endl;
    QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModel, true);
    BroadieKayaScheme broadieKayaSchemeQE(quadraticExponentialScheme);
    VarianceSwapsHestonMonteCarloPricer mcPricerBKQE(broadieKayaSchemeQE,nbSimulations);
    std::cout << mcPricerBKQE.price(varianceSwap) << std::endl;
}

int main()
{   
    //test_QEMC();
    //testThreeParametersSets();
    //testDiscretizationTimestep();
    //testNbOfObservations();
    //testNbOfSimulations();
    //Heston model parameters
    //testKappaParameter();
    testMaturityParameter();
    return 0;
}

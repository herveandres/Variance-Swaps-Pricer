#include "VarianceSwapsHestonMonteCarloPricer.h"
#include <iostream>

VarianceSwapsHestonMonteCarloPricer::VarianceSwapsHestonMonteCarloPricer
                                (const HestonLogSpotPathSimulator& hestonPathSimulator,
                                std::size_t nbSimulations):
            hestonPathSimulator_(hestonPathSimulator.clone()),
            nbSimulations_(nbSimulations)
{

}

VarianceSwapsHestonMonteCarloPricer::~VarianceSwapsHestonMonteCarloPricer()
{
    delete hestonPathSimulator_;
}

VarianceSwapsHestonMonteCarloPricer::VarianceSwapsHestonMonteCarloPricer(
                    const VarianceSwapsHestonMonteCarloPricer& mcPricer):
        hestonPathSimulator_(mcPricer.hestonPathSimulator_->clone()),
        nbSimulations_(mcPricer.nbSimulations_)
{

}
    
VarianceSwapsHestonMonteCarloPricer& VarianceSwapsHestonMonteCarloPricer::operator=(
                    const VarianceSwapsHestonMonteCarloPricer& mcPricer)
{
    if (this == &mcPricer)
		return *this;
	else
	{
		delete hestonPathSimulator_;												
		hestonPathSimulator_ = mcPricer.hestonPathSimulator_->clone();
	}
	return *this;
}

double VarianceSwapsHestonMonteCarloPricer::pathPrice(std::vector<double> path,
                                                double maturity) const
{
    double pathPrice = 0.0; 
    for(size_t i = 0; i < path.size()-1; i++)
    {   
        pathPrice += std::pow(path[i+1]-path[i],2); //Reminder : path[i] = log(S_ti) 
    }
    return 100*100*pathPrice/maturity;
}

double VarianceSwapsHestonMonteCarloPricer::price(const VarianceSwap& varianceSwap) const
{
    double price = 0.;
    std::vector<double> dates = varianceSwap.getDates();
    std::vector<double> simulationTimeSteps = hestonPathSimulator_->getTimePoints();

    //We look for the indexes of the simulated path corresponding to the dates 
    //of the variance swaps
    std::vector<double> indexes;
    std::size_t nbSimulationsBetweenDates = (simulationTimeSteps.size()+dates.size()-2)/(dates.size()-1);
    for(size_t i = 0; i < dates.size(); i++)
    {
        indexes.push_back((nbSimulationsBetweenDates-1)*i);
    }
    double maturity = dates.back();
    std::vector<double> simulatedPath;
    std::vector<double> pathForPricing;
	for (size_t simulationIdx = 0; simulationIdx < nbSimulations_; ++simulationIdx)
	{
		simulatedPath = hestonPathSimulator_->path();
        //"Slicing" of the path in order to keep only the logspot at the dates
        //of the variance swaps
        for(size_t i = 0; i < indexes.size(); i++)
        {
            pathForPricing.push_back(simulatedPath[indexes[i]]);
        }
		price += pathPrice(pathForPricing, maturity);
        pathForPricing.clear();
	}
	price /= nbSimulations_;
	return price;
}
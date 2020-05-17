#include "VarianceSwapsHestonMonteCarloPricer.h"

VarianceSwapsHestonMonteCarloPricer::VarianceSwapsHestonMonteCarloPricer
                                (const HestonModel& hestonModel,
                                const HestonLogSpotPathSimulator& hestonPathSimulator,
                                std::size_t nbSimulations):
            VarianceSwapsHestonPricer(hestonModel), hestonPathSimulator_(hestonPathSimulator.clone()),
            nbSimulations_(nbSimulations)
{

}

VarianceSwapsHestonMonteCarloPricer::~VarianceSwapsHestonMonteCarloPricer()
{
    delete hestonPathSimulator_;
}

double VarianceSwapsHestonMonteCarloPricer::pathPrice(std::vector<double> path,
                                                double maturity) const
{
    double pathPrice = 0.0; 
    for(size_t i = 0; i < path.size()-1; i++)
    {   
        pathPrice += std::pow(std::log(path[i+1]/path[i]),2); 
    }
    return 100*100*pathPrice/maturity;
}

double VarianceSwapsHestonMonteCarloPricer::price(const VarianceSwap& varianceSwap) const
{
    double price = 0.;
    std::vector<double> dates = varianceSwap.getDates();
    std::vector<double> simulationTimeSteps = hestonPathSimulator_->getTimePoints();
    std::vector<double> indexes{0};

    //We look for the indexes of the simulated path corresponding to the dates 
    //of the variance swaps
    for(size_t i = 1, j = 1; i < simulationTimeSteps.size(); i++)
    {
        if(abs(dates[j]-simulationTimeSteps[i]) <= std::pow(10,-8))
        {
            indexes.push_back(i); 
            j++;
        }
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
	}
	price /= nbSimulations_;
	return price;
}
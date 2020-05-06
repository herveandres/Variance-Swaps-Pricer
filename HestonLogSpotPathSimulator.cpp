#include <cmath>
#include "HestonLogSpotPathSimulator.h"

HestonLogSpotPathSimulator::HestonLogSpotPathSimulator(
                                const std::vector<double>& timePoints,
                                const HestonModel& hestonModel,
                                const HestonVariancePathSimulator& variancePathSimulator):
                                PathSimulator(log(hestonModel.getInitialAssetValue()),
                                              timePoints),
                                hestonModel_(new HestonModel(hestonModel)),
                                variancePathSimulator_(variancePathSimulator.clone())
{

}

HestonLogSpotPathSimulator::~HestonLogSpotPathSimulator()
{
    delete hestonModel_;
    delete variancePathSimulator_;
}

BroadieKayaScheme::BroadieKayaScheme(const std::vector<double>& timePoints,
                    const HestonModel& hestonModel,
                    const HestonVariancePathSimulator& variancePathSimulator):
        HestonLogSpotPathSimulator(timePoints, hestonModel, variancePathSimulator)
{

}

BroadieKayaScheme::BroadieKayaScheme(const BroadieKayaScheme& broadieKayaScheme):
        HestonLogSpotPathSimulator(broadieKayaScheme.timePoints_,
                                    *broadieKayaScheme.hestonModel_,
                                    *broadieKayaScheme.variancePathSimulator_)
{

}

BroadieKayaScheme* BroadieKayaScheme::clone() const{
    return new BroadieKayaScheme(*this);
}

double BroadieKayaScheme::nextStep(std::size_t currentIndex, double currentValue) const
{
    //A compléter
}

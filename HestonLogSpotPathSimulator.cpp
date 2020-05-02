#include "HestonLogSpotPathSimulator.h"

HestonLogSpotPathSimulator::HestonLogSpotPathSimulator(const HestonModel& hestonModel,
                                const HestonVariancePathSimulator& variancePathSimulator):
                                PathSimulator(hestonModel),
                                variancePathSimulator_(variancePathSimulator.clone())
{

}

HestonLogSpotPathSimulator::~HestonLogSpotPathSimulator()
{
    delete variancePathSimulator_;
}

BroadieKayaScheme::BroadieKayaScheme(const HestonModel& hestonModel,
                    const HestonVariancePathSimulator& variancePathSimulator):
                    HestonLogSpotPathSimulator(hestonModel, variancePathSimulator)
{

}

BroadieKayaScheme* BroadieKayaScheme::clone() const{
    return new BroadieKayaScheme(*this);
}


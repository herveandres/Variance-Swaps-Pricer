#ifndef HESTONLOGSPOTPATHSIMULATOR_H
#define HESTONLOGSPOTPATHSIMULATOR_H

#include "PathSimulator.h"
#include "HestonVariancePathSimulator.h"

class HestonLogSpotPathSimulator : public PathSimulator
{
protected:
    const HestonModel* hestonModel_;
    const HestonVariancePathSimulator* variancePathSimulator_;

    virtual double nextStep(std::size_t currentIndex, double currentValue) const = 0; 
public:
    HestonLogSpotPathSimulator(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel,
                                const HestonVariancePathSimulator& variancePathSimulator);
    virtual ~HestonLogSpotPathSimulator();
    virtual HestonLogSpotPathSimulator* clone() const =0;
};

class BroadieKayaScheme : public HestonLogSpotPathSimulator{
private:
    double nextStep(std::size_t currentIndex, double currentValue) const;
public:
    BroadieKayaScheme(const std::vector<double>& timePoints,
                    const HestonModel& hestonModel,
                    const HestonVariancePathSimulator& variancePathSimulator);
    BroadieKayaScheme(const BroadieKayaScheme& broadieKayaScheme);
    BroadieKayaScheme* clone() const;
};

/* Optional */ 

// class EulerScheme : public HestonLogSpotPathSimulator{
// private:
//     /* data */
// public:
//     EulerScheme(/* args */);
//     ~EulerScheme();
// };

// class KahlJackelScheme : public HestonLogSpotPathSimulator{
// private:
//     /* data */
// public:
//     KahlJackelScheme(/* args */);
//     ~KahlJackelScheme();
// };

#endif
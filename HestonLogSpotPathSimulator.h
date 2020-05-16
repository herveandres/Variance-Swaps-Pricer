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

class DiscretizationScheme : public HestonLogSpotPathSimulator
{
private:
    double nextStep(std::size_t currentIndex, double currentValue) const;
    double gamma1_;
    double gamma2_;

    std::vector<double>  k0_;
    std::vector<double>  k1_;
    std::vector<double>  k2_;
    std::vector<double>  k3_;
    std::vector<double>  k4_;
    void preComputations();
public:
    DiscretizationScheme(const std::vector<double>& timePoints,
                         const HestonModel& hestonModel,
                         const HestonVariancePathSimulator& variancePathSimulator);
    DiscretizationScheme(const DiscretizationScheme& discretizationScheme);
    DiscretizationScheme* clone() const;
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

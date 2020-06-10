#ifndef HESTONLOGSPOTPATHSIMULATOR_H
#define HESTONLOGSPOTPATHSIMULATOR_H

#include "PathSimulator.h"
#include "HestonVariancePathSimulator.h"

class HestonLogSpotPathSimulator : public PathSimulator
{
protected:
    const HestonVariancePathSimulator* variancePathSimulator_;
    virtual double nextStep(std::size_t currentIndex, double currentValue, const std::vector<double>& variancePath) const = 0;
public:
    HestonLogSpotPathSimulator(const HestonVariancePathSimulator& variancePathSimulator);
    virtual ~HestonLogSpotPathSimulator();
    virtual HestonLogSpotPathSimulator* clone() const =0;
    std::vector<double> path() const;
};

class BroadieKayaScheme : public HestonLogSpotPathSimulator{
private:
    double nextStep(std::size_t currentIndex, double currentValue, const std::vector<double>& variancePath) const;
    double gamma1_;
    double gamma2_;

    std::vector<double>  k0_;
    std::vector<double>  k1_;
    std::vector<double>  k2_;
    std::vector<double>  k3_;
    std::vector<double>  k4_;
    void preComputations();
public:
    BroadieKayaScheme(const HestonVariancePathSimulator& variancePathSimulator,
                         double gamma1 = 0.5,
                         double gamma2 = 0.5);
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

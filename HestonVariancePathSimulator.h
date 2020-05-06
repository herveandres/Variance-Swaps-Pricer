#ifndef HESTONVARIANCEPATHSIMULATOR_H
#define HESTONVARIANCEPATHSIMULATOR_H


#include "PathSimulator.h"

class HestonVariancePathSimulator : public PathSimulator
{
protected:
    const HestonModel* hestonModel_;

    virtual double nextStep(std::size_t currentIndex, double currentValue) const = 0;
public:
    HestonVariancePathSimulator(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel);
    virtual ~HestonVariancePathSimulator();
    virtual HestonVariancePathSimulator* clone() const = 0;
};

class TruncatedGaussianScheme : public HestonVariancePathSimulator
{
private:
    double nextStep(std::size_t currentIndex, double currentValue) const;
public:
    TruncatedGaussianScheme(const std::vector<double>& timePoints,
                            const HestonModel& hestonModel);
    TruncatedGaussianScheme(const TruncatedGaussianScheme& truncatedGaussianScheme);
    TruncatedGaussianScheme* clone() const;
};

class QuadraticExponentialScheme : public HestonVariancePathSimulator
{
private:
    double nextStep(std::size_t currentIndex, double currentValue) const;
public:
    QuadraticExponentialScheme(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel);
    QuadraticExponentialScheme(const QuadraticExponentialScheme& quadraticExponentialScheme);
    QuadraticExponentialScheme* clone() const;
};

#endif 
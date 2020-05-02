#ifndef HESTONVARIANCEPATHSIMULATOR_H
#define HESTONVARIANCEPATHSIMULATOR_H


#include "PathSimulator.h"

class HestonVariancePathSimulator : public PathSimulator
{
private:
    /* data */
public:
    HestonVariancePathSimulator(const HestonModel& hestonModel);
    virtual ~HestonVariancePathSimulator();
    virtual HestonVariancePathSimulator* clone() const = 0;
};

class TruncatedGaussianScheme : public HestonVariancePathSimulator
{
private:
    /* data */
public:
    TruncatedGaussianScheme(const HestonModel& hestonModel);
    TruncatedGaussianScheme* clone() const;
};

class QuadraticExponentialScheme : public HestonVariancePathSimulator
{
private:
    /* data */
public:
    QuadraticExponentialScheme(const HestonModel& hestonModel);
    QuadraticExponentialScheme* clone() const;
};

#endif 
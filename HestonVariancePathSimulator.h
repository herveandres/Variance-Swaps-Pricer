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

class TruncatedGaussian : public HestonVariancePathSimulator
{
private:
    /* data */
public:
    TruncatedGaussian(const HestonModel& hestonModel);
    TruncatedGaussian* clone() const;
};

class QuadraticExponential : public HestonVariancePathSimulator
{
private:
    /* data */
public:
    QuadraticExponential(const HestonModel& hestonModel);
    QuadraticExponential* clone() const;
};

#endif 
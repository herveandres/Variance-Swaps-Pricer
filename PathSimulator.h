#ifndef PATHSIMULATOR_H
#define PATHSIMULATOR_H

#include <vector>
#include "Model.h"

class PathSimulator
{
protected:
    double initialValue_;
    std::vector<double> timePoints_;
    virtual double nextStep(std::size_t currentIndex, double currentValue) const = 0; 
public:
    PathSimulator(double initialValue, const std::vector<double>& timePoints);
    virtual ~PathSimulator();
    virtual PathSimulator* clone() const = 0;
    std::vector<double> path() const;
};

#endif // !
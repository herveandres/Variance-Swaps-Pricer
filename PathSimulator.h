#ifndef PATHSIMULATOR_H
#define PATHSIMULATOR_H

#include <vector>
#include <random>
#include "Model.h"

class PathSimulator
{
protected:
    double initialValue_;
    std::vector<double> timePoints_;
public:
    PathSimulator(double initialValue, const std::vector<double>& timePoints);
    virtual ~PathSimulator();
    virtual PathSimulator* clone() const = 0;
    virtual std::vector<double> path() const = 0; 
    std::vector<double> getTimePoints() const;
};

#endif // !
#ifndef PATHSIMULATOR_H
#define PATHSIMULATOR_H

#include <vector>
#include <random>
#include "Model.h"

//Abstract class
class PathSimulator
{
protected:
    double initialValue_; // Each simulated path starts from the same initial value.
    std::vector<double> timePoints_; // Time interval which is dicretized in points.
public:
    PathSimulator(double initialValue, const std::vector<double>& timePoints);
    virtual ~PathSimulator();
    virtual PathSimulator* clone() const = 0;

    //Method simulating a random path
    virtual std::vector<double> path() const = 0; 
    std::vector<double> getTimePoints() const;
};

#endif // !
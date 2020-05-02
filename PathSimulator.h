#ifndef PATHSIMULATOR_H
#define PATHSIMULATOR_H

#include "Model.h"

class PathSimulator
{
private:
    Model* model_;
public:
    PathSimulator(const Model& model);
    virtual ~PathSimulator();
    virtual PathSimulator* clone() const = 0;
};

#endif // !
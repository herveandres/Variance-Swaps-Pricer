#include "PathSimulator.h"

 PathSimulator::PathSimulator(double initialValue, const std::vector<double>& timePoints):
    initialValue_(initialValue), timePoints_(timePoints)
 {

 }

PathSimulator::~PathSimulator()
{

}

std::vector<double> PathSimulator::path() const
{
	std::vector<double> path {initialValue_};
	for (std::size_t index = 0; index < timePoints_.size() - 1; ++index)
		path.push_back(nextStep(index, path[index]));

	return path;
}
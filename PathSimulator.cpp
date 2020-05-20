#include "PathSimulator.h"

 PathSimulator::PathSimulator(double initialValue, 
 							  const std::vector<double>& timePoints):
    initialValue_(initialValue), timePoints_(timePoints)
 {

 }

PathSimulator::~PathSimulator()
{

}

std::vector<double> PathSimulator::getTimePoints() const
{
	return timePoints_;
}
#include <iostream>
#include "HestonVariancePathSimulator.h"
#include "MathFunctions.h"

using namespace std;

int main()
{
    
    HestonModel hestonModel(0,0.5,0.04,1.0,-0.9,0.04,100);
    double maturity = 1;
    vector<double> timePoints = MathFunctions::buildLinearSpace(0,maturity,1000);

    cout << "Test of the variance simulation using the Truncated Gaussian method" << endl;
    TruncatedGaussianScheme truncatedGaussianScheme(timePoints,hestonModel);

    std::vector<double> pathTruncatedGaussian = truncatedGaussianScheme.path();
    // for(size_t i = 0; i < pathTruncatedGaussian.size(); i++)
    // {
    //     std::cout << pathTruncatedGaussian[i] << std::endl;
    // }

    cout << "Test of the variance simulation using the Quadratic Exponential method" << endl;
    QuadraticExponentialScheme quadraticExponentialScheme(timePoints,hestonModel);
    std::vector<double> pathQuadraticExponential = quadraticExponentialScheme.path();
    for(size_t i = 0; i < pathQuadraticExponential.size(); i++)
    {
        std::cout << pathQuadraticExponential[i] << std::endl;
    }

    return 0;
}
#include <iostream>
#include "HestonVariancePathSimulator.h"
#include "MathFunctions.h"

using namespace std;

int main()
{
    cout << "Test of the variance simulation using the Truncated Gaussian method" << endl;
    
    HestonModel hestonModel(0,0.5,0.04,1.0,-0.9,0.04,100);
    double maturity = 1;
    vector<double> timePoints = MathFunctions::buildLinearSpace(0,maturity,1000);
    TruncatedGaussianScheme truncatedGaussianScheme(timePoints,hestonModel);

    std::vector<double> pathTruncatedGaussian = truncatedGaussianScheme.path();
    for(size_t i = 0; i < pathTruncatedGaussian.size(); i++)
    {
        std::cout << pathTruncatedGaussian[i] << std::endl;
    }
    return 0;
}
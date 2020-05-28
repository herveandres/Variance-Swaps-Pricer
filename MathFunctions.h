#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H


#include <cmath>
#include <complex>
#include <functional>
#include <vector>
#include <random>
#include <iostream>


namespace MathFunctions
{
    double normalPDF(double x);

    double normalCDF(double x);

    double normalCDFInverse(double x);

    extern unsigned seed;
    extern std::default_random_engine generator;

    double simulateUniformRandomVariable();
    //Simulate a standard Gaussian variable using Box-Muller method
    double simulateGaussianRandomVariable();

    double newtonMethod(double initialGuess, std::function<double(double)> f, std::function<double(double)> fPrime, double precision = pow(10,-5));

    std::complex<double> differencesFinies(std::function<std::complex<double>(double,double)> f, double omega, double tau, double epsilon=pow(10,-2));

    //Function looking for the index i s.t. list[i] <= x < list[i+1] (list is sorted)
    std::size_t binarySearch(const std::vector<double>& list, double x);

    //Function building a linear interval [min,max] with N points
    std::vector<double> buildLinearSpace(double min, double max, std::size_t N);
}

#endif // !MATHFUNCTIONS_H

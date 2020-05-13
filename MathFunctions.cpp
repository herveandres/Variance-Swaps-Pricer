#include "MathFunctions.h"
#include <cstdlib>

namespace MathFunctions
{
    double normalPDF(double x)
    {
        return std::exp(-x*x/2)/std::sqrt(2*M_PI);
    }

    double normalCDF(double x)
    {
        return 0.5 * std::erfc(-x * M_SQRT1_2);
    }

    unsigned seed = 10;
    std::default_random_engine generator(seed);

    double simulateGaussianRandomVariable()
    {
        std::uniform_real_distribution<double> distribution(0,1);
        double u = distribution(generator), v = distribution(generator);
        return std::sqrt(-2*std::log(u))*std::sin(2*M_PI*v);
    }

    double newtonMethod(double initialGuess, std::function<double(double)> f, std::function<double(double)> fPrime, double precision)
    {
        double x = initialGuess;
        double fx = f(x);
        while(std::abs(fx) > precision)
        {
            x -= fx/fPrime(x);
            fx = f(x);
        }
        std::complex<double> j;
        return x;      
    }

    std::complex<double> differencesFinies(std::function<std::complex<double>(double,double)> f, double omega, double tau, double epsilon)
    {
        return (f(tau,omega+epsilon)-f(tau,omega-epsilon))/(2*epsilon);
    }

    //Function looking for the index i s.t. list[i] <= x < list[i+1] (list is sorted)
    std::size_t binarySearch(const std::vector<double>& list, double x)
    {
        if(list.back() <= x)
        {
            return list.size()-1;
        }
        else if (list.front() > x)
        {   
            return 0;
        }
        else
        {
            std::size_t begin = 0, end = list.size()-1;
            std::size_t idx = (begin+end)/2;
            while(!(list[idx] <= x && list[idx+1] > x))
            {
                if(list[idx] <= x)
                    begin = idx+1;
                else
                    end = idx;
                idx = (begin+end)/2;
            }
            return idx;
        }
        
    } 

    std::vector<double> buildLinearSpace(double min, double max, std::size_t N)
    {
        double delta = (max-min)/double(N-1);
        std::vector<double> linearSpace;

        for(std::size_t i = 0; i < N; i++)
        {
            linearSpace.push_back(min+i*delta);
        }  
        return linearSpace;
    }
}


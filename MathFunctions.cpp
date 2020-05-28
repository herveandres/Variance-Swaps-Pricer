#include "MathFunctions.h"
#include <cstdlib>
#include <random>

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


    double normalCDFInverse(double x)
    {
        static const double a0 = 2.50662823884;
        double a1 = -18.61500062529;
        double a2 = 41.39119773534;
        double a3 =-25.44106049637;
        double b0  = -8.47351093090;
        double b1 = 23.08336743743;
        double b2 =-21.06224101826;
        double b3 =  3.13082909833;
        double c0 = 0.3374754822726147;
        double c1 = 0.9761690190917186;
        double c2 = 0.1607979714918209;
        double c3 = 0.0276438810333863;
        double c4 = 0.0038405729373609;
        double c5 = 0.0003951896511919;
        double c6 =  0.0000321767881768;
        double c7 = 0.0000002888167364;
        double c8 = 0.0000003960315187;

        double result;
        double temp = x - 0.5;

        if (abs(temp)<0.42){
            result = temp*temp;
            result = temp*
                    (((a3*result+a2)*result+a1)*result+a0) /
                    ((((b3*result+b2)*result+b1)*result+b0)*result+1.0);
        } else{
            if (x<0.5)
                result = x;
            else
                result=1.0-x;
            result = std::log(-std::log(result));
            result = c0+result*(c1+result*(c2+result*(c3+result*
                                                        (c4+result*(c5+result*(c6+result*
                                                                                (c7+result*c8)))))));
            if (x<0.5)
                result=-result;
        }
        return result;
    }

    unsigned seed = 10;
    std::default_random_engine generator(seed);

    double simulateUniformRandomVariable()
    {
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        return distribution(generator);
    }
    double simulateGaussianRandomVariable()
    {
        double u = simulateUniformRandomVariable(), v = simulateUniformRandomVariable();
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


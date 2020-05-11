#include <cmath>
#include <functional>

namespace MathFunctions
{
    double normalPDF(double x)
    {
        return exp(-x*x/2)/sqrt(2*M_PI);
    }

    double normalCDF(double x)
    {
        return 0.5 * erfc(-x * M_SQRT1_2);
    }

    double newtonMethod(double initialGuess, std::function<double(double)> f, std::function<double(double)> fPrime, double precision = pow(10,-5))
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
};
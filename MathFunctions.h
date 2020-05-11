#include <cmath>

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
};
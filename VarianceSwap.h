#ifndef VARIANCESWAP_H
#define VARIANCESWAP_H

#include <vector>

using namespace std;

class VarianceSwap
{
private:
    double maturity_;
    size_t nbOfObservations_;
    vector<double> dates_;
public:
    VarianceSwap(double maturity, size_t nbOfObservations);
    ~VarianceSwap();
};

#endif
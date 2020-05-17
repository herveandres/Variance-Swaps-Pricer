#ifndef VARIANCESWAP_H
#define VARIANCESWAP_H

#include <vector>

class VarianceSwap
{
private:
    std::vector<double> dates_;
public:
    VarianceSwap(double maturity, std::size_t nbOfObservations);
    ~VarianceSwap();
    std::vector<double> getDates() const;
};

#endif
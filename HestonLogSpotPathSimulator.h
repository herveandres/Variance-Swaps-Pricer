#ifndef HESTONLOGSPOTPATHSIMULATOR_H
#define HESTONLOGSPOTPATHSIMULATOR_H

#include "PathSimulator.h"
#include "HestonVariancePathSimulator.h"

//Abstract class
class HestonLogSpotPathSimulator : public PathSimulator
{
protected:
    const HestonVariancePathSimulator* variancePathSimulator_;
    virtual double nextStep(std::size_t currentIndex, double currentValue, const std::vector<double>& variancePath) const = 0;
public:
    HestonLogSpotPathSimulator(const HestonVariancePathSimulator& variancePathSimulator);

    // Copy constructor, Assignement operator and Destructor are needed because one of the member variable is a pointer
    HestonLogSpotPathSimulator(const HestonLogSpotPathSimulator& logSpotPathSimulator); 
    virtual ~HestonLogSpotPathSimulator();
    HestonLogSpotPathSimulator& operator=(const HestonLogSpotPathSimulator& logSpotPathSimulator);

    virtual HestonLogSpotPathSimulator* clone() const =0;
    std::vector<double> path() const;
};

class BroadieKayaScheme : public HestonLogSpotPathSimulator{
private:
    double nextStep(std::size_t currentIndex, double currentValue, const std::vector<double>& variancePath) const;
    
    //Coefficients used for the approximation of the integral of V
    double gamma1_;
    double gamma2_;

    //Pre-computed coefficients of the diffusion that are path independent
    std::vector<double>  k0_;
    std::vector<double>  k1_;
    std::vector<double>  k2_;
    std::vector<double>  k3_;
    std::vector<double>  k4_;

    //Pre-computes k0_, k1_, k_2_, k3_ and k4_
    void preComputations();
public:
    BroadieKayaScheme(const HestonVariancePathSimulator& variancePathSimulator,
                         double gamma1 = 0.5,  //Default is central discretization
                         double gamma2 = 0.5); 
    BroadieKayaScheme(const BroadieKayaScheme& broadieKayaScheme);
    ~BroadieKayaScheme() = default;
    BroadieKayaScheme* clone() const;
};


#endif

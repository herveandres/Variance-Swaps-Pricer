#ifndef HESTONVARIANCEPATHSIMULATOR_H
#define HESTONVARIANCEPATHSIMULATOR_H


#include "PathSimulator.h"

class HestonVariancePathSimulator : public PathSimulator
{
protected:
    const HestonModel* hestonModel_;

    virtual double nextStep(std::size_t currentIndex, double currentValue) const = 0;
public:
    HestonVariancePathSimulator(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel);
    virtual ~HestonVariancePathSimulator();
    virtual HestonVariancePathSimulator* clone() const = 0;
};

class TruncatedGaussianScheme : public HestonVariancePathSimulator
{
private:
    /* Function that pre-computes some quantities that will be used in nextStep and caches them 
    NB : we allow the time grid to be non-equidistant so that the computed quantities are time dependent.*/
    void preComputations();
    double nextStep(std::size_t currentIndex, double currentValue) const;
    static double h(double r, double psi);
    static double hPrime(double r, double psi); 
    
    const double confidenceMultiplier_;
    std::vector<double> k1_;
    std::vector<double> k2_;
    std::vector<double> k3_;
    std::vector<double> k4_;
    std::vector<double> psiGrid_;
    double initialGuess_;
    std::vector<double> fmu_;
    std::vector<double> fsigma_;

public:
    TruncatedGaussianScheme(const std::vector<double>& timePoints,
                            const HestonModel& hestonModel, 
                            double confidenceMultiplier = 4,
                            std::size_t psiGridSize = 100,
                            double initialGuess = 1);
    //Copy constructor
    TruncatedGaussianScheme(const TruncatedGaussianScheme& truncatedGaussianScheme);
    TruncatedGaussianScheme* clone() const;
};

class QuadraticExponentialScheme : public HestonVariancePathSimulator
{
private:
    std::vector<double> k1_;
    std::vector<double> k2_;
    std::vector<double> k3_;
    std::vector<double> k4_;
    double psiC_;
    double kStar_;
    void preComputations();
    double nextStep(std::size_t currentIndex, double currentValue) const;

public:
    QuadraticExponentialScheme(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel, double psiC = 0.5);
    //Copy constructor
    QuadraticExponentialScheme(const QuadraticExponentialScheme& quadraticExponentialScheme);
    QuadraticExponentialScheme* clone() const;
};

#endif 

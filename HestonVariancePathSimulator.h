#ifndef HESTONVARIANCEPATHSIMULATOR_H
#define HESTONVARIANCEPATHSIMULATOR_H


#include "PathSimulator.h"

//Abstract class
class HestonVariancePathSimulator : public PathSimulator
{
protected:
    const HestonModel* hestonModel_;
    /* Function that pre-computes some quantities that will be used in nextStep and caches them */
    void preComputations();
    virtual double nextStep(std::size_t currentIndex, double currentValue) const = 0;
    std::vector<double> k1_;
    std::vector<double> k2_;
    std::vector<double> k3_;
    std::vector<double> k4_;
public:
    HestonVariancePathSimulator(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel);
    
    // Copy constructor, Assignement operator and Destructor are needed because one of the member variable is a pointer
    HestonVariancePathSimulator(const HestonVariancePathSimulator& variancePathSimulator);
    virtual ~HestonVariancePathSimulator();
    HestonVariancePathSimulator& operator=(
                        const HestonVariancePathSimulator& variancePathSimulator);
    
    virtual HestonVariancePathSimulator* clone() const = 0;
    std::vector<double> path() const;
    HestonModel getHestonModel() const;
};

class TruncatedGaussianScheme : public HestonVariancePathSimulator
{
private:
    //Pre-computation of f_mu and f_sigma on psiGrid
    void preComputationsTG();
    double nextStep(std::size_t currentIndex, double currentValue) const;

    /*Function that the r used to compute f_mu and f_sigma must nullifies. 
    It is declared as static since it doesn't need any attribute from the class*/
    static double h(double r, double psi);
    /*Function h's derivative with respect to r
    It is declared as static since it doesn't need any attribute from the class*/
    static double hPrime(double r, double psi); 
    
    /*Parameter alpha in the article that is used to compute the lower bound of
    the grid for psi */
    const double confidenceMultiplier_;

    //Grid for psi on which f_mu and f_sigma are computed
    std::vector<double> psiGrid_;

    //Initial guess for r inputed in Newton method
    double initialGuess_;

    std::vector<double> fmu_;
    std::vector<double> fsigma_;

public:
    TruncatedGaussianScheme(const std::vector<double>& timePoints,
                            const HestonModel& hestonModel, 
                            /*We chose to use a smaller confidence multiplier than suggested
                            in Andersen's article since we observed numerical instabilities 
                            for the first values of psi in the grid when we used a bigger 
                            confidence multiplier*/
                            double confidenceMultiplier = 2, 
                            //The two following default values are empirically chosen
                            std::size_t psiGridSize = 100,
                            double initialGuess = 1);

    TruncatedGaussianScheme(const TruncatedGaussianScheme& truncatedGaussianScheme);
    ~TruncatedGaussianScheme() = default;
    TruncatedGaussianScheme* clone() const;
};

class QuadraticExponentialScheme : public HestonVariancePathSimulator
{
protected:
    //Switching threshold
    double psiC_;

    double nextStep(std::size_t currentIndex, double currentValue) const;

public:
    QuadraticExponentialScheme(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel, double psiC = 1.5);

    QuadraticExponentialScheme(const QuadraticExponentialScheme& quadraticExponentialScheme);
    ~QuadraticExponentialScheme() = default;
    QuadraticExponentialScheme* clone() const;
};

class QuadraticExponentialMartingaleCorrectionScheme: public QuadraticExponentialScheme
{
private:
    mutable bool mCCase; // true = (psi <= psiC) and false = (psi> psiC)
    mutable double mCCoeff1; // will be used as a or p, depending on condition on A in LogSpotPathSimulator process
    mutable double mCCoeff2; // will be used as b or beta, depending on condition on A in LogSpotPathSimulator process
    double nextStep(std::size_t currentIndex, double currentValue) const;
public:
    QuadraticExponentialMartingaleCorrectionScheme(const std::vector<double>& timePoints,
                                const HestonModel& hestonModel, double psiC = 1.5);
    //Copy constructor
    QuadraticExponentialMartingaleCorrectionScheme(const QuadraticExponentialMartingaleCorrectionScheme& quadraticExponentialScheme);
    QuadraticExponentialMartingaleCorrectionScheme* clone() const;

    ~QuadraticExponentialMartingaleCorrectionScheme() = default;
    bool getMCCase() const;
    double getMCCoeff1() const;
    double getMCCoeff2() const;
};

#endif 

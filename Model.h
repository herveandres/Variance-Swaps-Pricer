#ifndef MODEL_H
#define MODEL_H


class HestonModel
{
private:
    double riskFreeRate_;
    double drift_;
    double meanReversionSpeed_;
    double meanReversionLevel_;
    double volOfVol_;
    double correlation_;
    double initialVolatility_;
    double initialAssetValue_;
public:
    HestonModel(double riskFreeRate,
                double drift,
                double meanReversionSpeed,
                double meanReversionLevel,
                double volOfVol,
                double correlation,
                double initialVolatility,
                double initialAssetValue);
    HestonModel(const HestonModel& hestonModel);
    double getRiskFreeRate() const;
    double getDrift() const;
    double getMeanReversionSpeed() const;
    double getMeanReversionLevel() const;
    double getVolOfVol() const;
    double getCorrelation() const;
    double getInitialVolatility() const;
    double getInitialAssetValue() const;
};

#endif // !MODEL_H
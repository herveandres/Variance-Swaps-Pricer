#include "Model.h"

HestonModel::HestonModel(double drift,
                        double meanReversionSpeed,
                        double meanReversionLevel,
                        double volOfVol,
                        double correlation,
                        double initialVolatility,
                        double initialAssetValue):
                        Model(),
                        drift_(drift),
                        meanReversionSpeed_(meanReversionSpeed),
                        meanReversionLevel_(meanReversionLevel),
                        volOfVol_(volOfVol),
                        correlation_(correlation),
                        initialVolatility_(initialVolatility),
                        initialAssetValue_(initialAssetValue)
{
    
}

HestonModel::HestonModel(const HestonModel& hestonModel):
    drift_(hestonModel.drift_), meanReversionSpeed_(hestonModel.meanReversionSpeed_),
    meanReversionLevel_(hestonModel.meanReversionLevel_), 
    volOfVol_(hestonModel.volOfVol_), correlation_(hestonModel.correlation_),
    initialVolatility_(hestonModel.initialVolatility_), 
    initialAssetValue_(hestonModel.initialAssetValue_)
{

}

double HestonModel::getInitialVolatility() const
{
    return initialVolatility_;
}

double HestonModel::getInitialAssetValue() const
{
    return initialAssetValue_;
}
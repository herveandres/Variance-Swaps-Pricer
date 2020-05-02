#include "Model.h"

Model::Model(){

}
Model::~Model(){
    
}

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

HestonModel* HestonModel::clone() const{
    return new HestonModel(*this);
}
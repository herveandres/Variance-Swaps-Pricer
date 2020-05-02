#ifndef MODEL_H
#define MODEL_H

class Model
{
private:
    /* data */
public:
    Model(/* args */);
    virtual ~Model();
    virtual Model* clone() const = 0;
};


class HestonModel : public Model
{
private:
    double drift_;
    double meanReversionSpeed_;
    double meanReversionLevel_;
    double volOfVol_;
    double correlation_;
    double initialVolatility_;
    double initialAssetValue_;
public:
    HestonModel(double drift,
                double meanReversionSpeed,
                double meanReversionLevel,
                double volOfVol,
                double correlation,
                double initialVolatility,
                double initialAssetValue);
    HestonModel* clone() const;
};

#endif // !MODEL_H
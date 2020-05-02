#ifndef HESTONLOGSPOTPATHSIMULATOR_H
#define HESTONLOGSPOTPATHSIMULATOR_H

#include "PathSimulator.h"
#include "HestonVariancePathSimulator.h"

class HestonLogSpotPathSimulator : public PathSimulator
{
protected:
    const HestonVariancePathSimulator* variancePathSimulator_;
public:
    HestonLogSpotPathSimulator(const HestonModel& hestonModel,
                                const HestonVariancePathSimulator& variancePathSimulator);
    virtual ~HestonLogSpotPathSimulator();
    virtual HestonLogSpotPathSimulator* clone() const =0;
};

class BroadieKayaScheme : public HestonLogSpotPathSimulator{
private:
    /* data */
public:
    BroadieKayaScheme(const HestonModel& hestonModel,
                    const HestonVariancePathSimulator& variancePathSimulator);
    BroadieKayaScheme* clone() const;
};

/* Optional */ 

// class EulerScheme : public HestonLogSpotPathSimulator{
// private:
//     /* data */
// public:
//     EulerScheme(/* args */);
//     ~EulerScheme();
// };

// class KahlJackelScheme : public HestonLogSpotPathSimulator{
// private:
//     /* data */
// public:
//     KahlJackelScheme(/* args */);
//     ~KahlJackelScheme();
// };

#endif
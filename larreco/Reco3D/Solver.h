// Christopher Backhouse - bckhouse@fnal.gov

#ifndef RECO3D_SOLVER_H
#define RECO3D_SOLVER_H

#include <vector>

#include "QuadExpr.h"

class InductionWireHit
{
public:
  InductionWireHit(int chan, double t, double q);

  //protected:
  int fChannel;
  double fTime;

  double fCharge;

  double fPred;
};

class SpaceCharge;
class Neighbour
{
public:
  Neighbour(SpaceCharge* sc, double dist, double coupling);

  SpaceCharge* fSC;
  double fDist;
  double fCoupling;
};

class CollectionWireHit;

class SpaceCharge
{
public:
  SpaceCharge(double t, double x, double y, double z,
              CollectionWireHit* cwire,
              InductionWireHit* wire1, InductionWireHit* wire2);

  //protected:
  double fTime, fX, fY, fZ;
  CollectionWireHit* fCWire;
  InductionWireHit *fWire1, *fWire2;

  std::vector<Neighbour> fNeighbours;

  double fPred;
  double fNeiPotential; ///< Neighbour-induced potential
};

class CollectionWireHit
{
public:
  CollectionWireHit(int chan, double t, double q, const std::vector<SpaceCharge*>& cross);
  ~CollectionWireHit();

  //protected:
  int fChannel;
  double fTime;

  double fCharge;

  std::vector<double> fWeights; ///< Distribution of charge over crossings
  std::vector<SpaceCharge*> fCrossings;
};

double Metric(const std::vector<SpaceCharge*>& scs, double alpha);
double Metric(const std::vector<CollectionWireHit*>& cwires, double alpha);
QuadExpr Metric(const SpaceCharge* sci, const SpaceCharge* scj, double Q, double alpha);

double SolvePair(CollectionWireHit* cwire,
                 SpaceCharge* sci, SpaceCharge* scj,
                 double xmin, double xmax,
                 double alpha);
void Iterate(CollectionWireHit* cwire,
             double alpha);
void Iterate(const std::vector<CollectionWireHit*>& cwires,
             double alpha);

#endif

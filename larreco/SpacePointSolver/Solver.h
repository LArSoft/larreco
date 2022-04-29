// Christopher Backhouse - bckhouse@fnal.gov

#ifndef RECO3D_SOLVER_H
#define RECO3D_SOLVER_H

#include <vector>

#include "larreco/SpacePointSolver/QuadExpr.h"

/// Allow InductionWireHit and CollectionWireHit to be put in the same maps
/// where necessary.
class WireHit
{
  // Do not add any virtual functions: will increase memory usage of all wires
};

class InductionWireHit: public WireHit
{
public:
  InductionWireHit(int chan, double q);

  //protected:
  int fChannel;

  double fCharge;

  double fPred;
};

class SpaceCharge;
class Neighbour
{
public:
  Neighbour(SpaceCharge* sc, double coupling);

  SpaceCharge* fSC;
  double fCoupling;
};

class CollectionWireHit;

class SpaceCharge
{
public:
  SpaceCharge(double x, double y, double z,
              CollectionWireHit* cwire,
              InductionWireHit* wire1, InductionWireHit* wire2);

  void AddCharge(double dq);

  //protected:
  double fX, fY, fZ;
  CollectionWireHit* fCWire;
  InductionWireHit *fWire1, *fWire2;

  std::vector<Neighbour> fNeighbours;

  double fPred;
  double fNeiPotential; ///< Neighbour-induced potential
};

class CollectionWireHit: public WireHit
{
public:
  CollectionWireHit(int chan, double q, const std::vector<SpaceCharge*>& cross);
  ~CollectionWireHit();

  //protected:
  int fChannel;

  double fCharge;

  std::vector<SpaceCharge*> fCrossings;
};

double Metric(const std::vector<SpaceCharge*>& scs, double alpha);
double Metric(const std::vector<CollectionWireHit*>& cwires, double alpha);
QuadExpr Metric(const SpaceCharge* sci, const SpaceCharge* scj, double alpha);

double SolvePair(CollectionWireHit* cwire,
                 SpaceCharge* sci, SpaceCharge* scj,
                 double xmin, double xmax,
                 double alpha);
void Iterate(CollectionWireHit* cwire, double alpha);
void Iterate(SpaceCharge* sc, double alpha);
void Iterate(const std::vector<CollectionWireHit*>& cwires,
             const std::vector<SpaceCharge*>& orphanSCs,
             double alpha);

#endif

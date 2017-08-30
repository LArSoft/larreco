#include "Solver.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <set>
#include <iostream>

template<class T> T sqr(T x){return x*x;}

// ---------------------------------------------------------------------------
InductionWireHit::InductionWireHit(int chan, double t, double q)
  : fChannel(chan), fTime(t), fCharge(q), fPred(0)
{
}

// ---------------------------------------------------------------------------
Neighbour::Neighbour(SpaceCharge* sc, double coupling)
  : fSC(sc), fCoupling(coupling)
{
}

// ---------------------------------------------------------------------------
SpaceCharge::SpaceCharge(double x, double y, double z,
                         CollectionWireHit* cwire,
                         InductionWireHit* wire1, InductionWireHit* wire2)
  : fX(x), fY(y), fZ(z),
    fCWire(cwire), fWire1(wire1), fWire2(wire2),
    fPred(0),
    fNeiPotential(0)
{
}

// ---------------------------------------------------------------------------
void SpaceCharge::AddCharge(double dq)
{
  fPred += dq;

  for(Neighbour& nei: fNeighbours)
    nei.fSC->fNeiPotential += dq * nei.fCoupling;

  fWire1->fPred += dq;
  fWire2->fPred += dq;
}

// ---------------------------------------------------------------------------
CollectionWireHit::CollectionWireHit(int chan, double t, double q,
                                     const std::vector<SpaceCharge*>& cross)
  : fChannel(chan), fTime(t), fCharge(q),
    fCrossings(cross)
{
  const double p = q/cross.size();

  for(SpaceCharge* iwires: cross){
    iwires->fPred += p;
    iwires->fWire1->fPred += p;
    iwires->fWire2->fPred += p;
  }
}

// ---------------------------------------------------------------------------
CollectionWireHit::~CollectionWireHit()
{
  // Each SpaceCharge can only be in one collection wire, so we own them. But
  // not SpaceCharge does not clean up its induction wires since they're
  // shared.
  for(SpaceCharge* sc: fCrossings) delete sc;
}

// ---------------------------------------------------------------------------
double Metric(double q, double p)
{
  return sqr(q-p);
}

// ---------------------------------------------------------------------------
QuadExpr Metric(double q, QuadExpr p)
{
  return sqr(q-p);
}

// ---------------------------------------------------------------------------
double Metric(const std::vector<SpaceCharge*>& scs, double alpha)
{
  double ret = 0;

  std::set<InductionWireHit*> iwires;
  for(const SpaceCharge* sc: scs){
    iwires.insert(sc->fWire1);
    iwires.insert(sc->fWire2);

    if(alpha != 0){
      ret -= alpha*sqr(sc->fPred);
      // "Double-counting" of the two ends of the connection is
      // intentional. Otherwise we'd have a half in the line above.
      ret -= alpha * sc->fPred * sc->fNeiPotential;
    }
  }

  for(const InductionWireHit* iwire: iwires){
    ret += Metric(iwire->fCharge, iwire->fPred);
  }

  return ret;
}

// ---------------------------------------------------------------------------
double Metric(const std::vector<CollectionWireHit*>& cwires, double alpha)
{
  std::vector<SpaceCharge*> scs;
  for(CollectionWireHit* cwire: cwires)
    scs.insert(scs.end(), cwire->fCrossings.begin(), cwire->fCrossings.end());
  return Metric(scs, alpha);
}

// ---------------------------------------------------------------------------
QuadExpr Metric(const SpaceCharge* sci, const SpaceCharge* scj, double alpha)
{
  QuadExpr ret = 0;

  // How much charge moves from scj to sci
  QuadExpr x = QuadExpr::X();

  const InductionWireHit* iwire1 = sci->fWire1;
  const InductionWireHit* iwire2 = sci->fWire2;
  const InductionWireHit* jwire1 = scj->fWire1;
  const InductionWireHit* jwire2 = scj->fWire2;

  const double scip = sci->fPred;
  const double scjp = scj->fPred;

  const double qi1 = iwire1->fCharge;
  const double qi2 = iwire2->fCharge;
  const double pi1 = iwire1->fPred;
  const double pi2 = iwire2->fPred;

  const double qj1 = jwire1->fCharge;
  const double qj2 = jwire2->fCharge;
  const double pj1 = jwire1->fPred;
  const double pj2 = jwire2->fPred;

  if(alpha != 0){
    // Self energy. SpaceCharges are never the same object
    ret -= alpha*sqr(scip + x);
    ret -= alpha*sqr(scjp - x);

    // Interaction. We're only seeing one end of the double-ended connection
    // here, so multiply by two.
    ret -= 2 * alpha * (scip + x) * sci->fNeiPotential;
    ret -= 2 * alpha * (scjp - x) * scj->fNeiPotential;
    // TODO - this miscounts if i and j are neighbours of each other?
  }

  if(iwire1 == jwire1){
    ret += Metric(qi1, pi1);
  }
  else{
    ret += Metric(qi1, pi1 + x);
    ret += Metric(qj1, pj1 - x);
  }

  if(iwire2 == jwire2){
    ret += Metric(qi2, pi2);
  }
  else{
    ret += Metric(qi2, pi2 + x);
    ret += Metric(qj2, pj2 - x);
  }

  return ret;
}

// ---------------------------------------------------------------------------
double SolvePair(CollectionWireHit* cwire,
                 SpaceCharge* sci, SpaceCharge* scj,
                 double alpha)
{
  const QuadExpr chisq = Metric(sci, scj, alpha);
  const double chisq0 = chisq.Eval(0);

  // Find the minimum of a quadratic expression
  double x = -chisq.Linear()/(2*chisq.Quadratic());

  // Don't allow either SpaceCharge to go negative
  const double xmin = -sci->fPred;
  const double xmax =  scj->fPred;

  // Clamp to allowed range
  x = std::min(xmax, x);
  x = std::max(xmin, x);

  const double chisq_new = chisq.Eval(x);

  // Should try these too, because the function might be convex not concave, so
  // d/dx=0 gives the max not the min, and the true min is at one extreme of
  // the range.
  const double chisq_p = chisq.Eval(xmax);
  const double chisq_n = chisq.Eval(xmin);

  if(std::min(std::min(chisq_p, chisq_n), chisq_new) > chisq0+1){
    for(double x = xmin; x < xmax; x += .01*(xmax-xmin)){
      std::cout << x << " " << chisq.Eval(x) << std::endl;
    }

    std::cout << chisq_new << " " << chisq0 << " " << chisq_p << " " << chisq_n << std::endl;
    abort();
  }

  if(std::min(chisq_n, chisq_p) < chisq_new){
    if(chisq_n < chisq_p) return xmin;
    return xmax;
  }

  return x;
}

// ---------------------------------------------------------------------------
void Iterate(CollectionWireHit* cwire, double alpha)
{
  // Consider all pairs of crossings
  const unsigned int N = cwire->fCrossings.size();

  for(unsigned int i = 0; i+1 < N; ++i){
    SpaceCharge* sci = cwire->fCrossings[i];

    InductionWireHit* iwire1 = sci->fWire1;
    InductionWireHit* iwire2 = sci->fWire2;

    for(unsigned int j = i+1; j < N; ++j){
      SpaceCharge* scj = cwire->fCrossings[j];

      InductionWireHit* jwire1 = scj->fWire1;
      InductionWireHit* jwire2 = scj->fWire2;

      // There can't be any cross-overs of U and V views like this
      if(iwire1 == jwire2 || iwire2 == jwire1) abort();

      // Driving all the same wires, no move can have any effect
      if(iwire1 == jwire1 && iwire2 == jwire2) continue;

      const double x = SolvePair(cwire, sci, scj, alpha);

      if(x == 0) continue;

      // Actually make the update
      sci->AddCharge(+x);
      scj->AddCharge(-x);
    } // end for j
  } // end for i
}

// ---------------------------------------------------------------------------
void Iterate(const std::vector<CollectionWireHit*>& cwires,
             double alpha)
{
  for(CollectionWireHit* cwire: cwires) Iterate(cwire, alpha);
}

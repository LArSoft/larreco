#include "Solver.h"

#include <algorithm>
#include <cstdlib>
#include <set>
#include <string>
#include <iostream>

template<class T> T sqr(T x){return x*x;}

// ---------------------------------------------------------------------------
InductionWireHit::InductionWireHit(int chan, double q)
  : fChannel(chan), fCharge(q), fPred(0)
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

  if(fWire1) fWire1->fPred += dq;
  if(fWire2) fWire2->fPred += dq;
}

// ---------------------------------------------------------------------------
CollectionWireHit::CollectionWireHit(int chan, double q,
                                     const std::vector<SpaceCharge*>& cross)
  : fChannel(chan), fCharge(q), fCrossings(cross)
{
  if(q < 0){
    std::cout << "Trying to construct collection wire with negative charge " << q << " this should never happen." << std::endl;
    abort();
  }

  const double p = q/cross.size();

  for(SpaceCharge* sc: cross) sc->AddCharge(p);
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
    if(sc->fWire1) iwires.insert(sc->fWire1);
    if(sc->fWire2) iwires.insert(sc->fWire2);

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

  if(alpha != 0){
    const double scip = sci->fPred;
    const double scjp = scj->fPred;

    // Self energy. SpaceCharges are never the same object
    ret -= alpha*sqr(scip + x);
    ret -= alpha*sqr(scjp - x);

    // Interaction. We're only seeing one end of the double-ended connection
    // here, so multiply by two.
    ret -= 2 * alpha * (scip + x) * sci->fNeiPotential;
    ret -= 2 * alpha * (scjp - x) * scj->fNeiPotential;

    // This miscounts if i and j are neighbours of each other
    for(const Neighbour& n: sci->fNeighbours){
      if(n.fSC == scj){
        // If we detect that case, remove the erroneous terms
        ret += 2 * alpha * (scip + x) * scjp * n.fCoupling;
        ret += 2 * alpha * (scjp - x) * scip * n.fCoupling;

        // And replace with the correct interaction terms
        ret -= 2 * alpha * (scip + x) * (scjp - x) * n.fCoupling;
        break;
      }
    }
  }


  const InductionWireHit* iwire1 = sci->fWire1;
  const InductionWireHit* jwire1 = scj->fWire1;

  const double qi1 = iwire1 ? iwire1->fCharge : 0;
  const double pi1 = iwire1 ? iwire1->fPred : 0;

  const double qj1 = jwire1 ? jwire1->fCharge : 0;
  const double pj1 = jwire1 ? jwire1->fPred : 0;

  if(iwire1 == jwire1){
    // Same wire means movement of charge cancels itself out
    if(iwire1) ret += Metric(qi1, pi1);
  }
  else{
    if(iwire1) ret += Metric(qi1, pi1 + x);
    if(jwire1) ret += Metric(qj1, pj1 - x);
  }

  const InductionWireHit* iwire2 = sci->fWire2;
  const InductionWireHit* jwire2 = scj->fWire2;

  const double qi2 = iwire2 ? iwire2->fCharge : 0;
  const double pi2 = iwire2 ? iwire2->fPred : 0;

  const double qj2 = jwire2 ? jwire2->fCharge : 0;
  const double pj2 = jwire2 ? jwire2->fPred : 0;

  if(iwire2 == jwire2){
    if(iwire2) ret += Metric(qi2, pi2);
  }
  else{
    if(iwire2) ret += Metric(qi2, pi2 + x);
    if(jwire2) ret += Metric(qj2, pj2 - x);
  }

  return ret;
}

// ---------------------------------------------------------------------------
QuadExpr Metric(const SpaceCharge* sc, double alpha)
{
  QuadExpr ret = 0;

  // How much charge is added to sc
  QuadExpr x = QuadExpr::X();

  if(alpha != 0){
    const double scp = sc->fPred;

    // Self energy
    ret -= alpha*sqr(scp + x);

    // Interaction. We're only seeing one end of the double-ended connection
    // here, so multiply by two.
    ret -= 2 * alpha * (scp + x) * sc->fNeiPotential;
  }

  // Prediction of the induction wires
  ret += Metric(sc->fWire1->fCharge, sc->fWire1->fPred + x);
  ret += Metric(sc->fWire2->fCharge, sc->fWire2->fPred + x);

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
    std::cout << "Solution at " << x << " is worse than current state! Scan from " << xmin << " to " << xmax << std::endl;
    for(double x = xmin; x < xmax; x += .01*(xmax-xmin)){
      std::cout << x << " " << chisq.Eval(x) << std::endl;
    }

    std::cout << "Soln, original, up edge, low edge:" << std::endl;
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

    for(unsigned int j = i+1; j < N; ++j){
      SpaceCharge* scj = cwire->fCrossings[j];

      const double x = SolvePair(cwire, sci, scj, alpha);

      if(x == 0) continue;

      // Actually make the update
      sci->AddCharge(+x);
      scj->AddCharge(-x);
    } // end for j
  } // end for i
}

// ---------------------------------------------------------------------------
void Iterate(SpaceCharge* sc, double alpha)
{
  const QuadExpr chisq = Metric(sc, alpha);

  // Find the minimum of a quadratic expression
  double x = -chisq.Linear()/(2*chisq.Quadratic());

  // Don't allow the SpaceCharge to go negative
  const double xmin = -sc->fPred;

  // Clamp to allowed range
  x = std::max(xmin, x);

  const double chisq_new = chisq.Eval(x);

  // Should try here too, because the function might be convex not concave, so
  // d/dx=0 gives the max not the min, and the true min is at one extreme of
  // the range.
  const double chisq_n = chisq.Eval(xmin);

  if(chisq_n < chisq_new)
    sc->AddCharge(xmin);
  else
    sc->AddCharge(x);
}

// ---------------------------------------------------------------------------
void Iterate(const std::vector<CollectionWireHit*>& cwires,
             const std::vector<SpaceCharge*>& orphanSCs,
             double alpha)
{
  // Visiting in a "random" order helps prevent local artefacts that are slow
  // to break up.
  unsigned int cwireIdx = 0;
  if (!cwires.empty()){
    do{
      Iterate(cwires[cwireIdx], alpha);

      const unsigned int prime = 1299827;
      cwireIdx = (cwireIdx+prime)%cwires.size();
    } while(cwireIdx != 0);
  }

  for(SpaceCharge* sc: orphanSCs) Iterate(sc, alpha);
}

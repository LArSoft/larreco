/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>

#include "TVirtualGeoTrack.h"
#include "larreco/Genfit/GFAbsRecoHit.h"
#include "larreco/Genfit/GFTrack.h"

genf::GFTrack::GFTrack(GFAbsTrackRep* defaultRep)
  : fTrackReps(NULL), fPDG(2112), fCardinal_rep(0), fNextHitToFit(0)
{
  addTrackRep(defaultRep);
}

genf::GFTrack::GFTrack() : fTrackReps(NULL), fCardinal_rep(0), fNextHitToFit(0)
{
  //trackReps = new TObjArray(defNumTrackReps);
}

genf::GFTrack::~GFTrack()
{
  if (fTrackReps != NULL) {
    for (unsigned int i = 0; i < getNumReps(); i++) {
      delete fTrackReps->At(i);
    }
    delete fTrackReps;
  }
  for (unsigned int i = 0; i < fHits.size(); i++) {
    delete fHits[i];
  }
  for (unsigned int i = 0; i < fBookkeeping.size(); ++i) {
    if (fBookkeeping.at(i) != NULL) delete fBookkeeping.at(i);
  }
}

genf::GFTrack::GFTrack(const GFTrack& _tr) : TObject()
{
  fCand = _tr.fCand;
  fCardinal_rep = _tr.fCardinal_rep;
  fNextHitToFit = _tr.fNextHitToFit;
  for (unsigned int i = 0; i < _tr.getNumHits(); i++) {
    fHits.push_back((_tr.getHit(i))->clone());
  }
  fTrackReps = NULL;
  for (unsigned int i = 0; i < _tr.getNumReps(); i++) {
    addTrackRep((_tr.getTrackRep(i))->clone());
  }
  for (unsigned int i = 0; i < fBookkeeping.size(); ++i)
    delete fBookkeeping[i];
  fBookkeeping.clear();

  for (unsigned int i = 0; i < _tr.fBookkeeping.size(); ++i) {
    fBookkeeping.push_back(new GFBookkeeping(*(_tr.fBookkeeping.at(i))));
  }
  fRepAtHit = _tr.fRepAtHit;
}

genf::GFTrack& genf::GFTrack::operator=(const GFTrack& _tr)
{
  if (fTrackReps != NULL) {
    for (unsigned int i = 0; i < getNumReps(); i++) {
      delete fTrackReps->At(i);
    }
    delete fTrackReps;
    fTrackReps = NULL;
  }
  for (unsigned int i = 0; i < fHits.size(); i++) {
    delete fHits[i];
  }
  for (unsigned int i = 0; i < fBookkeeping.size(); ++i) {
    if (fBookkeeping.at(i) != NULL) delete fBookkeeping.at(i);
  }

  for (unsigned int i = 0; i < _tr.getNumReps(); ++i) {
    addTrackRep(_tr.getTrackRep(i)->clone());
  }
  fCand = _tr.fCand;
  fCardinal_rep = _tr.fCardinal_rep;
  fNextHitToFit = _tr.fNextHitToFit;
  for (unsigned int i = 0; i < _tr.getNumHits(); i++) {
    fHits.push_back((_tr.getHit(i))->clone());
  }

  //clear the empty bookeeping objs made by addTrackRep and copy the others
  for (unsigned int i = 0; i < fBookkeeping.size(); ++i)
    delete fBookkeeping[i];
  fBookkeeping.clear();
  for (unsigned int i = 0; i < _tr.fBookkeeping.size(); ++i) {
    fBookkeeping.push_back(new GFBookkeeping(*(_tr.fBookkeeping.at(i))));
  }
  fRepAtHit = _tr.fRepAtHit;

  return *this;
}

void genf::GFTrack::reset()
{
  if (fTrackReps != NULL) {
    for (unsigned int i = 0; i < getNumReps(); i++) {
      delete fTrackReps->At(i);
      delete fBookkeeping.at(i);
    }
  }
  for (unsigned int i = 0; i < fHits.size(); i++) {
    delete fHits[i];
  }
  fHits.clear();
  fRepAtHit.clear();
}

void genf::GFTrack::mergeHits(GFTrack* trk)
{
  unsigned int nhits = trk->getNumHits();
  for (unsigned int i = 0; i < nhits; ++i) {
    unsigned int detId;
    unsigned int hitId;
    trk->getCand().getHit(i, detId, hitId);
    GFAbsRecoHit* hit = trk->getHit(i);
    addHit(hit, detId, hitId);
  }
  trk->fHits.clear();
}

void genf::GFTrack::setCandidate(const GFTrackCand& cand, bool doreset)
{
  fCand = cand;
  // reset fits
  if (doreset) {
    for (unsigned int i = 0; i < getNumReps(); i++) {
      ((GFAbsTrackRep*)fTrackReps->At(i))->reset();
    }
  }
}

void genf::GFTrack::fillGeoTrack(TVirtualGeoTrack* geotrk, unsigned int repid) const
{
  GFAbsTrackRep* rep = getTrackRep(repid);
  unsigned int n = fCand.getNHits();
  PrintROOTobject(std::cout, rep->getState());
  for (unsigned int i = 0; i < n; ++i) { // loop over hits
    GFDetPlane pl = fHits[i]->getDetPlane(rep);
    TVector3 pos = rep->getPos(pl);
#ifdef NDEBUG
    std::cout << pos.X() << "," << pos.Y() << "," << pos.Z() << std::endl;
#endif // NDEBUG
    geotrk->AddPoint(pos.X(), pos.Y(), pos.Z(), 0);
  } // end loop over hits
}

void genf::GFTrack::getResiduals(unsigned int detId, // which detector?
                                 unsigned int dim,   // which projection?
                                 unsigned int repid, // which trackrep ?
                                 std::vector<double>& result)
{
  unsigned int nhits = getNumHits();
  if (repid >= getNumReps()) return;
  GFAbsTrackRep* rep = getTrackRep(repid);      //->clone();
  for (unsigned int ih = 0; ih < nhits; ++ih) { // loop over hits
    unsigned int anid;
    unsigned int dummy;
    fCand.getHit(ih, anid, dummy); // check if this is a hit we want to look at
    if (anid == detId) {
      GFAbsRecoHit* hit = getHit(ih);
      // extrapolate trackrep there
      int repDim = rep->getDim();
      TMatrixT<Double_t> state(repDim, 1);
      GFDetPlane pl = hit->getDetPlane(rep);

      rep->extrapolate(pl, state);
      //rep->setState(state);
      //rep->setReferencePlane(pl);
      double res = hit->residualVector(rep, state, pl)[dim][0];

      //std::cout<<res<<std::endl;

      result.push_back(res);
    }
  }
}

void genf::GFTrack::printBookkeeping(std::ostream& out /* = std::cout */) const
{
  out << "GFTrack::printBookkeeping()" << std::endl;
  for (unsigned int i = 0; i < getNumReps(); ++i) {
    out << "trackRep " << i << ":" << std::endl;
    fBookkeeping.at(i)->Print(out);
  }
}

void genf::GFTrack::Print(std::ostream& out /* = std::cout */) const
{
  for (unsigned int i = 0; i < getNumReps(); ++i) {
    getTrackRep(i)->Print(out);
    fBookkeeping.at(i)->Print(out);
  }
  out << "GFTrack has " << getNumHits() << " detector hits." << std::endl;
}

void genf::GFTrack::getHitsByPlane(std::vector<std::vector<int>*>& retVal)
{
  for (int i = 0; retVal.size(); ++i) {
    delete retVal.at(i);
  }
  retVal.clear();
  //this method can only be called when all hits have been loaded
  if (fHits.size() != fCand.getNHits())
    throw GFException("genf::GFTrack::getResiduals(): inconsistent hits", __LINE__, __FILE__)
      .setFatal();
  if (fHits.size() < 2)
    throw GFException("genf::GFTrack::getResiduals(): less than 2 hits", __LINE__, __FILE__)
      .setFatal();
  unsigned int detId, hitId, planeId;
  fCand.getHitWithPlane(0, detId, hitId, planeId);
  //  std::cout << "$$$ " << 0 << " " << detId << " " << hitId << " " << planeId << std::endl;
  unsigned int lastPlane = planeId;
  retVal.push_back(new std::vector<int>);
  retVal.at(0)->push_back(0);
  for (unsigned int i = 1; i < fCand.getNHits(); ++i) {
    fCand.getHitWithPlane(i, detId, hitId, planeId);
    //std::cout << "$$$ " << i << " " << detId << " " << hitId << " " << planeId << std::endl;
    if (lastPlane == planeId) { retVal.at(retVal.size() - 1)->push_back(i); }
    else {
      lastPlane = planeId;
      retVal.push_back(new std::vector<int>);
      retVal.at(retVal.size() - 1)->push_back(i);
    }
  }
}

//ClassImp(GFTrack)

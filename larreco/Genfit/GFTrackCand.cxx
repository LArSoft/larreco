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
#include "larreco/Genfit/GFTrackCand.h"
#include "larreco/Genfit/GFException.h"

#include <algorithm>
#include <ostream>
#include <vector>

//ClassImp(GFTrackCand)

genf::GFTrackCand::GFTrackCand() : fCurv(0), fDip(0), fInv(false), fQoverpSeed(0.), fMcTrackId(-1)
{}

genf::GFTrackCand::~GFTrackCand() {}

genf::GFTrackCand::GFTrackCand(double curv,
                               double dip,
                               double inv,
                               std::vector<unsigned int> detIDs,
                               std::vector<unsigned int> hitIDs)
  : fDetId(detIDs)
  , fHitId(hitIDs)
  , fCurv(curv)
  , fDip(dip)
  , fInv(inv)
  , fQoverpSeed(0.)
  , fMcTrackId(-1)
{
  if (fDetId.size() != fHitId.size())
    throw GFException("genf::GFTrackCand::GFTrackCand(): hit/det size mismatch", __LINE__, __FILE__)
      .setFatal();
  fRho.resize(fDetId.size(), 0.);
}
genf::GFTrackCand::GFTrackCand(double curv,
                               double dip,
                               double inv,
                               std::vector<unsigned int> detIDs,
                               std::vector<unsigned int> hitIDs,
                               std::vector<double> rhos)
  : fDetId(detIDs)
  , fHitId(hitIDs)
  , fRho(rhos)
  , fCurv(curv)
  , fDip(dip)
  , fInv(inv)
  , fQoverpSeed(0.)
  , fMcTrackId(-1)
{
  if (fDetId.size() != fHitId.size())
    throw GFException("genf::GFTrackCand::GFTrackCand(): hit/det size mismatch", __LINE__, __FILE__)
      .setFatal();
  if (fDetId.size() != fHitId.size())
    throw GFException("genf::GFTrackCand::GFTrackCand(): rho/det size mismatch", __LINE__, __FILE__)
      .setFatal();
}

void genf::GFTrackCand::addHit(unsigned int detId,
                               unsigned int hitId,
                               double rho,
                               unsigned int planeId)
{
  fDetId.push_back(detId);
  fHitId.push_back(hitId);
  fPlaneId.push_back(planeId);
  fRho.push_back(rho);
}

std::vector<unsigned int> genf::GFTrackCand::GetHitIDs(int detId)
{
  if (detId < 0) { // return hits from all detectors
    return fHitId;
  }
  else {
    std::vector<unsigned int> result;
    unsigned int n = fHitId.size();
    for (unsigned int i = 0; i < n; ++i) {
      if (fDetId[i] == (unsigned int)detId) result.push_back(fHitId[i]);
    }
    return result;
  }
}

void genf::GFTrackCand::reset()
{
  fDetId.clear();
  fHitId.clear();
}

bool genf::GFTrackCand::HitInTrack(unsigned int detId, unsigned int hitId)
{
  for (unsigned int i = 0; i < fDetId.size(); i++) {
    if (detId == fDetId[i])
      if (hitId == fHitId[i]) return true;
  }
  return false;
}

bool genf::operator==(const GFTrackCand& lhs, const GFTrackCand& rhs)
{
  if (lhs.getNHits() != rhs.getNHits()) return false;
  bool result = std::equal(lhs.fDetId.begin(), lhs.fDetId.end(), rhs.fDetId.begin());
  result &= std::equal(lhs.fHitId.begin(), lhs.fHitId.end(), rhs.fHitId.begin());
  return result;
}

void genf::GFTrackCand::Print(std::ostream& out /* = std::cout */) const
{
  out << "======== GFTrackCand::print ========";
  if (fMcTrackId >= 0) out << "\nmcTrackId=" << fMcTrackId;
  out << "\nseed values for pos,direction, and q/p: " << std::endl;
  PrintROOTobject(out, fPosSeed);
  PrintROOTobject(out, fDirSeed);
  out << "q/p=" << fQoverpSeed << std::endl;
  if (fDetId.size() != fHitId.size())
    throw std::runtime_error("genf::GFTrackCand::GFTrackCand(): hit/det size mismatch");
  out << "detId|hitId|rho ";
  for (unsigned int i = 0; i < fDetId.size(); ++i)
    out << fDetId.at(i) << "|" << fHitId.at(i) << "|" << fRho.at(i) << " ";
  out << std::endl;
}

void genf::GFTrackCand::append(const GFTrackCand& rhs)
{
  unsigned int detId, hitId;
  double rho;
  for (unsigned int i = 0; i < rhs.getNHits(); ++i) {
    rhs.getHit(i, detId, hitId, rho);
    addHit(detId, hitId, rho);
  }
}

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
#include "Genfit/GFAbsRecoHit.h"

//ClassImp(GFAbsRecoHit)


genf::GFAbsRecoHit::GFAbsRecoHit(int NparHit) : fHitCoord(NparHit,1),
  fHitCov(NparHit,NparHit),
  fNparHit(NparHit)
{
}

genf::GFAbsRecoHit::GFAbsRecoHit() : fNparHit(-1) {
}

genf::GFAbsRecoHit::~GFAbsRecoHit(){;}

const std::string& genf::GFAbsRecoHit::getPolicyName(){
  std::cerr << "GFAbsRecoHit::getPolicyName() called for a reco hit, which wasnt derived from GFRecoHitIfc -> abort" << std::endl;
  throw;
}

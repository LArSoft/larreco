/* Copyright 2008-2009, Technische Universitaet Muenchen,
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
#include "GFGeoMatManager.h"

#include "math.h"
#include <cassert>

#include "TGeoManager.h"
#include "TGeoMaterial.h"

#include "GFException.h"

float MeanExcEnergy_get(int Z);
float MeanExcEnergy_get(TGeoMaterial*);

void genf::GFGeoMatManager::getMaterialParameters(double& matDensity,
                                                  double& matZ,
                                                  double& matA,
                                                  double& radiationLength,
                                                  double& mEE)
{
  if (!gGeoManager->GetCurrentVolume()->GetMedium())
    throw GFException(
      "genf::GFGeoMatManager::getMaterialParameters(): no medium in volume!", __LINE__, __FILE__)
      .setFatal();
  TGeoMaterial* mat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();
  //std::cout << "GFGeoMatManager::getMaterialParameters: CurrentVolume is " << std::endl;
  //gGeoManager->GetCurrentVolume()->Print();
  //std::cout << "GFGeoMatManager::getMaterialParameters: Material is " << std::endl;
  //gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial()->Print();

  matDensity = mat->GetDensity();
  matZ = mat->GetZ();
  matA = mat->GetA();
  radiationLength = mat->GetRadLen();
  mEE = MeanExcEnergy_get(mat);
}

void genf::GFGeoMatManager::initTrack(const double& posx,
                                      const double& posy,
                                      const double& posz,
                                      const double& dirx,
                                      const double& diry,
                                      const double& dirz)
{
  gGeoManager->InitTrack(posx, posy, posz, dirx, diry, dirz);
}

double genf::GFGeoMatManager::stepOrNextBoundary(const double& maxStep)
{
  gGeoManager->FindNextBoundaryAndStep(maxStep);
  return gGeoManager->GetStep();
}

//ClassImp(GFGeoMatManager)

/*
Reference for elemental mean excitation energies at:
http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
*/

const int MeanExcEnergy_NELEMENTS = 92;
const float MeanExcEnergy_vals[MeanExcEnergy_NELEMENTS] = {
  19.2,  41.8,  40.0,  63.7,  76.0,  78.,   82.0,  95.0,  115.0, 137.0, 149.0, 156.0, 166.0, 173.0,
  173.0, 180.0, 174.0, 188.0, 190.0, 191.0, 216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0,
  322.0, 330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0, 417.0, 424.0,
  428.0, 441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0, 488.0, 491.0,
  501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0, 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0,
  694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0, 823.0, 823.0, 830.0,
  825.0, 794.0, 827.0, 826.0, 841.0, 847.0, 878.0, 890.0};

float MeanExcEnergy_get(int Z)
{
  if ((Z <= 0) || (Z > MeanExcEnergy_NELEMENTS))
    throw GFException("MeanExcEnergy_get(): Z out of range", __LINE__, __FILE__).setFatal();
  return MeanExcEnergy_vals[Z - 1];
}

float MeanExcEnergy_get(TGeoMaterial* mat)
{
  if (mat->IsMixture()) {
    double logMEE = 0.;
    double denom = 0.;
    TGeoMixture* mix = (TGeoMixture*)mat;
    for (int i = 0; i < mix->GetNelements(); ++i) {
      int index = int(floor((mix->GetZmixt())[i]));
      //check whether the floor command worked
      assert(fabs(index - ((mix->GetZmixt())[i])) < 1.e-3);
      logMEE += 1. / (mix->GetAmixt())[i] * (mix->GetWmixt())[i] * (mix->GetZmixt())[i] *
                log(MeanExcEnergy_get(index));
      denom += (mix->GetWmixt())[i] * (mix->GetZmixt())[i] * 1. / (mix->GetAmixt())[i];
    }
    logMEE /= denom;
    return exp(logMEE);
  }
  else { // not a mixture
    int index = int(floor(mat->GetZ()));
    //check whether the floor command worked
    assert(fabs(index - mat->GetZ()) < 1.e-3);
    return MeanExcEnergy_get(index);
  }
}

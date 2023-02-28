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
#include "larreco/Genfit/GFAbsTrackRep.h"
#include "larreco/Genfit/GFException.h"
#include <iostream>

genf::GFAbsTrackRep::GFAbsTrackRep()
  : fDimension(5)
  , fState(5, 1)
  , fCov(5, 5)
  , fChiSqu(0)
  , fNdf(0)
  , fStatusFlag(0)
  , fInverted(false)
  , fFirstState(5, 1)
  , fFirstCov(5, 5)
  , fLastState(5, 1)
  , fLastCov(5, 5)
{}

genf::GFAbsTrackRep::GFAbsTrackRep(int dim)
  : fDimension(dim)
  , fState(dim, 1)
  , fCov(dim, dim)
  , fChiSqu(0)
  , fNdf(0)
  , fStatusFlag(0)
  , fInverted(false)
  , fFirstState(dim, 1)
  , fFirstCov(dim, dim)
  , fLastState(dim, 1)
  , fLastCov(dim, dim)
{}

genf::GFAbsTrackRep::~GFAbsTrackRep() {}

double genf::GFAbsTrackRep::extrapolate(const GFDetPlane& plane)
{
  TMatrixT<Double_t> statePred(fDimension, 1);
  TMatrixT<Double_t> covPred(fDimension, fDimension);
  double retVal = extrapolate(plane, statePred, covPred);
  setData(statePred, plane, &covPred);
  return retVal;
}

//default implentation might be overwritten, please see the doxy docu
double genf::GFAbsTrackRep::extrapolate(const GFDetPlane& plane, TMatrixT<Double_t>& statePred)
{
  TMatrixT<Double_t> cov(fDimension, fDimension);
  return extrapolate(plane, statePred, cov);
}

void genf::GFAbsTrackRep::Abort(std::string method)
{
  std::cerr << method << " as implemented in " << __FILE__
            << " was called. This means that this feature was used "
            << "in a track rep which didnt overwrite this method. " << std::endl
            << "C++ throw;" << std::endl;
  //system call abort
  throw GFException("genf::GFAbsTrackRep: " + method + "() not implemented", __LINE__, __FILE__)
    .setFatal();
}

void genf::GFAbsTrackRep::extrapolateToPoint(const TVector3& /* point */,
                                             TVector3& /* poca */,
                                             TVector3& /* normVec */)
{
  Abort("extrapolateToPoca()");
}

void genf::GFAbsTrackRep::extrapolateToLine(const TVector3& /* point1 */,
                                            const TVector3& /* point2 */,
                                            TVector3& /* poca */,
                                            TVector3& /* normVec */,
                                            TVector3& /* poca_onwire */)
{
  Abort("extrapolateToLine()");
}

void genf::GFAbsTrackRep::stepalong(double /* h */)
{
  Abort("stepalong()");
}

void genf::GFAbsTrackRep::getPosMomCov(const GFDetPlane& /* pl */,
                                       TVector3& /* pos */,
                                       TVector3& /* mom */,
                                       TMatrixT<Double_t>& /* cov */)
{
  Abort("getPosMomCov()");
}

void genf::GFAbsTrackRep::reset()
{
  std::cout << "GFAbsTrackRep::reset" << std::endl;
  TVector3 nullVec(0., 0., 0.);
  fRefPlane.set(nullVec, nullVec, nullVec);
  fState.Zero();
  fCov.Zero();
  fFirstState.Zero();
  fFirstCov.Zero();
  fLastState.Zero();
  fLastCov.Zero();
}

void genf::GFAbsTrackRep::Print(std::ostream& out /* = std::cout */) const
{
  out << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  out << "GFAbsTrackRep::Parameters at reference plane ";
  fRefPlane.Print(out);
  out << "GFAbsTrackRep::State" << std::endl;
  PrintROOTmatrix(out, fState);
  out << "GFAbsTrackRep::Covariances" << std::endl;
  PrintROOTmatrix(out, fCov);
  out << "GFAbsTrackRep::chi^2" << std::endl;
  out << fChiSqu << std::endl;
  out << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
}

//ClassImp(GFAbsTrackRep)

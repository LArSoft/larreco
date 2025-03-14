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

/*
  The DAF algorithm stems from the following sources:
  R. Fruehwirth & A. Strandlie,
  Computer Physics Communications 120 (199) 197-214
  &
  CERN thesis: Dissertation by Matthias Winkler
 */

#include "larreco/Genfit/GFDaf.h"
#include "larreco/Genfit/GFBookkeeping.h"
#include "larreco/Genfit/GFDetPlane.h"
#include "larreco/Genfit/GFException.h"

#include "TMath.h"
#include <math.h>

#include "larreco/Genfit/GFAbsRecoHit.h"
#include "larreco/Genfit/GFAbsTrackRep.h"
#include "larreco/Genfit/GFTrack.h"

#define COVEXC "cov_is_zero"

genf::GFDaf::GFDaf() : fBlowUpFactor(500.)
{
  setProbCut(0.01);
  setBetas(81, 8, 4, 1, 1, 1);
}

genf::GFDaf::~GFDaf()
{
  ;
}

void genf::GFDaf::processTrack(GFTrack* trk)
{
  trk->clearBookkeeping();

  unsigned int nreps = trk->getNumReps();
  unsigned int nhits = trk->getNumHits();

  for (unsigned int i = 0; i < nreps; ++i) {
    GFBookkeeping* bk = trk->getBK(i);
    bk->setNhits(nhits);
    bk->bookMatrices("fPredSt");
    bk->bookMatrices("fPredCov");
    bk->bookMatrices("bPredSt");
    bk->bookMatrices("bPredCov");
    bk->bookMatrices("H");
    bk->bookMatrices("m");
    bk->bookMatrices("V");
    bk->bookGFDetPlanes("pl");
    bk->bookMatrices("smooResid");
    bk->bookNumbers("p", 1.);
  }

  //plane grouping
  std::vector<std::vector<int>*> planes;
  trk->getHitsByPlane(planes);
  int nPlanes = planes.size();

  trk->clearFailedHits();
  for (unsigned int iBeta = 0; iBeta < fBeta.size(); iBeta++) {
    //std::cout << "@@ beta " << fBeta.at(iBeta) << std::endl;
    if (iBeta > 0) blowUpCovs(trk);
    for (int ipl = 0; ipl < nPlanes; ++ipl) {
      //std::cout << "@@ plane " << ipl << std::endl;
      std::vector<GFAbsRecoHit*> hits;
      unsigned int nhitsInPlane = planes.at(ipl)->size();
      for (unsigned int ihitInPl = 0; ihitInPl < nhitsInPlane; ++ihitInPl) {
        hits.push_back(trk->getHit(planes.at(ipl)->at(ihitInPl)));
      }
      //now hits contains pointers to all hits in the next plane
      //for non-planar detectors this will always be just one hit
      for (unsigned int irep = 0; irep < nreps; ++irep) {
        GFAbsTrackRep* rep = trk->getTrackRep(irep);
        if (rep->getStatusFlag() != 0) continue;
        //std::cout << "@@ rep " << irep << std::endl;
        GFBookkeeping* bk = trk->getBK(irep);
        if (bk->hitFailed(planes.at(ipl)->at(0)) > 0) continue;
        GFDetPlane pl;
        // get prototypes for matrices
        int repDim = rep->getDim();
        TMatrixT<Double_t> state(repDim, 1);
        TMatrixT<Double_t> cov(repDim, repDim);
        ;

        if (ipl > 0 || (ipl == 0 && iBeta == 0)) {
          // get the (virtual) detector plane
          //std::cout << "forward ### iBeta, ipl, irep " << iBeta<< " " << ipl << " " << irep << std::endl;
          //	  rep->getState().Print();
          try {
            pl = hits.at(0)->getDetPlane(rep);
            //do the extrapolation
            rep->extrapolate(pl, state, cov);

            if (cov[0][0] < 1.E-50) {
              GFException exc(COVEXC, __LINE__, __FILE__);
              throw exc;
            }
          }
          catch (GFException& exc) {
            std::cerr << exc.what();
            exc.info();
            for (unsigned int i = 0; i < planes.at(ipl)->size(); ++i)
              trk->addFailedHit(irep, planes.at(ipl)->at(i));
            if (exc.isFatal()) {
              rep->setStatusFlag(1);
              //abort();
            }
            continue; //go to next rep immediately
          }
        }
        else {
          pl = rep->getReferencePlane();
          state = rep->getState();
          cov = rep->getCov();
        }
        //the H matrix has to be identical for all hits in the plane
        TMatrixT<Double_t> H = hits.at(0)->getHMatrix(rep);
        bk->setMatrix("fPredSt", planes.at(ipl)->at(0), state);
        bk->setMatrix("fPredCov", planes.at(ipl)->at(0), cov);

        TMatrixT<Double_t> covInv;
        invertMatrix(cov, covInv);

        TMatrixT<Double_t> Htransp(H);
        Htransp.T();

        bk->setMatrix("H", planes.at(ipl)->at(0), H);
        bk->setDetPlane("pl", planes.at(ipl)->at(0), pl);

        TMatrixT<Double_t> stMod(state.GetNrows(), 1);

        double sumPk(0);
        for (unsigned int ihit = 0; ihit < nhitsInPlane; ++ihit) {
          double pki;
          bk->getNumber("p", planes.at(ipl)->at(ihit), pki);
          sumPk += pki;
        }
        TMatrixT<Double_t> V = hits.at(0)->getHitCov(pl);
        TMatrixT<Double_t> Vinv;
        invertMatrix(V, Vinv);
        bk->setMatrix("V", planes.at(ipl)->at(0), V);
        for (unsigned int ihit = 0; ihit < nhitsInPlane; ++ihit) {
          //std::cout << "%%%% forward hit " << planes.at(ipl)->at(ihit) << std::endl;

          TMatrixT<Double_t> m = hits.at(ihit)->getHitCoord(pl);
          bk->setMatrix("m", planes.at(ipl)->at(ihit), m);

          double pki;
          bk->getNumber("p", planes.at(ipl)->at(ihit), pki);
          TMatrixT<Double_t> Gain = calcGain(cov, V, H, sumPk);
          //std::cout << "using weight " << pki << std::endl;
          stMod += pki * Gain * (m - H * state);
        }
        state += stMod;
        invertMatrix(covInv + sumPk * Htransp * Vinv * H, cov);
        rep->setData(state, pl, &cov);

      } //loop over reps
    }   //loop over planes for forward fit
    blowUpCovs(trk);
    //now do the loop over the planes in backward direction
    for (int ipl = nPlanes - 1; ipl >= 0; --ipl) {
      //std::cout << "@bw@ plane " << ipl << std::endl;
      std::vector<GFAbsRecoHit*> hits;
      unsigned int nhitsInPlane = planes.at(ipl)->size();
      for (unsigned int ihitInPl = 0; ihitInPl < nhitsInPlane; ++ihitInPl) {
        hits.push_back(trk->getHit(planes.at(ipl)->at(ihitInPl)));
      }
      //now hits contains pointers to all hits in the next plane
      //for non-planar detectors this will always be just one hit
      for (unsigned int irep = 0; irep < nreps; ++irep) {
        GFAbsTrackRep* rep = trk->getTrackRep(irep);
        if (rep->getStatusFlag() != 0) continue;
        //std::cout << "@bw@ rep " << irep << std::endl;
        GFBookkeeping* bk = trk->getBK(irep);
        if (bk->hitFailed(planes.at(ipl)->at(0)) > 0) continue;
        // get prototypes for matrices
        int repDim = rep->getDim();
        TMatrixT<Double_t> state(repDim, 1);
        TMatrixT<Double_t> cov(repDim, repDim);
        ;

        TMatrixT<Double_t> H;
        bk->getMatrix("H", planes.at(ipl)->at(0), H);
        TMatrixT<Double_t> Htransp(H);
        Htransp.T();
        GFDetPlane pl;
        bk->getDetPlane("pl", planes.at(ipl)->at(0), pl);
        //std::cout << "backward ### iBeta, ipl, irep " << iBeta<< " " << ipl << " " << irep << std::endl;
        if (ipl < (nPlanes - 1)) {
          try {
            //do the extrapolation
            rep->extrapolate(pl, state, cov);
            if (cov[0][0] < 1.E-50) {
              GFException exc(COVEXC, __LINE__, __FILE__);
              throw exc;
            }
          }
          catch (GFException& exc) {
            //std::cout << __FILE__ << " " << __LINE__ << std::endl;
            std::cerr << exc.what();
            exc.info();
            //TAG
            for (unsigned int i = 0; i < planes.at(ipl)->size(); ++i)
              trk->addFailedHit(irep, planes.at(ipl)->at(i));
            if (exc.isFatal()) { rep->setStatusFlag(1); }
            continue; //go to next rep immediately
          }
        }
        else {
          state = rep->getState();
          cov = rep->getCov();
        }
        TMatrixT<Double_t> covInv;
        invertMatrix(cov, covInv);

        bk->setMatrix("bPredSt", planes.at(ipl)->at(0), state);
        bk->setMatrix("bPredCov", planes.at(ipl)->at(0), cov);

        TMatrixT<Double_t> stMod(state.GetNrows(), 1);

        double sumPk(0);
        for (unsigned int ihit = 0; ihit < nhitsInPlane; ++ihit) {
          double pki;
          bk->getNumber("p", planes.at(ipl)->at(ihit), pki);
          sumPk += pki;
        }

        TMatrixT<Double_t> V;
        bk->getMatrix("V", planes.at(ipl)->at(0), V);
        TMatrixT<Double_t> Vinv;
        invertMatrix(V, Vinv);

        for (unsigned int ihit = 0; ihit < nhitsInPlane; ++ihit) {
          //std::cout << "%%bw%% hit " << planes.at(ipl)->at(ihit) << std::endl;
          TMatrixT<Double_t> m;
          bk->getMatrix("m", planes.at(ipl)->at(ihit), m);

          double pki;
          bk->getNumber("p", planes.at(ipl)->at(ihit), pki);
          TMatrixT<Double_t> Gain = calcGain(cov, V, H, sumPk);
          //std::cout << "using pki " << pki << std::endl;

          stMod += pki * Gain * (m - H * state);
        }
        state += stMod;
        invertMatrix(covInv + sumPk * Htransp * Vinv * H, cov);

        rep->setData(state, pl, &cov);
      } //done with loop over reps
    }   //loop over planes for backward
    //calculate smoothed states and covs
    TMatrixT<Double_t> fSt, fCov, bSt, bCov, smooSt, smooCov, fCovInv, bCovInv, m, H, smooResid;
    for (int ipl = 0; ipl < nPlanes; ++ipl) {
      for (unsigned int irep = 0; irep < nreps; ++irep) {
        if (trk->getTrackRep(irep)->getStatusFlag() != 0) continue;
        GFBookkeeping* bk = trk->getBK(irep);
        if (bk->hitFailed(planes.at(ipl)->at(0)) > 0) continue;
        bk->getMatrix("fPredSt", planes.at(ipl)->at(0), fSt);
        bk->getMatrix("bPredSt", planes.at(ipl)->at(0), bSt);
        bk->getMatrix("fPredCov", planes.at(ipl)->at(0), fCov);
        bk->getMatrix("bPredCov", planes.at(ipl)->at(0), bCov);
        fCovInv.ResizeTo(fCov);
        bCovInv.ResizeTo(bCov);
        invertMatrix(fCov, fCovInv);
        invertMatrix(bCov, bCovInv);
        smooCov.ResizeTo(fCov);
        invertMatrix(fCovInv + bCovInv, smooCov);
        smooSt.ResizeTo(fSt.GetNrows(), 1);
        smooSt = smooCov * (fCovInv * fSt + bCovInv * bSt);
        //std::cout << "#@#@#@#@#@# smooResid ipl, irep " << ipl << " " << irep << std::endl;
        bk->getMatrix("H", planes.at(ipl)->at(0), H);
        for (unsigned int ihit = 0; ihit < planes.at(ipl)->size(); ++ihit) {
          //TAG
          bk->getMatrix("m", planes.at(ipl)->at(ihit), m);
          smooResid.ResizeTo(m.GetNrows(), 1);
          smooResid = m - H * smooSt;
          //std::cout << "%%%% hit " << planes.at(ipl)->at(ihit) << std::endl;
          //smooResid.Print();
          bk->setMatrix("smooResid", planes.at(ipl)->at(ihit), smooResid);
        }
      }
    } //end loop over planes for smoothed state calculation

    //calculate the new weights
    if (iBeta != fBeta.size() - 1) { //dont do it for the last beta, ause fitting will stop
      for (unsigned int irep = 0; irep < nreps; ++irep) {
        GFBookkeeping* bk = trk->getBK(irep);
        for (int ipl = 0; ipl < nPlanes; ++ipl) {
          if (trk->getTrackRep(irep)->getStatusFlag() != 0) continue;
          if (bk->hitFailed(planes.at(ipl)->at(0)) > 0) continue;
          double sumPhi = 0.;
          double sumPhiCut = 0.;
          std::vector<double> phi;
          TMatrixT<Double_t> V;
          bk->getMatrix("V", planes.at(ipl)->at(0), V);
          V *= fBeta.at(iBeta);
          TMatrixT<Double_t> Vinv;
          invertMatrix(V, Vinv);
          for (unsigned int ihit = 0; ihit < planes.at(ipl)->size(); ++ihit) {
            //TMatrixT<Double_t> smooResid;
            bk->getMatrix("smooResid", planes.at(ipl)->at(ihit), smooResid);
            TMatrixT<Double_t> smooResidTransp(smooResid);
            smooResidTransp.T();
            int dimV = V.GetNrows();
            double detV = V.Determinant();
            TMatrixT<Double_t> expArg = smooResidTransp * Vinv * smooResid;
            if (expArg.GetNrows() != 1)
              throw std::runtime_error("genf::GFDaf::processTrack(): wrong matrix dimensions!");
            double thisPhi =
              1. / (pow(2. * TMath::Pi(), dimV / 2.) * sqrt(detV)) * exp(-0.5 * expArg[0][0]);
            phi.push_back(thisPhi);
            sumPhi += thisPhi;

            double cutVal = chi2Cuts[dimV];
            if (cutVal <= 1.E-6)
              throw std::runtime_error("genf::GFDaf::processTrack(): chi2 cut too low");
            sumPhiCut += 1. / (2 * TMath::Pi() * sqrt(detV)) * exp(-0.5 * cutVal / fBeta.at(iBeta));
          }
          for (unsigned int ihit = 0; ihit < planes.at(ipl)->size(); ++ihit) {
            bk->setNumber("p", planes.at(ipl)->at(ihit), phi.at(ihit) / (sumPhi + sumPhiCut));
          }
        }
      }
    }

  } //loop over betas

  return;
}

void genf::GFDaf::blowUpCovs(GFTrack* trk)
{
  int nreps = trk->getNumReps();
  for (int irep = 0; irep < nreps; ++irep) {
    GFAbsTrackRep* arep = trk->getTrackRep(irep);
    //dont do it for already compromsied reps, since they wont be fitted anyway
    if (arep->getStatusFlag() == 0) {
      TMatrixT<Double_t> cov = arep->getCov();
      for (int i = 0; i < cov.GetNrows(); ++i) {
        for (int j = 0; j < cov.GetNcols(); ++j) {
          if (i != j) { //off diagonal
            cov[i][j] = 0.;
          }
          else { //diagonal
            cov[i][j] = cov[i][j] * fBlowUpFactor;
          }
        }
      }
      arep->setCov(cov);
    }
  }
}

void genf::GFDaf::invertMatrix(const TMatrixT<Double_t>& mat, TMatrixT<Double_t>& inv)
{
  inv.ResizeTo(mat);
  inv = (mat);
  double det = 0;
  inv.Invert(&det);
  if (TMath::IsNaN(det)) {
    GFException e("Daf Gain: det of matrix is nan", __LINE__, __FILE__);
    e.setFatal();
    throw e;
  }
  if (det == 0) {
    GFException e("cannot invert matrix in Daf - det=0", __LINE__, __FILE__);
    e.setFatal();
    throw e;
  }
}

TMatrixT<Double_t> genf::GFDaf::calcGain(const TMatrixT<Double_t>& C,
                                         const TMatrixT<Double_t>& V,
                                         const TMatrixT<Double_t>& H,
                                         const double& p)
{

  //get C^-1
  TMatrixT<Double_t> Cinv;
  invertMatrix(C, Cinv);
  TMatrixT<Double_t> Vinv;
  invertMatrix(V, Vinv);
  //get H^T
  TMatrixT<Double_t> Htransp(H);
  Htransp.T();

  TMatrixT<Double_t> covsum = Cinv + (p * Htransp * Vinv * H);
  TMatrixT<Double_t> covsumInv;
  invertMatrix(covsum, covsumInv);

  return (covsumInv * Htransp * Vinv);
}

void genf::GFDaf::setProbCut(double val)
{

  if (fabs(val - 0.01) < 1.E-10) {
    chi2Cuts[1] = 6.63;
    chi2Cuts[2] = 9.21;
    chi2Cuts[3] = 11.34;
    chi2Cuts[4] = 13.23;
    chi2Cuts[5] = 15.09;
  }
  else if (fabs(val - 0.005) < 1.E-10) {
    chi2Cuts[1] = 7.88;
    chi2Cuts[2] = 10.60;
    chi2Cuts[3] = 12.84;
    chi2Cuts[4] = 14.86;
    chi2Cuts[5] = 16.75;
  }
  else if (fabs(val - 0.001) < 1.E-10) {
    chi2Cuts[1] = 10.83;
    chi2Cuts[2] = 13.82;
    chi2Cuts[3] = 16.27;
    chi2Cuts[4] = 18.47;
    chi2Cuts[5] = 20.51;
  }
  else {
    throw GFException(
      "genf::GFTrackCand::GFDafsetProbCut(): value not supported", __LINE__, __FILE__)
      .setFatal()
      .setNumbers("value", {val});
  }
}

void genf::GFDaf::setBetas(double b1,
                           double b2,
                           double b3 /* = -1. */,
                           double b4 /* = -1. */,
                           double b5 /* = -1. */,
                           double b6 /* = -1. */,
                           double b7 /* = -1. */,
                           double b8 /* = -1. */,
                           double b9 /* = -1. */,
                           double b10 /* = -1. */
)
{

  /* // asserting version
  assert(b1>0);fBeta.push_back(b1);
  assert(b2>0 && b2<=b1);fBeta.push_back(b2);
  if(b3>=0.) {
    assert(b3<=b2);fBeta.push_back(b3);
    if(b4>=0.) {
      assert(b4<=b3);fBeta.push_back(b4);
      if(b5>=0.) {
        assert(b5<=b4);fBeta.push_back(b5);
        if(b6>=0.) {
          assert(b6<=b5);fBeta.push_back(b6);
          if(b7>=0.) {
            assert(b7<=b6);fBeta.push_back(b7);
            if(b8>=0.) {
              assert(b8<=b7);fBeta.push_back(b8);
              if(b9>=0.) {
                assert(b9<=b8);fBeta.push_back(b9);
                if(b10>=0.) {
                  assert(b10<=b9);fBeta.push_back(b10);
                }
              }
            }
          }
        }
      }
    }
  }
  */
  if (b1 <= 0) throw std::runtime_error("genf::GFDaf::setBetas(): b1 <= 0");
  fBeta.push_back(b1);

  if (b2 <= 0 || b2 > b1) throw std::runtime_error("genf::GFDaf::setBetas(): b2 in wrong range");
  fBeta.push_back(b2);

  if (b3 < 0.) return;
  if (b3 > b2) throw std::runtime_error("genf::GFDaf::setBetas(): b3 > b2");
  fBeta.push_back(b3);

  if (b4 < 0.) return;
  if (b4 > b3) throw std::runtime_error("genf::GFDaf::setBetas(): b4 > b3");
  fBeta.push_back(b4);

  if (b5 < 0.) return;
  if (b5 > b4) throw std::runtime_error("genf::GFDaf::setBetas(): b5 > b4");
  fBeta.push_back(b5);

  if (b6 < 0.) return;
  if (b6 > b5) throw std::runtime_error("genf::GFDaf::setBetas(): b6 > b5");
  fBeta.push_back(b6);

  if (b7 < 0.) return;
  if (b7 > b6) throw std::runtime_error("genf::GFDaf::setBetas(): b7 > b6");
  fBeta.push_back(b7);

  if (b8 < 0.) return;
  if (b8 > b7) throw std::runtime_error("genf::GFDaf::setBetas(): b8 > b7");
  fBeta.push_back(b8);

  if (b9 < 0.) return;
  if (b9 > b8) throw std::runtime_error("genf::GFDaf::setBetas(): b9 > b8");
  fBeta.push_back(b9);

  if (b10 < 0.) return;
  if (b10 > b9) throw std::runtime_error("genf::GFDaf::setBetas(): b10 > b9");
  fBeta.push_back(b10);
}

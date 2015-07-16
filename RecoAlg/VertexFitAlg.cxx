//////////////////////////////////////////////////////////////////////
///
/// VertexFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 3D vertex given a set of track hits
///
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>
#include <iostream>
#include <iomanip>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RecoAlg/VertexFitAlg.h"


namespace trkf{

  VertexFitMinuitStruct VertexFitAlg::fVtxFitMinStr;

  /////////////////////////////////////////
  void VertexFitAlg::fcnVtxPos(Int_t &, Double_t *, Double_t &fval, double *par, Int_t flag)
  {
    // Minuit function for fitting the vertex position and vertex track directions
    
    fval = 0;
    double vWire, DirX, DirY, DirZ, DirU, dX, dU, arg;
    unsigned short ipl, lastpl, indx;
    
    for(unsigned short itk = 0; itk < fVtxFitMinStr.HitX.size(); ++itk) {
      lastpl = 4;
      // index of the track Y direction vector. Z direction is the next one
      indx = 3 + 2 * itk;
      for(unsigned short iht = 0; iht < fVtxFitMinStr.HitX[itk].size(); ++iht) {
        ipl = fVtxFitMinStr.Plane[itk][iht];
        if(ipl != lastpl) {
          // get the vertex position in this plane
          // vertex wire number in the Detector coordinate system (equivalent to WireCoordinate)
          //vtx wir = vtx Y  * OrthY                + vtx Z  * OrthZ                    - wire offset
          vWire = par[1] * fVtxFitMinStr.OrthY[ipl] + par[2] * fVtxFitMinStr.OrthZ[ipl] - fVtxFitMinStr.FirstWire[ipl];
//          if(flag == 1) mf::LogVerbatim("VF")<<"fcn vtx "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" vWire "<<vWire<<" OrthY "<<fVtxFitMinStr.OrthY[ipl]<<" OrthZ "<<fVtxFitMinStr.OrthZ[ipl];
          lastpl = ipl;
        } // ipl != lastpl
        DirY = par[indx];
        DirZ = par[indx + 1];
        // rotate the track direction DirY, DirZ into the wire coordinate of this plane. The OrthVectors in ChannelMapStandardAlg
        // are divided by the wire pitch so we need to correct for that here
        DirU = fVtxFitMinStr.WirePitch * (DirY * fVtxFitMinStr.OrthY[ipl] + DirZ * fVtxFitMinStr.OrthZ[ipl]);
        // distance (cm) between the wire and the vertex in the wire coordinate system (U)
        dU = fVtxFitMinStr.WirePitch * (fVtxFitMinStr.Wire[itk][iht] - vWire);
        if(std::abs(DirU) < 1E-3 || std::abs(dU) < 1E-3) {
          // vertex is on the wire
          dX = par[0] - fVtxFitMinStr.HitX[itk][iht];
        } else {
          // project from vertex to the wire. We need to find dX/dU so first find DirX
          DirX = 1 - DirY * DirY - DirZ * DirZ;
          // DirX should be > 0 but the bounds on DirY and DirZ are +/- 1 so it is possible for a non-physical result.
          if(DirX < 0) DirX = 0;
          DirX = sqrt(DirX);
          // Get the DirX sign from the relative X position of the hit and the vertex
          if(fVtxFitMinStr.HitX[itk][iht] < par[0]) DirX = -DirX;
          dX = par[0] + (dU * DirX / DirU) - fVtxFitMinStr.HitX[itk][iht];
        }
        arg = dX / fVtxFitMinStr.HitXErr[itk][iht];
//        if(flag == 1) mf::LogVerbatim("VF")<<"fcn itk "<<itk<<" iht "<<iht<<" ipl "<<ipl<<" DirX "<<DirX<<" DirY "<<DirY<<" DirZ "<<DirZ
//        <<" DirU "<<DirU<<" W "<<fVtxFitMinStr.Wire[itk][iht]<<" X "<<fVtxFitMinStr.HitX[itk][iht]<<" dU "<<dU<<" dX "<<dX<<" arg "<<arg;
        fval += arg * arg;
      } // iht
    } //itk
    
    fval /= fVtxFitMinStr.DoF;
    
    // save the final chisq/dof in the struct on the last call
    if(flag == 3) fVtxFitMinStr.ChiDoF = fval;
    
  } // fcnVtxPos

  /////////////////////////////////////////
  
  VertexFitAlg::VertexFitAlg() { }

  VertexFitAlg::~VertexFitAlg() { }


  void VertexFitAlg::VertexFit(std::vector<std::vector<geo::WireID>>& hitWID,
                               std::vector<std::vector<double>>& hitX,
                               std::vector<std::vector<double>>& hitXErr,
                               TVector3& VtxPos, TVector3& VtxPosErr,
                               std::vector<TVector3>& TrkDir, std::vector<TVector3>& TrkDirErr,
                               float& ChiDOF)
  {
    // The passed set of hit WireIDs, X positions and X errors associated with a Track
    // are fitted to a vertex position VtxPos. The fitted track direction vectors trkDir, TrkDirErr
    // and ChiDOF are returned to the calling routine

    // assume failure
    ChiDOF = 9999;
    
    // need at least hits for two tracks
    if(hitX.size() < 2) return;
    if(hitX.size() != hitWID.size()) return;
    if(hitX.size() != hitXErr.size()) return;
    if(hitX.size() != TrkDir.size()) return;
    
    // number of variables = 3 for the vertex position + 2 * number of track directions
    const unsigned int ntrks = hitX.size();
    const unsigned int npars = 3 + 2 * ntrks;
    unsigned int npts = 0, itk;
    for(itk = 0; itk < ntrks; ++itk) npts += hitX[itk].size();
    
    if(npts < ntrks) return;
    
    // Get the cryostat and tpc from the first hit
    unsigned int cstat, tpc, nplanes, ipl, iht;
    cstat = hitWID[0][0].Cryostat;
    tpc = hitWID[0][0].TPC;
    nplanes = geom->Cryostat(cstat).TPC(tpc).Nplanes();
    
    fVtxFitMinStr.Cstat = cstat;
    fVtxFitMinStr.TPC = tpc;
    fVtxFitMinStr.NPlanes = nplanes;
    fVtxFitMinStr.WirePitch = geom->WirePitch(hitWID[0][0]);

    // Put geometry conversion factors into the struct
    for(ipl = 0; ipl < nplanes; ++ipl) {
      fVtxFitMinStr.FirstWire[ipl] = -geom->WireCoordinate(0, 0, ipl, tpc, cstat);
      fVtxFitMinStr.OrthY[ipl] = geom->WireCoordinate(1, 0, ipl, tpc, cstat) + fVtxFitMinStr.FirstWire[ipl];
      fVtxFitMinStr.OrthZ[ipl] = geom->WireCoordinate(0, 1, ipl, tpc, cstat) + fVtxFitMinStr.FirstWire[ipl];
    }
    // and the vertex starting position
    fVtxFitMinStr.VtxPos = VtxPos;

    // and the track direction and hits
    fVtxFitMinStr.HitX = hitX;
    fVtxFitMinStr.HitXErr = hitXErr;
    fVtxFitMinStr.Plane.resize(ntrks);
    fVtxFitMinStr.Wire.resize(ntrks);
    for(itk = 0; itk < ntrks; ++itk) {
      fVtxFitMinStr.Plane[itk].resize(hitX[itk].size());
      fVtxFitMinStr.Wire[itk].resize(hitX[itk].size());
      for(iht = 0; iht < hitWID[itk].size(); ++iht) {
        fVtxFitMinStr.Plane[itk][iht] = hitWID[itk][iht].Plane;
        fVtxFitMinStr.Wire[itk][iht] = hitWID[itk][iht].Wire;
      }
    } // itk
    fVtxFitMinStr.Dir = TrkDir;
    
    fVtxFitMinStr.DoF = npts - npars;
      
    // define the starting parameters
    std::vector<double> par(npars);
    std::vector<double> stp(npars);
    std::vector<double> parerr(npars);
    
    TMinuit *gMin = new TMinuit(npars);
    gMin->SetFCN(VertexFitAlg::fcnVtxPos);
    int errFlag = 0;
    double arglist[10];
    
    // print level (-1 = none, 1 = yes)
    arglist[0] = -1;
    //
    // ***** remember to delete gMin before returning on an error
    //
    gMin->mnexcm("SET PRINT", arglist, 1, errFlag);
    
    // the vertex position
    unsigned short ipar;
    for(ipar = 0; ipar < 3; ++ipar) {
      par[ipar] = fVtxFitMinStr.VtxPos[ipar]; // in cm
      stp[ipar] = 0.1;  // 1 mm initial step
      gMin->mnparm(ipar,"", par[ipar], stp[ipar], -1E6, 1E6, errFlag);
    }
    // use Y, Z track directions. There is no constraint that the direction vector is unit-normalized
    // since we are only passing two of the components. Minuit could violate this requirement when
    // fitting. fcnVtxPos prevents non-physical values and Minuit may throw error messages as a result.
    for(itk = 0; itk < ntrks; ++itk) {
      ipar = 3 + 2 * itk;
      par[ipar]     = fVtxFitMinStr.Dir[itk](1);
      stp[ipar]     = 0.03;
      gMin->mnparm(ipar,"", par[ipar], stp[ipar], -1.05, 1.05, errFlag);
      ++ipar;
      par[ipar] = fVtxFitMinStr.Dir[itk](2);
      stp[ipar] = 0.03;
      gMin->mnparm(ipar,"", par[ipar], stp[ipar], -1.05, 1.05, errFlag);
    } // itk

/*
    // Single call to fcnVtxPos for debugging it
    std::cout<<"Starting: par  ";
    for(unsigned short ip = 0; ip < par.size(); ++ip) std::cout<<" "<<std::fixed<<std::setprecision(2)<<par[ip];
    std::cout<<"\n";
    // call fcn with starting parameters
    arglist[0] = 1;
    gMin->mnexcm("CALL", arglist, 1, errFlag);
*/
    // set strategy 0 for faster Minuit fitting
    arglist[0] = 0.;
    gMin->mnexcm("SET STRATEGY", arglist, 1, errFlag);
    
    // execute Minuit command: Migrad, max calls, tolerance
    arglist[0] = 500; // max calls
    arglist[1] = 1.; // tolerance on fval in fcn
    gMin->mnexcm("MIGRAD", arglist, 2, errFlag);
    
    // call fcn to get final fit values
    arglist[0] = 3;
    gMin->mnexcm("CALL", arglist, 1, errFlag);
    ChiDOF = fVtxFitMinStr.ChiDoF;
    
    // get the parameters
    for(unsigned short ip = 0; ip < par.size(); ++ip) {
      gMin->GetParameter(ip, par[ip], parerr[ip]);
    }

    // return the vertex position and errors
    for(ipar = 0; ipar < 3; ++ipar) {
      VtxPos[ipar] = par[ipar];
      VtxPosErr[ipar] = parerr[ipar];
    }
    // return the track directions and the direction errors if applicable
    bool returnTrkDirErrs = (TrkDirErr.size() == TrkDir.size());
    for(itk = 0; itk < ntrks; ++itk) {
      ipar = 3 + 2 * itk;
      double arg = 1 - par[ipar] * par[ipar] - par[ipar + 1] * par[ipar + 1];
      if(arg < 0) arg = 0;
      TrkDir[itk](0) = sqrt(arg);
      TrkDir[itk](1) = par[ipar];
      TrkDir[itk](2) = par[ipar + 1];
      if(returnTrkDirErrs) {
        // TODO This needs to be checked if there is a need for trajectory errors
        double errY = parerr[ipar] / par[ipar];
        double errZ = parerr[ipar + 1] / par[ipar + 1];
        TrkDirErr[itk](0) = sqrt(arg * (errY * errY + errZ * errZ));
        TrkDirErr[itk](1) = parerr[ipar];
        TrkDirErr[itk](2) = parerr[ipar + 1];
      }
    } // itk
    
    delete gMin;

  } // VertexFit()

} // namespace trkf

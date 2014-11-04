//////////////////////////////////////////////////////////////////////
///
/// TrackTrajectoryAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm fitting a 3D trajectory through a set of hit pairs
///
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>
#include <iostream>
#include <iomanip>

#include "RecoAlg/TrackTrajectoryAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

struct PlnLen{
  unsigned short index;
  unsigned short length;
};

bool greaterThan (PlnLen p1, PlnLen p2) { return (p1.length > p2.length);}

namespace trkf{

  TrackTrajectoryAlg::TrackTrajectoryAlg() { }

  TrackTrajectoryAlg::~TrackTrajectoryAlg() { }

//------------------------------------------------------------------------------
  void TrackTrajectoryAlg::TrackTrajectory(
      std::array<std::vector<std::pair<double, geo::WireID>>,3>& trkXW,
      std::vector<TVector3>& TrajPos, std::vector<TVector3>& TrajDir,
      std::array<std::vector<double>,3> trkChg)
  {
    // Make a track trajectory (position, direction) and put it in the TrajPos
    // and TrajDir vectors
    
    // Track hits must be ordered such that the first trkHit (X, WireID) in each
    // plane is at the start (or end) of the track. The trajectory is returned
    // in the same order. The number of trajectory points is set so that there
    // are 8 - 20 hits in each trajectory point in the plane with the 2nd-most hits.
    
    // If the trkChg array is NOT supplied, equidistant trajectory points are made.
    // If the trkChg array IS supplied, trajectory points are made so that the 
    // normalized charge is the same in all planes for each point.
    
    unsigned short ip, ipl, iht, ii;
    
    TrajPos.clear();
    TrajDir.clear();
  
    // sort the planes by decreasing number of hits
    std::vector<PlnLen> spl;
    PlnLen plnlen;
    // the number of planes that have hits
    unsigned short nplanes = 0;
    for(ipl = 0; ipl < 3; ++ipl) {
      plnlen.index = ipl;
      plnlen.length = trkXW[ipl].size();
      if(plnlen.length > 0) ++nplanes;
      spl.push_back(plnlen);
    }
    std::sort (spl.begin(),spl.end(), greaterThan);

    // spl[0] has the most hits and spl[2] has the least
    if(spl[1].length < 3) {
      mf::LogWarning("TrackTrajectoryAlg")<<"Not enough hits in two planes\n";
      return;
    } // spl[1].length < 3

    std::vector< std::pair<double, geo::WireID> > fitHits;
    
    // TrackLineFit results will be passed to these variables
    TVector3 xyz, dir;
    float ChiDOF;
    
    // determine how many hits should be included in the trajectory fit for each
    // range and in each plane. 
    ipl = spl[1].index;
    unsigned short nhtraj = 5;
    if(trkXW[ipl].size() > 50) {
      nhtraj = 20;
    } else if(trkXW[ipl].size() > 20) {
      nhtraj = 8;
    }
    // number of trajectory points in the plane with the 2nd-most hits
    unsigned short ntp = trkXW[ipl].size() / nhtraj;
    
    // Handle short tracks. Two trajectory points
    if(ntp < 3) {
      double xmin = 9999, xmax = 0;
      for(ip = 0; ip < nplanes; ++ip) {
        ipl = spl[ip].index;
        for(iht = 0; iht < trkXW[ipl].size(); ++iht) {
          if(trkXW[ipl][iht].first < xmin) xmin = trkXW[ipl][iht].first;
          if(trkXW[ipl][iht].first > xmax) xmax = trkXW[ipl][iht].first;
          fitHits.push_back(trkXW[ipl][iht]);
        }
      } // ip
      fTrackLineFitAlg.TrkLineFit(fitHits, xmin, xyz, dir, ChiDOF);
      // should get the order correct here....
      TrajPos.push_back(xyz);
      TrajDir.push_back(dir);
      for(ii = 0; ii < 3; ++ii) {
        xyz(ii) = TrajPos[0](ii) + (xmax - TrajPos[0](0)) * TrajDir[0](ii) / TrajDir[0](0);
      }
      TrajPos.push_back(xyz);
      TrajDir.push_back(TrajDir[0]);
      return;
    }

    unsigned short firsthit;
    double xOrigin;

    if(trkChg[0].size() == 0) {
/////////////////////// No trkChg ////////////////////////////////////////////////
      std::array<unsigned short, 3> nhtp;
      for(ii = 0; ii < 3; ++ii) nhtp[ii] = nhtraj * trkXW[ii].size() / trkXW[ipl].size();

      for(unsigned short itp = 0; itp < ntp; ++itp) {
        fitHits.clear();
        xOrigin = 0;
        for(ip = 0; ip < nplanes; ++ip) {
          ipl = spl[ip].index;
          firsthit = itp * nhtp[ipl];
          xOrigin += trkXW[ipl][firsthit].first;
          for(ii = 0; ii < nhtp[ipl]; ++ii) {
            iht = firsthit + ii;
            fitHits.push_back(trkXW[ipl][iht]);
          } // iht
        } // ipl
        xOrigin /= (double)nplanes;
        fTrackLineFitAlg.TrkLineFit(fitHits, xOrigin, xyz, dir, ChiDOF);
        TrajPos.push_back(xyz);
        TrajDir.push_back(dir);
      } // itp
      
      // make the last traj point
      fitHits.clear();
      xOrigin = 0;
      for(ip = 0; ip < nplanes; ++ip) {
        ipl = spl[ip].index;
        firsthit = trkXW[ipl].size() - 1;
        xOrigin += trkXW[ipl][firsthit].first;
        for(ii = 0; ii < nhtp[ipl]; ++ii) {
          iht = firsthit - ii;
          fitHits.push_back(trkXW[ipl][iht]);
        } //
      } // ip
      xOrigin /= (double)nplanes;
      fTrackLineFitAlg.TrkLineFit(fitHits, xOrigin, xyz, dir, ChiDOF);
      TrajPos.push_back(xyz);
      TrajDir.push_back(dir);
    } else {
/////////////////////// Use trkChg ////////////////////////////////////////////////
      // sum the charge in each plane
      std::array<double, 3> totChg;
      for(ip = 0; ip < nplanes; ++ip) {
        ipl = spl[ip].index;
        if(trkChg[ipl].size() != trkXW[ipl].size()) {
          mf::LogWarning("TrackTrajectoryAlg")<<"Incompatible trkChg length ";
          return;
        }
        totChg[ipl] = 0;
        for(iht = 0; iht < trkChg[ipl].size(); ++iht) totChg[ipl] += trkChg[ipl][iht];
        // divide by the number of trajectory points - 1
        totChg[ipl] /= (double)(ntp - 1);
      } // ii

      double chgSum;
      bool hitTheEnd = false;
      // array of hit start indices
      std::array<unsigned short, 3> hstart = {0,0,0};
      for(unsigned short itp = 0; itp < ntp - 1; ++itp) {
        fitHits.clear();
        xOrigin = 0;
        for(ip = 0; ip < nplanes; ++ip) {
          ipl = spl[ip].index;
          chgSum = 0;
          firsthit = hstart[ipl];
          xOrigin += trkXW[ipl][firsthit].first;
          for(iht = firsthit; iht < trkChg[ipl].size(); ++iht) {
            fitHits.push_back(trkXW[ipl][iht]);
            chgSum += trkChg[ipl][iht];
            if(chgSum > totChg[ipl]) {
              hstart[ipl] = iht + 1;
              break;
            } // chgSum > totChg[ipl]
          } // iht
          // hit the end of trkChg without reaching totChg. Indicates an
          // incomplete set of hits. Break out and fit the last trajectory point
          if(iht == trkChg[ipl].size()) hitTheEnd = true;
          if(hitTheEnd) break;
        } // ip
        if(hitTheEnd) break;
        xOrigin /= (double)nplanes;
        fTrackLineFitAlg.TrkLineFit(fitHits, xOrigin, xyz, dir, ChiDOF);
        if(ChiDOF < 0) continue;
        TrajPos.push_back(xyz);
        TrajDir.push_back(dir);
      } // itp
      
      // make the last traj point
      fitHits.clear();
      xOrigin = 0;
      for(ip = 0; ip < nplanes; ++ip) {
        ipl = spl[ip].index;
        chgSum = 0;
        firsthit = trkXW[ipl].size() - 1;
        xOrigin += trkXW[ipl][firsthit].first;
        for(ii = 0; ii < firsthit; ++ii) {
          iht = firsthit - ii;
          fitHits.push_back(trkXW[ipl][iht]);
          chgSum += trkChg[ipl][iht];
          if(chgSum > totChg[ipl]) break;
        } //
      } // ip
      xOrigin /= (double)nplanes;
      fTrackLineFitAlg.TrkLineFit(fitHits, xOrigin, xyz, dir, ChiDOF);
      if(ChiDOF > 0) {
        TrajPos.push_back(xyz);
        TrajDir.push_back(dir);
      }
    } // use trkChg
    
    
    if(TrajPos.size() > 4) {
      // smooth the YZ points since they are less well known than the XZ points
      std::vector<double> ypts(TrajPos.size());
      ypts[0] = (3 * TrajPos[0](1) +  TrajPos[1](1)) / 4;
      ii = TrajPos.size() - 1;
      ypts[ii] = (3 * TrajPos[ii](1) + TrajPos[ii-1](1)) / 4;
      for(ii = 1; ii < ypts.size() - 1; ++ii) 
        ypts[ii] = (TrajPos[ii-1](1) + 2 * TrajPos[ii](1) + TrajPos[ii+1](1)) / 4;
      for(ii = 1; ii < ypts.size(); ++ii)  TrajPos[ii](1) = ypts[ii];
    } // TrajPos.size() > 4

  } // TrackTrajectoryAlg

} // trkf

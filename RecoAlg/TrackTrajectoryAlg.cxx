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
  float length;
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
    
    unsigned int cstat = trkXW[ipl][0].second.Cryostat;
    unsigned int tpc = trkXW[ipl][0].second.TPC;
    
    // x,y,z of the endpoints
    TVector3 begpt, endpt;
    // Find the first and last trajectory points using hits in two planes that have the most
    // similar times. Check for the case where there are hits in only two planes
    unsigned int wire1, wire2, plane1, plane2, iht1, iht2;
    if(spl[2].length == 0) {
      // 2 planes - first end
      plane1 = spl[0].index; 
      plane2 = spl[1].index;
      wire1 = trkXW[plane1][0].second.Wire;
      wire2 = trkXW[plane2][0].second.Wire;
      begpt(0) = 0.5 * (trkXW[plane1][0].first + trkXW[plane2][0].first);
      geom->IntersectionPoint(wire1, wire2, plane1, plane2, cstat, tpc, 
        begpt(1), begpt(2));
      // other end
      iht1 = trkXW[plane1].size() - 1;
      iht2 = trkXW[plane2].size() - 1;
      wire1 = trkXW[plane1][iht1].second.Wire;
      wire2 = trkXW[plane2][iht2].second.Wire;
      endpt(0) = 0.5 * (trkXW[plane1][iht1].first + trkXW[plane2][iht2].first);
      geom->IntersectionPoint(wire1, wire2, plane1, plane2, cstat, tpc, 
        endpt(1), endpt(2));
    } else {
      // sort the X positions of the first hit in each plane
      std::vector<PlnLen> xpl;
      for(ipl = 0; ipl < 3; ++ipl) {
        plnlen.index = ipl;
        plnlen.length = trkXW[ipl][0].first;
        xpl.push_back(plnlen);
      } // ipl
      std::sort (xpl.begin(),xpl.end(), greaterThan);
      // find which two hits have the most similar times
      if(fabs(xpl[xpl[0].index].length - xpl[xpl[1].index].length) < 
         fabs(xpl[xpl[1].index].length - xpl[xpl[2].index].length)) {
        plane1 = xpl[0].index;
        plane2 = xpl[1].index;
      } else {
        plane1 = xpl[1].index;
        plane2 = xpl[2].index;
      }
      wire1 = trkXW[plane1][0].second.Wire;
      wire2 = trkXW[plane2][0].second.Wire;
      begpt(0) = 0.5 * (trkXW[plane1][0].first + trkXW[plane2][0].first);
      geom->IntersectionPoint(wire1, wire2, plane1, plane2, cstat, tpc,
        begpt(1), begpt(2));
/*
  std::cout<<"begpt using P:W "<<plane1<<":"<<wire1<<" "<<plane2<<":"<<wire2
    <<" X1 "<<trkXW[plane1][0].first<<" X2 "<<trkXW[plane2][0].first
    <<" Y "<<begpt(1)<<" Z "<<begpt(2)<<"\n";
*/
      // now for the other end using the last hit
      xpl.clear();
      for(ipl = 0; ipl < 3; ++ipl) {
        plnlen.index = ipl;
        plnlen.length = trkXW[ipl][trkXW[ipl].size()-1].first;
        xpl.push_back(plnlen);
      } // ipl
      std::sort (xpl.begin(),xpl.end(), greaterThan);
      // find which two hits have the most similar times
      if(fabs(xpl[xpl[0].index].length - xpl[xpl[1].index].length) < 
         fabs(xpl[xpl[1].index].length - xpl[xpl[2].index].length)) {
        plane1 = xpl[0].index;
        plane2 = xpl[1].index;
      } else {
        plane1 = xpl[1].index;
        plane2 = xpl[2].index;
      }
      iht1 = trkXW[plane1].size() - 1;
      iht2 = trkXW[plane2].size() - 1;
      wire1 = trkXW[plane1][iht1].second.Wire;
      wire2 = trkXW[plane2][iht2].second.Wire;
      endpt(0) = 0.5 * (trkXW[plane1][iht1].first + trkXW[plane2][iht2].first);
      geom->IntersectionPoint(wire1, wire2, plane1, plane2, cstat, tpc,
        endpt(1), endpt(2));
/*
  std::cout<<"endpt using P:W "<<plane1<<":"<<wire1<<" "<<plane2<<":"<<wire2
    <<" X1 "<<trkXW[plane1][iht1].first<<" X2 "<<trkXW[plane2][iht2].first
    <<" Y "<<endpt(1)<<" Z "<<endpt(2)<<"\n";
*/
    } // spl[2].length

    // start the trajectory with the first end point
    TrajPos.push_back(begpt);

    // Handle short tracks. Two trajectory points
    if(ntp < 3) {
      // the other end point
      TrajPos.push_back(endpt);
      dir = endpt - begpt;
      dir = dir.Unit();
      TrajDir.push_back(dir);
      TrajDir.push_back(dir);
      return;
    }

    unsigned short firsthit;
    double xOrigin;

    if(trkChg[0].size() == 0) {
/////////////////////// No trkChg ////////////////////////////////////////////////
      ipl = spl[1].index;
      std::array<unsigned short, 3> nhtp;
      for(ii = 0; ii < 3; ++ii) nhtp[ii] = nhtraj * trkXW[ii].size() / trkXW[ipl].size();

      bool first = true;
      for(unsigned short itp = 1; itp < ntp; ++itp) {
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
        if(ChiDOF < 0) continue;
        if(first) {
          // put in the first traj point direction
          TVector3 tmp;
          tmp = xyz - begpt;
          tmp = tmp.Unit();
          TrajDir.push_back(tmp);
          first = false;
        }
        TrajPos.push_back(xyz);
        TrajDir.push_back(dir);
      } // itp
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
      bool first = true;
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
        if(first) {
          // put in the first traj point direction
          TVector3 tmp;
          tmp = TrajPos[TrajPos.size()-1] - begpt;
          tmp = tmp.Unit();
          TrajDir.push_back(tmp);
          first = false;
        }
        TrajPos.push_back(xyz);
        TrajDir.push_back(dir);
      } // itp
    } // use trkChg

    // make the last traj point
    TrajPos.push_back(endpt);
    dir = endpt - TrajPos[TrajPos.size()-1];
    dir = dir.Unit();
    TrajDir.push_back(dir);
    
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

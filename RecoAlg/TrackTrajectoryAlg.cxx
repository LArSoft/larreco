//////////////////////////////////////////////////////////////////////
///
/// TrackTrajectoryAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm fitting a 3D trajectory through a set of hits
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

namespace trkf{

  TrackTrajectoryAlg::TrackTrajectoryAlg() { }

  TrackTrajectoryAlg::~TrackTrajectoryAlg() { }

//------------------------------------------------------------------------------
  void TrackTrajectoryAlg::TrackTrajectory(std::array<std::vector<geo::WireID>,3> trkWID,
                                           std::array<std::vector<double>,3> trkX,
                                           std::array<std::vector<double>,3> trkXErr,
                                           std::vector<TVector3>& TrajPos, std::vector<TVector3>& TrajDir)
  {

    // Make a track trajectory (position, direction) and return it in the TrajPos
    // and TrajDir vectors. The track hits are received as 3 vectors (one vector per wire plane)
    // of Wire ID's, hit X and X position errors. The X position errors are used to 1) determine
    // the significance of the X difference between the beginning and end of the track trajectory
    // 2) determine the number of trajectory points and 3) weight the track line fit in TrackLineFitAlg.
    // This code assumes that hits at each end (e.g. trkXW[Plane][0]) of the vectors define the end
    // points of the trajectory. The ordering of planes in the array is irrelevant since the
    // plane number is extracted from the WireID. This algorithm will return with a failed condition
    // (TrajPos, TrajDir size = 0) if there are fewer than 2 planes with hits at each end that are less
    // than 5 * trkXErr apart. Valid and invalid conditions are shown schematically below where a . represents
    // hits that are separated by X values > 5 * trkXErr
    //
    //      minX            maxX         maxX            minX          minX            maxX
    // Pln0 ....................    Pln0 ....................     Pln0 ...................
    // Pln1 ...............         Pln1 ....................     Pln1
    // Pln2 ....................    Pln2 ................         Pln2 ...................
    // VALID                        VALID                         VALID - no hits in one plane is OK
    //
    //      minX            maxX
    // Pln0 .................
    // Pln1 ...............
    // Pln2 ....................
    // NOT VALID - Only one plane has a hit at MaxX
    

    TrajPos.clear();
    TrajDir.clear();
    
    bool prt = false;
    
    unsigned short minLen = 9999;
    unsigned short maxLen = 0;
    unsigned short nPlnsWithHits = 0;
    unsigned short ipl;
    for(ipl = 0; ipl < 3; ++ipl) {
      if(trkX[ipl].size() == 0) continue;
      ++nPlnsWithHits;
      if(trkX[ipl].size() < minLen) {
        minLen = trkX[ipl].size();
      }
      if(trkX[ipl].size() > maxLen) {
        maxLen = trkX[ipl].size();
      }
    } // ipl
    if(prt) std::cout<<"nPlnsWithHits "<<nPlnsWithHits<<"\n";
    if(nPlnsWithHits < 2) return;

    unsigned short iht;
    
    // reverse the order of the hit arrays if necessary so that the minimum X end is at trk[][0] end.
    // We will use posX to reverse the trajectory later if necessary
    bool posX = true;
    iht = trkX[0].size() - 1;
    if(trkX[0][0] > trkX[0][iht]) {
      posX = false;
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        std::reverse(trkWID[ipl].begin(), trkWID[ipl].end());
        std::reverse(trkX[ipl].begin(), trkX[ipl].end());
        std::reverse(trkXErr[ipl].begin(), trkXErr[ipl].end());
      } // ipl
      if(prt) std::cout<<"Swapped order\n";
    } // posX check
    
    if(prt) {
      std::cout<<"TrkXW end0";
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        std::cout<<" "<<trkX[ipl][0];
      }
      std::cout<<"\n";
      std::cout<<"TrkXW end1";
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        iht = trkX[ipl].size() - 1;
        std::cout<<" "<<trkX[ipl][iht];
      }
      std::cout<<"\n";
    }
    
    // find the min/max X
    minX = 1E6;
    minXPln = 0;
    maxX = -1E6;
    maxXPln = 0;
    
    for(unsigned short ipl = 0; ipl < 3; ++ipl) {
      if(trkX[ipl].size() == 0) continue;
      if(trkX[ipl][0] < minX) {
        minX = trkX[ipl][0];
        minXPln = ipl;
      }
      iht = trkX[ipl].size() - 1;
      if(trkX[ipl][iht] > maxX) {
        maxX = trkX[ipl][iht];
        maxXPln = ipl;
      }
    } // ipl
    
    if(prt) std::cout<<"minX "<<minX<<" in plane "<<minXPln<<" maxX "<<maxX<<" in plane "<<maxXPln<<"\n";
    
    // estimate the number of trajectory points we will want based on the delta T and the
    // hit errors
    double aveHitXErr = 0;
    unsigned short nHit = 0;
    for(ipl = 0; ipl < 3; ++ipl) {
      for(iht = 0; iht < trkXErr[ipl].size(); ++iht) {
        aveHitXErr += trkXErr[ipl][iht];
        ++nHit;
      }
    } // ipl
    aveHitXErr /= (double)nHit;
    
    unsigned short npt = (maxX - minX) / (1 * aveHitXErr);
    if(npt < 2) npt = 2;
    if(npt > maxLen) npt = maxLen;
    if(npt > 100) npt = 100;
    if(prt) std::cout<<" aveHitXErr "<<aveHitXErr<<" number of traj points "<<npt<<"\n";

    double maxBinX = (maxX - minX) / (double)(npt - 1);
    double binX = maxBinX;
    double xOrigin = minX;
    
    TVector3 xyz, dir;
    // working vectors passed to TrackLineFitAlg
    std::vector<geo::WireID> hitWID;
    std::vector<double> hitX;
    std::vector<double> hitXErr;
    float ChiDOF;

    if(maxLen < 4 || npt < 2) {
      ShortTrackTrajectory(trkWID, trkX, trkXErr, TrajPos, TrajDir);
      if(!posX) {
        // reverse everything
        std::reverse(TrajPos.begin(), TrajPos.end());
        std::reverse(TrajDir.begin(), TrajDir.end());
      }
      return;
    } // maxLen < 4
    
    // Start index of hits to include in the next fit
    std::array<unsigned short, 3> hStart;
    for(ipl = 0; ipl < 3; ++ipl) hStart[ipl] = 0;
    
    for(unsigned short ipt = 0; ipt < npt + 1; ++ipt) {
      hitWID.clear();
      hitX.clear();
      hitXErr.clear();
      for(ipl = 0; ipl < 3; ++ipl) {
        for(iht = hStart[ipl]; iht < trkX[ipl].size(); ++iht) {
          if(trkX[ipl][iht] < xOrigin - binX) continue;
          if(trkX[ipl][iht] > xOrigin + binX) break;
          hitWID.push_back(trkWID[ipl][iht]);
          hitX.push_back(trkX[ipl][iht]);
          hitXErr.push_back(trkXErr[ipl][iht]);
          hStart[ipl] =iht;
        } // iht
      } // ipl
      if(hitX.size() < 4) continue;
      fTrackLineFitAlg.TrkLineFit(hitWID, hitX, hitXErr, xOrigin, xyz, dir, ChiDOF);
      if(prt) std::cout<<" xyz "<<xyz(0)<<" "<<xyz(1)<<" "<<xyz(2)<<" dir "<<dir(0)<<" "<<dir(1)<<" "<<dir(2)<<" ChiDOF "<<ChiDOF<<" hitX size "<<hitX.size()<<"\n";
      // tweak xOrigin if we are close to the end
      if(xOrigin >= maxX) break;
      if(maxX - xOrigin < binX) {
        xOrigin = maxX;
      } else {
        xOrigin += binX;
      }
      if(ChiDOF < 0 || ChiDOF > 100) continue;
      TrajPos.push_back(xyz);
      TrajDir.push_back(dir);
    } // ipt
    if(TrajPos.size() < 2) {
      TrajPos.clear();
      TrajDir.clear();
      ShortTrackTrajectory(trkWID, trkX, trkXErr, TrajPos, TrajDir);
    }
    if(!posX) {
      // reverse everything
      std::reverse(TrajPos.begin(), TrajPos.end());
      std::reverse(TrajDir.begin(), TrajDir.end());
    } // !posX
/*
    // smooth trajectory Y points
    if(TrajPos.size() > 4) {
      unsigned short ii;
      std::vector<double> ypts(TrajPos.size());
      ypts[0] = (3 * TrajPos[0](1) +  TrajPos[1](1)) / 4;
      ii = TrajPos.size() - 1;
      ypts[ii] = (3 * TrajPos[ii](1) + TrajPos[ii-1](1)) / 4;
      for(ii = 1; ii < ypts.size() - 1; ++ii) ypts[ii] = (TrajPos[ii-1](1) + 2 * TrajPos[ii](1) + TrajPos[ii+1](1)) / 4;
      for(ii = 1; ii < ypts.size(); ++ii)  TrajPos[ii](1) = ypts[ii];
    }
*/
    
  } // TrackTrajectoryAlg

////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void TrackTrajectoryAlg::ShortTrackTrajectory(std::array<std::vector<geo::WireID>,3> trkWID,
                                                std::array<std::vector<double>,3> trkX,
                                                std::array<std::vector<double>,3> trkXErr,
                                                std::vector<TVector3>& TrajPos, std::vector<TVector3>& TrajDir)
  {
    // Make a short track trajectory composed of a space point at each end.
    
    std::vector<unsigned short> usePlns;
    double y, z, endX, endY, endZ;
    unsigned short ipl, iht, nend, jpl, iPlane, jPlane;
    unsigned int iWire, jWire;
    
    // do the min X end first
    float xCut = minX + 5 * trkXErr[minXPln][0];
    for(ipl = 0; ipl < 3; ++ipl) {
      if(trkX[ipl].size() == 0) continue;
      if(trkX[ipl][0] < xCut) usePlns.push_back(ipl);
    } // ipl
    // Not enough information to find a space point
//    std::cout<<"ShortTrack minX end "<<xCut<<" usePlns size "<<usePlns.size()<<"\n";
    if(usePlns.size() < 2) return;
    endY = 0; endZ = 0; nend = 0;
    iht = 0;
    ipl = usePlns[0];
    endX   = trkX[ipl][iht];
    iPlane = trkWID[ipl][iht].Plane;
    iWire  = trkWID[ipl][iht].Wire;
    unsigned int cstat = trkWID[ipl][iht].Cryostat;
    unsigned int tpc = trkWID[ipl][iht].TPC;
    for(unsigned short jj = 1; jj < usePlns.size(); ++jj) {
      jpl = usePlns[jj];
      endX  += trkX[jpl][iht];
      jPlane = trkWID[jpl][iht].Plane;
      jWire  = trkWID[jpl][iht].Wire;
      geom->IntersectionPoint(iWire, jWire, iPlane, jPlane, cstat, tpc, y, z);
      endY += y;
      endZ += z;
      ++nend;
    } // ii
    // nend is the number of y,z values. There are (nend + 1) X values
    TVector3 xyz;
    xyz(0) = endX / (float)(nend + 1);
    xyz(1) = endY / (float)nend;
    xyz(2) = endZ / (float)nend;
//    std::cout<<" xyz "<<xyz(0)<<" "<<xyz(1)<<" "<<xyz(2)<<"\n";
    TrajPos.push_back(xyz);
    
    // do the same for the other end
    xCut = maxX - 5 * trkXErr[maxXPln][0];
    usePlns.clear();
    for(ipl = 0; ipl < 3; ++ipl) {
      if(trkX[ipl].size() == 0) continue;
      if(trkX[ipl][0] > xCut) usePlns.push_back(ipl);
    } // ipl
    // Not enough information to find a space point
//    std::cout<<"ShortTrack maxX end "<<xCut<<" usePlns size "<<usePlns.size()<<"\n";
    if(usePlns.size() < 2) {
      TrajPos.clear();
      TrajDir.clear();
      return;
    }
    endY = 0; endZ = 0; nend = 0;
    ipl = usePlns[0];
    iht = trkX[ipl].size() - 1;
    endX = trkX[ipl][iht];
    iPlane = trkWID[ipl][iht].Plane;
    iWire  = trkWID[ipl][iht].Wire;
    for(unsigned short jj = 1; jj < usePlns.size(); ++jj) {
      jpl = usePlns[jj];
      iht = trkX[jpl].size() - 1;
      endX += trkX[jpl][iht];
      jPlane = trkWID[jpl][iht].Plane;
      jWire  = trkWID[jpl][iht].Wire;
      geom->IntersectionPoint(iWire, jWire, iPlane, jPlane, cstat, tpc, y, z);
      endY += y;
      endZ += z;
      ++nend;
    } // ii
    // nend is the number of y,z values. There are nend + 1 X values
    xyz(0) = endX / (float)(nend + 1);
    xyz(1) = endY / (float)nend;
    xyz(2) = endZ / (float)nend;
//    std::cout<<" xyz "<<xyz(0)<<" "<<xyz(1)<<" "<<xyz(2)<<"\n";
    TrajPos.push_back(xyz);
    TVector3 dir = TrajPos[1] - TrajPos[0];
    dir = dir.Unit();
    TrajDir.push_back(dir);
    TrajDir.push_back(dir);
  
//    std::cout<<">>>> Short track ("<<TrajPos[0](0)<<", "<<TrajPos[0](1)<<", "<<TrajPos[0](2)<<") to ("<<TrajPos[1](0)<<", "<<TrajPos[1](1)<<", "<<TrajPos[1](2)<<")\n";
  }

  
} // trkf

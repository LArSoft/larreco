//////////////////////////////////////////////////////////////////////
///
/// TrackTrajectoryAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm fitting a 3D trajectory through a set of hits
///
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iomanip>

#include "TVector3.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackTrajectoryAlg.h"

namespace trkf{

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
    // than fHitWidthFactor * trkXErr apart. Valid and invalid conditions are shown schematically below
    // where a . represents hits that are separated by X values > fHitWidthFactor * trkXErr
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

    prt = false;

    unsigned short minLen = 9999;
    unsigned short maxLen = 0;
    unsigned short nPlnsWithHits = 0;
    unsigned short ipl, aPlane = 3;
    for(ipl = 0; ipl < 3; ++ipl) {
      if(trkX[ipl].size() == 0) continue;
      ++nPlnsWithHits;
      if(trkX[ipl].size() < minLen) {
        minLen = trkX[ipl].size();
      }
      if(trkX[ipl].size() > maxLen) {
        maxLen = trkX[ipl].size();
        aPlane = ipl;
      }
    } // ipl
    if(prt) mf::LogVerbatim("TTA")<<"trkX sizes "<<trkX[0].size()<<" "<<trkX[1].size()<<" "<<trkX[2].size()<<" "<<" nPlnsWithHits "<<nPlnsWithHits;
    if(nPlnsWithHits < 2) return;
    if(aPlane > 2) return;

    fMaxTrajPoints = 100;
    fHitWidthFactor = 5.;

    unsigned short iht;

    // reverse the order of the hit arrays if necessary so that the minimum X end is at trk[][0] end.
    // We will use posX to reverse the trajectory later if necessary
    bool posX = true;
    iht = trkX[aPlane].size() - 1;
    if(trkX[aPlane][0] > trkX[aPlane][iht]) {
      posX = false;
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        std::reverse(trkWID[ipl].begin(), trkWID[ipl].end());
        std::reverse(trkX[ipl].begin(), trkX[ipl].end());
        std::reverse(trkXErr[ipl].begin(), trkXErr[ipl].end());
      } // ipl
      if(prt) mf::LogVerbatim("TTA")<<"Swapped order";
    } // posX check

    if(prt) {
      mf::LogVerbatim myprt("TTA");
      myprt<<"TrkXW end0";
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        myprt<<" "<<trkX[ipl][0];
      }
      myprt<<"\n";
      myprt<<"TrkXW end1";
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        iht = trkX[ipl].size() - 1;
        myprt<<" "<<trkX[ipl][iht];
      }
      myprt<<"\n";
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

    if(prt) mf::LogVerbatim("TTA")<<"minX "<<minX<<" in plane "<<minXPln<<" maxX "<<maxX<<" in plane "<<maxXPln;

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
    if(npt > fMaxTrajPoints) npt = fMaxTrajPoints;
    if(prt) mf::LogVerbatim("TTA")<<" aveHitXErr "<<aveHitXErr<<" number of traj points ";

    double maxBinX = (maxX - minX) / (double)(npt - 1);
    double binX = maxBinX;
    double xOrigin = minX;

    TVector3 xyz, dir;
    // working vectors passed to TrackLineFitAlg
    std::vector<geo::WireID> hitWID;
    std::vector<double> hitX;
    std::vector<double> hitXErr;
    float ChiDOF;

    // make a short track trajectory (end points only) to use in case an error occurs
    std::vector<TVector3> STPos;
    std::vector<TVector3> STDir;
    ShortTrackTrajectory(trkWID, trkX, trkXErr, STPos, STDir);

    if(STPos.size() != 2 || STDir.size() != STPos.size()) {
      TrajPos.clear();
      TrajDir.clear();
      if(posX) {
        for(ipl = 0; ipl < 3; ++ipl) {
          if(trkX[ipl].size() == 0) continue;
          std::reverse(trkWID[ipl].begin(), trkWID[ipl].end());
          std::reverse(trkX[ipl].begin(), trkX[ipl].end());
          std::reverse(trkXErr[ipl].begin(), trkXErr[ipl].end());
        } // ipl
      } // posX
      return;
    } // bad STPos, STDir

    if(prt) {
      mf::LogVerbatim("TTA")<<"STPos";
      for(unsigned short ii = 0; ii < STPos.size(); ++ii) mf::LogVerbatim("TTA")<<ii<<" "<<std::fixed<<std::setprecision(1)<<STPos[ii](0)<<" "<<STPos[ii](1)<<" "<<STPos[ii](2);
    }

    if(maxLen < 4 || npt < 2) {
      TrajPos = STPos;
      TrajDir = STDir;
      if(!posX) {
        // reverse everything
        std::reverse(TrajPos.begin(), TrajPos.end());
        std::reverse(TrajDir.begin(), TrajDir.end());
        for(ipl = 0; ipl < 3; ++ipl) {
          if(trkX[ipl].size() == 0) continue;
          std::reverse(trkWID[ipl].begin(), trkWID[ipl].end());
          std::reverse(trkX[ipl].begin(), trkX[ipl].end());
          std::reverse(trkXErr[ipl].begin(), trkXErr[ipl].end());
        } // ipl
      }
      return;
    } // maxLen < 4

    // Start index of hits to include in the next fit
    std::array<unsigned short, 3> hStart;
    for(ipl = 0; ipl < 3; ++ipl) hStart[ipl] = 0;

    bool gotLastPoint = false;
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
      if(prt) mf::LogVerbatim("TTA")<<"ipt "<<ipt<<" xOrigin "<<xOrigin<<" binX "<<binX<<" hitX size "<<hitX.size();
      if(hitX.size() > 3) {
        fTrackLineFitAlg.TrkLineFit(hitWID, hitX, hitXErr, xOrigin, xyz, dir, ChiDOF);
        if(prt) mf::LogVerbatim("TTA")<<" xyz "<<xyz(0)<<" "<<xyz(1)<<" "<<xyz(2)<<" dir "<<dir(0)<<" "<<dir(1)<<" "<<dir(2)<<" ChiDOF "<<ChiDOF<<" hitX size "<<hitX.size();
      } else if(ipt == 0 && STPos.size() == 2) {
        // failure on the first traj point. Use STPos
        xyz = STPos[0];
        dir = STDir[0];
      }
      if(prt && hitX.size() < 4) mf::LogVerbatim("TTA")<<"\n";
      if(xOrigin >= maxX) break;
      // tweak xOrigin if we are close to the end
      if(maxX - xOrigin < binX) {
        xOrigin = maxX;
      } else {
        xOrigin += binX;
      }
      if(ChiDOF < 0 || ChiDOF > 100) continue;
      TrajPos.push_back(xyz);
      TrajDir.push_back(dir);
      if(ipt == npt) gotLastPoint = true;
    } // ipt
    if(prt) {
      mf::LogVerbatim("TTA")<<"gotLastPoint "<<gotLastPoint<<" TTA Traj \n";
      for(unsigned short ii = 0; ii < TrajPos.size(); ++ii) mf::LogVerbatim("TTA")<<ii<<" "<<std::fixed<<std::setprecision(1)<<TrajPos[ii](0)<<" "<<TrajPos[ii](1)<<" "<<TrajPos[ii](2);
    }

    if(!gotLastPoint && STPos.size() == 2) {
      // failure on the last point. Use STPos, etc
      TrajPos.push_back(STPos[1]);
      TrajDir.push_back(STDir[1]);
    }

    if(TrajPos.size() < 2) {
      TrajPos = STPos;
      TrajDir = STDir;
    }

    if(!posX) {
      // reverse everything
      std::reverse(TrajPos.begin(), TrajPos.end());
      std::reverse(TrajDir.begin(), TrajDir.end());
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        std::reverse(trkWID[ipl].begin(), trkWID[ipl].end());
        std::reverse(trkX[ipl].begin(), trkX[ipl].end());
        std::reverse(trkXErr[ipl].begin(), trkXErr[ipl].end());
      } // ipl

    } // !posX
    if(prt) {
      mf::LogVerbatim("TTA")<<"TTA Traj2\n";
      for(unsigned short ii = 0; ii < TrajPos.size(); ++ii) mf::LogVerbatim("TTA")<<ii<<" "<<std::fixed<<std::setprecision(1)<<TrajPos[ii](0)<<" "<<TrajPos[ii](1)<<" "<<TrajPos[ii](2);
    }

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
    float xCut;
    for(unsigned short nit = 0; nit < 5; ++nit) {
      xCut = minX + fHitWidthFactor * trkXErr[minXPln][0];
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        if(trkX[ipl][0] < xCut) usePlns.push_back(ipl);
      } // ipl
      if(usePlns.size() >= 2) break;
      fHitWidthFactor += 5;
    }
    // Not enough information to find a space point
    if(prt) mf::LogVerbatim("TTA")<<"ShortTrack minX end "<<xCut<<" usePlns size "<<usePlns.size();
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
    if(prt) mf::LogVerbatim("TTA")<<" xyz "<<xyz(0)<<" "<<xyz(1)<<" "<<xyz(2);
    TrajPos.push_back(xyz);

    // do the same for the other end
    fHitWidthFactor = 5;
//    xCut = maxX - fHitWidthFactor * trkXErr[maxXPln][0];
    usePlns.clear();
    for(unsigned short nit = 0; nit < 5; ++nit) {
      xCut = maxX - fHitWidthFactor * trkXErr[maxXPln][0];
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkX[ipl].size() == 0) continue;
        iht = trkX[ipl].size() - 1;
        if(trkX[ipl][iht] > xCut) usePlns.push_back(ipl);
      } // ipl
      if(usePlns.size() >= 2) break;
      fHitWidthFactor += 5;
    }
    // Not enough information to find a space point
    if(prt) mf::LogVerbatim("TTA")<<"ShortTrack maxX end "<<xCut<<" usePlns size "<<usePlns.size();
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
    if(prt) mf::LogVerbatim("TTA")<<" xyz "<<xyz(0)<<" "<<xyz(1)<<" "<<xyz(2);
    TrajPos.push_back(xyz);
    TVector3 dir = TrajPos[1] - TrajPos[0];
    dir = dir.Unit();
    TrajDir.push_back(dir);
    TrajDir.push_back(dir);

    if(prt) mf::LogVerbatim("TTA")<<">>>> Short track ("<<TrajPos[0](0)<<", "<<TrajPos[0](1)<<", "<<TrajPos[0](2)<<") to ("<<TrajPos[1](0)<<", "<<TrajPos[1](1)<<", "<<TrajPos[1](2)<<")";
  }


} // trkf

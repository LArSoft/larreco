#include "larreco/RecoAlg/TCAlg/TCCR.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <float.h>
#include <map>
#include <math.h>
#include <utility>
#include <vector>

namespace tca {

  ////////////////////////////////////////////////
  void SaveCRInfo(detinfo::DetectorClocksData const& clockData,
                  TCSlice& slc,
                  PFPStruct& pfp,
                  bool prt,
                  bool fIsRealData)
  {

    //Check the origin of pfp
    if (tcc.modes[kSaveCRTree]) {
      if (fIsRealData) { slc.crt.cr_origin.push_back(-1); }
      else {
        slc.crt.cr_origin.push_back(GetOrigin(clockData, slc, pfp));
      }
    }

    // save the xmin and xmax of each pfp
    auto& startPos = pfp.TP3Ds[0].Pos;
    auto& endPos = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    slc.crt.cr_pfpxmin.push_back(std::min(startPos[0], endPos[0]));
    slc.crt.cr_pfpxmax.push_back(std::max(startPos[0], endPos[0]));

    //find max
    const geo::TPCGeo& tpc = tcc.geom->TPC(0);
    float mindis0 = FLT_MAX;
    float mindis1 = FLT_MAX;
    if (std::abs(startPos[1] - tpc.MinY()) < mindis0) mindis0 = std::abs(startPos[1] - tpc.MinY());
    if (std::abs(startPos[1] - tpc.MaxY()) < mindis0) mindis0 = std::abs(startPos[1] - tpc.MaxY());
    if (std::abs(startPos[2] - tpc.MinZ()) < mindis0) mindis0 = std::abs(startPos[2] - tpc.MinZ());
    if (std::abs(startPos[2] - tpc.MaxZ()) < mindis0) mindis0 = std::abs(startPos[2] - tpc.MaxZ());
    if (std::abs(endPos[1] - tpc.MinY()) < mindis1) mindis1 = std::abs(endPos[1] - tpc.MinY());
    if (std::abs(endPos[1] - tpc.MaxY()) < mindis1) mindis1 = std::abs(endPos[1] - tpc.MaxY());
    if (std::abs(endPos[2] - tpc.MinZ()) < mindis1) mindis1 = std::abs(endPos[2] - tpc.MinZ());
    if (std::abs(endPos[2] - tpc.MaxZ()) < mindis1) mindis1 = std::abs(endPos[2] - tpc.MaxZ());
    //std::cout<<startPos[1]<<" "<<startPos[2]<<" "<<endPos[1]<<" "<<endPos[2]<<" "<<tpc.MinY()<<" "<<tpc.MaxY()<<" "<<tpc.MinZ()<<" "<<tpc.MaxZ()<<" "<<mindis0<<" "<<mindis1<<" "<<mindis0+mindis1<<std::endl;
    slc.crt.cr_pfpyzmindis.push_back(mindis0 + mindis1);

    if (slc.crt.cr_pfpxmin.back() < -2 || slc.crt.cr_pfpxmax.back() > 260 ||
        slc.crt.cr_pfpyzmindis.back() < 30) {
      pfp.CosmicScore = 1.;
    }
    else
      pfp.CosmicScore = 0;
  }

  ////////////////////////////////////////////////
  int GetOrigin(detinfo::DetectorClocksData const& clockData, TCSlice& slc, PFPStruct& pfp)
  {

    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;

    std::map<int, float> omap; //<origin, energy>

    for (auto& tjID : pfp.TjIDs) {

      Trajectory& tj = slc.tjs[tjID - 1];
      for (auto& tp : tj.Pts) {
        for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if (!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          TCHit& slhit = slc.slHits[iht];
          auto& hit = (*evt.allHits)[slhit.allHitsIndex];
          raw::ChannelID_t channel = tcc.geom->PlaneWireToChannel((int)hit.WireID().Plane,
                                                                  (int)hit.WireID().Wire,
                                                                  (int)hit.WireID().TPC,
                                                                  (int)hit.WireID().Cryostat);
          double startTick = hit.PeakTime() - hit.RMS();
          double endTick = hit.PeakTime() + hit.RMS();
          // get a list of track IDEs that are close to this hit
          std::vector<sim::TrackIDE> tides;
          tides = bt_serv->ChannelToTrackIDEs(clockData, channel, startTick, endTick);
          for (auto itide = tides.begin(); itide != tides.end(); ++itide) {
            omap[pi_serv->TrackIdToMCTruth_P(itide->trackID)->Origin()] += itide->energy;
          }
        }
      }
    }

    float maxe = -1;
    int origin = 0;
    for (auto& i : omap) {
      if (i.second > maxe) {
        maxe = i.second;
        origin = i.first;
      }
    }
    return origin;
  }

  ////////////////////////////////////////////////
  void ClearCRInfo(TCSlice& slc)
  {
    slc.crt.cr_origin.clear();
    slc.crt.cr_pfpxmin.clear();
    slc.crt.cr_pfpxmax.clear();
    slc.crt.cr_pfpyzmindis.clear();
  }
}

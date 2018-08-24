#include "larreco/RecoAlg/TCAlg/TCCR.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace tca {

  ////////////////////////////////////////////////
  void SaveCRInfo(TCSlice& slc, PFPStruct& pfp, bool prt, bool fIsRealData){

    //Check the origin of pfp
    if (tcc.modes[kSaveCRTree]){
      if (fIsRealData){
        slc.crt.cr_origin.push_back(-1);
      }
      else{
        slc.crt.cr_origin.push_back(GetOrigin(slc, pfp));
      }
    }

    // save the xmin and xmax of each pfp
    slc.crt.cr_pfpxmin.push_back(std::min(pfp.XYZ[0][0], pfp.XYZ[1][0]));
    slc.crt.cr_pfpxmax.push_back(std::max(pfp.XYZ[0][0], pfp.XYZ[1][0]));

    //find max 
    const geo::TPCGeo &tpc = tcc.geom->TPC(0);
    float mindis0 = FLT_MAX;
    float mindis1 = FLT_MAX;
    if (std::abs(pfp.XYZ[0][1] - tpc.MinY())<mindis0) mindis0 = std::abs(pfp.XYZ[0][1] - tpc.MinY());
    if (std::abs(pfp.XYZ[0][1] - tpc.MaxY())<mindis0) mindis0 = std::abs(pfp.XYZ[0][1] - tpc.MaxY());
    if (std::abs(pfp.XYZ[0][2] - tpc.MinZ())<mindis0) mindis0 = std::abs(pfp.XYZ[0][2] - tpc.MinZ());
    if (std::abs(pfp.XYZ[0][2] - tpc.MaxZ())<mindis0) mindis0 = std::abs(pfp.XYZ[0][2] - tpc.MaxZ());
    if (std::abs(pfp.XYZ[1][1] - tpc.MinY())<mindis1) mindis1 = std::abs(pfp.XYZ[1][1] - tpc.MinY());
    if (std::abs(pfp.XYZ[1][1] - tpc.MaxY())<mindis1) mindis1 = std::abs(pfp.XYZ[1][1] - tpc.MaxY());
    if (std::abs(pfp.XYZ[1][2] - tpc.MinZ())<mindis1) mindis1 = std::abs(pfp.XYZ[1][2] - tpc.MinZ());
    if (std::abs(pfp.XYZ[1][2] - tpc.MaxZ())<mindis1) mindis1 = std::abs(pfp.XYZ[1][2] - tpc.MaxZ());
    //std::cout<<pfp.XYZ[0][1]<<" "<<pfp.XYZ[0][2]<<" "<<pfp.XYZ[1][1]<<" "<<pfp.XYZ[1][2]<<" "<<tpc.MinY()<<" "<<tpc.MaxY()<<" "<<tpc.MinZ()<<" "<<tpc.MaxZ()<<" "<<mindis0<<" "<<mindis1<<" "<<mindis0+mindis1<<std::endl;
    slc.crt.cr_pfpyzmindis.push_back(mindis0+mindis1);

    if (slc.crt.cr_pfpxmin.back()<-2||
        slc.crt.cr_pfpxmax.back()>260||
        slc.crt.cr_pfpyzmindis.back()<30){
      pfp.CosmicScore = 1.;
    }
    else pfp.CosmicScore = 0;

    /*
    float minx = FLT_MAX;
    float maxx = FLT_MIN;

    for(auto& tjID : pfp.TjIDs) {
      Trajectory& tj = slc.allTraj[tjID - 1];
      TrajPoint& beginPoint = tj.Pts[tj.EndPt[0]];
      TrajPoint& endPoint = tj.Pts[tj.EndPt[1]];
      if (beginPoint.Pos[1]/slc.unitsPerTick<minx){
        minx = beginPoint.Pos[1]/slc.unitsPerTick;
      }
      if (beginPoint.Pos[1]/slc.unitsPerTick>maxx){
        maxx = beginPoint.Pos[1]/slc.unitsPerTick;
      }
      if (endPoint.Pos[1]/slc.unitsPerTick<minx){
        minx = endPoint.Pos[1]/slc.unitsPerTick;
      }
      if (endPoint.Pos[1]/slc.unitsPerTick>maxx){
        maxx = endPoint.Pos[1]/slc.unitsPerTick;
      }
    } // tjID

    slc.crt.cr_pfpmintick.push_back(minx);
    slc.crt.cr_pfpmaxtick.push_back(maxx);
    */
//    std::cout<<pfp.mcpListIndex<<std::endl;
//    std::cout<<pfp.XYZ[0][0]<<" "<<pfp.XYZ[1][0]<<std::endl;

  }

  ////////////////////////////////////////////////
  int GetOrigin(TCSlice& slc, PFPStruct& pfp){

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;    
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;    

    std::map<int, float> omap; //<origin, energy>

    for(auto& tjID : pfp.TjIDs) {
      
      Trajectory& tj = slc.tjs[tjID - 1];
      for(auto& tp : tj.Pts) {
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          TCHit& slhit = slc.slHits[iht];
          auto& hit = (*evt.allHits)[slhit.allHitsIndex];
          raw::ChannelID_t channel = tcc.geom->PlaneWireToChannel((int)hit.WireID().Plane, (int)hit.WireID().Wire, (int)hit.WireID().TPC, (int)hit.WireID().Cryostat);
          double startTick = hit.PeakTime() - hit.RMS();
          double endTick = hit.PeakTime() + hit.RMS();
          // get a list of track IDEs that are close to this hit
          std::vector<sim::TrackIDE> tides;
          tides = bt_serv->ChannelToTrackIDEs(channel, startTick, endTick);          
          for(auto itide = tides.begin(); itide != tides.end(); ++itide) {
            omap[pi_serv->TrackIdToMCTruth_P(itide->trackID)->Origin()] += itide->energy;
          }
        }
      }
    }
    
    float maxe = -1;
    int origin = 0;
    for (auto & i : omap){
      if (i.second > maxe){
        maxe = i.second;
        origin = i.first;
      }
    }
    return origin;
  }

  ////////////////////////////////////////////////
  void ClearCRInfo(TCSlice& slc){
    slc.crt.cr_origin.clear();
    slc.crt.cr_pfpxmin.clear();
    slc.crt.cr_pfpxmax.clear();
    slc.crt.cr_pfpyzmindis.clear();
  }
}

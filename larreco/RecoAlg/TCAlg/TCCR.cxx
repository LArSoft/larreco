#include "larreco/RecoAlg/TCAlg/TCCR.h"

namespace tca {

  ////////////////////////////////////////////////
  void SaveCRInfo(TjStuff& tjs, MatchStruct& ms, bool prt){

    tjs.crt.cr_pfpxmin.push_back(std::min(ms.XYZ[0][0], ms.XYZ[1][0]));
    tjs.crt.cr_pfpxmax.push_back(std::max(ms.XYZ[0][0], ms.XYZ[1][0]));
    /*
    float minx = FLT_MAX;
    float maxx = FLT_MIN;

    for(auto& tjID : ms.TjIDs) {
      Trajectory& tj = tjs.allTraj[tjID - 1];
      TrajPoint& beginPoint = tj.Pts[tj.EndPt[0]];
      TrajPoint& endPoint = tj.Pts[tj.EndPt[1]];
      if (beginPoint.Pos[1]/tjs.UnitsPerTick<minx){
        minx = beginPoint.Pos[1]/tjs.UnitsPerTick;
      }
      if (beginPoint.Pos[1]/tjs.UnitsPerTick>maxx){
        maxx = beginPoint.Pos[1]/tjs.UnitsPerTick;
      }
      if (endPoint.Pos[1]/tjs.UnitsPerTick<minx){
        minx = endPoint.Pos[1]/tjs.UnitsPerTick;
      }
      if (endPoint.Pos[1]/tjs.UnitsPerTick>maxx){
        maxx = endPoint.Pos[1]/tjs.UnitsPerTick;
      }
    } // tjID

    tjs.crt.cr_pfpmintick.push_back(minx);
    tjs.crt.cr_pfpmaxtick.push_back(maxx);
    */
//    std::cout<<ms.MCPartListIndex<<std::endl;
//    std::cout<<ms.XYZ[0][0]<<" "<<ms.XYZ[1][0]<<std::endl;

  }

  ////////////////////////////////////////////////
  void ClearCRInfo(TjStuff& tjs){
    tjs.crt.cr_pfpxmin.clear();
    tjs.crt.cr_pfpxmax.clear();
    tjs.crt.cr_origin.clear();
  }
}

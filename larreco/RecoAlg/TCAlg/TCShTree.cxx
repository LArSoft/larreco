#include "larreco/RecoAlg/TCAlg/TCShTree.h"

namespace tca {

  void SaveTjInfo(TjStuff& tjs,  const CTP_t& inCTP, std::vector<std::vector<int>>& tjList,
                  std::string stageName) {
    
    int stageNum = GetStageNum(tjs.stv, stageName);

    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      Trajectory& tj1 = tjs.allTraj[it1];
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) continue;

      SaveTjInfoStuff(tjs, inCTP, tj1,  stageNum, stageName);

      int trajID = tj1.ID;
      bool inShower = false;

      for (size_t l1 = 0; l1 < tjList.size(); ++l1) {
        if (inShower) break;
        for (size_t l2 = 0; l2 < tjList[l1].size(); ++l2) {

          if (trajID == tjList[l1][l2]) {
            tjs.stv.ShowerID.back() = l1;
            inShower = true;
            break;
          }
        } // end list loop 2 
      } // end list loop 1
    } // end tjs loop
    // add meaningless envelope to list for counting purposes
    // envelopes are defined once DefineShower is called
    // fill four times, one for each side of polygon
    for (int i = 0; i < 8; i++) {
      tjs.stv.Envelope.push_back(-999);
      tjs.stv.EnvStage.push_back(stageNum);
      tjs.stv.EnvPlane.push_back(-1);
      tjs.stv.EnvShowerID.push_back(-1);
    }

  } // SaveTjInfo (tjlist) 

  void SaveTjInfo(TjStuff& tjs,  const CTP_t& inCTP, const unsigned short& cotIndex,
                  std::string stageName) {

    int stageNum = GetStageNum(tjs.stv, stageName);

    ShowerStruct& ss = tjs.cots[cotIndex];

    bool noMatch = true;

    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {

      Trajectory& tj1 = tjs.allTraj[it1];

      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) {
	continue;
      }

      int trajID = tj1.ID;

      // check if this tj has already been added to the list
      // for this particular stage and plane
      int tjIndex = -1;
      bool isShowerTj = false;
      for (size_t i = 0; i < tjs.stv.TjID.size(); ++i) {
	if (tjs.stv.StageNum.at(i) != (int)stageNum) continue;
	if (tjs.stv.PlaneNum.at(i) != (short)DecodeCTP(ss.CTP).Plane) continue;

	if (tjs.stv.TjID.at(i) == trajID) {
	  tjIndex = i;
	  if (tjs.stv.IsShowerTj.at(tjIndex) == 1) isShowerTj = true;
	  //beenDoneBefore = true;
	  break;
	}
      }

      if (isShowerTj) continue;
      if (tjIndex == -1) SaveTjInfoStuff(tjs, inCTP, tj1, stageNum, stageName);

      for (size_t i = 0; i < ss.TjIDs.size(); ++i) {
        if (trajID == ss.TjIDs[i]) {
	  noMatch = false;
          if (tjIndex == -1) tjs.stv.ShowerID.back() = cotIndex;
	  else tjs.stv.ShowerID.at(tjIndex) = cotIndex;
        }

	if (it1 == (ss.ShowerTjID - 1)) tjs.stv.IsShowerTj.back() = 1;

        // check if tj is shower parent. if so, add to ttree
        // and mark parent flag	
        if (trajID == ss.ParentID) {
          if (tjIndex == -1) {
	    tjs.stv.ShowerID.back() = cotIndex;
	    tjs.stv.IsShowerParent.back() = 1;
	  }
	  else {
	    tjs.stv.ShowerID.at(tjIndex) = cotIndex;
	    tjs.stv.IsShowerParent.at(tjIndex) = 1;
	  }
	  break;

        }
      } // ss TjID loop
    } // end tjs loop
    
    if (noMatch) return;

    // add envelope information to showertreevars 
    geo::PlaneID iPlnID = DecodeCTP(ss.CTP);

    for (int i = 0; i < 8; i++) {
      tjs.stv.EnvStage.push_back(stageNum);
      tjs.stv.EnvPlane.push_back(iPlnID.Plane);
      tjs.stv.EnvShowerID.push_back(cotIndex);
    }

    tjs.stv.Envelope.push_back(ss.Envelope[0][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[0][1]/tjs.UnitsPerTick);
    tjs.stv.Envelope.push_back(ss.Envelope[1][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[1][1]/tjs.UnitsPerTick);
    tjs.stv.Envelope.push_back(ss.Envelope[2][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[2][1]/tjs.UnitsPerTick);
    tjs.stv.Envelope.push_back(ss.Envelope[3][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[3][1]/tjs.UnitsPerTick);

  } // SaveTjInfo (cots)

  void SaveTjInfoStuff(TjStuff& tjs,  const CTP_t& inCTP, Trajectory& tj, int stageNum, std::string stageName) {

    TrajPoint& beginPoint = tj.Pts[tj.EndPt[0]];
    TrajPoint& endPoint = tj.Pts[tj.EndPt[1]];

    tjs.stv.BeginWir.push_back(std::nearbyint(beginPoint.Pos[0]));
    tjs.stv.BeginTim.push_back(std::nearbyint(beginPoint.Pos[1]/tjs.UnitsPerTick));
    tjs.stv.BeginAng.push_back(beginPoint.Ang);
    tjs.stv.BeginChg.push_back(beginPoint.Chg);
    tjs.stv.BeginVtx.push_back(tj.VtxID[0]);

    tjs.stv.EndWir.push_back(std::nearbyint(endPoint.Pos[0]));
    tjs.stv.EndTim.push_back(std::nearbyint(endPoint.Pos[1]/tjs.UnitsPerTick));
    tjs.stv.EndAng.push_back(endPoint.Ang);
    tjs.stv.EndChg.push_back(endPoint.Chg);
    tjs.stv.EndVtx.push_back(tj.VtxID[1]);

    tjs.stv.MCSMom.push_back(tj.MCSMom);
    tjs.stv.TjID.push_back(tj.ID);
    tjs.stv.IsShowerTj.push_back(-1);

    tjs.stv.ShowerID.push_back(-1);
    tjs.stv.IsShowerParent.push_back(-1);
    tjs.stv.StageNum.push_back(stageNum);
    tjs.stv.nStages = stageNum;
    geo::PlaneID iPlnID = DecodeCTP(tj.CTP);
    tjs.stv.PlaneNum.push_back(iPlnID.Plane);

    tjs.stv.nPlanes = tjs.NumPlanes;

  } // SaveTjInfoStuff  

  ////////////////////////////////////////////////
  void SaveAllCots(TjStuff& tjs, const CTP_t& inCTP, std::string someText)
  {
    if(!tjs.SaveShowerTree) return;
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      auto& ss = tjs.cots[cotIndex];
      if (ss.CTP != inCTP) continue;
      if(ss.ID == 0) continue;
      SaveTjInfo(tjs, ss.CTP, cotIndex, someText);
    } // cotIndex
  } // SaveAllCots

  int GetStageNum(ShowerTreeVars& stv, std::string stageName) {
    int stageNum;
    bool existingStage = false;
    for (unsigned short i = 0; i < stv.StageName.size(); ++i) {
      if (stv.StageName.at(i) == stageName) {
        existingStage = true;
        stageNum = i+1;
      }
    }

    if (!existingStage) {
      stv.StageName.push_back(stageName);
      stageNum = stv.StageName.size();
    }

    return stageNum;
  }


  void ClearShowerTree(ShowerTreeVars& stv) {
    stv.BeginWir.clear();
    stv.BeginTim.clear();
    stv.BeginAng.clear();
    stv.BeginChg.clear();
    stv.BeginVtx.clear();
    stv.EndWir.clear();
    stv.EndTim.clear();
    stv.EndAng.clear();
    stv.EndChg.clear();
    stv.EndVtx.clear();
    stv.MCSMom.clear();
    stv.PlaneNum.clear();
    stv.TjID.clear();
    stv.IsShowerTj.clear();
    stv.ShowerID.clear();
    stv.IsShowerParent.clear();
    stv.StageNum.clear();
    stv.Envelope.clear();
    stv.EnvPlane.clear();
    stv.EnvStage.clear();
    stv.EnvShowerID.clear();

    return;

  } // ClearShowerTree

} // namespace tca

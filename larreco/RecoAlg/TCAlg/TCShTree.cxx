#include "larreco/RecoAlg/TCAlg/TCShTree.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h"

#include <array>
#include <bitset>
#include <math.h>
#include <stddef.h>
#include <string>
#include <vector>

namespace tca {

  void SaveTjInfo(TCSlice& slc, std::vector<std::vector<int>>& tjList, std::string stageName)
  {
    if (!tcc.modes[kSaveShowerTree]) return;
    if (tjList.empty()) return;
    int stageNum = GetStageNum(stv, stageName);

    // get the CTP from the first tj
    CTP_t inCTP = slc.tjs[tjList[0][0] - 1].CTP;
    for (unsigned short it1 = 0; it1 < slc.tjs.size(); ++it1) {
      Trajectory& tj1 = slc.tjs[it1];
      if (tj1.CTP != inCTP) continue;
      if (tj1.AlgMod[kKilled]) continue;

      SaveTjInfoStuff(slc, tj1, stageNum, stageName);

      int trajID = tj1.ID;
      bool inShower = false;

      for (size_t l1 = 0; l1 < tjList.size(); ++l1) {
        if (inShower) break;
        for (size_t l2 = 0; l2 < tjList[l1].size(); ++l2) {

          if (trajID == tjList[l1][l2]) {
            stv.ShowerID.back() = l1;
            inShower = true;
            break;
          }
        } // end list loop 2
      }   // end list loop 1
    }     // end tjs loop
    // add meaningless envelope to list for counting purposes
    // envelopes are defined once DefineShower is called
    // fill four times, one for each side of polygon
    for (int i = 0; i < 8; i++) {
      stv.Envelope.push_back(-999);
      stv.EnvStage.push_back(stageNum);
      stv.EnvPlane.push_back(-1);
      stv.EnvShowerID.push_back(-1);
    }

  } // SaveTjInfo (tjlist)

  void SaveTjInfo(TCSlice& slc, const ShowerStruct& ss, std::string stageName)
  {
    if (!tcc.modes[kSaveShowerTree]) return;
    int stageNum = GetStageNum(stv, stageName);

    // killed shower?
    if (ss.ID == 0) return;

    bool noMatch = true;

    for (unsigned short it1 = 0; it1 < slc.tjs.size(); ++it1) {

      Trajectory& tj1 = slc.tjs[it1];

      if (tj1.AlgMod[kKilled]) continue;

      int trajID = tj1.ID;

      // check if this tj has already been added to the list
      // for this particular stage and plane
      int tjIndex = -1;
      bool isShowerTj = false;
      for (size_t i = 0; i < stv.TjID.size(); ++i) {
        if (stv.StageNum.at(i) != (int)stageNum) continue;
        if (stv.PlaneNum.at(i) != (short)DecodeCTP(ss.CTP).Plane) continue;

        if (stv.TjID.at(i) == trajID) {
          tjIndex = i;
          if (stv.IsShowerTj.at(tjIndex) == 1) isShowerTj = true;
          //beenDoneBefore = true;
          break;
        }
      }

      if (isShowerTj) continue;
      if (tjIndex == -1) SaveTjInfoStuff(slc, tj1, stageNum, stageName);

      for (size_t i = 0; i < ss.TjIDs.size(); ++i) {
        if (trajID == ss.TjIDs[i]) {
          noMatch = false;
          if (tjIndex == -1)
            stv.ShowerID.back() = ss.ID;
          else
            stv.ShowerID.at(tjIndex) = ss.ID;
        }

        if (it1 == (ss.ShowerTjID - 1))
          stv.IsShowerTj.back() = 1;
        else if (tj1.AlgMod[kShowerTj])
          stv.IsShowerTj.back() = 1; // this is a better check
        // check if tj is shower parent. if so, add to ttree
        // and mark parent flag
        if (trajID == ss.ParentID) {
          if (tjIndex == -1) {
            stv.ShowerID.back() = ss.ID;
            stv.IsShowerParent.back() = 1;
          }
          else {
            stv.ShowerID.at(tjIndex) = ss.ID;
            stv.IsShowerParent.at(tjIndex) = 1;
          }
          break;
        }
      } // ss TjID loop
    }   // end tjs loop

    if (noMatch) return;

    // add envelope information to showertreevars
    geo::PlaneID iPlnID = DecodeCTP(ss.CTP);

    for (int i = 0; i < 8; i++) {
      stv.EnvStage.push_back(stageNum);
      stv.EnvPlane.push_back(iPlnID.Plane);
      stv.EnvShowerID.push_back(ss.ID);
    }

    stv.Envelope.push_back(ss.Envelope[0][0]);
    stv.Envelope.push_back(ss.Envelope[0][1] / tcc.unitsPerTick);
    stv.Envelope.push_back(ss.Envelope[1][0]);
    stv.Envelope.push_back(ss.Envelope[1][1] / tcc.unitsPerTick);
    stv.Envelope.push_back(ss.Envelope[2][0]);
    stv.Envelope.push_back(ss.Envelope[2][1] / tcc.unitsPerTick);
    stv.Envelope.push_back(ss.Envelope[3][0]);
    stv.Envelope.push_back(ss.Envelope[3][1] / tcc.unitsPerTick);

  } // SaveTjInfo (cots)

  void SaveTjInfoStuff(TCSlice& slc, Trajectory& tj, int stageNum, std::string stageName)
  {
    if (!tcc.modes[kSaveShowerTree]) return;

    TrajPoint& beginPoint = tj.Pts[tj.EndPt[0]];
    TrajPoint& endPoint = tj.Pts[tj.EndPt[1]];

    stv.BeginWir.push_back(std::nearbyint(beginPoint.Pos[0]));
    stv.BeginTim.push_back(std::nearbyint(beginPoint.Pos[1] / tcc.unitsPerTick));
    stv.BeginAng.push_back(beginPoint.Ang);
    stv.BeginChg.push_back(beginPoint.Chg);
    stv.BeginVtx.push_back(tj.VtxID[0]);

    stv.EndWir.push_back(std::nearbyint(endPoint.Pos[0]));
    stv.EndTim.push_back(std::nearbyint(endPoint.Pos[1] / tcc.unitsPerTick));
    stv.EndAng.push_back(endPoint.Ang);
    stv.EndChg.push_back(endPoint.Chg);
    stv.EndVtx.push_back(tj.VtxID[1]);

    stv.MCSMom.push_back(tj.MCSMom);
    stv.TjID.push_back(tj.ID);
    stv.IsShowerTj.push_back(-1);

    stv.ShowerID.push_back(-1);
    stv.IsShowerParent.push_back(-1);
    stv.StageNum.push_back(stageNum);
    stv.nStages = stageNum;
    geo::PlaneID iPlnID = DecodeCTP(tj.CTP);
    stv.PlaneNum.push_back(iPlnID.Plane);

    stv.nPlanes = slc.nPlanes;

  } // SaveTjInfoStuff

  ////////////////////////////////////////////////
  void SaveAllCots(TCSlice& slc, const CTP_t& inCTP, std::string someText)
  {
    if (!tcc.modes[kSaveShowerTree]) return;
    for (unsigned short cotIndex = 0; cotIndex < slc.cots.size(); ++cotIndex) {
      auto& ss = slc.cots[cotIndex];
      if (ss.CTP != inCTP) continue;
      if (ss.ID == 0) continue;
      SaveTjInfo(slc, ss, someText);
    } // cotIndex
  }   // SaveAllCots

  void SaveAllCots(TCSlice& slc, std::string someText)
  {
    if (!tcc.modes[kSaveShowerTree]) return;
    for (unsigned short cotIndex = 0; cotIndex < slc.cots.size(); ++cotIndex) {
      auto& ss = slc.cots[cotIndex];
      if (ss.ID == 0) continue;
      SaveTjInfo(slc, ss, someText);
    } // cotIndex
  }

  int GetStageNum(ShowerTreeVars& stv, std::string stageName)
  {
    int stageNum;
    bool existingStage = false;
    for (unsigned short i = 0; i < stv.StageName.size(); ++i) {
      if (stv.StageName.at(i) == stageName) {
        existingStage = true;
        stageNum = i + 1;
      }
    }

    if (!existingStage) {
      stv.StageName.push_back(stageName);
      stageNum = stv.StageName.size();
    }

    return stageNum;
  }

  void ClearShowerTree(ShowerTreeVars& stv)
  {
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

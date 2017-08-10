#include "larreco/RecoAlg/TCAlg/TCShTree.h"

namespace tca {

  void SaveTjInfo(TjStuff& tjs,  const CTP_t& inCTP, std::vector<std::vector<int>>& tjList,
                  unsigned int stageNum) {

    //    std::cout << "RUNNING FOR STAGE NUM: " << stageNum << std::endl;

    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      Trajectory& tj1 = tjs.allTraj[it1];
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) continue;


      SaveTjInfoStuff(tjs, inCTP, tj1,  stageNum);

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
                  unsigned int stageNum) {

    // std::cout << "RUNNING FOR STAGE NUM: " << stageNum << std::endl;

    ShowerStruct& ss = tjs.cots[cotIndex];

    bool noMatch = true;

    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {

      //std::cout << tjs.allTraj.size() << " "  << it1 << std::endl << std::flush;
   
      // right now I exclude the shower trajectory
      //      if (it1 == (ss.ShowerTjID - 1)) continue;

      Trajectory& tj1 = tjs.allTraj[it1];

      //std::cout << "HERE" << std::endl << std::flush;

      //      std::cout << tj1.CTP << std::endl << std::flush;
      //      std::cout << tj1.AlgMod[kKilled] << std::endl << std::flush;
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) {
	//std::cout << "INSIDE KILLED LOOP" << std::endl << std::flush;
	continue;
      }

      //      std::cout << "HERE2" << std::endl << std::flush;

      int trajID = tj1.ID;

      //      std::cout << "HERE3" << std::endl << std::flush;

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

      //      std::cout << "TJINDEX: " << tjIndex << std::endl << std::flush;
      //      if (tjIndex != -1) std::cout << "is shower tree? " << tjs.stv.IsShowerTj.at(tjIndex) << std::endl << std::flush;

      //      if (isShowerTj) std::cout << "THIS IS A SHOWER TRAJECTORY" << std::endl << std::flush;

      if (isShowerTj) continue;
      if (tjIndex == -1) SaveTjInfoStuff(tjs, inCTP, tj1, stageNum);

      for (size_t i = 0; i < ss.TjIDs.size(); ++i) {
        if (trajID == ss.TjIDs[i]) {
	  //	  std::cout << "TJID MATCH\n" << std::flush;
	  noMatch = false;
          if (tjIndex == -1) tjs.stv.ShowerID.back() = cotIndex;
	  else tjs.stv.ShowerID.at(tjIndex) = cotIndex;
	  //break;
        }

	if (it1 == (ss.ShowerTjID - 1)) tjs.stv.IsShowerTj.back() = 1;
        // check if tj is shower parent. if so, add to ttree
        // and mark parent flag
	
        if (trajID == ss.ParentID) {

	  //	  std::cout << "FOUND A PARENT!!" << std::endl;
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

    //    std::cout << "OUTSIDE LOOP" << std::endl << std::flush;
    // add envelope information to showertreevars 
    geo::PlaneID iPlnID = DecodeCTP(ss.CTP);

    //    std::cout << "HAVE PLANE" << std::endl << std::flush;
    for (int i = 0; i < 8; i++) {
      tjs.stv.EnvStage.push_back(stageNum);
      tjs.stv.EnvPlane.push_back(iPlnID.Plane);
      tjs.stv.EnvShowerID.push_back(cotIndex);
    }

    //    std::cout << "ADDED ENVELOPE VARIABLES" << std::endl << std::flush;

    tjs.stv.Envelope.push_back(ss.Envelope[0][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[0][1]/tjs.UnitsPerTick);
    tjs.stv.Envelope.push_back(ss.Envelope[1][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[1][1]/tjs.UnitsPerTick);
    tjs.stv.Envelope.push_back(ss.Envelope[2][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[2][1]/tjs.UnitsPerTick);
    tjs.stv.Envelope.push_back(ss.Envelope[3][0]);
    tjs.stv.Envelope.push_back(ss.Envelope[3][1]/tjs.UnitsPerTick);

    //    std::cout << "ADDED MORE ENVELOPE VARIABLES" << std::endl << std::flush;
  } // SaveTjInfo (cots)

  void SaveTjInfoStuff(TjStuff& tjs,  const CTP_t& inCTP, Trajectory& tj,  unsigned int stageNum) {

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

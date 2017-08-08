#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  const std::vector<std::string> AlgBitNames {
    "MaskHits",
    "MaskBadTPs",
    "CTKink",
    "CTStepChk",
    "TryNextPass",
    "RvPrp",
    "CHMH",
    "Split",
    "Comp3DVx",
    "Comp3DVxIG",
    "HED",
    "HamVx",
    "HamVx2",
    "JunkTj",
    "Killed",
    "Merge",
    "TEP",
    "CHMEH",
    "FillGap",
    "Ghost",
    "ChkInTraj",
    "StopBadFits",
    "FixBegin",
    "FixEnd",
    "UUH",
    "MisdVxTj",
    "VtxTj",
    "RefVtx",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "FTBRvProp",
    "StopAtTj",
    "M3D",
    "Mat3DMerge",
    "TjHiVx3Score",
    "VtxHitsSwap",
    "SplitHiChgHits",
    "InShower",
    "ShowerTj",
    "ShwrParent",
    "MergeOverlap",
    "MergeSubShowers",
    "MergeNrShowers"
  };

  const std::vector<std::string> StopFlagNames {
    "Signal",
    "AtKink",
    "AtVtx",
    "Bragg",
    "AtTj"
  };
  
  const std::vector<std::string> VtxBitNames {
    "VtxTrjTried",
    "Fixed",
    "OnDeadWire",
    "VtxRefined",
    "HiVx3Score",
    "VtxTruMatch",
    "VtxMerged"
  } ;
  
  geo::PlaneID DecodeCTP(CTP_t CTP) {
    geo::PlaneID tmp;
    tmp.Cryostat = CTP / Cpad;
    tmp.TPC = (CTP - tmp.Cryostat * Cpad) / Tpad;
    tmp.Plane = (CTP % 10);
    return tmp;
  }
  
} // namespace tca


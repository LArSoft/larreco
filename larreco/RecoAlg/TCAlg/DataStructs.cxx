#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  const std::vector<std::string> AlgBitNames {
    "MaskHits",
    "MaskBadTPs",
    "CTKink",
    "CTStepChk",
    "TryNextPass",
    "RevProp",
    "CHMH",
    "SplitTraj",
    "Comp3DVx",
    "Comp3DVxIG",
    "HED",
    "HamVx",
    "HamVx2",
    "JunkTj",
    "Killed",
    "EndMerge",
    "TrimEndPts",
    "CHMEH",
    "FillGap",
    "Ghost",
    "ChkInTraj",
    "StopBadFits",
    "FixBegin",
    "FixEnd",
    "UseUnusedHits",
    "VtxTj",
    "RefVtx",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkAllStop",
    "FTBRevProp",
    "StopAtTj",
    "Match3D",
    "VtxHitsSwap",
    "SplitHiChgHits",
    "InShower",
    "ShowerParent",
    "ShowerTj"
  };

  const std::vector<std::string> StopFlagNames {
    "Signal",
    "AtKink",
    "AtVtx",
    "Bragg",
    "AtTj"
  };
  
  const std::vector<std::string> VtxBitNames {
    "Fixed",
    "VtxTrjTried",
    "OnDeadWire",
    "VtxRefined",
    "NiceVtx",
    "kInShower"
  } ;
  
  geo::PlaneID DecodeCTP(CTP_t CTP) {
    geo::PlaneID tmp;
    tmp.Cryostat = CTP / Cpad;
    tmp.TPC = (CTP - tmp.Cryostat * Cpad) / Tpad;
    tmp.Plane = (CTP % 10);
    return tmp;
  }
  
} // namespace tca


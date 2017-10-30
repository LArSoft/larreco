#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  const std::vector<std::string> AlgBitNames {
    "MaskHits",
    "MaskBadTPs",
    "Michel",
    "DeltaRay",
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
    "MrgGhost",
    "ChkInTraj",
    "StopBadFits",
    "FixBegin",
    "BeginChg",
    "FixEnd",
    "UUH",
    "MisdVxTj",
    "VtxTj",
    "ChkVxTj",
    "RefVtx",
    "VxMerge",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkStopEP",
    "FTBRvProp",
    "StopAtTj",
    "Mat3D",
    "Mat3DMerge",
    "TjHiVx3Score",
    "VtxHitsSwap",
    "SplitHiChgHits",
    "InShower",
    "ShowerTj",
    "ShwrParent",
    "MergeOverlap",
    "MergeSubShowers",
    "MergeNrShowers",
    "MergeShChain",
    "SplitTjCVx",
    "SetDir"
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
    auto const cryo = (CTP / Cpad);
    return geo::PlaneID(
      /* Cryostat */ cryo,
           /* TPC */ (CTP - cryo * Cpad) / Tpad,
         /* Plane */ (CTP % 10)
      );
  }
  
} // namespace tca


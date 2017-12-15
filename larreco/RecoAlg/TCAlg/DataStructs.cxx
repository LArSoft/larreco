#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  const std::vector<std::string> AlgBitNames {
    "HitsOrdered",
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
    "JunkVx",
    "JunkTj",
    "Killed",
    "Merge",
    "MergeChain",
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
    "Photon",
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
    "KillInShowerVx",
    "ShowerTj",
    "ShwrParent",
    "ChkShwrParEnd",  // Ensure that the end of a shower parent already inside a shower has an end near a shower end
    "KillShwrNuPFP",  // Kill neutrino PFP particles with a vertex inside a shower
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


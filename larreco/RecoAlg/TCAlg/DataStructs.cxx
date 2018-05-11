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
    "CHMUH",
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
    "TEP",
    "CHMEH",
    "FillGap",
    "Ghost",
    "MrgGhost",
    "ChkInTraj",
    "StopBadFits",
    "FixBegin",
    "FTBChg",
    "BeginChg",
    "FixEnd",
    "UUH",
    "MisdVxTj",
    "VtxTj",
    "ChkVxTj",
    "Photon",
    "NoFitToVx",
    "VxMerge",
    "VxNeutral",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkStopEP",
    "ChkChgAsym",
    "FTBRvProp",
    "StopAtTj",
    "Mat3D",
    "Mat3DMerge",
    "Split3DKink",
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
    "AtTj",
    "OutFV"
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


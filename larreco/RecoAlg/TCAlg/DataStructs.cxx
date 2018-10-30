#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  
  TCEvent evt;
  TCConfig tcc;
  std::vector<TjForecast> tjfs;
  ShowerTreeVars stv;
  // vector of hits, tjs, etc in each slice
  std::vector<TCSlice> slices;
  //    TruthMatcher tm{tjs};

  const std::vector<std::string> AlgBitNames {
    "Stiff",
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
    "VtxTj",
    "ChkVxTj",
    "MisdVxTj",
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
    "M3DVxTj",
    "Mat3DMerge",
    "Split3DKink",
    "TjHiVx3Score",
    "VtxHitsSwap",
    "SplitHiChgHits",
    "ShowerLike",
    "KillInShowerVx",
    "ShowerTj",
    "ShwrParent",
    "MergeOverlap",
    "MergeSubShowers",
    "MergeSubShowersTj",
    "MergeNrShowers",
    "MergeShChain",
    "CompleteShower",
    "SplitTjCVx",
    "NewStpCuts",
    "NewVtxCuts"
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
    "VtxMerged",
    "VtxIndPlnNoChg"
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


#include "larreco/RecoAlg/TCAlg/DataStructs.h"

#include <vector>

namespace tca {

  TCEvent evt;
  TCConfig tcc;
  std::vector<TjForecast> tjfs;
  ShowerTreeVars stv;
  // vector of hits, tjs, etc in each slice
  std::vector<TCSlice> slices;
  std::vector<TrajPoint> seeds;

  const std::vector<std::string> AlgBitNames {
    "FillGaps3D",
    "Kink3D",
    "TEP3D",
    "Junk3D",
    "RTPs3D",
    "SmallAng3D",
    "KillBadPts3D",
    "Mat3D",
    "MaskHits",
    "MaskBadTPs",
    "Michel",
    "DeltaRay",
    "CTStepChk",
    "RvPrp",
    "CHMUH",
    "Split",
    "Comp3DVx",
    "Comp3DVxIG",
    "DecayVx",
    "HamVx",
    "HamVx2",
    "JunkVx",
    "ShrtLong2V",
    "JunkTj",
    "Killed",
    "Merge",
    "LastEndMerge",
    "TEP",
    "EndPtFit",
    "THiQEP",
    "CHMEH",
    "FillGaps",
    "UseGhostHits",
    "MrgGhost",
    "ChkInTraj",
    "StopBadFits",
    "FixBegin",
    "FTBChg",
    "BeginChg",
    "BraggSplit",
    "UUH",
    "VtxTj",
    "ChkVxTj",
    "Photon",
    "HaloTj",
    "NoFitToVx",
    "VxMerge",
    "VxNeutral",
    "KillOrphan2V",
    "NoKinkChk",
    "ChkStop",
    "ChkStopEP",
    "ChkChgAsym",
    "FTBRvProp",
    "TjHiVx3Score",
    "VxEndSwap",
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
    "MakePFPTjs",
    "StopShort",
    "Reconcile2Vs",
    "FTBMod",
    "NewCuts"
  };

  const std::vector<std::string> EndFlagNames {
    "HitsAfterEnd",
    "EndKink",
    "EndBragg",
    "EndOutFV",
    "EndNoFitVx"
  };

  const std::vector<std::string> VtxBitNames {
    "VxTrjTried",
    "Fixed",
    "OnDeadWire",
    "HiVx3Score",
    "VxTruMatch",
    "VxMerged",
    "VxIndPlnNoChg",
    "VxEnvOK"
  };

  const std::vector<std::string> StrategyBitNames {
    "Normal",
    "StiffEl",
    "StiffMu",
    "Slowing"
  };

  geo::PlaneID DecodeCTP(CTP_t CTP) {
    auto const cryo = (CTP / Cpad);
    return geo::PlaneID(
      /* Cryostat */ cryo,
           /* TPC */ (CTP - cryo * Cpad) / Tpad,
         /* Plane */ (CTP % 10)
      );
  }

} // namespace tca

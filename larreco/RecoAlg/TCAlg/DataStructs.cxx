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

  const std::vector<std::string> AlgBitNames{"FillGaps3D",
                                             "Kink3D",
                                             "TEP3D",
                                             "Junk3D",
                                             "RTPs3D",
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
                                             "HamBragg",
                                             "HamVx",
                                             "HamVx2",
                                             "JunkVx",
                                             "JunkTj",
                                             "Killed",
                                             "Merge",
                                             "LastEndMerge",
                                             "TEP",
                                             "THCEP",   // TrimHiChgEndPoints
                                             "EndKink", // CheckEndKink
                                             "CHMEH",
                                             "FillGaps",
                                             "UseGhostHits",
                                             "MrgGhost",
                                             "ChkInTraj",
                                             "StopBadFits",
                                             "FixBegin",
                                             "FTBChg",
                                             "BeginChg",
                                             "FixEnd",
                                             "BraggSplit",
                                             "UUH",
                                             "VtxTj",
                                             "ChkVxTj",
                                             "Photon",
                                             "HaloTj",
                                             "NoFitToVx",
                                             "VxMerge",
                                             "VxNeutral",
                                             "NoKinkChk",
                                             "ChkStop",
                                             "ChkStopEP",
                                             "ChkChgAsym",
                                             "FTBRvProp",
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
                                             "MakePFPTjs",
                                             "StopShort",
                                             "Reconcile2Vs",
                                             "FTBMod"};

  const std::vector<std::string>
    EndFlagNames{"Signal", "AtKink", "AtVtx", "Bragg", "AtTj", "OutFV", "NoFitVx"};

  const std::vector<std::string> VtxBitNames{"VxTrjTried",
                                             "Fixed",
                                             "OnDeadWire",
                                             "HiVx3Score",
                                             "VxTruMatch",
                                             "VxMerged",
                                             "VxIndPlnNoChg",
                                             "VxEnvOK"};

  const std::vector<std::string> StrategyBitNames{"Normal", "StiffEl", "StiffMu", "Slowing"};

  geo::PlaneID DecodeCTP(CTP_t CTP)
  {
    auto const cryo = (CTP / Cpad);
    return geo::PlaneID(
      /* Cryostat */ cryo,
      /* TPC */ (CTP - cryo * Cpad) / Tpad,
      /* Plane */ (CTP % 10));
  }

} // namespace tca

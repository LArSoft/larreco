#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  const std::vector<std::string> AlgBitNames {
    "MaskHits",
    "UnMaskHits",
    "Kink",
    "CTKink",
    "CTStepChk",
    "TryNextPass",
    "RevProp",
    "ChkHiMultHits",
    "SplitTraj",
    "Comp3DVx",
    "HiEndDelta",
    "HammerVx",
    "HammerVx2",
    "JunkTj",
    "Killed",
    "StopAtVtx",
    "EndMerge",
    "TrimHits",
    "ChkHiMultEndHits",
    "FillGap",
    "UseGhostHits",
    "ChkInTraj",
    "FixEnd",
    "UseUnusedHits",
    "VtxTj",
    "RefineVtx",
    "MaskBadTPs",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkAllStop",
    "FTBRevProp"
  };
  
  const std::vector<std::string> VtxBitNames {
    "Fixed",
    "VtxTrjTried",
    "VtxRefined"
  } ;
} // namespace tca


#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  const std::vector<std::string> AlgBitNames {
    "MaskHits",
    "UnMaskHits",
    "CTKink",
    "CTStepChk",
    "TryNextPass",
    "RevProp",
    "CHMH",
    "SplitTraj",
    "Comp3DVx",
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
    "FixBegin",
    "FixEnd",
    "UseUnusedHits",
    "VtxTj",
    "RefVtx",
    "MBadTPs",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkAllStop",
    "FTBRevProp",
    "MHOK2",
    "FUH2",
    "StopAtTj",
  };

  const std::vector<std::string> StopFlagNames {
    "Signal",
    "AtKink",
    "AtVtx",
    "Bragg",
    "RvPrp",
    "AtTj"
  };
  
  const std::vector<std::string> VtxBitNames {
    "Fixed",
    "VtxTrjTried",
    "VtxRefined"
  } ;
} // namespace tca


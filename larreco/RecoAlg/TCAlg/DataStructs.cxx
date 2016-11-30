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
    "TrimHits",
    "CHMEH",
    "FillGap",
    "Ghost",
    "CIT",
    "FixEnd",
    "UUH",
    "VtxTj",
    "RefVtx",
    "MBadTPs",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkAllStop",
    "FTBRP"
  };

  const std::vector<std::string> StopFlagNames {
    "Signal",
    "AtKink",
    "AtVtx",
    "Bragg"
  };
  
  const std::vector<std::string> VtxBitNames {
    "Fixed",
    "VtxTrjTried",
    "VtxRefined"
  } ;
} // namespace tca


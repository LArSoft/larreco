#include "larreco/RecoAlg/TCAlg/TCCR.h"

namespace tca {

  ////////////////////////////////////////////////
  void SaveCRInfo(TjStuff& tjs, MatchStruct& ms, bool prt){

    tjs.crt.cr_pfpx0.push_back(std::min(ms.XYZ[0][0], ms.XYZ[1][0]));
    tjs.crt.cr_pfpx1.push_back(std::max(ms.XYZ[0][0], ms.XYZ[1][0]));

  }
}

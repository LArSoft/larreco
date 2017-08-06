#include "larreco/RecoAlg/TCAlg/TCCR.h"

namespace tca {

  ////////////////////////////////////////////////
  void SaveCRInfo(TjStuff& tjs, MatchStruct& ms, bool prt){

    tjs.crh.cr_pfpx0->Fill(std::min(ms.XYZ[0][0], ms.XYZ[1][0]));
    tjs.crh.cr_pfpx1->Fill(std::min(ms.XYZ[0][0], ms.XYZ[1][0]));

  }
}

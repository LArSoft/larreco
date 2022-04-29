#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larreco/RecoAlg/TCAlg/TCHist.h"

#include <cmath>
#include "TH1F.h"
#include "TProfile.h"

namespace  tca {

  void HistStuff::CreateHists(art::ServiceHandle<art::TFileService const>& tfs)
  {

    fTruT[0] = tfs->make<TH1F>("TruT_El","True KE (MeV) - Electrons", 100, 0, 1000);
    fTruT[1] = tfs->make<TH1F>("TruT_Mu","True KE (MeV) - Muons", 100, 0, 10000);
    fTruT[2] = tfs->make<TH1F>("TruT_Pi","True KE (MeV) - Pions", 100, 0, 1000);
    fTruT[3] = tfs->make<TH1F>("TruT_Ka","True KE (MeV) - Kaons", 100, 0, 10000);
    fTruT[4] = tfs->make<TH1F>("TruT_Pr","True KE (MeV) - Protons", 100, 0, 1000);

    fEff_T[0] = tfs->make<TProfile>("Eff_T_El","Eff vs T(MeV) - Electrons", 20, 0, 1000);
    fEff_T[1] = tfs->make<TProfile>("Eff_T_Mu","Eff vs T(MeV) - Muons", 20, 0, 10000);
    fEff_T[2] = tfs->make<TProfile>("Eff_T_Pi","Eff vs T(MeV) - Pions", 20, 0, 1000);
    fEff_T[3] = tfs->make<TProfile>("Eff_T_Ka","Eff vs T(MeV) - Kaons", 20, 0, 1000);
    fEff_T[4] = tfs->make<TProfile>("Eff_T_Pr","Eff vs T(MeV) - Protons", 20, 0, 1000);

    fPur_T[0] = tfs->make<TProfile>("Pur_T_El","Pur vs T(MeV) - Electrons", 20, 0, 1000);
    fPur_T[1] = tfs->make<TProfile>("Pur_T_Mu","Pur vs T(MeV) - Muons", 20, 0, 10000);
    fPur_T[2] = tfs->make<TProfile>("Pur_T_Pi","Pur vs T(MeV) - Pions", 20, 0, 1000);
    fPur_T[3] = tfs->make<TProfile>("Pur_T_Ka","Pur vs T(MeV) - Kaons", 20, 0, 1000);
    fPur_T[4] = tfs->make<TProfile>("Pur_T_Pr","Pur vs T(MeV) - Protons", 20, 0, 1000);

  }


} // tca

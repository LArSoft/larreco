#include "art_root_io/TFileService.h"
#include "larreco/RecoAlg/TCAlg/TCHist.h"

#include <math.h>
#include "TH1F.h"
#include "TH2.h"
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


/* TTree used to develop the TMVA showerParentReader
    fShowerParentSig = tfs->make<TTree>("shwr_parent_tree_sig", "shwr_parent_tree_sig");
    fShowerParentSig->Branch("fShEnergy", &fShEnergy, "fShEnergy/F");
    fShowerParentSig->Branch("fPfpEnergy", &fPfpEnergy, "fPfpEnergy/F");
    fShowerParentSig->Branch("fMCSMom", &fMCSMom, "fMCSMom/F");
    fShowerParentSig->Branch("fPfpLen", &fPfpLen, "fPfpLen/F");
    fShowerParentSig->Branch("fSep", &fSep, "fSep/F");
    fShowerParentSig->Branch("fDang1", &fDang1, "fDang1/F");
    fShowerParentSig->Branch("fDang2", &fDang2, "fDang2/F");
    fShowerParentSig->Branch("fChgFrac", &fChgFrac, "fChgFrac/F");
    fShowerParentSig->Branch("fAlong", &fAlong, "fAlong/F");
    fShowerParentSig->Branch("fTrans", &fTrans, "fTrans/F");
    fShowerParentSig->Branch("fInShwrProb", &fInShwrProb, "fInShwrProb/F");

    fShowerParentBkg = tfs->make<TTree>("shwr_parent_tree_bkg", "shwr_parent_tree_bkg");
    fShowerParentBkg->Branch("fShEnergy", &fShEnergy, "fShEnergy/F");
    fShowerParentBkg->Branch("fPfpEnergy", &fPfpEnergy, "fPfpEnergy/F");
    fShowerParentBkg->Branch("fMCSMom", &fMCSMom, "fMCSMom/F");
    fShowerParentBkg->Branch("fPfpLen", &fPfpLen, "fPfpLen/F");
    fShowerParentBkg->Branch("fSep", &fSep, "fSep/F");
    fShowerParentBkg->Branch("fDang1", &fDang1, "fDang1/F");
    fShowerParentBkg->Branch("fDang2", &fDang2, "fDang2/F");
    fShowerParentBkg->Branch("fChgFrac", &fChgFrac, "fChgFrac/F");
    fShowerParentBkg->Branch("fAlong", &fAlong, "fAlong/F");
    fShowerParentBkg->Branch("fTrans", &fTrans, "fTrans/F");
    fShowerParentBkg->Branch("fInShwrProb", &fInShwrProb, "fInShwrProb/F");
*/
  }


} // tca

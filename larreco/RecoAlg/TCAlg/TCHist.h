////////////////////////////////////////////////////////////////////////
//
//
// TCAlg debug struct
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGHISTSTRUCT_H
#define TRAJCLUSTERALGHISTSTRUCT_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

class TH1F;
class TH2F;
class TProfile;

namespace tca {

  struct HistStuff {
    void CreateHists(art::ServiceHandle<art::TFileService const>& tfs);

    // True kinetic energy (MeV)
    TH1F* fTruT[5];
    TProfile* fEff_T[5];
    TProfile* fPur_T[5];
  };
} // namespace tca

#endif // ifndef TRAJCLUSTERALGHISTSTRUCT_H

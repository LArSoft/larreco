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

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TH1F.h"
#include "TProfile.h"

namespace art { class TFileService; }
class TH1F;
class TH2F;
class TProfile;

namespace tca {

  struct HistStuff {
    void CreateHists(art::ServiceHandle<art::TFileService const>& tfs);

    // True kinetic energy (MeV)
    TH1F *fTruT[5];
    TProfile* fEff_T[5];
    TProfile* fPur_T[5];


/*
    // study electrons
    TH1F *fChgRMS[5];
    TH1F *fMomAsym[5];
    TH1F *fElectronLike[5];
    TH2F *fElectronLike_Len[5];

    TH1F *fChgToMeV[3];
    TProfile *fChgToMeV_Etru;
    TH2F *AlongTrans1;
    TH2F *AlongTrans5;
    TH2F *AlongTrans9;

    TH1F *fUnMatchedHitFrac;


    TH2F *fMCSMom_TruMom_e;
    TH2F *fMCSMom_TruMom_mu;
    TH2F *fMCSMom_TruMom_pi;
    TH2F *fMCSMom_TruMom_p;

    TH2F *fMCSMomEP_TruMom_e;

    // Reco-MC vertex position difference
    TH1F* fNuVtx_dx;
    TH1F* fNuVtx_dy;
    TH1F* fNuVtx_dz;

    // Vertex score for 2D vertices that are near the neutrino interaction vertex
    TH1F* fNuVx3Score;
    TH1F* fNuVx2Score;
    TH1F* fNuVx3ScoreDiff;
    TH1F* fVxTopoMat;
    TH1F* fVxTopoNoMat;
    // Vertex score for 2D and 3D vertices
    TH1F* fVx2Score;
    TH1F* fVx3Score;

    // Reco-MC stopping wire difference for different MC Particles
    TH1F* fdWire[5];
    // EP vs KE for different MC Particles
    TProfile* fEP_T[5];

    // fraction of TPs that have the kEnvNearTj bit set
    TProfile* fNearTj[5];

    // PFParticle PDGCode vs true PDG code
    TH2F* PDGCode_reco_true;
    TH1F* fPFPStartEnd;
    TH1F* fPFPStartdX[5];
    TH1F* fPFPStartdY[5];
    TH1F* fPFPStartdZ[5];
    TH1F* fPFPStartAngDiff[5];

    TH1F* fEff;
    TH1F* fPur;
*/
//    TTree* fShowerParentSig;
//    TTree* fShowerParentBkg;

//    float fShEnergy, fPfpEnergy, fMCSMom, fPfpLen, fSep, fDang1, fDang2, fChgFrac, fAlong, fTrans, fInShwrProb;

  };
} // namespace tca

#endif // ifndef TRAJCLUSTERALGHISTSTRUCT_H

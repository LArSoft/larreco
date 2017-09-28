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

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

namespace tca {
  
  struct HistStuff {
    void CreateHists(art::TFileService& tfs);
    
    TH1F *fUnMatchedHitFrac;
    
    // 
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
    
    // PFParticle PDGCode vs true PDG code
    TH2F* PDGCode_reco_true;

    //Test KF Fit
    TH1F* dXkf_match;
    TH1F* dXtc_match;
    TH1F* dYkf_match;
    TH1F* dYtc_match;
    TH1F* dZkf_match;
    TH1F* dZtc_match;
    TH1F* dUXkf_match;
    TH1F* dUXtc_match;
    TH1F* dUYkf_match;
    TH1F* dUYtc_match;
    TH1F* dUZkf_match;
    TH1F* dUZtc_match;
    TH1F* dXpull_match;
    TH1F* dYpull_match;
    TH1F* dZpull_match;
    TH1F* dUXpull_match;
    TH1F* dUYpull_match;
    TH1F* dUZpull_match;
    TH1F* covtrace_match;
    TH1F* hasfit_match;
    TH1F* nchi2_match;
    TH1F* nvalidpoints_match;
    //
    TH1F* dXkf_okid;
    TH1F* dXtc_okid;
    TH1F* dYkf_okid;
    TH1F* dYtc_okid;
    TH1F* dZkf_okid;
    TH1F* dZtc_okid;
    TH1F* dUXkf_okid;
    TH1F* dUXtc_okid;
    TH1F* dUYkf_okid;
    TH1F* dUYtc_okid;
    TH1F* dUZkf_okid;
    TH1F* dUZtc_okid;
    TH1F* dXpull_okid;
    TH1F* dYpull_okid;
    TH1F* dZpull_okid;
    TH1F* dUXpull_okid;
    TH1F* dUYpull_okid;
    TH1F* dUZpull_okid;
    TH1F* hasfit_okid;
    TH1F* nchi2_okid;
    //
    TH1F* dXkf_wrongid;
    TH1F* dXtc_wrongid;
    TH1F* dYkf_wrongid;
    TH1F* dYtc_wrongid;
    TH1F* dZkf_wrongid;
    TH1F* dZtc_wrongid;
    TH1F* dUXkf_wrongid;
    TH1F* dUXtc_wrongid;
    TH1F* dUYkf_wrongid;
    TH1F* dUYtc_wrongid;
    TH1F* dUZkf_wrongid;
    TH1F* dUZtc_wrongid;
    TH1F* hasfit_wrongid;
    TH1F* nchi2_wrongid;
    //
    TH1F* covtrace_nomatch;
    TH1F* hasfit_nomatch;
    TH1F* nchi2_nomatch;
    TH1F* nvalidpoints_nomatch;
  };
} // namespace tca

#endif // ifndef TRAJCLUSTERALGHISTSTRUCT_H

#include "larreco/RecoAlg/TCAlg/TCHist.h"

namespace  tca {
  
  void HistStuff::CreateHists(art::TFileService& tfs)
  {
    
    // True - Reco unmatched hit fraction
    fUnMatchedHitFrac = tfs.make<TH1F>("UnMatchedHitFrac","True-Reco Unmatched Hit Fraction", 200, 0, 0.2);
    
    // True - Reco vertex difference
    fNuVtx_dx = tfs.make<TH1F>("Vtx dx","Vtx dx",80,-10,10);
    fNuVtx_dy = tfs.make<TH1F>("Vtx dy","Vtx dy",80,-10,10);
    fNuVtx_dz = tfs.make<TH1F>("Vtx dz","Vtx dz",80,-10,10);
    
    // Score of all 2D and 3D vertices
    fVx2Score = tfs.make<TH1F>("Vx2_Score","Vx2 Score",80, 0, 100);
    fVx3Score = tfs.make<TH1F>("Vx3_Score","Vx2 Score",100, 0, 100);
    // Score of vertices that are matched to a true neutrino interaction vertex
    fNuVx3Score = tfs.make<TH1F>("NuVx3Score","Vx3Score - TruNuVx matched",80, 0, 80);
    fNuVx2Score = tfs.make<TH1F>("NuVx2Score","Vx2Score - TruNuVx matched",80, 0, 80);
    // Difference in score between all others and the one matched to the true neutrino vertex
    // Score(True matched) - Score(all others). Ideally the one matched to the true neutrino
    // has the highest score resulting in all entries begin negative in this histogram
    fNuVx3ScoreDiff = tfs.make<TH1F>("NuVx3ScoreDiff","Vx3Score - TruNuVx matched", 100, -100, 100);
    
    fVxTopoMat = tfs.make<TH1F>("VxTopoMat","VxTopo - True Vtx match",10, 0, 10);
    fVxTopoNoMat = tfs.make<TH1F>("VxTopoNoMat","VxTopo - No MC match",10, 0, 10);
    
    fdWire[0] = tfs.make<TH1F>("dWireEl","dWire - Electrons",21,-10,10);
    fdWire[1] = tfs.make<TH1F>("dWireMu","dWire - Muons",21,-10,10);
    fdWire[2] = tfs.make<TH1F>("dWirePi","dWire - Pions",21,-10,10);
    fdWire[3] = tfs.make<TH1F>("dWireKa","dWire - Kaons",21,-10,10);
    fdWire[4] = tfs.make<TH1F>("dWirePr","dWire - Protons",21,-10,10);
    
    fEP_T[0] = tfs.make<TProfile>("EP_T_El","EP vs T(MeV) - Electrons", 20, 0, 100);
    fEP_T[1] = tfs.make<TProfile>("EP_T_Mu","EP vs T(MeV) - Muons", 20, 0, 1000);
    fEP_T[2] = tfs.make<TProfile>("EP_T_Pi","EP vs T(MeV) - Pions", 20, 0, 1000);
    fEP_T[3] = tfs.make<TProfile>("EP_T_Ka","EP vs T(MeV) - Kaons", 20, 0, 1000);
    fEP_T[4] = tfs.make<TProfile>("EP_T_Pr","EP vs T(MeV) - Protons", 20, 0, 1000);
    
    
    fMCSMom_TruMom_e = tfs.make<TH2F>("MCSMom_TruMom_e","MCSMom vs Tru Mom electrons", 50, 0, 100, 50, 0, 1000);
    fMCSMom_TruMom_mu = tfs.make<TH2F>("MCSMom_TruMom_mu","MCSMom vs Tru Mom electrons", 50, 0, 1000, 50, 0, 1000);
    fMCSMom_TruMom_pi = tfs.make<TH2F>("MCSMom_TruMom_pi","MCSMom vs Tru Mom electrons", 50, 0, 1000, 50, 0, 1000);
    fMCSMom_TruMom_p = tfs.make<TH2F>("MCSMom_TruMom_p","MCSMom vs Tru Mom electrons", 50, 0, 1000, 50, 0, 1000);
    
    // Same as above but with good Efficiency * Purity
    fMCSMomEP_TruMom_e = tfs.make<TH2F>("MCSMomEP_TruMom_e","MCSMom vs Tru Mom electrons", 50, 0, 100, 50, 0, 1000);
    
    PDGCode_reco_true = tfs.make<TH2F>("PDGCode_reco_true", "PDGCode Reco vs True", 5, -0.5, 4.5, 5, -0.5, 4.5);

    art::TFileDirectory kffitdir = tfs.mkdir( "kffit" );
    dXkf_match = kffitdir.make<TH1F>("dXkf_match","dXkf_match",50,-10,10);
    dXtc_match = kffitdir.make<TH1F>("dXtc_match","dXtc_match",50,-10,10);
    dYkf_match = kffitdir.make<TH1F>("dYkf_match","dYkf_match",50,-10,10);
    dYtc_match = kffitdir.make<TH1F>("dYtc_match","dYtc_match",50,-10,10);
    dZkf_match = kffitdir.make<TH1F>("dZkf_match","dZkf_match",50,-10,10);
    dZtc_match = kffitdir.make<TH1F>("dZtc_match","dZtc_match",50,-10,10);
    dUXkf_match = kffitdir.make<TH1F>("dUXkf_match","dUXkf_match",50,-10,10);
    dUXtc_match = kffitdir.make<TH1F>("dUXtc_match","dUXtc_match",50,-10,10);
    dUYkf_match = kffitdir.make<TH1F>("dUYkf_match","dUYkf_match",50,-10,10);
    dUYtc_match = kffitdir.make<TH1F>("dUYtc_match","dUYtc_match",50,-10,10);
    dUZkf_match = kffitdir.make<TH1F>("dUZkf_match","dUZkf_match",50,-10,10);
    dUZtc_match = kffitdir.make<TH1F>("dUZtc_match","dUZtc_match",50,-10,10);
    dXpull_match = kffitdir.make<TH1F>("dXpull_match","dXpull_match",50,-10,10);
    dYpull_match = kffitdir.make<TH1F>("dYpull_match","dYpull_match",50,-10,10);
    dZpull_match = kffitdir.make<TH1F>("dZpull_match","dZpull_match",50,-10,10);
    dUXpull_match = kffitdir.make<TH1F>("dUXpull_match","dUXpull_match",50,-10,10);
    dUYpull_match = kffitdir.make<TH1F>("dUYpull_match","dUYpull_match",50,-10,10);
    dUZpull_match = kffitdir.make<TH1F>("dUZpull_match","dUZpull_match",50,-10,10);
    covtrace_match = kffitdir.make<TH1F>("covtrace_match", "covtrace_match", 100, 0, 50);
    hasfit_match = kffitdir.make<TH1F>("hasfit_match", "hasfit_match", 3, 0, 3);
    nchi2_match = kffitdir.make<TH1F>("nchi2_match", "nchi2_match", 50, 0, 100);
    nvalidpoints_match = kffitdir.make<TH1F>("nvalidpoints_match", "nvalidpoints_match", 100, 0, 1000);
    nrejectpoints_match = kffitdir.make<TH1F>("nrejectpoints_match", "nrejectpoints_match", 100, 0, 100);
    fracreject_match = kffitdir.make<TH1F>("fracreject_match", "fracreject_match", 100, 0, 1);
    //
    dXkf_okid = kffitdir.make<TH1F>("dXkf_okid","dXkf_okid",50,-10,10);
    dXtc_okid = kffitdir.make<TH1F>("dXtc_okid","dXtc_okid",50,-10,10);
    dYkf_okid = kffitdir.make<TH1F>("dYkf_okid","dYkf_okid",50,-10,10);
    dYtc_okid = kffitdir.make<TH1F>("dYtc_okid","dYtc_okid",50,-10,10);
    dZkf_okid = kffitdir.make<TH1F>("dZkf_okid","dZkf_okid",50,-10,10);
    dZtc_okid = kffitdir.make<TH1F>("dZtc_okid","dZtc_okid",50,-10,10);
    dUXkf_okid = kffitdir.make<TH1F>("dUXkf_okid","dUXkf_okid",50,-10,10);
    dUXtc_okid = kffitdir.make<TH1F>("dUXtc_okid","dUXtc_okid",50,-10,10);
    dUYkf_okid = kffitdir.make<TH1F>("dUYkf_okid","dUYkf_okid",50,-10,10);
    dUYtc_okid = kffitdir.make<TH1F>("dUYtc_okid","dUYtc_okid",50,-10,10);
    dUZkf_okid = kffitdir.make<TH1F>("dUZkf_okid","dUZkf_okid",50,-10,10);
    dUZtc_okid = kffitdir.make<TH1F>("dUZtc_okid","dUZtc_okid",50,-10,10);
    dXpull_okid = kffitdir.make<TH1F>("dXpull_okid","dXpull_okid",50,-10,10);
    dYpull_okid = kffitdir.make<TH1F>("dYpull_okid","dYpull_okid",50,-10,10);
    dZpull_okid = kffitdir.make<TH1F>("dZpull_okid","dZpull_okid",50,-10,10);
    dUXpull_okid = kffitdir.make<TH1F>("dUXpull_okid","dUXpull_okid",50,-10,10);
    dUYpull_okid = kffitdir.make<TH1F>("dUYpull_okid","dUYpull_okid",50,-10,10);
    dUZpull_okid = kffitdir.make<TH1F>("dUZpull_okid","dUZpull_okid",50,-10,10);
    hasfit_okid = kffitdir.make<TH1F>("hasfit_okid", "hasfit_okid", 3, 0, 3);
    nchi2_okid = kffitdir.make<TH1F>("nchi2_okid", "nchi2_okid", 50, 0, 100);
    //
    dXkf_wrongid = kffitdir.make<TH1F>("dXkf_wrongid","dXkf_wrongid",50,-10,10);
    dXtc_wrongid = kffitdir.make<TH1F>("dXtc_wrongid","dXtc_wrongid",50,-10,10);
    dYkf_wrongid = kffitdir.make<TH1F>("dYkf_wrongid","dYkf_wrongid",50,-10,10);
    dYtc_wrongid = kffitdir.make<TH1F>("dYtc_wrongid","dYtc_wrongid",50,-10,10);
    dZkf_wrongid = kffitdir.make<TH1F>("dZkf_wrongid","dZkf_wrongid",50,-10,10);
    dZtc_wrongid = kffitdir.make<TH1F>("dZtc_wrongid","dZtc_wrongid",50,-10,10);
    dUXkf_wrongid = kffitdir.make<TH1F>("dUXkf_wrongid","dUXkf_wrongid",50,-10,10);
    dUXtc_wrongid = kffitdir.make<TH1F>("dUXtc_wrongid","dUXtc_wrongid",50,-10,10);
    dUYkf_wrongid = kffitdir.make<TH1F>("dUYkf_wrongid","dUYkf_wrongid",50,-10,10);
    dUYtc_wrongid = kffitdir.make<TH1F>("dUYtc_wrongid","dUYtc_wrongid",50,-10,10);
    dUZkf_wrongid = kffitdir.make<TH1F>("dUZkf_wrongid","dUZkf_wrongid",50,-10,10);
    dUZtc_wrongid = kffitdir.make<TH1F>("dUZtc_wrongid","dUZtc_wrongid",50,-10,10);
    hasfit_wrongid = kffitdir.make<TH1F>("hasfit_wrongid", "hasfit_wrongid", 3, 0, 3);
    nchi2_wrongid = kffitdir.make<TH1F>("nchi2_wrongid", "nchi2_wrongid", 50, 0, 100);
    pdgid_wrongid = kffitdir.make<TH1F>("pdgid_wrongid", "pdgid_wrongid", 10, 0, 10);
    //
    covtrace_nomatch = kffitdir.make<TH1F>("covtrace_nomatch", "covtrace_nomatch", 100, 0, 50);
    hasfit_nomatch = kffitdir.make<TH1F>("hasfit_nomatch", "hasfit_nomatch", 3, 0, 3);
    nchi2_nomatch = kffitdir.make<TH1F>("nchi2_nomatch", "nchi2_nomatch", 50, 0, 100);
    nvalidpoints_nomatch = kffitdir.make<TH1F>("nvalidpoints_nomatch", "nvalidpoints_nomatch", 100, 0, 1000);
    nrejectpoints_nomatch = kffitdir.make<TH1F>("nrejectpoints_nomatch", "nrejectpoints_nomatch", 100, 0, 100);
    fracreject_nomatch = kffitdir.make<TH1F>("fracreject_nomatch", "fracreject_nomatch", 100, 0, 1);
  }
  
  
} // tca

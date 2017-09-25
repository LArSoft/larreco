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

  }
  
  
} // tca
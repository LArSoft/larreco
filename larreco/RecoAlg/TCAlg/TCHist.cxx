#include "larreco/RecoAlg/TCAlg/TCHist.h"

namespace  tca {
  
  void HistStuff::CreateHists(art::TFileService& tfs)
  {
    
    // study electrons
    fChgRMS[0] = tfs.make<TH1F>("ChgRMS0","ChgRMS - Electrons", 50, 0, 1);
    fChgRMS[1] = tfs.make<TH1F>("ChgRMS1","ChgRMS - Muons", 50, 0, 1);
    fChgRMS[2] = tfs.make<TH1F>("ChgRMS2","ChgRMS - Pions", 50, 0, 1);
    fChgRMS[3] = tfs.make<TH1F>("ChgRMS3","ChgRMS - Kaons", 50, 0, 1);
    fChgRMS[4] = tfs.make<TH1F>("ChgRMS4","ChgRMS - Protons", 50, 0, 1);
    
    fMomAsym[0] = tfs.make<TH1F>("MomAsym0","MomAsym - Electrons", 50, 0, 1);
    fMomAsym[1] = tfs.make<TH1F>("MomAsym1","MomAsym - Muons", 50, 0, 1);
    fMomAsym[2] = tfs.make<TH1F>("MomAsym2","MomAsym - Pions", 50, 0, 1);
    fMomAsym[3] = tfs.make<TH1F>("MomAsym3","MomAsym - Kaons", 50, 0, 1);
    fMomAsym[4] = tfs.make<TH1F>("MomAsym4","MomAsym - Protons", 50, 0, 1);
    
    fElectronLike[0] = tfs.make<TH1F>("ElectronLike0","ElectronLike - Electrons", 50, 0, 1);
    fElectronLike[1] = tfs.make<TH1F>("ElectronLike1","ElectronLike - Muons", 50, 0, 1);
    fElectronLike[2] = tfs.make<TH1F>("ElectronLike2","ElectronLike - Pions", 50, 0, 1);
    fElectronLike[3] = tfs.make<TH1F>("ElectronLike3","ElectronLike - Kaons", 50, 0, 1);
    fElectronLike[4] = tfs.make<TH1F>("ElectronLike4","ElectronLike - Protons", 50, 0, 1);
    
    fElectronLike_Len[0] = tfs.make<TH2F>("ElectronLike_Len0","ElectronLike vs length - Electrons", 100, 0, 100, 100, 0, 1);
    fElectronLike_Len[1] = tfs.make<TH2F>("ElectronLike_Len1","ElectronLike vs length - Muons", 100, 0, 100, 100, 0, 1);
    fElectronLike_Len[2] = tfs.make<TH2F>("ElectronLike_Len2","ElectronLike vs length - Pions", 100, 0, 100, 100, 0, 1);
    fElectronLike_Len[3] = tfs.make<TH2F>("ElectronLike_Len3","ElectronLike vs length - Kaons", 100, 0, 100, 100, 0, 1);
    fElectronLike_Len[4] = tfs.make<TH2F>("ElectronLike_Len4","ElectronLike vs length - Protons", 100, 0, 100, 100, 0, 1);
    
    // study showers
    AlongTrans1 = tfs.make<TH2F>("AlongTrans1","Trans vs Along 100 MeV",400, -50, 150, 40, 0, 20);
    AlongTrans5 = tfs.make<TH2F>("AlongTrans5","Trans vs Along 500 MeV",400, -50, 150, 40, 0, 20);
    AlongTrans9 = tfs.make<TH2F>("AlongTrans9","Trans vs Along 900 MeV",400, -50, 150, 40, 0, 20);
    
    fChgToMeV[0] = tfs.make<TH1F>("ChgToMeV0","ChgToMeV Pln 0", 80, 0, 0.04);
    fChgToMeV[1] = tfs.make<TH1F>("ChgToMeV1","ChgToMeV Pln 1", 80, 0, 0.04);
    fChgToMeV[2] = tfs.make<TH1F>("ChgToMeV2","ChgToMeV Pln 2", 80, 0, 0.04);
    fChgToMeV_Etru = tfs.make<TProfile>("ChgToMeV_Etru","ChgToMeV vs Etru Pln 0", 10, 0, 1000);
    
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
    
    fTruT[0] = tfs.make<TH1F>("TruT_El","True KE (MeV) - Electrons", 100, 0, 1000);
    fTruT[1] = tfs.make<TH1F>("TruT_Mu","True KE (MeV) - Muons", 100, 0, 1000);
    fTruT[2] = tfs.make<TH1F>("TruT_Pi","True KE (MeV) - Pions", 100, 0, 1000);
    fTruT[3] = tfs.make<TH1F>("TruT_Ka","True KE (MeV) - Kaons", 100, 0, 1000);
    fTruT[4] = tfs.make<TH1F>("TruT_Pr","True KE (MeV) - Protons", 100, 0, 1000);
    
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
    
    fNearTj[0] = tfs.make<TProfile>("fNearTj0","Frac Tps NearTj - Electrons", 5, 0, 500);
    fNearTj[1] = tfs.make<TProfile>("fNearTj1","Frac Tps NearTj - Muons", 5, 0, 500);
    fNearTj[2] = tfs.make<TProfile>("fNearTj2","Frac Tps NearTj - Pions", 5, 0, 500);
    fNearTj[3] = tfs.make<TProfile>("fNearTj3","Frac Tps NearTj - Kaons", 5, 0, 500);
    fNearTj[4] = tfs.make<TProfile>("fNearTj4","Frac Tps NearTj - Protons", 5, 0, 500);

    
    PDGCode_reco_true = tfs.make<TH2F>("PDGCode_reco_true", "PDGCode Reco vs True", 5, -0.5, 4.5, 5, -0.5, 4.5);
    // PFParticle start position reco - true
    fPFPStartdX[0] = tfs.make<TH1F>("PFPStartdX0","PFP Start dX Reco-MC - Electrons", 100, -10, 10);
    fPFPStartdX[1] = tfs.make<TH1F>("PFPStartdX1","PFP Start dX Reco-MC - Muons", 100, -10, 10);
    fPFPStartdX[2] = tfs.make<TH1F>("PFPStartdX2","PFP Start dX Reco-MC - Pions", 100, -10, 10);
    fPFPStartdX[3] = tfs.make<TH1F>("PFPStartdX3","PFP Start dX Reco-MC - Kaons", 100, -10, 10);
    fPFPStartdX[4] = tfs.make<TH1F>("PFPStartdX4","PFP Start dX Reco-MC - Protons", 100, -10, 10);
    fPFPStartdY[0] = tfs.make<TH1F>("PFPStartdY0","PFP Start dY Reco-MC - Electrons", 100, -10, 10);
    fPFPStartdY[1] = tfs.make<TH1F>("PFPStartdY1","PFP Start dY Reco-MC - Muons", 100, -10, 10);
    fPFPStartdY[2] = tfs.make<TH1F>("PFPStartdY2","PFP Start dY Reco-MC - Pions", 100, -10, 10);
    fPFPStartdY[3] = tfs.make<TH1F>("PFPStartdY3","PFP Start dY Reco-MC - Kaons", 100, -10, 10);
    fPFPStartdY[4] = tfs.make<TH1F>("PFPStartdY4","PFP Start dY Reco-MC - Protons", 100, -10, 10);
    fPFPStartdZ[0] = tfs.make<TH1F>("PFPStartdZ0","PFP Start dZ Reco-MC - Electrons", 100, -10, 10);
    fPFPStartdZ[1] = tfs.make<TH1F>("PFPStartdZ1","PFP Start dZ Reco-MC - Muons", 100, -10, 10);
    fPFPStartdZ[2] = tfs.make<TH1F>("PFPStartdZ2","PFP Start dZ Reco-MC - Pions", 100, -10, 10);
    fPFPStartdZ[3] = tfs.make<TH1F>("PFPStartdZ3","PFP Start dZ Reco-MC - Kaons", 100, -10, 10);
    fPFPStartdZ[4] = tfs.make<TH1F>("PFPStartdZ4","PFP Start dZ Reco-MC - Protons", 100, -10, 10);
    fPFPStartEnd = tfs.make<TH1F>("PFPStartEnd","PFP Start End", 2, -0.001, 1.001);
    // PFParticle start direction reco - true
    fPFPStartAngDiff[0] = tfs.make<TH1F>("PFPStartAngDiff0","PFP Start Ang Reco-MC - Electrons", 100, 0, M_PI);
    fPFPStartAngDiff[1] = tfs.make<TH1F>("PFPStartAngDiff1","PFP Start Ang Reco-MC - Muons", 100, 0, M_PI);
    fPFPStartAngDiff[2] = tfs.make<TH1F>("PFPStartAngDiff2","PFP Start Ang Reco-MC - Pions", 100, 0, M_PI);
    fPFPStartAngDiff[3] = tfs.make<TH1F>("PFPStartAngDiff3","PFP Start Ang Reco-MC - Kaons", 100, 0, M_PI);
    fPFPStartAngDiff[4] = tfs.make<TH1F>("PFPStartAngDiff4","PFP Start Ang Reco-MC - Protons", 100, 0, M_PI);
    
    fEff = tfs.make<TH1F>("Eff","Efficiency", 50, 0, 1);
    fPur = tfs.make<TH1F>("Pur","Purity", 50, 0, 1);
/* TTree used to develop the TMVA showerParentReader
    fShowerParentSig = tfs.make<TTree>("shwr_parent_tree_sig", "shwr_parent_tree_sig");
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
    
    fShowerParentBkg = tfs.make<TTree>("shwr_parent_tree_bkg", "shwr_parent_tree_bkg");
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

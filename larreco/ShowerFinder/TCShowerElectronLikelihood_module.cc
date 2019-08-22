// -------------------------------------------------
// makes shower profile templates
//
// Author: Rory Fitzpatrick (roryfitz@umich.edu)
// Created: 8/3/18
// -------------------------------------------------

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace shower {

  class TCShowerElectronLikelihood : public art::EDAnalyzer {

  public:

    explicit TCShowerElectronLikelihood(fhicl::ParameterSet const& pset);

    void beginJob();
    void analyze(const art::Event& evt);

    void getShowerProfile(std::vector< art::Ptr<recob::Hit> > showerhits, TVector3 shwvtx, TVector3 shwdir);
    void findEnergyBin();
    void getLongLikelihood();
    void getTranLikelihood();

  protected:

  private:
    void resetProfiles();

    std::string fTemplateFile;
    std::string fROOTfile;

    TH3F* longTemplate;
    TH3F* tranTemplate;
    TH3F* tranTemplate_1;
    TH3F* tranTemplate_2;
    TH3F* tranTemplate_3;
    TH3F* tranTemplate_4;
    TH3F* tranTemplate_5;
    TProfile2D* longTemplateProf2D;
    TProfile2D* tranTemplateProf2D;
    TProfile2D* tranTemplateProf2D_1;
    TProfile2D* tranTemplateProf2D_2;
    TProfile2D* tranTemplateProf2D_3;
    TProfile2D* tranTemplateProf2D_4;
    TProfile2D* tranTemplateProf2D_5;

    //TTree* fTree;
    TH1F* energyDist;
    TH1F* longLikelihoodHist;
    TH1F* tranLikelihoodHist;
    // just for printing purposes
    TH1F* longProfHist;
    TH1F* tranProfHist_1;
    TH1F* tranProfHist_2;
    TH1F* tranProfHist_3;
    TH1F* tranProfHist_4;
    TH1F* tranProfHist_5;

    TH1F* longProfile;
    TH1F* tranProfile;
    TH1F* tranProfile_1;
    TH1F* tranProfile_2;
    TH1F* tranProfile_3;
    TH1F* tranProfile_4;
    TH1F* tranProfile_5;

    int energyGuess;
    double energyChi2;
    int maxt;

    double longLikelihood;
    double tranLikelihood;
    double tranLikelihood_1;
    double tranLikelihood_2;
    double tranLikelihood_3;
    double tranLikelihood_4;
    double tranLikelihood_5;

    const int LBINS = 20;
    const int LMIN = 0;
    const int LMAX = 5;

    const int TBINS = 20;
    const int TMIN = -5;
    const int TMAX = 5;

    /*
    const int EBINS = 10;
    const double EMIN = 0.5;
    const double EMAX = 10.5;
    */

    // unused const int EBINS = 20;
    // unused const double EMIN = 0.5;
    // unused const double EMAX = 20.5;

    const double X0 = 14;

    std::string fHitModuleLabel;
    std::string fShowerModuleLabel;
    std::string fTemplateModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fDigitModuleLabel;

    calo::CalorimetryAlg fCalorimetryAlg;

  }; // class TCShowerElectronLikelihood

} // shower

// -------------------------------------------------

shower::TCShowerElectronLikelihood::TCShowerElectronLikelihood(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fHitModuleLabel           (pset.get< std::string >("HitModuleLabel", "trajcluster" ) ),
  fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel", "tcshower" ) ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel", "generator") ),
  fCalorimetryAlg           (pset.get< fhicl::ParameterSet >("CalorimetryAlg") ) {
  fTemplateFile           = pset.get< std::string >("TemplateFile");
  cet::search_path sp("FW_SEARCH_PATH");
  if( !sp.find_file(fTemplateFile, fROOTfile) )
    throw cet::exception("TCShowerElectronLikelihood") << "cannot find the root template file: \n"
						       << fTemplateFile
						       << "\n bail ungracefully.\n";

  TFile *file = TFile::Open(fROOTfile.c_str());

  longTemplate = (TH3F*)file->Get("tcshowertemplate/fLongitudinal");
  tranTemplate = (TH3F*)file->Get("tcshowertemplate/fTransverse");
  tranTemplate_1 = (TH3F*)file->Get("tcshowertemplate/fTransverse_1");
  tranTemplate_2 = (TH3F*)file->Get("tcshowertemplate/fTransverse_2");
  tranTemplate_3 = (TH3F*)file->Get("tcshowertemplate/fTransverse_3");
  tranTemplate_4 = (TH3F*)file->Get("tcshowertemplate/fTransverse_4");
  tranTemplate_5 = (TH3F*)file->Get("tcshowertemplate/fTransverse_5");

  longTemplateProf2D = (TProfile2D*)file->Get("tcshowertemplate/fShowerProfileRecoLong2D");
  tranTemplateProf2D = (TProfile2D*)file->Get("tcshowertemplate/fShowerProfileRecoTrans2D");
  tranTemplateProf2D_1 = (TProfile2D*)file->Get("tcshowertemplate/fShowerProfileRecoTrans2D_1");
  tranTemplateProf2D_2 = (TProfile2D*)file->Get("tcshowertemplate/fShowerProfileRecoTrans2D_2");
  tranTemplateProf2D_3 = (TProfile2D*)file->Get("tcshowertemplate/fShowerProfileRecoTrans2D_3");
  tranTemplateProf2D_4 = (TProfile2D*)file->Get("tcshowertemplate/fShowerProfileRecoTrans2D_4");
  tranTemplateProf2D_5 = (TProfile2D*)file->Get("tcshowertemplate/fShowerProfileRecoTrans2D_5");

  longProfile = new TH1F("longProfile", "longitudinal shower profile;t;Q", LBINS, LMIN, LMAX);
  tranProfile = new TH1F("tranProfile", "transverse shower profile;dist (cm);Q", TBINS, TMIN, TMAX);;

  tranProfile_1 = new TH1F("tranProfile_1", "transverse shower profile [0 <= t < 1];dist (cm);Q", TBINS, TMIN, TMAX);;
  tranProfile_2 = new TH1F("tranProfile_2", "transverse shower profile [1 <= t < 2];dist (cm);Q", TBINS, TMIN, TMAX);;
  tranProfile_3 = new TH1F("tranProfile_3", "transverse shower profile [2 <= t < 3];dist (cm);Q", TBINS, TMIN, TMAX);;
  tranProfile_4 = new TH1F("tranProfile_4", "transverse shower profile [3 <= t < 4];dist (cm);Q", TBINS, TMIN, TMAX);;
  tranProfile_5 = new TH1F("tranProfile_5", "transverse shower profile [4 <= t < 5];dist (cm);Q", TBINS, TMIN, TMAX);;
} // TCShowerElectronLikelihood

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::beginJob() {

  art::ServiceHandle<art::TFileService const> tfs;
  //fTree = tfs->make<TTree>("elikelihood", "elikelihood");

  energyDist = tfs->make<TH1F>("energyDist", "true energy - guess energy", 41, -20.5, 20.5);
  longLikelihoodHist = tfs->make<TH1F>("longLikelihoodHist", "longitudinal likelihood", 20, 0, 2);
  tranLikelihoodHist = tfs->make<TH1F>("tranLikelihoodHist", "transverse likelihood", 20, 0, 3);

  // just for printing purposes
  longProfHist = tfs->make<TH1F>("longProfHist", "longitudinal e- profile (reco);t;Q", LBINS, LMIN, LMAX);
  tranProfHist_1 = tfs->make<TH1F>("tranProfHist", "transverse e- profile (reco) [0 <= t < 1];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_2 = tfs->make<TH1F>("tranProfHist", "transverse e- profile (reco) [1 <= t < 2];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_3 = tfs->make<TH1F>("tranProfHist", "transverse e- profile (reco) [2 <= t < 3];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_4 = tfs->make<TH1F>("tranProfHist", "transverse e- profile (reco) [3 <= t < 4];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_5 = tfs->make<TH1F>("tranProfHist", "transverse e- profile (reco) [4 <= t < 5];dist (cm);Q", TBINS, TMIN, TMAX);

} // beginJob

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::analyze(const art::Event& evt) {

  resetProfiles();

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if (evt.getByLabel(fShowerModuleLabel,showerListHandle))
    art::fill_ptr_vector(showerlist, showerListHandle);

  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  art::FindManyP<recob::Hit> shwfm(showerListHandle, evt, fShowerModuleLabel);

  if (showerlist.size()) {
    std::vector< art::Ptr<recob::Hit> > showerhits = shwfm.at(0);
    getShowerProfile(showerhits, showerlist[0]->ShowerStart(), showerlist[0]->Direction());

    maxt = std::ceil((90 - showerlist[0]->ShowerStart().Z())/X0);

    findEnergyBin();
    getLongLikelihood();
    getTranLikelihood();

    longLikelihoodHist->Fill(longLikelihood);
    tranLikelihoodHist->Fill(tranLikelihood);

    // check true shower energy
    if (mclist.size()) {
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
      if (mctruth->NeutrinoSet()) {
	if (std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 12 && mctruth->GetNeutrino().CCNC() == 0) {
	  double elep =  mctruth->GetNeutrino().Lepton().E();
	  //	  std::cout << "true shower energy: " << elep << std::endl;
	  //	  std::cout << "energy guess:       " << bestE << " (" << bestchi2 << ")" << std::endl;
	  energyDist->Fill(elep - energyGuess);
	} // cc nue
      } // neutrinoset
    } // mclist
  } // showerlist

  longProfHist = longProfile;
  tranProfHist_1 = tranProfile_1;
  tranProfHist_2 = tranProfile_2;
  tranProfHist_3 = tranProfile_3;
  tranProfHist_4 = tranProfile_4;
  tranProfHist_5 = tranProfile_5;

  //fTree->Fill();

} // analyze

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::resetProfiles() {

  longProfile->Reset();
  tranProfile->Reset();
  tranProfile_1->Reset();
  tranProfile_2->Reset();
  tranProfile_3->Reset();
  tranProfile_4->Reset();
  tranProfile_5->Reset();

  energyGuess = -9999;
  energyChi2 = -9999;
  maxt = -9999;

  longLikelihood = -9999;
  tranLikelihood = -9999;
  tranLikelihood_1 = -9999;
  tranLikelihood_2 = -9999;
  tranLikelihood_3 = -9999;
  tranLikelihood_4 = -9999;
  tranLikelihood_5 = -9999;

  return;

} // resetProfiles

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::getShowerProfile(std::vector< art::Ptr<recob::Hit> > showerhits, TVector3 shwvtx, TVector3 shwdir) {

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry const> geom;

  auto collectionPlane = geo::PlaneID(0, 0, 1);

  double shwVtxTime = detprop->ConvertXToTicks(shwvtx[0], collectionPlane);
  double shwVtxWire = geom->WireCoordinate(shwvtx[1], shwvtx[2], collectionPlane);

  double shwTwoTime = detprop->ConvertXToTicks(shwvtx[0]+shwdir[0], collectionPlane);
  double shwTwoWire = geom->WireCoordinate(shwvtx[1]+shwdir[1], shwvtx[2]+shwdir[2], collectionPlane);

  for (size_t i = 0; i < showerhits.size(); ++i) {
    if (showerhits[i]->WireID().Plane != collectionPlane.Plane) continue;

    double wirePitch = geom->WirePitch(showerhits[i]->WireID());
    double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns

    double xvtx = shwVtxTime * tickToDist;
    double yvtx = shwVtxWire * wirePitch;

    double xtwo = shwTwoTime * tickToDist;
    double ytwo = shwTwoWire * wirePitch;

    double xtwoorth = (ytwo - yvtx) + xvtx;
    double ytwoorth = -(xtwo - xvtx) + yvtx;

    double xhit = showerhits[i]->PeakTime() * tickToDist;
    double yhit = showerhits[i]->WireID().Wire * wirePitch;

    double ldist = std::abs((ytwoorth-yvtx)*xhit - (xtwoorth-xvtx)*yhit + xtwoorth*yvtx - ytwoorth*xvtx)/std::sqrt( pow((ytwoorth-yvtx), 2) + pow((xtwoorth-xvtx), 2) );
    double tdist = ((ytwo-yvtx)*xhit - (xtwo-xvtx)*yhit + xtwo*yvtx - ytwo*xvtx)/std::sqrt( pow((ytwo-yvtx), 2) + pow((xtwo-xvtx), 2) );

    double to3D = 1. / sqrt( pow(xvtx-xtwo,2) + pow(yvtx-ytwo,2) ) ; // distance between two points in 3D space is one
    ldist *= to3D;
    tdist *= to3D;
    double t = ldist/X0;

    double Q = showerhits[i]->Integral() * fCalorimetryAlg.LifetimeCorrection(showerhits[i]->PeakTime());

    longProfile->Fill(t, Q);
    tranProfile->Fill(tdist, Q);

    if (t < 1) tranProfile_1->Fill(tdist, Q);
    else if (t < 2) tranProfile_2->Fill(tdist, Q);
    else if (t < 3) tranProfile_3->Fill(tdist, Q);
    else if (t < 4) tranProfile_4->Fill(tdist, Q);
    else if (t < 5) tranProfile_5->Fill(tdist, Q);

  } // loop through showerhits

  return;

} // getShowerProfile

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::findEnergyBin() {

  if (longProfile->GetNbinsX() != longTemplate->GetNbinsX())
    throw cet::exception("TCShowerElectronLikelihood") << "Bin mismatch in longitudinal profile template \n";

  if (tranProfile->GetNbinsX() != tranTemplate->GetNbinsX())
    throw cet::exception("TCShowerElectronLikelihood") << "Bin mismatch in transverse profile template \n";

  double chi2min = 999999;
  double bestbin = -1;

  int ebins = longTemplate->GetNbinsY();
  int lbins = longTemplate->GetNbinsX();
  int tbins = tranTemplate->GetNbinsX();

  //  lbins = floor(lbins/2); // only use the first half of the bins

  TProfile* ltemp;
  TProfile* ttemp_1;
  TProfile* ttemp_2;
  TProfile* ttemp_3;
  TProfile* ttemp_4;
  TProfile* ttemp_5;

  for (int i = 0; i < ebins; ++i) {
    double thischi2 = 0;

    ltemp = (TProfile*)longTemplateProf2D->ProfileX("_x", i+1, i+1);
    ttemp_1 = (TProfile*)tranTemplateProf2D_1->ProfileX("_x_1", i+1, i+1);
    ttemp_2 = (TProfile*)tranTemplateProf2D_2->ProfileX("_x_2", i+1, i+1);
    ttemp_3 = (TProfile*)tranTemplateProf2D_3->ProfileX("_x_3", i+1, i+1);
    ttemp_4 = (TProfile*)tranTemplateProf2D_4->ProfileX("_x_4", i+1, i+1);
    ttemp_5 = (TProfile*)tranTemplateProf2D_5->ProfileX("_x_5", i+1, i+1);

    int nlbins = 0;
    int ntbins = 0;

    for (int j = 0; j < lbins; ++j) {
      double obs = longProfile->GetBinContent(j+1);
      double exp = ltemp->GetBinContent(j+1);
      if (obs != 0) {
	thischi2 += pow(obs - exp, 2) / exp;
	++nlbins;
      }
    } // loop through longitudinal bins

    for (int j = 0; j < tbins; ++j) {
      double obs = tranProfile_1->GetBinContent(j+1);
      double exp = ttemp_1->GetBinContent(j+1);
      if (obs != 0) {
	thischi2 += pow(obs - exp, 2) / exp;
	++ntbins;
      }

      obs = tranProfile_2->GetBinContent(j+1);
      exp = ttemp_2->GetBinContent(j+1);
      if (obs != 0) {
	thischi2 += pow(obs - exp, 2) / exp;
	++ntbins;
      }

      obs = tranProfile_3->GetBinContent(j+1);
      exp = ttemp_3->GetBinContent(j+1);
      if (obs != 0) {
	thischi2 += pow(obs - exp, 2) / exp;
	++ntbins;
      }

      obs = tranProfile_4->GetBinContent(j+1);
      exp = ttemp_4->GetBinContent(j+1);
      if (obs != 0) {
	thischi2 += pow(obs - exp, 2) / exp;
	++ntbins;
      }

      obs = tranProfile_5->GetBinContent(j+1);
      exp = ttemp_5->GetBinContent(j+1);
      if (obs != 0) {
	thischi2 += pow(obs - exp, 2) / exp;
	++ntbins;
      }
    } // loop through longitudinal bins

    thischi2 /= (nlbins+ntbins);

    if (thischi2 < chi2min) {
      chi2min = thischi2;
      bestbin = i;
    }

  } // loop through energy bins

  energyChi2 = chi2min;
  energyGuess = bestbin+1;

  return;

} // findEnergyBin

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::getLongLikelihood() {

  if (energyGuess < 0) return;
  int energyBin = energyGuess;

  longLikelihood = 0;
  int nbins = 0;

  for (int i = 0; i < LBINS; ++i) {
    double qval = longProfile->GetBinContent(i+1);
    int qbin = longTemplate->GetZaxis()->FindBin(qval);
    int binentries = longTemplate->GetBinContent(i+1, energyBin, qbin);
    int totentries = longTemplate->Integral(i+1, i+1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries/totentries * 100;
      if (binentries > 0) longLikelihood += log(prob);
    }
  } // loop through

  longLikelihood /= nbins;

  std::cout << longLikelihood << std::endl;

  return;

} // getLongLikelihood

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::getTranLikelihood() {

  if (energyGuess < 0) return;
  int energyBin = energyGuess;

  tranLikelihood_1 = 0;
  tranLikelihood_2 = 0;
  tranLikelihood_3 = 0;
  tranLikelihood_4 = 0;
  tranLikelihood_5 = 0;

  double qval;
  int qbin, binentries, totentries;

  int nbins = 0;

  for (int i = 0; i < TBINS; ++i) {
    qval = tranProfile_1->GetBinContent(i+1);
    qbin = tranTemplate_1->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_1->GetBinContent(i+1, energyBin, qbin);
    totentries = tranTemplate_1->Integral(i+1, i+1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries/totentries * 100;
      if (binentries > 0) tranLikelihood_1 += log(prob);
    }

    qval = tranProfile_2->GetBinContent(i+1);
    qbin = tranTemplate_2->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_2->GetBinContent(i+1, energyBin, qbin);
    totentries = tranTemplate_2->Integral(i+1, i+1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries/totentries * 100;
      if (binentries > 0) tranLikelihood_2 += log(prob);
    }

    qval = tranProfile_3->GetBinContent(i+1);
    qbin = tranTemplate_3->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_3->GetBinContent(i+1, energyBin, qbin);
    totentries = tranTemplate_3->Integral(i+1, i+1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries/totentries * 100;
      if (binentries > 0) tranLikelihood_3 += log(prob);
    }

    qval = tranProfile_4->GetBinContent(i+1);
    qbin = tranTemplate_4->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_4->GetBinContent(i+1, energyBin, qbin);
    totentries = tranTemplate_4->Integral(i+1, i+1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries/totentries * 100;
      if (binentries > 0) tranLikelihood_4 += log(prob);
    }

    qval = tranProfile_5->GetBinContent(i+1);
    qbin = tranTemplate_5->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_5->GetBinContent(i+1, energyBin, qbin);
    totentries = tranTemplate_5->Integral(i+1, i+1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries/totentries * 100;
      if (binentries > 0) tranLikelihood_5 += log(prob);
    }

  } // loop through

  /*
  std::cout << tranLikelihood_1 << std::endl;
  std::cout << tranLikelihood_2 << std::endl;
  std::cout << tranLikelihood_3 << std::endl;
  std::cout << tranLikelihood_4 << std::endl;
  std::cout << tranLikelihood_5 << std::endl;
  */

  tranLikelihood = tranLikelihood_1 + tranLikelihood_2 + tranLikelihood_3 + tranLikelihood_4 + tranLikelihood_5;

  tranLikelihood /= nbins;

  std::cout << tranLikelihood << std::endl;

  return;

} // getTranLikelihood

// -------------------------------------------------

DEFINE_ART_MODULE(shower::TCShowerElectronLikelihood)

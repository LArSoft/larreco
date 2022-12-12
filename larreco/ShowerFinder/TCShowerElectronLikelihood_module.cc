// -------------------------------------------------
// makes shower profile templates
//
// Author: Rory Fitzpatrick (roryfitz@umich.edu)
// Created: 8/3/18
// -------------------------------------------------

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "cetlib/pow.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace shower {

  class TCShowerElectronLikelihood : public art::EDAnalyzer {
  public:
    explicit TCShowerElectronLikelihood(fhicl::ParameterSet const& pset);

  private:
    void beginJob() override;
    void analyze(const art::Event& evt) override;

    void getShowerProfile(detinfo::DetectorClocksData const& clockdata,
                          detinfo::DetectorPropertiesData const& detProp,
                          std::vector<art::Ptr<recob::Hit>> showerhits,
                          TVector3 shwvtx,
                          TVector3 shwdir);
    void findEnergyBin();
    double getLongLikelihood();
    double getTranLikelihood();

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

    const int LBINS = 20;
    const int LMIN = 0;
    const int LMAX = 5;

    const int TBINS = 20;
    const int TMIN = -5;
    const int TMAX = 5;

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

shower::TCShowerElectronLikelihood::TCShowerElectronLikelihood(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
  , fHitModuleLabel(pset.get<std::string>("HitModuleLabel", "trajcluster"))
  , fShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel", "tcshower"))
  , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel", "generator"))
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  fTemplateFile = pset.get<std::string>("TemplateFile");
  cet::search_path sp("FW_SEARCH_PATH");
  if (!sp.find_file(fTemplateFile, fROOTfile))
    throw cet::exception("TCShowerElectronLikelihood")
      << "cannot find the root template file: \n"
      << fTemplateFile << "\n bail ungracefully.\n";

  std::unique_ptr<TFile> file{TFile::Open(fROOTfile.c_str())};

  longTemplate = file->Get<TH3F>("tcshowertemplate/fLongitudinal");
  tranTemplate = file->Get<TH3F>("tcshowertemplate/fTransverse");
  tranTemplate_1 = file->Get<TH3F>("tcshowertemplate/fTransverse_1");
  tranTemplate_2 = file->Get<TH3F>("tcshowertemplate/fTransverse_2");
  tranTemplate_3 = file->Get<TH3F>("tcshowertemplate/fTransverse_3");
  tranTemplate_4 = file->Get<TH3F>("tcshowertemplate/fTransverse_4");
  tranTemplate_5 = file->Get<TH3F>("tcshowertemplate/fTransverse_5");

  longTemplateProf2D = file->Get<TProfile2D>("tcshowertemplate/fShowerProfileRecoLong2D");
  tranTemplateProf2D = file->Get<TProfile2D>("tcshowertemplate/fShowerProfileRecoTrans2D");
  tranTemplateProf2D_1 = file->Get<TProfile2D>("tcshowertemplate/fShowerProfileRecoTrans2D_1");
  tranTemplateProf2D_2 = file->Get<TProfile2D>("tcshowertemplate/fShowerProfileRecoTrans2D_2");
  tranTemplateProf2D_3 = file->Get<TProfile2D>("tcshowertemplate/fShowerProfileRecoTrans2D_3");
  tranTemplateProf2D_4 = file->Get<TProfile2D>("tcshowertemplate/fShowerProfileRecoTrans2D_4");
  tranTemplateProf2D_5 = file->Get<TProfile2D>("tcshowertemplate/fShowerProfileRecoTrans2D_5");

  longProfile = new TH1F("longProfile", "longitudinal shower profile;t;Q", LBINS, LMIN, LMAX);
  tranProfile = new TH1F("tranProfile", "transverse shower profile;dist (cm);Q", TBINS, TMIN, TMAX);

  tranProfile_1 = new TH1F(
    "tranProfile_1", "transverse shower profile [0 <= t < 1];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfile_2 = new TH1F(
    "tranProfile_2", "transverse shower profile [1 <= t < 2];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfile_3 = new TH1F(
    "tranProfile_3", "transverse shower profile [2 <= t < 3];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfile_4 = new TH1F(
    "tranProfile_4", "transverse shower profile [3 <= t < 4];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfile_5 = new TH1F(
    "tranProfile_5", "transverse shower profile [4 <= t < 5];dist (cm);Q", TBINS, TMIN, TMAX);
} // TCShowerElectronLikelihood

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::beginJob()
{

  art::ServiceHandle<art::TFileService const> tfs;

  energyDist = tfs->make<TH1F>("energyDist", "true energy - guess energy", 41, -20.5, 20.5);
  longLikelihoodHist = tfs->make<TH1F>("longLikelihoodHist", "longitudinal likelihood", 20, 0, 2);
  tranLikelihoodHist = tfs->make<TH1F>("tranLikelihoodHist", "transverse likelihood", 20, 0, 3);

  // just for printing purposes
  longProfHist =
    tfs->make<TH1F>("longProfHist", "longitudinal e- profile (reco);t;Q", LBINS, LMIN, LMAX);
  tranProfHist_1 = tfs->make<TH1F>(
    "tranProfHist", "transverse e- profile (reco) [0 <= t < 1];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_2 = tfs->make<TH1F>(
    "tranProfHist", "transverse e- profile (reco) [1 <= t < 2];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_3 = tfs->make<TH1F>(
    "tranProfHist", "transverse e- profile (reco) [2 <= t < 3];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_4 = tfs->make<TH1F>(
    "tranProfHist", "transverse e- profile (reco) [3 <= t < 4];dist (cm);Q", TBINS, TMIN, TMAX);
  tranProfHist_5 = tfs->make<TH1F>(
    "tranProfHist", "transverse e- profile (reco) [4 <= t < 5];dist (cm);Q", TBINS, TMIN, TMAX);

} // beginJob

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::analyze(const art::Event& evt)
{
  resetProfiles();

  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (evt.getByLabel(fHitModuleLabel, hitListHandle)) art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle<std::vector<recob::Shower>> showerListHandle;
  std::vector<art::Ptr<recob::Shower>> showerlist;
  if (evt.getByLabel(fShowerModuleLabel, showerListHandle))
    art::fill_ptr_vector(showerlist, showerListHandle);

  art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth>> mclist;
  if (evt.getByLabel(fGenieGenModuleLabel, mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  art::FindManyP<recob::Hit> shwfm(showerListHandle, evt, fShowerModuleLabel);

  if (showerlist.size()) {
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);
    std::vector<art::Ptr<recob::Hit>> showerhits = shwfm.at(0);
    getShowerProfile(
      clock_data, det_prop, showerhits, showerlist[0]->ShowerStart(), showerlist[0]->Direction());

    maxt = std::ceil((90 - showerlist[0]->ShowerStart().Z()) / X0);

    findEnergyBin();

    longLikelihoodHist->Fill(getLongLikelihood());
    tranLikelihoodHist->Fill(getTranLikelihood());

    // check true shower energy
    if (mclist.size()) {
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
      if (mctruth->NeutrinoSet()) {
        if (std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 12 &&
            mctruth->GetNeutrino().CCNC() == 0) {
          double elep = mctruth->GetNeutrino().Lepton().E();
          energyDist->Fill(elep - energyGuess);
        } // cc nue
      }   // neutrinoset
    }     // mclist
  }       // showerlist

  longProfHist = longProfile;
  tranProfHist_1 = tranProfile_1;
  tranProfHist_2 = tranProfile_2;
  tranProfHist_3 = tranProfile_3;
  tranProfHist_4 = tranProfile_4;
  tranProfHist_5 = tranProfile_5;
} // analyze

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::resetProfiles()
{

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
} // resetProfiles

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::getShowerProfile(
  detinfo::DetectorClocksData const& clockdata,
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::Hit>> showerhits,
  TVector3 shwvtx,
  TVector3 shwdir)
{
  art::ServiceHandle<geo::Geometry const> geom;

  auto collectionPlaneID = geo::PlaneID(0, 0, 1);
  auto const& collectionPlane = geom->Plane(collectionPlaneID);

  using namespace geo::vect;
  auto const shwvtx_p = toPoint(shwvtx);
  auto const shwvtx_p2 = shwvtx_p + toVector(shwdir);

  double shwVtxTime = detProp.ConvertXToTicks(shwvtx_p.X(), collectionPlaneID);
  double shwVtxWire = collectionPlane.WireCoordinate(shwvtx_p);

  double shwTwoTime = detProp.ConvertXToTicks(shwvtx_p2.X(), collectionPlaneID);
  double shwTwoWire = collectionPlane.WireCoordinate(shwvtx_p2);

  for (size_t i = 0; i < showerhits.size(); ++i) {
    if (showerhits[i]->WireID().Plane != collectionPlaneID.Plane) continue;

    double wirePitch = geom->Plane(showerhits[i]->WireID()).WirePitch();
    double tickToDist = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
    tickToDist *= 1.e-3 * sampling_rate(clockdata); // 1e-3 is conversion of 1/us to 1/ns

    double xvtx = shwVtxTime * tickToDist;
    double yvtx = shwVtxWire * wirePitch;

    double xtwo = shwTwoTime * tickToDist;
    double ytwo = shwTwoWire * wirePitch;

    double xtwoorth = (ytwo - yvtx) + xvtx;
    double ytwoorth = -(xtwo - xvtx) + yvtx;

    double xhit = showerhits[i]->PeakTime() * tickToDist;
    double yhit = showerhits[i]->WireID().Wire * wirePitch;

    double ldist = std::abs((ytwoorth - yvtx) * xhit - (xtwoorth - xvtx) * yhit + xtwoorth * yvtx -
                            ytwoorth * xvtx) /
                   std::hypot(ytwoorth - yvtx, xtwoorth - xvtx);
    double tdist = ((ytwo - yvtx) * xhit - (xtwo - xvtx) * yhit + xtwo * yvtx - ytwo * xvtx) /
                   std::hypot(ytwo - yvtx, xtwo - xvtx);

    double to3D = 1. / std::hypot(xvtx - xtwo,
                                  yvtx - ytwo); // distance between two points in 3D space is one
    ldist *= to3D;
    tdist *= to3D;
    double t = ldist / X0;

    double Q = showerhits[i]->Integral() *
               fCalorimetryAlg.LifetimeCorrection(clockdata, detProp, showerhits[i]->PeakTime());

    longProfile->Fill(t, Q);
    tranProfile->Fill(tdist, Q);

    if (t < 1)
      tranProfile_1->Fill(tdist, Q);
    else if (t < 2)
      tranProfile_2->Fill(tdist, Q);
    else if (t < 3)
      tranProfile_3->Fill(tdist, Q);
    else if (t < 4)
      tranProfile_4->Fill(tdist, Q);
    else if (t < 5)
      tranProfile_5->Fill(tdist, Q);

  } // loop through showerhits

  return;

} // getShowerProfile

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::findEnergyBin()
{

  if (longProfile->GetNbinsX() != longTemplate->GetNbinsX())
    throw cet::exception("TCShowerElectronLikelihood")
      << "Bin mismatch in longitudinal profile template \n";

  if (tranProfile->GetNbinsX() != tranTemplate->GetNbinsX())
    throw cet::exception("TCShowerElectronLikelihood")
      << "Bin mismatch in transverse profile template \n";

  double chi2min = 999999;
  double bestbin = -1;

  int ebins = longTemplate->GetNbinsY();
  int lbins = longTemplate->GetNbinsX();
  int tbins = tranTemplate->GetNbinsX();

  for (int i = 0; i < ebins; ++i) {
    double thischi2 = 0;

    TProfile* ltemp = longTemplateProf2D->ProfileX("_x", i + 1, i + 1);
    TProfile* ttemp_1 = tranTemplateProf2D_1->ProfileX("_x_1", i + 1, i + 1);
    TProfile* ttemp_2 = tranTemplateProf2D_2->ProfileX("_x_2", i + 1, i + 1);
    TProfile* ttemp_3 = tranTemplateProf2D_3->ProfileX("_x_3", i + 1, i + 1);
    TProfile* ttemp_4 = tranTemplateProf2D_4->ProfileX("_x_4", i + 1, i + 1);
    TProfile* ttemp_5 = tranTemplateProf2D_5->ProfileX("_x_5", i + 1, i + 1);

    int nlbins = 0;
    int ntbins = 0;

    for (int j = 0; j < lbins; ++j) {
      double obs = longProfile->GetBinContent(j + 1);
      double exp = ltemp->GetBinContent(j + 1);
      if (obs != 0) {
        thischi2 += cet::square(obs - exp) / exp;
        ++nlbins;
      }
    } // loop through longitudinal bins

    for (int j = 0; j < tbins; ++j) {
      double obs = tranProfile_1->GetBinContent(j + 1);
      double exp = ttemp_1->GetBinContent(j + 1);
      if (obs != 0) {
        thischi2 += cet::square(obs - exp) / exp;
        ++ntbins;
      }

      obs = tranProfile_2->GetBinContent(j + 1);
      exp = ttemp_2->GetBinContent(j + 1);
      if (obs != 0) {
        thischi2 += cet::square(obs - exp) / exp;
        ++ntbins;
      }

      obs = tranProfile_3->GetBinContent(j + 1);
      exp = ttemp_3->GetBinContent(j + 1);
      if (obs != 0) {
        thischi2 += cet::square(obs - exp) / exp;
        ++ntbins;
      }

      obs = tranProfile_4->GetBinContent(j + 1);
      exp = ttemp_4->GetBinContent(j + 1);
      if (obs != 0) {
        thischi2 += cet::square(obs - exp) / exp;
        ++ntbins;
      }

      obs = tranProfile_5->GetBinContent(j + 1);
      exp = ttemp_5->GetBinContent(j + 1);
      if (obs != 0) {
        thischi2 += cet::square(obs - exp) / exp;
        ++ntbins;
      }
    } // loop through longitudinal bins

    thischi2 /= (nlbins + ntbins);

    if (thischi2 < chi2min) {
      chi2min = thischi2;
      bestbin = i;
    }

  } // loop through energy bins

  energyChi2 = chi2min;
  energyGuess = bestbin + 1;
} // findEnergyBin

// -------------------------------------------------

double shower::TCShowerElectronLikelihood::getLongLikelihood()
{
  if (energyGuess < 0) return -9999.;
  int energyBin = energyGuess;

  double longLikelihood = 0.;
  int nbins = 0;

  for (int i = 0; i < LBINS; ++i) {
    double qval = longProfile->GetBinContent(i + 1);
    int qbin = longTemplate->GetZaxis()->FindBin(qval);
    int binentries = longTemplate->GetBinContent(i + 1, energyBin, qbin);
    int totentries = longTemplate->Integral(i + 1, i + 1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries / totentries * 100;
      if (binentries > 0) longLikelihood += log(prob);
    }
  } // loop through

  longLikelihood /= nbins;

  std::cout << longLikelihood << std::endl;
  return longLikelihood;
} // getLongLikelihood

// -------------------------------------------------

double shower::TCShowerElectronLikelihood::getTranLikelihood()
{
  if (energyGuess < 0) return -9999.;
  int energyBin = energyGuess;

  double tranLikelihood_1 = 0;
  double tranLikelihood_2 = 0;
  double tranLikelihood_3 = 0;
  double tranLikelihood_4 = 0;
  double tranLikelihood_5 = 0;

  double qval;
  int qbin, binentries, totentries;

  int nbins = 0;

  for (int i = 0; i < TBINS; ++i) {
    qval = tranProfile_1->GetBinContent(i + 1);
    qbin = tranTemplate_1->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_1->GetBinContent(i + 1, energyBin, qbin);
    totentries = tranTemplate_1->Integral(i + 1, i + 1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries / totentries * 100;
      if (binentries > 0) tranLikelihood_1 += log(prob);
    }

    qval = tranProfile_2->GetBinContent(i + 1);
    qbin = tranTemplate_2->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_2->GetBinContent(i + 1, energyBin, qbin);
    totentries = tranTemplate_2->Integral(i + 1, i + 1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries / totentries * 100;
      if (binentries > 0) tranLikelihood_2 += log(prob);
    }

    qval = tranProfile_3->GetBinContent(i + 1);
    qbin = tranTemplate_3->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_3->GetBinContent(i + 1, energyBin, qbin);
    totentries = tranTemplate_3->Integral(i + 1, i + 1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries / totentries * 100;
      if (binentries > 0) tranLikelihood_3 += log(prob);
    }

    qval = tranProfile_4->GetBinContent(i + 1);
    qbin = tranTemplate_4->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_4->GetBinContent(i + 1, energyBin, qbin);
    totentries = tranTemplate_4->Integral(i + 1, i + 1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries / totentries * 100;
      if (binentries > 0) tranLikelihood_4 += log(prob);
    }

    qval = tranProfile_5->GetBinContent(i + 1);
    qbin = tranTemplate_5->GetZaxis()->FindBin(qval);
    binentries = tranTemplate_5->GetBinContent(i + 1, energyBin, qbin);
    totentries = tranTemplate_5->Integral(i + 1, i + 1, energyBin, energyBin, 0, 100);
    if (qval > 0) {
      ++nbins;
      double prob = (double)binentries / totentries * 100;
      if (binentries > 0) tranLikelihood_5 += log(prob);
    }

  } // loop through

  double const tranLikelihood =
    (tranLikelihood_1 + tranLikelihood_2 + tranLikelihood_3 + tranLikelihood_4 + tranLikelihood_5) /
    nbins;
  std::cout << tranLikelihood << std::endl;
  return tranLikelihood;
} // getTranLikelihood

// -------------------------------------------------

DEFINE_ART_MODULE(shower::TCShowerElectronLikelihood)

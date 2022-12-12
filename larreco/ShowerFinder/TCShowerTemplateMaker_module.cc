// -------------------------------------------------
// makes shower profile templates
//
// Author: Rory Fitzpatrick (roryfitz@umich.edu)
// Created: 7/16/18
// -------------------------------------------------

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/MCBase/MCDataHolder.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace {
  constexpr geo::PlaneID collectionPlaneID{0, 0, 1};
}

namespace shower {

  class TCShowerTemplateMaker : public art::EDAnalyzer {
  public:
    explicit TCShowerTemplateMaker(fhicl::ParameterSet const& pset);

  private:
    void beginJob() override;
    void analyze(const art::Event& evt) override;

    void showerProfile(detinfo::DetectorClocksData const& clockData,
                       detinfo::DetectorPropertiesData const& detProp,
                       std::vector<art::Ptr<recob::Hit>> showerhits,
                       TVector3 shwvtx,
                       TVector3 shwdir,
                       double elep);
    void showerProfileTrue(detinfo::DetectorClocksData const& clockData,
                           detinfo::DetectorPropertiesData const& detProp,
                           std::vector<art::Ptr<recob::Hit>> allhits,
                           double elep);
    void showerProfileTrue(std::vector<art::Ptr<sim::SimChannel>> allchan,
                           simb::MCParticle electron);

    TProfile* fShowerProfileSimLong;
    TProfile* fShowerProfileHitLong;
    TProfile* fShowerProfileRecoLong;

    TProfile2D* fShowerProfileSimLong2D;
    TProfile2D* fShowerProfileHitLong2D;
    TProfile2D* fShowerProfileRecoLong2D;

    TProfile* fShowerProfileSimTrans;
    TProfile* fShowerProfileHitTrans;
    TProfile* fShowerProfileRecoTrans;

    TProfile2D* fShowerProfileSimTrans2D;
    TProfile2D* fShowerProfileHitTrans2D;
    TProfile2D* fShowerProfileRecoTrans2D;

    TProfile2D* fShowerProfileSimTrans2D_1;
    TProfile2D* fShowerProfileHitTrans2D_1;
    TProfile2D* fShowerProfileRecoTrans2D_1;

    TProfile2D* fShowerProfileSimTrans2D_2;
    TProfile2D* fShowerProfileHitTrans2D_2;
    TProfile2D* fShowerProfileRecoTrans2D_2;

    TProfile2D* fShowerProfileSimTrans2D_3;
    TProfile2D* fShowerProfileHitTrans2D_3;
    TProfile2D* fShowerProfileRecoTrans2D_3;

    TProfile2D* fShowerProfileSimTrans2D_4;
    TProfile2D* fShowerProfileHitTrans2D_4;
    TProfile2D* fShowerProfileRecoTrans2D_4;

    TProfile2D* fShowerProfileSimTrans2D_5;
    TProfile2D* fShowerProfileHitTrans2D_5;
    TProfile2D* fShowerProfileRecoTrans2D_5;

    TH3F* fLongitudinal;
    TH3F* fTransverse;
    TH3F* fTransverse_1;
    TH3F* fTransverse_2;
    TH3F* fTransverse_3;
    TH3F* fTransverse_4;
    TH3F* fTransverse_5;

    // electrons
    TProfile* fLongitudinal_electron;
    TProfile* fTransverse1_electron;
    TProfile* fTransverse2_electron;
    TProfile* fTransverse3_electron;
    TProfile* fTransverse4_electron;
    TProfile* fTransverse5_electron;
    // photons
    TProfile* fLongitudinal_photon;
    TProfile* fTransverse1_photon;
    TProfile* fTransverse2_photon;
    TProfile* fTransverse3_photon;
    TProfile* fTransverse4_photon;
    TProfile* fTransverse5_photon;
    // other
    TProfile* fLongitudinal_other;
    TProfile* fTransverse1_other;
    TProfile* fTransverse2_other;
    TProfile* fTransverse3_other;
    TProfile* fTransverse4_other;
    TProfile* fTransverse5_other;

    const int LBINS = 20;
    const int LMIN = 0;
    const int LMAX = 5;

    const int TBINS = 20;
    const int TMIN = -5;
    const int TMAX = 5;

    const int EBINS = 20;
    const double EMIN = 0.5;
    const double EMAX = 20.5;

    const double X0 = 14;

    std::string fHitModuleLabel;
    std::string fShowerModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fDigitModuleLabel;

    calo::CalorimetryAlg fCalorimetryAlg;

  }; // class TCShowerAnalysis

} // shower

// -------------------------------------------------

shower::TCShowerTemplateMaker::TCShowerTemplateMaker(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
  , fHitModuleLabel(pset.get<std::string>("HitModuleLabel", "trajcluster"))
  , fShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel", "tcshower"))
  , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel", "generator"))
  , fDigitModuleLabel(pset.get<std::string>("DigitModuleLabel", "largeant"))
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{} // TCShowerTemplateMaker

// -------------------------------------------------

void shower::TCShowerTemplateMaker::beginJob()
{

  art::ServiceHandle<art::TFileService const> tfs;
  fShowerProfileSimLong =
    tfs->make<TProfile>("fShowerProfileSimLong",
                        "longitudinal e- profile (true, simchannel);t;E (MeV)",
                        LBINS,
                        LMIN,
                        LMAX);
  fShowerProfileHitLong = tfs->make<TProfile>(
    "fShowerProfileHitLong", "longitudinal e- profile (true, hit);t;E (MeV)", LBINS, LMIN, LMAX);
  fShowerProfileRecoLong = tfs->make<TProfile>(
    "fShowerProfileRecoLong", "longitudinal e- profile (reco);t;Q", LBINS, LMIN, LMAX);

  fShowerProfileSimLong2D = tfs->make<TProfile2D>(
    "fShowerProfileSimLong2D",
    "longitudinal e- profile (true, simchannel);t;electron energy (MeV);E (MeV)",
    LBINS,
    LMIN,
    LMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileHitLong2D =
    tfs->make<TProfile2D>("fShowerProfileHitLong2D",
                          "longitudinal e- profile (true, hit);t;electron energy (MeV);E (MeV)",
                          LBINS,
                          LMIN,
                          LMAX,
                          EBINS,
                          EMIN,
                          EMAX);
  fShowerProfileRecoLong2D =
    tfs->make<TProfile2D>("fShowerProfileRecoLong2D",
                          "longitudinal e- profile (reco);t;electron energy (MeV);Q",
                          LBINS,
                          LMIN,
                          LMAX,
                          EBINS,
                          EMIN,
                          EMAX);

  fShowerProfileSimTrans =
    tfs->make<TProfile>("fShowerProfileSimTrans",
                        "transverse e- profile (true, simchannel);dist (cm);E (MeV)",
                        TBINS,
                        TMIN,
                        TMAX);
  fShowerProfileHitTrans =
    tfs->make<TProfile>("fShowerProfileHitTrans",
                        "transverse e- profile (true, hit);dist (cm);E (MeV)",
                        TBINS,
                        TMIN,
                        TMAX);
  fShowerProfileRecoTrans = tfs->make<TProfile>(
    "fShowerProfileRecoTrans", "transverse e- profile (reco);dist (cm);Q", TBINS, TMIN, TMAX);

  fShowerProfileSimTrans2D = tfs->make<TProfile2D>(
    "fShowerProfileSimTrans2D",
    "transverse e- profile (true, simchannel);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileHitTrans2D =
    tfs->make<TProfile2D>("fShowerProfileHitTrans2D",
                          "transverse e- profile (true, hit);t;electron energy (MeV);E (MeV)",
                          TBINS,
                          TMIN,
                          TMAX,
                          EBINS,
                          EMIN,
                          EMAX);
  fShowerProfileRecoTrans2D =
    tfs->make<TProfile2D>("fShowerProfileRecoTrans2D",
                          "transverse e- profile (reco);t;electron energy (MeV);Q",
                          TBINS,
                          TMIN,
                          TMAX,
                          EBINS,
                          EMIN,
                          EMAX);

  fShowerProfileSimTrans2D_1 = tfs->make<TProfile2D>(
    "fShowerProfileSimTrans2D_1",
    "transverse e- profile [0 <= t < 1] (true, simchannel);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileHitTrans2D_1 = tfs->make<TProfile2D>(
    "fShowerProfileHitTrans2D_1",
    "transverse e- profile [0 <= t < 1] (true, hit);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileRecoTrans2D_1 =
    tfs->make<TProfile2D>("fShowerProfileRecoTrans2D_1",
                          "transverse e- profile [0 <= t < 1] (reco);t;electron energy (MeV);Q",
                          TBINS,
                          TMIN,
                          TMAX,
                          EBINS,
                          EMIN,
                          EMAX);

  fShowerProfileSimTrans2D_2 = tfs->make<TProfile2D>(
    "fShowerProfileSimTrans2D_2",
    "transverse e- profile [1 <= t < 2] (true, simchannel);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileHitTrans2D_2 = tfs->make<TProfile2D>(
    "fShowerProfileHitTrans2D_2",
    "transverse e- profile [1 <= t < 2] (true, hit);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileRecoTrans2D_2 =
    tfs->make<TProfile2D>("fShowerProfileRecoTrans2D_2",
                          "transverse e- profile [1 <= t < 2] (reco);t;electron energy (MeV);Q",
                          TBINS,
                          TMIN,
                          TMAX,
                          EBINS,
                          EMIN,
                          EMAX);

  fShowerProfileSimTrans2D_3 = tfs->make<TProfile2D>(
    "fShowerProfileSimTrans2D_3",
    "transverse e- profile [2 <= t < 3] (true, simchannel);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileHitTrans2D_3 = tfs->make<TProfile2D>(
    "fShowerProfileHitTrans2D_3",
    "transverse e- profile [2 <= t < 3] (true, hit);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileRecoTrans2D_3 =
    tfs->make<TProfile2D>("fShowerProfileRecoTrans2D_3",
                          "transverse e- profile [2 <= t < 3] (reco);t;electron energy (MeV);Q",
                          TBINS,
                          TMIN,
                          TMAX,
                          EBINS,
                          EMIN,
                          EMAX);

  fShowerProfileSimTrans2D_4 = tfs->make<TProfile2D>(
    "fShowerProfileSimTrans2D_4",
    "transverse e- profile [3 <= t < 4] (true, simchannel);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileHitTrans2D_4 = tfs->make<TProfile2D>(
    "fShowerProfileHitTrans2D_4",
    "transverse e- profile [3 <= t < 4] (true, hit);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileRecoTrans2D_4 =
    tfs->make<TProfile2D>("fShowerProfileRecoTrans2D_4",
                          "transverse e- profile [3 <= t < 4] (reco);t;electron energy (MeV);Q",
                          TBINS,
                          TMIN,
                          TMAX,
                          EBINS,
                          EMIN,
                          EMAX);

  fShowerProfileSimTrans2D_5 = tfs->make<TProfile2D>(
    "fShowerProfileSimTrans2D_5",
    "transverse e- profile [4 <= t < 5] (true, simchannel);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileHitTrans2D_5 = tfs->make<TProfile2D>(
    "fShowerProfileHitTrans2D_5",
    "transverse e- profile [4 <= t < 5] (true, hit);t;electron energy (MeV);E (MeV)",
    TBINS,
    TMIN,
    TMAX,
    EBINS,
    EMIN,
    EMAX);
  fShowerProfileRecoTrans2D_5 =
    tfs->make<TProfile2D>("fShowerProfileRecoTrans2D_5",
                          "transverse e- profile [4 <= t < 5] (reco);t;electron energy (MeV);Q",
                          TBINS,
                          TMIN,
                          TMAX,
                          EBINS,
                          EMIN,
                          EMAX);

  fLongitudinal = tfs->make<TH3F>("fLongitudinal",
                                  "longitudinal e- profile;t;electron energy (MeV);Q",
                                  LBINS,
                                  LMIN,
                                  LMAX,
                                  EBINS,
                                  EMIN,
                                  EMAX,
                                  100,
                                  0,
                                  150000);
  fTransverse = tfs->make<TH3F>("fTransverse",
                                "transverse e- profile;dist (cm);electron energy (MeV);Q",
                                TBINS,
                                TMIN,
                                TMAX,
                                EBINS,
                                EMIN,
                                EMAX,
                                100,
                                0,
                                150000);
  fTransverse_1 =
    tfs->make<TH3F>("fTransverse_1",
                    "transverse e- profile [0 <= t < 1];dist (cm);electron energy (MeV);Q",
                    TBINS,
                    TMIN,
                    TMAX,
                    EBINS,
                    EMIN,
                    EMAX,
                    100,
                    0,
                    100000);
  fTransverse_2 =
    tfs->make<TH3F>("fTransverse_2",
                    "transverse e- profile [1 <= t < 2];dist (cm);electron energy (MeV);Q",
                    TBINS,
                    TMIN,
                    TMAX,
                    EBINS,
                    EMIN,
                    EMAX,
                    100,
                    0,
                    100000);
  fTransverse_3 =
    tfs->make<TH3F>("fTransverse_3",
                    "transverse e- profile [2 <= t < 3];dist (cm);electron energy (MeV);Q",
                    TBINS,
                    TMIN,
                    TMAX,
                    EBINS,
                    EMIN,
                    EMAX,
                    100,
                    0,
                    100000);
  fTransverse_4 =
    tfs->make<TH3F>("fTransverse_4",
                    "transverse e- profile [3 <= t < 4];dist (cm);electron energy (MeV);Q",
                    TBINS,
                    TMIN,
                    TMAX,
                    EBINS,
                    EMIN,
                    EMAX,
                    100,
                    0,
                    100000);
  fTransverse_5 =
    tfs->make<TH3F>("fTransverse_5",
                    "transverse e- profile [4 <= t < 5];dist (cm);electron energy (MeV);Q",
                    TBINS,
                    TMIN,
                    TMAX,
                    EBINS,
                    EMIN,
                    EMAX,
                    100,
                    0,
                    100000);

  // electrons
  fLongitudinal_electron = tfs->make<TProfile>(
    "fLongitudinal_electron", "longitudinal e- profile (reco);t;Q", LBINS, LMIN, LMAX);
  fTransverse1_electron =
    tfs->make<TProfile>("fTransverse1_electron",
                        "transverse e- profile [0 <= t < 1] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse2_electron =
    tfs->make<TProfile>("fTransverse2_electron",
                        "transverse e- profile [1 <= t < 2] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse3_electron =
    tfs->make<TProfile>("fTransverse3_electron",
                        "transverse e- profile [2 <= t < 3] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse4_electron =
    tfs->make<TProfile>("fTransverse4_electron",
                        "transverse e- profile [3 <= t < 4] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse5_electron =
    tfs->make<TProfile>("fTransverse5_electron",
                        "transverse e- profile [4 <= t < 5] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);

  // photons
  fLongitudinal_photon = tfs->make<TProfile>(
    "fLongitudinal_photon", "longitudinal photon profile (reco);t;Q", LBINS, LMIN, LMAX);
  fTransverse1_photon =
    tfs->make<TProfile>("fTransverse1_photon",
                        "transverse photon profile [0 <= t < 1] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse2_photon =
    tfs->make<TProfile>("fTransverse2_photon",
                        "transverse photon profile [1 <= t < 2] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse3_photon =
    tfs->make<TProfile>("fTransverse3_photon",
                        "transverse photon profile [2 <= t < 3] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse4_photon =
    tfs->make<TProfile>("fTransverse4_photon",
                        "transverse photon profile [3 <= t < 4] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse5_photon =
    tfs->make<TProfile>("fTransverse5_photon",
                        "transverse photon profile [4 <= t < 5] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);

  // other
  fLongitudinal_other = tfs->make<TProfile>(
    "fLongitudinal_other", "longitudinal other profile (reco);t;Q", LBINS, LMIN, LMAX);
  fTransverse1_other =
    tfs->make<TProfile>("fTransverse1_other",
                        "transverse other profile [0 <= t < 1] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse2_other =
    tfs->make<TProfile>("fTransverse2_other",
                        "transverse other profile [1 <= t < 2] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse3_other =
    tfs->make<TProfile>("fTransverse3_other",
                        "transverse other profile [2 <= t < 3] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse4_other =
    tfs->make<TProfile>("fTransverse4_other",
                        "transverse other profile [3 <= t < 4] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);
  fTransverse5_other =
    tfs->make<TProfile>("fTransverse5_other",
                        "transverse other profile [4 <= t < 5] (reco);dist (cm);Q",
                        TBINS,
                        TMIN,
                        TMAX);

} // beginJob

// -------------------------------------------------

void shower::TCShowerTemplateMaker::analyze(const art::Event& evt)
{

  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (evt.getByLabel(fHitModuleLabel, hitListHandle)) art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle<std::vector<sim::SimChannel>> scListHandle;
  std::vector<art::Ptr<sim::SimChannel>> simchanlist;
  if (evt.getByLabel(fDigitModuleLabel, scListHandle))
    art::fill_ptr_vector(simchanlist, scListHandle);

  art::Handle<std::vector<recob::Shower>> showerListHandle;
  std::vector<art::Ptr<recob::Shower>> showerlist;
  if (evt.getByLabel(fShowerModuleLabel, showerListHandle))
    art::fill_ptr_vector(showerlist, showerListHandle);

  art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth>> mclist;
  if (evt.getByLabel(fGenieGenModuleLabel, mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  art::FindManyP<recob::Hit> shwfm(showerListHandle, evt, fShowerModuleLabel);

  if (empty(mclist)) return;

  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const det_prop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);

  art::Ptr<simb::MCTruth> mctruth = mclist[0];
  if (mctruth->NeutrinoSet()) {
    if (std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 12 &&
        mctruth->GetNeutrino().CCNC() == 0) {
      double elep = mctruth->GetNeutrino().Lepton().E();
      if (showerlist.size()) {
        std::vector<art::Ptr<recob::Hit>> showerhits = shwfm.at(0);
        showerProfile(clock_data,
                      det_prop,
                      showerhits,
                      showerlist[0]->ShowerStart(),
                      showerlist[0]->Direction(),
                      elep);
      }
      showerProfileTrue(clock_data, det_prop, hitlist, elep);
      showerProfileTrue(simchanlist, mctruth->GetNeutrino().Lepton());
    }
  }
} // analyze

// -------------------------------------------------

void shower::TCShowerTemplateMaker::showerProfile(detinfo::DetectorClocksData const& clockData,
                                                  detinfo::DetectorPropertiesData const& detProp,
                                                  std::vector<art::Ptr<recob::Hit>> showerhits,
                                                  TVector3 shwvtx,
                                                  TVector3 shwdir,
                                                  double elep)
{
  art::ServiceHandle<geo::Geometry const> geom;
  auto const& collectionPlane = geom->Plane(collectionPlaneID);

  double shwVtxTime = detProp.ConvertXToTicks(shwvtx[0], collectionPlaneID);
  using geo::vect::toPoint;
  double shwVtxWire = collectionPlane.WireCoordinate(toPoint(shwvtx));

  double shwTwoTime = detProp.ConvertXToTicks(shwvtx[0] + shwdir[0], collectionPlaneID);
  double shwTwoWire = collectionPlane.WireCoordinate(toPoint(shwvtx + shwdir));

  TH1F* ltemp = new TH1F("ltemp", "ltemp", LBINS, LMIN, LMAX);
  TH1F* ttemp = new TH1F("ttemp", "ttemp", TBINS, TMIN, TMAX);

  TH1F* ttemp_1 = new TH1F("ttemp_1", "ttemp_1", TBINS, TMIN, TMAX);
  TH1F* ttemp_2 = new TH1F("ttemp_2", "ttemp_2", TBINS, TMIN, TMAX);
  TH1F* ttemp_3 = new TH1F("ttemp_3", "ttemp_3", TBINS, TMIN, TMAX);
  TH1F* ttemp_4 = new TH1F("ttemp_4", "ttemp_4", TBINS, TMIN, TMAX);
  TH1F* ttemp_5 = new TH1F("ttemp_5", "ttemp_5", TBINS, TMIN, TMAX);

  for (size_t i = 0; i < showerhits.size(); ++i) {
    if (showerhits[i]->WireID().Plane != collectionPlaneID.Plane) continue;

    double wirePitch = geom->Plane(showerhits[i]->WireID()).WirePitch();
    double tickToDist = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
    tickToDist *= 1.e-3 * sampling_rate(clockData); // 1e-3 is conversion of 1/us to 1/ns

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
               fCalorimetryAlg.LifetimeCorrection(clockData, detProp, showerhits[i]->PeakTime());

    ltemp->Fill(t, Q);
    ttemp->Fill(tdist, Q);

    if (t < 1)
      ttemp_1->Fill(tdist, Q);
    else if (t < 2)
      ttemp_2->Fill(tdist, Q);
    else if (t < 3)
      ttemp_3->Fill(tdist, Q);
    else if (t < 4)
      ttemp_4->Fill(tdist, Q);
    else if (t < 5)
      ttemp_5->Fill(tdist, Q);

  } // loop through showerhits

  for (int i = 0; i < LBINS; ++i) {
    if (ltemp->GetBinContent(i + 1) == 0) continue;
    fShowerProfileRecoLong->Fill(ltemp->GetBinCenter(i + 1), ltemp->GetBinContent(i + 1));
    fShowerProfileRecoLong2D->Fill(ltemp->GetBinCenter(i + 1), elep, ltemp->GetBinContent(i + 1));
    fLongitudinal->Fill(ltemp->GetBinCenter(i + 1), elep, ltemp->GetBinContent(i + 1));
  }

  for (int i = 0; i < TBINS; ++i) {
    if (ttemp->GetBinContent(i + 1) == 0) continue;
    fShowerProfileRecoTrans->Fill(ttemp->GetBinCenter(i + 1), ttemp->GetBinContent(i + 1));
    fShowerProfileRecoTrans2D->Fill(ttemp->GetBinCenter(i + 1), elep, ttemp->GetBinContent(i + 1));
    fTransverse->Fill(ttemp->GetBinCenter(i + 1), elep, ttemp->GetBinContent(i + 1));

    fShowerProfileRecoTrans2D_1->Fill(
      ttemp_1->GetBinCenter(i + 1), elep, ttemp_1->GetBinContent(i + 1));
    fShowerProfileRecoTrans2D_2->Fill(
      ttemp_2->GetBinCenter(i + 1), elep, ttemp_2->GetBinContent(i + 1));
    fShowerProfileRecoTrans2D_3->Fill(
      ttemp_3->GetBinCenter(i + 1), elep, ttemp_3->GetBinContent(i + 1));
    fShowerProfileRecoTrans2D_4->Fill(
      ttemp_4->GetBinCenter(i + 1), elep, ttemp_4->GetBinContent(i + 1));
    fShowerProfileRecoTrans2D_5->Fill(
      ttemp_5->GetBinCenter(i + 1), elep, ttemp_5->GetBinContent(i + 1));

    fTransverse_1->Fill(ttemp_1->GetBinCenter(i + 1), elep, ttemp_1->GetBinContent(i + 1));
    fTransverse_2->Fill(ttemp_2->GetBinCenter(i + 1), elep, ttemp_2->GetBinContent(i + 1));
    fTransverse_3->Fill(ttemp_3->GetBinCenter(i + 1), elep, ttemp_3->GetBinContent(i + 1));
    fTransverse_4->Fill(ttemp_4->GetBinCenter(i + 1), elep, ttemp_4->GetBinContent(i + 1));
    fTransverse_5->Fill(ttemp_5->GetBinCenter(i + 1), elep, ttemp_5->GetBinContent(i + 1));
  }
} // showerProfile

// -------------------------------------------------

void shower::TCShowerTemplateMaker::showerProfileTrue(
  detinfo::DetectorClocksData const& clockData,
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::Hit>> allhits,
  double elep)
{
  art::ServiceHandle<geo::Geometry const> geom;
  auto const& collectionPlane = geom->Plane(collectionPlaneID);
  art::ServiceHandle<cheat::BackTrackerService const> btserv;
  art::ServiceHandle<cheat::ParticleInventoryService const> piserv;
  std::map<int, double> trkID_E;

  TH1F* ltemp = new TH1F("ltemp", "ltemp", LBINS, LMIN, LMAX);
  TH1F* ttemp = new TH1F("ttemp", "ttemp", TBINS, TMIN, TMAX);

  TH1F* ttemp_1 = new TH1F("ttemp_1", "ttemp_1", TBINS, TMIN, TMAX);
  TH1F* ttemp_2 = new TH1F("ttemp_2", "ttemp_2", TBINS, TMIN, TMAX);
  TH1F* ttemp_3 = new TH1F("ttemp_3", "ttemp_3", TBINS, TMIN, TMAX);
  TH1F* ttemp_4 = new TH1F("ttemp_4", "ttemp_4", TBINS, TMIN, TMAX);
  TH1F* ttemp_5 = new TH1F("ttemp_5", "ttemp_5", TBINS, TMIN, TMAX);

  double xvtx = -999;
  double yvtx = -999;
  double zvtx = -999;
  double xtwo = -999;
  double ytwo = -999;
  double ztwo = -999;
  double shwvtxT = -999;
  double shwvtxW = -999;
  double shwtwoT = -999;
  double shwtwoW = -999;

  double shwvtxx = -999;
  double shwvtxy = -999;
  double shwtwox = -999;
  double shwtwoy = -999;
  double xtwoorth = -999;
  double ytwoorth = -999;

  double wirePitch = -999;
  double tickToDist = -999;

  bool foundParent = false;

  for (auto const& hit : allhits) {
    if (hit->WireID().Plane != collectionPlaneID.Plane) continue;

    // art::Ptr<recob::Hit> hit = allhits[i];
    std::vector<sim::TrackIDE> trackIDs = btserv->HitToEveTrackIDEs(clockData, hit);

    for (size_t j = 0; j < trackIDs.size(); ++j) {
      // only want energy associated with the electron and electron must have neutrino mother
      if (std::abs((piserv->TrackIdToParticle_P(trackIDs[j].trackID))->PdgCode()) != 11) continue;
      if (std::abs((piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Mother()) != 0) continue;

      if (!foundParent) {
        xvtx = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Vx();
        yvtx = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Vy();
        zvtx = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Vz();

        xtwo = xvtx + (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Px();
        ytwo = yvtx + (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Py();
        ztwo = zvtx + (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Pz();

        shwvtxT = detProp.ConvertXToTicks(xvtx, collectionPlaneID);
        geo::Point_t const vtx{xvtx, yvtx, zvtx};
        shwvtxW = collectionPlane.WireCoordinate(vtx);

        shwtwoT = detProp.ConvertXToTicks(xtwo, collectionPlaneID);
        geo::Point_t const vtwo{xtwo, ytwo, ztwo};
        shwtwoW = collectionPlane.WireCoordinate(vtwo);

        wirePitch = geom->Plane(hit->WireID()).WirePitch();
        tickToDist = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
        tickToDist *= 1.e-3 * sampling_rate(clockData); // 1e-3 is conversion of 1/us to 1/ns

        shwvtxx = shwvtxT * tickToDist;
        shwvtxy = shwvtxW * wirePitch;

        shwtwox = shwtwoT * tickToDist;
        shwtwoy = shwtwoW * wirePitch;

        xtwoorth = (shwtwoy - shwvtxy) + shwvtxx;
        ytwoorth = -(shwtwox - shwvtxx) + shwvtxy;

        foundParent = true;
      }
      double xhit = hit->PeakTime() * tickToDist;
      double yhit = hit->WireID().Wire * wirePitch;

      double ldist = abs((ytwoorth - shwvtxy) * xhit - (xtwoorth - shwvtxx) * yhit +
                         xtwoorth * shwvtxy - ytwoorth * shwvtxx) /
                     std::hypot(ytwoorth - shwvtxy, xtwoorth - shwvtxx);
      double tdist = ((shwtwoy - shwvtxy) * xhit - (shwtwox - shwvtxx) * yhit + shwtwox * shwvtxy -
                      shwtwoy * shwvtxx) /
                     std::hypot(shwtwoy - shwvtxy, shwtwox - shwvtxx);

      double to3D = std::hypot(xvtx - xtwo, yvtx - ytwo, zvtx - ztwo) /
                    std::hypot(shwvtxx - shwtwox,
                               shwvtxy - shwtwoy); // distance between two points in 3D space is one
      ldist *= to3D;
      tdist *= to3D;
      double t = ldist / X0;

      double energy = trackIDs[j].energy;

      ltemp->Fill(t, energy);
      ttemp->Fill(tdist, energy);

      if (t < 1)
        ttemp_1->Fill(tdist, energy);
      else if (t < 2)
        ttemp_2->Fill(tdist, energy);
      else if (t < 3)
        ttemp_3->Fill(tdist, energy);
      else if (t < 4)
        ttemp_4->Fill(tdist, energy);
      else if (t < 5)
        ttemp_5->Fill(tdist, energy);

    } // loop through track IDE

  } // loop through all hits

  for (int i = 0; i < LBINS; ++i) {
    if (ltemp->GetBinContent(i + 1) == 0) continue;
    fShowerProfileHitLong->Fill(ltemp->GetBinCenter(i + 1), ltemp->GetBinContent(i + 1));
    fShowerProfileHitLong2D->Fill(ltemp->GetBinCenter(i + 1), elep, ltemp->GetBinContent(i + 1));
  }

  for (int i = 0; i < TBINS; ++i) {
    if (ttemp->GetBinContent(i + 1) == 0) continue;
    fShowerProfileHitTrans->Fill(ttemp->GetBinCenter(i + 1), ttemp->GetBinContent(i + 1));
    fShowerProfileHitTrans2D->Fill(ttemp->GetBinCenter(i + 1), elep, ttemp->GetBinContent(i + 1));

    fShowerProfileHitTrans2D_1->Fill(
      ttemp_1->GetBinCenter(i + 1), elep, ttemp_1->GetBinContent(i + 1));
    fShowerProfileHitTrans2D_2->Fill(
      ttemp_2->GetBinCenter(i + 1), elep, ttemp_2->GetBinContent(i + 1));
    fShowerProfileHitTrans2D_3->Fill(
      ttemp_3->GetBinCenter(i + 1), elep, ttemp_3->GetBinContent(i + 1));
    fShowerProfileHitTrans2D_4->Fill(
      ttemp_4->GetBinCenter(i + 1), elep, ttemp_4->GetBinContent(i + 1));
    fShowerProfileHitTrans2D_5->Fill(
      ttemp_5->GetBinCenter(i + 1), elep, ttemp_5->GetBinContent(i + 1));
  }
} // showerProfileTrue

// -------------------------------------------------

void shower::TCShowerTemplateMaker::showerProfileTrue(
  std::vector<art::Ptr<sim::SimChannel>> allchan,
  simb::MCParticle electron)
{
  art::ServiceHandle<cheat::ParticleInventoryService const> piserv;
  auto const* channelMapAlg =
    art::ServiceHandle<geo::ExptGeoHelperInterface const>()->ChannelMapAlgPtr();

  std::vector<sim::MCEnDep> alledep;

  TH1F* ltemp = new TH1F("ltemp", "ltemp", LBINS, LMIN, LMAX);
  TH1F* ttemp = new TH1F("ttemp", "ttemp", TBINS, TMIN, TMAX);

  TH1F* ttemp_1 = new TH1F("ttemp_1", "ttemp_1", TBINS, TMIN, TMAX);
  TH1F* ttemp_2 = new TH1F("ttemp_2", "ttemp_2", TBINS, TMIN, TMAX);
  TH1F* ttemp_3 = new TH1F("ttemp_3", "ttemp_3", TBINS, TMIN, TMAX);
  TH1F* ttemp_4 = new TH1F("ttemp_4", "ttemp_4", TBINS, TMIN, TMAX);
  TH1F* ttemp_5 = new TH1F("ttemp_5", "ttemp_5", TBINS, TMIN, TMAX);

  // get electron energy depositions
  for (size_t i = 0; i < allchan.size(); ++i) {
    art::Ptr<sim::SimChannel> simchan = allchan[i];
    if (channelMapAlg->View(simchan->Channel()) != geo::kV) continue;
    auto tdc_ide_map = simchan->TDCIDEMap();

    for (auto const& tdc_ide_pair : tdc_ide_map) {
      for (auto const& ide : tdc_ide_pair.second) {
        if (piserv->TrackIdToMotherParticle_P(ide.trackID) == nullptr) continue;
        if (std::abs(piserv->TrackIdToMotherParticle_P(ide.trackID)->PdgCode()) != 11) continue;

        sim::MCEnDep edep;
        edep.SetVertex(ide.x, ide.y, ide.z);
        edep.SetEnergy(ide.energy);
        edep.SetTrackId(ide.trackID);

        alledep.push_back(edep);

      } // loop through ide
    }   // loop through tdc_ide

  } // loop through channels

  double x0 = electron.Vx();
  double y0 = electron.Vy();
  double z0 = electron.Vz();

  double x2 = electron.Px();
  double y2 = electron.Py();
  double z2 = electron.Pz();

  TVector3 v0(x2, y2, z2);
  v0 = v0.Unit();

  for (size_t i = 0; i < alledep.size(); ++i) {
    double x = (double)alledep[i].Vertex()[0];
    double y = (double)alledep[i].Vertex()[1];
    double z = (double)alledep[i].Vertex()[2];

    TVector3 v1(x - x0, y - y0, z - z0);

    double ldist = v0.Dot(v1);
    double t = ldist / X0;
    double tdist = (v0.Orthogonal()).Dot(v1);

    double energy = alledep[i].Energy();

    ltemp->Fill(t, energy);
    ttemp->Fill(tdist, energy);

    if (t < 1)
      ttemp_1->Fill(tdist, energy);
    else if (t < 2)
      ttemp_2->Fill(tdist, energy);
    else if (t < 3)
      ttemp_3->Fill(tdist, energy);
    else if (t < 4)
      ttemp_4->Fill(tdist, energy);
    else if (t < 5)
      ttemp_5->Fill(tdist, energy);
  }

  for (int i = 0; i < LBINS; ++i) {
    if (ltemp->GetBinContent(i + 1) == 0) continue;
    fShowerProfileSimLong->Fill(ltemp->GetBinCenter(i + 1), ltemp->GetBinContent(i + 1));
    fShowerProfileSimLong2D->Fill(
      ltemp->GetBinCenter(i + 1), electron.E(), ltemp->GetBinContent(i + 1));
  }

  for (int i = 0; i < TBINS; ++i) {
    if (ttemp->GetBinContent(i + 1) == 0) continue;
    fShowerProfileSimTrans->Fill(ttemp->GetBinCenter(i + 1), ttemp->GetBinContent(i + 1));
    fShowerProfileSimTrans2D->Fill(
      ttemp->GetBinCenter(i + 1), electron.E(), ttemp->GetBinContent(i + 1));

    fShowerProfileSimTrans2D_1->Fill(
      ttemp_1->GetBinCenter(i + 1), electron.E(), ttemp_1->GetBinContent(i + 1));
    fShowerProfileSimTrans2D_2->Fill(
      ttemp_2->GetBinCenter(i + 1), electron.E(), ttemp_2->GetBinContent(i + 1));
    fShowerProfileSimTrans2D_3->Fill(
      ttemp_3->GetBinCenter(i + 1), electron.E(), ttemp_3->GetBinContent(i + 1));
    fShowerProfileSimTrans2D_4->Fill(
      ttemp_4->GetBinCenter(i + 1), electron.E(), ttemp_4->GetBinContent(i + 1));
    fShowerProfileSimTrans2D_5->Fill(
      ttemp_5->GetBinCenter(i + 1), electron.E(), ttemp_5->GetBinContent(i + 1));
  }
} // showerProfileTrue

// -------------------------------------------------

DEFINE_ART_MODULE(shower::TCShowerTemplateMaker)

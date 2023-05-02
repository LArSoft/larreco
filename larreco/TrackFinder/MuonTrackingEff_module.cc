//
// **Muon Tracking Efficiency module**
// This module is based on NeutrinoTrackingEff_module.cc from A. Higuera.
// I have changed the module so that it only calculates the Muon Tracking
// Efficiency by checking the Completeness and Purity (see function: "truthMatcher")
// and track length (see function: "truthLength") of the leading reconstructed muon
// track (This is the one with the highest Completeness).
// In case the leading muon track failed and there is more than one reconstructed
// muon track (e.g. due tue to track splitting in the reconstruction bc of a kink
// after multiple scattering), I check the efficiency criteria for the sum of the
// leading and second reconstructed muon track.
//
// Christoph Alt
// christoph.alt@cern.ch

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TVector3.h"

#define MAX_TRACKS 1000
using namespace std;

//========================================================================

namespace DUNE {

  class MuonTrackingEff : public art::EDAnalyzer {
  public:
    explicit MuonTrackingEff(fhicl::ParameterSet const& pset);

  private:
    void beginJob() override;
    void endJob() override;
    void beginRun(const art::Run& run) override;
    void analyze(const art::Event& evt) override;

    void processEff(const art::Event& evt, bool& isFiducial);

    void truthMatcher(detinfo::DetectorClocksData const& clockData,
                      std::vector<recob::Hit> const& AllHits,
                      std::vector<art::Ptr<recob::Hit>> track_hits,
                      const simb::MCParticle*& MCparticle,
                      double& Purity,
                      double& Completeness,
                      double& TotalRecoEnergy);

    void FuncDistanceAndAngleBetweenTracks(art::Ptr<recob::Track> Track1,
                                           art::Ptr<recob::Track> Track2,
                                           double& TempDistanceBetweenTracks,
                                           double& TempAngleBetweenTracks,
                                           double& TempCriteriaTwoTracks);

    void FuncDistanceAndAngleBetweenTruthAndRecoTrack(const simb::MCParticle*& MCparticle,
                                                      art::Ptr<recob::Track> Track,
                                                      double& TempDistanceBetweenTruthAndRecoTrack,
                                                      double& TempAngleBeetweenTruthAndRecoTrack);

    double truthLength(const simb::MCParticle* MCparticle);

    bool insideFV(double vertex[4]);

    void doEfficiencies();

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fTrackModuleLabel;
    int fMuonPDGCode;

    //  int    MC_isCC;
    //  int    MC_incoming_PDG;
    //  double MC_incoming_P[4];
    double MCTruthMuonVertex[4];
    //  double MCTruthMuonStartMomentum[4];

    int MCTruthMuonID;
    double MCTruthMuonMomentum;

    //  double MCTruthMuonMomentumXZ=0;
    //  double MCTruthMuonMomentumYZ=0;
    double MCTruthMuonThetaXZ = 0;
    double MCTruthMuonThetaYZ = 0;

    //Counts to calculate efficiency and failed criteria
    int EventCounter = 0;

    int CountMCTruthMuon = 0;
    int CountRecoMuon = 0;

    int CountGoodLeadingMuonTrack = 0;
    int CountNoRecoTracks = 0;
    int CountNoMuonTracks = 0;
    int CountBadLeadingMuonTrack = 0;
    int CountCompleteness = 0;
    int CountPurity = 0;
    int CountTrackLengthTooShort = 0;
    int CountTrackLengthTooLong = 0;

    int CountBadLeadingMuonTrackAndOnlyOneMuonTrack = 0;
    int CountBadLeadingMuonTrackButLeadingPlusSecondGood = 0;
    int CountBadLeadingMuonTrackAndLeadingPlusSecondBad = 0;

    int CountBadLeadingMuonTrackAndLeadingPlusSecondBadCompleteness = 0;
    int CountBadLeadingMuonTrackAndLeadingPlusSecondBadPurity = 0;
    int CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooShort = 0;
    int CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooLong = 0;

    int CountBadLeadingMuonTrackAndOnlyOneMuonTrackCompleteness = 0;
    int CountBadLeadingMuonTrackAndOnlyOneMuonTrackPurity = 0;
    int CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooShort = 0;
    int CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooLong = 0;

    double Criteria;
    //  int NMuonTracksTooShort=0;

    int GoodEvents1MuonTrack = 0;
    int GoodEvents2MuonTrack = 0;
    int GoodEvents3MuonTrack = 0;
    int GoodEvents4OrMoreMuonTrack = 0;

    int BadEvents0MuonTrack = 0;
    int BadEvents1MuonTrack = 0;
    int BadEvents2MuonTrack = 0;
    int BadEvents3MuonTrack = 0;
    int BadEvents4OrMoreMuonTrack = 0;

    //TH1Ds
    //Single Criteria and total reco energy
    TH1D* h_Purity;
    TH1D* h_Completeness;
    TH1D* h_TrackRes;
    TH1D* h_TotalRecoEnergy;
    TH1D* h_TruthLength;
    TH1D* h_VertexRes;
    TH1D* h_DirectionRes;

    //Efficiency ThetaXZ
    TH1D* h_Efficiency_ThetaXZ;
    TH1D* h_ThetaXZ_den;
    TH1D* h_ThetaXZ_num;

    //Efficiency ThetaYZ
    TH1D* h_Efficiency_ThetaYZ;
    TH1D* h_ThetaYZ_den;
    TH1D* h_ThetaYZ_num;

    //Efficiency SinThetaYZ
    TH1D* h_Efficiency_SinThetaYZ;
    TH1D* h_SinThetaYZ_den;
    TH1D* h_SinThetaYZ_num;

    //TH2Ds
    //Efficiency ThetaXZ vs. ThetaYZ
    TH2D* h_Efficiency_ThetaXZ_ThetaYZ;
    TH2D* h_ThetaXZ_ThetaYZ_den;
    TH2D* h_ThetaXZ_ThetaYZ_num;
    TH2D* h_FailedReconstruction_ThetaXZ_ThetaYZ;

    //Efficiency ThetaXZ vs. SinThetaYZ
    TH2D* h_Efficiency_ThetaXZ_SinThetaYZ;
    TH2D* h_ThetaXZ_SinThetaYZ_den;
    TH2D* h_ThetaXZ_SinThetaYZ_num;
    TH2D* h_FailedReconstruction_ThetaXZ_SinThetaYZ;

    //Efficiency ThetaXZ vs. ThetaYZ after summing up leading plus second
    TH2D* h_Efficiency_ThetaXZ_ThetaYZ_LeadingPlusSecond;
    TH2D* h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk;
    TH2D* h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk_num;

    //Difference in efficiency before and after summing up leading plus second: ThetaXZ vs. ThetaYZ
    TH2D* h_Efficiency_ThetaXZ_ThetaYZ_DifferenceLeadingAndLeadingPlusSecond;

    //Efficiency ThetaXZ vs. SinThetaYZ after summing up leading plus second
    TH2D* h_Efficiency_ThetaXZ_SinThetaYZ_LeadingPlusSecond;
    TH2D* h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk;
    TH2D* h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk_num;

    //Failed Criteria
    TH2D* h_NoRecoTrackAtAll_ThetaXZ_SinThetaYZ;
    TH2D* h_NoMuonTrack_ThetaXZ_SinThetaYZ;
    TH2D* h_TrackTooShort_ThetaXZ_SinThetaYZ;
    TH2D* h_TrackTooLong_ThetaXZ_SinThetaYZ;
    TH2D* h_Completeness_ThetaXZ_SinThetaYZ;
    TH2D* h_Purity_ThetaXZ_SinThetaYZ;

    //Criteria vs. NRecoTrack
    TH2D* h_Criteria_NRecoTrack;
    TH2D* h_Criteria_NRecoTrack_den;
    TH2D* h_Criteria_NRecoTrack_num;

    //Criteria vs. NMuonTrack
    TH2D* h_Criteria_NMuonTrack;
    TH2D* h_Criteria_NMuonTrack_den;
    TH2D* h_Criteria_NMuonTrack_num;

    //NoMuonTrack: Max length of no muon track vs. PDG code of that track (MC truth)
    TH2D* h_NoMuonTrack_MaxTrackLength_PDGCode;

    //Stitching variables: all events
    TH2D* h_MuonTrackStitching_TrackRes_Completeness;
    TH2D* h_MuonTrackStitching_TrackResLeading_TrackResSecond;
    TH2D* h_MuonTrackStitching_Distance_Angle;
    TH2D* h_MuonTrackStitching_TrackResSecondMuon_Angle;
    TH2D* h_MuonTrackStitching_CompletenessSecondMuon_Angle;
    TH2D* h_MuonTrackStitching_CriteriaTwoTracks_Angle;

    //Stitching variables: bad events
    TH2D* h_MuonTrackStitching_FailedCriteria_TrackRes_Completeness;
    TH2D* h_MuonTrackStitching_FailedCriteria_TrackResLeading_TrackResSecond;
    TH2D* h_MuonTrackStitching_FailedCriteria_CompletenessLeading_CompletenessSecond;
    TH2D* h_MuonTrackStitching_FailedCriteria_Distance_Angle;
    TH2D* h_MuonTrackStitching_FailedCriteria_TrackResSecondMuon_Angle;
    TH2D* h_MuonTrackStitching_FailedCriteria_CompletenessSecondMuon_Angle;
    TH2D* h_MuonTrackStitching_FailedCriteria_CriteriaTwoTracks_Angle;

    //Stitching variables: good events
    TH2D* h_MuonTrackStitching_MatchedCriteria_TrackRes_Completeness;
    TH2D* h_MuonTrackStitching_MatchedCriteria_TrackResLeading_TrackResSecond;
    TH2D* h_MuonTrackStitching_MatchedCriteria_Distance_Angle;
    TH2D* h_MuonTrackStitching_MatchedCriteria_CriteriaTwoTracks_Angle;

    //Stitching variables: bad events but leading plus second ok
    TH2D* h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackRes_Completeness;
    TH2D*
      h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackResLeading_TrackResSecond;
    TH2D* h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_Distance_Angle;

    //Stitching variables: bad events and leading plus second not ok
    TH2D* h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackRes_Completeness;
    TH2D* h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackResLeading_TrackResSecond;
    TH2D* h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_Distance_Angle;

    float fFidVolCutX;
    float fFidVolCutY;
    float fFidVolCutZ;

    float fFidVolXmin;
    float fFidVolXmax;
    float fFidVolYmin;
    float fFidVolYmax;
    float fFidVolZmin;
    float fFidVolZmax;

    art::ServiceHandle<geo::Geometry const> geom;

    //My histograms
    int NThetaXZBins = 36;
    int ThetaXZBinMin = 0;
    int ThetaXZBinMax = 360;

    int NThetaYZBins = 18;
    int ThetaYZBinMin = -90;
    int ThetaYZBinMax = 90;

    int NSinThetaYZBins = 18;
    int SinThetaYZBinMin = -1;
    int SinThetaYZBinMax = 1;

    int NCriteriaBins = 13;
    double CriteriaBinMin = -0.25;
    double CriteriaBinMax = 6.25;

    int NRecoTracksBins = 19;
    double RecoTracksBinMin = -0.25;
    double RecoTracksBinMax = 9.25;

  }; // class MuonTrackingEff

  //========================================================================
  MuonTrackingEff::MuonTrackingEff(fhicl::ParameterSet const& p) : EDAnalyzer(p)
  {
    fMCTruthModuleLabel = p.get<std::string>("MCTruthModuleLabel");
    fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
    fMuonPDGCode = p.get<int>("MuonPDGCode");
    fFidVolCutX = p.get<float>("FidVolCutX");
    fFidVolCutY = p.get<float>("FidVolCutY");
    fFidVolCutZ = p.get<float>("FidVolCutZ");
  }
  //========================================================================
  void MuonTrackingEff::beginJob()
  {
    std::cout << "job begin..." << std::endl;
    // Get geometry.
    auto const* geo = lar::providerFrom<geo::Geometry>();
    // Define histogram boundaries (cm).
    // For now only draw cryostat=0.
    double minx = 1e9;
    double maxx = -1e9;
    double miny = 1e9;
    double maxy = -1e9;
    double minz = 1e9;
    double maxz = -1e9;
    for (auto const& tpc : geo->Iterate<geo::TPCGeo>(geo::CryostatID{0})) {
      auto const world = tpc.GetCenter();
      if (minx > world.X() - tpc.HalfWidth()) minx = world.X() - tpc.HalfWidth();
      if (maxx < world.X() + tpc.HalfWidth()) maxx = world.X() + tpc.HalfWidth();
      if (miny > world.Y() - tpc.HalfHeight()) miny = world.Y() - tpc.HalfHeight();
      if (maxy < world.Y() + tpc.HalfHeight()) maxy = world.Y() + tpc.HalfHeight();
      if (minz > world.Z() - tpc.Length() / 2.) minz = world.Z() - tpc.Length() / 2.;
      if (maxz < world.Z() + tpc.Length() / 2.) maxz = world.Z() + tpc.Length() / 2.;
    }

    fFidVolXmin = minx + fFidVolCutX;
    fFidVolXmax = maxx - fFidVolCutX;
    fFidVolYmin = miny + fFidVolCutY;
    fFidVolYmax = maxy - fFidVolCutY;
    fFidVolZmin = minz + fFidVolCutZ;
    fFidVolZmax = maxz - fFidVolCutZ;

    std::cout << "Fiducial volume:"
              << "\n"
              << fFidVolXmin << "\t< x <\t" << fFidVolXmax << "\n"
              << fFidVolYmin << "\t< y <\t" << fFidVolYmax << "\n"
              << fFidVolZmin << "\t< z <\t" << fFidVolZmax << "\n";

    art::ServiceHandle<art::TFileService const> tfs;

    //TH1D's
    gStyle->SetTitleOffset(1.3, "Y");

    //Single Criteria and total reco energy
    h_Purity =
      tfs->make<TH1D>("h_Purity", "All events: Purity vs. # events; Purity; # events", 60, 0, 1.2);

    h_Completeness = tfs->make<TH1D>(
      "h_Completeness", "All events: Completeness vs # events; Completeness; # events", 60, 0, 1.2);
    h_Completeness->SetLineColor(kBlue);

    h_TrackRes =
      tfs->make<TH1D>("h_TrackRes",
                      "All events: L_{reco}/L_{truth} vs. # events; L_{reco}/L_{truth}; # events;",
                      75,
                      0,
                      1.5);
    h_TrackRes->SetLineColor(kRed);

    h_TotalRecoEnergy = tfs->make<TH1D>("h_TotalRecoEnergy",
                                        "All events: Total reco energy (sum of all hits in all "
                                        "tracks) vs. # events; E_{reco., tot.} [MeV]; # events",
                                        100,
                                        0,
                                        1000);

    h_TruthLength = tfs->make<TH1D>(
      "h_TruthLength",
      "All events: truth muon length vs. # events; truth muon length [cm]; # events",
      100,
      0,
      300);

    h_VertexRes = tfs->make<TH1D>(
      "h_VertexRes",
      "All events: Vertex residuals vs. # events; #Delta vertex_{truth-teco} [cm]; # events",
      300,
      0,
      300);

    h_DirectionRes = tfs->make<TH1D>(
      "h_DirectionRes",
      "All events: Angular residuals vs. # events; #Delta#theta_{truth-reco} [#circ]; # events",
      180,
      0,
      180);

    //Efficiency ThetaXZ
    h_Efficiency_ThetaXZ =
      tfs->make<TH1D>("h_Efficiency_ThetaXZ",
                      "Muon reco efficiency vs. #theta_{XZ}; #theta_{XZ} [#circ]; Efficiency",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax);
    h_ThetaXZ_den =
      tfs->make<TH1D>("h_ThetaXZ_den",
                      "# generated muons vs. #theta_{XZ}; #theta_{XZ} [#circ]; # generated muons",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax);
    h_ThetaXZ_num = tfs->make<TH1D>(
      "h_ThetaXZ_num",
      "# reconstructed muons vs. #theta_{XZ}; #theta_{XZ} [#circ]; # reconstructed muons",
      NThetaXZBins,
      ThetaXZBinMin,
      ThetaXZBinMax);

    //Efficiency ThetaYZ
    h_Efficiency_ThetaYZ = tfs->make<TH1D>(
      "h_Efficiency_ThetaYZ",
      "Muon reco efficiency vs. #theta_{YZ}; #theta_{YZ} [#circ]; Muon reco efficiency",
      NThetaYZBins,
      ThetaYZBinMin,
      ThetaYZBinMax);
    ;
    h_ThetaYZ_den =
      tfs->make<TH1D>("h_ThetaYZ_den",
                      "# generated muons vs. #theta_{YZ}; #theta_{YZ} [#circ]; # generated muons",
                      NThetaYZBins,
                      ThetaYZBinMin,
                      ThetaYZBinMax);
    h_ThetaYZ_num = tfs->make<TH1D>(
      "h_ThetaYZ_num",
      "# reconstructed muons vs. #theta_{YZ}; #theta_{YZ} [#circ]; # reconstructed muons",
      NThetaYZBins,
      ThetaYZBinMin,
      ThetaYZBinMax);

    //Efficiency SinThetaYZ
    h_Efficiency_SinThetaYZ = tfs->make<TH1D>(
      "h_Efficiency_SinThetaYZ",
      "Muon reco efficiency vs. sin(#theta_{YZ}); sin(#theta_{YZ}); Muon reco efficiency",
      NSinThetaYZBins,
      SinThetaYZBinMin,
      SinThetaYZBinMax);
    ;
    h_SinThetaYZ_den =
      tfs->make<TH1D>("h_SinThetaYZ_den",
                      "# generated muons vs. sin(#theta_{YZ}); sin(#theta_{YZ}); # generated muons",
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_SinThetaYZ_num = tfs->make<TH1D>(
      "h_SinThetaYZ_num",
      "# reconstructed muons vs. sin(#theta_{YZ}); sin(#theta_{YZ}); # reconstructed muons",
      NSinThetaYZBins,
      SinThetaYZBinMin,
      SinThetaYZBinMax);

    h_Purity->Sumw2();
    h_Completeness->Sumw2();
    h_TrackRes->Sumw2();
    h_TotalRecoEnergy->Sumw2();
    h_TruthLength->Sumw2();
    h_VertexRes->Sumw2();
    h_DirectionRes->Sumw2();

    h_Efficiency_SinThetaYZ->Sumw2();
    h_SinThetaYZ_den->Sumw2();
    h_SinThetaYZ_num->Sumw2();

    h_Efficiency_ThetaXZ->Sumw2();
    h_ThetaXZ_den->Sumw2();
    h_ThetaXZ_num->Sumw2();

    h_Efficiency_ThetaYZ->Sumw2();
    h_ThetaYZ_den->Sumw2();
    h_ThetaYZ_num->Sumw2();

    //TH2D's
    //Efficiency and Failed Reconstruction ThetaXZ vs. ThetaYZ
    h_Efficiency_ThetaXZ_ThetaYZ = tfs->make<TH2D>(
      "h_Efficiency_ThetaXZ_ThetaYZ",
      "Muon reco efficiency: #theta_{XZ} vs. #theta_{YZ}; #theta_{XZ} [#circ]; #theta_{YZ} [#circ]",
      NThetaXZBins,
      ThetaXZBinMin,
      ThetaXZBinMax,
      NThetaYZBins,
      ThetaYZBinMin,
      ThetaYZBinMax);
    h_Efficiency_ThetaXZ_ThetaYZ->SetOption("colz");

    h_ThetaXZ_ThetaYZ_den = tfs->make<TH2D>(
      "h_ThetaXZ_ThetaYZ_den",
      "# generated muons: #theta_{XZ} vs. #theta_{YZ}; #theta_{XZ} [#circ]; #theta_{YZ} [#circ]",
      NThetaXZBins,
      ThetaXZBinMin,
      ThetaXZBinMax,
      NThetaYZBins,
      ThetaYZBinMin,
      ThetaYZBinMax);
    h_ThetaXZ_ThetaYZ_den->SetOption("colz");

    h_ThetaXZ_ThetaYZ_num = tfs->make<TH2D>("h_ThetaXZ_ThetaYZ_num",
                                            "# reconstructed muons: #theta_{XZ} vs. #theta_{YZ}; "
                                            "#theta_{XZ} [#circ]; #theta_{YZ} [#circ]",
                                            NThetaXZBins,
                                            ThetaXZBinMin,
                                            ThetaXZBinMax,
                                            NThetaYZBins,
                                            ThetaYZBinMin,
                                            ThetaYZBinMax);
    h_ThetaXZ_ThetaYZ_num->SetOption("colz");

    h_FailedReconstruction_ThetaXZ_ThetaYZ =
      tfs->make<TH2D>("h_FailedReconstruction_ThetaXZ_ThetaYZ",
                      "# failed reconstructions: #theta_{XZ} vs. #theta_{YZ}; #theta_{XZ} [#circ]; "
                      "#theta_{YZ} [#circ]",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NThetaYZBins,
                      ThetaYZBinMin,
                      ThetaYZBinMax);
    h_FailedReconstruction_ThetaXZ_ThetaYZ->SetOption("colz");

    //Efficiency and Failed Reconstruction ThetaXZ vs. SinThetaYZ
    h_Efficiency_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_Efficiency_ThetaXZ_SinThetaYZ",
                      "Muon reco efficiency: #theta_{XZ} vs. sin(#theta_{YZ}); #theta_{XZ} "
                      "[#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_Efficiency_ThetaXZ_SinThetaYZ->SetOption("colz");

    h_ThetaXZ_SinThetaYZ_den = tfs->make<TH2D>(
      "h_ThetaXZ_SinThetaYZ_den",
      "# generated muons: #theta_{XZ} vs. sin(#theta_{YZ}); #theta_{XZ} [#circ]; sin(#theta_{YZ})",
      NThetaXZBins,
      ThetaXZBinMin,
      ThetaXZBinMax,
      NSinThetaYZBins,
      SinThetaYZBinMin,
      SinThetaYZBinMax);
    h_ThetaXZ_SinThetaYZ_den->SetOption("colz");

    h_ThetaXZ_SinThetaYZ_num =
      tfs->make<TH2D>("h_ThetaXZ_SinThetaYZ_num",
                      "# reconstructed muons: #theta_{XZ} vs. sin(#theta_{YZ}); #theta_{XZ} "
                      "[#circ]; sin(#theta_{YZ}); #theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_ThetaXZ_SinThetaYZ_num->SetOption("colz");

    h_FailedReconstruction_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_FailedReconstruction_ThetaXZ_SinThetaYZ",
                      "# failed reconstructions: #theta_{XZ} vs. sin(#theta_{YZ}); #theta_{XZ} "
                      "[#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_FailedReconstruction_ThetaXZ_SinThetaYZ->SetOption("colz");

    //Efficiency ThetaXZ vs. ThetaYZ after summing up leading plus second
    h_Efficiency_ThetaXZ_ThetaYZ_LeadingPlusSecond =
      tfs->make<TH2D>("h_Efficiency_ThetaXZ_ThetaYZ_LeadingPlusSecond",
                      "Muon reco efficiency after stitching: #theta_{XZ} vs. #theta_{YZ}; "
                      "#theta_{XZ} [#circ]; #theta_{YZ} [#circ]",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NThetaYZBins,
                      ThetaYZBinMin,
                      ThetaYZBinMax);
    h_Efficiency_ThetaXZ_ThetaYZ_LeadingPlusSecond->SetOption("colz");

    h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk =
      tfs->make<TH2D>("h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk",
                      "# reconstructed muons after stitching (failed before stitching): "
                      "#theta_{XZ} vs #theta_{YZ}; #theta_{XZ} [#circ]; #theta_{YZ} [#circ]",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NThetaYZBins,
                      ThetaYZBinMin,
                      ThetaYZBinMax);
    h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk->SetOption("colz");

    h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk_num =
      tfs->make<TH2D>("h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk_num",
                      "# reconstructed muons after stitching: #theta_{XZ} vs. #theta_{YZ}; "
                      "#theta_{XZ} [#circ]; #theta_{YZ} [#circ]",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NThetaYZBins,
                      ThetaYZBinMin,
                      ThetaYZBinMax);
    h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk_num->SetOption("colz");

    //Efficiency ThetaXZ vs. SinThetaYZ after summing up leading plus second
    h_Efficiency_ThetaXZ_SinThetaYZ_LeadingPlusSecond =
      tfs->make<TH2D>("h_Efficiency_ThetaXZ_SinThetaYZ_LeadingPlusSecond",
                      "Muon reco efficiency after stitching: #theta_{XZ} vs. sin(#theta_{YZ}); "
                      "#theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_Efficiency_ThetaXZ_SinThetaYZ_LeadingPlusSecond->SetOption("colz");

    h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk =
      tfs->make<TH2D>("h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk",
                      "# reconstructed muons after stitching (failed before stitching): "
                      "#theta_{XZ} vs. sin(#theta_{YZ}); #theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk->SetOption("colz");

    h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk_num =
      tfs->make<TH2D>("h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk_num",
                      "# reconstructed muons after stitching: #theta_{XZ} vs. sin(#theta_{YZ}); "
                      "#theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk_num->SetOption("colz");

    //Difference in efficiency before and after summing up leading plus second: ThetaXZ vs. ThetaYZ
    h_Efficiency_ThetaXZ_ThetaYZ_DifferenceLeadingAndLeadingPlusSecond =
      tfs->make<TH2D>("h_Efficiency_ThetaXZ_ThetaYZ_DifferenceLeadingAndLeadingPlusSecond",
                      "Muon reco efficiency: difference before and after stitching: #theta_{XZ} "
                      "vs. #theta_{YZ}; #theta_{XZ} [#circ]; #theta_{YZ} [#circ]",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NThetaYZBins,
                      ThetaYZBinMin,
                      ThetaYZBinMax);
    h_Efficiency_ThetaXZ_ThetaYZ_DifferenceLeadingAndLeadingPlusSecond->SetOption("colz");

    //Failed Criteria
    h_NoRecoTrackAtAll_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_NoRecoTrackAtAll_ThetaXZ_SinThetaYZ",
                      "# events with no reco track at all: #theta_{XZ} vs. sin(#theta_{YZ}); "
                      "#theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_NoRecoTrackAtAll_ThetaXZ_SinThetaYZ->SetOption("colz");

    h_NoMuonTrack_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_NoMuonTrack_ThetaXZ_SinThetaYZ",
                      "# events with no muon track: #theta_{XZ} vs. sin(#theta_{YZ}); #theta_{XZ} "
                      "[#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_NoMuonTrack_ThetaXZ_SinThetaYZ->SetOption("colz");

    h_TrackTooShort_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_TrackTooShort_ThetaXZ_SinThetaYZ",
                      "# events with L_{reco}/L_{truth} < 75%: #theta_{XZ} vs. sin(#theta_{YZ}); "
                      "#theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_TrackTooShort_ThetaXZ_SinThetaYZ->SetOption("colz");

    h_TrackTooLong_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_TrackTooLong_ThetaXZ_SinThetaYZ",
                      "#events with L_{reco}/L_{truth} > 125%: #theta_{XZ} vs. sin(#theta_{YZ}); "
                      "#theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_TrackTooLong_ThetaXZ_SinThetaYZ->SetOption("colz");

    h_Completeness_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_Completeness_ThetaXZ_SinThetaYZ",
                      "# events with Completeness < 50%: #theta_{XZ} vs. sin(#theta_{YZ}); "
                      "#theta_{XZ} [#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_Completeness_ThetaXZ_SinThetaYZ->SetOption("colz");

    h_Purity_ThetaXZ_SinThetaYZ =
      tfs->make<TH2D>("h_Purity_ThetaXZ_SinThetaYZ",
                      "# events with Purity < 50%: #theta_{XZ} vs. sin(#theta_{YZ}); #theta_{XZ} "
                      "[#circ]; sin(#theta_{YZ})",
                      NThetaXZBins,
                      ThetaXZBinMin,
                      ThetaXZBinMax,
                      NSinThetaYZBins,
                      SinThetaYZBinMin,
                      SinThetaYZBinMax);
    h_Purity_ThetaXZ_SinThetaYZ->SetOption("colz");

    //Criteria vs. NRecoTrack
    h_Criteria_NRecoTrack =
      tfs->make<TH2D>("h_Criteria_NRecoTrack",
                      "Ratio: criteria vs. # reco tracks; Criteria; # reco tracks",
                      NCriteriaBins,
                      CriteriaBinMin,
                      CriteriaBinMax,
                      NRecoTracksBins,
                      RecoTracksBinMin,
                      RecoTracksBinMax);
    h_Criteria_NRecoTrack->SetOption("colz");

    h_Criteria_NRecoTrack_num =
      tfs->make<TH2D>("h_Criteria_NRecoTrack_num",
                      "# events: criteria vs. # reco tracks; Criteria; # reco tracks",
                      NCriteriaBins,
                      CriteriaBinMin,
                      CriteriaBinMax,
                      NRecoTracksBins,
                      RecoTracksBinMin,
                      RecoTracksBinMax);
    h_Criteria_NRecoTrack_num->SetOption("colz");

    h_Criteria_NRecoTrack_den = tfs->make<TH2D>("h_Criteria_NRecoTrack_den",
                                                "Divider histogram; Criteria; # reco tracks",
                                                NCriteriaBins,
                                                CriteriaBinMin,
                                                CriteriaBinMax,
                                                NRecoTracksBins,
                                                RecoTracksBinMin,
                                                RecoTracksBinMax);
    h_Criteria_NRecoTrack_den->SetOption("colz");

    //Criteria vs. NMuonTrack
    h_Criteria_NMuonTrack =
      tfs->make<TH2D>("h_Criteria_NMuonTrack",
                      "Ratio: criteria vs. # muon tracks; Criteria; # muon tracks",
                      NCriteriaBins,
                      CriteriaBinMin,
                      CriteriaBinMax,
                      NRecoTracksBins,
                      RecoTracksBinMin,
                      RecoTracksBinMax);
    h_Criteria_NMuonTrack->SetOption("colz");

    h_Criteria_NMuonTrack_num =
      tfs->make<TH2D>("h_Criteria_NMuonTrack_num",
                      "# events: criteria vs. # muon tracks; Criteria; # muon tracks",
                      NCriteriaBins,
                      CriteriaBinMin,
                      CriteriaBinMax,
                      NRecoTracksBins,
                      RecoTracksBinMin,
                      RecoTracksBinMax);
    h_Criteria_NMuonTrack_num->SetOption("colz");

    h_Criteria_NMuonTrack_den = tfs->make<TH2D>("h_Criteria_NMuonTrack_den",
                                                "Divider histogram; Criteria; # muon tracks",
                                                NCriteriaBins,
                                                CriteriaBinMin,
                                                CriteriaBinMax,
                                                NRecoTracksBins,
                                                RecoTracksBinMin,
                                                RecoTracksBinMax);
    h_Criteria_NMuonTrack_den->SetOption("colz");

    //NoMuonTrack: Max length of no muon track vs. PDG code of that track (MC truth)
    h_NoMuonTrack_MaxTrackLength_PDGCode = tfs->make<TH2D>(
      "h_NoMuonTrack_MaxTrackLength_PDGCode",
      "Events with no muon track: L_{reco, max} vs. PDG Code; L_{reco} [cm]; PDG Code",
      100,
      0.,
      200.,
      200,
      -100.,
      100.);
    h_NoMuonTrack_MaxTrackLength_PDGCode->SetOption("colz");

    //Stitching variables: all events
    h_MuonTrackStitching_TrackRes_Completeness =
      tfs->make<TH2D>("h_MuonTrackStitching_TrackRes_Completeness",
                      "All events: L_{reco}/L_{truth} (leading) vs. Completeness (leading); "
                      "L_{reco}/L_{truth} (leading); Completeness (leading)",
                      150,
                      0.,
                      1.5,
                      150,
                      0.,
                      1.5);
    h_MuonTrackStitching_TrackRes_Completeness->SetOption("colz");

    h_MuonTrackStitching_TrackResLeading_TrackResSecond =
      tfs->make<TH2D>("h_MuonTrackStitching_TrackResLeading_TrackResSecond",
                      "All events: L_{reco}/L_{truth} (leading) vs. L_{reco}/L_{truth} (second); "
                      "L_{reco}/L_{truth} (leading); L_{reco}/L_{truth} (second)",
                      150,
                      0.,
                      1.5,
                      150,
                      0.,
                      1.5);
    h_MuonTrackStitching_TrackResLeading_TrackResSecond->SetOption("colz");

    h_MuonTrackStitching_Distance_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_Distance_Angle",
                      "All events: distance vs. angle b/w leading and second muon track; Distance "
                      "[cm]; Angle [#circ]",
                      100,
                      0.,
                      100.,
                      100,
                      0.,
                      180.);
    h_MuonTrackStitching_Distance_Angle->SetOption("colz");

    h_MuonTrackStitching_TrackResSecondMuon_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_TrackResSecondMuon_Angle",
                      "All events: L_{reco}/L_{truth} (second) vs. angle; L_{reco}/L_{truth} "
                      "(second); Angle [#circ]",
                      150,
                      0.,
                      1.5,
                      180,
                      0,
                      180.);
    h_MuonTrackStitching_TrackResSecondMuon_Angle->SetOption("colz");

    h_MuonTrackStitching_CompletenessSecondMuon_Angle = tfs->make<TH2D>(
      "h_MuonTrackStitching_CompletenessSecondMuon_Angle",
      "All events: Completeness (second) vs. angle; Completeness (second); Angle [#circ]",
      120,
      0.,
      1.2,
      180,
      0,
      180.);
    h_MuonTrackStitching_CompletenessSecondMuon_Angle->SetOption("colz");

    h_MuonTrackStitching_CriteriaTwoTracks_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_CriteriaTwoTracks_Angle",
                      "All events: CriteriaTwoTracks vs. angle; Criteria; Angle [#circ]",
                      7,
                      0.75,
                      4.25,
                      180,
                      0,
                      180.);
    h_MuonTrackStitching_CriteriaTwoTracks_Angle->SetOption("colz");

    //Stitching variables: bad events
    h_MuonTrackStitching_FailedCriteria_TrackRes_Completeness =
      tfs->make<TH2D>("h_MuonTrackStitching_FailedCriteria_TrackRes_Completeness",
                      "Bad events: L_{reco}/L_{truth}  (leading) vs. Completeness (leading); "
                      "L_{reco}/L_{truth} (leading); Completeness (leading)",
                      150,
                      0.,
                      1.5,
                      150,
                      0.,
                      1.5);
    h_MuonTrackStitching_FailedCriteria_TrackRes_Completeness->SetOption("colz");

    h_MuonTrackStitching_FailedCriteria_Distance_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_FailedCriteria_Distance_Angle",
                      "Bad events: distance vs. angle b/w leading and second muon track; Distance "
                      "[cm]; Angle [#circ]",
                      100,
                      0.,
                      100.,
                      100,
                      0.,
                      180.);
    h_MuonTrackStitching_FailedCriteria_Distance_Angle->SetOption("colz");

    h_MuonTrackStitching_FailedCriteria_TrackResLeading_TrackResSecond =
      tfs->make<TH2D>("h_MuonTrackStitching_FailedCriteria_TrackResLeading_TrackResSecond",
                      "Bad events: L_{reco}/L_{truth} (leading) vs. L_{reco}/L_{truth} (second); "
                      "L_{reco}/L_{truth} (leading); L_{reco}/L_{truth} (second)",
                      150,
                      0.,
                      1.5,
                      150,
                      0.,
                      1.5);
    h_MuonTrackStitching_FailedCriteria_TrackResLeading_TrackResSecond->SetOption("colz");

    h_MuonTrackStitching_FailedCriteria_CompletenessLeading_CompletenessSecond =
      tfs->make<TH2D>("*h_MuonTrackStitching_FailedCriteria_CompletenessLeading_CompletenessSecond",
                      "Bad events: Completeness (leading) vs. Completeness (second); Completeness "
                      "(leading); Completeness (second)",
                      150,
                      0.,
                      1.5,
                      150,
                      0.,
                      1.5);
    h_MuonTrackStitching_FailedCriteria_CompletenessLeading_CompletenessSecond->SetOption("colz");

    h_MuonTrackStitching_FailedCriteria_TrackResSecondMuon_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_FailedCriteria_TrackResSecondMuon_Angle",
                      "Bad events: L_{reco}/L_{truth} (second) vs. angle; L_{reco}/L_{truth} "
                      "(second); Angle [#circ]",
                      150,
                      0.,
                      1.5,
                      180,
                      0,
                      180.);
    h_MuonTrackStitching_FailedCriteria_TrackResSecondMuon_Angle->SetOption("colz");

    h_MuonTrackStitching_FailedCriteria_CompletenessSecondMuon_Angle = tfs->make<TH2D>(
      "h_MuonTrackStitching_FailedCriteria_CompletenessSecondMuon_Angle",
      "Bad events: Completeness (second) vs. angle; Completeness (second); Angle [#circ]",
      120,
      0.,
      1.2,
      180,
      0,
      180.);
    h_MuonTrackStitching_FailedCriteria_CompletenessSecondMuon_Angle->SetOption("colz");

    h_MuonTrackStitching_FailedCriteria_CriteriaTwoTracks_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_FailedCriteria_CriteriaTwoTracks_Angle",
                      "Bad events: CriteriaTwoTracks vs. angle; Criteria; Angle [#circ]",
                      7,
                      0.75,
                      4.25,
                      180,
                      0,
                      180.);
    h_MuonTrackStitching_FailedCriteria_CriteriaTwoTracks_Angle->SetOption("colz");

    //Stitching variables: good events
    h_MuonTrackStitching_MatchedCriteria_TrackRes_Completeness =
      tfs->make<TH2D>("h_MuonTrackStitching_MatchedCriteria_TrackRes_Completeness",
                      "Good events: L_{reco}/L_{truth} (leading) vs. Completeness (leading); "
                      "L_{reco}/L_{truth} (leading); Completeness (leading)",
                      150,
                      0.,
                      1.5,
                      150,
                      0.,
                      1.5);
    h_MuonTrackStitching_MatchedCriteria_TrackRes_Completeness->SetOption("colz");

    h_MuonTrackStitching_MatchedCriteria_Distance_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_MatchedCriteria_Distance_Angle",
                      "Good events: distance vs. angle b/w leading and second muon track; Distance "
                      "[cm]; Angle [#circ]",
                      100,
                      0.,
                      100.,
                      100,
                      0.,
                      180.);
    h_MuonTrackStitching_MatchedCriteria_Distance_Angle->SetOption("colz");

    h_MuonTrackStitching_MatchedCriteria_TrackResLeading_TrackResSecond =
      tfs->make<TH2D>("h_MuonTrackStitching_MatchedCriteria_TrackResLeading_TrackResSecond",
                      "Good events: L_{reco}/L_{truth} (leading) vs. L_{reco}/L_{truth} (second); "
                      "L_{reco}/L_{truth} (leading); L_{reco}/L_{truth} (second)",
                      150,
                      0.,
                      1.5,
                      150,
                      0.,
                      1.5);
    h_MuonTrackStitching_MatchedCriteria_TrackResLeading_TrackResSecond->SetOption("colz");

    h_MuonTrackStitching_MatchedCriteria_CriteriaTwoTracks_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_MatchedCriteria_CriteriaTwoTracks_Angle",
                      "Good events: CriteriaTwoTracks vs. angle b/w leading and second muon track; "
                      "Criteria; Angle [#circ]",
                      7,
                      0.75,
                      4.25,
                      100,
                      0.,
                      180.);
    h_MuonTrackStitching_MatchedCriteria_CriteriaTwoTracks_Angle->SetOption("colz");

    //Stitching variables: bad events but leading plus second ok
    h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackRes_Completeness =
      tfs->make<TH2D>(
        "h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackRes_Completeness",
        "Bad events but leading + second is good: L_{reco}/L_{truth}  (leading) vs. Completeness "
        "(leading); L_{reco}/L_{truth} (leading); Completeness (leading)",
        150,
        0.,
        1.5,
        150,
        0.,
        1.5);
    h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackRes_Completeness->SetOption(
      "colz");

    h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_Distance_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_Distance_Angle",
                      "Bad events but leading + second is good: distance vs. angle b/w leading and "
                      "second muon track; Distance [cm]; Angle [#circ]",
                      100,
                      0.,
                      100.,
                      100,
                      0.,
                      180.);
    h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_Distance_Angle->SetOption("colz");

    h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackResLeading_TrackResSecond =
      tfs->make<TH2D>(
        "h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackResLeading_"
        "TrackResSecond",
        "Bad events but leading + second is good: L_{reco}/L_{truth} (leading) vs. "
        "L_{reco}/L_{truth} (second); L_{reco}/L_{truth} (leading); L_{reco}/L_{truth} (second)",
        150,
        0.,
        1.5,
        150,
        0.,
        1.5);
    h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackResLeading_TrackResSecond
      ->SetOption("colz");

    //Stitching variables: bad events and leading plus second not ok
    h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackRes_Completeness =
      tfs->make<TH2D>(
        "h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackRes_Completeness",
        "Bad events and leading + second is bad: L_{reco}/L_{truth}  (leading) vs. Completeness "
        "(leading); L_{reco}/L_{truth} (leading); Completeness (leading)",
        150,
        0.,
        1.5,
        150,
        0.,
        1.5);
    h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackRes_Completeness->SetOption(
      "colz");

    h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_Distance_Angle =
      tfs->make<TH2D>("h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_Distance_Angle",
                      "Bad events and leading + second is bad: distance vs. angle b/w leading and "
                      "second muon track; Distance [cm]; Angle [#circ]",
                      100,
                      0.,
                      100.,
                      100,
                      0.,
                      180.);
    h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_Distance_Angle->SetOption("colz");

    h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackResLeading_TrackResSecond =
      tfs->make<TH2D>(
        "h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackResLeading_TrackResSecond",
        "Bad events and leading + second is bad: L_{reco}/L_{truth} (leading) vs. "
        "L_{reco}/L_{truth} (second); L_{reco}/L_{truth} (leading); L_{reco}/L_{truth} (second)",
        150,
        0.,
        1.5,
        150,
        0.,
        1.5);
    h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackResLeading_TrackResSecond
      ->SetOption("colz");
  }
  //========================================================================
  void MuonTrackingEff::endJob() { doEfficiencies(); }
  //========================================================================
  void MuonTrackingEff::beginRun(const art::Run& /*run*/)
  {
    mf::LogInfo("MuonTrackingEff") << "begin run..." << std::endl;
  }
  //========================================================================
  void MuonTrackingEff::analyze(const art::Event& event)
  {
    if (event.isRealData()) return;

    bool isFiducial = false;
    processEff(event, isFiducial);
  }
  //========================================================================
  void MuonTrackingEff::processEff(const art::Event& event, bool& isFiducial)
  {

    EventCounter++;
    simb::MCParticle* MCTruthMuonParticle = NULL;

    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    simb::MCParticle* particle = 0;

    for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar) {
      particle = ipar->second;

      if (particle->Mother() ==
          0) { //0=particle has no mother particle, 1=particle has a mother particle
        const TLorentzVector& positionStart = particle->Position(0);
        positionStart.GetXYZT(
          MCTruthMuonVertex); //MCTruthMuonVertex[0-2]=truth vertex, MCTruthMuonVertex[3]=t=0
      }

      if (particle->PdgCode() == fMuonPDGCode) { // Particle cannon muon
        MCTruthMuonParticle = particle;
        MCTruthMuonID = particle->TrackId();
        MCTruthMuonMomentum =
          sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
               pow(particle->Momentum().Pz(), 2));

        if (particle->Momentum().Pz() >= 0 && particle->Momentum().Px() >= 0) {
          MCTruthMuonThetaXZ =
            (180.0 / 3.14159) * atan(particle->Momentum().Px() / particle->Momentum().Pz());
        }
        else if (particle->Momentum().Pz() < 0 && particle->Momentum().Px() >= 0) {
          MCTruthMuonThetaXZ =
            180.0 + (180.0 / 3.14159) * atan(particle->Momentum().Px() / particle->Momentum().Pz());
        }
        else if (particle->Momentum().Pz() < 0 && particle->Momentum().Px() < 0) {
          MCTruthMuonThetaXZ =
            180.0 + (180.0 / 3.14159) * atan(particle->Momentum().Px() / particle->Momentum().Pz());
        }
        else if (particle->Momentum().Pz() >= 0 && particle->Momentum().Px() < 0) {
          MCTruthMuonThetaXZ =
            360.0 + (180.0 / 3.14159) * atan(particle->Momentum().Px() / particle->Momentum().Pz());
        }

        MCTruthMuonThetaYZ =
          (180.0 / 3.14159) * asin(particle->Momentum().Py() / MCTruthMuonMomentum);
      }
    }
    double MCTruthLengthMuon = truthLength(MCTruthMuonParticle);
    h_TruthLength->Fill(MCTruthLengthMuon);

    //===================================================================
    //Saving denominator histograms
    //===================================================================
    isFiducial = insideFV(MCTruthMuonVertex);
    if (!isFiducial) return;

    //save events for Nucleon decay and particle cannon
    if (MCTruthMuonParticle) {
      h_ThetaXZ_den->Fill(MCTruthMuonThetaXZ);
      h_ThetaYZ_den->Fill(MCTruthMuonThetaYZ);
      h_SinThetaYZ_den->Fill(sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
      h_ThetaXZ_ThetaYZ_den->Fill(MCTruthMuonThetaXZ, MCTruthMuonThetaYZ);
      h_ThetaXZ_SinThetaYZ_den->Fill(MCTruthMuonThetaXZ,
                                     sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
      CountMCTruthMuon++;
    }

    //========================================================================
    // Reco stuff, once we have selected a MC Particle let's find out if there is a track associated
    //========================================================================

    int NMuonTracks = 0;

    art::Handle<std::vector<recob::Track>> TrackListHandle;
    if (!event.getByLabel(fTrackModuleLabel, TrackListHandle)) return;
    std::vector<art::Ptr<recob::Track>> TrackList;
    art::fill_ptr_vector(TrackList, TrackListHandle);
    int NRecoTracks = TrackList.size();
    art::FindManyP<recob::Hit> track_hits(TrackListHandle, event, fTrackModuleLabel);
    if (NRecoTracks == 0) {
      MF_LOG_DEBUG("MuonTrackingEff") << "There are no reco tracks... bye";
      std::cout << "There are no reco tracks! MCTruthMuonThetaXZ: " << std::endl;

      Criteria = 0.;
      h_FailedReconstruction_ThetaXZ_ThetaYZ->Fill(MCTruthMuonThetaXZ, MCTruthMuonThetaYZ);
      h_FailedReconstruction_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                                      sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
      h_NoRecoTrackAtAll_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                                  sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
      h_Criteria_NRecoTrack_num->Fill(Criteria, static_cast<double>(NRecoTracks));
      h_Criteria_NMuonTrack_num->Fill(Criteria, static_cast<double>(NMuonTracks));
      CountNoRecoTracks++;
      return;
    }
    MF_LOG_DEBUG("MuonTrackingEff") << "Found this many reco tracks " << NRecoTracks;

    //std::cout << "NRecoTracks: " << NRecoTracks << std::endl;

    double PurityLeadingMuon = 0.;
    double CompletenessLeadingMuon = 0.;
    double RecoLengthLeadingMuon = 0.;
    art::Ptr<recob::Track> TrackLeadingMuon;

    double RecoLengthSecondMuon = 0.;
    double CompletenessSecondMuon = 0.;
    double PuritySecondMuon = 0.;
    art::Ptr<recob::Track> TrackSecondMuon;

    double TrackLengthMuonSum = 0.;
    double tmpTotalRecoEnergy = 0.;

    double MaxLengthNoRecoMuon = 0;
    int PDGCodeMaxLengthNoRecoMuon = 0;

    const simb::MCParticle* RecoMuonParticle = NULL;

    std::vector<art::Ptr<recob::Hit>> const& tmp_TrackHits = track_hits.at(0);
    std::vector<recob::Hit> const& AllHits = tmp_TrackHits[0].parentAs<std::vector>();

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);

    // Loop over reco tracks
    for (int i = 0; i < NRecoTracks; i++) {
      art::Ptr<recob::Track> track = TrackList[i];
      std::vector<art::Ptr<recob::Hit>> TrackHits = track_hits.at(i);
      double tmpPurity = 0.;
      double tmpCompleteness = 0.;
      const simb::MCParticle* particle;

      truthMatcher(
        clockData, AllHits, TrackHits, particle, tmpPurity, tmpCompleteness, tmpTotalRecoEnergy);

      if (!particle) {
        std::cout << "ERROR: Truth matcher didn't find a particle!" << std::endl;
        continue;
      }

      if (track->Length() > MaxLengthNoRecoMuon && particle->PdgCode() != fMuonPDGCode &&
          particle->TrackId() != MCTruthMuonID) {
        MaxLengthNoRecoMuon = track->Length();
        PDGCodeMaxLengthNoRecoMuon = particle->PdgCode();
      }
      //some muon tracks have Completeness=0 and Purity=0 but are still considered as muon tracks in function truthmatcher. Getting rid of these tracks by asking tmpCompleteness > 0 && tmpPurity > 0
      if ((particle->PdgCode() == fMuonPDGCode) && (particle->TrackId() == MCTruthMuonID) &&
          tmpCompleteness > 0 && tmpPurity > 0) {

        NMuonTracks++;
        TrackLengthMuonSum += track->Length();

        if (NMuonTracks == 1) {
          CompletenessLeadingMuon = tmpCompleteness;
          PurityLeadingMuon = tmpPurity;
          RecoLengthLeadingMuon = track->Length();
          TrackLeadingMuon = track;

          RecoMuonParticle = particle;
        }

        if (NMuonTracks >= 2 && tmpCompleteness > CompletenessLeadingMuon) {

          CompletenessSecondMuon = CompletenessLeadingMuon;
          PuritySecondMuon = PurityLeadingMuon;
          RecoLengthSecondMuon = RecoLengthLeadingMuon;
          TrackSecondMuon = TrackLeadingMuon;

          CompletenessLeadingMuon = tmpCompleteness;
          PurityLeadingMuon = tmpPurity;
          RecoLengthLeadingMuon = track->Length();
          TrackLeadingMuon = track;

          RecoMuonParticle = particle;
        }

        else if (NMuonTracks >= 2 && tmpCompleteness < CompletenessLeadingMuon &&
                 tmpCompleteness > CompletenessSecondMuon) {
          CompletenessSecondMuon = tmpCompleteness;
          PuritySecondMuon = tmpPurity;
          RecoLengthSecondMuon = track->Length();
          TrackSecondMuon = track;
        }
      }
    }

    h_TotalRecoEnergy->Fill(tmpTotalRecoEnergy);

    //Muon events
    //=========================================================

    if (RecoMuonParticle && MCTruthMuonParticle) {
      CountRecoMuon++;
      h_Purity->Fill(PurityLeadingMuon);
      h_Completeness->Fill(CompletenessLeadingMuon);
      h_TrackRes->Fill(RecoLengthLeadingMuon / MCTruthLengthMuon);

      std::cout << "TrackLeadingMuon->Vertex().X(): " << TrackLeadingMuon->Vertex().X()
                << std::endl;
      std::cout << "MCTruthMuonParticle->Vz(): " << MCTruthMuonParticle->Vz() << std::endl;
      std::cout << " " << std::endl;

      double DistanceBetweenTruthAndRecoTrack;
      double AngleBetweenTruthAndRecoTrack;
      FuncDistanceAndAngleBetweenTruthAndRecoTrack(RecoMuonParticle,
                                                   TrackLeadingMuon,
                                                   DistanceBetweenTruthAndRecoTrack,
                                                   AngleBetweenTruthAndRecoTrack);

      h_VertexRes->Fill(DistanceBetweenTruthAndRecoTrack);
      h_DirectionRes->Fill(AngleBetweenTruthAndRecoTrack);

      h_MuonTrackStitching_TrackRes_Completeness->Fill(RecoLengthLeadingMuon / MCTruthLengthMuon,
                                                       CompletenessLeadingMuon);
      if (NMuonTracks >= 2) {
        double DistanceBetweenTracks;
        double AngleBetweenTracks;
        double CriteriaTwoTracks;

        FuncDistanceAndAngleBetweenTracks(TrackSecondMuon,
                                          TrackLeadingMuon,
                                          DistanceBetweenTracks,
                                          AngleBetweenTracks,
                                          CriteriaTwoTracks);

        h_MuonTrackStitching_TrackResLeading_TrackResSecond->Fill(
          RecoLengthLeadingMuon / MCTruthLengthMuon, RecoLengthSecondMuon / MCTruthLengthMuon);
        h_MuonTrackStitching_Distance_Angle->Fill(DistanceBetweenTracks, AngleBetweenTracks);
        h_MuonTrackStitching_TrackResSecondMuon_Angle->Fill(
          RecoLengthSecondMuon / MCTruthLengthMuon, AngleBetweenTracks);
        h_MuonTrackStitching_CompletenessSecondMuon_Angle->Fill(CompletenessSecondMuon,
                                                                AngleBetweenTracks);
        h_MuonTrackStitching_CriteriaTwoTracks_Angle->Fill(CriteriaTwoTracks, AngleBetweenTracks);
      }

      //Completeness
      if (CompletenessLeadingMuon < 0.5) {
        CountCompleteness++;
        Criteria = 2.;
        h_Completeness_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                                sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
        h_Criteria_NRecoTrack_num->Fill(Criteria, static_cast<double>(NRecoTracks));
        h_Criteria_NMuonTrack_num->Fill(Criteria, static_cast<double>(NMuonTracks));
      }
      //Purity
      if (PurityLeadingMuon < 0.5) {
        CountPurity++;
        Criteria = 3.;
        h_Purity_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                          sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
        h_Criteria_NRecoTrack_num->Fill(Criteria, static_cast<double>(NRecoTracks));
        h_Criteria_NMuonTrack_num->Fill(Criteria, static_cast<double>(NMuonTracks));
      }
      //Track too short
      if (RecoLengthLeadingMuon / MCTruthLengthMuon < 0.75) {
        CountTrackLengthTooShort++;
        Criteria = 4.;
        h_TrackTooShort_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                                 sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
        h_Criteria_NRecoTrack_num->Fill(Criteria, static_cast<double>(NRecoTracks));
        h_Criteria_NMuonTrack_num->Fill(Criteria, static_cast<double>(NMuonTracks));
      }
      //Track too long
      if (RecoLengthLeadingMuon / MCTruthLengthMuon > 1.25) {
        CountTrackLengthTooLong++;
        Criteria = 5.;
        h_TrackTooLong_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                                sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
        h_Criteria_NRecoTrack_num->Fill(Criteria, static_cast<double>(NRecoTracks));
        h_Criteria_NMuonTrack_num->Fill(Criteria, static_cast<double>(NMuonTracks));

        /*std::cout << "Track too long!" << std::endl;
                        art::ServiceHandle<cheat::BackTracker const> bt2;
                        const sim::ParticleList& plist2 = bt2->ParticleList();
                        simb::MCParticle *particle2=0;

                        std::cout << "EventCounter: " << EventCounter << std::endl;
                        for( sim::ParticleList::const_iterator ipar2 = plist2.begin(); ipar2!=plist2.end(); ++ipar2)
                        {
                         particle2 = ipar2->second;
                        std::cout << "particle2->TrackId(): " << particle2->TrackId() << std::endl;
                        std::cout << "particle2->PdgCode(): " << particle2->PdgCode() << std::endl;
                        std::cout << "particle2->Momentum().P(): " << particle2->Momentum().P() << std::endl;
                        std::cout << "tuthLength(particle2): " << truthLength(particle2) << std::endl;
                        std::cout << "particle2->Position(): (x,y,z): " << particle2->Vx() << "\t" << particle2->Vy() << "\t" << particle2->Vz() << std::endl;
                        std::cout << "particle2->Momentum(): (Px,Py,Pz): " << particle2->Momentum().Px() << "\t" << particle2->Momentum().Py() << "\t" << particle2->Momentum().Pz() << std::endl;
                        std::cout << "particle2->Position().T(): " << particle2->Position().T() << std::endl;
                        std::cout << "" << std::endl;
                        }*/
      }
      //Reco failed at least one of the above criteria
      if (CompletenessLeadingMuon < 0.5 || PurityLeadingMuon < 0.5 ||
          RecoLengthLeadingMuon / MCTruthLengthMuon < 0.75 ||
          RecoLengthLeadingMuon / MCTruthLengthMuon > 1.25) {
        CountBadLeadingMuonTrack++;
        h_FailedReconstruction_ThetaXZ_ThetaYZ->Fill(MCTruthMuonThetaXZ, MCTruthMuonThetaYZ);
        h_FailedReconstruction_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                                        sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
        h_MuonTrackStitching_FailedCriteria_TrackRes_Completeness->Fill(
          RecoLengthLeadingMuon / MCTruthLengthMuon, CompletenessLeadingMuon);

        if (NMuonTracks == 1) { BadEvents1MuonTrack++; }
        if (NMuonTracks == 2) { BadEvents2MuonTrack++; }
        if (NMuonTracks == 3) { BadEvents3MuonTrack++; }
        if (NMuonTracks == 4) { BadEvents4OrMoreMuonTrack++; }

        if (NMuonTracks >= 2) {
          double AngleBetweenTracks;
          double DistanceBetweenTracks;
          double CriteriaTwoTracks;
          FuncDistanceAndAngleBetweenTracks(TrackSecondMuon,
                                            TrackLeadingMuon,
                                            DistanceBetweenTracks,
                                            AngleBetweenTracks,
                                            CriteriaTwoTracks);

          if (AngleBetweenTracks > 160.) {
            std::cout << "EventCounter: " << EventCounter << std::endl;
            std::cout << "Angle b/w tracks: " << AngleBetweenTracks << std::endl;
            std::cout << "CriteriaTwoTracks: " << CriteriaTwoTracks << std::endl;
            std::cout << "CompletenessLeadingMuon: " << CompletenessLeadingMuon << std::endl;
            std::cout << "CompletenessSecondMuon: " << CompletenessSecondMuon << std::endl;
            std::cout << "PurityLeadingMuon: " << PurityLeadingMuon << std::endl;
            std::cout << "PuritySecondMuon: " << PuritySecondMuon << std::endl;
            std::cout << "TrackLeadingMuon: " << RecoLengthLeadingMuon / MCTruthLengthMuon
                      << std::endl;
            std::cout << "TrackResSecondMuon: " << RecoLengthSecondMuon / MCTruthLengthMuon
                      << std::endl;
            std::cout << "TrackID leading: " << TrackLeadingMuon->ID() << std::endl;
            std::cout << "TrackID second: " << TrackSecondMuon->ID() << std::endl;
          }

          h_MuonTrackStitching_FailedCriteria_Distance_Angle->Fill(DistanceBetweenTracks,
                                                                   AngleBetweenTracks);
          h_MuonTrackStitching_FailedCriteria_TrackResLeading_TrackResSecond->Fill(
            RecoLengthLeadingMuon / MCTruthLengthMuon, RecoLengthSecondMuon / MCTruthLengthMuon);
          h_MuonTrackStitching_FailedCriteria_CompletenessLeading_CompletenessSecond->Fill(
            CompletenessLeadingMuon, CompletenessSecondMuon);
          h_MuonTrackStitching_FailedCriteria_TrackResSecondMuon_Angle->Fill(
            RecoLengthSecondMuon / MCTruthLengthMuon, AngleBetweenTracks);
          h_MuonTrackStitching_FailedCriteria_CompletenessSecondMuon_Angle->Fill(
            CompletenessSecondMuon, AngleBetweenTracks);
          h_MuonTrackStitching_FailedCriteria_CriteriaTwoTracks_Angle->Fill(CriteriaTwoTracks,
                                                                            AngleBetweenTracks);
          if ((CompletenessLeadingMuon + CompletenessSecondMuon) >= 0.5 &&
              PurityLeadingMuon >= 0.5 && PuritySecondMuon >= 0.5 &&
              (RecoLengthLeadingMuon + RecoLengthSecondMuon) / MCTruthLengthMuon >= 0.75 &&
              (RecoLengthLeadingMuon + RecoLengthSecondMuon) / MCTruthLengthMuon <= 1.25) {
            CountBadLeadingMuonTrackButLeadingPlusSecondGood++;
            h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk->Fill(MCTruthMuonThetaXZ, MCTruthMuonThetaYZ);
            h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk->Fill(
              MCTruthMuonThetaXZ, sin((3.14159 / 180.) * MCTruthMuonThetaYZ));

            h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackRes_Completeness->Fill(
              RecoLengthLeadingMuon / MCTruthLengthMuon, CompletenessLeadingMuon);
            h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_TrackResLeading_TrackResSecond
              ->Fill(RecoLengthLeadingMuon / MCTruthLengthMuon,
                     RecoLengthSecondMuon / MCTruthLengthMuon);
            h_MuonTrackStitching_FailedCriteriaButLeadingPlusSecondGood_Distance_Angle->Fill(
              DistanceBetweenTracks, AngleBetweenTracks);
          }
          if ((CompletenessLeadingMuon + CompletenessSecondMuon) < 0.5 || PurityLeadingMuon < 0.5 ||
              PuritySecondMuon < 0.5 ||
              (RecoLengthLeadingMuon + RecoLengthSecondMuon) / MCTruthLengthMuon < 0.75 ||
              (RecoLengthLeadingMuon + RecoLengthSecondMuon) / MCTruthLengthMuon > 1.25) {
            CountBadLeadingMuonTrackAndLeadingPlusSecondBad++;

            h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackRes_Completeness->Fill(
              RecoLengthLeadingMuon / MCTruthLengthMuon, CompletenessLeadingMuon);
            h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_TrackResLeading_TrackResSecond
              ->Fill(RecoLengthLeadingMuon / MCTruthLengthMuon,
                     RecoLengthSecondMuon / MCTruthLengthMuon);
            h_MuonTrackStitching_FailedCriteriaAndLeadingPlusSecondBad_Distance_Angle->Fill(
              DistanceBetweenTracks, AngleBetweenTracks);
          }
          if ((CompletenessLeadingMuon + CompletenessSecondMuon) < 0.5) {
            CountBadLeadingMuonTrackAndLeadingPlusSecondBadCompleteness++;
          }
          if (PurityLeadingMuon < 0.5 || PuritySecondMuon < 0.5) {
            CountBadLeadingMuonTrackAndLeadingPlusSecondBadPurity++;
          }
          if ((RecoLengthLeadingMuon + RecoLengthSecondMuon) / MCTruthLengthMuon < 0.75) {
            CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooShort++;
          }
          if ((RecoLengthLeadingMuon + RecoLengthSecondMuon) / MCTruthLengthMuon > 1.25) {
            CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooLong++;
          }
        }
        else if (NMuonTracks == 1) {
          CountBadLeadingMuonTrackAndOnlyOneMuonTrack++;
          if (CompletenessLeadingMuon < 0.5) {
            CountBadLeadingMuonTrackAndOnlyOneMuonTrackCompleteness++;
          }
          if (PurityLeadingMuon < 0.5) { CountBadLeadingMuonTrackAndOnlyOneMuonTrackPurity++; }
          if (RecoLengthLeadingMuon / MCTruthLengthMuon < 0.75) {
            CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooShort++;
          }
          if (RecoLengthLeadingMuon / MCTruthLengthMuon > 1.25) {
            CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooLong++;
          }
        }
      }
      //Reco succesful according to the above criteria
      if (CompletenessLeadingMuon >= 0.5 && PurityLeadingMuon >= 0.5 &&
          RecoLengthLeadingMuon / MCTruthLengthMuon >= 0.75 &&
          RecoLengthLeadingMuon / MCTruthLengthMuon <= 1.25) {
        CountGoodLeadingMuonTrack++;
        Criteria = 6.;
        h_ThetaXZ_num->Fill(MCTruthMuonThetaXZ);
        h_ThetaYZ_num->Fill(MCTruthMuonThetaYZ);
        h_SinThetaYZ_num->Fill(sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
        h_ThetaXZ_ThetaYZ_num->Fill(MCTruthMuonThetaXZ, MCTruthMuonThetaYZ);
        h_ThetaXZ_SinThetaYZ_num->Fill(MCTruthMuonThetaXZ,
                                       sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
        h_Criteria_NRecoTrack_num->Fill(Criteria, static_cast<double>(NRecoTracks));
        h_Criteria_NMuonTrack_num->Fill(Criteria, static_cast<double>(NMuonTracks));

        h_MuonTrackStitching_MatchedCriteria_TrackRes_Completeness->Fill(
          RecoLengthLeadingMuon / MCTruthLengthMuon, CompletenessLeadingMuon);

        if (NMuonTracks == 1) { GoodEvents1MuonTrack++; }
        if (NMuonTracks == 2) { GoodEvents2MuonTrack++; }
        if (NMuonTracks == 3) { GoodEvents3MuonTrack++; }
        if (NMuonTracks == 4) { GoodEvents4OrMoreMuonTrack++; }

        if (NMuonTracks >= 2) {
          h_MuonTrackStitching_MatchedCriteria_TrackResLeading_TrackResSecond->Fill(
            RecoLengthLeadingMuon / MCTruthLengthMuon, RecoLengthSecondMuon / MCTruthLengthMuon);

          double AngleBetweenTracks;
          double DistanceBetweenTracks;
          double CriteriaTwoTracks;
          FuncDistanceAndAngleBetweenTracks(TrackSecondMuon,
                                            TrackLeadingMuon,
                                            DistanceBetweenTracks,
                                            AngleBetweenTracks,
                                            CriteriaTwoTracks);
          h_MuonTrackStitching_MatchedCriteria_Distance_Angle->Fill(DistanceBetweenTracks,
                                                                    AngleBetweenTracks);
          h_MuonTrackStitching_MatchedCriteria_CriteriaTwoTracks_Angle->Fill(CriteriaTwoTracks,
                                                                             AngleBetweenTracks);
        }
      }
    }
    //No muon track
    if (!RecoMuonParticle && MCTruthMuonParticle) {
      CountNoMuonTracks++;
      BadEvents0MuonTrack++;
      Criteria = 1.;
      h_NoMuonTrack_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                             sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
      h_FailedReconstruction_ThetaXZ_ThetaYZ->Fill(MCTruthMuonThetaXZ, MCTruthMuonThetaYZ);
      h_FailedReconstruction_ThetaXZ_SinThetaYZ->Fill(MCTruthMuonThetaXZ,
                                                      sin((3.14159 / 180.) * MCTruthMuonThetaYZ));
      h_Criteria_NRecoTrack_num->Fill(Criteria, static_cast<double>(NRecoTracks));
      h_Criteria_NMuonTrack_num->Fill(Criteria, static_cast<double>(NMuonTracks));
      h_NoMuonTrack_MaxTrackLength_PDGCode->Fill(MaxLengthNoRecoMuon,
                                                 static_cast<double>(PDGCodeMaxLengthNoRecoMuon));
    }
  }
  //========================================================================
  void MuonTrackingEff::truthMatcher(detinfo::DetectorClocksData const& clockData,
                                     std::vector<recob::Hit> const& AllHits,
                                     std::vector<art::Ptr<recob::Hit>> track_hits,
                                     const simb::MCParticle*& MCparticle,
                                     double& Purity,
                                     double& Completeness,
                                     double& TotalRecoEnergy)
  {
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    std::map<int, double> trkID_E; // map that connects TrackID and energy for
                                   // each hit <trackID, energy>
    for (size_t j = 0; j < track_hits.size(); ++j) {
      art::Ptr<recob::Hit> hit = track_hits[j];
      std::vector<sim::TrackIDE> TrackIDs =
        bt_serv->HitToTrackIDEs(clockData,
                                hit); // TrackIDE contains TrackID, energy and energyFrac. A hit can
      // have several TrackIDs (so this hit is associated with multiple
      // MC truth track IDs (EM shower IDs are negative). If a hit ahs
      // multiple trackIDs, "energyFrac" contains the fraction of the
      // energy of for each ID compared to the total energy of the hit.
      // "energy" contains only the energy associated with the specific
      // ID in that case. This requires MC truth info!
      for (size_t k = 0; k < TrackIDs.size(); k++) {
        trkID_E[TrackIDs[k].trackID] +=
          TrackIDs[k].energy; // sum up the energy for each TrackID and store
                              // <TrackID, energy> in "TrkID_E"
      }
    }

    double E_em = 0.0;
    double max_E = -999.0;
    double TotalEnergyTrack = 0.0;
    int TrackID = -999;
    double PartialEnergyTrackID =
      0.0; // amount of energy deposited by the particle that deposited more
           // energy... tomato potato... blabla
    //! if the collection of hits have more than one particle associate save the
    //! particle w/ the highest energy deposition since we are looking for
    //! muons/pions/protons this should be enough
    if (!trkID_E.size()) {
      MCparticle = 0;
      return; // Ghost track???
    }
    for (std::map<int, double>::iterator ii = trkID_E.begin(); ii != trkID_E.end();
         ++ii) {                      // trkID_E contains the trackID (first) and corresponding
                                      // energy (second) for a specific track, summed up over all
                                      // events. here looping over all trekID_E's
      TotalEnergyTrack += ii->second; // and summing up the energy of all hits
                                      // in the track (TotalEnergyTrack)
      if ((ii->second) > max_E) {     // looking for the trakID with the highest energy in the
                                      // track.  this is PartialEnergyTrackID and max_E then
        PartialEnergyTrackID = ii->second;
        max_E = ii->second;
        TrackID = ii->first;                 // saving trackID of the ID with the highest energy
                                             // contribution in the track to assign it to
                                             // MCparticle later
        if (TrackID < 0) E_em += ii->second; // IDs of em shower particles are negative
      }
    }
    MCparticle = pi_serv->TrackIdToParticle_P(TrackID);
    // In the current simulation, we do not save EM Shower daughters in GEANT.
    // But we do save the energy deposition in TrackIDEs. If the energy
    // deposition is from a particle that is the daughter of an EM particle, the
    // negative of the parent track ID is saved in TrackIDE for the daughter
    // particle we don't want to track gammas or any other EM activity
    if (TrackID < 0) { return; }

    // Purity = (PartialEnergyTrackID+E_em)/TotalEnergyTrack;
    Purity = PartialEnergyTrackID / TotalEnergyTrack;

    // completeness
    TotalRecoEnergy = 0;
    for (recob::Hit const& hit : AllHits) {
      std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
      for (size_t l = 0; l < TrackIDs.size(); ++l) { // and over all track IDs of the hits
        if (TrackIDs[l].trackID == TrackID)
          TotalRecoEnergy += TrackIDs[l].energy; // and sum up the energy fraction of all hits
                                                 // that correspond ot the saved trackID
      }
    }
    Completeness = PartialEnergyTrackID / TotalRecoEnergy;
  }

  void MuonTrackingEff::FuncDistanceAndAngleBetweenTracks(art::Ptr<recob::Track> Track1,
                                                          art::Ptr<recob::Track> Track2,
                                                          double& TempDistanceBetweenTracks,
                                                          double& TempAngleBetweenTracks,
                                                          double& TempCriteriaTwoTracks)
  {

    TempDistanceBetweenTracks = sqrt(pow(Track1->End().X() - Track2->Vertex().X(), 2) +
                                     pow(Track1->End().Y() - Track2->Vertex().Y(), 2) +
                                     pow(Track1->End().Z() - Track2->Vertex().Z(), 2));
    TempAngleBetweenTracks = (180.0 / 3.14159) * Track1->EndDirection<TVector3>().Angle(
                                                   Track2->VertexDirection<TVector3>());
    TempCriteriaTwoTracks = 1.;

    if (TempDistanceBetweenTracks > sqrt(pow(Track1->End().X() - Track2->End().X(), 2) +
                                         pow(Track1->End().Y() - Track2->End().Y(), 2) +
                                         pow(Track1->End().Z() - Track2->End().Z(), 2))) {
      TempDistanceBetweenTracks = sqrt(pow(Track1->End().X() - Track2->End().X(), 2) +
                                       pow(Track1->End().Y() - Track2->End().Y(), 2) +
                                       pow(Track1->End().Z() - Track2->End().Z(), 2));
      TempAngleBetweenTracks = 180. - (180.0 / 3.14159) * Track1->EndDirection<TVector3>().Angle(
                                                            Track2->EndDirection<TVector3>());
      TempCriteriaTwoTracks = 2.;
    }

    if (TempDistanceBetweenTracks > sqrt(pow(Track1->Vertex().X() - Track2->End().X(), 2) +
                                         pow(Track1->Vertex().Y() - Track2->End().Y(), 2) +
                                         pow(Track1->Vertex().Z() - Track2->End().Z(), 2))) {
      TempDistanceBetweenTracks = sqrt(pow(Track1->Vertex().X() - Track2->End().X(), 2) +
                                       pow(Track1->Vertex().Y() - Track2->End().Y(), 2) +
                                       pow(Track1->Vertex().Z() - Track2->End().Z(), 2));
      TempAngleBetweenTracks = (180.0 / 3.14159) * Track1->VertexDirection<TVector3>().Angle(
                                                     Track2->EndDirection<TVector3>());
      TempCriteriaTwoTracks = 3.;
    }

    if (TempDistanceBetweenTracks > sqrt(pow(Track1->Vertex().X() - Track2->Vertex().X(), 2) +
                                         pow(Track1->Vertex().Y() - Track2->Vertex().Y(), 2) +
                                         pow(Track1->Vertex().Z() - Track2->Vertex().Z(), 2))) {
      TempDistanceBetweenTracks = sqrt(pow(Track1->Vertex().X() - Track2->Vertex().X(), 2) +
                                       pow(Track1->Vertex().Y() - Track2->Vertex().Y(), 2) +
                                       pow(Track1->Vertex().Z() - Track2->Vertex().Z(), 2));
      TempAngleBetweenTracks = 180. - (180.0 / 3.14159) * Track1->VertexDirection<TVector3>().Angle(
                                                            Track2->VertexDirection<TVector3>());
      TempCriteriaTwoTracks = 4.;
    }
  }

  void MuonTrackingEff::FuncDistanceAndAngleBetweenTruthAndRecoTrack(
    const simb::MCParticle*& MCparticle,
    art::Ptr<recob::Track> Track,
    double& TempDistanceBetweenTruthAndRecoTrack,
    double& TempAngleBeetweenTruthAndRecoTrack)
  {
    TempDistanceBetweenTruthAndRecoTrack = sqrt(pow(Track->Vertex().X() - MCparticle->Vx(), 2) +
                                                pow(Track->Vertex().Y() - MCparticle->Vy(), 2) +
                                                pow(Track->Vertex().Z() - MCparticle->Vz(), 2));

    TempAngleBeetweenTruthAndRecoTrack = 0;
  }

  //========================================================================
  double MuonTrackingEff::truthLength(const simb::MCParticle* MCparticle)
  {
    //calculate the truth length considering only the part that is inside the TPC
    //Base on a piece of code from dune/TrackingAna/TrackingEfficiency_module.cc

    if (!MCparticle) return -999.0;
    int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
    std::vector<double> TPCLengthHits(numberTrajectoryPoints, 0);
    int FirstHit = 0, LastHit = 0;
    double TPCLength = 0.0;
    bool BeenInVolume = false;

    for (unsigned int MCHit = 0; MCHit < TPCLengthHits.size(); ++MCHit) {
      if (MCHit != 0)
        TPCLengthHits[MCHit] = sqrt(pow((MCparticle->Vx(MCHit - 1) - MCparticle->Vx(MCHit)), 2) +
                                    pow((MCparticle->Vy(MCHit - 1) - MCparticle->Vy(MCHit)), 2) +
                                    pow((MCparticle->Vz(MCHit - 1) - MCparticle->Vz(MCHit)), 2));
      geo::TPCID tpcid =
        geom->FindTPCAtPosition(geo::vect::toPoint(MCparticle->Position(MCHit).Vect()));
      if (tpcid.isValid) {
        LastHit = MCHit;
        if (!BeenInVolume) {
          BeenInVolume = true;
          FirstHit = MCHit;
        }
      }
    }
    for (int Hit = FirstHit + 1; Hit <= LastHit; ++Hit)
      TPCLength += TPCLengthHits[Hit];
    return TPCLength;
  }
  //========================================================================
  bool MuonTrackingEff::insideFV(double vertex[4])
  {
    double x = vertex[0];
    double y = vertex[1];
    double z = vertex[2];

    return x > fFidVolXmin && x < fFidVolXmax && y > fFidVolYmin && y < fFidVolYmax &&
           z > fFidVolZmin && z < fFidVolZmax;
  }
  //========================================================================
  void MuonTrackingEff::doEfficiencies()
  {
    std::cout << std::endl;

    std::cout << "EventCounter: "
              << "\t" << EventCounter << std::endl;

    std::cout << "CountMCTruthMuon: "
              << "\t" << CountMCTruthMuon << " = "
              << 100 * static_cast<double>(CountMCTruthMuon) / static_cast<double>(EventCounter)
              << "%" << std::endl;

    std::cout << "CountGoodLeadingMuonTrack (=good events): "
              << "\t" << CountGoodLeadingMuonTrack << "/" << CountMCTruthMuon << " = "
              << 100 * static_cast<double>(CountGoodLeadingMuonTrack) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountBadLeadingMuonTrack+CountNoRecoTracks+CountNoMuonTracks (=bad events): "
              << "\t" << CountBadLeadingMuonTrack + CountNoRecoTracks + CountNoMuonTracks << " = "
              << 100 *
                   static_cast<double>(CountBadLeadingMuonTrack + CountNoRecoTracks +
                                       CountNoMuonTracks) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountNoRecoTracks+CountNoMuonTracks: "
              << "\t" << CountNoRecoTracks + CountNoMuonTracks << " = "
              << 100 * static_cast<double>(CountNoRecoTracks + CountNoMuonTracks) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountTrackLengthTooShort: "
              << "\t" << CountTrackLengthTooShort << " = "
              << 100 * static_cast<double>(CountTrackLengthTooShort) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountCompleteness: "
              << "\t" << CountCompleteness << " = "
              << 100 * static_cast<double>(CountCompleteness) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountTrackLengthTooLong: "
              << "\t" << CountTrackLengthTooLong << " = "
              << 100 * static_cast<double>(CountTrackLengthTooLong) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountPurity: "
              << "\t" << CountPurity << " = "
              << 100 * static_cast<double>(CountPurity) / static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << std::endl;

    std::cout << "GoodLeadingMuonTrack+CountBadLeadingMuonTrackButLeadingPlusSecondGood (=good "
                 "events after stitching): "
              << "\t"
              << CountGoodLeadingMuonTrack + CountBadLeadingMuonTrackButLeadingPlusSecondGood << "/"
              << CountMCTruthMuon << " = "
              << 100 *
                   static_cast<double>(CountGoodLeadingMuonTrack +
                                       CountBadLeadingMuonTrackButLeadingPlusSecondGood) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountBadLeadingMuonTrack+CountNoRecoTracks+CountNoMuonTracks-"
                 "CountBadLeadingMuonTrackButLeadingPlusSecondGood (=bad events after stitching) : "
              << "\t"
              << CountBadLeadingMuonTrack + CountNoRecoTracks + CountNoMuonTracks -
                   CountBadLeadingMuonTrackButLeadingPlusSecondGood
              << " = "
              << 100 *
                   static_cast<double>(CountBadLeadingMuonTrack + CountNoRecoTracks +
                                       CountNoMuonTracks -
                                       CountBadLeadingMuonTrackButLeadingPlusSecondGood) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << std::endl;

    std::cout << "CountBadLeadingMuonTrackButLeadingPlusSecondGood: "
              << "\t" << CountBadLeadingMuonTrackButLeadingPlusSecondGood << " = "
              << 100 * static_cast<double>(CountBadLeadingMuonTrackButLeadingPlusSecondGood) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountBadLeadingMuonTrackAndOnlyOneMuonTrack: "
              << "\t" << CountBadLeadingMuonTrackAndOnlyOneMuonTrack << " = "
              << 100 * static_cast<double>(CountBadLeadingMuonTrackAndOnlyOneMuonTrack) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountBadLeadingMuonTrackAndLeadingPlusSecondBad: "
              << "\t" << CountBadLeadingMuonTrackAndLeadingPlusSecondBad << " = "
              << 100 * static_cast<double>(CountBadLeadingMuonTrackAndLeadingPlusSecondBad) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountNoRecoTracks: "
              << "\t" << CountNoRecoTracks << " = "
              << 100 * static_cast<double>(CountNoRecoTracks) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountNoMuonTracks: "
              << "\t" << CountNoMuonTracks << " = "
              << 100 * static_cast<double>(CountNoMuonTracks) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << std::endl;

    std::cout << "CountBadLeadingMuonTrackAndOnlyOneMuonTrackCompleteness: "
              << CountBadLeadingMuonTrackAndOnlyOneMuonTrackCompleteness << std::endl;
    std::cout << "CountBadLeadingMuonTrackAndOnlyOneMuonTrackPurity: "
              << CountBadLeadingMuonTrackAndOnlyOneMuonTrackPurity << std::endl;
    std::cout << "CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooShort: "
              << CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooShort << std::endl;
    std::cout << "CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooLong: "
              << CountBadLeadingMuonTrackAndOnlyOneMuonTrackTrackTooLong << std::endl;

    std::cout << std::endl;

    std::cout << "CountBadLeadingMuonTrackAndLeadingPlusSecondBadCompleteness: "
              << "\t" << CountBadLeadingMuonTrackAndLeadingPlusSecondBadCompleteness << " = "
              << 100 *
                   static_cast<double>(
                     CountBadLeadingMuonTrackAndLeadingPlusSecondBadCompleteness) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountBadLeadingMuonTrackAndLeadingPlusSecondBadPurity: "
              << "\t" << CountBadLeadingMuonTrackAndLeadingPlusSecondBadPurity << " = "
              << 100 * static_cast<double>(CountBadLeadingMuonTrackAndLeadingPlusSecondBadPurity) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooShort: "
              << "\t" << CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooShort << " = "
              << 100 *
                   static_cast<double>(
                     CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooShort) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << "CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooLong: "
              << "\t" << CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooLong << " = "
              << 100 *
                   static_cast<double>(
                     CountBadLeadingMuonTrackAndLeadingPlusSecondBadTrackTooLong) /
                   static_cast<double>(CountMCTruthMuon)
              << "%" << std::endl;

    std::cout << std::endl;

    std::cout << "GoodEvents1MuonTrack: " << GoodEvents1MuonTrack << std::endl;
    std::cout << "GoodEvents2MuonTrack: " << GoodEvents2MuonTrack << std::endl;
    std::cout << "GoodEvents3MuonTrack: " << GoodEvents3MuonTrack << std::endl;
    std::cout << "GoodEvents4OrMoreMuonTrack: " << GoodEvents4OrMoreMuonTrack << std::endl;

    std::cout << "BadEvents0MuonTrack: " << BadEvents0MuonTrack << std::endl;
    std::cout << "BadEvents1MuonTrack: " << BadEvents1MuonTrack << std::endl;
    std::cout << "BadEvents2MuonTrack: " << BadEvents2MuonTrack << std::endl;
    std::cout << "BadEvents3MuonTrack: " << BadEvents3MuonTrack << std::endl;
    std::cout << "BadEvents4OrMoreMuonTrack: " << BadEvents4OrMoreMuonTrack << std::endl;

    art::ServiceHandle<art::TFileService const> tfs;

    for (int i = 0; i < (NCriteriaBins + 1) / 2; ++i) {
      for (int j = 0; j < (NRecoTracksBins + 1) / 2; ++j) {
        if (i == 0) {
          h_Criteria_NRecoTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountNoRecoTracks);
          h_Criteria_NMuonTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountNoRecoTracks);
        }
        else if (i == 1) {
          h_Criteria_NRecoTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountNoMuonTracks);
          h_Criteria_NMuonTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountNoMuonTracks);
        }
        else if (i == 2) {
          h_Criteria_NRecoTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountCompleteness);
          h_Criteria_NMuonTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountCompleteness);
        }
        else if (i == 3) {
          h_Criteria_NRecoTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountPurity);
          h_Criteria_NMuonTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountPurity);
        }
        else if (i == 4) {
          h_Criteria_NRecoTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountTrackLengthTooShort);
          h_Criteria_NMuonTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountTrackLengthTooShort);
        }
        else if (i == 5) {
          h_Criteria_NRecoTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountTrackLengthTooLong);
          h_Criteria_NMuonTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountTrackLengthTooLong);
        }
        else if (i == 6) {
          h_Criteria_NRecoTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountRecoMuon);
          h_Criteria_NMuonTrack_den->SetBinContent(1 + 2 * i, 1 + 2 * j, CountRecoMuon);
        }
      }
    }

    h_Efficiency_ThetaXZ->Divide(h_ThetaXZ_num, h_ThetaXZ_den, 1, 1, "");
    h_Efficiency_ThetaYZ->Divide(h_ThetaYZ_num, h_ThetaYZ_den, 1, 1, "");
    h_Efficiency_SinThetaYZ->Divide(h_SinThetaYZ_num, h_SinThetaYZ_den, 1, 1, "");

    h_Efficiency_ThetaXZ_ThetaYZ->Divide(h_ThetaXZ_ThetaYZ_num, h_ThetaXZ_ThetaYZ_den, 1, 1, "");
    h_Efficiency_ThetaXZ_SinThetaYZ->Divide(
      h_ThetaXZ_SinThetaYZ_num, h_ThetaXZ_SinThetaYZ_den, 1, 1, "");

    h_Criteria_NRecoTrack->Divide(h_Criteria_NRecoTrack_num, h_Criteria_NRecoTrack_den, 1, 1, "");
    h_Criteria_NMuonTrack->Divide(h_Criteria_NMuonTrack_num, h_Criteria_NMuonTrack_den, 1, 1, "");

    h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk_num->Add(h_ThetaXZ_ThetaYZ_num, 1);
    h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk_num->Add(h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk, 1);
    h_Efficiency_ThetaXZ_ThetaYZ_LeadingPlusSecond->Divide(
      h_ThetaXZ_ThetaYZ_LeadingPlusSecondOk_num, h_ThetaXZ_ThetaYZ_den, 1, 1, "");

    h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk_num->Add(h_ThetaXZ_SinThetaYZ_num, 1);
    h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk_num->Add(h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk, 1);
    h_Efficiency_ThetaXZ_SinThetaYZ_LeadingPlusSecond->Divide(
      h_ThetaXZ_SinThetaYZ_LeadingPlusSecondOk_num, h_ThetaXZ_SinThetaYZ_den, 1, 1, "");

    h_Efficiency_ThetaXZ_ThetaYZ_DifferenceLeadingAndLeadingPlusSecond->Add(
      h_Efficiency_ThetaXZ_ThetaYZ_LeadingPlusSecond, 1);
    h_Efficiency_ThetaXZ_ThetaYZ_DifferenceLeadingAndLeadingPlusSecond->Add(
      h_Efficiency_ThetaXZ_ThetaYZ, -1);
  }
  //========================================================================
  DEFINE_ART_MODULE(MuonTrackingEff)

}

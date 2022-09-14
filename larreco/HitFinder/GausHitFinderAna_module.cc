////////////////////////////////////////////////////////////////////////
//
// Gaus(s)HitFinder class designed to analyze signal on a wire in the TPC
//
// jaasaadi@syr.edu
//
// Note: This is a rework of the original hit finder ana module
//       The only histograms that are saved are ones that can be used
//	 to make sure the hit finder is functioning...the rest is
//       outputted to a TTree for offline analysis.
////////////////////////////////////////////////////////////////////////

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/MCCheater/BackTrackerService.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// C++ includes
#include <string>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

constexpr int kMaxHits = 20000;

namespace hit {

  /// Base class for creation of raw signals on wires.
  class GausHitFinderAna : public art::EDAnalyzer {
  public:
    explicit GausHitFinderAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;
    void beginJob() override;

    std::string fHitFinderModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel;

    TH1F* fHitResidualAll;
    TH1F* fHitResidualAllAlt;
    TH1F* fNumberOfHitsPerEvent;
    TH2F* fPeakTimeVsWire;

    // ### TTree for offline analysis ###
    TTree* fHTree;

    // === Event Information ===
    Int_t fRun; // Run Number
    Int_t fEvt; // Event Number

    // === Wire Information ====
    Float_t fWireTotalCharge; // Charge on all wires

    // === Hit Information ===
    Int_t fnHits;                      // Number of Hits in the Event
    Int_t fWire[kMaxHits];             // Wire Number
    Float_t fStartTime[kMaxHits];      // Start Time
    Float_t fEndTime[kMaxHits];        // End Time
    Float_t fPeakTime[kMaxHits];       // Peak Time
    Float_t fPeakTimeUncert[kMaxHits]; // Peak Time Uncertainty
    Float_t fCharge[kMaxHits];         // Charge of this hit
    Float_t fChargeUncert[kMaxHits];   // Charge Uncertainty of this hit
    Int_t fMultiplicity[kMaxHits];     // Hit pulse multiplicity
    Float_t fGOF[kMaxHits];            // Goodness of Fit (Chi2/NDF)

    // === Total Hit Information ===
    Float_t fTotalHitChargePerEvent; //Total charge recorded in each event

    // === Truth Hit Info from BackTrackerService ===
    Float_t fTruePeakPos[kMaxHits]; // Truth Time Tick info from BackTrackerService

  }; // class GausHitFinderAna

  //-------------------------------------------------
  GausHitFinderAna::GausHitFinderAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {
    fHitFinderModuleLabel = pset.get<std::string>("HitsModuleLabel");
    fLArG4ModuleLabel = pset.get<std::string>("LArGeantModuleLabel");
    fCalDataModuleLabel = pset.get<std::string>("CalDataModuleLabel");
  }

  //-------------------------------------------------
  void GausHitFinderAna::beginJob()
  {
    art::ServiceHandle<art::TFileService const> tfs;
    fHitResidualAll = tfs->make<TH1F>("fHitResidualAll", "Hit Residual All", 1600, -400, 400);
    fHitResidualAllAlt = tfs->make<TH1F>("fHitResidualAllAlt", "Hit Residual All", 1600, -400, 400);
    fNumberOfHitsPerEvent =
      tfs->make<TH1F>("fNumberOfHitsPerEvent", "Number of Hits in Each Event", 10000, 0, 10000);
    fPeakTimeVsWire =
      tfs->make<TH2F>("fPeakTimeVsWire", "Peak Time vs Wire Number", 3200, 0, 3200, 9500, 0, 9500);

    fHTree = tfs->make<TTree>("HTree", "HTree");
    fHTree->Branch("Evt", &fEvt, "Evt/I");
    fHTree->Branch("Run", &fRun, "Run/I");
    fHTree->Branch("WireTotalCharge", &fWireTotalCharge, "WireTotalCharge/F");

    // === Hit Info ===
    fHTree->Branch("nHits", &fnHits, "nHits/I");
    fHTree->Branch("Wire", &fWire, "Wire[nHits]/I");
    fHTree->Branch("StartTime", &fStartTime, "fStartTime[nHits]/F");
    fHTree->Branch("EndTime", &fEndTime, "fEndTime[nHits]/F");
    fHTree->Branch("PeakTime", &fPeakTime, "fPeakTime[nHits]/F");
    fHTree->Branch("PeakTimeUncert", &fPeakTimeUncert, "fPeakTimeUncert[nHits]/F");
    fHTree->Branch("Charge", &fCharge, "fCharge[nHits]/F");
    fHTree->Branch("ChargeUncert", &fChargeUncert, "fChargeUncert[nHits]/F");
    fHTree->Branch("Multiplicity", &fMultiplicity, "fMultiplicity[nHits]/I");
    fHTree->Branch("GOF", &fGOF, "fGOF[nHits]/F");

    // === Total Hit Information ===
    fHTree->Branch("TotalHitChargePerEvent", &fTotalHitChargePerEvent, "TotalHitChargePerEvent/F");

    // === Truth Hit Information from BackTrackerService ===
    fHTree->Branch("TruePeakPos", &fTruePeakPos, "fTruePeakPos[nHits]/F");
  }

  //-------------------------------------------------
  void GausHitFinderAna::analyze(const art::Event& evt)
  {
    // ### TTree Run/Event ###
    fEvt = evt.id().event();
    fRun = evt.run();

    art::ServiceHandle<geo::Geometry const> geom;
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);

    art::Handle<std::vector<recob::Wire>> wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel, wireVecHandle);

    // Charge directly from wire info
    float TotWireCharge = 0;

    for (size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++) {
      art::Ptr<recob::Wire> wire(wireVecHandle, wireIter);
      std::vector<float> signal(wire->Signal());

      for (auto timeIter = signal.begin(); timeIter + 2 < signal.end(); timeIter++) {

        if (*timeIter < 2) { continue; }

        TotWireCharge += *timeIter;
      }
    }

    fWireTotalCharge = TotWireCharge;

    // Reconstructed hit information
    art::Handle<std::vector<recob::Hit>> hitHandle;
    evt.getByLabel(fHitFinderModuleLabel, hitHandle);

    std::vector<art::Ptr<recob::Hit>> hits;
    art::fill_ptr_vector(hits, hitHandle);

    float TotCharge = 0;
    int hitCount = 0;
    fnHits = hitHandle->size();
    fNumberOfHitsPerEvent->Fill(hitHandle->size());

    for (size_t numHit = 0; numHit < hitHandle->size(); ++numHit) {
      // === Finding Channel associated with the hit ===
      art::Ptr<recob::Hit> hit(hitHandle, numHit);

      fWire[hitCount] = hit->WireID().Wire;
      fStartTime[hitCount] = hit->PeakTimeMinusRMS();
      fEndTime[hitCount] = hit->PeakTimePlusRMS();
      fPeakTime[hitCount] = hit->PeakTime();
      fPeakTimeUncert[hitCount] = hit->SigmaPeakTime();
      fCharge[hitCount] = hit->Integral();
      fChargeUncert[hitCount] = hit->SigmaIntegral();
      fMultiplicity[hitCount] = hit->Multiplicity();
      fGOF[hitCount] = hit->GoodnessOfFit();

      hitCount++;
      TotCharge += hit->Integral();

      fPeakTimeVsWire->Fill(hit->WireID().Wire, hit->PeakTime());
    } //<---End numHit
    fTotalHitChargePerEvent = TotCharge;

    // Truth hit info from BackTracker

    Float_t TruthHitTime = 0, TruthHitCalculated = 0;
    int count = 0;

    double time_tick = sampling_rate(clock_data) / 1000.;
    double drift_velocity = det_prop.DriftVelocity(det_prop.Efield(), det_prop.Temperature());

    for (size_t nh = 0; nh < hitHandle->size(); nh++) {
      // === Finding Channel associated with the hit ===
      art::Ptr<recob::Hit> hitPoint(hitHandle, nh);
      auto const& planeID = hitPoint->WireID().asPlaneID();

      // ===================================================================
      // Using Track IDE's to locate the XYZ location from truth information
      // ===================================================================
      std::vector<sim::TrackIDE> trackides;
      std::vector<double> xyz;
      try {
        art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
        trackides = bt_serv->HitToTrackIDEs(clock_data, hitPoint);
        xyz = bt_serv->HitToXYZ(clock_data, hitPoint);
      }
      catch (cet::exception const&) {
        mf::LogWarning("GausHitFinderAna") << "BackTrackerService Failed";
        continue;
      }

      // ==============================================================
      // Calculating the truth tick position of the hit using 2 methods
      // Method 1: ConvertXtoTicks from the detector properties package
      // Method 2: Actually do the calculation myself to double check things
      // ==============================================================

      // ### Method 1 ###
      TruthHitTime = det_prop.ConvertXToTicks(xyz[0], planeID);

      // ### Method 2 ###
      // ================================================
      // Establishing the x-position of the current plane
      // ================================================
      auto const pos = geom->Plane(planeID).GetBoxCenter();
      double planePos_timeCorr = (pos.X() / drift_velocity) * (1. / time_tick) + 60;
      //<---x position of plane / drift velocity + 60 (Trigger offset)

      TruthHitCalculated = ((xyz[0]) / (drift_velocity * time_tick)) + planePos_timeCorr;

      fTruePeakPos[count] = TruthHitTime;
      count++;
      double hitresid = ((TruthHitTime - hitPoint->PeakTime()) / hitPoint->SigmaPeakTime());
      fHitResidualAll->Fill(hitresid);

      double hitresidAlt =
        ((TruthHitCalculated - hitPoint->PeakTime()) / hitPoint->SigmaPeakTime());
      fHitResidualAllAlt->Fill(hitresidAlt);

    } //<---End nh loop

    fHTree->Fill();
  } // end analyze method

  DEFINE_ART_MODULE(GausHitFinderAna)

} // end of hit namespace

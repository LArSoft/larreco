////////////////////////////////////////////////////////////////////////
//
// MagDriftAna class
//
// dmckee@phys.ksu.edu
//
////////////////////////////////////////////////////////////////////////

// C++ std library includes
#include <algorithm>
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

// Root Includes
#include "TH2.h"
#include "TLine.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"

/// Detector simulation of raw signals on wires
namespace hit {

  /// Base class for creation of raw signals on wires.
  class MagDriftAna : public art::EDAnalyzer {
  public:
    explicit MagDriftAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;
    void endJob() override;

    // intilize the histograms
    //
    // Can't be done in Begin job because I want to use LArProperties
    // which used the database, so I test and run on each
    // event. Wasteful and silly, but at least it *works*.
    void ensureHists(art::Event const& evt, detinfo::DetectorClocksData const& clockData);

    std::string fFFTHitFinderModuleLabel;
    std::string fLArG4ModuleLabel;

    // Flag for initialization done, because we set up histograms the
    // first time through beginRun() so that we can use the
    // database...
    bool initDone{false};

    // Drift properties
    double fDirCosY{0.};
    double fDirCosZ{0.};

    TH1D* fChargeXpos{nullptr}; // << position of the MC Truth charge deposition
    TH1D* fChargeYpos{nullptr};
    TH1D* fChargeZpos{nullptr};
    TH1D* fHitZpos{nullptr}; // << Z position of the recorded hit (from the
                             //    z-sensitive wire)

    TH1D* fDriftDeltaZ{nullptr}; // << Difference in MC charge Z and recorded hit Z
    TH1D* fDeltaZoverX{nullptr}; // << Delta Z as a function of drift distance
    TH2D* fDeltaZvsX{nullptr};

    // Same as above, but only for long drift distances (greater than
    // 4/5 of the detector)
    TH1D* fDriftDeltaZAway{nullptr}; // << Difference in MC charge Z and recorded hit Z
    TH1D* fDeltaZoverXAway{nullptr}; // << Delta Z as a function of drift distance

  }; // class MagdriftAna

  //-------------------------------------------------
  MagDriftAna::MagDriftAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fFFTHitFinderModuleLabel{pset.get<std::string>("HitsModuleLabel")}
    , fLArG4ModuleLabel{pset.get<std::string>("LArGeantModuleLabel")}
  {}

  //-------------------------------------------------
  void MagDriftAna::ensureHists(art::Event const& evt, detinfo::DetectorClocksData const& clockData)
  {
    if (initDone) return; // Bail if we've already done this.
    initDone = true;      // Insure that we bail later on

    art::ServiceHandle<art::TFileService const> tfs;

    // Find magnetic field related corrections
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    // art::ServiceHandle<mag::MagneticField const> MagField;
    art::ServiceHandle<mag::MagneticFieldService const> MagFieldHandle;
    auto const* MagField = MagFieldHandle->provider();
    double Efield = detProp.Efield();
    double Temperature = detProp.Temperature();
    double DriftVelocity = detProp.DriftVelocity(Efield, Temperature) / 1000.;

    // MagneticField::FieldAtPoint() returns (0, 0, 0) if there is no
    // field at the requested point, so these direction cosines are
    // 0 if there is no field
    fDirCosY = -DriftVelocity * MagField->FieldAtPoint().z() / Efield;
    fDirCosZ = +DriftVelocity * MagField->FieldAtPoint().y() / Efield;
    MF_LOG_VERBATIM("MagDriftAna") << "Drift ratios: "
                                   << "dY/dX = " << fDirCosY << ", "
                                   << "dZ/dX = " << fDirCosZ;

    // geometry data.
    art::ServiceHandle<geo::Geometry const> geom;
    // assumes all TPCs are the same
    auto const& tpc_0 = geom->TPC();
    double width = 2 * tpc_0.HalfWidth();
    double halfHeight = tpc_0.HalfHeight();
    double length = tpc_0.Length();

    double zScale = std::max(fDirCosZ / 2.0, 4e-4);

    // Assumes microboone dimensions. Ideally we'd fix this later...
    fChargeXpos =
      tfs->make<TH1D>("hChargeXpos", "MC X charge depositions; X (cm); Events", 101, 0.0, width);
    fChargeYpos = tfs->make<TH1D>(
      "hChargeYpos", "MC Y charge depositions; Y (cm); Events", 101, -halfHeight, halfHeight);
    fChargeZpos =
      tfs->make<TH1D>("hChargeZpos", "MC Z charge depositions; Z (cm); Events", 101, 0.0, length);
    fHitZpos = tfs->make<TH1D>("hHitZpos", "Z charge collection; Z (cm); Events", 101, 0.0, length);

    fDriftDeltaZ = tfs->make<TH1D>("hDriftDeltaZ",
                                   "Z drift of charge; delta Z (cm); Events",
                                   101,
                                   -5 * zScale * width,
                                   5 * zScale * width);
    fDeltaZoverX = tfs->make<TH1D>(
      "hDeltaZoverX", "Z drift of charge; delta Z/X; Events", 51, -10 * zScale, 10 * zScale);
    fDeltaZvsX = tfs->make<TH2D>("hDeltaZvsX",
                                 "delta Z vs X; X (cm); delta Z (cm), Events",
                                 51,
                                 0.0,
                                 width,
                                 51,
                                 -20 * zScale,
                                 20 * zScale);

    // Some stats only when the Xdrift is large (more than 3/4)
    fDriftDeltaZAway = tfs->make<TH1D>("hDriftDeltaZAway",
                                       "Z drift of charge (long drift); delta Z (cm); Events",
                                       101,
                                       -5 * zScale * width,
                                       5 * zScale * width);
    fDeltaZoverXAway = tfs->make<TH1D>("hDeltaZoverXAway",
                                       "Z drift of charge (long drift); delta Z/X; Events",
                                       51,
                                       -10 * zScale,
                                       10 * zScale);
  }

  //-------------------------------------------------
  void MagDriftAna::endJob()
  {
    // Add a line on the deltaZ/X graph to denote the calculated value
    // of the drift ration
    TLine* l = new TLine(fDirCosZ, 0, fDirCosZ, 1.05 * fDeltaZoverX->GetMaximum());
    l->SetLineColor(kRed);
    l->SetLineStyle(kDotted);
    fDeltaZoverX->GetListOfFunctions()->Add(l);

    // I know this looks like a memory leak, but each historgram needs
    // it's own copy of the line to prevent double freeing by the
    // framework...
    l = new TLine(fDirCosZ, 0, fDirCosZ, 1.05 * fDeltaZoverX->GetMaximum());

    l->SetLineColor(kRed);
    l->SetLineStyle(kDotted);
    fDeltaZoverXAway->GetListOfFunctions()->Add(l);
  }

  //-------------------------------------------------
  void MagDriftAna::analyze(const art::Event& evt)
  {
    if (evt.isRealData()) {
      throw cet::exception("MagDriftAna: ") << "Not for use on Data yet...\n";
    }

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    ensureHists(evt, clockData);

    art::Handle<std::vector<recob::Hit>> hitHandle;
    evt.getByLabel(fFFTHitFinderModuleLabel, hitHandle);

    art::ServiceHandle<geo::Geometry const> geom;

    // We're going to want to compare the reconstructed Z with the
    // simulted Z. For that purpose we use the simultion backtracking.

    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

    std::vector<art::Ptr<recob::Hit>> hits;
    art::fill_ptr_vector(hits, hitHandle);

    geo::WireID hitWireID;

    //++++++++++
    // Loop over the hits (data) and fill histos
    //++++++++++
    for (auto itr : hits) {

      hitWireID = itr->WireID();
      // By assumption the drift occurs only in the z-direction, so
      // we can get all the info we need from the z-measug plane.
      if (hitWireID.Plane != (geom->Nplanes() - 1)) continue;

      // Charge collected at the wire
      //
      // Exactly once for each recob::Hit
      auto const w0pos = geom->Plane(hitWireID).Wire(0).GetCenter();
      double HitZpos = w0pos.Z() + hitWireID.Wire * geom->TPC(hitWireID).WirePitch();
      double Charge = itr->Integral();
      fHitZpos->Fill(HitZpos, Charge);

      // Charge deposition in the detector
      std::vector<double> xyz = bt_serv->HitToXYZ(clockData, itr);
      fChargeXpos->Fill(xyz[0], Charge);
      fChargeYpos->Fill(xyz[1], Charge);
      double ChargeZpos = xyz[2];
      fChargeZpos->Fill(ChargeZpos, Charge);

      // Delta-Z
      //
      // Compares the collected z position from the wire to the
      // simulated z position
      double DeltaZ = HitZpos - ChargeZpos;
      fDriftDeltaZ->Fill(DeltaZ, Charge);
      // Delta Z correlation with X
      fDeltaZoverX->Fill(DeltaZ / xyz[0], Charge);
      fDeltaZvsX->Fill(xyz[0], DeltaZ, Charge);
      // The X related histograms know the dimensions of the
      // detector, so we use them to set the "away" limit
      if (xyz[0] > (fChargeYpos->GetXaxis()->GetXmax() * 0.80)) {
        fDriftDeltaZAway->Fill(DeltaZ, Charge);
        fDeltaZoverXAway->Fill(DeltaZ / xyz[0], Charge);
      }

    } // loop on Hits

    return;
  } // end analyze method

  DEFINE_ART_MODULE(MagDriftAna)

} // end of hit namespace

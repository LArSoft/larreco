//
// Name: SpacePointAna_module.cc
//
// Purpose: Module SpacePointAna.
//
//
// Configuration parameters.
//
//  HitModuleLabel:        Hit module label.
//  UseClusterHits:        If true, use clustered hits (otherwise all hits).
//  ClusterModuleLabel:    Cluster module label.
//  UseMC:                 Use MC truth information.
//  SpacePointAlgTime:     SpacePointAlg configuration (loose time cuts).
//  SpacePointAlgSep:      SpacePointAlg configuration (loose separation cuts).
//  SpacePointAlgDefault:  SpacePointAlg configuration (default cuts).
//
// Created: 2-Aug-2011  H. Greenlee
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/ArtDataHelper/ToElement.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "TH1F.h"
#include "TH2F.h"

#include "range/v3/view.hpp"

using lar::to_element;
using ranges::views::transform;

namespace trkf {

  class SpacePointAna : public art::EDAnalyzer {
  public:
    explicit SpacePointAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;

    void bookHistograms(bool mc);

    // Fcl Attributes.

    const SpacePointAlg fSptalgTime;    // Algorithm object (increased time cut).
    const SpacePointAlg fSptalgSep;     // Algorithm object (increased sepataion cut).
    const SpacePointAlg fSptalgDefault; // Algorithm object (default cuts).
    std::string fHitModuleLabel;
    bool fUseClusterHits;
    std::string fClusterModuleLabel;
    bool fUseMC;
    double fMinX; // Minimum x.
    double fMaxX; // Maximum x.
    double fMinY; // Minimum y.
    double fMaxY; // Maximum y.
    double fMinZ; // Minimum z.
    double fMaxZ; // Maximum z.

    // Histograms.

    bool fBooked;    // Have histograms been booked yet?
    TH1F* fHDTUE;    // U-drift electrons time difference.
    TH1F* fHDTVE;    // V-drift electrons time difference.
    TH1F* fHDTWE;    // W-drift electrons time difference.
    TH1F* fHDTUPull; // U-drift electrons time pull.
    TH1F* fHDTVPull; // V-drift electrons time pull.
    TH1F* fHDTWPull; // W-drift electrons time pull.
    TH1F* fHDTUV;    // U-V time difference.
    TH1F* fHDTVW;    // V-W time difference.
    TH1F* fHDTWU;    // W-U time difference.
    TH2F* fHDTUVU;   // U-V time difference vs. U.
    TH2F* fHDTUVV;   // U-V time difference vs. V.
    TH2F* fHDTVWV;   // V-W time difference vs. V.
    TH2F* fHDTVWW;   // V-W time difference vs. W.
    TH2F* fHDTWUW;   // W-U time difference vs. W.
    TH2F* fHDTWUU;   // W-U time difference vs. U.
    TH1F* fHS;       // Spatial separation.
    TH1F* fHchisq;   // Space point chisquare.
    TH1F* fHx;       // X position.
    TH1F* fHy;       // Y position.
    TH1F* fHz;       // Z position.
    TH1F* fHAmpU;    // U hit amplitude.
    TH1F* fHAmpV;    // V hit amplitude.
    TH1F* fHAmpW;    // W hit amplitude.
    TH1F* fHAreaU;   // U hit area.
    TH1F* fHAreaV;   // V hit area.
    TH1F* fHAreaW;   // W hit area.
    TH1F* fHSumU;    // U hit sum ADC.
    TH1F* fHSumV;    // V hit sum ADC.
    TH1F* fHSumW;    // W hit sum ADC.
    TH1F* fHMCdx;    // X residual (reco vs. mc truth).
    TH1F* fHMCdy;    // Y residual (reco vs. mc truth).
    TH1F* fHMCdz;    // Z residual (reco vs. mc truth).
    TH1F* fHMCxpull; // X pull (reco vs. mc truth).
    TH1F* fHMCypull; // Y pull (reco vs. mc truth).
    TH1F* fHMCzpull; // Z pull (reco vs. mc truth).

    // Statistics.

    int fNumEvent;
  };

  DEFINE_ART_MODULE(SpacePointAna)

  SpacePointAna::SpacePointAna(const fhicl::ParameterSet& pset)
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    : EDAnalyzer(pset)
    , fSptalgTime(pset.get<fhicl::ParameterSet>("SpacePointAlgTime"))
    , fSptalgSep(pset.get<fhicl::ParameterSet>("SpacePointAlgSep"))
    , fSptalgDefault(pset.get<fhicl::ParameterSet>("SpacePointAlgDefault"))
    , fHitModuleLabel(pset.get<std::string>("HitModuleLabel"))
    , fUseClusterHits(pset.get<bool>("UseClusterHits"))
    , fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel"))
    , fUseMC(pset.get<bool>("UseMC"))
    , fMinX(pset.get<double>("MinX", -1.e10))
    , fMaxX(pset.get<double>("MaxX", 1.e10))
    , fMinY(pset.get<double>("MinY", -1.e10))
    , fMaxY(pset.get<double>("MaxY", 1.e10))
    , fMinZ(pset.get<double>("MinZ", -1.e10))
    , fMaxZ(pset.get<double>("MaxZ", 1.e10))
    , fBooked(false)
    , fHDTUE(0)
    , fHDTVE(0)
    , fHDTWE(0)
    , fHDTUPull(0)
    , fHDTVPull(0)
    , fHDTWPull(0)
    , fHDTUV(0)
    , fHDTVW(0)
    , fHDTWU(0)
    , fHDTUVU(0)
    , fHDTUVV(0)
    , fHDTVWV(0)
    , fHDTVWW(0)
    , fHDTWUW(0)
    , fHDTWUU(0)
    , fHS(0)
    , fHchisq(0)
    , fHx(0)
    , fHy(0)
    , fHz(0)
    , fHAmpU(0)
    , fHAmpV(0)
    , fHAmpW(0)
    , fHAreaU(0)
    , fHAreaV(0)
    , fHAreaW(0)
    , fHSumU(0)
    , fHSumV(0)
    , fHSumW(0)
    , fHMCdx(0)
    , fHMCdy(0)
    , fHMCdz(0)
    , fHMCxpull(0)
    , fHMCypull(0)
    , fHMCzpull(0)
    , fNumEvent(0)
  {
    mf::LogInfo("SpacePointAna") << "SpacePointAna configured with the following parameters:\n"
                                 << "  HitModuleLabel = " << fHitModuleLabel << "\n"
                                 << "  UseClusterHits = " << fUseClusterHits << "\n"
                                 << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
                                 << "  UseMC = " << fUseMC;
  }

  void SpacePointAna::bookHistograms(bool mc)
  {
    if (!fBooked) {
      fBooked = true;

      art::ServiceHandle<geo::Geometry const> geom;
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory dir = tfs->mkdir("sptana", "SpacePointAna histograms");

      unsigned int nwiresU = 0, nwiresV = 0, nwiresW = 0;

      // Figure out the number of wires in U, V, and W planes.

      // Loop over cryostats, tpcs, and planes.

      for (auto const& pgeom : geom->Iterate<geo::PlaneGeo>()) {
        unsigned int const nwires = pgeom.Nwires();
        switch (pgeom.View()) {
        case geo::kU: nwiresU = nwires; break;
        case geo::kV: nwiresV = nwires; break;
        case geo::kZ: nwiresW = nwires; break;
        default: {}
        }
      }

      if (mc && fUseMC) {
        fHDTUE = dir.make<TH1F>("MCDTUE", "U-Drift Electrons Time Difference", 100, -5., 5.);
        fHDTVE = dir.make<TH1F>("MCDTVE", "V-Drift Electrons Time Difference", 100, -5., 5.);
        fHDTWE = dir.make<TH1F>("MCDTWE", "W-Drift Electrons Time Difference", 100, -5., 5.);
        fHDTUPull = dir.make<TH1F>("MCDTUPull", "U-Drift Electrons Time Pull", 100, -50., 50.);
        fHDTVPull = dir.make<TH1F>("MCDTVPull", "V-Drift Electrons Time Pull", 100, -50., 50.);
        fHDTWPull = dir.make<TH1F>("MCDTWPull", "W-Drift Electrons Time Pull", 100, -50., 50.);
      }
      if (!fSptalgTime.merge()) {
        fHDTUV = dir.make<TH1F>("DTUV", "U-V time difference", 100, -20., 20.);
        fHDTVW = dir.make<TH1F>("DTVW", "V-W time difference", 100, -20., 20.);
        fHDTWU = dir.make<TH1F>("DTWU", "W-U time difference", 100, -20., 20.);
        fHDTUVU = dir.make<TH2F>(
          "DTUVU", "U-V time difference vs. U", 100, 0., double(nwiresU), 100, -20., 20.);
        fHDTUVV = dir.make<TH2F>(
          "DTUVV", "U-V time difference vs. V", 100, 0., double(nwiresV), 100, -20., 20.);
        fHDTVWV = dir.make<TH2F>(
          "DTVWV", "V-W time difference vs. V", 100, 0., double(nwiresV), 100, -20., 20.);
        fHDTVWW = dir.make<TH2F>(
          "DTVWW", "V-W time difference vs. W", 100, 0., double(nwiresW), 100, -20., 20.);
        fHDTWUW = dir.make<TH2F>(
          "DTWUW", "W-U time difference vs. W", 100, 0., double(nwiresW), 100, -20., 20.);
        fHDTWUU = dir.make<TH2F>(
          "DTWUU", "W-U time difference vs. U", 100, 0., double(nwiresU), 100, -20., 20.);
        fHS = dir.make<TH1F>("DS", "Spatial Separation", 100, -2., 2.);
      }
      fHchisq = dir.make<TH1F>("chisq", "Chisquare", 100, 0., 20.);

      fHx = dir.make<TH1F>("xpos", "X Position", 100, 0., 2. * geom->DetHalfWidth());
      fHy =
        dir.make<TH1F>("ypos", "Y Position", 100, -geom->DetHalfHeight(), geom->DetHalfHeight());
      fHz = dir.make<TH1F>("zpos", "Z Position", 100, 0., geom->DetLength());
      fHAmpU = dir.make<TH1F>("ampU", "U Hit Amplitude", 50, 0., 50.);
      fHAmpV = dir.make<TH1F>("ampV", "V Hit Amplitude", 50, 0., 50.);
      fHAmpW = dir.make<TH1F>("ampW", "W Hit Amplitude", 50, 0., 50.);
      fHAreaU = dir.make<TH1F>("areaU", "U Hit Area", 100, 0., 500.);
      fHAreaV = dir.make<TH1F>("areaV", "V Hit Area", 100, 0., 500.);
      fHAreaW = dir.make<TH1F>("areaW", "W Hit Area", 100, 0., 500.);
      fHSumU = dir.make<TH1F>("sumU", "U Hit Sum ADC", 100, 0., 500.);
      fHSumV = dir.make<TH1F>("sumV", "V Hit Sum ADC", 100, 0., 500.);
      fHSumW = dir.make<TH1F>("sumW", "W Hit Sum ADC", 100, 0., 500.);
      if (mc && fUseMC) {
        fHMCdx = dir.make<TH1F>("MCdx", "X MC Residual", 100, -1., 1.);
        fHMCdy = dir.make<TH1F>("MCdy", "Y MC Residual", 100, -1., 1.);
        fHMCdz = dir.make<TH1F>("MCdz", "Z MC Residual", 100, -1., 1.);
        fHMCxpull = dir.make<TH1F>("MCxpull", "X MC Pull", 100, -50., 50.);
        fHMCypull = dir.make<TH1F>("MCypull", "Y MC Pull", 100, -50., 50.);
        fHMCzpull = dir.make<TH1F>("MCzpull", "Z MC Pull", 100, -50., 50.);
      }
    }
  }

  void SpacePointAna::analyze(const art::Event& evt)
  {
    ++fNumEvent;

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();
    bookHistograms(mc);

    // Get Services.

    art::ServiceHandle<geo::Geometry const> geom;

    // Get Hits.

    art::PtrVector<recob::Hit> hits;

    if (fUseClusterHits) {

      // Get clusters.

      art::Handle<std::vector<recob::Cluster>> clusterh;
      evt.getByLabel(fClusterModuleLabel, clusterh);

      // Get hits from all clusters.
      art::FindManyP<recob::Hit> fm(clusterh, evt, fClusterModuleLabel);

      if (clusterh.isValid()) {
        int nclus = clusterh->size();

        for (int i = 0; i < nclus; ++i) {
          art::Ptr<recob::Cluster> pclus(clusterh, i);
          std::vector<art::Ptr<recob::Hit>> clushits = fm.at(i);
          int nhits = clushits.size();
          hits.reserve(hits.size() + nhits);

          for (std::vector<art::Ptr<recob::Hit>>::const_iterator ihit = clushits.begin();
               ihit != clushits.end();
               ++ihit) {
            hits.push_back(*ihit);
          }
        }
      }
    }
    else {

      // Get unclustered hits.

      art::Handle<std::vector<recob::Hit>> hith;
      evt.getByLabel(fHitModuleLabel, hith);
      if (hith.isValid()) {
        int nhits = hith->size();
        hits.reserve(nhits);

        for (int i = 0; i < nhits; ++i)
          hits.push_back(art::Ptr<recob::Hit>(hith, i));
      }
    }

    // Fill histograms that don't depend on space points.

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    if (mc && fUseMC) {

      art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

      // Loop over hits and fill hit-electron time difference histogram.

      for (auto const& hitPtr : hits) {
        const recob::Hit& hit = *hitPtr;

        double tpeak = hit.PeakTime();
        double terr = hit.SigmaPeakTime();

        bool tav_ok = true;
        double tav = 0.;
        try {
          std::vector<double> hitxyz = bt_serv->HitToXYZ(clockData, hit);
          tav = detProp.ConvertXToTicks(
            hitxyz[0], hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat);
        }
        catch (cet::exception& x) {
          tav_ok = false;
        }
        if (tav_ok) {
          if (hit.View() == geo::kU) {
            fHDTUE->Fill(tpeak - tav);
            fHDTUPull->Fill((tpeak - tav) / terr);
          }
          else if (hit.View() == geo::kV) {
            fHDTVE->Fill(tpeak - tav);
            fHDTVPull->Fill((tpeak - tav) / terr);
          }
          else if (hit.View() == geo::kZ) {
            fHDTWE->Fill(tpeak - tav);
            fHDTWPull->Fill((tpeak - tav) / terr);
          }
          else
            throw cet::exception("SpacePointAna") << "Bad view = " << hit.View() << "\n";
        }
      }
    }

    std::vector<recob::SpacePoint> spts1; // For time histograms.
    std::vector<recob::SpacePoint> spts2; // For separation histogram.
    std::vector<recob::SpacePoint> spts3; // Default cuts.

    // If nonzero time cut is specified, make space points using that
    // time cut (for time histograms).

    if (!fSptalgTime.merge()) {
      if (mc && fUseMC)
        fSptalgTime.makeMCTruthSpacePoints(clockData, detProp, hits, spts1);
      else
        fSptalgTime.makeSpacePoints(clockData, detProp, hits, spts1);

      // Report number of space points.

      MF_LOG_DEBUG("SpacePointAna")
        << "Found " << spts1.size() << " space points using special time cut.";
    }

    // If nonzero separation cut is specified, make space points using that
    // separation cut (for separation histogram).

    if (!fSptalgSep.merge()) {
      if (mc && fUseMC)
        fSptalgSep.makeMCTruthSpacePoints(clockData, detProp, hits, spts2);
      else
        fSptalgSep.makeSpacePoints(clockData, detProp, hits, spts2);

      // Report number of space points.

      MF_LOG_DEBUG("SpacePointAna")
        << "Found " << spts2.size() << " space points using special seperation cut.";
    }

    // Make space points using default cuts.

    if (mc && fUseMC)
      fSptalgDefault.makeMCTruthSpacePoints(clockData, detProp, hits, spts3);
    else
      fSptalgDefault.makeSpacePoints(clockData, detProp, hits, spts3);

    // Report number of space points.

    MF_LOG_DEBUG("SpacePointAna") << "Found " << spts3.size()
                                  << " space points using default cuts.";

    if (!fSptalgTime.merge()) {

      // Loop over space points and fill time histograms.

      for (auto const& spt : spts1) {
        if (spt.XYZ()[0] < fMinX || spt.XYZ()[0] > fMaxX || spt.XYZ()[1] < fMinY ||
            spt.XYZ()[1] > fMaxY || spt.XYZ()[2] < fMinZ || spt.XYZ()[2] > fMaxZ)
          continue;

        // Get hits associated with this SpacePoint.

        const art::PtrVector<recob::Hit>& spthits = fSptalgTime.getAssociatedHits(spt);

        // Make a double loop over hits and fill hit time difference histograms.

        for (auto const& hit1 : spthits | transform(to_element)) {
          geo::WireID hit1WireID = hit1.WireID();
          unsigned int tpc1 = hit1WireID.TPC;
          unsigned int plane1 = hit1WireID.Plane;
          unsigned int wire1 = hit1WireID.Wire;

          geo::View_t view1 = hit1.View();
          double t1 = fSptalgTime.correctedTime(detProp, hit1);

          for (auto const& hit2 : spthits | transform(to_element)) {
            geo::WireID hit2WireID = hit2.WireID();
            unsigned int tpc2 = hit2WireID.TPC;
            unsigned int plane2 = hit2WireID.Plane;
            unsigned int wire2 = hit2WireID.Wire;

            // Require same tpc, different view.

            if (tpc1 == tpc2 && plane1 != plane2) {

              geo::View_t view2 = hit2.View();
              double t2 = fSptalgTime.correctedTime(detProp, hit2);

              if (view1 == geo::kU) {
                if (view2 == geo::kV) {
                  fHDTUV->Fill(t1 - t2);
                  fHDTUVU->Fill(double(wire1), t1 - t2);
                  fHDTUVV->Fill(double(wire2), t1 - t2);
                }
                if (view2 == geo::kZ) {
                  fHDTWU->Fill(t2 - t1);
                  fHDTWUW->Fill(double(wire2), t2 - t1);
                  fHDTWUU->Fill(double(wire1), t2 - t1);
                }
              }
              if (view1 == geo::kV) {
                if (view2 == geo::kZ) {
                  fHDTVW->Fill(t1 - t2);
                  fHDTVWV->Fill(double(wire1), t1 - t2);
                  fHDTVWW->Fill(double(wire2), t1 - t2);
                }
                if (view2 == geo::kU) {
                  fHDTUV->Fill(t2 - t1);
                  fHDTUVU->Fill(double(wire2), t2 - t1);
                  fHDTUVV->Fill(double(wire1), t2 - t1);
                }
              }
              if (view1 == geo::kZ) {
                if (view2 == geo::kU) {
                  fHDTWU->Fill(t1 - t2);
                  fHDTWUW->Fill(double(wire1), t1 - t2);
                  fHDTWUU->Fill(double(wire2), t1 - t2);
                }
                if (view2 == geo::kV) {
                  fHDTVW->Fill(t2 - t1);
                  fHDTVWV->Fill(double(wire2), t2 - t1);
                  fHDTVWW->Fill(double(wire1), t2 - t1);
                }
              }
            }
          }
        }
      }

      // Loop over space points and fill seperation histograms.

      for (auto const& spt : spts2) {
        if (spt.XYZ()[0] < fMinX || spt.XYZ()[0] > fMaxX || spt.XYZ()[1] < fMinY ||
            spt.XYZ()[1] > fMaxY || spt.XYZ()[2] < fMinZ || spt.XYZ()[2] > fMaxZ)
          continue;

        // Get hits associated with this SpacePoint.

        const art::PtrVector<recob::Hit>& spthits = fSptalgSep.getAssociatedHits(spt);

        // Fill separation histogram.

        double sep = fSptalgSep.separation(spthits);
        fHS->Fill(sep);
      }
    }

    // Loop over default space points and fill histograms.

    for (auto const& spt : spts3) {
      if (spt.XYZ()[0] < fMinX || spt.XYZ()[0] > fMaxX || spt.XYZ()[1] < fMinY ||
          spt.XYZ()[1] > fMaxY || spt.XYZ()[2] < fMinZ || spt.XYZ()[2] > fMaxZ)
        continue;

      fHchisq->Fill(spt.Chisq());
      fHx->Fill(spt.XYZ()[0]);
      fHy->Fill(spt.XYZ()[1]);
      fHz->Fill(spt.XYZ()[2]);

      // Get hits associated with this SpacePoint.

      std::vector<art::Ptr<recob::Hit>> spthits;
      const art::PtrVector<recob::Hit>& av_spthits = fSptalgDefault.getAssociatedHits(spt);
      for (auto const& ptr : av_spthits) {
        spthits.push_back(ptr);
      }

      // Fill single hit histograms.

      for (auto const& hitPtr : spthits) {
        const recob::Hit& hit = *hitPtr;

        geo::View_t view = hit.View();

        if (view == geo::kU) {
          fHAmpU->Fill(hit.PeakAmplitude());
          fHAreaU->Fill(hit.Integral());
          fHSumU->Fill(hit.ROISummedADC());
        }
        if (view == geo::kV) {
          fHAmpV->Fill(hit.PeakAmplitude());
          fHAreaV->Fill(hit.Integral());
          fHSumV->Fill(hit.ROISummedADC());
        }
        if (view == geo::kZ) {
          fHAmpW->Fill(hit.PeakAmplitude());
          fHAreaW->Fill(hit.Integral());
          fHSumW->Fill(hit.ROISummedADC());
        }
      }

      if (mc && fUseMC) {

        art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

        try {
          std::vector<double> mcxyz = bt_serv->SpacePointHitsToWeightedXYZ(clockData, spthits);
          fHMCdx->Fill(spt.XYZ()[0] - mcxyz[0]);
          fHMCdy->Fill(spt.XYZ()[1] - mcxyz[1]);
          fHMCdz->Fill(spt.XYZ()[2] - mcxyz[2]);
          if (spt.ErrXYZ()[0] > 0.)
            fHMCxpull->Fill((spt.XYZ()[0] - mcxyz[0]) / std::sqrt(spt.ErrXYZ()[0]));
          if (spt.ErrXYZ()[2] > 0.)
            fHMCypull->Fill((spt.XYZ()[1] - mcxyz[1]) / std::sqrt(spt.ErrXYZ()[2]));
          if (spt.ErrXYZ()[5] > 0.)
            fHMCzpull->Fill((spt.XYZ()[2] - mcxyz[2]) / std::sqrt(spt.ErrXYZ()[5]));
        }
        catch (cet::exception&) {
        }
      }
    }
  }
}

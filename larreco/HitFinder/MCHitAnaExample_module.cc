////////////////////////////////////////////////////////////////////////
// Class:       MCHitAnaExample
// Module Type: analyzer
// File:        MCHitAnaExample_module.cc
//
// Generated at Tue Jun 24 22:54:57 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "lardataobj/RecoBase/Hit.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TStopwatch.h>

namespace hit {

  class MCHitAnaExample;

  class MCHitAnaExample : public art::EDAnalyzer {
  public:
    explicit MCHitAnaExample(fhicl::ParameterSet const& p);
    virtual ~MCHitAnaExample();

  private:
    void analyze(art::Event const& e) override;

    //
    // Module names
    //

    /// Producer module for recob::Hit
    std::string fRecoHitModuleName;

    /// Producer module for MCHit
    std::string fMCHitModuleName;

    //
    // Stopwatch for timings
    //
    TStopwatch fAnaWatch, fSearchWatch, fReadWatch;

    //
    // Histogram for timing record
    //

    /// Time to read recob::Hit from disk
    TH1D* hRecoHitReadTime;

    /// Time to read MCHit from disk
    TH1D* hMCHitReadTime;

    /// Time to run analyze() function
    TH1D* hAnalysisTime;

    /// Time to search for corresponding MCHit per RecoHit
    TH1D* hMCHitSearchTime;

    /// Time to search for corresponding MCHit per RecoHit summed over all RecoHits (per event)
    TH1D* hMCHitSearchTimeSum;

    //
    // Histograms purely from MCHit
    //

    /// Histogram (per plane) for all MCHit charge
    std::vector<TH1D*> hMCHitQ_v;

    /// Histogram (per plane) for MCHit multiplicity in each event
    std::vector<TH1D*> hMCHitMult_v;

    //
    // Histograms purely from RecoHit
    //

    /// Histogram (per plane) for RecoHit charge
    std::vector<TH1D*> hRecoHitQ_v;

    /// Histogram (per plane) for RecoHit multiplicity in each event
    std::vector<TH1D*> hRecoHitMult_v;

    //
    // Histograms from MC+Reco combined/compared information
    //

    /// Histogram (per plane) for # of MCHit vs. # of RecoHit
    std::vector<TH2D*> hCorrMult_v;

    /// Histogram (per plane) for charge of some RecoHit that has no corresponding MCHit
    std::vector<TH1D*> hVoidHitQ_v;

    /// Histogram (per plane) for charge of one MCHit that is closest (in time) to one RecoHit
    std::vector<TH2D*> hCorrQ_v;

    /// Histogram (per plane) for charge sum of multiple MCHit found within start=>end time of one RecoHit
    std::vector<TH2D*> hCorrQSum_v;

    /// Histogram (per plane) for a ratio of RecoHit charge to the closest (in time) RecoHit charge
    std::vector<TH1D*> hQRatio_v;

    /**
       Histogram (per plane) for a ratio of charge sum of multiple MCHit to one RecoHit. Multiple
       MCHits are found within start=>end time of RecoHit.
    */
    std::vector<TH1D*> hQSumRatio_v;

    /// Histogram (per plane) for time-diff between one Reco hit and the closest (in time) RecoHit
    std::vector<TH1D*> hDT_v;

    /// Histogram (per plane) for MCHit multiplicity within start=>end time region of RecoHit
    std::vector<TH1D*> hMCHitMultPerRecoHit_v;
  };

  MCHitAnaExample::MCHitAnaExample(fhicl::ParameterSet const& p) : EDAnalyzer(p) // ,

  {
    fRecoHitModuleName = p.get<std::string>("RecoHitModuleName");
    fMCHitModuleName = p.get<std::string>("MCHitModuleName");

    art::ServiceHandle<art::TFileService const> fs;

    art::ServiceHandle<geo::Geometry const> geo;

    for (unsigned char plane = 0; plane < geo->Nplanes(); ++plane) {

      hMCHitQ_v.push_back(
        fs->make<TH1D>(Form("hMCHitQ_%d", plane),
                       Form("MCHit Charge (Plane %d); Charge [ADC]; # MCHit", plane),
                       250,
                       0,
                       2000));

      hMCHitMult_v.push_back(
        fs->make<TH1D>(Form("hMCHitMult_%d", plane),
                       Form("MCHit Multiplicity (Plane %d); Multiplicity; # MCHit", plane),
                       250,
                       0,
                       5000));

      hRecoHitQ_v.push_back(
        fs->make<TH1D>(Form("hRecoHitQ_%d", plane),
                       Form("RecoHit Charge (Plane %d); Charge [ADC]; # RecoHits", plane),
                       250,
                       0,
                       2000));

      hRecoHitMult_v.push_back(
        fs->make<TH1D>(Form("hRecoHitMult_%d", plane),
                       Form("RecoHit Multiplicity (Plane %d); Multiplicity; # RecoHits", plane),
                       200,
                       0,
                       2000));

      hCorrMult_v.push_back(
        fs->make<TH2D>(Form("hCorrMult_%d", plane),
                       Form("# MCHits vs. # RecoHits (Plane %d); # MCHits; # RecoHits", plane),
                       250,
                       0,
                       5000,
                       200,
                       0,
                       2000));

      hVoidHitQ_v.push_back(fs->make<TH1D>(
        Form("hVoidHitQ_%d", plane),
        Form("RecoHit Charge (No Corresponding MCHit)) (Plane %d); Charge [ADC]; # RecoHits",
             plane),
        250,
        0,
        2000));

      hCorrQ_v.push_back(fs->make<TH2D>(
        Form("hCorrQ_%d", plane),
        Form("RecoHit Q vs. Closest MCHit Q (Plane %d); MCHit Charge [ADC]; RecoHit Charge [ADC]",
             plane),
        250,
        0,
        2000,
        250,
        0,
        2000));

      hCorrQSum_v.push_back(
        fs->make<TH2D>(Form("hCorrQSum_%d", plane),
                       Form("RecoHit Q vs. MCHit QSum Within Start=>End  (Plane %d); MCHit Charge "
                            "[ADC]; RecoHit Charge [ADC]",
                            plane),
                       250,
                       0,
                       2000,
                       250,
                       0,
                       2000));

      hQRatio_v.push_back(
        fs->make<TH1D>(Form("hQRatio_%d", plane),
                       Form("Reco/MCHit Charge (Plane %d); Ratio; # RecoHits", plane),
                       200,
                       0,
                       2));

      hQSumRatio_v.push_back(
        fs->make<TH1D>(Form("hQSumRatio_%d", plane),
                       Form("Reco/MCHit Charge (Plane %d); Ratio; # RecoHits", plane),
                       200,
                       0,
                       2));

      hDT_v.push_back(fs->make<TH1D>(
        Form("hDT_%d", plane),
        Form("#DeltaT btw Reco and Closest MCHit (Plane %d); #DeltaT; # RecoHits", plane),
        190,
        -9.5,
        9.5));

      hMCHitMultPerRecoHit_v.push_back(fs->make<TH1D>(
        Form("hMCHitMultPerRecoHit_%d", plane),
        Form("MCHit Multiplicity per RecoHit (Plane %d); Multiplicity; # RecoHits", plane),
        21,
        -1.5,
        19.5));
    }

    hRecoHitReadTime =
      fs->make<TH1D>("hRecoHitReadTime",
                     "Real Time to Read recob::Hit From Disk; Real Time [ms]; Events",
                     200,
                     0,
                     1000);

    hMCHitReadTime =
      fs->make<TH1D>("hMCHitReadTime",
                     "Real Time to Read sim::MCHit From Disk; Real Time [ms]; Events",
                     200,
                     0,
                     1000);

    hAnalysisTime =
      fs->make<TH1D>("hAnalysisTime",
                     "Analysis Time to Run analyze() Function; Real Time [ms]; Events",
                     200,
                     0,
                     1000);

    hMCHitSearchTime =
      fs->make<TH1D>("hMCHitSearchTime",
                     "Time to Search Multiple sim::MCHit Per RecoHit; RealTime [us]; # RecoHit",
                     200,
                     0,
                     200);

    hMCHitSearchTimeSum = fs->make<TH1D>(
      "hMCHitSearchTimeSum",
      "Time to Search Multiple sim::MCHit for All RecoHit in Event; RealTime [ms]; # Event",
      200,
      0,
      200);
  }

  MCHitAnaExample::~MCHitAnaExample() {}

  void
  MCHitAnaExample::analyze(art::Event const& e)
  {
    fAnaWatch.Start();

    art::ServiceHandle<geo::Geometry const> geo;
    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

    fReadWatch.Start();
    auto const& mchits_v = *e.getValidHandle<std::vector<sim::MCHitCollection>>(fMCHitModuleName);
    hMCHitReadTime->Fill(fReadWatch.RealTime() * 1.e3);

    fReadWatch.Start();
    auto const& recohits = *e.getValidHandle<std::vector<recob::Hit>>(fRecoHitModuleName);
    hRecoHitReadTime->Fill(fReadWatch.RealTime() * 1.e3);

    //
    // Work on purely MCHit info
    //
    std::vector<size_t> mchit_mult(geo->Nplanes(), 0);
    for (auto const& mchits : mchits_v) {

      auto plane = geo->ChannelToWire(mchits.Channel()).at(0).Plane;

      mchit_mult.at(plane) += mchits.size();

      for (auto const& mchit : mchits)

        hMCHitQ_v.at(plane)->Fill(mchit.Charge(true));
    }

    std::vector<size_t> recohit_mult(geo->Nplanes(), 0);

    double search_time_sum = 0;

    // Loop over RecoHit
    for (size_t hit_index = 0; hit_index < recohits.size(); ++hit_index) {

      auto const& hit = recohits.at(hit_index);

      auto const& wire_id = hit.WireID();

      // Fill purly RecoHit info
      hRecoHitQ_v.at(wire_id.Plane)->Fill(hit.PeakAmplitude());
      recohit_mult.at(wire_id.Plane) += 1;

      // Figure out channel & retrieve MCHitCollection for this channel
      auto ch = geo->PlaneWireToChannel(wire_id.Plane, wire_id.Wire, wire_id.TPC, wire_id.Cryostat);

      if (mchits_v.size() <= ch)
        throw cet::exception(__PRETTY_FUNCTION__)
          << "Channel " << ch << " not found in MCHit vector!" << std::endl;

      // Found MCHitCollection
      auto const& mchits = mchits_v.at(ch);

      if (ch != mchits.Channel())
        throw cet::exception(__PRETTY_FUNCTION__)
          << "MCHit channel & vector index mismatch: " << mchits.Channel() << " != " << ch
          << std::endl;

      // Locate corresponding MCHit(s) to this RecoHit
      sim::MCHit start_time, end_time, peak_time;

      start_time.SetTime(clock_data.TPCTick2TDC(hit.PeakTimeMinusRMS()), 0);
      peak_time.SetTime(clock_data.TPCTick2TDC(hit.PeakTime()), 0);
      end_time.SetTime(clock_data.TPCTick2TDC(hit.PeakTimePlusRMS()), 0);

      fSearchWatch.Start();
      auto start_iter = std::lower_bound(mchits.begin(), mchits.end(), start_time);
      auto end_iter = std::upper_bound(mchits.begin(), mchits.end(), end_time);
      search_time_sum += fSearchWatch.RealTime() * 1.e3;
      hMCHitSearchTime->Fill(fSearchWatch.RealTime() * 1.e6);

      double reco_q = hit.PeakAmplitude();
      double mc_qsum = 0;
      double mc_q = 0;
      double mult = 0;
      double dt_min = 0;
      double abs_dt_min = 1e9;

      // Loop over MCHit(s) that reside in start=>end time of MCHit
      while (start_iter != end_iter) {
        mc_qsum += (*start_iter).Charge(true);
        ++mult;
        double dt = (*start_iter).PeakTime() - peak_time.PeakTime();
        double abs_dt = dt;
        if (abs_dt < 0) abs_dt *= -1;

        if (abs_dt < abs_dt_min) {
          abs_dt_min = abs_dt;
          dt_min = dt;
          mc_q = (*start_iter).Charge(true);
        }
        ++start_iter;
      }

      // Fill histograms
      hMCHitMultPerRecoHit_v.at(wire_id.Plane)->Fill(mult);

      if (!mc_qsum)
        hVoidHitQ_v.at(wire_id.Plane)->Fill(reco_q);

      else {

        hRecoHitQ_v.at(wire_id.Plane)->Fill(reco_q);
        hCorrQ_v.at(wire_id.Plane)->Fill(mc_q, reco_q);
        hCorrQSum_v.at(wire_id.Plane)->Fill(mc_qsum, reco_q);
        hQSumRatio_v.at(wire_id.Plane)->Fill(reco_q / mc_qsum);
        hQRatio_v.at(wire_id.Plane)->Fill(reco_q / mc_q);
        hDT_v.at(wire_id.Plane)->Fill(dt_min);
      }

    } // end looping over hits

    // Fill purely-MCHit and purely-RecoHit multiplicity histograms
    for (unsigned char plane = 0; plane < geo->Nplanes(); ++plane) {
      std::cout << mchit_mult.at(plane) << std::endl;
      hMCHitMult_v.at(plane)->Fill(mchit_mult.at(plane));
      hRecoHitMult_v.at(plane)->Fill(recohit_mult.at(plane));
      hCorrMult_v.at(plane)->Fill(mchit_mult.at(plane), recohit_mult.at(plane));
    }

    hAnalysisTime->Fill(fAnaWatch.RealTime() * 1.e3);
    hMCHitSearchTimeSum->Fill(search_time_sum);
  }
}
DEFINE_ART_MODULE(hit::MCHitAnaExample)

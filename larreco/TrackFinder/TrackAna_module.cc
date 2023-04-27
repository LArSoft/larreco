//
// Name: TrackAna_module.cc
//
// Purpose: Module TrackAna.
//
// Configuration parameters.
//
//  TrackModuleLabel:   Track module label.
//  MCTrackModuleLabel: MCTrack module label.
//  MinMCKE:            Minimum MC particle kinetic energy.
//  MatchColinearity:   Minimum colinearity for mc-track matching.
//  MatchDisp:          Maximum uv displacement for mc-track matching.
//  WMatchDisp:         maximum w displacement for mc-track matching.
//  MatchLength:        Minimum length fraction for mc-track matching.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/Simulation/sim.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixD.h"

namespace {

  // Calculate distance to boundary.
  //----------------------------------------------------------------------------
  double bdist(const TVector3& pos, unsigned int /*tpc*/ = 0, unsigned int /*cstat*/ = 0)
  {
    // Get geometry.

    art::ServiceHandle<geo::Geometry const> geom;

    double d1 = pos.X();                             // Distance to right side (wires).
    double d2 = 2. * geom->DetHalfWidth() - pos.X(); // Distance to left side (cathode).
    double d3 = pos.Y() + geom->DetHalfHeight();     // Distance to bottom.
    double d4 = geom->DetHalfHeight() - pos.Y();     // Distance to top.
    double d5 = pos.Z();                             // Distance to front.
    double d6 = geom->DetLength() - pos.Z();         // Distance to back.

    return std::min({d1, d2, d3, d4, d5, d6});
  }

  // Length of reconstructed track.
  //----------------------------------------------------------------------------
  double length(const recob::Track& track) { return track.Length(); }

  // Length of MC particle.
  //----------------------------------------------------------------------------
  double length(detinfo::DetectorClocksData const& clockData,
                detinfo::DetectorPropertiesData const& detProp,
                const simb::MCParticle& part,
                double dx,
                geo::Point_t& start,
                geo::Point_t& end,
                TVector3& startmom,
                TVector3& endmom)
  {
    art::ServiceHandle<geo::Geometry const> geom;

    // Get fiducial volume boundary.
    double xmin = 0.;
    double xmax = 2. * geom->DetHalfWidth();
    double ymin = -geom->DetHalfHeight();
    double ymax = geom->DetHalfHeight();
    double zmin = 0.;
    double zmax = geom->DetLength();
    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;

    for (int i = 0; i < n; ++i) {
      TVector3 pos = part.Position(i).Vect();

      // Make fiducial cuts.  Require the particle to be within the physical volume of
      // the tpc, and also require the apparent x position to be within the expanded
      // readout frame.

      if (pos.X() >= xmin && pos.X() <= xmax && pos.Y() >= ymin && pos.Y() <= ymax &&
          pos.Z() >= zmin && pos.Z() <= zmax) {
        pos[0] += dx;
        double ticks = detProp.ConvertXToTicks(pos[0], 0, 0, 0);
        if (ticks >= 0. && ticks < detProp.ReadOutWindowSize()) {
          if (first) {
            start = pos;
            startmom = part.Momentum(i).Vect();
          }
          else {
            disp -= pos;
            result += disp.Mag();
          }
          first = false;
          disp = pos;
          end = pos;
          endmom = part.Momentum(i).Vect();
        }
      }
    }

    return result;
  }

  // Length of MCTrack.
  // In this function, the extracted start and end momenta are converted to GeV
  // (MCTrack stores momenta in Mev).
  //----------------------------------------------------------------------------
  double length(detinfo::DetectorClocksData const& clockData,
                detinfo::DetectorPropertiesData const& detProp,
                const sim::MCTrack& mctrk,
                double dx,
                geo::Point_t& start,
                geo::Point_t& end,
                TVector3& startmom,
                TVector3& endmom)
  {
    art::ServiceHandle<geo::Geometry const> geom;

    // Get fiducial volume boundary.

    double xmin = 0.;
    double xmax = 2. * geom->DetHalfWidth();
    double ymin = -geom->DetHalfHeight();
    double ymax = geom->DetHalfHeight();
    double zmin = 0.;
    double zmax = geom->DetLength();
    double result = 0.;
    TVector3 disp;
    int n = mctrk.size();
    bool first = true;

    for (int i = 0; i < n; ++i) {
      TVector3 pos = mctrk[i].Position().Vect();

      // Make fiducial cuts.  Require the particle to be within the physical volume of
      // the tpc, and also require the apparent x position to be within the expanded
      // readout frame.

      if (pos.X() >= xmin && pos.X() <= xmax && pos.Y() >= ymin && pos.Y() <= ymax &&
          pos.Z() >= zmin && pos.Z() <= zmax) {
        pos[0] += dx;
        double ticks = detProp.ConvertXToTicks(pos[0], 0, 0, 0);
        if (ticks >= 0. && ticks < detProp.ReadOutWindowSize()) {
          if (first) {
            start = pos;
            startmom = 0.001 * mctrk[i].Momentum().Vect();
          }
          else {
            disp -= pos;
            result += disp.Mag();
          }
          first = false;
          disp = pos;
          end = pos;
          endmom = 0.001 * mctrk[i].Momentum().Vect();
        }
      }
    }

    return result;
  }

  // Fill efficiency histogram assuming binomial errors.

  void effcalc(const TH1* hnum, const TH1* hden, TH1* heff)
  {
    int nbins = hnum->GetNbinsX();
    if (nbins != hden->GetNbinsX())
      throw cet::exception("TrackAna") << "effcalc[" __FILE__ "]: incompatible histograms (I)\n";
    if (nbins != heff->GetNbinsX())
      throw cet::exception("TrackAna") << "effcalc[" __FILE__ "]: incompatible histograms (II)\n";

    // Loop over bins, including underflow and overflow.

    for (int ibin = 0; ibin <= nbins + 1; ++ibin) {
      double num = hnum->GetBinContent(ibin);
      double den = hden->GetBinContent(ibin);
      if (den == 0.) {
        heff->SetBinContent(ibin, 0.);
        heff->SetBinError(ibin, 0.);
      }
      else {
        double eff = num / den;
        if (eff < 0.) eff = 0.;
        if (eff > 1.) eff = 1.;
        double err = std::sqrt(eff * (1. - eff) / den);
        heff->SetBinContent(ibin, eff);
        heff->SetBinError(ibin, err);
      }
    }

    heff->SetMinimum(0.);
    heff->SetMaximum(1.05);
    heff->SetMarkerStyle(20);
  }

  class flattener : public std::vector<unsigned int> {
  public:
    flattener() = default;

    flattener(const std::vector<std::vector<unsigned int>>& input) { Convert(input); }

    void Convert(const std::vector<std::vector<unsigned int>>& input)
    {
      clear();
      size_t length = 0;
      for (auto const& vec : input)
        length += vec.size();
      reserve(length);

      for (auto const& vec : input)
        for (auto const& value : vec)
          push_back(value);
    }
  }; // end class flattener

} // end namespace

namespace trkf {

  class TrackAna : public art::EDAnalyzer {
  public:
    // Embedded structs.

    // Struct for histograms that depend on reco track only.

    struct RecoHists {
      explicit RecoHists(const std::string& subdir);

      // Pure reco track histograms.

      TH1F* fHstartx{nullptr};   // Starting x position.
      TH1F* fHstarty{nullptr};   // Starting y position.
      TH1F* fHstartz{nullptr};   // Starting z position.
      TH1F* fHstartd{nullptr};   // Starting distance to boundary.
      TH1F* fHendx{nullptr};     // Ending x position.
      TH1F* fHendy{nullptr};     // Ending y position.
      TH1F* fHendz{nullptr};     // Ending z position.
      TH1F* fHendd{nullptr};     // Ending distance to boundary.
      TH1F* fHtheta{nullptr};    // Theta.
      TH1F* fHphi{nullptr};      // Phi.
      TH1F* fHtheta_xz{nullptr}; // Theta_xz.
      TH1F* fHtheta_yz{nullptr}; // Theta_yz.
      TH1F* fHmom{nullptr};      // Momentum.
      TH1F* fHmoml{nullptr};     // Momentum (low momentum).
      TH1F* fHlen{nullptr};      // Length.
      TH1F* fHlens{nullptr};     // Length (short tracks).

      // Histograms for the consituent Hits

      TH1F* fHHitChg{nullptr};    // hit charge
      TH1F* fHHitWidth{nullptr};  // hit width
      TH1F* fHHitPdg{nullptr};    // Pdg primarily responsible.
      TH1F* fHHitTrkId{nullptr};  // TrkId
      TH1F* fModeFrac{nullptr};   // mode fraction
      TH1F* fNTrkIdTrks{nullptr}; // # of stitched tracks in which unique TrkId appears
      TH2F* fNTrkIdTrks2{nullptr};
      TH2F* fNTrkIdTrks3{nullptr};
    };

    // Struct for mc particles and mc-matched tracks.

    struct MCHists {
      explicit MCHists(const std::string& subdir);

      // Reco-MC matching.

      TH2F* fHduvcosth{nullptr};  // 2D mc vs. data matching, duv vs. cos(theta).
      TH1F* fHcosth{nullptr};     // 1D direction matching, cos(theta).
      TH1F* fHmcu{nullptr};       // 1D endpoint truth u.
      TH1F* fHmcv{nullptr};       // 1D endpoint truth v.
      TH1F* fHmcw{nullptr};       // 1D endpoint truth w.
      TH1F* fHupull{nullptr};     // 1D endpoint u pull.
      TH1F* fHvpull{nullptr};     // 1D endpoint v pull.
      TH1F* fHmcdudw{nullptr};    // Truth du/dw.
      TH1F* fHmcdvdw{nullptr};    // Truth dv/dw.
      TH1F* fHdudwpull{nullptr};  // du/dw pull.
      TH1F* fHdvdwpull{nullptr};  // dv/dw pull.
      TH1F* fHHitEff{nullptr};    // Hit efficiency.
      TH1F* fHHitPurity{nullptr}; // Hit purity.

      // Histograms for matched tracks.

      TH1F* fHstartdx{nullptr}; // Start dx.
      TH1F* fHstartdy{nullptr}; // Start dy.
      TH1F* fHstartdz{nullptr}; // Start dz.
      TH1F* fHenddx{nullptr};   // End dx.
      TH1F* fHenddy{nullptr};   // End dy.
      TH1F* fHenddz{nullptr};   // End dz.
      TH2F* fHlvsl{nullptr};    // MC vs. reco length.
      TH1F* fHdl{nullptr};      // Delta(length).
      TH2F* fHpvsp{nullptr};    // MC vs. reco momentum.
      TH2F* fHpvspc{nullptr};   // MC vs. reco momentum (contained tracks).
      TH1F* fHdp{nullptr};      // Momentum difference.
      TH1F* fHdpc{nullptr};     // Momentum difference (contained tracks).
      TH1F* fHppull{nullptr};   // Momentum pull.
      TH1F* fHppullc{nullptr};  // Momentum pull (contained tracks).

      // Pure MC particle histograms (efficiency denominator).

      TH1F* fHmcstartx{nullptr};   // Starting x position.
      TH1F* fHmcstarty{nullptr};   // Starting y position.
      TH1F* fHmcstartz{nullptr};   // Starting z position.
      TH1F* fHmcendx{nullptr};     // Ending x position.
      TH1F* fHmcendy{nullptr};     // Ending y position.
      TH1F* fHmcendz{nullptr};     // Ending z position.
      TH1F* fHmctheta{nullptr};    // Theta.
      TH1F* fHmcphi{nullptr};      // Phi.
      TH1F* fHmctheta_xz{nullptr}; // Theta_xz.
      TH1F* fHmctheta_yz{nullptr}; // Theta_yz.
      TH1F* fHmcmom{nullptr};      // Momentum.
      TH1F* fHmcmoml{nullptr};     // Momentum (low momentum).
      TH1F* fHmcke{nullptr};       // Kinetic energy.
      TH1F* fHmckel{nullptr};      // Kinetic energy (low energy).
      TH1F* fHmclen{nullptr};      // Length.
      TH1F* fHmclens{nullptr};     // Length (short tracks).

      // Histograms for well-reconstructed matched tracks (efficiency numerator).

      TH1F* fHgstartx{nullptr};   // Starting x position.
      TH1F* fHgstarty{nullptr};   // Starting y position.
      TH1F* fHgstartz{nullptr};   // Starting z position.
      TH1F* fHgendx{nullptr};     // Ending x position.
      TH1F* fHgendy{nullptr};     // Ending y position.
      TH1F* fHgendz{nullptr};     // Ending z position.
      TH1F* fHgtheta{nullptr};    // Theta.
      TH1F* fHgphi{nullptr};      // Phi.
      TH1F* fHgtheta_xz{nullptr}; // Theta_xz.
      TH1F* fHgtheta_yz{nullptr}; // Theta_yz.
      TH1F* fHgmom{nullptr};      // Momentum.
      TH1F* fHgmoml{nullptr};     // Momentum (low momentum).
      TH1F* fHgke{nullptr};       // Kinetic energy.
      TH1F* fHgkel{nullptr};      // Kinetic energy (low momentum).
      TH1F* fHglen{nullptr};      // Length.
      TH1F* fHglens{nullptr};     // Length (short tracks).

      // Efficiency histograms.

      TH1F* fHestartx{nullptr};   // Starting x position.
      TH1F* fHestarty{nullptr};   // Starting y position.
      TH1F* fHestartz{nullptr};   // Starting z position.
      TH1F* fHeendx{nullptr};     // Ending x position.
      TH1F* fHeendy{nullptr};     // Ending y position.
      TH1F* fHeendz{nullptr};     // Ending z position.
      TH1F* fHetheta{nullptr};    // Theta.
      TH1F* fHephi{nullptr};      // Phi.
      TH1F* fHetheta_xz{nullptr}; // Theta_xz.
      TH1F* fHetheta_yz{nullptr}; // Theta_yz.
      TH1F* fHemom{nullptr};      // Momentum.
      TH1F* fHemoml{nullptr};     // Momentum (low momentum).
      TH1F* fHeke{nullptr};       // Kinetic energy.
      TH1F* fHekel{nullptr};      // Kinetic energy (low momentum).
      TH1F* fHelen{nullptr};      // Length.
      TH1F* fHelens{nullptr};     // Length (short tracks).
    };

    // Constructors, destructor

    explicit TrackAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;
    void endJob() override;

    void anaStitch(const art::Event& evt,
                   detinfo::DetectorClocksData const& clockData,
                   detinfo::DetectorPropertiesData const& detProp);

    std::vector<size_t> fsort_indexes(const std::vector<double>& v);

    // Fcl Attributes.

    std::string fTrackModuleLabel;
    std::string fMCTrackModuleLabel;
    std::string fSpacepointModuleLabel;
    std::string fStitchModuleLabel;
    std::string fTrkSpptAssocModuleLabel;
    std::string fHitSpptAssocModuleLabel;
    std::string fHitModuleLabel;

    int fDump;                // Number of events to dump to debug message facility.
    double fMinMCKE;          // Minimum MC particle kinetic energy (GeV).
    double fMinMCLen;         // Minimum MC particle length in tpc (cm).
    double fMatchColinearity; // Minimum matching colinearity.
    double fMatchDisp;        // Maximum matching displacement.
    double fWMatchDisp;       // Maximum matching displacement in the w direction.
    double fMatchLength;      // Minimum length fraction.
    bool fIgnoreSign;         // Ignore sign of mc particle if true.
    bool fStitchedAnalysis;   // if true, do the whole drill-down from stitched track to assd hits

    std::string fOrigin;
    bool fCheckOrigin;
    simb::Origin_t fOriginValue;
    int fPrintLevel; // 0 = none, 1 = event summary, 2 = track detail

    // Histograms.

    std::map<int, MCHists> fMCHistMap;     // Indexed by pdg id.
    std::map<int, RecoHists> fRecoHistMap; // Indexed by pdg id.

    // Statistics.

    int fNumEvent;
  };

  DEFINE_ART_MODULE(TrackAna)

  // RecoHists methods.

  TrackAna::RecoHists::RecoHists(const std::string& subdir)
  {
    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHstartx = dir.make<TH1F>(
      "xstart", "X Start Position", 100, -2. * geom->DetHalfWidth(), 4. * geom->DetHalfWidth());
    fHstarty = dir.make<TH1F>(
      "ystart", "Y Start Position", 100, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHstartz = dir.make<TH1F>("zstart", "Z Start Position", 100, 0., geom->DetLength());
    fHstartd = dir.make<TH1F>(
      "dstart", "Start Position Distance to Boundary", 100, -10., geom->DetHalfWidth());
    fHendx = dir.make<TH1F>(
      "xend", "X End Position", 100, -2. * geom->DetHalfWidth(), 4. * geom->DetHalfWidth());
    fHendy =
      dir.make<TH1F>("yend", "Y End Position", 100, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHendz = dir.make<TH1F>("zend", "Z End Position", 100, 0., geom->DetLength());
    fHendd =
      dir.make<TH1F>("dend", "End Position Distance to Boundary", 100, -10., geom->DetHalfWidth());
    fHtheta = dir.make<TH1F>("theta", "Theta", 100, 0., 3.142);
    fHphi = dir.make<TH1F>("phi", "Phi", 100, -3.142, 3.142);
    fHtheta_xz = dir.make<TH1F>("theta_xz", "Theta_xz", 100, -3.142, 3.142);
    fHtheta_yz = dir.make<TH1F>("theta_yz", "Theta_yz", 100, -3.142, 3.142);
    fHmom = dir.make<TH1F>("mom", "Momentum", 100, 0., 10.);
    fHmoml = dir.make<TH1F>("moml", "Momentum", 100, 0., 1.);
    fHlen = dir.make<TH1F>("len", "Track Length", 100, 0., 1.1 * geom->DetLength());
    fHlens = dir.make<TH1F>("lens", "Track Length", 100, 0., 0.1 * geom->DetLength());
    fHHitChg = dir.make<TH1F>("hchg", "Hit Charge (ADC counts)", 100, 0., 4000.);
    fHHitWidth = dir.make<TH1F>("hwid", "Hit Width (ticks)", 40, 0., 20.);
    fHHitPdg = dir.make<TH1F>("hpdg", "Hit Pdg code", 5001, -2500.5, +2500.5);
    fHHitTrkId = dir.make<TH1F>("htrkid", "Hit Track ID", 401, -200.5, +200.5);
    fModeFrac =
      dir.make<TH1F>("hmodefrac",
                     "quasi-Purity: Fraction of component tracks with the Track mode value",
                     20,
                     0.01,
                     1.01);
    fNTrkIdTrks =
      dir.make<TH1F>("hntrkids",
                     "quasi-Efficiency: Number of stitched tracks in which TrkId appears",
                     20,
                     0.,
                     +10.0);
    fNTrkIdTrks2 = dir.make<TH2F>("hntrkids2",
                                  "Number of stitched tracks in which TrkId appears vs KE [GeV]",
                                  20,
                                  0.,
                                  +10.0,
                                  20,
                                  0.0,
                                  1.5);
    fNTrkIdTrks3 = dir.make<TH2F>("hntrkids3",
                                  "MC Track vs Reco Track, wtd by nhits on Collection Plane",
                                  10,
                                  -0.5,
                                  9.5,
                                  10,
                                  -0.5,
                                  9.5);
    fNTrkIdTrks3->GetXaxis()->SetTitle("Sorted-by-Descending-CPlane-Hits outer Track Number");
    fNTrkIdTrks3->GetYaxis()->SetTitle("Sorted-by-Descending-True-Length G4Track");
  }

  // MCHists methods.

  TrackAna::MCHists::MCHists(const std::string& subdir)
  {
    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHduvcosth =
      dir.make<TH2F>("duvcosth", "Delta(uv) vs. Colinearity", 100, 0.95, 1., 100, 0., 1.);
    fHcosth = dir.make<TH1F>("colin", "Colinearity", 100, 0.95, 1.);
    fHmcu = dir.make<TH1F>("mcu", "MC Truth U", 100, -5., 5.);
    fHmcv = dir.make<TH1F>("mcv", "MC Truth V", 100, -5., 5.);
    fHmcw = dir.make<TH1F>("mcw", "MC Truth W", 100, -20., 20.);
    fHupull = dir.make<TH1F>("dupull", "U Pull", 100, -20., 20.);
    fHvpull = dir.make<TH1F>("dvpull", "V Pull", 100, -20., 20.);
    fHmcdudw = dir.make<TH1F>("mcdudw", "MC Truth U Slope", 100, -0.2, 0.2);
    fHmcdvdw = dir.make<TH1F>("mcdvdw", "MV Truth V Slope", 100, -0.2, 0.2);
    fHdudwpull = dir.make<TH1F>("dudwpull", "U Slope Pull", 100, -10., 10.);
    fHdvdwpull = dir.make<TH1F>("dvdwpull", "V Slope Pull", 100, -10., 10.);
    fHHitEff = dir.make<TH1F>("hiteff", "MC Hit Efficiency", 100, 0., 1.0001);
    fHHitPurity = dir.make<TH1F>("hitpurity", "MC Hit Purity", 100, 0., 1.0001);
    fHstartdx = dir.make<TH1F>("startdx", "Start Delta x", 100, -10., 10.);
    fHstartdy = dir.make<TH1F>("startdy", "Start Delta y", 100, -10., 10.);
    fHstartdz = dir.make<TH1F>("startdz", "Start Delta z", 100, -10., 10.);
    fHenddx = dir.make<TH1F>("enddx", "End Delta x", 100, -10., 10.);
    fHenddy = dir.make<TH1F>("enddy", "End Delta y", 100, -10., 10.);
    fHenddz = dir.make<TH1F>("enddz", "End Delta z", 100, -10., 10.);
    fHlvsl = dir.make<TH2F>("lvsl",
                            "Reco Length vs. MC Truth Length",
                            100,
                            0.,
                            1.1 * geom->DetLength(),
                            100,
                            0.,
                            1.1 * geom->DetLength());
    fHdl = dir.make<TH1F>("dl", "Track Length Minus MC Particle Length", 100, -50., 50.);
    fHpvsp =
      dir.make<TH2F>("pvsp", "Reco Momentum vs. MC Truth Momentum", 100, 0., 5., 100, 0., 5.);
    fHpvspc = dir.make<TH2F>(
      "pvspc", "Reco Momentum vs. MC Truth Momentum (Contained Tracks)", 100, 0., 5., 100, 0., 5.);
    fHdp = dir.make<TH1F>("dp", "Reco-MC Momentum Difference", 100, -5., 5.);
    fHdpc = dir.make<TH1F>("dpc", "Reco-MC Momentum Difference (Contained Tracks)", 100, -5., 5.);
    fHppull = dir.make<TH1F>("ppull", "Momentum Pull", 100, -10., 10.);
    fHppullc = dir.make<TH1F>("ppullc", "Momentum Pull (Contained Tracks)", 100, -10., 10.);

    fHmcstartx = dir.make<TH1F>(
      "mcxstart", "MC X Start Position", 10, -2. * geom->DetHalfWidth(), 4. * geom->DetHalfWidth());
    fHmcstarty = dir.make<TH1F>(
      "mcystart", "MC Y Start Position", 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmcstartz = dir.make<TH1F>("mczstart", "MC Z Start Position", 10, 0., geom->DetLength());
    fHmcendx = dir.make<TH1F>(
      "mcxend", "MC X End Position", 10, -2. * geom->DetHalfWidth(), 4. * geom->DetHalfWidth());
    fHmcendy = dir.make<TH1F>(
      "mcyend", "MC Y End Position", 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmcendz = dir.make<TH1F>("mczend", "MC Z End Position", 10, 0., geom->DetLength());
    fHmctheta = dir.make<TH1F>("mctheta", "MC Theta", 20, 0., 3.142);
    fHmcphi = dir.make<TH1F>("mcphi", "MC Phi", 10, -3.142, 3.142);
    fHmctheta_xz = dir.make<TH1F>("mctheta_xz", "MC Theta_xz", 40, -3.142, 3.142);
    fHmctheta_yz = dir.make<TH1F>("mctheta_yz", "MC Theta_yz", 40, -3.142, 3.142);
    fHmcmom = dir.make<TH1F>("mcmom", "MC Momentum", 10, 0., 10.);
    fHmcmoml = dir.make<TH1F>("mcmoml", "MC Momentum", 10, 0., 1.);
    fHmcke = dir.make<TH1F>("mcke", "MC Kinetic Energy", 10, 0., 10.);
    fHmckel = dir.make<TH1F>("mckel", "MC Kinetic Energy", 10, 0., 1.);
    fHmclen = dir.make<TH1F>("mclen", "MC Particle Length", 10, 0., 1.1 * geom->DetLength());
    fHmclens = dir.make<TH1F>("mclens", "MC Particle Length", 10, 0., 0.1 * geom->DetLength());

    fHgstartx = dir.make<TH1F>("gxstart",
                               "Good X Start Position",
                               10,
                               -2. * geom->DetHalfWidth(),
                               4. * geom->DetHalfWidth());
    fHgstarty = dir.make<TH1F>(
      "gystart", "Good Y Start Position", 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHgstartz = dir.make<TH1F>("gzstart", "Good Z Start Position", 10, 0., geom->DetLength());
    fHgendx = dir.make<TH1F>(
      "gxend", "Good X End Position", 10, -2. * geom->DetHalfWidth(), 4. * geom->DetHalfWidth());
    fHgendy = dir.make<TH1F>(
      "gyend", "Good Y End Position", 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHgendz = dir.make<TH1F>("gzend", "Good Z End Position", 10, 0., geom->DetLength());
    fHgtheta = dir.make<TH1F>("gtheta", "Good Theta", 20, 0., 3.142);
    fHgphi = dir.make<TH1F>("gphi", "Good Phi", 10, -3.142, 3.142);
    fHgtheta_xz = dir.make<TH1F>("gtheta_xz", "Good Theta_xz", 40, -3.142, 3.142);
    fHgtheta_yz = dir.make<TH1F>("gtheta_yz", "Good Theta_yz", 40, -3.142, 3.142);
    fHgmom = dir.make<TH1F>("gmom", "Good Momentum", 10, 0., 10.);
    fHgmoml = dir.make<TH1F>("gmoml", "Good Momentum", 10, 0., 1.);
    fHgke = dir.make<TH1F>("gke", "Good Kinetic Energy", 10, 0., 10.);
    fHgkel = dir.make<TH1F>("gkel", "Good Kinetic Energy", 10, 0., 1.);
    fHglen = dir.make<TH1F>("glen", "Good Particle Length", 10, 0., 1.1 * geom->DetLength());
    fHglens = dir.make<TH1F>("glens", "Good Particle Length", 10, 0., 0.1 * geom->DetLength());

    fHestartx = dir.make<TH1F>("exstart",
                               "Efficiency vs. X Start Position",
                               10,
                               -2. * geom->DetHalfWidth(),
                               4. * geom->DetHalfWidth());
    fHestarty = dir.make<TH1F>("eystart",
                               "Efficiency vs. Y Start Position",
                               10,
                               -geom->DetHalfHeight(),
                               geom->DetHalfHeight());
    fHestartz =
      dir.make<TH1F>("ezstart", "Efficiency vs. Z Start Position", 10, 0., geom->DetLength());
    fHeendx = dir.make<TH1F>("exend",
                             "Efficiency vs. X End Position",
                             10,
                             -2. * geom->DetHalfWidth(),
                             4. * geom->DetHalfWidth());
    fHeendy = dir.make<TH1F>(
      "eyend", "Efficiency vs. Y End Position", 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHeendz = dir.make<TH1F>("ezend", "Efficiency vs. Z End Position", 10, 0., geom->DetLength());
    fHetheta = dir.make<TH1F>("etheta", "Efficiency vs. Theta", 20, 0., 3.142);
    fHephi = dir.make<TH1F>("ephi", "Efficiency vs. Phi", 10, -3.142, 3.142);
    fHetheta_xz = dir.make<TH1F>("etheta_xz", "Efficiency vs. Theta_xz", 40, -3.142, 3.142);
    fHetheta_yz = dir.make<TH1F>("etheta_yz", "Efficiency vs. Theta_yz", 40, -3.142, 3.142);
    fHemom = dir.make<TH1F>("emom", "Efficiency vs. Momentum", 10, 0., 10.);
    fHemoml = dir.make<TH1F>("emoml", "Efficiency vs. Momentum", 10, 0., 1.);
    fHeke = dir.make<TH1F>("eke", "Efficiency vs. Kinetic Energy", 10, 0., 10.);
    fHekel = dir.make<TH1F>("ekel", "Efficiency vs. Kinetic Energy", 10, 0., 1.);
    fHelen =
      dir.make<TH1F>("elen", "Efficiency vs. Particle Length", 10, 0., 1.1 * geom->DetLength());
    fHelens =
      dir.make<TH1F>("elens", "Efficiency vs. Particle Length", 10, 0., 0.1 * geom->DetLength());
  }

  TrackAna::TrackAna(const fhicl::ParameterSet& pset)
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    : EDAnalyzer(pset)
    , fTrackModuleLabel(pset.get<std::string>("TrackModuleLabel"))
    , fMCTrackModuleLabel(pset.get<std::string>("MCTrackModuleLabel"))
    , fSpacepointModuleLabel(pset.get<std::string>("SpacepointModuleLabel"))
    , fStitchModuleLabel(pset.get<std::string>("StitchModuleLabel"))
    , fTrkSpptAssocModuleLabel(pset.get<std::string>("TrkSpptAssocModuleLabel"))
    , fHitSpptAssocModuleLabel(pset.get<std::string>("HitSpptAssocModuleLabel"))
    , fHitModuleLabel(pset.get<std::string>("HitModuleLabel"))
    , fDump(pset.get<int>("Dump"))
    , fMinMCKE(pset.get<double>("MinMCKE"))
    , fMinMCLen(pset.get<double>("MinMCLen"))
    , fMatchColinearity(pset.get<double>("MatchColinearity"))
    , fMatchDisp(pset.get<double>("MatchDisp"))
    , fWMatchDisp(pset.get<double>("WMatchDisp"))
    , fMatchLength(pset.get<double>("MatchLength"))
    , fIgnoreSign(pset.get<bool>("IgnoreSign"))
    , fStitchedAnalysis(pset.get<bool>("StitchedAnalysis", false))
    , fOrigin(pset.get<std::string>("MCTrackOrigin", "Any"))
    , fPrintLevel(pset.get<int>("PrintLevel", 0))
    , fNumEvent(0)
  {

    // Decide whether to check MCTrack origin
    fCheckOrigin = false;
    fOriginValue = simb::kUnknown;
    if (fOrigin.find("Beam") != std::string::npos) {
      fCheckOrigin = true;
      fOriginValue = simb::kBeamNeutrino;
    }
    else if (fOrigin.find("Cosmic") != std::string::npos) {
      fCheckOrigin = true;
      fOriginValue = simb::kCosmicRay;
    }
    else if (fOrigin.find("Super") != std::string::npos) {
      fCheckOrigin = true;
      fOriginValue = simb::kSuperNovaNeutrino;
    }
    else if (fOrigin.find("Single") != std::string::npos) {
      fCheckOrigin = true;
      fOriginValue = simb::kSingleParticle;
    }

    // Report.

    mf::LogInfo("TrackAna") << "TrackAna configured with the following parameters:\n"
                            << "  TrackModuleLabel = " << fTrackModuleLabel << "\n"
                            << "  MCTrackModuleLabel = " << fMCTrackModuleLabel << "\n"
                            << "  StitchModuleLabel = " << fStitchModuleLabel << "\n"
                            << "  TrkSpptAssocModuleLabel = " << fTrkSpptAssocModuleLabel << "\n"
                            << "  HitSpptAssocModuleLabel = " << fHitSpptAssocModuleLabel << "\n"
                            << "  HitModuleLabel = " << fHitModuleLabel << "\n"
                            << "  Dump = " << fDump << "\n"
                            << "  MinMCKE = " << fMinMCKE << "\n"
                            << "  MinMCLen = " << fMinMCLen << "  Origin = " << fOrigin
                            << " Origin value " << fOriginValue;
  }

  void TrackAna::analyze(const art::Event& evt)
  {
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    art::ServiceHandle<geo::Geometry const> geom;

    ++fNumEvent;

    // Optional dump stream.

    std::unique_ptr<mf::LogInfo> pdump;
    if (fDump > 0) {
      --fDump;
      pdump = std::make_unique<mf::LogInfo>("TrackAna");
    }

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();

    // Get MCTracks.

    art::Handle<std::vector<sim::MCTrack>> mctrackh;
    evt.getByLabel(fMCTrackModuleLabel, mctrackh);
    // pair of MCTrack and index of matched reco track
    std::vector<std::pair<const sim::MCTrack*, int>> selected_mctracks;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    if (mc && mctrackh.isValid()) {

      selected_mctracks.reserve(mctrackh->size());

      // Dump MCTracks.

      if (pdump) {
        *pdump << "MC Tracks\n"
               << "       Id   pdg           x         y         z          dx        dy        dz "
                  "          p\n"
               << "--------------------------------------------------------------------------------"
                  "-----------\n";
      }

      // Loop over mc tracks, and fill histograms that depend only
      // on mc information.  Also, fill a secondary list of mc tracks
      // that pass various selection criteria.

      for (std::vector<sim::MCTrack>::const_iterator imctrk = mctrackh->begin();
           imctrk != mctrackh->end();
           ++imctrk) {
        const sim::MCTrack& mctrk = *imctrk;
        int pdg = mctrk.PdgCode();
        if (fIgnoreSign) pdg = std::abs(pdg);

        // Ignore everything except stable charged nonshowering particles.

        int apdg = std::abs(pdg);
        if (apdg == 13 ||   // Muon
            apdg == 211 ||  // Charged pion
            apdg == 321 ||  // Charged kaon
            apdg == 2212) { // (Anti)proton

          // check MC track origin?
          if (fCheckOrigin && mctrk.Origin() != fOriginValue) continue;

          // Apply minimum energy cut.

          if (mctrk.Start().E() >= mctrk.Start().Momentum().Mag() + 1000. * fMinMCKE) {

            // Calculate the x offset due to nonzero mc particle time.

            double mctime = mctrk.Start().T();                      // nsec
            double mcdx = mctime * 1.e-3 * detProp.DriftVelocity(); // cm

            // Calculate the length of this mc particle inside the fiducial volume.

            geo::Point_t mcstart;
            geo::Point_t mcend;
            TVector3 mcstartmom;
            TVector3 mcendmom;
            double plen =
              length(clockData, detProp, mctrk, mcdx, mcstart, mcend, mcstartmom, mcendmom);

            // Apply minimum fiducial length cut.  Always reject particles that have
            // zero length in the tpc regardless of the configured cut.

            if (plen > 0. && plen > fMinMCLen) {

              // This is a good mc particle (capable of making a track).

              selected_mctracks.push_back(std::make_pair(&mctrk, -1));

              // Dump MC particle information here.

              if (pdump) {
                double pstart = mcstartmom.Mag();
                double pend = mcendmom.Mag();
                *pdump << "\nOffset" << std::setw(3) << mctrk.TrackID() << std::setw(6)
                       << mctrk.PdgCode() << "  " << std::fixed << std::setprecision(2)
                       << std::setw(10) << mcdx << "\nStart " << std::setw(3) << mctrk.TrackID()
                       << std::setw(6) << mctrk.PdgCode() << "  " << std::fixed
                       << std::setprecision(2) << std::setw(10) << mcstart.X() << std::setw(10)
                       << mcstart.Y() << std::setw(10) << mcstart.Z();
                if (pstart > 0.) {
                  *pdump << "  " << std::fixed << std::setprecision(3) << std::setw(10)
                         << mcstartmom[0] / pstart << std::setw(10) << mcstartmom[1] / pstart
                         << std::setw(10) << mcstartmom[2] / pstart;
                }
                else
                  *pdump << std::setw(32) << " ";
                *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pstart;
                *pdump << "\nEnd   " << std::setw(3) << mctrk.TrackID() << std::setw(6)
                       << mctrk.PdgCode() << "  " << std::fixed << std::setprecision(2)
                       << std::setw(10) << mcend.X() << std::setw(10) << mcend.Y() << std::setw(10)
                       << mcend.Z();
                if (pend > 0.01) {
                  *pdump << "  " << std::fixed << std::setprecision(3) << std::setw(10)
                         << mcendmom[0] / pend << std::setw(10) << mcendmom[1] / pend
                         << std::setw(10) << mcendmom[2] / pend;
                }
                else
                  *pdump << std::setw(32) << " ";
                *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pend
                       << "\nLength: " << plen << "\n";
              }

              // Fill histograms.

              if (fMCHistMap.count(pdg) == 0) {
                std::ostringstream ostr;
                ostr << "MC" << (fIgnoreSign ? "All" : (pdg > 0 ? "Pos" : "Neg")) << std::abs(pdg);
                fMCHistMap.emplace(pdg, MCHists{ostr.str()});
              }
              const MCHists& mchists = fMCHistMap.at(pdg);

              double mctheta_xz = std::atan2(mcstartmom.X(), mcstartmom.Z());
              double mctheta_yz = std::atan2(mcstartmom.Y(), mcstartmom.Z());
              double mcmom = mcstartmom.Mag();
              double mcmass = 0.001 * mctrk.Start().Momentum().Mag();
              double mcke = mcmom * mcmom / (std::hypot(mcmom, mcmass) + mcmass);

              mchists.fHmcstartx->Fill(mcstart.X());
              mchists.fHmcstarty->Fill(mcstart.Y());
              mchists.fHmcstartz->Fill(mcstart.Z());
              mchists.fHmcendx->Fill(mcend.X());
              mchists.fHmcendy->Fill(mcend.Y());
              mchists.fHmcendz->Fill(mcend.Z());
              mchists.fHmctheta->Fill(mcstartmom.Theta());
              mchists.fHmcphi->Fill(mcstartmom.Phi());
              mchists.fHmctheta_xz->Fill(mctheta_xz);
              mchists.fHmctheta_yz->Fill(mctheta_yz);
              mchists.fHmcmom->Fill(mcmom);
              mchists.fHmcmoml->Fill(mcmom);
              mchists.fHmcke->Fill(mcke);
              mchists.fHmckel->Fill(mcke);
              mchists.fHmclen->Fill(plen);
              mchists.fHmclens->Fill(plen);
            }
          }
        }
      }
    } //mc

    // Get tracks and spacepoints and hits
    art::Handle<std::vector<recob::Track>> trackh;
    art::Handle<std::vector<art::PtrVector<recob::Track>>> trackvh;
    art::Handle<std::vector<recob::Hit>> hith;

    evt.getByLabel(fTrackModuleLabel, trackh);
    evt.getByLabel(fStitchModuleLabel, trackvh);
    evt.getByLabel(fHitModuleLabel, hith);

    // Extract all hits into a vector of art::Ptrs (the format used by back tracker).

    std::vector<art::Ptr<recob::Hit>> allhits;
    if (hith.isValid()) {
      allhits.reserve(hith->size());
      for (unsigned int i = 0; i < hith->size(); ++i) {
        allhits.emplace_back(hith, i);
      }
    }

    // Construct FindManyP object to be used for finding track-hit associations.

    art::FindManyP<recob::Hit> tkhit_find(trackh, evt, fTrackModuleLabel);

    // This new top part of TrackAna between two long lines of ************s
    // is particular to analyzing Stitched Tracks.
    // *******************************************************************//

    if (trackvh.isValid() && fStitchedAnalysis) {
      mf::LogDebug("TrackAna") << "TrackAna read " << trackvh->size()
                               << "  vectors of Stitched PtrVectorsof tracks.";
      anaStitch(evt, clockData, detProp);
    }

    if (trackh.isValid()) {

      if (pdump) {
        *pdump << "\nReconstructed Tracks\n"
               << "       Id  MCid           x         y         z          dx        dy        dz "
                  "          p\n"
               << "--------------------------------------------------------------------------------"
                  "-----------\n";
      }

      // Loop over tracks.

      int ntrack = trackh->size();
      for (int i = 0; i < ntrack; ++i) {
        art::Ptr<recob::Track> ptrack(trackh, i);
        const recob::Track& track = *ptrack;

        // Extract hits associated with this track.

        std::vector<art::Ptr<recob::Hit>> trackhits;
        tkhit_find.get(i, trackhits);

        // Calculate the x offset due to nonzero reconstructed time.

        double trackdx = 0;

        // Fill histograms involving reco tracks only.

        int ntraj = track.NumberTrajectoryPoints();
        if (ntraj > 0) {
          TVector3 pos = track.Vertex<TVector3>();
          TVector3 dir = track.VertexDirection<TVector3>();
          TVector3 end = track.End<TVector3>();
          pos[0] += trackdx;
          end[0] += trackdx;

          double dpos = bdist(pos);
          double dend = bdist(end);
          double tlen = length(track);
          double theta_xz = std::atan2(dir.X(), dir.Z());
          double theta_yz = std::atan2(dir.Y(), dir.Z());

          if (fRecoHistMap.count(0) == 0) fRecoHistMap.emplace(0, RecoHists{"Reco"});
          const RecoHists& rhists = fRecoHistMap.at(0);

          rhists.fHstartx->Fill(pos.X());
          rhists.fHstarty->Fill(pos.Y());
          rhists.fHstartz->Fill(pos.Z());
          rhists.fHstartd->Fill(dpos);
          rhists.fHendx->Fill(end.X());
          rhists.fHendy->Fill(end.Y());
          rhists.fHendz->Fill(end.Z());
          rhists.fHendd->Fill(dend);
          rhists.fHtheta->Fill(dir.Theta());
          rhists.fHphi->Fill(dir.Phi());
          rhists.fHtheta_xz->Fill(theta_xz);
          rhists.fHtheta_yz->Fill(theta_yz);

          double mom = 0.;
          if (track.HasMomentum()) mom = track.VertexMomentum();
          rhists.fHmom->Fill(mom);
          rhists.fHmoml->Fill(mom);
          rhists.fHlen->Fill(tlen);
          rhists.fHlens->Fill(tlen);

          // Id of matching mc particle.

          int mcid = -1;

          // Loop over direction.

          for (int swap = 0; swap < 2; ++swap) {

            // Analyze reversed tracks only if start momentum = end momentum.

            if (swap != 0 && track.HasMomentum() &&
                std::abs(track.VertexMomentum() - track.EndMomentum()) > 1.e-3)
              continue;

            // Calculate the global-to-local rotation matrix.

            int start_point = (swap == 0 ? 0 : ntraj - 1);
            TMatrixD rot = track.GlobalToLocalRotationAtPoint<TMatrixD>(start_point);

            // Update track data for reversed case.

            if (swap != 0) {
              rot(1, 0) = -rot(1, 0);
              rot(2, 0) = -rot(2, 0);
              rot(1, 1) = -rot(1, 1);
              rot(2, 1) = -rot(2, 1);
              rot(1, 2) = -rot(1, 2);
              rot(2, 2) = -rot(2, 2);

              pos = track.End<TVector3>();
              dir = -track.EndDirection<TVector3>();
              end = track.Vertex<TVector3>();
              pos[0] += trackdx;
              end[0] += trackdx;

              dpos = bdist(pos);
              dend = bdist(end);
              theta_xz = std::atan2(dir.X(), dir.Z());
              theta_yz = std::atan2(dir.Y(), dir.Z());

              if (track.HasMomentum()) mom = track.EndMomentum();
            }

            // Get covariance matrix.

            const auto& cov = (swap == 0 ? track.VertexCovariance() : track.EndCovariance());

            // Loop over track-like mc particles.

            for (unsigned int imc = 0; imc < selected_mctracks.size(); ++imc) {
              const sim::MCTrack& mctrk = *selected_mctracks[imc].first;
              int pdg = mctrk.PdgCode();
              if (fIgnoreSign) pdg = std::abs(pdg);
              auto iMCHistMap = fMCHistMap.find(pdg);
              if (iMCHistMap == fMCHistMap.end())
                throw cet::exception("TrackAna") << "no particle with ID=" << pdg << "\n";
              const MCHists& mchists = iMCHistMap->second;

              // Calculate the x offset due to nonzero mc particle time.

              double mctime = mctrk.Start().T();                      // nsec
              double mcdx = mctime * 1.e-3 * detProp.DriftVelocity(); // cm

              // Calculate the points where this mc particle enters and leaves the
              // fiducial volume, and the length in the fiducial volume.

              geo::Point_t mcstart;
              geo::Point_t mcend;
              TVector3 mcstartmom;
              TVector3 mcendmom;
              double plen =
                length(clockData, detProp, mctrk, mcdx, mcstart, mcend, mcstartmom, mcendmom);

              // Get the displacement of this mc particle in the global coordinate system.

              TVector3 mcpos{mcstart.X() - pos[0], mcstart.Y() - pos[1], mcstart.Z() - pos[2]};

              // Rotate the momentum and position to the
              // track-local coordinate system.

              TVector3 mcmoml = rot * mcstartmom;
              TVector3 mcposl = rot * mcpos;

              double colinearity = mcmoml.Z() / mcmoml.Mag();

              double u = mcposl.X();
              double v = mcposl.Y();
              double w = mcposl.Z();

              double pu = mcmoml.X();
              double pv = mcmoml.Y();
              double pw = mcmoml.Z();

              double dudw = pu / pw;
              double dvdw = pv / pw;

              double u0 = u - w * dudw;
              double v0 = v - w * dvdw;
              double uv0 = std::hypot(u0, v0);

              mchists.fHduvcosth->Fill(colinearity, uv0);

              if (std::abs(uv0) < fMatchDisp) {

                // Fill slope matching histograms.

                mchists.fHmcdudw->Fill(dudw);
                mchists.fHmcdvdw->Fill(dvdw);
                mchists.fHdudwpull->Fill(dudw / std::sqrt(cov(2, 2)));
                mchists.fHdvdwpull->Fill(dvdw / std::sqrt(cov(3, 3)));
              }
              mchists.fHcosth->Fill(colinearity);
              if (colinearity > fMatchColinearity) {

                // Fill displacement matching histograms.

                mchists.fHmcu->Fill(u0);
                mchists.fHmcv->Fill(v0);
                mchists.fHmcw->Fill(w);
                mchists.fHupull->Fill(u0 / std::sqrt(cov(0, 0)));
                mchists.fHvpull->Fill(v0 / std::sqrt(cov(1, 1)));

                if (std::abs(uv0) < fMatchDisp) {

                  // Fill matching histograms.

                  double mctheta_xz = std::atan2(mcstartmom.X(), mcstartmom.Z());
                  double mctheta_yz = std::atan2(mcstartmom.Y(), mcstartmom.Z());
                  double mcmom = mcstartmom.Mag();
                  double mcmass = 0.001 * mctrk.Start().Momentum().Mag();
                  double mcke = mcmom * mcmom / (std::hypot(mcmom, mcmass) + mcmass);

                  mchists.fHstartdx->Fill(pos.X() - mcstart.X());
                  mchists.fHstartdy->Fill(pos.Y() - mcstart.Y());
                  mchists.fHstartdz->Fill(pos.Z() - mcstart.Z());
                  mchists.fHenddx->Fill(end.X() - mcend.X());
                  mchists.fHenddy->Fill(end.Y() - mcend.Y());
                  mchists.fHenddz->Fill(end.Z() - mcend.Z());
                  mchists.fHlvsl->Fill(plen, tlen);
                  mchists.fHdl->Fill(tlen - plen);
                  mchists.fHpvsp->Fill(mcmom, mom);
                  double dp = mom - mcmom;
                  mchists.fHdp->Fill(dp);
                  mchists.fHppull->Fill(dp / std::sqrt(cov(4, 4)));
                  if (std::abs(dpos) >= 5. && std::abs(dend) >= 5.) {
                    mchists.fHpvspc->Fill(mcmom, mom);
                    mchists.fHdpc->Fill(dp);
                    mchists.fHppullc->Fill(dp / std::sqrt(cov(4, 4)));
                  }

                  // Count this track as well-reconstructed if it is matched to an
                  // mc particle (which it is if get here), and if in addition the
                  // starting position (w) matches and the reconstructed track length
                  // is more than 0.5 of the mc particle trajectory length.

                  bool good = std::abs(w) <= fWMatchDisp && tlen > fMatchLength * plen;
                  if (good) {
                    mcid = mctrk.TrackID();

                    // Calculate and fill hit efficiency and purity.

                    std::set<int> tkidset;
                    tkidset.insert(mcid);
                    double hiteff = bt_serv->HitCollectionEfficiency(
                      clockData, tkidset, trackhits, allhits, geo::k3D);
                    double hitpurity = bt_serv->HitCollectionPurity(clockData, tkidset, trackhits);
                    mchists.fHHitEff->Fill(hiteff);
                    mchists.fHHitPurity->Fill(hitpurity);

                    // Fill efficiency numerator histograms.

                    mchists.fHgstartx->Fill(mcstart.X());
                    mchists.fHgstarty->Fill(mcstart.Y());
                    mchists.fHgstartz->Fill(mcstart.Z());
                    mchists.fHgendx->Fill(mcend.X());
                    mchists.fHgendy->Fill(mcend.Y());
                    mchists.fHgendz->Fill(mcend.Z());
                    mchists.fHgtheta->Fill(mcstartmom.Theta());
                    mchists.fHgphi->Fill(mcstartmom.Phi());
                    mchists.fHgtheta_xz->Fill(mctheta_xz);
                    mchists.fHgtheta_yz->Fill(mctheta_yz);
                    mchists.fHgmom->Fill(mcmom);
                    mchists.fHgmoml->Fill(mcmom);
                    mchists.fHgke->Fill(mcke);
                    mchists.fHgkel->Fill(mcke);
                    mchists.fHglen->Fill(plen);
                    mchists.fHglens->Fill(plen);

                    // set the match flag
                    selected_mctracks[imc].second = i;

                    if (fPrintLevel > 0) {
                      const simb::MCParticle* ptkl = pi_serv->TrackIdToParticle_P(mcid);
                      float KE = ptkl->E() - ptkl->Mass();
                      std::string KEUnits = " GeV";
                      if (mctrk.Origin() != simb::kCosmicRay) {
                        // MeV for low energy particles
                        KE *= 1000;
                        KEUnits = " MeV";
                      }
                      mf::LogVerbatim("TrackAna")
                        << evt.run() << "." << evt.event() << " Match MCTkID " << std::setw(6)
                        << mctrk.TrackID() << " Origin " << mctrk.Origin() << " PDG" << std::setw(5)
                        << mctrk.PdgCode() << " KE" << std::setw(4) << (int)KE << KEUnits
                        << " RecoTrkID " << track.ID() << " hitEff " << std::setprecision(2)
                        << hiteff << " hitPur " << hitpurity;
                      // this won't work for DUNE
                      for (unsigned short ipl = 0; ipl < geom->Nplanes(); ++ipl) {
                        geo::PlaneID const planeID{0, 0, ipl};
                        auto const sWire = geom->NearestWireID(mcstart, planeID).Wire;
                        auto const sTick = detProp.ConvertXToTicks(mcstart.X(), planeID);
                        auto const eWire = geom->NearestWireID(mcend, planeID).Wire;
                        auto const eTick = detProp.ConvertXToTicks(mcend.X(), planeID);
                        mf::LogVerbatim("TrackAna")
                          << "   Wire:Tick in Pln " << ipl << " W:T " << sWire << ":" << sTick
                          << " - " << eWire << ":" << eTick;
                      } // ipl
                    }   // fPrintLevel == 2
                  }     // good
                }
              }
            }
          }

          // Dump track information here.

          if (pdump) {
            TVector3 pos = track.Vertex<TVector3>();
            TVector3 dir = track.VertexDirection<TVector3>();
            TVector3 end = track.End<TVector3>();
            pos[0] += trackdx;
            end[0] += trackdx;
            TVector3 enddir = track.EndDirection<TVector3>();
            double pstart = track.VertexMomentum();
            double pend = track.EndMomentum();
            *pdump << "\nOffset" << std::setw(3) << track.ID() << std::setw(6) << mcid << "  "
                   << std::fixed << std::setprecision(2) << std::setw(10) << trackdx << "\nStart "
                   << std::setw(3) << track.ID() << std::setw(6) << mcid << "  " << std::fixed
                   << std::setprecision(2) << std::setw(10) << pos[0] << std::setw(10) << pos[1]
                   << std::setw(10) << pos[2];
            if (pstart > 0.) {
              *pdump << "  " << std::fixed << std::setprecision(3) << std::setw(10) << dir[0]
                     << std::setw(10) << dir[1] << std::setw(10) << dir[2];
            }
            else
              *pdump << std::setw(32) << " ";
            *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pstart;
            *pdump << "\nEnd   " << std::setw(3) << track.ID() << std::setw(6) << mcid << "  "
                   << std::fixed << std::setprecision(2) << std::setw(10) << end[0] << std::setw(10)
                   << end[1] << std::setw(10) << end[2];
            if (pend > 0.01) {
              *pdump << "  " << std::fixed << std::setprecision(3) << std::setw(10) << enddir[0]
                     << std::setw(10) << enddir[1] << std::setw(10) << enddir[2];
            }
            else
              *pdump << std::setw(32) << " ";
            *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pend
                   << "\nLength: " << tlen << "\n";
          }
        }
      }
    } // i

    // print out un-matched MC tracks
    if (fPrintLevel > 0) {
      for (unsigned int imc = 0; imc < selected_mctracks.size(); ++imc) {
        if (selected_mctracks[imc].second >= 0) continue;
        const sim::MCTrack& mctrk = *selected_mctracks[imc].first;
        const simb::MCParticle* ptkl = pi_serv->TrackIdToParticle_P(mctrk.TrackID());
        float KE = ptkl->E() - ptkl->Mass();
        std::string KEUnits = " GeV";
        if (mctrk.Origin() != simb::kCosmicRay) {
          // MeV for low energy particles
          KE *= 1000;
          KEUnits = " MeV";
        }
        // find the start/end wire:time in each plane
        geo::Point_t mcstart, mcend;
        TVector3 mcstartmom, mcendmom;
        double mcdx = mctrk.Start().T() * 1.e-3 * detProp.DriftVelocity(); // cm
        double plen = length(clockData, detProp, mctrk, mcdx, mcstart, mcend, mcstartmom, mcendmom);
        mf::LogVerbatim("TrackAna")
          << evt.run() << "." << evt.event() << " NoMat MCTkID " << std::setw(6) << mctrk.TrackID()
          << " Origin " << mctrk.Origin() << " PDG" << std::setw(5) << mctrk.PdgCode() << " KE"
          << std::setw(4) << (int)KE << KEUnits << " Length " << std::fixed << std::setprecision(1)
          << plen << " cm";
        if (fPrintLevel > 1) {
          // this won't work for DUNE
          for (unsigned short ipl = 0; ipl < geom->Nplanes(); ++ipl) {
            geo::PlaneID const& planeID{0, 0, ipl};
            auto const sWire = geom->NearestWireID(mcstart, planeID).Wire;
            auto const sTick = detProp.ConvertXToTicks(mcstart.X(), planeID);
            auto const eWire = geom->NearestWireID(mcend, planeID).Wire;
            auto const eTick = detProp.ConvertXToTicks(mcend.X(), planeID);
            mf::LogVerbatim("TrackAna") << "   Wire:Tick in Pln " << ipl << " W:T " << sWire << ":"
                                        << sTick << " - " << eWire << ":" << eTick;
          } // ipl
        }   // fPrintLevel > 1
      }     // imc
    }       // fPrintLevel > 0
  }

  void TrackAna::anaStitch(const art::Event& evt,
                           detinfo::DetectorClocksData const& clockData,
                           detinfo::DetectorPropertiesData const& detProp)
  {
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    art::ServiceHandle<geo::Geometry const> geom;

    std::map<int, std::map<int, art::PtrVector<recob::Hit>>> hitmap; // trkID, otrk, hitvec
    std::map<int, int> KEmap; // length traveled in det [cm]?, trkID want to sort by KE
    bool mc = !evt.isRealData();
    art::Handle<std::vector<recob::Track>> trackh;
    art::Handle<std::vector<recob::SpacePoint>> sppth;
    art::Handle<std::vector<art::PtrVector<recob::Track>>> trackvh;

    evt.getByLabel(fTrackModuleLabel, trackh);
    evt.getByLabel(fSpacepointModuleLabel, sppth);
    evt.getByLabel(fStitchModuleLabel, trackvh);
    int ntv(trackvh->size());

    std::vector<art::PtrVector<recob::Track>>::const_iterator cti = trackvh->begin();

    if (trackh.isValid()) {
      art::FindManyP<recob::SpacePoint> fswhole(trackh, evt, fTrkSpptAssocModuleLabel);
      int nsppts_assnwhole = fswhole.size();
      std::cout << "TrackAna: Number of clumps of Spacepoints from Assn for all Tracks: "
                << nsppts_assnwhole << std::endl;
    }

    if (fRecoHistMap.count(0) == 0) {
      fRecoHistMap.emplace(0, RecoHists{"Reco"});
      std::cout << "\n"
                << "\t\t  TrkAna: Fresh fRecoHistMap[0] ******* \n"
                << std::endl;
    }
    const RecoHists& rhistsStitched = fRecoHistMap.at(0);

    std::vector<std::vector<unsigned int>> NtrkIdsAll;
    std::vector<double> ntvsorted;
    hitmap.clear();
    KEmap.clear();

    // Look at the components of the stitched tracks. Grab their sppts/hits from Assns.
    for (int o = 0; o < ntv; ++o) // o for outer
    {
      const art::PtrVector<recob::Track> pvtrack(*(cti++));
      int ntrack = pvtrack.size();
      std::vector<std::vector<unsigned int>> NtrkId_Hit; // hit IDs in inner tracks
      std::vector<unsigned int> vecMode;
      art::FindManyP<recob::SpacePoint> fs(pvtrack, evt, fTrkSpptAssocModuleLabel);

      for (int i = 0; i < ntrack; ++i) {

        // From gdb> ptype fs, the vector of Ptr<SpacePoint>s it appears is
        // grabbed after fs.at(0)
        bool assns(true);
        try {
          // Got Spacepoints from this Track; now get Hits from those
          // Spacepoints.

          int nsppts_assn = fs.at(i).size();

          const auto& sppt = fs.at(i); //.at(0);
          // since we're in a try and worried about failure, we won't pull the following
          // FindManyP out of the loop.
          art::FindManyP<recob::Hit> fh(sppt, evt, fHitSpptAssocModuleLabel);

          // Importantly, loop on all sppts, though they don't all contribute to the track.
          // As opposed to looping on the trajectory pts, which is a lower number.
          // Also, important, in job in whch this runs I set TrackKal3DSPS parameter MaxPass=1,
          // cuz I don't want merely the sparse set of sppts as follows from the uncontained
          // MS-measurement in 2nd pass.
          std::vector<unsigned int> vecNtrkIds;
          for (int is = 0; is < nsppts_assn; ++is) {
            int nhits = fh.at(is).size(); // should be 2 or 3: number of planes.
            for (int ih = 0; ih < nhits; ++ih) {
              const auto& hit = fh.at(is).at(ih); // Our vector is after the .at(is) this time.
              if (hit->SignalType() != geo::kCollection) continue;
              rhistsStitched.fHHitChg->Fill(hit->Integral());
              rhistsStitched.fHHitWidth->Fill(2. * hit->RMS());
              if (mc) {
                std::vector<sim::TrackIDE> tids = bt_serv->HitToTrackIDEs(clockData, hit);
                // more here.
                // Loop over track ids.
                bool justOne(true); // Only take first trk that contributed to this hit
                for (std::vector<sim::TrackIDE>::const_iterator itid = tids.begin();
                     itid != tids.end();
                     ++itid) {
                  int trackID = std::abs(itid->trackID);
                  hitmap[trackID][o].push_back(hit);

                  if (justOne) {
                    vecNtrkIds.push_back(trackID);
                    justOne = false;
                  }
                  // Add hit to PtrVector corresponding to this track id.
                  rhistsStitched.fHHitTrkId->Fill(trackID);
                  const simb::MCParticle* part = pi_serv->TrackIdToParticle_P(trackID);
                  if (!part) break;

                  rhistsStitched.fHHitPdg->Fill(part->PdgCode());
                  // This really needs to be indexed as KE deposited in volTPC, not just KE. EC, 24-July-2014.

                  geo::Point_t mcstart;
                  geo::Point_t mcend;
                  TVector3 mcstartmom;
                  TVector3 mcendmom;
                  double mctime = part->T();                              // nsec
                  double mcdx = mctime * 1.e-3 * detProp.DriftVelocity(); // cm

                  double plen =
                    length(clockData, detProp, *part, mcdx, mcstart, mcend, mcstartmom, mcendmom);

                  KEmap[(int)(1e6 * plen)] = trackID; // multiple assignment but always the same, so
                                                      // fine.
                }

              } // mc
            }   //  hits

          } //    spacepoints

          if (mc) {
            NtrkId_Hit.push_back(vecNtrkIds);
            // Find the trkID mode for this i^th track
            unsigned int ii(1);
            int max(-12), n(1), ind(0);
            std::sort(vecNtrkIds.begin(), vecNtrkIds.end());
            std::vector<unsigned int> strkIds(vecNtrkIds);
            while (ii < vecNtrkIds.size()) {
              if (strkIds.at(ii) != strkIds.at(ii - 1)) { n = 1; }
              else {
                n++;
              }
              if (n > max) {
                max = n;
                ind = ii;
              }
              ii++;
            }
            unsigned int mode(sim::NoParticleId);
            if (strkIds.begin() != strkIds.end()) mode = strkIds.at(ind);
            vecMode.push_back(mode);

            if (strkIds.size() != 0)
              rhistsStitched.fModeFrac->Fill((double)max / (double)strkIds.size());
            else
              rhistsStitched.fModeFrac->Fill(-1.0);
          } // mc

        } // end try
        catch (cet::exception& x) {
          assns = false;
        }
        if (!assns) throw cet::exception("TrackAna") << "Bad Associations. \n";

      } // i

      if (mc) {
        // one vector per o trk, for all modes of stitched i trks
        NtrkIdsAll.push_back(vecMode);

        std::unique(NtrkIdsAll.back().begin(), NtrkIdsAll.back().end());
        double sum(0.0);
        for (auto const val : NtrkIdsAll.back()) {
          sum += hitmap[val][o].size();
        }
        ntvsorted.push_back(sum);
      }

      //
    } // o

    int vtmp(0);
    // get KEmap indices by most energetic first, least last.
    for (auto it = KEmap.rbegin(); it != KEmap.rend(); ++it) {
      vtmp++;
    }

    // get o trk indices by their hits. Most populous track first, least last.
    for (auto const o : fsort_indexes(ntvsorted)) {
      int v(0);
      // get KEmap indices by longest trajectory first, least last.
      for (auto it = KEmap.rbegin(); it != KEmap.rend(); ++it) {
        int val = it->second; // grab trkIDs in order, since they're sorted by KE
        rhistsStitched.fNTrkIdTrks3->Fill(o, v, hitmap[val][o].size());
        v++;
      }
    }

    // In how many o tracks did each trkId appear? Histo it. Would like it to be
    // precisely 1. Histo it vs. particle KE.
    flattener flat(NtrkIdsAll);
    std::vector<unsigned int>& v = flat;
    //  auto const it ( std::unique(v.begin(),v.end()) ); // never use this it,
    //  perhaps.
    for (auto const val : v) {
      if (val != (unsigned int)sim::NoParticleId) {
        const simb::MCParticle* part = pi_serv->TrackIdToParticle_P(val);
        double T(part->E() - 0.001 * part->Mass());
        rhistsStitched.fNTrkIdTrks->Fill(std::count(v.begin(), v.end(), val));
        rhistsStitched.fNTrkIdTrks2->Fill(std::count(v.begin(), v.end(), val), T);
      }
      else {
        rhistsStitched.fNTrkIdTrks2->Fill(-1.0, 0.0);
      }
    }
  }

  void TrackAna::endJob()
  {
    // Print summary.

    mf::LogInfo("TrackAna") << "TrackAna statistics:\n"
                            << "  Number of events = " << fNumEvent;

    // Fill efficiency histograms.

    for (std::map<int, MCHists>::const_iterator i = fMCHistMap.begin(); i != fMCHistMap.end();
         ++i) {
      const MCHists& mchists = i->second;
      effcalc(mchists.fHgstartx, mchists.fHmcstartx, mchists.fHestartx);
      effcalc(mchists.fHgstarty, mchists.fHmcstarty, mchists.fHestarty);
      effcalc(mchists.fHgstartz, mchists.fHmcstartz, mchists.fHestartz);
      effcalc(mchists.fHgendx, mchists.fHmcendx, mchists.fHeendx);
      effcalc(mchists.fHgendy, mchists.fHmcendy, mchists.fHeendy);
      effcalc(mchists.fHgendz, mchists.fHmcendz, mchists.fHeendz);
      effcalc(mchists.fHgtheta, mchists.fHmctheta, mchists.fHetheta);
      effcalc(mchists.fHgphi, mchists.fHmcphi, mchists.fHephi);
      effcalc(mchists.fHgtheta_xz, mchists.fHmctheta_xz, mchists.fHetheta_xz);
      effcalc(mchists.fHgtheta_yz, mchists.fHmctheta_yz, mchists.fHetheta_yz);
      effcalc(mchists.fHgmom, mchists.fHmcmom, mchists.fHemom);
      effcalc(mchists.fHgmoml, mchists.fHmcmoml, mchists.fHemoml);
      effcalc(mchists.fHgke, mchists.fHmcke, mchists.fHeke);
      effcalc(mchists.fHgkel, mchists.fHmckel, mchists.fHekel);
      effcalc(mchists.fHglen, mchists.fHmclen, mchists.fHelen);
      effcalc(mchists.fHglens, mchists.fHmclens, mchists.fHelens);
    }
  }

  // Stole this from online. Returns indices sorted by corresponding vector values.
  std::vector<size_t> TrackAna::fsort_indexes(const std::vector<double>& v)
  {
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i)
      idx[i] = i;
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
      return v[i1] > v[i2];
    }); // Want most occupied trks first. was <, EC, 21-July.
    return idx;
  }

}

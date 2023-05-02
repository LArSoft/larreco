//
// Name: SeedAna_module.cc
//
// Purpose: Module SeedAna.
//
// Configuration parameters.
//
//  SeedModuleLabel:    Seed module label.
//  MCTrackModuleLabel: MCTrack module label.
//  MinMCKE:            Minimum MC particle kinetic energy.
//  MatchColinearity:   Minimum colinearity for mc-seed matching.
//  MatchDisp:          Maximum uv displacement for mc-seed matching.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Seed.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixD.h"

namespace {

  // Calculate distance to boundary.
  //----------------------------------------------------------------------------
  double bdist(const TVector3& pos)
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

  // Find the closest matching mc trajectory point (sim::MCStep) for a given seed.
  // Returned value is index of the trajectory point.
  // Return -1 in case of no match.
  int mcmatch(detinfo::DetectorPropertiesData const& detProp,
              const sim::MCTrack& mctrk,
              const recob::Seed& seed)
  {
    // Get seed point.

    double pos[3];
    double err[3];
    seed.GetPoint(pos, err);

    // Calculate the x offset due to nonzero mc particle time.
    double mctime = mctrk.Start().T();                      // nsec
    double mcdx = mctime * 1.e-3 * detProp.DriftVelocity(); // cm

    // Loop over trajectory points.

    int best_traj = -1;
    double max_dist = 0.;
    int ntraj = mctrk.size();
    for (int itraj = 0; itraj < ntraj; ++itraj) {
      const TLorentzVector& vec = mctrk[itraj].Position();
      double dx = pos[0] - vec.X() - mcdx;
      double dy = pos[1] - vec.Y();
      double dz = pos[2] - vec.Z();
      double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (best_traj < 0 || dist < max_dist) {
        best_traj = itraj;
        max_dist = dist;
      }
    }
    return best_traj;
  }

  // Length of MCTrack.
  // In this function, the extracted start and end momenta are converted to GeV
  // (MCTrack stores momenta in Mev).
  //----------------------------------------------------------------------------
  double length(detinfo::DetectorPropertiesData const& detProp,
                const sim::MCTrack& mctrk,
                double dx,
                TVector3& start,
                TVector3& end,
                TVector3& startmom,
                TVector3& endmom,
                unsigned int /*tpc*/ = 0,
                unsigned int /*cstat*/ = 0)
  {
    // Get services.

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
      throw cet::exception("SeedAna") << "effcalc[" __FILE__ "]: incompatible histograms (I)\n";
    if (nbins != heff->GetNbinsX())
      throw cet::exception("SeedAna") << "effcalc[" __FILE__ "]: incompatible histograms (II)\n";

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

  // Fill multiplicity histogram.

  void mulcalc(const TH1* hnum, const TH1* hden, TH1* hmul)
  {
    int nbins = hnum->GetNbinsX();
    if (nbins != hden->GetNbinsX())
      throw cet::exception("SeedAna") << "mulcalc[" __FILE__ "]: incompatible histograms (I)\n";
    if (nbins != hmul->GetNbinsX())
      throw cet::exception("SeedAna") << "mulcalc[" __FILE__ "]: incompatible histograms (II)\n";

    // Loop over bins, including underflow and overflow.

    for (int ibin = 0; ibin <= nbins + 1; ++ibin) {
      double num = hnum->GetBinContent(ibin);
      double den = hden->GetBinContent(ibin);
      if (den == 0.) {
        hmul->SetBinContent(ibin, 0.);
        hmul->SetBinError(ibin, 0.);
      }
      else {
        double mul = num / den;
        if (mul < 0.) mul = 0.;
        double err = std::sqrt((1. + mul) * mul / den);
        hmul->SetBinContent(ibin, mul);
        hmul->SetBinError(ibin, err);
      }
    }

    hmul->SetMinimum(0.);
    hmul->SetMarkerStyle(20);
  }
}

namespace trkf {

  class SeedAna : public art::EDAnalyzer {
  public:
    // Embedded structs.

    // Struct for histograms that depend on seeds only.

    struct RecoHists {
      RecoHists(const std::string& subdir);

      // Pure reco seed histograms.

      TH1F* fHx{nullptr};        // Seed x position.
      TH1F* fHy{nullptr};        // Seed y position.
      TH1F* fHz{nullptr};        // Seed z position.
      TH1F* fHdist{nullptr};     // Seed distance to boundary.
      TH1F* fHtheta{nullptr};    // Theta.
      TH1F* fHphi{nullptr};      // Phi.
      TH1F* fHtheta_xz{nullptr}; // Theta_xz.
      TH1F* fHtheta_yz{nullptr}; // Theta_yz.
    };

    // Struct for mc particles and mc-matched tracks.

    struct MCHists {
      MCHists(const std::string& subdir);

      // Reco-MC matching.

      TH2F* fHduvcosth{nullptr}; // 2D mc vs. data matching, duv vs. cos(theta).
      TH1F* fHcosth{nullptr};    // 1D direction matching, cos(theta).
      TH1F* fHmcu{nullptr};      // 1D endpoint truth u.
      TH1F* fHmcv{nullptr};      // 1D endpoint truth v.
      TH1F* fHmcw{nullptr};      // 1D endpoint truth w.
      TH1F* fHmcdudw{nullptr};   // Truth du/dw.
      TH1F* fHmcdvdw{nullptr};   // Truth dv/dw.

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
      TH1F* fHmclen{nullptr};      // Length.

      // Matched seed histograms (multiplicity numerator).

      TH1F* fHmstartx{nullptr};   // Starting x position.
      TH1F* fHmstarty{nullptr};   // Starting y position.
      TH1F* fHmstartz{nullptr};   // Starting z position.
      TH1F* fHmendx{nullptr};     // Ending x position.
      TH1F* fHmendy{nullptr};     // Ending y position.
      TH1F* fHmendz{nullptr};     // Ending z position.
      TH1F* fHmtheta{nullptr};    // Theta.
      TH1F* fHmphi{nullptr};      // Phi.
      TH1F* fHmtheta_xz{nullptr}; // Theta_xz.
      TH1F* fHmtheta_yz{nullptr}; // Theta_yz.
      TH1F* fHmmom{nullptr};      // Momentum.
      TH1F* fHmlen{nullptr};      // Length.

      // Matched seed histograms (efficiency numerator).

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
      TH1F* fHglen{nullptr};      // Length.

      // Multiplicity histograms.

      TH1F* fHmulstartx{nullptr};   // Starting x position.
      TH1F* fHmulstarty{nullptr};   // Starting y position.
      TH1F* fHmulstartz{nullptr};   // Starting z position.
      TH1F* fHmulendx{nullptr};     // Ending x position.
      TH1F* fHmulendy{nullptr};     // Ending y position.
      TH1F* fHmulendz{nullptr};     // Ending z position.
      TH1F* fHmultheta{nullptr};    // Theta.
      TH1F* fHmulphi{nullptr};      // Phi.
      TH1F* fHmultheta_xz{nullptr}; // Theta_xz.
      TH1F* fHmultheta_yz{nullptr}; // Theta_yz.
      TH1F* fHmulmom{nullptr};      // Momentum.
      TH1F* fHmullen{nullptr};      // Length.

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
      TH1F* fHelen{nullptr};      // Length.
    };

    explicit SeedAna(fhicl::ParameterSet const& pset);

  private:
    void analyze(const art::Event& evt) override;
    void endJob() override;

    // Fcl Attributes.

    std::string fSeedModuleLabel;
    std::string fMCTrackModuleLabel;
    int fDump;                // Number of events to dump to debug message facility.
    double fMinMCKE;          // Minimum MC particle kinetic energy (GeV).
    double fMinMCLen;         // Minimum MC particle length in tpc (cm).
    double fMatchColinearity; // Minimum matching colinearity.
    double fMatchDisp;        // Maximum matching displacement.
    bool fIgnoreSign;         // Ignore sign of mc particle if true.

    // Histograms.

    std::map<int, MCHists> fMCHistMap;     // Indexed by pdg id.
    std::map<int, RecoHists> fRecoHistMap; // Indexed by pdg id.

    // Statistics.

    int fNumEvent;
  };

  DEFINE_ART_MODULE(SeedAna)

  // RecoHists methods.

  SeedAna::RecoHists::RecoHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("seedana", "SeedAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHx =
      dir.make<TH1F>("x", "X Position", 100, -2. * geom->DetHalfWidth(), 4. * geom->DetHalfWidth());
    fHy = dir.make<TH1F>("y", "Y Position", 100, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHz = dir.make<TH1F>("z", "Z Position", 100, 0., geom->DetLength());
    fHdist =
      dir.make<TH1F>("dist", "Position Distance to Boundary", 100, -10., geom->DetHalfWidth());
    fHtheta = dir.make<TH1F>("theta", "Theta", 100, 0., 3.142);
    fHphi = dir.make<TH1F>("phi", "Phi", 100, -3.142, 3.142);
    fHtheta_xz = dir.make<TH1F>("theta_xz", "Theta_xz", 100, -3.142, 3.142);
    fHtheta_yz = dir.make<TH1F>("theta_yz", "Theta_yz", 100, -3.142, 3.142);
  }

  // MCHists methods.

  SeedAna::MCHists::MCHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("seedana", "SeedAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHduvcosth =
      dir.make<TH2F>("duvcosth", "Delta(uv) vs. Colinearity", 100, 0.95, 1., 100, 0., 1.);
    fHcosth = dir.make<TH1F>("colin", "Colinearity", 100, 0.95, 1.);
    fHmcu = dir.make<TH1F>("mcu", "MC Truth U", 100, -5., 5.);
    fHmcv = dir.make<TH1F>("mcv", "MC Truth V", 100, -5., 5.);
    fHmcw = dir.make<TH1F>("mcw", "MC Truth W", 100, -20., 20.);
    fHmcdudw = dir.make<TH1F>("mcdudw", "MC Truth U Slope", 100, -0.2, 0.2);
    fHmcdvdw = dir.make<TH1F>("mcdvdw", "MV Truth V Slope", 100, -0.2, 0.2);

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
    fHmclen = dir.make<TH1F>("mclen", "MC Particle Length", 10, 0., 1.1 * geom->DetLength());

    fHmstartx = dir.make<TH1F>("mxstart",
                               "Matched X Start Position",
                               10,
                               -2. * geom->DetHalfWidth(),
                               4. * geom->DetHalfWidth());
    fHmstarty = dir.make<TH1F>(
      "mystart", "Matched Y Start Position", 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmstartz = dir.make<TH1F>("mzstart", "Matched Z Start Position", 10, 0., geom->DetLength());
    fHmendx = dir.make<TH1F>(
      "mxend", "Matched X End Position", 10, -2. * geom->DetHalfWidth(), 4. * geom->DetHalfWidth());
    fHmendy = dir.make<TH1F>(
      "myend", "Matched Y End Position", 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmendz = dir.make<TH1F>("mzend", "Matched Z End Position", 10, 0., geom->DetLength());
    fHmtheta = dir.make<TH1F>("mtheta", "Matched Theta", 20, 0., 3.142);
    fHmphi = dir.make<TH1F>("mphi", "Matched Phi", 10, -3.142, 3.142);
    fHmtheta_xz = dir.make<TH1F>("mtheta_xz", "Matched Theta_xz", 40, -3.142, 3.142);
    fHmtheta_yz = dir.make<TH1F>("mtheta_yz", "Matched Theta_yz", 40, -3.142, 3.142);
    fHmmom = dir.make<TH1F>("mmom", "Matched Momentum", 10, 0., 10.);
    fHmlen = dir.make<TH1F>("mlen", "Matched Particle Length", 10, 0., 1.1 * geom->DetLength());

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
    fHglen = dir.make<TH1F>("glen", "Good Particle Length", 10, 0., 1.1 * geom->DetLength());

    fHmulstartx = dir.make<TH1F>("mulxstart",
                                 "Multiplicity vs. X Start Position",
                                 10,
                                 -2. * geom->DetHalfWidth(),
                                 4. * geom->DetHalfWidth());
    fHmulstarty = dir.make<TH1F>("mulystart",
                                 "Multiplicity vs. Y Start Position",
                                 10,
                                 -geom->DetHalfHeight(),
                                 geom->DetHalfHeight());
    fHmulstartz =
      dir.make<TH1F>("mulzstart", "Multiplicity vs. Z Start Position", 10, 0., geom->DetLength());
    fHmulendx = dir.make<TH1F>("mulxend",
                               "Multiplicity vs. X End Position",
                               10,
                               -2. * geom->DetHalfWidth(),
                               4. * geom->DetHalfWidth());
    fHmulendy = dir.make<TH1F>("mulyend",
                               "Multiplicity vs. Y End Position",
                               10,
                               -geom->DetHalfHeight(),
                               geom->DetHalfHeight());
    fHmulendz =
      dir.make<TH1F>("mulzend", "Multiplicity vs. Z End Position", 10, 0., geom->DetLength());
    fHmultheta = dir.make<TH1F>("multheta", "Multiplicity vs. Theta", 20, 0., 3.142);
    fHmulphi = dir.make<TH1F>("mulphi", "Multiplicity vs. Phi", 10, -3.142, 3.142);
    fHmultheta_xz = dir.make<TH1F>("multheta_xz", "Multiplicity vs. Theta_xz", 40, -3.142, 3.142);
    fHmultheta_yz = dir.make<TH1F>("multheta_yz", "Multiplicity vs. Theta_yz", 40, -3.142, 3.142);
    fHmulmom = dir.make<TH1F>("mulmom", "Multiplicity vs. Momentum", 10, 0., 10.);
    fHmullen =
      dir.make<TH1F>("mullen", "Multiplicity vs. Particle Length", 10, 0., 1.1 * geom->DetLength());

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
    fHelen =
      dir.make<TH1F>("elen", "Efficiency vs. Particle Length", 10, 0., 1.1 * geom->DetLength());
  }

  SeedAna::SeedAna(const fhicl::ParameterSet& pset)
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    : EDAnalyzer(pset)
    , fSeedModuleLabel(pset.get<std::string>("SeedModuleLabel"))
    , fMCTrackModuleLabel(pset.get<std::string>("MCTrackModuleLabel"))
    , fDump(pset.get<int>("Dump"))
    , fMinMCKE(pset.get<double>("MinMCKE"))
    , fMinMCLen(pset.get<double>("MinMCLen"))
    , fMatchColinearity(pset.get<double>("MatchColinearity"))
    , fMatchDisp(pset.get<double>("MatchDisp"))
    , fIgnoreSign(pset.get<bool>("IgnoreSign"))
    , fNumEvent(0)
  {
    mf::LogInfo("SeedAna") << "SeedAna configured with the following parameters:\n"
                           << "  SeedModuleLabel = " << fSeedModuleLabel << "\n"
                           << "  MCTrackModuleLabel = " << fMCTrackModuleLabel << "\n"
                           << "  Dump = " << fDump << "\n"
                           << "  MinMCKE = " << fMinMCKE << "\n"
                           << "  MinMCLen = " << fMinMCLen;
  }

  void SeedAna::analyze(const art::Event& evt)
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Optional dump stream.

    std::unique_ptr<mf::LogInfo> pdump;
    if (fDump > 0) {
      --fDump;
      pdump = std::unique_ptr<mf::LogInfo>(new mf::LogInfo("TrackAna"));
    }

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();

    // Get seed handle.

    art::Handle<std::vector<recob::Seed>> seedh;
    evt.getByLabel(fSeedModuleLabel, seedh);

    // Seed->mc track matching map.

    std::map<const recob::Seed*, int> seedmap;

    if (mc) {

      // Get MCTracks.

      art::Handle<std::vector<sim::MCTrack>> mctrackh;
      evt.getByLabel(fMCTrackModuleLabel, mctrackh);

      // Dump MCTracks.

      if (pdump) {
        *pdump << "MC Tracks\n"
               << "       Id   pdg           x         y         z          dx        dy        dz "
                  "          p\n"
               << "--------------------------------------------------------------------------------"
                  "-----------\n";
      }

      // Loop over mc tracks, and fill histograms that depend only
      // on mc particles.
      auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

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

          // Apply minimum energy cut.

          if (mctrk.Start().E() >= mctrk.Start().Momentum().Mag() + 1000. * fMinMCKE) {

            // Calculate the x offset due to nonzero mc particle time.

            double mctime = mctrk.Start().T();                      // nsec
            double mcdx = mctime * 1.e-3 * detProp.DriftVelocity(); // cm

            // Calculate the length of this mc particle inside the fiducial volume.

            TVector3 mcstart;
            TVector3 mcend;
            TVector3 mcstartmom;
            TVector3 mcendmom;
            double plen = length(detProp, mctrk, mcdx, mcstart, mcend, mcstartmom, mcendmom);

            // Apply minimum fiducial length cut.  Always reject particles that have
            // zero length in the tpc regardless of the configured cut.

            if (plen > 0. && plen > fMinMCLen) {

              // Dump MC particle information here.

              if (pdump) {
                double pstart = mcstartmom.Mag();
                double pend = mcendmom.Mag();
                *pdump << "\nOffset" << std::setw(3) << mctrk.TrackID() << std::setw(6)
                       << mctrk.PdgCode() << "  " << std::fixed << std::setprecision(2)
                       << std::setw(10) << mcdx << "\nStart " << std::setw(3) << mctrk.TrackID()
                       << std::setw(6) << mctrk.PdgCode() << "  " << std::fixed
                       << std::setprecision(2) << std::setw(10) << mcstart[0] << std::setw(10)
                       << mcstart[1] << std::setw(10) << mcstart[2];
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
                       << std::setw(10) << mcend[0] << std::setw(10) << mcend[1] << std::setw(10)
                       << mcend[2];
                if (pend > 0.01) {
                  *pdump << "  " << std::fixed << std::setprecision(3) << std::setw(10)
                         << mcendmom[0] / pend << std::setw(10) << mcendmom[1] / pend
                         << std::setw(10) << mcendmom[2] / pend;
                }
                else
                  *pdump << std::setw(32) << " ";
                *pdump << std::setw(12) << std::fixed << std::setprecision(3) << pend << "\n";
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
              mchists.fHmcmom->Fill(mcstartmom.Mag());
              mchists.fHmclen->Fill(plen);

              // Loop over seeds and do matching.

              int nmatch = 0;
              if (seedh.isValid()) {

                // Loop over seeds.

                int nseed = seedh->size();
                for (int i = 0; i < nseed; ++i) {
                  art::Ptr<recob::Seed> pseed(seedh, i);
                  const recob::Seed& seed = *pseed;
                  if (seedmap.count(&seed) == 0) seedmap[&seed] = -1;

                  // Get parameters of this seed.

                  TVector3 pos;
                  TVector3 dir;
                  double err[3];
                  seed.GetPoint(&pos[0], err);
                  seed.GetDirection(&dir[0], err);

                  // Calculate the global-to-local rotation matrix.
                  // Copied from Track.cxx.

                  TMatrixD rot(3, 3);
                  double dirmag = dir.Mag();
                  double diryz = std::sqrt(dir.Y() * dir.Y() + dir.Z() * dir.Z());

                  double sinth = dir.X() / dirmag;
                  double costh = diryz / dirmag;
                  double sinphi = 0.;
                  double cosphi = 1.;
                  if (diryz != 0) {
                    sinphi = -dir.Y() / diryz;
                    cosphi = dir.Z() / diryz;
                  }
                  rot(0, 0) = costh;
                  rot(1, 0) = 0.;
                  rot(2, 0) = sinth;
                  rot(0, 1) = sinth * sinphi;
                  rot(1, 1) = cosphi;
                  rot(2, 1) = -costh * sinphi;
                  rot(0, 2) = -sinth * cosphi;
                  rot(1, 2) = sinphi;
                  rot(2, 2) = costh * cosphi;

                  // Get best matching mc trajectory point.

                  int itraj = mcmatch(detProp, mctrk, seed);
                  if (itraj >= 0) {

                    // Get mc relative position and direction at matching trajectory point.

                    TVector3 mcpos = mctrk[itraj].Position().Vect() - pos;
                    TVector3 mcmom = mctrk[itraj].Momentum().Vect();
                    mcpos[0] += mcdx;

                    // Rotate the momentum and position to the
                    // seed-local coordinate system.

                    TVector3 mcmoml = rot * mcmom;
                    TVector3 mcposl = rot * mcpos;

                    if (mcmoml.Z() < 0.) mcmoml = -mcmoml;
                    double costh = mcmoml.Z() / mcmoml.Mag();

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
                    double uv0 = std::sqrt(u0 * u0 + v0 * v0);

                    // Fill matching histograms.

                    mchists.fHduvcosth->Fill(costh, uv0);
                    if (std::abs(uv0) < fMatchDisp) {

                      // Fill slope matching histograms.

                      mchists.fHmcdudw->Fill(dudw);
                      mchists.fHmcdvdw->Fill(dvdw);
                    }
                    mchists.fHcosth->Fill(costh);
                    if (costh > fMatchColinearity) {

                      // Fill displacement matching histograms.

                      mchists.fHmcu->Fill(u0);
                      mchists.fHmcv->Fill(v0);
                      mchists.fHmcw->Fill(w);

                      if (std::abs(uv0) < fMatchDisp) {

                        // Now we have passed all matching cuts and we have a matching
                        // mc particle + seed pair.

                        ++nmatch;
                        seedmap[&seed] = mctrk.TrackID();

                        // Fill matched seed histograms (seed multiplicity).

                        mchists.fHmstartx->Fill(mcstart.X());
                        mchists.fHmstarty->Fill(mcstart.Y());
                        mchists.fHmstartz->Fill(mcstart.Z());
                        mchists.fHmendx->Fill(mcend.X());
                        mchists.fHmendy->Fill(mcend.Y());
                        mchists.fHmendz->Fill(mcend.Z());
                        mchists.fHmtheta->Fill(mcstartmom.Theta());
                        mchists.fHmphi->Fill(mcstartmom.Phi());
                        mchists.fHmtheta_xz->Fill(mctheta_xz);
                        mchists.fHmtheta_yz->Fill(mctheta_yz);
                        mchists.fHmmom->Fill(mcstartmom.Mag());
                        mchists.fHmlen->Fill(plen);
                      }
                    }
                  }
                }

                // If we found at least one matched seed, fill good
                // particle histograms.

                if (nmatch > 0) {
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
                  mchists.fHgmom->Fill(mcstartmom.Mag());
                  mchists.fHglen->Fill(plen);
                }
              }
            }
          }
        }
      }
    }

    // Loop over seeds and fill reco-only seed histograms.

    if (seedh.isValid()) {

      // Loop over seeds.

      int nseed = seedh->size();

      if (nseed > 0 && pdump != 0) {
        *pdump << "\nReconstructed Seeds\n"
               << "           MCid           x         y         z          dx        dy        dz "
                  "          p\n"
               << "--------------------------------------------------------------------------------"
                  "-----------\n";
      }

      for (int i = 0; i < nseed; ++i) {
        art::Ptr<recob::Seed> pseed(seedh, i);
        const recob::Seed& seed = *pseed;

        // Fill histograms involving reco seeds only.

        TVector3 pos;
        TVector3 dir;
        double err[3];
        seed.GetPoint(&pos[0], err);
        seed.GetDirection(&dir[0], err);
        double mdir = dir.Mag();
        if (mdir != 0.) { dir *= (1. / mdir); }

        double dpos = bdist(pos);
        double theta_xz = std::atan2(dir.X(), dir.Z());
        double theta_yz = std::atan2(dir.Y(), dir.Z());

        // Dump seed information here.

        if (pdump) {
          int mcid = seedmap[&seed];
          *pdump << std::setw(15) << mcid << "  " << std::fixed << std::setprecision(2)
                 << std::setw(10) << pos[0] << std::setw(10) << pos[1] << std::setw(10) << pos[2]
                 << "  " << std::fixed << std::setprecision(3) << std::setw(10) << dir[0]
                 << std::setw(10) << dir[1] << std::setw(10) << dir[2] << "\n";
        }

        // Fill histograms.

        if (fRecoHistMap.count(0) == 0) fRecoHistMap.emplace(0, RecoHists{"Reco"});
        const RecoHists& rhists = fRecoHistMap.at(0);

        rhists.fHx->Fill(pos.X());
        rhists.fHy->Fill(pos.Y());
        rhists.fHz->Fill(pos.Z());
        rhists.fHdist->Fill(dpos);
        rhists.fHtheta->Fill(dir.Theta());
        rhists.fHphi->Fill(dir.Phi());
        rhists.fHtheta_xz->Fill(theta_xz);
        rhists.fHtheta_yz->Fill(theta_yz);
      }
    }
  }

  void SeedAna::endJob()
  //
  // Purpose: End of job.
  //
  {
    // Print summary.

    mf::LogInfo("SeedAna") << "SeedAna statistics:\n"
                           << "  Number of events = " << fNumEvent;

    // Fill multiplicity histograms.

    for (std::map<int, MCHists>::const_iterator i = fMCHistMap.begin(); i != fMCHistMap.end();
         ++i) {
      const MCHists& mchists = i->second;
      mulcalc(mchists.fHmstartx, mchists.fHmcstartx, mchists.fHmulstartx);
      mulcalc(mchists.fHmstarty, mchists.fHmcstarty, mchists.fHmulstarty);
      mulcalc(mchists.fHmstartz, mchists.fHmcstartz, mchists.fHmulstartz);
      mulcalc(mchists.fHmendx, mchists.fHmcendx, mchists.fHmulendx);
      mulcalc(mchists.fHmendy, mchists.fHmcendy, mchists.fHmulendy);
      mulcalc(mchists.fHmendz, mchists.fHmcendz, mchists.fHmulendz);
      mulcalc(mchists.fHmtheta, mchists.fHmctheta, mchists.fHmultheta);
      mulcalc(mchists.fHmphi, mchists.fHmcphi, mchists.fHmulphi);
      mulcalc(mchists.fHmtheta_xz, mchists.fHmctheta_xz, mchists.fHmultheta_xz);
      mulcalc(mchists.fHmtheta_yz, mchists.fHmctheta_yz, mchists.fHmultheta_yz);
      mulcalc(mchists.fHmmom, mchists.fHmcmom, mchists.fHmulmom);
      mulcalc(mchists.fHmlen, mchists.fHmclen, mchists.fHmullen);
    }
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
      effcalc(mchists.fHglen, mchists.fHmclen, mchists.fHelen);
    }
  }
}

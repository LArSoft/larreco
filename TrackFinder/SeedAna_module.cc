//
// Name: SeedAna_module.cc
//
// Purpose: Module SeedAna.
//
// Configuration parameters.
//
//  SeedModuleLabel:   Seed module label.
//  MinMCKE:           Minimum MC particle kinetic energy.
//  MatchColinearity:  Minimum colinearity for mc-seed matching.
//  MatchDisp:         Maximum uv displacement for mc-seed matching.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include <map>
#include <iostream>
#include <sstream>
#include <cmath>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "RecoBase/Seed.h"
#include "MCCheater/BackTracker.h"
#include "SimulationBase/MCParticle.h"

#include "TH2F.h"
#include "TFile.h"
#include "TMatrixD.h"

namespace {

  // Local functions.

  // Calculate distance to boundary.
  //----------------------------------------------------------------------------
  double bdist(const TVector3& pos, unsigned int tpc = 0, unsigned int cstat = 0)
  {
    // Get geometry.

    art::ServiceHandle<geo::Geometry> geom;

    double d1 = pos.X();                             // Distance to right side (wires).
    double d2 = 2.*geom->DetHalfWidth() - pos.X();   // Distance to left side (cathode).
    double d3 = pos.Y() + geom->DetHalfHeight();     // Distance to bottom.
    double d4 = geom->DetHalfHeight() - pos.Y();     // Distance to top.
    double d5 = pos.Z();                             // Distance to front.
    double d6 = geom->DetLength() - pos.Z();         // Distance to back.

    double result = std::min(std::min(std::min(std::min(std::min(d1, d2), d3), d4), d5), d6);
    return result;
  }

  // Find the closest matching mc trajectory point for a given seed.
  // Returned value is index of the trajectory point.
  // Return -1 in case of no match.
  int mcmatch(const simb::MCParticle& part, const recob::Seed& seed)
  {
    // Get seed point.

    double pos[3];
    double err[3];
    seed.GetPoint(pos, err);

    // Loop over trajectory points.

    int best_traj = -1;
    double max_dist = 0.;
    int ntraj = part.NumberTrajectoryPoints();
    for(int itraj = 0; itraj < ntraj; ++itraj) {
      const TLorentzVector& vec = part.Position(itraj);
      double dx = pos[0] - vec.X();
      double dy = pos[1] - vec.Y();
      double dz = pos[2] - vec.Z();
      double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
      if(best_traj < 0 || dist < max_dist) {
	best_traj = itraj;
	max_dist = dist;
      }
    }
    return best_traj;
  }

  // Length of MC particle.
  //----------------------------------------------------------------------------
  double length(const simb::MCParticle& part, 
		TVector3& start, TVector3& end,
		unsigned int tpc = 0, unsigned int cstat = 0)
  {
    // Get geometry.

    art::ServiceHandle<geo::Geometry> geom;

    // Get fiducial volume boundary.

    double xmin = 0.;
    double xmax = 2.*geom->DetHalfWidth();
    double ymin = -geom->DetHalfHeight();
    double ymax = geom->DetHalfHeight();
    double zmin = 0.;
    double zmax = geom->DetLength();

    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;

    for(int i = 0; i < n; ++i) {
      const TVector3& pos = part.Position(i).Vect();
      if(pos.X() >= xmin &&
	 pos.X() <= xmax &&
	 pos.Y() >= ymin &&
	 pos.Y() <= ymax &&
	 pos.Z() >= zmin &&
	 pos.Z() <= zmax) {
	if(first)
	  start = pos;
	else {
	  disp -= pos;
	  result += disp.Mag();
	}
	first = false;
	disp = pos;
	end = pos;
      }
    }

    return result;
  }

  // Fill efficiency histogram assuming binomial errors.

  void effcalc(const TH1* hnum, const TH1* hden, TH1* heff)
  {
    int nbins = hnum->GetNbinsX();
    assert(nbins == hden->GetNbinsX());
    assert(nbins == heff->GetNbinsX());

    // Loop over bins, including underflow and overflow.

    for(int ibin = 0; ibin <= nbins+1; ++ibin) {
      double num = hnum->GetBinContent(ibin);
      double den = hden->GetBinContent(ibin);
      if(den == 0.) {
	heff->SetBinContent(ibin, 0.);
	heff->SetBinError(ibin, 0.);
      }
      else {
	double eff = num / den;
	if(eff < 0.)
	  eff = 0.;
	if(eff > 1.)
	  eff = 1.;
	double err = std::sqrt(eff * (1.-eff) / den);
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
    assert(nbins == hden->GetNbinsX());
    assert(nbins == hmul->GetNbinsX());

    // Loop over bins, including underflow and overflow.

    for(int ibin = 0; ibin <= nbins+1; ++ibin) {
      double num = hnum->GetBinContent(ibin);
      double den = hden->GetBinContent(ibin);
      if(den == 0.) {
	hmul->SetBinContent(ibin, 0.);
	hmul->SetBinError(ibin, 0.);
      }
      else {
	double mul = num / den;
	if(mul < 0.)
	  mul = 0.;
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

  class SeedAna : public art::EDAnalyzer
  {
  public:

    // Embedded structs.

    // Struct for histograms that depend on seeds only.

    struct RecoHists
    {
      // Constructors.

      RecoHists();
      RecoHists(const std::string& subdir);

      // Pure reco seed histograms.

      TH1F* fHx;           // Seed x position.
      TH1F* fHy;           // Seed y position.
      TH1F* fHz;           // Seed z position.
      TH1F* fHdist;        // Seed distance to boundary.
      TH1F* fHtheta;       // Theta.
      TH1F* fHphi;         // Phi.
      TH1F* fHtheta_xz;    // Theta_xz.
      TH1F* fHtheta_yz;    // Theta_yz.
    };

    // Struct for mc particles and mc-matched tracks.

    struct MCHists
    {
      // Constructors.

      MCHists();
      MCHists(const std::string& subdir);

      // Reco-MC matching.

      TH2F* fHduvcosth;    // 2D mc vs. data matching, duv vs. cos(theta).
      TH1F* fHcosth;       // 1D direction matching, cos(theta).
      TH1F* fHmcu;         // 1D endpoint truth u.
      TH1F* fHmcv;         // 1D endpoint truth v.
      TH1F* fHmcw;         // 1D endpoint truth w.
      TH1F* fHmcdudw;      // Truth du/dw.
      TH1F* fHmcdvdw;      // Truth dv/dw.

      // Pure MC particle histograms (efficiency denominator).

      TH1F* fHmcstartx;    // Starting x position.
      TH1F* fHmcstarty;    // Starting y position.
      TH1F* fHmcstartz;    // Starting z position.
      TH1F* fHmcendx;      // Ending x position.
      TH1F* fHmcendy;      // Ending y position.
      TH1F* fHmcendz;      // Ending z position.
      TH1F* fHmctheta;     // Theta.
      TH1F* fHmcphi;       // Phi.
      TH1F* fHmctheta_xz;  // Theta_xz.
      TH1F* fHmctheta_yz;  // Theta_yz.
      TH1F* fHmcmom;       // Momentum.
      TH1F* fHmclen;       // Length.

      // Matched seed histograms (multiplicity numerator).

      TH1F* fHmstartx;     // Starting x position.
      TH1F* fHmstarty;     // Starting y position.
      TH1F* fHmstartz;     // Starting z position.
      TH1F* fHmendx;       // Ending x position.
      TH1F* fHmendy;       // Ending y position.
      TH1F* fHmendz;       // Ending z position.
      TH1F* fHmtheta;      // Theta.
      TH1F* fHmphi;        // Phi.
      TH1F* fHmtheta_xz;   // Theta_xz.
      TH1F* fHmtheta_yz;   // Theta_yz.
      TH1F* fHmmom;        // Momentum.
      TH1F* fHmlen;        // Length.

      // Matched seed histograms (efficiency numerator).

      TH1F* fHgstartx;     // Starting x position.
      TH1F* fHgstarty;     // Starting y position.
      TH1F* fHgstartz;     // Starting z position.
      TH1F* fHgendx;       // Ending x position.
      TH1F* fHgendy;       // Ending y position.
      TH1F* fHgendz;       // Ending z position.
      TH1F* fHgtheta;      // Theta.
      TH1F* fHgphi;        // Phi.
      TH1F* fHgtheta_xz;   // Theta_xz.
      TH1F* fHgtheta_yz;   // Theta_yz.
      TH1F* fHgmom;        // Momentum.
      TH1F* fHglen;        // Length.

      // Multiplicity histograms.

      TH1F* fHmulstartx;     // Starting x position.
      TH1F* fHmulstarty;     // Starting y position.
      TH1F* fHmulstartz;     // Starting z position.
      TH1F* fHmulendx;       // Ending x position.
      TH1F* fHmulendy;       // Ending y position.
      TH1F* fHmulendz;       // Ending z position.
      TH1F* fHmultheta;      // Theta.
      TH1F* fHmulphi;        // Phi.
      TH1F* fHmultheta_xz;   // Theta_xz.
      TH1F* fHmultheta_yz;   // Theta_yz.
      TH1F* fHmulmom;        // Momentum.
      TH1F* fHmullen;        // Length.

      // Efficiency histograms.

      TH1F* fHestartx;     // Starting x position.
      TH1F* fHestarty;     // Starting y position.
      TH1F* fHestartz;     // Starting z position.
      TH1F* fHeendx;       // Ending x position.
      TH1F* fHeendy;       // Ending y position.
      TH1F* fHeendz;       // Ending z position.
      TH1F* fHetheta;      // Theta.
      TH1F* fHephi;        // Phi.
      TH1F* fHetheta_xz;   // Theta_xz.
      TH1F* fHetheta_yz;   // Theta_yz.
      TH1F* fHemom;        // Momentum.
      TH1F* fHelen;        // Length.
    };

    // Constructors, destructor

    explicit SeedAna(fhicl::ParameterSet const& pset);
    virtual ~SeedAna();

    // Overrides.

    void analyze(const art::Event& evt);
    void endJob();

  private:

    // Fcl Attributes.

    std::string fSeedModuleLabel;
    double fMinMCKE;           // Minimum MC particle kinetic energy (GeV).
    double fMatchColinearity;  // Minimum matching colinearity.
    double fMatchDisp;         // Maximum matching displacement.
    bool fIgnoreSign;          // Ignore sign of mc particle if true.

    // Histograms.

    std::map<int, MCHists> fMCHistMap;       // Indexed by pdg id.
    std::map<int, RecoHists> fRecoHistMap;   // Indexed by pdg id.

    // Statistics.

    int fNumEvent;
  };

  DEFINE_ART_MODULE(SeedAna);

  // RecoHists methods.

  SeedAna::RecoHists::RecoHists() :
    //
    // Purpose: Default constructor.
    //
    fHx(0),
    fHy(0),
    fHz(0),
    fHdist(0),
    fHtheta(0),
    fHphi(0),
    fHtheta_xz(0),
    fHtheta_yz(0)
  {}

  SeedAna::RecoHists::RecoHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.

    *this = RecoHists();

    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("seedana", "SeedAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHx = dir.make<TH1F>("x", "X Position", 100, 0., 2.*geom->DetHalfWidth());
    fHy = dir.make<TH1F>("y", "Y Position", 100, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHz = dir.make<TH1F>("z", "Z Position", 100, 0., geom->DetLength());
    fHdist = dir.make<TH1F>("dist", "Position Distance to Boundary", 
			    100, -10., geom->DetHalfWidth());
    fHtheta = dir.make<TH1F>("theta", "Theta", 100, 0., 3.142);
    fHphi = dir.make<TH1F>("phi", "Phi", 100, -3.142, 3.142);
    fHtheta_xz = dir.make<TH1F>("theta_xz", "Theta_xz", 100, -3.142, 3.142);
    fHtheta_yz = dir.make<TH1F>("theta_yz", "Theta_yz", 100, -3.142, 3.142);
  }

  // MCHists methods.

  SeedAna::MCHists::MCHists() :
    //
    // Purpose: Default constructor.
    //
    fHduvcosth(0),
    fHcosth(0),
    fHmcu(0),
    fHmcv(0),
    fHmcw(0),
    fHmcdudw(0),
    fHmcdvdw(0),
    fHmcstartx(0),
    fHmcstarty(0),
    fHmcstartz(0),
    fHmcendx(0),
    fHmcendy(0),
    fHmcendz(0),
    fHmctheta(0),
    fHmcphi(0),
    fHmctheta_xz(0),
    fHmctheta_yz(0),
    fHmcmom(0),
    fHmclen(0),
    fHmstartx(0),
    fHmstarty(0),
    fHmstartz(0),
    fHmendx(0),
    fHmendy(0),
    fHmendz(0),
    fHmtheta(0),
    fHmphi(0),
    fHmtheta_xz(0),
    fHmtheta_yz(0),
    fHmmom(0),
    fHmlen(0),
    fHgstartx(0),
    fHgstarty(0),
    fHgstartz(0),
    fHgendx(0),
    fHgendy(0),
    fHgendz(0),
    fHgtheta(0),
    fHgphi(0),
    fHgtheta_xz(0),
    fHgtheta_yz(0),
    fHgmom(0),
    fHglen(0),
    fHmulstartx(0),
    fHmulstarty(0),
    fHmulstartz(0),
    fHmulendx(0),
    fHmulendy(0),
    fHmulendz(0),
    fHmultheta(0),
    fHmulphi(0),
    fHmultheta_xz(0),
    fHmultheta_yz(0),
    fHmulmom(0),
    fHmullen(0),
    fHestartx(0),
    fHestarty(0),
    fHestartz(0),
    fHeendx(0),
    fHeendy(0),
    fHeendz(0),
    fHetheta(0),
    fHephi(0),
    fHetheta_xz(0),
    fHetheta_yz(0),
    fHemom(0),
    fHelen(0)
  {}

  SeedAna::MCHists::MCHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.

    *this = MCHists();

    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("seedana", "SeedAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHduvcosth = dir.make<TH2F>("duvcosth", "Delta(uv) vs. Colinearity", 
				100, 0.95, 1., 100, 0., 1.);
    fHcosth = dir.make<TH1F>("colin", "Colinearity", 100, 0.95, 1.);
    fHmcu = dir.make<TH1F>("mcu", "MC Truth U", 100, -1., 1.);
    fHmcv = dir.make<TH1F>("mcv", "MC Truth V", 100, -1., 1.);
    fHmcw = dir.make<TH1F>("mcw", "MC Truth W", 100, -10., 10.);
    fHmcdudw = dir.make<TH1F>("mcdudw", "MC Truth U Slope", 100, -0.2, 0.2);
    fHmcdvdw = dir.make<TH1F>("mcdvdw", "MV Truth V Slope", 100, -0.2, 0.2);

    fHmcstartx = dir.make<TH1F>("mcxstart", "MC X Start Position",
				10, 0., 2.*geom->DetHalfWidth());
    fHmcstarty = dir.make<TH1F>("mcystart", "MC Y Start Position",
				10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmcstartz = dir.make<TH1F>("mczstart", "MC Z Start Position",
				10, 0., geom->DetLength());
    fHmcendx = dir.make<TH1F>("mcxend", "MC X End Position",
			      10, 0., 2.*geom->DetHalfWidth());
    fHmcendy = dir.make<TH1F>("mcyend", "MC Y End Position",
			      10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmcendz = dir.make<TH1F>("mczend", "MC Z End Position",
			      10, 0., geom->DetLength());
    fHmctheta = dir.make<TH1F>("mctheta", "MC Theta", 20, 0., 3.142);
    fHmcphi = dir.make<TH1F>("mcphi", "MC Phi", 10, -3.142, 3.142);
    fHmctheta_xz = dir.make<TH1F>("mctheta_xz", "MC Theta_xz", 40, -3.142, 3.142);
    fHmctheta_yz = dir.make<TH1F>("mctheta_yz", "MC Theta_yz", 40, -3.142, 3.142);
    fHmcmom = dir.make<TH1F>("mcmom", "MC Momentum", 10, 0., 10.);
    fHmclen = dir.make<TH1F>("mclen", "MC Particle Length", 100, 0., 1.1 * geom->DetLength());

    fHmstartx = dir.make<TH1F>("mxstart", "Matched X Start Position",
			       10, 0., 2.*geom->DetHalfWidth());
    fHmstarty = dir.make<TH1F>("mystart", "Matched Y Start Position",
			       10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmstartz = dir.make<TH1F>("mzstart", "Matched Z Start Position",
			       10, 0., geom->DetLength());
    fHmendx = dir.make<TH1F>("mxend", "Matched X End Position",
			     10, 0., 2.*geom->DetHalfWidth());
    fHmendy = dir.make<TH1F>("myend", "Matched Y End Position",
			     10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmendz = dir.make<TH1F>("mzend", "Matched Z End Position",
			     10, 0., geom->DetLength());
    fHmtheta = dir.make<TH1F>("mtheta", "Matched Theta", 20, 0., 3.142);
    fHmphi = dir.make<TH1F>("mphi", "Matched Phi", 10, -3.142, 3.142);
    fHmtheta_xz = dir.make<TH1F>("mtheta_xz", "Matched Theta_xz", 40, -3.142, 3.142);
    fHmtheta_yz = dir.make<TH1F>("mtheta_yz", "Matched Theta_yz", 40, -3.142, 3.142);
    fHmmom = dir.make<TH1F>("mmom", "Matched Momentum", 10, 0., 10.);
    fHmlen = dir.make<TH1F>("mlen", "Matched Particle Length", 100, 0., 1.1 * geom->DetLength());

    fHgstartx = dir.make<TH1F>("gxstart", "Good X Start Position",
			       10, 0., 2.*geom->DetHalfWidth());
    fHgstarty = dir.make<TH1F>("gystart", "Good Y Start Position",
			       10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHgstartz = dir.make<TH1F>("gzstart", "Good Z Start Position",
			       10, 0., geom->DetLength());
    fHgendx = dir.make<TH1F>("gxend", "Good X End Position",
			     10, 0., 2.*geom->DetHalfWidth());
    fHgendy = dir.make<TH1F>("gyend", "Good Y End Position",
			     10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHgendz = dir.make<TH1F>("gzend", "Good Z End Position",
			     10, 0., geom->DetLength());
    fHgtheta = dir.make<TH1F>("gtheta", "Good Theta", 20, 0., 3.142);
    fHgphi = dir.make<TH1F>("gphi", "Good Phi", 10, -3.142, 3.142);
    fHgtheta_xz = dir.make<TH1F>("gtheta_xz", "Good Theta_xz", 40, -3.142, 3.142);
    fHgtheta_yz = dir.make<TH1F>("gtheta_yz", "Good Theta_yz", 40, -3.142, 3.142);
    fHgmom = dir.make<TH1F>("gmom", "Good Momentum", 10, 0., 10.);
    fHglen = dir.make<TH1F>("glen", "Good Particle Length", 100, 0., 1.1 * geom->DetLength());

    fHmulstartx = dir.make<TH1F>("mulxstart", "Multiplicity vs. X Start Position",
				 10, 0., 2.*geom->DetHalfWidth());
    fHmulstarty = dir.make<TH1F>("mulystart", "Multiplicity vs. Y Start Position",
				 10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmulstartz = dir.make<TH1F>("mulzstart", "Multiplicity vs. Z Start Position",
				 10, 0., geom->DetLength());
    fHmulendx = dir.make<TH1F>("mulxend", "Multiplicity vs. X End Position",
			       10, 0., 2.*geom->DetHalfWidth());
    fHmulendy = dir.make<TH1F>("mulyend", "Multiplicity vs. Y End Position",
			       10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmulendz = dir.make<TH1F>("mulzend", "Multiplicity vs. Z End Position",
			       10, 0., geom->DetLength());
    fHmultheta = dir.make<TH1F>("multheta", "Multiplicity vs. Theta", 20, 0., 3.142);
    fHmulphi = dir.make<TH1F>("mulphi", "Multiplicity vs. Phi", 10, -3.142, 3.142);
    fHmultheta_xz = dir.make<TH1F>("multheta_xz", "Multiplicity vs. Theta_xz", 40, -3.142, 3.142);
    fHmultheta_yz = dir.make<TH1F>("multheta_yz", "Multiplicity vs. Theta_yz", 40, -3.142, 3.142);
    fHmulmom = dir.make<TH1F>("mulmom", "Multiplicity vs. Momentum", 10, 0., 10.);
    fHmullen = dir.make<TH1F>("mullen", "Multiplicity vs. Particle Length",
			      100, 0., 1.1 * geom->DetLength());

    fHestartx = dir.make<TH1F>("exstart", "Efficiency vs. X Start Position",
			       10, 0., 2.*geom->DetHalfWidth());
    fHestarty = dir.make<TH1F>("eystart", "Efficiency vs. Y Start Position",
			       10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHestartz = dir.make<TH1F>("ezstart", "Efficiency vs. Z Start Position",
			       10, 0., geom->DetLength());
    fHeendx = dir.make<TH1F>("exend", "Efficiency vs. X End Position",
			     10, 0., 2.*geom->DetHalfWidth());
    fHeendy = dir.make<TH1F>("eyend", "Efficiency vs. Y End Position",
			     10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHeendz = dir.make<TH1F>("ezend", "Efficiency vs. Z End Position",
			     10, 0., geom->DetLength());
    fHetheta = dir.make<TH1F>("etheta", "Efficiency vs. Theta", 20, 0., 3.142);
    fHephi = dir.make<TH1F>("ephi", "Efficiency vs. Phi", 10, -3.142, 3.142);
    fHetheta_xz = dir.make<TH1F>("etheta_xz", "Efficiency vs. Theta_xz", 40, -3.142, 3.142);
    fHetheta_yz = dir.make<TH1F>("etheta_yz", "Efficiency vs. Theta_yz", 40, -3.142, 3.142);
    fHemom = dir.make<TH1F>("emom", "Efficiency vs. Momentum", 10, 0., 10.);
    fHelen = dir.make<TH1F>("elen", "Efficiency vs. Particle Length",
			    100, 0., 1.1 * geom->DetLength());
  }

  SeedAna::SeedAna(const fhicl::ParameterSet& pset)
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    : EDAnalyzer(pset)
    , fSeedModuleLabel(pset.get<std::string>("SeedModuleLabel"))
    , fMinMCKE(pset.get<double>("MinMCKE"))
    , fMatchColinearity(pset.get<double>("MatchColinearity"))
    , fMatchDisp(pset.get<double>("MatchDisp"))
    , fIgnoreSign(pset.get<bool>("IgnoreSign"))
    , fNumEvent(0)
  {

    // Report.

    mf::LogInfo("SeedAna") 
      << "SeedAna configured with the following parameters:\n"
      << "  SeedModuleLabel = " << fSeedModuleLabel << "\n"
      << "  MinMCKE = " << fMinMCKE;
  }

  SeedAna::~SeedAna()
  //
  // Purpose: Destructor.
  //
  {}

  void SeedAna::analyze(const art::Event& evt)
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();

    // Get seed handle.

    art::Handle< std::vector<recob::Seed> > seedh;
    evt.getByLabel(fSeedModuleLabel, seedh);

    // Get mc particles.

    sim::ParticleList plist;
    std::vector<const simb::MCParticle*> plist2;

    if(mc) {
      art::ServiceHandle<cheat::BackTracker> bt;
      plist = bt->ParticleList();

      // Loop over mc particles, and fill histograms that depend only
      // on mc particles.

      for(sim::ParticleList::const_iterator ipart = plist.begin();
	  ipart != plist.end(); ++ipart) {
	const simb::MCParticle* part = (*ipart).second;
	assert(part != 0);
	int pdg = part->PdgCode();
	if(fIgnoreSign)
	  pdg = std::abs(pdg);

	// Ignore everything except stable charged nonshowering particles.

	int apdg = std::abs(pdg);
	if(apdg == 13 ||     // Muon
	   apdg == 211 ||    // Charged pion
	   apdg == 321 ||    // Charged kaon
	   apdg == 2212) {   // (Anti)proton

	  // Apply minimum energy cut.

	  if(part->E() >= 0.001*part->Mass() + fMinMCKE) {

	    // Fill mc histograms (efficiency denominstors).

	    if(fMCHistMap.count(pdg) == 0) {
	      std::ostringstream ostr;
	      ostr << "MC" << (fIgnoreSign ? "All" : (pdg > 0 ? "Pos" : "Neg")) << std::abs(pdg);
	      fMCHistMap[pdg] = MCHists(ostr.str());
	    }
	    const MCHists& mchists = fMCHistMap[pdg];

	    TVector3 mcstart;
	    TVector3 mcend;
	    double plen = length(*part, mcstart, mcend);
	    double mctheta_xz = std::atan2(part->Px(), part->Pz());
	    double mctheta_yz = std::atan2(part->Py(), part->Pz());

	    mchists.fHmcstartx->Fill(mcstart.X());
	    mchists.fHmcstarty->Fill(mcstart.Y());
	    mchists.fHmcstartz->Fill(mcstart.Z());
	    mchists.fHmcendx->Fill(mcend.X());
	    mchists.fHmcendy->Fill(mcend.Y());
	    mchists.fHmcendz->Fill(mcend.Z());
	    mchists.fHmctheta->Fill(part->Momentum().Theta());
	    mchists.fHmcphi->Fill(part->Momentum().Phi());
	    mchists.fHmctheta_xz->Fill(mctheta_xz);
	    mchists.fHmctheta_yz->Fill(mctheta_yz);
	    mchists.fHmcmom->Fill(part->Momentum().Vect().Mag());
	    mchists.fHmclen->Fill(plen);

	    // Loop over seeds and do matching.

	    int nmatch = 0;
	    if(seedh.isValid()) {

	      // Loop over seeds.

	      int nseed = seedh->size();
	      for(int i = 0; i < nseed; ++i) {
		art::Ptr<recob::Seed> pseed(seedh, i);
		const recob::Seed& seed = *pseed;

		// Get parameters of this seed.

		TVector3 pos;
		TVector3 dir;
		double err[3];
		seed.GetPoint(&pos[0], err);
		seed.GetDirection(&dir[0], err);

		// Calculate the global-to-local rotation matrix.
		// Copied from Track.cxx.

		TMatrixD rot(3,3);
		double dirmag = dir.Mag();
		double diryz = std::sqrt(dir.Y()*dir.Y() + dir.Z()*dir.Z());

		double sinth = dir.X() / dirmag;
		double costh = diryz / dirmag;
		double sinphi = 0.;
		double cosphi = 1.;
		if(diryz != 0) {
		  sinphi = -dir.Y() / diryz;
		  cosphi = dir.Z() / diryz;
		}
		rot(0,0) = costh;
		rot(1,0) = 0.;
		rot(2,0) = sinth;
		rot(0,1) = sinth * sinphi;
		rot(1,1) = cosphi;
		rot(2,1) = -costh * sinphi;
		rot(0,2) = -sinth * cosphi;
		rot(1,2) = sinphi;
		rot(2,2) = costh * cosphi;

		// Get best matching mc trajectory point.

		int itraj = mcmatch(*part, seed);
		if(itraj >= 0) {

		  // Get mc relative position and direction at matching trajectory point.

		  TVector3 mcpos = part->Position(itraj).Vect() - pos;
		  TVector3 mcmom = part->Momentum(itraj).Vect();

		  // Rotate the momentum and position to the
		  // seed-local coordinate system.

		  TVector3 mcmoml = rot * mcmom;
		  TVector3 mcposl = rot * mcpos;

		  if(mcmoml.Z() < 0.)
		    mcmoml = -mcmoml;
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
		  double uv0 = std::sqrt(u0*u0 + v0*v0);

		  // Fill matching histograms.

		  mchists.fHduvcosth->Fill(costh, uv0);
		  if(std::abs(uv0) < fMatchDisp) {

		    // Fill slope matching histograms.

		    mchists.fHmcdudw->Fill(dudw);
		    mchists.fHmcdvdw->Fill(dvdw);
		  }
		  mchists.fHcosth->Fill(costh);
		  if(costh > fMatchColinearity) {

		    // Fill displacement matching histograms.

		    mchists.fHmcu->Fill(u0);
		    mchists.fHmcv->Fill(v0);
		    mchists.fHmcw->Fill(w);

		    if(std::abs(uv0) < fMatchDisp) {

		      // Now we have passed all matching cuts and we have a matching
		      // mc particle + seed pair.

		      ++nmatch;

		      // Fill matched seed histograms (seed multiplicity).

		      mchists.fHmstartx->Fill(mcstart.X());
		      mchists.fHmstarty->Fill(mcstart.Y());
		      mchists.fHmstartz->Fill(mcstart.Z());
		      mchists.fHmendx->Fill(mcend.X());
		      mchists.fHmendy->Fill(mcend.Y());
		      mchists.fHmendz->Fill(mcend.Z());
		      mchists.fHmtheta->Fill(part->Momentum().Theta());
		      mchists.fHmphi->Fill(part->Momentum().Phi());
		      mchists.fHmtheta_xz->Fill(mctheta_xz);
		      mchists.fHmtheta_yz->Fill(mctheta_yz);
		      mchists.fHmmom->Fill(part->Momentum().Vect().Mag());
		      mchists.fHmlen->Fill(plen);
		    }
		  }
		}
	      }

	      // If we found at least one matched seed, fill good
	      // particle histograms.

	      if(nmatch > 0) {
		mchists.fHgstartx->Fill(mcstart.X());
		mchists.fHgstarty->Fill(mcstart.Y());
		mchists.fHgstartz->Fill(mcstart.Z());
		mchists.fHgendx->Fill(mcend.X());
		mchists.fHgendy->Fill(mcend.Y());
		mchists.fHgendz->Fill(mcend.Z());
		mchists.fHgtheta->Fill(part->Momentum().Theta());
		mchists.fHgphi->Fill(part->Momentum().Phi());
		mchists.fHgtheta_xz->Fill(mctheta_xz);
		mchists.fHgtheta_yz->Fill(mctheta_yz);
		mchists.fHgmom->Fill(part->Momentum().Vect().Mag());
		mchists.fHglen->Fill(plen);
	      }
	    }
	  }
	}
      }
    }

    // Loop over seeds and fill reco-only seed histograms.

    if(seedh.isValid()) {

      // Loop over seeds.

      int nseed = seedh->size();
      for(int i = 0; i < nseed; ++i) {
	art::Ptr<recob::Seed> pseed(seedh, i);
	const recob::Seed& seed = *pseed;

	// Fill histograms involving reco seeds only.

	TVector3 pos;
	TVector3 dir;
	double err[3];
	seed.GetPoint(&pos[0], err);
	seed.GetDirection(&dir[0], err);

	double dpos = bdist(pos);
	double theta_xz = std::atan2(dir.X(), dir.Z());
	double theta_yz = std::atan2(dir.Y(), dir.Z());

	if(fRecoHistMap.count(0) == 0)
	  fRecoHistMap[0] = RecoHists("Reco");
	const RecoHists& rhists = fRecoHistMap[0];

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

    mf::LogInfo("SeedAna") 
      << "SeedAna statistics:\n"
      << "  Number of events = " << fNumEvent;

    // Fill multiplicity histograms.

    for(std::map<int, MCHists>::const_iterator i = fMCHistMap.begin();
	i != fMCHistMap.end(); ++i) {
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

    for(std::map<int, MCHists>::const_iterator i = fMCHistMap.begin();
	i != fMCHistMap.end(); ++i) {
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

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

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RecoAlg/SpacePointAlg.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "MCCheater/BackTracker.h"

#include "TH1F.h"


namespace trkf {

  class SpacePointAna : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit SpacePointAna(fhicl::ParameterSet const& pset);
    virtual ~SpacePointAna();

    // Book histograms.

    void bookHistograms(bool mc);

    // Overrides.

    void analyze(const art::Event& evt);

  private:

    // Fcl Attributes.

    const SpacePointAlg fSptalgTime;    // Algorithm object (increased time cut).
    const SpacePointAlg fSptalgSep;     // Algorithm object (increased sepataion cut).
    const SpacePointAlg fSptalgDefault; // Algorithm object (default cuts).
    std::string fHitModuleLabel;
    bool fUseClusterHits;
    std::string fClusterModuleLabel;
    bool fUseMC;

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
    TH1F* fHS;       // Spatial separation.
    TH1F* fHchisq;   // Space point chisquare.
    TH1F* fHx;       // X position.
    TH1F* fHy;       // Y position.
    TH1F* fHz;       // Z position.
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
   , fHS(0)
   , fHchisq(0)
   , fHx(0)
   , fHy(0)
   , fHz(0)
   , fHMCdx(0)
   , fHMCdy(0)
   , fHMCdz(0)
   , fHMCxpull(0)
   , fHMCypull(0)
   , fHMCzpull(0)
   , fNumEvent(0)
  {

    // Report.

    mf::LogInfo("SpacePointAna") 
      << "SpacePointAna configured with the following parameters:\n"
      << "  HitModuleLabel = " << fHitModuleLabel << "\n"
      << "  UseClusterHits = " << fUseClusterHits << "\n"
      << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
      << "  UseMC = " << fUseMC;
  }

  SpacePointAna::~SpacePointAna()
  //
  // Purpose: Destructor.
  //
  {}

  void SpacePointAna::bookHistograms(bool mc)
  //
  // Purpose: Book histograms.
  //
  {
    if(!fBooked) {
      fBooked = true;

      art::ServiceHandle<geo::Geometry> geom;
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory dir = tfs->mkdir("sptana", "SpacePointAna histograms");

      if(mc) {
	fHDTUE = dir.make<TH1F>("MCDTUE", "U-Drift Electrons Time Difference", 100, -2., 2.);
	fHDTVE = dir.make<TH1F>("MCDTVE", "V-Drift Electrons Time Difference", 100, -2., 2.);
	fHDTWE = dir.make<TH1F>("MCDTWE", "W-Drift Electrons Time Difference", 100, -2., 2.);
	fHDTUPull = dir.make<TH1F>("MCDTUPull", "U-Drift Electrons Time Pull", 100, -50., 50.);
	fHDTVPull = dir.make<TH1F>("MCDTVPull", "V-Drift Electrons Time Pull", 100, -50., 50.);
	fHDTWPull = dir.make<TH1F>("MCDTWPull", "W-Drift Electrons Time Pull", 100, -50., 50.);
      }
      if(!fSptalgTime.merge()) {
	fHDTUV = dir.make<TH1F>("DTUV", "U-V time difference", 100, -20., 20.);
	fHDTVW = dir.make<TH1F>("DTVW", "V-W time difference", 100, -20., 20.);
	fHDTWU = dir.make<TH1F>("DTWU", "W-U time difference", 100, -20., 20.);
	fHS = dir.make<TH1F>("DS", "Spatial Separatoin", 100, -2., 2.);
      }
      fHchisq = dir.make<TH1F>("chisq", "Chisquare", 100, 0., 20.);

      fHx = dir.make<TH1F>("xpos", "X Position",
			   100, 0., 2.*geom->DetHalfWidth());
      fHy = dir.make<TH1F>("ypos", "Y Position",
			   100, -geom->DetHalfHeight(), geom->DetHalfHeight());
      fHz = dir.make<TH1F>("zpos", "Z Position",
			   100, 0., geom->DetLength());
      if(mc) {
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
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    art::ServiceHandle<cheat::BackTracker> bt;

    ++fNumEvent;

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();
    bookHistograms(mc);

    // Get Services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<util::LArProperties> larprop;

    // Get Hits.

    art::PtrVector<recob::Hit> hits;

    if(fUseClusterHits) {

      // Get clusters.

      art::Handle< std::vector<recob::Cluster> > clusterh;
      evt.getByLabel(fClusterModuleLabel, clusterh);

      // Get hits from all clusters.
      art::FindManyP<recob::Hit> fm(clusterh, evt, fClusterModuleLabel);

      if(clusterh.isValid()) {
	int nclus = clusterh->size();

	for(int i = 0; i < nclus; ++i) {
	  art::Ptr<recob::Cluster> pclus(clusterh, i);
	  std::vector< art::Ptr<recob::Hit> > clushits = fm.at(i);
	  int nhits = clushits.size();
	  hits.reserve(hits.size() + nhits);

	  for(std::vector< art::Ptr<recob::Hit> >::const_iterator ihit = clushits.begin(); ihit != clushits.end(); ++ihit) {
	    hits.push_back(*ihit);
	  }
	}
      }
    }
    else {

      // Get unclustered hits.

      art::Handle< std::vector<recob::Hit> > hith;
      evt.getByLabel(fHitModuleLabel, hith);
      if(hith.isValid()) {
	int nhits = hith->size();
	hits.reserve(nhits);

	for(int i = 0; i < nhits; ++i)
	  hits.push_back(art::Ptr<recob::Hit>(hith, i));
      }
    }

    // Fill histograms that don't depend on space points.

    if(mc) {

      // Loop over hits and fill hit-electron time difference histogram.

      for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	  ihit != hits.end(); ++ihit) {
	const recob::Hit& hit = **ihit;

	//unsigned int channel = hit.Channel();
	//geo::View_t geo_view = geom->View(channel);
	//geo::View_t hit_view = hit.View();
	//assert(geo_view == hit_view);
	double tpeak = hit.PeakTime();
	double terr = hit.SigmaPeakTime();

	//assert(channel == hit.Channel());

	// Loop over electrons associated with this hit/channel and fill
	// hit-electron time difference histograms.

	// loop over the map of TDC to sim::IDE to get the TDC for each energy dep
	// Find the average time in ticks for this hit.

	//double sumw = 0.;
	//double sumt = 0.;

	//const std::map<unsigned short, std::vector<sim::IDE> > &idemap = simchan.TDCIDEMap();
	//std::map<unsigned short, std::vector<sim::IDE> >::const_iterator mitr = idemap.begin();
	//for(mitr = idemap.begin(); mitr != idemap.end(); mitr++) {
	//  double tdc = double((*mitr).first);
	//  if(tdc >= tstart && tdc <= tend) {
	//    const std::vector<sim::IDE>& idevec = (*mitr).second;
	//    for(std::vector<sim::IDE>::const_iterator iide=idevec.begin();
	//	iide != idevec.end(); ++iide) {
	//      const sim::IDE& ide = *iide;
	//      double w = ide.numElectrons;
	//      sumw += w;
	//      sumt += w*tdc;
	//    }
	//  }
	//}
	//double tav = 0.;
	//if(sumw != 0.)
	//  tav = sumt / sumw;

	bool tav_ok = true;
	double tav = 0.;
	try {
	  std::vector<double> hitxyz = bt->HitToXYZ(*ihit);
	  tav = detprop->ConvertXToTicks(hitxyz[0], (*ihit)->WireID().Plane, (*ihit)->WireID().TPC, (*ihit)->WireID().Cryostat);
	}
	catch(cet::exception& x) {
	  tav_ok = false;
	}
	if(tav_ok) {
	  if((*ihit)->View() == geo::kU) {
	    fHDTUE->Fill(tpeak - tav);
	    fHDTUPull->Fill((tpeak - tav) / terr);
	  }
	  else if((*ihit)->View() == geo::kV) {
	    fHDTVE->Fill(tpeak - tav);
	    fHDTVPull->Fill((tpeak - tav) / terr);
	  }
	  else if((*ihit)->View() == geo::kZ) {
	    fHDTWE->Fill(tpeak - tav);
	    fHDTWPull->Fill((tpeak - tav) / terr);
	  }
	  else
	    throw cet::exception("SpacePointAna") << "Bad view = " << (*ihit)->View() << "\n";
	}
      }
    }
    
    std::vector<recob::SpacePoint> spts1;  // For time histograms.
    std::vector<recob::SpacePoint> spts2;  // For separation histogram.
    std::vector<recob::SpacePoint> spts3;  // Default cuts.

    // If nonzero time cut is specified, make space points using that
    // time cut (for time histograms).

    if(!fSptalgTime.merge()) {
      if(mc && fUseMC)
	fSptalgTime.makeMCTruthSpacePoints(hits, spts1);
      else
	fSptalgTime.makeSpacePoints(hits, spts1);

      // Report number of space points.

      mf::LogDebug("SpacePointAna") << "Found " << spts1.size() 
				    << " space points using special time cut.";
    }

    // If nonzero separation cut is specified, make space points using that 
    // separation cut (for separation histogram).

    if(!fSptalgSep.merge()) {
      if(mc && fUseMC)
	fSptalgSep.makeMCTruthSpacePoints(hits, spts2);
      else
	fSptalgSep.makeSpacePoints(hits, spts2);

      // Report number of space points.

      mf::LogDebug("SpacePointAna") << "Found " << spts2.size() 
				    << " space points using special seperation cut.";
    }

    // Make space points using default cuts.

    if(mc && fUseMC)
      fSptalgDefault.makeMCTruthSpacePoints(hits, spts3);
    else
      fSptalgDefault.makeSpacePoints(hits, spts3);

    // Report number of space points.

    mf::LogDebug("SpacePointAna") << "Found " << spts3.size() 
				  << " space points using default cuts.";

    if(!fSptalgTime.merge()) {

      // Loop over space points and fill time histograms.

      for(std::vector<recob::SpacePoint>::const_iterator i = spts1.begin(); 
	  i != spts1.end(); ++i) {
	const recob::SpacePoint& spt = *i;

	// Get hits associated with this SpacePoint.

	const art::PtrVector<recob::Hit>& spthits = fSptalgTime.getAssociatedHits(spt);

	// Make a double loop over hits and fill hit time difference histograms.

	for(art::PtrVector<recob::Hit>::const_iterator ihit = spthits.begin();
	    ihit != spthits.end(); ++ihit) {
	  const recob::Hit& hit1 = **ihit;

	  geo::WireID hit1WireID = hit1.WireID();
	  unsigned int tpc1, plane1;
	  tpc1 = hit1WireID.TPC;
	  plane1 = hit1WireID.Plane;
		
	  geo::View_t view1 = hit1.View();
	  double t1 = fSptalgTime.correctedTime(hit1);

	  for(art::PtrVector<recob::Hit>::const_iterator jhit = spthits.begin();
	      jhit != spthits.end(); ++jhit) {
	    const recob::Hit& hit2 = **jhit;

	    geo::WireID hit2WireID = hit2.WireID();
	    unsigned int tpc2, plane2;
	    tpc2 = hit2WireID.TPC;
	    plane2 = hit2WireID.Plane;

	    // Require same tpc, different view.

	    if(tpc1 == tpc2 && plane1 != plane2) {

	      geo::View_t view2 = hit2.View();
	      double t2 = fSptalgTime.correctedTime(hit2);

	      if(view1 == geo::kU) {
		if(view2 == geo::kV)
		  fHDTUV->Fill(t1-t2);
		if(view2 == geo::kZ)
		  fHDTWU->Fill(t2-t1);
	      }
	      if(view1 == geo::kV) {
		if(view2 == geo::kZ)
		  fHDTVW->Fill(t1-t2);
		if(view2 == geo::kU)
		  fHDTUV->Fill(t2-t1);
	      }
	      if(view1 == geo::kZ) {
		if(view2 == geo::kU)
		  fHDTWU->Fill(t1-t2);
		if(view2 == geo::kV)
		  fHDTVW->Fill(t2-t1);
	      }
	    }
	  }
	}
      }

      // Loop over space points and fill seperation histograms.

      for(std::vector<recob::SpacePoint>::const_iterator i = spts2.begin(); 
	  i != spts2.end(); ++i) {
	const recob::SpacePoint& spt = *i;

	// Get hits associated with this SpacePoint.

	const art::PtrVector<recob::Hit>& spthits = fSptalgSep.getAssociatedHits(spt);

	// Fill separation histogram.

	double sep = fSptalgSep.separation(spthits);
	fHS->Fill(sep);
      }
    }

    // Loop over default space points and fill histograms.

    for(std::vector<recob::SpacePoint>::const_iterator i = spts3.begin();
	i != spts3.end(); ++i) {
      const recob::SpacePoint& spt = *i;

      fHchisq->Fill(spt.Chisq());
      fHx->Fill(spt.XYZ()[0]);
      fHy->Fill(spt.XYZ()[1]);
      fHz->Fill(spt.XYZ()[2]);
      if(mc) {
	try {
	  std::vector<double> mcxyz = bt->SpacePointHitsToXYZ(fSptalgDefault.getAssociatedHits(spt));
	  fHMCdx->Fill(spt.XYZ()[0] - mcxyz[0]);
	  fHMCdy->Fill(spt.XYZ()[1] - mcxyz[1]);
	  fHMCdz->Fill(spt.XYZ()[2] - mcxyz[2]);
	  if(spt.ErrXYZ()[0] > 0.)
	    fHMCxpull->Fill((spt.XYZ()[0] - mcxyz[0]) / std::sqrt(spt.ErrXYZ()[0]));
	  if(spt.ErrXYZ()[2] > 0.)
	    fHMCypull->Fill((spt.XYZ()[1] - mcxyz[1]) / std::sqrt(spt.ErrXYZ()[2]));
	  if(spt.ErrXYZ()[5] > 0.)
	    fHMCzpull->Fill((spt.XYZ()[2] - mcxyz[2]) / std::sqrt(spt.ErrXYZ()[5]));
	}
	catch(cet::exception& x) {}
      }
    }
  }
}

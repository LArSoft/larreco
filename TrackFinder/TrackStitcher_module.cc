////////////////////////////////////////////////////////////////////////
//
// \file TrackStitcher
//
// \author echurch@fnal.gov
//
//  This algorithm is designed to join tracks that point in roughly same direction
//  and whose endpoints are suitably close.
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <tuple>

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Principal/View.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Track.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoAlg/StitchAlg.h"
#include "Simulation/sim.h"
#include "SimulationBase/MCTruth.h"
#include "Utilities/AssociationUtil.h"

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include <TTree.h>
#include <TMatrixT.h>


#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <vector>
#include <string>


// ROOT includes
#include "TVectorD.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TDatabasePDG.h"
#include "Utilities/AssociationUtil.h"

class StitchAlg;

namespace trkf {

  class TrackStitcher : public art::EDProducer {
    
  public:
    
    explicit TrackStitcher(fhicl::ParameterSet const& pset);
    virtual ~TrackStitcher();
    
    //////////////////////////////////////////////////////////
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:

    const art::PtrVector<recob::Hit> GetHitsFromComponentTracks(const art::PtrVector<recob::Track> &, const art::Event& evt);
    const art::PtrVector<recob::SpacePoint> GetSpacePointsFromComponentTracks(const art::PtrVector<recob::Track> &, const art::Event& evt);
    const art::PtrVector<recob::Hit> GetHitsFromAssdSpacePoints(const art::PtrVector<recob::SpacePoint> &, const art::Event& evt);
    std::string     fTrackModuleLabel;// label for input collection
    std::string     fSpptModuleLabel;// label for input collection
    bool            fStizatch; // CommonComponentStitch
    StitchAlg fStitchAlg;

  }; // class TrackStitcher

} // end namespace for declarations

namespace trkf {

//-------------------------------------------------
  TrackStitcher::TrackStitcher(fhicl::ParameterSet const& pset) : 
    fStitchAlg(pset.get< fhicl::ParameterSet >("StitchAlg"))
  {
    
    this->reconfigure(pset);
    
    produces< std::vector<recob::Track>                  >();
    produces<std::vector<art::PtrVector<recob::Track> >  >(); 
    produces<art::Assns<recob::Track, recob::Hit>        >();
    produces<art::Assns<recob::Track, recob::SpacePoint> >();
    produces<art::Assns<recob::SpacePoint, recob::Hit>   >();
    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine( seed );

  }

//-------------------------------------------------
  void TrackStitcher::reconfigure(fhicl::ParameterSet const& pset) 
  {
    fTrackModuleLabel    = pset.get< std::string >("TrackModuleLabel");
    fSpptModuleLabel     = pset.get< std::string >("SpptModuleLabel"); 
    fStizatch            = pset.get< bool >       ("CommonComponentStitch",true); 
    fStitchAlg.reconfigure(pset.get< fhicl::ParameterSet >("StitchAlg"));
  }
  
  //-------------------------------------------------
  TrackStitcher::~TrackStitcher()
  {
  }


//-------------------------------------------------
  void TrackStitcher::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
  }

//-------------------------------------------------
  void TrackStitcher::endJob()
  {
  }


//------------------------------------------------------------------------------------//
  void TrackStitcher::produce(art::Event& evt)
  { 

    // get services
    art::ServiceHandle<geo::Geometry> geom;

    //////////////////////////////////////////////////////
    // Make a std::unique_ptr<> for the thing you want to put into the event
    //////////////////////////////////////////////////////
    // tcol is the collection of new tracks
    std::unique_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);
    std::unique_ptr<std::vector<recob::SpacePoint> > scol(new std::vector<recob::SpacePoint>);
    // tvcol is the collection of vectors that comprise each tcol
    std::unique_ptr<std::vector< art::PtrVector<recob::Track> > > tvcol(new std::vector< art::PtrVector<recob::Track> >);
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> > thassn(new art::Assns<recob::Track, recob::Hit>);     
    std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint> > tsptassn(new art::Assns<recob::Track, recob::SpacePoint>);     
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit > > spthassn(new art::Assns<recob::SpacePoint, recob::Hit>);     


    // Get the original Spacepoints.
    art::Handle< std::vector<recob::SpacePoint> > sppth;
    evt.getByLabel(fSpptModuleLabel, sppth);
    for (size_t ii=0; ii<sppth->size() ;ii++)
      scol->push_back(sppth->at(ii));


    // Find the best match for each track's Head and Tail.
    fStitchAlg.FindHeadsAndTails(evt,fTrackModuleLabel);
    // walk through each vertex of one track to its match on another, and so on and stitch 'em.
    fStitchAlg.WalkStitch();
    // search composite tracks and stitch further if there are components in common. Do it until all are stitched.
    bool stizatch(fStizatch);
    while (stizatch)
      {
	stizatch = fStitchAlg.CommonComponentStitch();
      }
    mf::LogVerbatim("TrackStitcher.beginning") << "There are " <<  fStitchAlg.ftListHandle->size() << " Tracks in this event before stitching.";

    fStitchAlg.GetTracks(*tcol);    
    fStitchAlg.GetTrackComposites(*tvcol);

    if (tcol->size()!=tvcol->size())
      throw cet::exception("TrackStitcher") << "Tracks and TrackComposites do not match: "<<tcol->size()<<" vs "<<tvcol->size()<<"\n";      

    for (size_t ii=0; ii<tvcol->size(); ii++)
      {
	const art::PtrVector<recob::Hit>& hits(GetHitsFromComponentTracks(tvcol->at(ii), evt));
	// Now make the Assns of relevant Hits to stitched Track
	util::CreateAssn(*this, evt, *tcol, hits, *thassn, ii);
	const art::PtrVector<recob::SpacePoint>& sppts(GetSpacePointsFromComponentTracks(tvcol->at(ii), evt));
	// Now make the Assns of relevant Sppts to stitched Track
	util::CreateAssn(*this, evt, *tcol, sppts, *tsptassn, ii);
	// Now Assns of sppts to hits

	const art::PtrVector<recob::Hit>& hitsFromSppts(GetHitsFromAssdSpacePoints(sppts, evt));
	for ( size_t jj=0; jj<sppts.size(); jj++ ) 
	  {
	    // find jjth sppt in the list of scol. Meaning, find kkth element of sppth.
	    size_t ll(2e5);
	    for ( size_t kk=0; kk<scol->size(); kk++ ) 
	      {
		art::Ptr<recob::SpacePoint> sppthOne (sppth,kk);
		if ( sppthOne != sppts.at(jj)) continue;
		ll = kk;
	      }
	    if (ll!=2e5)
	      util::CreateAssn(*this, evt, *scol, hitsFromSppts, *spthassn, ll);
	  }
      }


    mf::LogVerbatim("TrackStitcher.end") << "There are " <<  tvcol->size() << " Tracks in this event after stitching.";
    evt.put(std::move(tcol)); 
    evt.put(std::move(tvcol));
    // Add Hit-to-Track and Sppt-to-Track Assns.
    evt.put(std::move(thassn));
    evt.put(std::move(tsptassn));
    evt.put(std::move(spthassn));

  }
  
  const art::PtrVector<recob::Hit> TrackStitcher::GetHitsFromComponentTracks(const art::PtrVector<recob::Track> &tcomp, const art::Event& evtGHFCT) 
  {

    art::PtrVector<recob::Hit> hits;
    art::FindManyP<recob::Hit> hitAssns(tcomp, evtGHFCT, fTrackModuleLabel); 

    for (unsigned int ii=0; ii < tcomp.size(); ++ii )
      {
	 hits.insert(hits.end(),hitAssns.at(ii).begin(), hitAssns.at(ii).end());
      }



    //    const art::PtrVector<recob::Hit> chits(hits);
    return hits;
  }

  const art::PtrVector<recob::SpacePoint> TrackStitcher::GetSpacePointsFromComponentTracks(const art::PtrVector<recob::Track> &tcomp, const art::Event& evtGHFCT) 
  {

    art::PtrVector<recob::SpacePoint> sppts;
    art::FindManyP<recob::SpacePoint> spptAssns(tcomp, evtGHFCT, fTrackModuleLabel); 
    for (unsigned int ii=0; ii < tcomp.size(); ++ii )
      {
	sppts.insert(sppts.end(),spptAssns.at(ii).begin(), spptAssns.at(ii).end());
      }
    
    //    const art::PtrVector<recob::Hit> chits(hits);
    return sppts;
  }

  const art::PtrVector<recob::Hit> TrackStitcher::GetHitsFromAssdSpacePoints(const art::PtrVector<recob::SpacePoint> &sppts, const art::Event& evtGHFCT) 
  {

    art::PtrVector<recob::Hit> hits;
    art::FindManyP<recob::Hit> hitAssns(sppts, evtGHFCT, fSpptModuleLabel); 

    for (unsigned int ii=0; ii < sppts.size(); ++ii )
      {
	 hits.insert(hits.end(),hitAssns.at(ii).begin(), hitAssns.at(ii).end());
      }

    return hits;
  }

  DEFINE_ART_MODULE(TrackStitcher)


}  // end namespace

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

#include "Simulation/sim.h"
#include "SimulationBase/MCTruth.h"
#include "Utilities/AssociationUtil.h"

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
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

// Leave the sorting functions here, in case we want them later.
/*
static bool sp_sort_3dz(const art::Ptr<recob::SpacePoint>& h1, const art::Ptr<recob::SpacePoint>& h2)
{
  const double* xyz1 = h1->XYZ();
  const double* xyz2 = h2->XYZ();
  return xyz1[2] < xyz2[2];
}
*/

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

    const recob::Track Stitch(const art::PtrVector<recob::Track> &);
    const art::PtrVector<recob::Hit> GetHitsFromComponentTracks(const art::PtrVector<recob::Track> &, const art::Event& evt);
    const art::PtrVector<recob::SpacePoint> GetSpacePointsFromComponentTracks(const art::PtrVector<recob::Track> &, const art::Event& evt);
    std::string     fTrackModuleLabel;// label for input collection
    std::string     fSpptModuleLabel;// label for input collection
    double fCosAngTol;
    double fSepTol;

    int ftNo;

  }; // class TrackStitcher

} // end namespace for declarations

namespace trkf {

//-------------------------------------------------
  TrackStitcher::TrackStitcher(fhicl::ParameterSet const& pset) 
    : ftNo(0)
  {
    
    this->reconfigure(pset);
    
    produces< std::vector<recob::Track>                  >();
    produces<std::vector<art::PtrVector<recob::Track> >  >(); 
    produces<art::Assns<recob::Track, recob::Hit>        >();
    produces<art::Assns<recob::Track, recob::SpacePoint> >();
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
    fCosAngTol           = pset.get< double >("CosAngTolerance", 0.95); 
    fSepTol              = pset.get< double >("SpptSepTolerance", 10.0); //cm 
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
    // tvcol is the collection of vectors that comprise each tcol
    std::unique_ptr<std::vector< art::PtrVector<recob::Track> > > tvcol(new std::vector< art::PtrVector<recob::Track> >);
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> > thassn(new art::Assns<recob::Track, recob::Hit>);     
    std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint> > tsptassn(new art::Assns<recob::Track, recob::SpacePoint>);     

    art::PtrVector <recob::Track> tcolTmp;
    // define TPC parameters
    TString tpcName = geom->GetLArTPCVolumeName();

    art::Handle< std::vector< recob::Track > > tListHandle;
    evt.getByLabel(fTrackModuleLabel,tListHandle);

    //    mf::LogVerbatim("TrackStitcher.beginning") << "There are " <<  tListHandle->size() << " Tracks in this event before stitching.";
    
    int ntrack = tListHandle->size();
    for(int ii = 0; ii < ntrack; ++ii) {
      art::Ptr<recob::Track> ptrack1(tListHandle, ii);
            const recob::Track& track1 = *ptrack1;
      const TVector3 start1(track1.Vertex());
      const TVector3 end1(track1.End());
      const TVector3 start1Dir(track1.VertexDirection());
      const TVector3 end1Dir(track1.EndDirection());

      // If we don't already have track1, start a new tcolTmp.
      // determination of which requires a loop over all tracks in all vectors
      // and also the current list of tracks, not yet put in a vector.
      bool alreadyHaveIt(false);
      for (std::vector< art::PtrVector<recob::Track> >::const_iterator  itv=tvcol->begin(); itv<tvcol->end() && !alreadyHaveIt; ++itv ) 
	{
	  art::PtrVector<recob::Track> tvcolTmpInner(*itv);
	  for (std::vector< art::Ptr<recob::Track> >::const_iterator  it=tvcolTmpInner.begin(); it<tvcolTmpInner.end() && !alreadyHaveIt; ++it ) 
	    {
	      if ((*it).get() == ptrack1.get()) 
		{
		  alreadyHaveIt=true;
		}
	    }
	}

      for (std::vector< art::Ptr<recob::Track> >::const_iterator  it=tcolTmp.begin(); it<tcolTmp.end() && !alreadyHaveIt; ++it ) 
	{
	  if ((*it).get() == ptrack1.get()) 
	    {
	      alreadyHaveIt=true;
	    }
	}
      

      if (!alreadyHaveIt)
	{
	  if (tcolTmp.size()!=0) 
	    {
	      // I've finished building a vector. push it back onto tvcol, 
	      // make a Stitched vector.
	      tvcol->push_back(tcolTmp);
	      const recob::Track& t = Stitch(tcolTmp);
	      // also make the cumulative track
	      tcol->push_back(t);

	      const art::PtrVector<recob::Hit>& hits(GetHitsFromComponentTracks(tcolTmp, evt));
	      // Now make the Assns of relevant Hits to stitched Track
	      util::CreateAssn(*this, evt, *tcol, hits, *thassn, tcol->size()-1);

	      const art::PtrVector<recob::SpacePoint>& sppts(GetSpacePointsFromComponentTracks(tcolTmp, evt));
	      // Now make the Assns of relevant Spacepoints to stitched Track
	      util::CreateAssn(*this, evt, *tcol, sppts, *tsptassn, tcol->size()-1);

	      tcolTmp.erase(tcolTmp.begin(),tcolTmp.end());

	    }
	  // start with 1st elment of a new vector of tracks.
	  tcolTmp.push_back(ptrack1);
	}


      for(int jj = ii+1; jj < ntrack; ++jj) {
	art::Ptr<recob::Track> ptrack2(tListHandle, jj);
	const recob::Track& track2 = *ptrack2;
	const TVector3& start2(track2.Vertex());
	const TVector3& end2(track2.End());
	const TVector3& start2Dir(track2.VertexDirection());
	const TVector3& end2Dir(track2.EndDirection());
	/*
	std::cout << "abs(start1Dir.Dot(end2Dir)) " << std::abs(start1Dir.Dot(end2Dir)) << ", start1-end2.Mag(): " << (start1-end2).Mag() << std::endl;
	std::cout << "abs(end1Dir.Dot(start2Dir)) " << std::abs(end1Dir.Dot(start2Dir)) << ", start2-end1.Mag(): " << (start2-end1).Mag() << std::endl;
	std::cout << "abs(start1Dir.Dot(start2Dir)) " << std::abs(start1Dir.Dot(start2Dir)) << ", start1-start2.Mag(): " << (start1-start2).Mag() << std::endl;
	std::cout << "abs(end1Dir.Dot(end2Dir)) " << std::abs(end1Dir.Dot(end2Dir)) << ", end1-end2.Mag(): " << (end1-end2).Mag() << std::endl;
	*/

	if (
	    (std::abs(start1Dir.Dot(end2Dir))>fCosAngTol && ((start1-end2).Mag())<fSepTol ) ||
	    (std::abs(end1Dir.Dot(start2Dir))>fCosAngTol && ((start2-end1).Mag())<fSepTol ) ||
	    (std::abs(start1Dir.Dot(start2Dir))>fCosAngTol && ((start1-start2).Mag())<fSepTol ) ||
	    (std::abs(end1Dir.Dot(end2Dir))>fCosAngTol && ((end1-end2).Mag())<fSepTol ) 
	    )
	  {
	    tcolTmp.push_back(ptrack2);
	    break;
	  }
	
      } // jj
      
    } // ii

    if (tcolTmp.size()!=0) 
      {
	tvcol->push_back(tcolTmp);
	const recob::Track& t = Stitch(tcolTmp);
	tcol->push_back(t);
	const art::PtrVector<recob::Hit>& hits(GetHitsFromComponentTracks(tcolTmp, evt));
	// Now make the Assns of relevant Hits to stitched Track
	util::CreateAssn(*this, evt, *tcol, hits, *thassn, tcol->size()-1);
	const art::PtrVector<recob::SpacePoint>& sppts(GetSpacePointsFromComponentTracks(tcolTmp, evt));
	// Now make the Assns of relevant Sppts to stitched Track
	util::CreateAssn(*this, evt, *tcol, sppts, *tsptassn, tcol->size()-1);

      }
    
    //        mf::LogVerbatim("TrackStitcher.end") << "There are " <<  tvcol->size() << " Tracks in this event after stitching.";
    evt.put(std::move(tcol)); 
    evt.put(std::move(tvcol));
    // Add Hit-to-Track and Sppt-to-Track Assns.
    evt.put(std::move(thassn));
    evt.put(std::move(tsptassn));

  }
  

  const recob::Track TrackStitcher::Stitch(const art::PtrVector<recob::Track> &tv)
  {

    // take the vector of tracks, walk through each track's vectors of xyz, dxdydz, etc 
    // and concatenate them into longer vectors. Use those to instantiate one new 
    // Stitched-together track.
    std::vector<TVector3> xyz;
    std::vector<TVector3> dxdydz;
    std::vector<TMatrixT<double> > cov;
    std::vector<double> mom;
    std::vector< std::vector <double> > dQdx;
    //art::PtrVector<recob::Track>::const_iterator
    for (auto it = tv.begin(); it!=tv.end(); ++it)
      {
	for (size_t pt = 0; pt!=(*it).get()->NumberTrajectoryPoints(); pt++)
	  {
	    try 
	      { 
		xyz.push_back((*it).get()->LocationAtPoint(pt));
		dxdydz.push_back((*it).get()->DirectionAtPoint(pt));
		TMatrixT<double>  dumc(5,5); 
		if (pt<(*it).get()->NumberCovariance())
		  dumc = (*it).get()->CovarianceAtPoint(pt);
		cov.push_back(dumc);
		double dumm(0.0); 
		if (pt<(*it).get()->NumberFitMomentum())
		  dumm = (*it).get()->MomentumAtPoint(pt);
		mom.push_back(dumm);
		std::vector <double> dum; 
		if (pt<(*it).get()->NumberdQdx(geo::kZ))
		  dum.push_back((*it).get()->DQdxAtPoint(pt,geo::kZ));
		else
		  dum.push_back(0.0);
		dQdx.push_back(dum);

	      }
	    catch (cet::exception &e)
	      {
		mf::LogVerbatim("TrackStitcher bailing. ") << " One or more of xyz, dxdydz, cov, mom, dQdx elements from original Track is out of range..." << e.what() << __LINE__;
		break;
	      }
	  }
      }
   
    /// TODO: sort according to spacepoint distance separation.
    /// As is, we're not sure we're not forming a stitched track with a (some) 
    /// jump(s) and a reversal(s) of direction in it.

    /// But for now we'll instantiate the stitched track faithfully.

    const recob::Track t(xyz,dxdydz,cov,dQdx,mom,ftNo++);
    return t;
  }

  const art::PtrVector<recob::Hit> TrackStitcher::GetHitsFromComponentTracks(const art::PtrVector<recob::Track> &tcomp, const art::Event& evtGHFCT) 
  {

    art::PtrVector<recob::Hit> hits;
    for (unsigned int ii=0; ii < tcomp.size(); ++ii )
      {
	// From the component tracks, get the Hits from the Assns
	const art::Ptr<recob::Track> ptrack( tcomp.at(ii) );
	auto p { ptrack };
	art::FindManyP<recob::Hit> hitAssns(p, evtGHFCT, fTrackModuleLabel); 
	for (unsigned int jj=0; jj < hitAssns.at(0).size(); ++jj)
	  hits.push_back(hitAssns.at(0).at(jj));
      }
    
    //    const art::PtrVector<recob::Hit> chits(hits);
    return hits;
  }

  const art::PtrVector<recob::SpacePoint> TrackStitcher::GetSpacePointsFromComponentTracks(const art::PtrVector<recob::Track> &tcomp, const art::Event& evtGHFCT) 
  {

    art::PtrVector<recob::SpacePoint> sppts;
    for (unsigned int ii=0; ii < tcomp.size(); ++ii )
      {
	// From the component tracks, get the Hits from the Assns
	const art::Ptr<recob::Track> ptrack( tcomp.at(ii) );
	auto p { ptrack };
	art::FindManyP<recob::SpacePoint> spptAssns(p, evtGHFCT, fTrackModuleLabel); 
	for (unsigned int jj=0; jj < spptAssns.at(0).size(); ++jj)
	  sppts.push_back(spptAssns.at(0).at(jj));
      }
    
    //    const art::PtrVector<recob::Hit> chits(hits);
    return sppts;
  }

  DEFINE_ART_MODULE(TrackStitcher)


}  // end namespace

////////////////////////////////////////////////////////////////////////
//
// TrackCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

#include <string>


// ROOT includes
#include "TVector3.h"

// LArSoft includes
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace trkf {
  class TrackCheater : public art::EDProducer {
  public:
    explicit TrackCheater(fhicl::ParameterSet const& pset);

 private:
    void produce(art::Event& evt) override;


    std::string fCheatedClusterLabel; ///< label for module creating recob::Cluster objects
    std::string fG4ModuleLabel;       ///< label for module running G4 and making particles, etc

  };
}

namespace trkf{

  //--------------------------------------------------------------------
  TrackCheater::TrackCheater(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
  {
    fCheatedClusterLabel = pset.get< std::string >("CheatedClusterLabel", "cluster" );
    fG4ModuleLabel       = pset.get< std::string >("G4ModuleLabel",       "largeant");

    produces< std::vector<recob::Track>                   >();
    produces< std::vector<recob::SpacePoint>              >();
    produces< art::Assns<recob::Track, recob::Cluster>    >();
    produces< art::Assns<recob::Track, recob::SpacePoint> >();
    produces< art::Assns<recob::Track, recob::Hit>        >();
  }

  //--------------------------------------------------------------------
  void TrackCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<geo::Geometry const>            geo;

    // grab the clusters that have been reconstructed
    art::Handle< std::vector<recob::Cluster> > clustercol;
    evt.getByLabel(fCheatedClusterLabel, clustercol);

    art::FindManyP<recob::Hit> fmh(clustercol, evt, fCheatedClusterLabel);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Cluster> > clusters;
    art::fill_ptr_vector(clusters, clustercol);

    // loop over the clusters and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Cluster> >::iterator itr = clusters.begin();

    // make a map of vectors of art::Ptrs keyed by eveID values
    std::map< int, std::vector< std::pair<size_t, art::Ptr<recob::Cluster> > > > eveClusterMap;

    // loop over all clusters and fill in the map
    for(size_t c = 0; c < clusters.size(); ++c){

      std::pair<size_t, art::Ptr<recob::Cluster> > idxPtr(c, clusters[c]);

      // in the ClusterCheater module we set the cluster ID to be
      // the eve particle track ID*1000 + plane*100 + tpc*10 + cryostat number.  The
      // floor function on the cluster ID / 1000 will give us
      // the eve track ID
      int eveID = floor((*itr)->ID()/1000.);

      eveClusterMap[eveID].push_back(idxPtr);
      ++itr;
    }// end loop over clusters

    // loop over the map and make prongs
    std::unique_ptr< std::vector<recob::Track> >                   trackcol(new std::vector<recob::Track>);
    std::unique_ptr< std::vector<recob::SpacePoint> >              spcol  (new std::vector<recob::SpacePoint>);
    std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
    std::unique_ptr< art::Assns<recob::Track, recob::Cluster> >    tcassn (new art::Assns<recob::Track, recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> >        thassn (new art::Assns<recob::Track, recob::Hit>);

    for(auto const& clusterMapItr : eveClusterMap){

      // separate out the hits for each particle into the different views
      std::vector< std::pair<size_t, art::Ptr<recob::Cluster> > > const& eveClusters = clusterMapItr.second;

      simb::MCParticle *part = pi_serv->ParticleList()[clusterMapItr.first];

      // is this prong electro-magnetic in nature or
      // hadronic/muonic?  EM --> shower, everything else is a track
      if( abs(part->PdgCode()) != 11  &&
	  abs(part->PdgCode()) != 22  &&
	  abs(part->PdgCode()) != 111 ){

	// vectors to hold the positions and directions of the track
	std::vector<TVector3> points;
	std::vector<TVector3> moms;

	// size_t nviews = geo->Nviews();
	// std::vector< std::vector<double> > dQdx(nviews);

	mf::LogInfo("TrackCheater") << "G4 id " << clusterMapItr.first
				    << " is a track with pdg code "
				    << part->PdgCode();

	art::PtrVector<recob::Cluster> ptrvs;
	std::vector<size_t> idxs;
	for(auto const& idxPtr : eveClusters){
	  ptrvs.push_back(idxPtr.second);
	  idxs.push_back(idxPtr.first);
	}

	// grab the hits associated with the collection plane
	std::vector< art::Ptr<recob::Hit> > hits;
	for(size_t p = 0; p < ptrvs.size(); ++p){
	  std::vector< art::Ptr<recob::Hit> > chits = fmh.at(idxs[p]);
	  if (!chits.size()) continue;
	  if (chits[0]->SignalType() != geo::kCollection) continue;
	  hits.insert(hits.end(), chits.begin(), chits.end());
	}

	// need at least 2 hits to make a track
	if(hits.size() < 2) continue;

	// loop over the hits to get the positions and directions
	size_t spStart = spcol->size();
	for(size_t t = 0; t < hits.size(); ++t){
	  std::vector<double> xyz = bt_serv->HitToXYZ(hits[t]);
	  TVector3 point(xyz[0], xyz[1], xyz[2]);
	  points.push_back(point);

	  std::vector<double> xyz1;
	  //double charge = hits[t]->Integral();
	  double dx     = 0.;
	  double sign   = 1.;

	  if(t < hits.size()-1){
	    xyz1 = bt_serv->HitToXYZ(hits[t+1]);
	  }
	  else{
	    xyz1 = bt_serv->HitToXYZ(hits[t-1]);
	    sign = -1.;
	  }

	  // dx is always positive
	  dx = std::sqrt(std::pow(xyz1[0] - xyz[0], 2) +
			 std::pow(xyz1[1] - xyz[1], 2) +
			 std::pow(xyz1[2] - xyz[2], 2));

	  // figure out momentum
	  double mom = 0;
	  double drmin = std::numeric_limits<double>::max();
	  for (unsigned int itp = 0; itp<part->NumberTrajectoryPoints(); itp++) {
	    TVector3 p(part->Vx(itp), part->Vy(itp), part->Vz(itp));
	    double dr = (p-point).Mag();
	    if ( dr<drmin ) {
	      mom = part->P(itp);
	      drmin = dr;
	    }
	  }

	  // direction is always from the previous point along track to
	  // the next point, that is why sign is there
	  moms.push_back(TVector3(mom*sign*(xyz1[0] - xyz[0])/dx,
				  mom*sign*(xyz1[1] - xyz[1])/dx,
				  mom*sign*(xyz1[2] - xyz[2])/dx));

	  /*************************************************************/
	  /*                          WARNING                          */
	  /*************************************************************/
	  /* The dQdx information in recob::Track has been deprecated  */
	  /* since 2016 and in 11/2018 the recob::Track interface was  */
	  /* changed and DQdxAtPoint and NumberdQdx were removed.      */
	  /* Therefore the code below is now commented out             */
	  /* (note that it was most likely not functional anyways).    */
	  /* For any issue please contact: larsoft-team@fnal.gov       */
	  /*************************************************************/
	  /*
	    dQdx[0].push_back(charge/dx);
	    dQdx[1].push_back(charge/dx);
	    dQdx[2].push_back(charge/dx);
	  */
	  /*************************************************************/

	  // make the space point and set its ID and XYZ
	  double xyzerr[6] = {1.e-3};
	  double chisqr    = 0.9;
	  recob::SpacePoint sp(&xyz[0], xyzerr, chisqr, clusterMapItr.first*10000 + t);
	  spcol->push_back(sp);
	}

	size_t spEnd = spcol->size();

	// add a track to the collection.  Make the track
	// ID be the same as the track ID for the eve particle
	trackcol->push_back(recob::Track(recob::TrackTrajectory(recob::tracking::convertCollToPoint(points),
								recob::tracking::convertCollToVector(moms),
								recob::Track::Flags_t(points.size()), true),
					 0, -1., 0, recob::tracking::SMatrixSym55(), recob::tracking::SMatrixSym55(), clusterMapItr.first));

	// associate the track with its clusters
	util::CreateAssn(*this, evt, *trackcol, ptrvs, *tcassn);

	// assume the input tracks were previously associated with hits
	hits.clear();
	for(size_t p = 0; p < ptrvs.size(); ++p){
	  hits  = fmh.at(idxs[p]);
	  util::CreateAssn(*this, evt, *trackcol, hits, *thassn);
	}

	// associate the track to the space points
	util::CreateAssn(*this, evt, *trackcol, *spcol, *tspassn, spStart, spEnd);

	mf::LogInfo("TrackCheater") << "adding track: \n"
				     << trackcol->back()
				     << "\nto collection.";

      }//end if this is a track

    } // end loop over the map

    evt.put(std::move(trackcol));
    evt.put(std::move(spcol));
    evt.put(std::move(tcassn));
    evt.put(std::move(thassn));
    evt.put(std::move(tspassn));

    return;

  } // end produce

} // end namespace

namespace trkf{

  DEFINE_ART_MODULE(TrackCheater)

}

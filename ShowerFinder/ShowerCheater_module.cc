////////////////////////////////////////////////////////////////////////
// $Id: ShowerCheater_module.cc Exp $
//
// ShowerCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef SHWF_SHOWERCHEATER_H
#define SHWF_SHOWERCHEATER_H
#include <string>
#include <vector>

// ROOT includes

// LArSoft includes
#include "Geometry/Geometry.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/AssociationUtil.h"
#include "SimpleTypesAndConstants/PhysicalConstants.h"


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace shwf {
  class ShowerCheater : public art::EDProducer {
  public:
    explicit ShowerCheater(fhicl::ParameterSet const& pset);
    virtual ~ShowerCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fCheatedClusterLabel; ///< label for module creating recob::Cluster objects	   
    std::string fG4ModuleLabel;       ///< label for module running G4 and making particles, etc

  };
}

namespace shwf{

  //--------------------------------------------------------------------
  ShowerCheater::ShowerCheater(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Shower> >();
    produces< std::vector<recob::SpacePoint> >();
    produces< art::Assns<recob::Shower, recob::Cluster> >();
    produces< art::Assns<recob::Shower, recob::SpacePoint> >();
    produces< art::Assns<recob::Shower, recob::Hit> >();
    produces< art::Assns<recob::Hit, recob::SpacePoint> >();
  }

  //--------------------------------------------------------------------
  ShowerCheater::~ShowerCheater()
  {
  }

  //--------------------------------------------------------------------
  void ShowerCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCheatedClusterLabel = pset.get< std::string >("CheatedClusterLabel", "cluster" );
    fG4ModuleLabel       = pset.get< std::string >("G4ModuleLabel",       "largeant");

    return;
  }

  //--------------------------------------------------------------------
  void ShowerCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<cheat::BackTracker> bt;
    art::ServiceHandle<geo::Geometry> geo;

    // grab the clusters that have been reconstructed
    art::Handle< std::vector<recob::Cluster> > clustercol;
    evt.getByLabel(fCheatedClusterLabel, clustercol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Cluster> > clusters;
    art::fill_ptr_vector(clusters, clustercol);
    
    // loop over the clusters and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Cluster> >::iterator itr = clusters.begin();

    // make a map of vectors of art::Ptrs keyed by eveID values
    std::map< int, std::vector< art::Ptr<recob::Cluster> > > eveClusterMap;
    std::map< int, std::vector< art::Ptr<recob::Cluster> > >::iterator clusterMapItr = eveClusterMap.begin();

    // loop over all clusters and fill in the map
    while( itr != clusters.end() ){

      // in the ClusterCheater module we set the cluster ID to be 
      // the eve particle track ID*1000 + plane*100 + tpc*10 + cryostat number.  The
      // floor function on the cluster ID / 1000 will give us
      // the eve track ID
      int eveID = floor((*itr)->ID()/1000.);

      clusterMapItr = eveClusterMap.find(eveID);
	
      // is this id already in the map, if so extend the collection 
      // by one hit, otherwise make a new collection and put it in
      // the map
      if( clusterMapItr != eveClusterMap.end() ){
	  ((*clusterMapItr).second).push_back((*itr));
      }
      else{
	std::vector< art::Ptr<recob::Cluster> > clustervec;
	clustervec.push_back(*itr);
	eveClusterMap[eveID] = clustervec;
      }

      itr++;
    }// end loop over clusters

    // loop over the map and make prongs
    std::unique_ptr< std::vector<recob::Shower> > showercol(new std::vector<recob::Shower>);
    std::unique_ptr< std::vector<recob::SpacePoint> > spcol(new std::vector<recob::SpacePoint>);
    std::unique_ptr< art::Assns<recob::Shower, recob::Cluster> > scassn(new art::Assns<recob::Shower, recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Shower, recob::Hit> > shassn(new art::Assns<recob::Shower, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Shower, recob::SpacePoint> > sspassn(new art::Assns<recob::Shower, recob::SpacePoint>);
    std::unique_ptr< art::Assns<recob::Hit, recob::SpacePoint> > sphassn(new art::Assns<recob::Hit, recob::SpacePoint>);

    for(clusterMapItr = eveClusterMap.begin(); clusterMapItr != eveClusterMap.end(); clusterMapItr++){

      // separate out the hits for each particle into the different views
      std::vector< art::Ptr<recob::Cluster> > eveClusters( (*clusterMapItr).second );

      art::PtrVector<recob::Cluster> ptrvs;

      size_t startSPIndx = spcol->size();

      double totalCharge = 0.;

      for(size_t c = 0; c < eveClusters.size(); ++c){
	ptrvs.push_back(eveClusters[c]);
	size_t cindx = ptrvs.size() - 1;
	// need to make the space points for this prong
	// loop over the hits for this cluster and make 
	// a space point for each one
	// set the SpacePoint ID to be the cluster ID*10000 
	// + the hit index in the cluster PtrVector of hits
	art::FindManyP<recob::Hit> fmh(ptrvs, evt, fCheatedClusterLabel);
	std::vector< art::Ptr<recob::Hit> > hits = fmh.at(cindx);
	
	for(size_t h = 0; h < hits.size(); ++h){
	  art::Ptr<recob::Hit> hit = hits[h];
	  // add up the charge from the hits on the collection plane
	  if(hit->SignalType() == geo::kCollection) totalCharge += hit->Charge();
	  std::vector<double> xyz = bt->HitToXYZ(hit);
	  double sperr[6] = {0.01, 0.01, 0.1, 0.001, 0.001, 0.001};

	  // make the space point and set its ID and XYZ
	  // the errors and chi^2 are set to "good" values as we know the information perfectly
	  recob::SpacePoint sp(&xyz[0],
			       sperr,
			       0.9,
			       eveClusters[c]->ID()*10000 + h);
	  spcol->push_back(sp);

	  // associate the space point to the hit
	  util::CreateAssn(*this, evt, *spcol, hit, *sphassn);

	}
      }
      
      size_t endSPIndx = spcol->size();

      // is this prong electro-magnetic in nature or 
      // hadronic/muonic?  EM --> shower, everything else is a track
      if( abs(bt->ParticleList()[(*clusterMapItr).first]->PdgCode()) == 11  ||
	  abs(bt->ParticleList()[(*clusterMapItr).first]->PdgCode()) == 22  ||
	  abs(bt->ParticleList()[(*clusterMapItr).first]->PdgCode()) == 111 ){

	mf::LogInfo("ShowerCheater") << "prong of " << (*clusterMapItr).first 
				    << " is a shower with pdg code "
				    << bt->ParticleList()[(*clusterMapItr).first]->PdgCode();

	// get the direction cosine for the eve ID particle
	// just use the same for both the start and end of the prong
	const TLorentzVector initmom = bt->ParticleList()[(*clusterMapItr).first]->Momentum();
	double dcos[3] = { initmom.Px()/initmom.Mag(),
			   initmom.Py()/initmom.Mag(),
			   initmom.Pz()/initmom.Mag() };
	double dcosErr[3] = { 1.e-3, 1.e-3, 1.e-3 };
	
	/// \todo figure out the max transverse width of the shower in the x and y directions
	double maxTransWidth[2] = { util::kBogusD, util::kBogusD };
	double distanceMaxWidth = util::kBogusD;

	// add a prong to the collection.  Make the prong
	// ID be the same as the track ID for the eve particle
	showercol->push_back(recob::Shower(dcos, dcosErr, maxTransWidth, 
					   distanceMaxWidth, totalCharge, (*clusterMapItr).first));

	// associate the shower with its clusters
	util::CreateAssn(*this, evt, *showercol, ptrvs, *scassn);

	art::FindManyP<recob::Hit> fmh(ptrvs, evt, fCheatedClusterLabel);

	// get the hits associated with each cluster and associate those with the shower
	for(size_t p = 0; p < ptrvs.size(); ++p){
	  std::vector< art::Ptr<recob::Hit> > hits = fmh.at(p);
	  util::CreateAssn(*this, evt, *showercol, hits, *shassn);
	}

	// associate the shower with the space points
	util::CreateAssn(*this, evt, *showercol, *spcol, *sspassn, startSPIndx, endSPIndx);

	mf::LogInfo("ShowerCheater") << "adding shower: \n" 
				     << showercol->back()
				     << "\nto collection.";

      }// end if this is a shower
    } // end loop over the map

    evt.put(std::move(showercol));
    evt.put(std::move(spcol));
    evt.put(std::move(scassn));
    evt.put(std::move(shassn));
    evt.put(std::move(sspassn));
    evt.put(std::move(sphassn));

    return;

  } // end produce

} // end namespace

namespace shwf{

  DEFINE_ART_MODULE(ShowerCheater)

}

#endif

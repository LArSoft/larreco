////////////////////////////////////////////////////////////////////////
// $Id: VertexCheater_module.cc Exp $
//
// VertexCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef VERTEX_VERTEXCHEATER_H
#define VERTEX_VERTEXCHEATER_H
#include <string>
#include <vector>

// ROOT includes

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "Simulation/ParticleList.h"

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

namespace vertex {
  class VertexCheater : public art::EDProducer {
  public:
    explicit VertexCheater(fhicl::ParameterSet const& pset);
    virtual ~VertexCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fCheatedTrackLabel;  ///< label for module creating recob::Track objects
    std::string fCheatedShowerLabel; ///< label for module creating recob::Shower objects
    std::string fG4ModuleLabel;      ///< label for module running G4 and making particles, etc

  };
}

namespace vertex{

  //--------------------------------------------------------------------
  VertexCheater::VertexCheater(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
    produces< art::Assns<recob::Vertex, recob::Track>  >();
    produces< art::Assns<recob::Vertex, recob::Hit>    >();

  }

  //--------------------------------------------------------------------
  VertexCheater::~VertexCheater()
  {
  }

  //--------------------------------------------------------------------
  void VertexCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCheatedTrackLabel  = pset.get< std::string >("CheatedTrackLabel",  "track"   );
    fCheatedShowerLabel = pset.get< std::string >("CheatedShowerLabel", "shower"  );
    fG4ModuleLabel      = pset.get< std::string >("G4ModuleLabel",      "largeant");

    return;
  }

  //--------------------------------------------------------------------
  void VertexCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<cheat::BackTracker> bt;

    // grab the sim::ParticleList
    const sim::ParticleList& plist = bt->ParticleList();

    // grab the showers that have been reconstructed
    art::Handle< std::vector<recob::Shower> > showercol;
    evt.getByLabel(fCheatedShowerLabel, showercol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Shower> > showers;
    try{
      art::fill_ptr_vector(showers, showercol);
    }
    catch(cet::exception &e){
      mf::LogWarning("VertexCheater") << "showers: " << e;
    }

    // grab the tracks that have been reconstructed
    art::Handle< std::vector<recob::Track> > trackcol;
    evt.getByLabel(fCheatedTrackLabel, trackcol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Track> > tracks;
    try{
      art::fill_ptr_vector(tracks, trackcol);
    }
    catch(cet::exception &e){
      mf::LogWarning("VertexCheater") << "tracks: " << e;
    }

    // loop over the prongs and figure out which primaries they are associated with
    std::vector< art::Ptr<recob::Shower> >::iterator shwitr = showers.begin();
    std::vector< art::Ptr<recob::Track> >::iterator trkitr  = tracks.begin();

    // protect against events where there are either no tracks or no showers
    if(tracks.size()  < 1) trkitr = tracks.end();
    if(showers.size() < 1) shwitr = showers.end();

    // make a map of eve id's to collections of prongs
    std::map<int, std::vector< art::Ptr<recob::Shower> > > eveShowerMap;
    std::map<int, std::vector< art::Ptr<recob::Shower> > >::iterator showerMapItr = eveShowerMap.begin();

    std::map<int, std::vector< art::Ptr<recob::Track> > > eveTrackMap;
    std::map<int, std::vector< art::Ptr<recob::Track> > >::iterator trackMapItr = eveTrackMap.begin();
    
    // make a map of eve id's to the number of prongs corresponding to that id
    std::vector<int> eveIDs;

    // loop over all showers
    while( shwitr != showers.end() ){

      // in the ProngCheater module we set the prong ID to be 
      // the particle track ID creating the energy depositions
      // of the prong
      int prongID = (*shwitr)->ID();

      // now get the mother particle of this prong if it exists
      // set the eveID to the mother particle track ID, or to the
      // ID of this prong if it is primary and mother ID = 0
      int eveID = plist[prongID]->Mother();
      if( eveID < 1 ) eveID = prongID;

      if(std::find(eveIDs.begin(), eveIDs.end(), eveID) == eveIDs.end())
	eveIDs.push_back(eveID);

      // now we want to associate all prongs having the same 
      // eve ID, so look to see if there are other prongs with that
      // ID
      showerMapItr = eveShowerMap.find(eveID);

      mf::LogInfo("VertexCheater") << "shower: " << prongID << " has mother " << eveID;
	
      // is this id already in the map, if so extend the collection 
      // by one hit, otherwise make a new collection and put it in
      // the map
      if( showerMapItr != eveShowerMap.end() ){
	  ((*showerMapItr).second).push_back((*shwitr));
      }
      else{
	std::vector< art::Ptr<recob::Shower> > showervec;
	showervec.push_back(*shwitr);
	eveShowerMap[eveID] = showervec;
      }

      shwitr++;
    }// end loop over showers

    // loop over all tracks
    while( trkitr != tracks.end() ){

      // in the ProngCheater module we set the prong ID to be 
      // the particle track ID creating the energy depositions
      // of the prong
      int prongID = (*trkitr)->ID();

      // now get the mother particle of this prong if it exists
      // set the eveID to the mother particle track ID, or to the
      // ID of this prong if it is primary and mother ID = 0
      int eveID = plist[prongID]->Mother();
      if( eveID < 1 ) eveID = prongID;

      if(std::find(eveIDs.begin(), eveIDs.end(), eveID) == eveIDs.end())
	eveIDs.push_back(eveID);

      // now we want to associate all prongs having the same 
      // eve ID, so look to see if there are other prongs with that
      // ID
      trackMapItr = eveTrackMap.find(eveID);

      mf::LogInfo("VertexCheater") << "track: " << prongID << " has mother " << eveID;
	
      // is this id already in the map, if so extend the collection 
      // by one hit, otherwise make a new collection and put it in
      // the map
      if( trackMapItr != eveTrackMap.end() ){
	  ((*trackMapItr).second).push_back((*trkitr));
      }
      else{
	std::vector< art::Ptr<recob::Track> > trackvec;
	trackvec.push_back(*trkitr);
	eveTrackMap[eveID] = trackvec;
      }

      trkitr++;
    }// end loop over tracks

    std::unique_ptr< std::vector<recob::Vertex> > vertexcol(new std::vector<recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> > vsassn(new art::Assns<recob::Vertex, recob::Shower>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track>  > vtassn(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit>    > vhassn(new art::Assns<recob::Vertex, recob::Hit>);

    // loop over the eve ID values and make Vertexs
    for(std::vector<int>::iterator eItr = eveIDs.begin(); eItr != eveIDs.end(); eItr++){

      int eveID = *eItr;
      
      // Vertex objects require PtrVectors of showers and tracks as well
      // as a vertex position for their constructor
      art::PtrVector<recob::Shower> ptrvshw;
      art::PtrVector<recob::Track>  ptrvtrk;

      // first get the showers
      if(eveShowerMap.find(eveID) != eveShowerMap.end()){
	std::vector< art::Ptr<recob::Shower> > eveShowers( eveShowerMap[eveID] );
	for(size_t s = 0; s < eveShowers.size(); ++s){
	  ptrvshw.push_back(eveShowers[s]);
	} // end loop over showers for this particle
      } // end find showers for this particle

      // now the tracks
      if(eveTrackMap.find(eveID) != eveTrackMap.end()){
	std::vector< art::Ptr<recob::Track> > eveTracks( eveTrackMap[eveID] );
	for(size_t t = 0; t < eveTracks.size(); ++t){
	  ptrvtrk.push_back(eveTracks[t]);
	} // end loop over tracks for this particle
      } // end find tracks for this particle
 
      double xyz[3] = { plist[eveID]->Vx(),
			plist[eveID]->Vy(),
			plist[eveID]->Vz() };

      // add a vector to the collection.  
      vertexcol->push_back(recob::Vertex(xyz, eveID));

      // associate the vertex with its showers and tracks

      if( ptrvtrk.size() > 0 ){
	util::CreateAssn(*this, evt, *vertexcol, ptrvtrk, *vtassn);

	// get the hits associated with each track and associate those with the vertex
	art::FindManyP<recob::Hit> fmh(ptrvtrk, evt, fCheatedTrackLabel);
	for(size_t p = 0; p < ptrvtrk.size(); ++p){
	  std::vector< art::Ptr<recob::Hit> > hits = fmh.at(p);
	  util::CreateAssn(*this, evt, *vertexcol, hits, *vhassn);
	}
      }
      
      if( ptrvshw.size() > 0 ){
	util::CreateAssn(*this, evt, *vertexcol, ptrvshw, *vsassn);
	// get the hits associated with each shower and associate those with the vertex
	art::FindManyP<recob::Hit> fmh(ptrvshw, evt, fCheatedShowerLabel);
	for(size_t p = 0; p < ptrvshw.size(); ++p){
	  std::vector< art::Ptr<recob::Hit> > hits = fmh.at(p);
	  util::CreateAssn(*this, evt, *vertexcol, hits, *vhassn);
	}
      }

      mf::LogInfo("VertexCheater") << "adding vertex: \n" 
				   << vertexcol->back()
				   << "\nto collection.";

    } // end loop over the eve ID values

    evt.put(std::move(vertexcol));
    evt.put(std::move(vsassn));
    evt.put(std::move(vtassn));
    evt.put(std::move(vhassn));

    return;

  } // end produce

} // end namespace


namespace vertex{

  DEFINE_ART_MODULE(VertexCheater)

}

#endif

////////////////////////////////////////////////////////////////////////
//
// VertexCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>

// ROOT includes

// LArSoft includes
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nutools/ParticleNavigation/ParticleList.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace vertex {
  class VertexCheater : public art::EDProducer {
  public:
    explicit VertexCheater(fhicl::ParameterSet const& pset);

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
    : EDProducer{pset}
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
    produces< art::Assns<recob::Vertex, recob::Track>  >();
    produces< art::Assns<recob::Vertex, recob::Hit>    >();

  }

  //--------------------------------------------------------------------
  void VertexCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCheatedTrackLabel  = pset.get< std::string >("CheatedTrackLabel",  "track"   );
    fCheatedShowerLabel = pset.get< std::string >("CheatedShowerLabel", "shower"  );
    fG4ModuleLabel      = pset.get< std::string >("G4ModuleLabel",      "largeant");
  }

  //--------------------------------------------------------------------
  void VertexCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // grab the sim::ParticleList
    const sim::ParticleList& plist = pi_serv->ParticleList();

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

    art::FindManyP<recob::Hit> fmhs(showers, evt, fCheatedShowerLabel);
    art::FindManyP<recob::Hit> fmht(tracks, evt, fCheatedTrackLabel);

    // loop over the prongs and figure out which primaries they are associated with
    std::vector< art::Ptr<recob::Shower> >::iterator shwitr = showers.begin();
    std::vector< art::Ptr<recob::Track>  >::iterator trkitr = tracks.begin();

    // protect against events where there are either no tracks or no showers
    if(tracks.size()  < 1) trkitr = tracks.end();
    if(showers.size() < 1) shwitr = showers.end();

    // make a map of eve id's to collections of prongs
    std::map<int, std::vector< std::pair<size_t, art::Ptr<recob::Shower> > > > eveShowerMap;
    std::map<int, std::vector< std::pair<size_t, art::Ptr<recob::Track>  > > > eveTrackMap;
    
    // make a map of eve id's to the number of prongs corresponding to that id
    std::vector<int> eveIDs;

    // loop over all showers
    for(size_t s = 0; s < showers.size(); ++s){

      std::pair<size_t, art::Ptr<recob::Shower> > idxShw(s, showers[s]);

      // in the ProngCheater module we set the prong ID to be 
      // the particle track ID creating the energy depositions
      // of the prong
      int prongID = showers[s]->ID();

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
      eveShowerMap[eveID].push_back(idxShw);

      mf::LogInfo("VertexCheater") << "shower: " << prongID << " has mother " << eveID;
	
    }// end loop over showers

    // loop over all tracks
    for(size_t t = 0; t < tracks.size(); ++t){

      std::pair<size_t, art::Ptr<recob::Track> > idxTrk(t, tracks[t]);

      // in the ProngCheater module we set the prong ID to be 
      // the particle track ID creating the energy depositions
      // of the prong
      int prongID = tracks[t]->ID();

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
      eveTrackMap[eveID].push_back(idxTrk);

      mf::LogInfo("VertexCheater") << "track: " << prongID << " has mother " << eveID;
	
    }// end loop over tracks

    std::unique_ptr< std::vector<recob::Vertex> > vertexcol(new std::vector<recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> > vsassn(new art::Assns<recob::Vertex, recob::Shower>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track>  > vtassn(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit>    > vhassn(new art::Assns<recob::Vertex, recob::Hit>);

    // loop over the eve ID values and make Vertexs
    for(auto const& eveID : eveIDs){

      // Vertex objects require PtrVectors of showers and tracks as well
      // as a vertex position for their constructor
      art::PtrVector<recob::Shower> ptrvshw;
      art::PtrVector<recob::Track>  ptrvtrk;
      std::vector<size_t> idxShw;
      std::vector<size_t> idxTrk;

      // first get the showers
      if(eveShowerMap.find(eveID) != eveShowerMap.end()){
	auto const& eveShowers = eveShowerMap[eveID];
	for(auto const& is : eveShowers){
	  ptrvshw.push_back(is.second);
	  idxShw.push_back(is.first);
	} // end loop over showers for this particle
      } // end find showers for this particle

      // now the tracks
      if(eveTrackMap.find(eveID) != eveTrackMap.end()){
	auto const& eveTracks = eveTrackMap[eveID];
	for(auto const& it : eveTracks){
	  ptrvtrk.push_back(it.second);
	  idxTrk.push_back(it.first);
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
	for(auto const& i : idxTrk){
	  std::vector< art::Ptr<recob::Hit> > hits = fmht.at(i);
	  util::CreateAssn(*this, evt, *vertexcol, hits, *vhassn);
	}
      }
      
      if( ptrvshw.size() > 0 ){
	util::CreateAssn(*this, evt, *vertexcol, ptrvshw, *vsassn);
	// get the hits associated with each shower and associate those with the vertex
	for(auto const& i : idxShw){
	  std::vector< art::Ptr<recob::Hit> > hits = fmhs.at(i);
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

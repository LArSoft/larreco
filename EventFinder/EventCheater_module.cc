////////////////////////////////////////////////////////////////////////
// $Id: EventCheater_module.cc Exp $
//
// EventCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef EVENT_EVENTCHEATER_H
#define EVENT_EVENTCHEATER_H
#include <string>
#include <vector>

// ROOT includes

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Event.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Hit.h"
#include "SimulationBase/MCTruth.h"
#include "Simulation/ParticleList.h"
#include "Utilities/AssociationUtil.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/FindOne.h"
#include "art/Framework/Core/EDProducer.h"

///Event finding and building
namespace event {
  class EventCheater : public art::EDProducer {
  public:
    explicit EventCheater(fhicl::ParameterSet const& pset);
    virtual ~EventCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

  private:

    std::string fCheatedVertexLabel; ///< label for module creating recob::Vertex objects	   
    std::string fG4ModuleLabel;      ///< label for module running G4 and making particles, etc

  };
}

namespace event{

  //--------------------------------------------------------------------
  EventCheater::EventCheater(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Event> >();
    produces< art::Assns<recob::Event, recob::Vertex> >();
    produces< art::Assns<recob::Event, recob::Hit> >();
  }

  //--------------------------------------------------------------------
  EventCheater::~EventCheater()
  {
  }

  //--------------------------------------------------------------------
  void EventCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCheatedVertexLabel = pset.get< std::string >("CheatedVertexLabel", "prong" );
    fG4ModuleLabel      = pset.get< std::string >("G4ModuleLabel",      "largeant");

    return;
  }

  //--------------------------------------------------------------------
  void EventCheater::produce(art::Event& evt)
  {

    art::View<simb::MCParticle> pcol;
    evt.getView(fG4ModuleLabel, pcol);

    art::FindOneP<simb::MCTruth> fo(pcol, evt, fG4ModuleLabel);

    // make a map of the track id for each sim::Particle to its entry in the 
    // collection of sim::Particles
    std::map<int, int> trackIDToPColEntry;
    for(size_t p = 0; p < pcol.vals().size(); ++p) trackIDToPColEntry[pcol.vals().at(p)->TrackId()] = p;

    // grab the vertices that have been reconstructed
    art::Handle< std::vector<recob::Vertex> > vertexcol;
    evt.getByLabel(fCheatedVertexLabel, vertexcol);

    art::FindManyP<recob::Hit> fm(vertexcol, evt, fCheatedVertexLabel);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Vertex> > vertices;
    art::fill_ptr_vector(vertices, vertexcol);

    // loop over the vertices and figure out which primaries they are associated with
    std::vector< art::Ptr<recob::Vertex> >::iterator vertexitr = vertices.begin();

    // make a map of primary product id's to collections of vertices
    std::map<art::Ptr<simb::MCTruth>, std::vector< art::Ptr<recob::Vertex> > > vertexMap;
    std::map<art::Ptr<simb::MCTruth>, std::vector< art::Ptr<recob::Vertex> > >::iterator vertexMapItr = vertexMap.begin();

    // loop over all prongs
    while( vertexitr != vertices.end() ){

      size_t pcolEntry = trackIDToPColEntry.find((*vertexitr)->ID())->second;
      const art::Ptr<simb::MCTruth> primary = fo.at(pcolEntry);
      
      vertexMap[primary].push_back(*vertexitr);

      vertexitr++;
    }// end loop over vertices

    std::unique_ptr< std::vector<recob::Event> >               eventcol(new std::vector<recob::Event>);
    std::unique_ptr< art::Assns<recob::Event, recob::Vertex> > evassn(new art::Assns<recob::Event, recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Event, recob::Hit> >    ehassn(new art::Assns<recob::Event, recob::Hit>);

    // loop over the map and associate all vertex objects with an event
    for(vertexMapItr = vertexMap.begin(); vertexMapItr != vertexMap.end(); vertexMapItr++){

      art::PtrVector<recob::Vertex> ptrvs;

      std::vector< art::Ptr<recob::Vertex> > verts( (*vertexMapItr).second );

      // add an event to the collection.  
      eventcol->push_back(recob::Event(eventcol->size()-1));

      // associate the event with its vertices
      util::CreateAssn(*this, evt, *eventcol, verts, *evassn);
      
      // get the hits associated with each vertex and associate those with the event
      for(size_t p = 0; p < ptrvs.size(); ++p){
	std::vector< art::Ptr<recob::Hit> > hits = fm.at(p);
	util::CreateAssn(*this, evt, *eventcol, hits, *ehassn);
      }


      mf::LogInfo("EventCheater") << "adding event: \n" 
				  << eventcol->back()
				  << "\nto collection";

    } // end loop over the map

    evt.put(std::move(eventcol));
    evt.put(std::move(evassn));
    evt.put(std::move(ehassn));

    return;

  } // end produce

} // end namespace

namespace event{

  DEFINE_ART_MODULE(EventCheater)

}

#endif

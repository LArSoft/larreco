#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/Event.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"

namespace event {
  class EventMaker;
}

class event::EventMaker : public art::EDProducer {
public:
  explicit EventMaker(fhicl::ParameterSet const &p);

private:
  void produce(art::Event &e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  std::string fVertexModuleLabel;  ///< label of the module making the recob::Vertex objects
  double      fProximity;          ///< how close a vertex needs to be to another to be from the same event

};

//--------------------------------------------------------
event::EventMaker::EventMaker(fhicl::ParameterSet const &p)
  : EDProducer{p}
{
  this->reconfigure(p);
  
  produces< std::vector<recob::Event> >();
  produces< art::Assns<recob::Event, recob::Vertex> >();
  produces< art::Assns<recob::Event, recob::Hit> >();
}

//--------------------------------------------------------
void event::EventMaker::produce(art::Event &e) 
{

  // make the unique_ptr of the vector for the recob::Events
  std::unique_ptr< std::vector<recob::Event> > eventcol(new std::vector<recob::Event>);
  std::unique_ptr< art::Assns<recob::Event, recob::Vertex> > evassn(new art::Assns<recob::Event, recob::Vertex>);
  std::unique_ptr< art::Assns<recob::Event, recob::Hit> >    ehassn(new art::Assns<recob::Event, recob::Hit>);

  // first get the recob::Vertex objects out of the event
  art::Handle< std::vector<recob::Vertex> > vtxHandle;
  e.getByLabel(fVertexModuleLabel, vtxHandle);

  art::FindManyP<recob::Hit> fmhv(vtxHandle, e, fVertexModuleLabel);

  // get a collection of art::Ptrs
  std::vector< art::Ptr<recob::Vertex> > vtxs;
  art::fill_ptr_vector(vtxs, vtxHandle);

  // if only 1 vertex in the event, life is easy
  if(vtxs.size() == 1){
    art::PtrVector<recob::Vertex> vtxvec;
    vtxvec.push_back(*(vtxs.begin()));
    recob::Event evt(0);
    eventcol->push_back(evt);

    // associate the event with its vertex
    util::CreateAssn(*this, e, *eventcol, vtxvec, *evassn);
    
    // get the hits associated with the vertex and associate those with the event
    std::vector< art::Ptr<recob::Hit> > hits = fmhv.at(0);
    util::CreateAssn(*this, e, *eventcol, hits, *ehassn);

    e.put(std::move(eventcol));
    e.put(std::move(ehassn));
    e.put(std::move(evassn));

    return;
  }

  // now the hard part, we have multiple vertex objects
  // things to consider:
  // 1. all particles coming from a common vertex location are in the same event
  // 2. particles coming from decays of another particle, ie the electron from a muon decay
  // 3. pi-zero decay vertex will be offset from the interaction vertex
  int evtctr = 0;
  std::vector< art::Ptr<recob::Vertex> >::iterator itr  = vtxs.begin();
  std::vector< art::Ptr<recob::Vertex> >::iterator itr2 = vtxs.begin();

  // make a map of the Ptr to the index in the original vector
  std::map<art::Ptr<recob::Vertex>, size_t> ptrToIdx;
  for(size_t v = 0; v < vtxs.size(); ++v) ptrToIdx[ vtxs[v] ] = v; 

  while( itr != vtxs.end() ){
    
    art::PtrVector< recob::Vertex > vtxVec;

    // get the current vertex object and put it into the vector
    art::Ptr<recob::Vertex> curvtx = *itr;
    vtxVec.push_back(curvtx);

    double curxyz[3] = {0.};
    curvtx->XYZ(curxyz);

    // make itr2 the next one in the list, also remove this one from the list
    // as it is being used
    itr2 = vtxs.erase(itr);
    while( itr2 != vtxs.end() ){
      art::Ptr<recob::Vertex> nexvtx = *itr2;
      
      // get the xyz location of the vertex to compare to the current vertex
      double nextxyz[3] = {0.};
      nexvtx->XYZ(nextxyz);
      if( std::sqrt((curxyz[0]-nextxyz[0])*(curxyz[0]-nextxyz[0]) +
		    (curxyz[1]-nextxyz[1])*(curxyz[1]-nextxyz[1]) +
		    (curxyz[2]-nextxyz[2])*(curxyz[2]-nextxyz[2]) ) <= fProximity){
	
	// add this one to the vector and remove it from the list as it is being used
	vtxVec.push_back(nexvtx);
	itr2 = vtxs.erase(itr2);
      }// end if vertices are close enough to be considered from the same event
      else itr2++;
    }// end inner loop
    
    // make an event from these vertex objects and add them to the collection
    recob::Event evt(++evtctr);
    eventcol->push_back(evt);

    // associate the event with its vertices
    util::CreateAssn(*this, e, *eventcol, vtxVec, *evassn);
    
    // get the hits associated with each vertex and associate those with the event
    for(size_t p = 0; p < vtxVec.size(); ++p){
      std::vector< art::Ptr<recob::Hit> > hits = fmhv.at( ptrToIdx[ vtxVec[p] ] );
      util::CreateAssn(*this, e, *eventcol, hits, *ehassn);
    }
    
    // move the initial iterator forward
    itr++;

  }

  // put the collection of events in the art::Event
  e.put(std::move(eventcol));
  e.put(std::move(ehassn));
  e.put(std::move(evassn));

  return;

}

//--------------------------------------------------------
void event::EventMaker::reconfigure(fhicl::ParameterSet const & p) 
{
  fVertexModuleLabel = p.get< std::string >("VertexModuleLabel");
  fProximity         = p.get< double      >("Proximity");
}



DEFINE_ART_MODULE(event::EventMaker)

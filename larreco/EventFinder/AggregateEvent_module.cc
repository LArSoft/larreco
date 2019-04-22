////////////////////////////////////////////////////////////////////////
//
// AggregateEvent module
//
// echurch@fnal.gov
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <fstream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "lardataobj/RecoBase/Event.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcore/Geometry/Geometry.h"

// ROOT includes
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TGraph.h"

#include "art/Framework/Core/EDProducer.h"

#include "TH2F.h"
#include "TF1.h"
#include <string>

namespace recob{ class Vertex; }

namespace event {

  class AggregateEvent : public art::EDProducer {

  public:

    explicit AggregateEvent(fhicl::ParameterSet const& );

  private:
    void produce(art::Event& evt) override;

    std::string fVertexModuleLabel;

    art::PtrVector<recob::Vertex> fvertexlist;

  }; // class AggregateEvent

} // namespace


namespace event {

  //-------------------------------------------------
  AggregateEvent::AggregateEvent(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    fVertexModuleLabel(pset.get< std::string >("VertexModuleLabel"))
  {

    produces< std::vector<recob::Event> >();

  }

  //------------------------------------------------------------------------------------//
  void AggregateEvent::produce(art::Event& evt)
  {

    std::unique_ptr<std::vector<recob::Event> > ecol(new std::vector<recob::Event>);

    // get the geometry
    art::ServiceHandle<geo::Geometry const> geom;

    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);
    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii){
      art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      fvertexlist.push_back(vertex); // class member
    }

  }

  DEFINE_ART_MODULE(AggregateEvent)

} // end namespace

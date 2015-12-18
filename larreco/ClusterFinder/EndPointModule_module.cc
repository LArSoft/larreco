////////////////////////////////////////////////////////////////////////
//
// EndPointModule class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find (weak) vertices from hits after deconvolution and hit finding. 
//  A weak vertex is a vertex that has been found using a dedicated vertex finding algorithm only. A 
//  strong vertex is a vertex that has been found using a dedicated vertex finding algorithm and matched 
//  to a crossing of two or more HoughLineFinder lines. The VertexMatch module finds strong vertices.
////////////////////////////////////////////////////////////////////////
/// The algorithm is based on:
///C. Harris and M. Stephens (1988). "A combined corner and edge detector". Proceedings of the 4th Alvey 
///Vision Conference. pp. 147-151.
///B. Morgan (2010). "Interest Point Detection for Reconstruction in High Granularity Tracking Detectors". 
///arXiv:1006.3012v1 [physics.ins-det]
//Thanks to B. Morgan of U. of Warwick for comments and suggestions



#include <iostream>
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"
#include "TMath.h"
#include <vector>
#include <string>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/EndPointAlg.h"

//#ifndef EndPointModule_H
//#define EndPointModule_H



/// 2D end point reconstruction
namespace cluster {

  ///module to find 2D end points 
 class EndPointModule :  public art::EDProducer {
    
  public:
    
    explicit EndPointModule(fhicl::ParameterSet const& pset); 
    virtual ~EndPointModule();        

    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt);
    
  private:

    std::string fDBScanModuleLabel;

    EndPointAlg fEPAlg;            ///< object that contains the end point finding algorithm
   
  };
    
}

//#endif // EndPointModule_H



namespace cluster {


  //-----------------------------------------------------------------------------
  EndPointModule::EndPointModule(fhicl::ParameterSet const& pset)
    : fEPAlg(pset.get< fhicl::ParameterSet >("EndPointAlg"))
  {
    this->reconfigure(pset);
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::EndPoint2D, recob::Hit> >();
  }
  
  //-----------------------------------------------------------------------------
  EndPointModule::~EndPointModule()
  {
  }
  
  //-----------------------------------------------------------------------------
  void EndPointModule::reconfigure(fhicl::ParameterSet const& p)
  {
    fDBScanModuleLabel = p.get<std::string>("DBScanModuleLabel");  
    fEPAlg.reconfigure(p.get< fhicl::ParameterSet >("EndPointAlg"));
  }
  
  //-----------------------------------------------------------------------------
  
  void EndPointModule::produce(art::Event& evt)
  {
  
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
    //Point to a collection of vertices to output.
    
    //.......................................
    art::PtrVector<recob::Cluster> clusIn;
    for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
      {
        art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
        clusIn.push_back(cluster);
      }
   
    // make a std::vector<recob::Cluster> for the output of the 
    // Hough Transform
    std::vector<recob::EndPoint2D> vtxOut;
    std::vector< art::PtrVector<recob::Hit> > vtxHitsOut;
    size_t numvtx = fEPAlg.EndPoint(clusIn, vtxOut, vtxHitsOut, evt, fDBScanModuleLabel);
  
    LOG_DEBUG("Vertex") << "found " << numvtx << "vertices with VertexService";
  
    //Point to a collection of vertices to output.
    std::unique_ptr<std::vector<recob::EndPoint2D> > vtxcol(new std::vector<recob::EndPoint2D>(vtxOut));
    std::unique_ptr< art::Assns<recob::EndPoint2D, recob::Hit> > assn(new art::Assns<recob::EndPoint2D, recob::Hit>);
  
    for(size_t v = 0; v < vtxcol->size(); ++v)
      util::CreateAssn(*this, evt, *(vtxcol.get()), vtxHitsOut[v], *(assn.get()), v);
    
    evt.put(std::move(vtxcol));   
    evt.put(std::move(assn));
  }
  

} // end namespace







namespace cluster{

  DEFINE_ART_MODULE(EndPointModule)
  
} 


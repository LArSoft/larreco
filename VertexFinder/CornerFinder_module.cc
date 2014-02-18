#ifndef CORNERFINDER_H
#define CORNDERFINDER_H
/*!
 * Title:   CornerFinder class
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Wire (calibrated)
 * Outputs: recob::EndPoint2D
 *
 * Description:
 * This module, is designed to produce EndPoint2D's using the CornerFinderAlg.
 * It's inlikely this module would actually be used on its own for vertex finding, but
 * is meant to be used as an input into a different module, or as a test piece.
 */


//Basic includes
#include <vector>
#include <string>

//Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

//LArSoft includes
#include "RecoBase/EndPoint2D.h"
#include "RecoAlg/CornerFinderAlg.h"

//header info for CornerFinder class

namespace vertex {
   
  class CornerFinder :  public art::EDProducer {
    
  public:
    
    explicit CornerFinder(fhicl::ParameterSet const& pset); 
    virtual ~CornerFinder();        
    void reconfigure(fhicl::ParameterSet const& p);

    
    void produce(art::Event& evt);
    
  private:
    cluster::CornerFinderAlg  fCornerAlg;

  }; //class CornerFinder
  
  

  //-----------------------------------------------------------------------------
  CornerFinder::CornerFinder(fhicl::ParameterSet const& pset):
    fCornerAlg(pset.get<fhicl::ParameterSet>("CornerAlgParamSet"))
  {  
    this->reconfigure(pset);    
    produces< std::vector<recob::EndPoint2D> >();
  }
  
  //-----------------------------------------------------------------------------
  CornerFinder::~CornerFinder(){}

  //---------------------------------------------------------------------------
  void CornerFinder::reconfigure(fhicl::ParameterSet const& pset) {
    fCornerAlg.reconfigure(pset.get<fhicl::ParameterSet>("CornerAlgParamSet"));
    return;
  }

  //-----------------------------------------------------------------------------
  void CornerFinder::produce(art::Event& evt){
  
    //We need do very little here, as it's all handled by the corner finder.

    //First, have it process the "raw" data.
    fCornerAlg.TakeInRaw(evt);

    //now, make a vector of recob::EndPoint2Ds, and hand that to CornerAlg to fill out
    std::unique_ptr< std::vector<recob::EndPoint2D> > corner_vector(new std::vector<recob::EndPoint2D>);
    fCornerAlg.get_feature_points(*corner_vector);

    mf::LogVerbatim("CornerFinderModule") << "CornerFinderAlg finished, and returned " 
					  << corner_vector->size() << " endpoints.";

    for(std::vector<recob::EndPoint2D>::iterator iter=corner_vector->begin(); iter!=corner_vector->end(); iter++){
      geo::WireID wid = iter->WireID();
      mf::LogVerbatim("CornerFinderModule") << "Endpoint found: (plane,wire,time,strength)=(" 
					    << wid.Plane << "," 
					    << wid.Wire << "," 
					    << iter->DriftTime() << ","
					    << iter->Strength() << ")";
    }

    //and now, put this on the event.
    evt.put(std::move(corner_vector));

    //Done!

  } // end of produce


  DEFINE_ART_MODULE(CornerFinder)

}
#endif // CORNERFINDER_H

#ifndef TTSPACEPOINTFINDER_H
#define TTSPACEPOINTFINDER_H
/*!
 * Title:   TTSpacePointFinder class
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Hit
 * Outputs: recob::SpacePoint
 *
 * Description:
 * This module, TimeTickSpacePointFinder (or TTSpacePointFinder for short) is 
 * designed to produce a spacepoint object based on hits from TTHitFinder.
 * There is intention to allow for a significant number of ghost spacepoints, 
 * with some downstream application dealing with the results.
 */
#include <string>
#include <math.h>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Core/EDProducer.h" 

// LArSoft Includes
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"

namespace sppt{

  class TTSpacePointFinder : public art::EDProducer {
    
  public:
    
    explicit TTSpacePointFinder(fhicl::ParameterSet const& pset); 
    virtual ~TTSpacePointFinder();
         
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& p);                

  private:
        
    std::string    fHitModuleLabel; /// Input hit module name
    float          fTimeDiffMax;    /// Maximum allowed time difference
   
  protected:     
  
  }; // class TTSpacePointFinder  
  
  //-------------------------------------------------
  TTSpacePointFinder::TTSpacePointFinder(fhicl::ParameterSet const& pset) {
    this->reconfigure(pset);
    produces< std::vector<recob::SpacePoint> >();
    produces<art::Assns<recob::SpacePoint, recob::Hit>       >();
  }

  //-------------------------------------------------
  TTSpacePointFinder::~TTSpacePointFinder(){}

  //-------------------------------------------------
  void TTSpacePointFinder::reconfigure(fhicl::ParameterSet const& p) {
    fHitModuleLabel = p.get< std::string >("HitModuleLabel");
    fTimeDiffMax    = p.get< float       >("TimeDiffMax");

    //enforce a minimum time diff
    if(fTimeDiffMax<0){
      mf::LogError("TTSpacePointFinder") << "Time difference must be greater than zero.";
      return 0;
    }

  }

  //-------------------------------------------------
  void TTSpacePointFinder::beginJob(){}

  //-------------------------------------------------
  void TTSpacePointFinder::endJob(){}

  //-------------------------------------------------
  void TTSpacePointFinder::produce(art::Event& evt)
  { 
    

  } // End of produce()  
    
  DEFINE_ART_MODULE(TTSpacePointFinder)

} // end of hit namespace


#endif // TTSPACEPOINTFINDER_H

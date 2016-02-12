//////////////////////////////////////////////////
//
// BezierCalorimetry produces a calo object based on
// a bezier trajectory and associated hit information
//
// bjpjones@mit.edu, June 2013
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <vector>
#include <string>
#include <algorithm>

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/AnalysisAlg/CalorimetryAlg.h"
#include "lardata/RecoObjects/BezierTrack.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

///calorimetry
namespace calo {
   
    class BezierCalorimetry : public art::EDProducer {
    
    public:
    
      explicit BezierCalorimetry(fhicl::ParameterSet const& pset); 
      virtual ~BezierCalorimetry();
    
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt);

    private:
        
      std::string    fBTrackModuleLabel; ///< module creating the track objects and assns to hits
      art::ServiceHandle<geo::Geometry> fGeo;

      CalorimetryAlg caloAlg;


    }; // class BezierCalorimetry
  
}

//-------------------------------------------------
calo::BezierCalorimetry::BezierCalorimetry(fhicl::ParameterSet const& pset):
  caloAlg(pset.get< fhicl::ParameterSet >("CaloAlg"))
{
  
  this->reconfigure(pset);

  produces< std::vector<anab::Calorimetry>              >();
  produces< art::Assns<recob::Track, anab::Calorimetry> >();
}

//-------------------------------------------------
calo::BezierCalorimetry::~BezierCalorimetry()
{
  
}

//------------------------------------------------------------------------------------//
void calo::BezierCalorimetry::reconfigure(fhicl::ParameterSet const& pset)
{
  fBTrackModuleLabel = pset.get< std::string >("BTrackModuleLabel");
  return;
}

//------------------------------------------------------------------------------------//
void calo::BezierCalorimetry::produce(art::Event& evt)
{ 
  art::Handle< std::vector<recob::Track> > trackHandle;
  evt.getByLabel(fBTrackModuleLabel, "bezierformat", trackHandle);
  std::vector< art::Ptr<recob::Track> > tracks;
  art::fill_ptr_vector(tracks, trackHandle);
  
  //create anab::Calorimetry objects and make association with recob::Track
  std::unique_ptr< std::vector<anab::Calorimetry> > calorimetrycol(new std::vector<anab::Calorimetry>);
  std::unique_ptr< art::Assns<recob::Track, anab::Calorimetry> > assn(new art::Assns<recob::Track, anab::Calorimetry>);

  art::FindManyP<recob::Hit> fmht(trackHandle, evt, fBTrackModuleLabel);
  
  
  // loop over the tracks
  for(size_t t = 0; t < tracks.size(); ++t)
    {
     
      std::vector<art::Ptr<recob::Hit> > hits = fmht.at(t);
      
      art::Ptr<recob::Track> trk = tracks.at(t);
      trkf::BezierTrack BTrack(*trk);
      
      calorimetrycol->push_back(BTrack.GetCalorimetryObject(hits, geo::kCollection, caloAlg)); 
      util::CreateAssn(*this, evt, *calorimetrycol, trk, *assn);
    }

 
  if (calorimetrycol->size()>0){
    evt.put(std::move(calorimetrycol));
    evt.put(std::move(assn));
  }
  
  return;
}

namespace calo{
  
  DEFINE_ART_MODULE(BezierCalorimetry)
  
} // end namespace 


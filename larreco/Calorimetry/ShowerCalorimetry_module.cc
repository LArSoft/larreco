////////////////////////////////////////////////////////////////////////
// Class:       ShowerCalorimetry
// Plugin Type: producer (art v3_02_06)
// File:        ShowerCalorimetry_module.cc
//
// Generated at Fri Jul 12 14:14:46 2019 by Jacob Calcutt using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>

namespace calo{
  class ShowerCalorimetry;
} 

class calo::ShowerCalorimetry : public art::EDProducer {
public:
  explicit ShowerCalorimetry(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerCalorimetry(ShowerCalorimetry const&) = delete;
  ShowerCalorimetry(ShowerCalorimetry&&) = delete;
  ShowerCalorimetry& operator=(ShowerCalorimetry const&) = delete;
  ShowerCalorimetry& operator=(ShowerCalorimetry&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  int GetShowerIndex( const recob::Shower & shower, art::Event const & evt ) const;

private:

  std::string fShowerTag;
};


calo::ShowerCalorimetry::ShowerCalorimetry(fhicl::ParameterSet const& p):
  EDProducer{p},
  fShowerTag( p.get< std::string >( "ShowerTag" ) )
{
  produces< std::vector< anab::Calorimetry > >();
  //produces< art::Assns< recob::Shower, anab::Calorimetry > >();
}

void calo::ShowerCalorimetry::produce(art::Event& e) {

  //Make the container for the calo product to put onto the event.
  std::unique_ptr< std::vector<anab::Calorimetry> > caloPtr(new std::vector<anab::Calorimetry>);
  //std::vector< anab::Calorimetry > & caloVector(*caloPtr);

/*Do this later
  //Make a container for the track<-->calo associations.
  //One entry per track, with entry equal to index in calorimetry collection of associated object.
  std::vector< size_t > assnShowerCaloVector;
  std::unique_ptr< art::Assns< recob::Shower,anab::Calorimetry> > associationPtr( new  art::Assns< recob::Shower, anab::Calorimetry > );

*/
  //Make the associations for ART 
  /*for( size_t i = 0; i < assnTrackCaloVector.size(); i++ ){
    if( assnTrackCaloVector[i] == std::numeric_limits< size_t >::max() ) continue;

    art::Ptr<recob::Track> trk_ptr(trackHandle,assnTrackCaloVector[i]);
    util::CreateAssn(*this, e, caloVector, trk_ptr, *assnTrackCaloPtr, i);
  }*/



  //Get the shower handle
  auto showerHandle = e.getValidHandle< std::vector< recob::Shower > >(fShowerTag);

  //Turn it into a vector of art pointers
  std::vector< art::Ptr< recob::Shower > > recoShowers;
  art::fill_ptr_vector( recoShowers, showerHandle );

  //Also get the hits from all the showers
  art::FindManyP<recob::Hit> findHitsFromShowers(showerHandle,e,fShowerTag);

  for( size_t i = 0; i < recoShowers.size(); ++i ){
    const recob::Shower & shower = *(recoShowers.at(i));
    std::cout << GetShowerIndex(shower, e) << std::endl;
  }


  //Finish up: Put the objects into the event
  e.put( std::move( caloPtr ) );
  //e.put( std::move( associationPtr ) );
  

}

int calo::ShowerCalorimetry::GetShowerIndex( const recob::Shower & shower, art::Event const & evt ) const{
  if(shower.ID() != -999) return shower.ID();

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(fShowerTag);

  // Iterate through all showers to find the matching one to our shower
  int actualIndex = shower.ID();
  if(shower.ID() < 0){
    for(unsigned int s = 0; s < recoShowers->size(); ++s){
      const recob::Shower thisShower = (*recoShowers)[s];
      // Can't compare actual objects so look at a property
      if(fabs(thisShower.Length() - shower.Length()) < 1.e-5){
        actualIndex = s;
        continue;
      }
    }
  }

  return actualIndex;
}

DEFINE_ART_MODULE(calo::ShowerCalorimetry)

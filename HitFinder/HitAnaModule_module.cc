////////////////////////////////////////////////////////////////////////
// Class:       HitAnaModule
// Module Type: analyzer
// File:        HitAnaModule_module.cc
//
// Generated at Sun Jun 22 08:40:08 2014 by Wesley Ketchum using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

/*!
 * Title:   HitAnaModule
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Wire (calibrated), recob::Hit, Assns<recob::Wire, recob::Hit>
 * Outputs: validation histograms
 *
 * Description:
 * This module is intended to be yet another hit analyzer module. Its intention is
 *   (1) to compare hit-finding modules against each other, and eventually
 *   (2) to compare those to truth
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <vector>
#include <string>

#include "HitAnaAlg.h"

namespace hit {
  class HitAnaModule;
}

class hit::HitAnaModule : public art::EDAnalyzer {
public:
  explicit HitAnaModule(fhicl::ParameterSet const & p);
  virtual ~HitAnaModule();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.

  std::vector<std::string> fHitModuleLabels;
  std::vector<std::string> fWireModuleLabels;

};


hit::HitAnaModule::HitAnaModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{ this->reconfigure(p); }

hit::HitAnaModule::~HitAnaModule()
{
  // Clean up dynamic memory and other resources here.
}

void hit::HitAnaModule::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

void hit::HitAnaModule::beginJob()
{
  // Implementation of optional member function here.
}

void hit::HitAnaModule::endJob()
{
  // Implementation of optional member function here.
}

void hit::HitAnaModule::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.

  fHitModuleLabels = p.get< std::vector<std::string> >("HitModuleLabels");
  size_t NHitModules = fHitModuleLabels.size();

  fWireModuleLabels = p.get< std::vector<std::string> >("WireModuleLabels");
  size_t NWireModules = fWireModuleLabels.size();

  if( NWireModules==NHitModules ) return;

  if( NWireModules < NHitModules ){
    mf::LogWarning("HitAnaModule") 
      << "Number of wire module labels provided is less than hit module labels provided.\n"
      << "Using last wire module label provided (" << fWireModuleLabels[NWireModules-1] << ")"
      << " for remaining hit module labels.";

    for(size_t i=NWireModules; i<NHitModules; i++)
      fWireModuleLabels.push_back( fWireModuleLabels[NWireModules-1] );

    return;
  }

  if( NWireModules > NHitModules ){
    mf::LogWarning("HitAnaModule") 
      << "Number of hit module labels provided is less than wire module labels provided.\n"
      << "Using last hit module label provided (" << fHitModuleLabels[NHitModules-1] << ")"
      << " for remaining wire module labels.";

    for(size_t i=NHitModules; i<NWireModules; i++)
      fHitModuleLabels.push_back( fHitModuleLabels[NHitModules-1] );
    
    return;
  }

}

DEFINE_ART_MODULE(hit::HitAnaModule)

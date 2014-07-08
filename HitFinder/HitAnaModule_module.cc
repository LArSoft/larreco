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

#include "art/Framework/Services/Optional/TFileService.h"

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
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  void createAssocVector(std::vector< art::Ptr<recob::Hit> > const&,
			 std::vector< std::vector<int> >&);

  void createMCAssocVector( std::vector<recob::Wire> const&,
			    std::vector<sim::MCHitCollection> const&,
			    std::vector< std::vector<int> >&);

  // Declare member data here.
  std::vector<std::string> fHitModuleLabels;
  std::string fMCHitModuleLabel;
  std::string fWireModuleLabel;

  TTree *wireDataTree;

  HitAnaAlg analysisAlg;

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

void hit::HitAnaModule::createAssocVector( std::vector< art::Ptr<recob::Hit> > const& hitPtrs,
					   std::vector< std::vector<int> > & WireHitAssocVector){

  for(auto const& hitptr : hitPtrs)
    {
      size_t wire_key = ( hitptr->Wire() ).key();
      WireHitAssocVector.at(wire_key).push_back(hitptr.key());
    }

}

void hit::HitAnaModule::createMCAssocVector( std::vector<recob::Wire> const& wireVector,
					     std::vector<sim::MCHitCollection> const& mcHitVector,
					     std::vector< std::vector<int> > & WireMCHitAssocVector){

  WireMCHitAssocVector.clear(); WireMCHitAssocVector.resize(wireVector.size());

  //first, store all the MCHitCollection indices in a map keyed on channel
  //then, loop through wires, and lookup mchitcollections based on the wire's channel

  std::map<unsigned int,std::vector<int> > mcHitIndicesByChannel;
  for(int icol=0; icol<mcHitVector.size(); icol++)
    mcHitIndicesByChannel[mcHitVector[icol].Channel()].push_back(icol);
  

  for(int iwire=0; iwire<wireVector.size(); iwire++)
    WireMCHitAssocVector[iwire].insert(WireMCHitAssocVector[iwire].end(),
				       mcHitIndicesByChannel[wireVector[iwire].Channel()].begin(),
				       mcHitIndicesByChannel[wireVector[iwire].Channel()].end());

}

void hit::HitAnaModule::analyze(art::Event const & e)
{
  //get event and run numbers
  unsigned int eventNumber = e.id().event();
  unsigned int runNumber   = e.run();

  analysisAlg.ClearHitModules();

  //get the wire data
  art::Handle< std::vector<recob::Wire> > wireHandle;
  e.getByLabel(fWireModuleLabel,wireHandle);
  std::vector<recob::Wire> const& wireVector(*wireHandle);

  //get the MC hits
  art::Handle< std::vector<sim::MCHitCollection> > mcHitHandle;
  e.getByLabel(fMCHitModuleLabel,mcHitHandle);
  std::vector<sim::MCHitCollection> const& mcHitVector(*mcHitHandle);

  //make the association vector. First index is wire index, second is mcHitCollection index
  std::vector< std::vector<int> > WireMCHitAssocVector;
  createMCAssocVector(wireVector,mcHitVector,WireMCHitAssocVector);

  //get the hit data
  size_t nHitModules = fHitModuleLabels.size();
  std::vector< art::Handle< std::vector<recob::Hit> > > hitHandles(nHitModules);
  std::vector< std::vector< std::vector<int> > > WireHitAssocVectors(nHitModules);
  for(size_t iter=0; iter < nHitModules; iter++){

    e.getByLabel(fHitModuleLabels[iter],hitHandles[iter]);

    //create association vectors by hand for now
    std::vector< art::Ptr<recob::Hit> > hitPtrs;
    art::fill_ptr_vector(hitPtrs,hitHandles[iter]);    
    WireHitAssocVectors[iter].resize(wireVector.size());
    createAssocVector(hitPtrs,WireHitAssocVectors[iter]);

    //load in this hit/assoc pair
    analysisAlg.LoadHitAssocPair( *(hitHandles[iter]),
				  WireHitAssocVectors[iter],
				  fHitModuleLabels[iter]);
    
  }

  //run the analzyer alg
  analysisAlg.AnalyzeWires(wireVector,eventNumber,runNumber);

}

void hit::HitAnaModule::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  wireDataTree = tfs->make<TTree>("wireDataTree","WireDataTree");
  analysisAlg.SetWireDataTree(wireDataTree);
}

void hit::HitAnaModule::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.

  fHitModuleLabels  = p.get< std::vector<std::string> >("HitModuleLabels");
  fWireModuleLabel  = p.get< std::string              >("WireModuleLabel");
  fMCHitModuleLabel = p.get< std::string              >("MCHitModuleLabel");

}

DEFINE_ART_MODULE(hit::HitAnaModule)

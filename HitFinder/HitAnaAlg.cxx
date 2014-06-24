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

#include "HitAnaAlg.h"
#include <functional>
#include <iostream>

hit::HitAnaAlg::HitAnaAlg(){
  NHitModules=0;
  wireData.NHitModules = 0;
}

void hit::HitAnaAlg::LoadHitAssocPair( std::vector<recob::Hit> const& HitVector,
				       std::vector< std::vector<int> > const& AssocVector,
				       std::string const& HitModuleLabel){
  
  HitProcessingQueue.push_back( std::make_pair( std::cref(HitVector), std::cref(AssocVector)) );
  HitModuleLabels.push_back(HitModuleLabel);

  if(HitProcessingQueue.size()!=HitModuleLabels.size())
    throw hitanaalgexception;

  NHitModules = HitProcessingQueue.size();
  
}

void hit::HitAnaAlg::AnalyzeWires(std::vector<recob::Wire> const& WireVector){
  
  InitWireData();
  for(auto const& wire : WireVector)
    FillWireInfo(wire);

}

void hit::HitAnaAlg::InitWireData(){

  wireData.NHitModules = NHitModules;
  wireData.HitModuleLabels = HitModuleLabels;
  wireData.NHits.resize(NHitModules);
  wireData.Hits.resize(NHitModules);

}

void hit::HitAnaAlg::FillWireInfo(recob::Wire const& wire){

  for( auto const& range : wire.SignalROI().get_ranges() )
    ProcessROI(range);

}

void hit::HitAnaAlg::ProcessROI(lar::sparse_vector<float>::datarange_t const& range){
}

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

hit::HitAnaAlg::HitAnaAlg(){
  wireData.NHitModules = 0;
}

void hit::HitAnaAlg::SetWireDataTree(TTree *wdt){
  wireDataTree = wdt;
  SetupWireDataTree();
}

void hit::HitAnaAlg::SetupWireDataTree(){
  wireDataTree->Branch("event", &wireData.event, "event/i");
  wireDataTree->Branch("run", &wireData.run, "run/i");
  wireDataTree->Branch("channel", &wireData.channel, "channel/i");
  wireDataTree->Branch("roi_index", &wireData.range_index, "roi_index/i");
  wireDataTree->Branch("roi_start", &wireData.range_start, "roi_start/i");
  wireDataTree->Branch("roi_size", &wireData.range_size, "roi_size/i");
  wireDataTree->Branch("roi_charge", &wireData.integrated_charge, "roi_charge/F");
  wireDataTree->Branch("nHitModules", &wireData.NHitModules, "nHitModules/I");
  wireDataTree->Branch("HitModuleLabels",&wireData.HitModuleLabels);
  wireDataTree->Branch("NHits",&wireData.NHits);
  wireDataTree->Branch("Hits_IntegratedCharge",&wireData.Hits_IntegratedCharge);
}

void hit::HitAnaAlg::ClearHitModules(){
  HitModuleLabels.clear();
  HitProcessingQueue.clear();
  wireData.NHitModules = 0;
}

void hit::HitAnaAlg::LoadHitAssocPair( std::vector<recob::Hit> const& HitVector,
				       std::vector< std::vector<int> > const& AssocVector,
				       std::string const& HitModuleLabel){
  
  HitProcessingQueue.push_back( std::make_pair( std::cref(HitVector), std::cref(AssocVector)) );
  HitModuleLabels.push_back(HitModuleLabel);

  if(HitProcessingQueue.size()!=HitModuleLabels.size())
    throw hitanaalgexception;

}

void hit::HitAnaAlg::AnalyzeWires(std::vector<recob::Wire> const& WireVector,
				  unsigned int event, unsigned int run){
  
  InitWireData(event,run);
  for(size_t iwire=0 ; iwire < WireVector.size(); iwire++)
    FillWireInfo(WireVector[iwire], iwire);

}

void hit::HitAnaAlg::InitWireData(unsigned int event, unsigned int run){

  wireData.event = event;
  wireData.run = run;
  wireData.NHitModules = HitModuleLabels.size();
  wireData.HitModuleLabels = HitModuleLabels;

}

void hit::HitAnaAlg::ClearWireDataHitInfo(){
    wireData.NHits.assign(wireData.NHitModules,0);
    wireData.Hits_IntegratedCharge.assign(wireData.NHitModules,0);
    wireData.Hits.clear(); wireData.Hits.resize(wireData.NHitModules);
}

void hit::HitAnaAlg::FillWireInfo(recob::Wire const& wire, int WireIndex){

  wireData.channel = wire.Channel();
  unsigned int range_index = 0;

  for( auto const& range : wire.SignalROI().get_ranges() ){

    wireData.range_index = range_index;
    wireData.range_start = range.begin_index();
    wireData.range_size = range.size();

    ClearWireDataHitInfo();

    ProcessROI(range, WireIndex);
    range_index++;

  }//end loop over roi ranges

}

float hit::HitAnaAlg::ROIIntegral(lar::sparse_vector<float>::datarange_t const& range){

  float charge_sum=0;
  for(auto const& value : range)
    charge_sum += value;

  return charge_sum;
}

void hit::HitAnaAlg::ProcessROI(lar::sparse_vector<float>::datarange_t const& range, int WireIndex){

  wireData.integrated_charge = ROIIntegral(range);

  for(size_t iter = 0; iter < HitProcessingQueue.size(); iter++)
    FindAndStoreHitsInRange(HitProcessingQueue[iter].first, 
			    HitProcessingQueue[iter].second.at(WireIndex),
			    iter,
			    range.begin_index(),
			    range.begin_index()+range.size());

  wireDataTree->Fill();
}

void hit::HitAnaAlg::FindAndStoreHitsInRange( std::vector<recob::Hit> const& HitVector,
					      std::vector<int> const& HitsOnWire,
					      size_t hitmodule_iter,
					      size_t begin_wire_tdc,
					      size_t end_wire_tdc){

  for( auto const& hit_index : HitsOnWire){
    recob::Hit const& thishit = HitVector.at(hit_index);

    //check if this hit is on this ROI
    if( thishit.StartTime() < begin_wire_tdc ||
	thishit.EndTime() > end_wire_tdc)
      continue;

    FillHitInfo(thishit,wireData.Hits[hitmodule_iter]);
    wireData.NHits[hitmodule_iter]++;
    wireData.Hits_IntegratedCharge[hitmodule_iter] += thishit.Charge();

  }

}

void hit::HitAnaAlg::FillHitInfo(recob::Hit const& hit, std::vector<HitInfo>& HitInfoVector){
  HitInfoVector.emplace_back(hit.PeakTime(),
			     hit.SigmaPeakTime(),
			     hit.StartTime(),
			     hit.SigmaStartTime(),
			     hit.EndTime(),
			     hit.SigmaEndTime(),
			     hit.Charge(),
			     hit.SigmaCharge(),
			     hit.Charge(true),
			     hit.SigmaCharge(true),
			     hit.GoodnessOfFit());
}

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
#include <unordered_map>

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
  wireDataTree->Branch("plane", &wireData.plane, "plane/i");
  wireDataTree->Branch("roi_index", &wireData.range_index, "roi_index/i");
  wireDataTree->Branch("roi_start", &wireData.range_start, "roi_start/i");
  wireDataTree->Branch("roi_size", &wireData.range_size, "roi_size/i");
  wireDataTree->Branch("roi_charge", &wireData.integrated_charge, "roi_charge/F");
  wireDataTree->Branch("nHitModules", &wireData.NHitModules, "nHitModules/I");
  wireDataTree->Branch("HitModuleLabels",&wireData.HitModuleLabels);
  wireDataTree->Branch("NHits",&wireData.NHits);
  wireDataTree->Branch("Hits_IntegratedCharge",&wireData.Hits_IntegratedCharge);
  wireDataTree->Branch("Hits_AverageCharge",&wireData.Hits_AverageCharge);
  wireDataTree->Branch("Hits_wAverageCharge",&wireData.Hits_wAverageCharge);
  wireDataTree->Branch("Hits_wAverageTime",&wireData.Hits_wAverageTime);
  wireDataTree->Branch("NMCHits",&wireData.NMCHits);
  wireDataTree->Branch("MCHits_IntegratedCharge",&wireData.MCHits_IntegratedCharge);
  wireDataTree->Branch("MCHits_AverageCharge",&wireData.MCHits_AverageCharge);
  wireDataTree->Branch("MCHits_wAverageCharge",&wireData.MCHits_wAverageCharge);
  wireDataTree->Branch("MCHits_wAverageTime",&wireData.MCHits_wAverageTime);
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
				  std::vector<sim::MCHitCollection> const& MCHitCollectionVector,
				  std::vector< std::vector<int> > const& AssocVector,
				  util::TimeService const& TimeService,
				  unsigned int event, unsigned int run){
  
  InitWireData(event,run);
  for(size_t iwire=0 ; iwire < WireVector.size(); iwire++)
    FillWireInfo(WireVector[iwire], iwire, MCHitCollectionVector, AssocVector[iwire], TimeService);

}

void hit::HitAnaAlg::InitWireData(unsigned int event, unsigned int run){

  wireData.event = event;
  wireData.run = run;
  wireData.NHitModules = HitModuleLabels.size();
  wireData.HitModuleLabels = HitModuleLabels;

}

void hit::HitAnaAlg::ClearWireDataHitInfo(){
  wireData.NMCHits = 0;
  wireData.MCHits_IntegratedCharge = 0;
  wireData.MCHits_AverageCharge = 0;
  wireData.MCHits_wAverageCharge = 0;
  wireData.MCHits_wAverageTime = 0;
  
  wireData.NHits.assign(wireData.NHitModules,0);
  wireData.Hits_IntegratedCharge.assign(wireData.NHitModules,0);
  wireData.Hits_AverageCharge.assign(wireData.NHitModules,0);
  wireData.Hits_wAverageCharge.assign(wireData.NHitModules,0);
  wireData.Hits_wAverageTime.assign(wireData.NHitModules,0);
  wireData.Hits.clear(); wireData.Hits.resize(wireData.NHitModules);
}

void hit::HitAnaAlg::FillWireInfo(recob::Wire const& wire, 
				  int WireIndex,
				  std::vector<sim::MCHitCollection> const& MCHitCollectionVector,
				  std::vector<int> const& thisAssocVector,
				  util::TimeService const& TimeService){

  wireData.channel = wire.Channel();
  wireData.plane = wire.View();
  unsigned int range_index = 0;

  for( auto const& range : wire.SignalROI().get_ranges() ){

    wireData.range_index = range_index;
    wireData.range_start = range.begin_index();
    wireData.range_size = range.size();

    ClearWireDataHitInfo();

    ProcessROI(range, WireIndex, MCHitCollectionVector, thisAssocVector, TimeService);
    range_index++;

  }//end loop over roi ranges

}

float hit::HitAnaAlg::ROIIntegral(lar::sparse_vector<float>::datarange_t const& range){

  float charge_sum=0;
  for(auto const& value : range)
    charge_sum += value;

  return charge_sum;
}

void hit::HitAnaAlg::ProcessROI(lar::sparse_vector<float>::datarange_t const& range, 
				int WireIndex,
				std::vector<sim::MCHitCollection> const& MCHitCollectionVector,
				std::vector<int> const& thisAssocVector,
				util::TimeService const& TimeService){

  wireData.integrated_charge = ROIIntegral(range);

  //std::cout << "----------------------------------------------------------------" << std::endl;
  //std::cout << "WireIndex = " << WireIndex << std::endl;
  //std::cout << "\tRange begin: " << range.begin_index() << std::endl;
  //std::cout << "\tRange end: " << range.begin_index()+range.size() << std::endl;

  for(size_t iter = 0; iter < HitProcessingQueue.size(); iter++)
    FindAndStoreHitsInRange(HitProcessingQueue[iter].first, 
			    HitProcessingQueue[iter].second.at(WireIndex),
			    iter,
			    range.begin_index(),
			    range.begin_index()+range.size());

  FindAndStoreMCHitsInRange(MCHitCollectionVector,
			    thisAssocVector,
			    range.begin_index(),
			    range.begin_index()+range.size(),
			    TimeService);

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

    wireData.Hits_wAverageCharge[hitmodule_iter] += thishit.Charge()*thishit.Charge();
    wireData.Hits_wAverageTime[hitmodule_iter]   += thishit.Charge()*thishit.PeakTime();
  }

  wireData.Hits_AverageCharge[hitmodule_iter] = 
    wireData.Hits_IntegratedCharge[hitmodule_iter]/wireData.NHits[hitmodule_iter];
  wireData.Hits_wAverageCharge[hitmodule_iter] = 
    wireData.Hits_wAverageCharge[hitmodule_iter]/wireData.Hits_IntegratedCharge[hitmodule_iter];
  wireData.Hits_wAverageTime[hitmodule_iter] = 
    wireData.Hits_wAverageTime[hitmodule_iter]/wireData.Hits_IntegratedCharge[hitmodule_iter];

}

void hit::HitAnaAlg::FindAndStoreMCHitsInRange( std::vector<sim::MCHitCollection> const& MCHitCollectionVector,
						std::vector<int> const& HitsOnWire,
						size_t begin_wire_tdc,
						size_t end_wire_tdc,
						util::TimeService const& TimeService){

  for( auto const& hit_index : HitsOnWire){
    sim::MCHitCollection const& thismchitcol = MCHitCollectionVector.at(hit_index);

    //let's have a map to keep track of the number of total particles
    std::unordered_map<int,unsigned int> nmchits_per_trackID_map;
    for( auto const& thishit : thismchitcol){

      //std::cout << "\t************************************************************" << std::endl;
      //std::cout << "\t\tMCHit begin: " << TimeService.TPCTDC2Tick( thishit.PeakTime()-thishit.PeakWidth() ) << std::endl;
      //std::cout << "\t\tMCHit end: " << TimeService.TPCTDC2Tick( thishit.PeakTime()+thishit.PeakWidth() ) << std::endl;

      //check if this hit is on this ROI
      if( TimeService.TPCTDC2Tick( thishit.PeakTime()-thishit.PeakWidth() ) < begin_wire_tdc ||
	  TimeService.TPCTDC2Tick( thishit.PeakTime()+thishit.PeakWidth() ) > end_wire_tdc )
	continue;

      nmchits_per_trackID_map[thishit.PartTrackId()] += 1;
      wireData.MCHits_IntegratedCharge += thishit.Charge();
      
      wireData.MCHits_wAverageCharge += thishit.Charge()*thishit.Charge();
      wireData.MCHits_wAverageTime   += thishit.Charge()*TimeService.TPCTDC2Tick(thishit.PeakTime());
    }
    
    wireData.NMCHits = nmchits_per_trackID_map.size();

    wireData.MCHits_AverageCharge = 
      wireData.MCHits_IntegratedCharge/wireData.NMCHits;
    wireData.MCHits_wAverageCharge = 
      wireData.MCHits_wAverageCharge/wireData.MCHits_IntegratedCharge;
    wireData.MCHits_wAverageTime = 
      wireData.MCHits_wAverageTime/wireData.MCHits_IntegratedCharge;
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

#ifndef HITANAALG_H
#define HITANAALG_H

/*!
 * Title:   HitAnaModule
 * Author:  wketchum@lanl.gov. 
 * Inputs:  recob::Wire (calibrated), recob::Hit, Assns<recob::Wire, recob::Hit>
 * Outputs: validation histograms for wire aggregated hits.
 *
 * Embellished by echurch@fnal.gov for just Hits
 *
 * Description:
 * This module is intended to be yet another hit analyzer module. Its intention is
 *   (1) to compare hit-finding modules against each other, and eventually
 *   (2) to compare those to truth
 */

#include <vector>
#include <string>
#include <exception>

#include "MCBase/MCHitCollection.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "Utilities/IDetectorClocksService.h"

#include "TTree.h"

namespace hit{

  struct HitInfo{
    
    //need a constructor here
    HitInfo(float pt, float pt_s,
	    float w,
	    int st, int et,
	    float c, float c_s,
	    float mc, float mc_s,
	    float gof)
      : peaktime(pt)
      , peaktime_sigma(pt_s)
      , rms(w)
      , starttick(st)
      , endtick(et)
      , charge(c)
      , charge_sigma(c_s)
      , maxcharge(mc)
      , maxcharge_sigma(mc_s)
      , goodness_of_fit(gof)
    {}

    float peaktime;
    float peaktime_sigma;
    float rms;
    int starttick;
    int endtick;
    float charge;
    float charge_sigma;
    float maxcharge;
    float maxcharge_sigma;
    float goodness_of_fit;
  };

  struct WireROIInfo{
    unsigned int event;
    unsigned int run;
    unsigned int channel;
    unsigned int plane;
    unsigned int range_index;
    unsigned int range_start;
    size_t range_size;
    float integrated_charge;
    float peak_charge;
    float peak_time;
    int NHitModules;
    std::vector<std::string> HitModuleLabels;
    std::vector<int> NHits;
    std::vector<float> Hits_IntegratedCharge;
    std::vector<float> Hits_AverageCharge;
    std::vector<float> Hits_PeakCharge;
    std::vector<float> Hits_PeakTime;
    std::vector<float> Hits_wAverageCharge;
    std::vector<float> Hits_wAverageTime;
    std::vector<float> Hits_MeanMultiplicity;
    std::vector< std::vector<HitInfo> > Hits;
    int NMCHits;
    float MCHits_IntegratedCharge;
    float MCHits_AverageCharge;
    float MCHits_PeakCharge;
    float MCHits_PeakTime;
    float MCHits_wAverageCharge;
    float MCHits_wAverageTime;
  };


  class HitAnaAlgException : public std::exception{
    virtual const char* what() const throw(){
      return "HitAnaAlg Exception";
    }
  } hitanaalgexception;

  class HitAnaAlg{

    typedef std::pair< const std::vector<recob::Hit>& , const std::vector< std::vector<int> >& > HitAssocPair;

  public:

    HitAnaAlg();
    
    void SetWireDataTree(TTree*);

    void SetHitDataTree(std::vector<TTree*>& trees);

    void AnalyzeWires(std::vector<recob::Wire> const&,
		      std::vector<sim::MCHitCollection> const&,
		      std::vector< std::vector<int> > const&,
		      const dataprov::IDetectorClocks *,
		      unsigned int,
		      unsigned int);

    void LoadHitAssocPair( std::vector<recob::Hit> const&, 
			   std::vector< std::vector<int> > const&,
			   std::string const&);

    void ClearHitModules();
    
  private:
    
    void InitWireData(unsigned int, unsigned int);
    void ClearWireDataHitInfo();

    void FillHitInfo(recob::Hit const&, std::vector<HitInfo>&);

    void FillWireInfo(recob::Wire const&, 
		      int,
		      std::vector<sim::MCHitCollection> const&,
		      std::vector<int> const&,
		      const dataprov::IDetectorClocks *);

    void ProcessROI(lar::sparse_vector<float>::datarange_t const&, int,
		    std::vector<sim::MCHitCollection> const&,
		    std::vector<int> const&,
		    const dataprov::IDetectorClocks *);

    void ROIInfo(lar::sparse_vector<float>::datarange_t const&,
		 float&,float&,float&);

    void FindAndStoreHitsInRange(std::vector<recob::Hit> const&,
				 std::vector<int> const&,
				 size_t,size_t,size_t);
    void FindAndStoreMCHitsInRange(std::vector<sim::MCHitCollection> const&,
				   std::vector<int> const&,
				   size_t,size_t,
				   const dataprov::IDetectorClocks *);
    
    WireROIInfo wireData;
    std::vector<recob::Hit*> hitData;

    std::vector<std::string> HitModuleLabels;
    std::vector< HitAssocPair > HitProcessingQueue;


    void SetupWireDataTree();
    TTree* wireDataTree;

    std::vector<TTree*> hitDataTree;

    //this is for unit testing...class has no other purpose
    friend class HitAnaAlgTest;

  };

}//end namespace hit


#endif

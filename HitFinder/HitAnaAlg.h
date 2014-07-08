#ifndef HITANAALG_H
#define HITANAALG_H

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

#include <vector>
#include <string>
#include <exception>

#include "MCBase/MCHitCollection.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "Utilities/TimeService.h"

#include "TTree.h"

namespace hit{

  struct HitInfo{
    
    //need a constructor here
    HitInfo(float pt, float pt_s,
	    float st, float st_s,
	    float et, float et_s,
	    float c, float c_s,
	    float mc, float mc_s,
	    float gof)
    {
      peaktime = pt; peaktime_sigma = pt_s;
      starttime = st; starttime_sigma = st_s;
      endtime = et; endtime_sigma = et_s;
      charge = c; charge_sigma = c_s;
      maxcharge = mc; maxcharge = mc_s;
      goodness_of_fit = gof;
    }

    float peaktime;
    float peaktime_sigma;
    float starttime;
    float starttime_sigma;
    float endtime;
    float endtime_sigma;
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
    unsigned int range_index;
    unsigned int range_start;
    size_t range_size;
    float integrated_charge;
    int NHitModules;
    std::vector<std::string> HitModuleLabels;
    std::vector<int> NHits;
    std::vector<float> Hits_IntegratedCharge;
    std::vector<float> Hits_AverageCharge;
    std::vector<float> Hits_wAverageCharge;
    std::vector<float> Hits_wAverageTime;
    std::vector< std::vector<HitInfo> > Hits;
    int NMCHits;
    float MCHits_IntegratedCharge;
    float MCHits_AverageCharge;
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

    void AnalyzeWires(std::vector<recob::Wire> const&,
		      std::vector<sim::MCHitCollection> const&,
		      std::vector< std::vector<int> > const&,
		      util::TimeService const&,
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
		      util::TimeService const&);

    void ProcessROI(lar::sparse_vector<float>::datarange_t const&, int,
		    std::vector<sim::MCHitCollection> const&,
		    std::vector<int> const&,
		    util::TimeService const&);

    float ROIIntegral(lar::sparse_vector<float>::datarange_t const&);

    void FindAndStoreHitsInRange(std::vector<recob::Hit> const&,
				 std::vector<int> const&,
				 size_t,size_t,size_t);
    void FindAndStoreMCHitsInRange(std::vector<sim::MCHitCollection> const&,
				   std::vector<int> const&,
				   size_t,size_t,
				   util::TimeService const&);
    
    WireROIInfo wireData;
    std::vector<std::string> HitModuleLabels;
    std::vector< HitAssocPair > HitProcessingQueue;


    void SetupWireDataTree();
    TTree* wireDataTree;

    //this is for unit testing...class has no other purpose
    friend class HitAnaAlgTest;

  };

}//end namespace hit


#endif

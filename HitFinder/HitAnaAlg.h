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

#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"

namespace hit{

  struct HitInfo{
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
    unsigned int channel;
    unsigned int range_index;
    float integrated_charge;
    int NHitModules;
    std::vector<std::string> HitModuleLabels;
    std::vector<int> NHits;
    std::vector< std::vector<HitInfo> > Hits;
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
    
    void AnalyzeWires(std::vector<recob::Wire> const&);

    void LoadHitAssocPair( std::vector<recob::Hit> const&, 
			   std::vector< std::vector<int> > const&,
			   std::string const&);
    
  private:
    
    void InitWireData();

    void FillHitInfo(recob::Hit const&, HitInfo &);
    void FillWireInfo(recob::Wire const&);
    void ProcessROI(lar::sparse_vector<float>::datarange_t const&);
    
    WireROIInfo wireData;
    int NHitModules;
    std::vector<std::string> HitModuleLabels;
    std::vector< HitAssocPair > HitProcessingQueue;

  };

}

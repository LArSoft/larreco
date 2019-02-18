/**
 *  @file   HitsStandard_tool.cc
 * 
 *  @brief  art Tool for comparing clusters and merging those that are consistent
 * 
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/SpacePointSolver/HitReaders/IHitReader.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace reco3d {
    
class HitsStandard : virtual public IHitReader
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit HitsStandard(const fhicl::ParameterSet&);
    
    /**
     *  @brief  Destructor
     */
    ~HitsStandard();
    
    void configure(fhicl::ParameterSet const &pset) override;
    
    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    bool readHits(const std::vector<art::Ptr<recob::Hit>>&,            // input hits
                  std::vector<art::Ptr<recob::Hit>>&,                  // output hits plane 0
                  std::vector<art::Ptr<recob::Hit>>&,                  // output hits plane 1
                  std::vector<art::Ptr<recob::Hit>>&) const override;  // output hits plane 2

};

HitsStandard::HitsStandard(fhicl::ParameterSet const &pset)
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitsStandard::~HitsStandard()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitsStandard::configure(fhicl::ParameterSet const &pset)
{
//    m_enableMonitoring   = pset.get<bool>  ("EnableMonitoring",  true  );
    
    return;
}

bool HitsStandard::readHits(const std::vector<art::Ptr<recob::Hit>>& inputHits,   // input hits
                            std::vector<art::Ptr<recob::Hit>>&       xhits,       // output hits plane 0
                            std::vector<art::Ptr<recob::Hit>>&       uhits,       // output hits plane 1
                            std::vector<art::Ptr<recob::Hit>>&       vhits) const // output hits plane 2
{

    bool is2view = false;

    for(auto& hit: inputHits)
    {
        if(hit->Integral() < 0 || isnan(hit->Integral()) || isinf(hit->Integral()))
        {
            mf::LogWarning("HitsStandard") << "WARNING: bad recob::Hit::Integral() = "
            << hit->Integral()
            << ". Skipping." << std::endl;
            continue;
        }
            
        if(hit->SignalType() == geo::kCollection){
            // For DualPhase, both view are collection. Arbitrarily map V to the main
            // "X" view. For Argoneut and Lariat, collection=V is also the right
            // convention.
            if(hit->View() == geo::kZ){
                xhits.push_back(hit);
            }
            if(hit->View() == geo::kV){
                xhits.push_back(hit);
                is2view = true;
            }
            if(hit->View() == geo::kU || hit->View() == geo::kY){
                uhits.push_back(hit);
                is2view = true;
            }
        }
        else{
            if(hit->View() == geo::kU) uhits.push_back(hit);
            if(hit->View() == geo::kV) vhits.push_back(hit);
        }
    } // end for hit

    mf::LogDebug("HitsStandard") << ">>>>> Reading hits done" << std::endl;
    
    return is2view;
}
    
DEFINE_ART_CLASS_TOOL(HitsStandard)
} // namespace lar_cluster3d

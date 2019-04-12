/**
 *  @file   HitsICARUS_tool.cc
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
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace reco3d {

class HitsICARUS : virtual public IHitReader
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit HitsICARUS(const fhicl::ParameterSet&);

    /**
     *  @brief  Destructor
     */
    ~HitsICARUS();

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

HitsICARUS::HitsICARUS(fhicl::ParameterSet const &pset)
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitsICARUS::~HitsICARUS()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitsICARUS::configure(fhicl::ParameterSet const &pset)
{
    //    m_enableMonitoring   = pset.get<bool>  ("EnableMonitoring",  true  );

    return;
}

bool HitsICARUS::readHits(const std::vector<art::Ptr<recob::Hit>>& inputHits,           // input hits
                          std::vector<art::Ptr<recob::Hit>>&       collectionHits,      // output hits plane 0
                          std::vector<art::Ptr<recob::Hit>>&       firstIndHits,        // output hits plane 1
                          std::vector<art::Ptr<recob::Hit>>&       secondIndHits) const // output hits plane 2
{
    for(auto& hit: inputHits)
    {
        if(hit->Integral() < 0 || isnan(hit->Integral()) || isinf(hit->Integral()))
        {
            mf::LogWarning("Hits_ICARUS") << "WARNING: bad recob::Hit::Integral() = "
            << hit->Integral()
            << ". Skipping." << std::endl;
            continue;
        }

        if      (hit->WireID().Plane == 0) firstIndHits.push_back(hit);
        else if (hit->WireID().Plane == 1) secondIndHits.push_back(hit);
        else                               collectionHits.push_back(hit);
    } // end for hit

    mf::LogDebug("Hits_ICARUS") << ">>>>> Reading hits done" << std::endl;

    return false;
}

DEFINE_ART_CLASS_TOOL(HitsICARUS)
} // namespace lar_cluster3d

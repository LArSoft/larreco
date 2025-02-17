/**
 *  @file   HitsICARUS_tool.cc
 *
 *  @brief  art Tool for comparing clusters and merging those that are consistent
 *
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/SpacePointSolver/HitReaders/IHitReader.h"

// std includes
#include <cmath>
#include <ostream>
#include <vector>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace reco3d {

  class HitsICARUS : virtual public IHitReader {
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

    void configure(fhicl::ParameterSet const& pset) override;

    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    bool readHits(const std::vector<art::Ptr<recob::Hit>>&,           // input hits
                  std::vector<art::Ptr<recob::Hit>>&,                 // output hits plane 0
                  std::vector<art::Ptr<recob::Hit>>&,                 // output hits plane 1
                  std::vector<art::Ptr<recob::Hit>>&) const override; // output hits plane 2
  };

  HitsICARUS::HitsICARUS(fhicl::ParameterSet const& pset)
  {
    this->configure(pset);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  HitsICARUS::~HitsICARUS() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void HitsICARUS::configure(fhicl::ParameterSet const&)
  {
    //    m_enableMonitoring   = pset.get<bool>  ("EnableMonitoring",  true  );

    return;
  }

  bool HitsICARUS::readHits(
    const std::vector<art::Ptr<recob::Hit>>& inputHits,     // input hits
    std::vector<art::Ptr<recob::Hit>>& collectionHits,      // output hits plane 0
    std::vector<art::Ptr<recob::Hit>>& firstIndHits,        // output hits plane 1
    std::vector<art::Ptr<recob::Hit>>& secondIndHits) const // output hits plane 2
  {
    for (auto& hit : inputHits) {
      if (hit->Integral() < 0 || std::isnan(hit->Integral()) || std::isinf(hit->Integral())) {
        mf::LogWarning("Hits_ICARUS") << "WARNING: bad recob::Hit::Integral() = " << hit->Integral()
                                      << ". Skipping." << std::endl;
        continue;
      }

      if (hit->WireID().Plane == 0)
        firstIndHits.push_back(hit);
      else if (hit->WireID().Plane == 1)
        secondIndHits.push_back(hit);
      else
        collectionHits.push_back(hit);
    } // end for hit

    mf::LogDebug("Hits_ICARUS") << ">>>>> Reading hits done" << std::endl;

    return false;
  }

  DEFINE_ART_CLASS_TOOL(HitsICARUS)
} // namespace lar_cluster3d

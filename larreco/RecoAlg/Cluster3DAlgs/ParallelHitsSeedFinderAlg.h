/**
 *  @file   ParallelHitsSeedFinderAlg.h
 *
 *  @brief  This is an algorithm for finding recob::Seed objects in 3D clusters
 *
 */
#ifndef ParallelHitsSeedFinderAlg_h
#define ParallelHitsSeedFinderAlg_h

// Framework Includes
#include "fhiclcpp/fwd.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/SeedFinderAlgBase.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d {

  /**
 *  @brief  ParallelHitsSeedFinderAlg class
 */
  class ParallelHitsSeedFinderAlg final : public SeedFinderAlgBase {
  public:
    /**
   *  @brief  Constructor
   *
   *  @param  pset
   */
    explicit ParallelHitsSeedFinderAlg(fhicl::ParameterSet const& pset);

    /**
     *  @brief Given the list of hits this will search for candidate Seed objects and return them
     */
    bool findTrackSeeds(reco::HitPairListPtr& hitPairListPtr,
                        reco::PrincipalComponents& inputPCA,
                        SeedHitPairListPairVec& seedHitMap) const override;

  private:
    size_t m_maxNumEdgeHits; ///< Maximum number hits each end of PCA axis
    double m_gapDistance;    ///< Maximum allowed distance between hits
    size_t m_numSeed2DHits;  ///< Number 2D seed hits desired

    PrincipalComponentsAlg m_pcaAlg; // For running Principal Components Analysis
  };

} // namespace lar_cluster3d
#endif

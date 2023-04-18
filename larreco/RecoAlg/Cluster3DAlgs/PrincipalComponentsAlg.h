/**
 *  @file   PrincipalComponentsAlg.h
 *
 *  @brief  This header file defines the interface to a principal components analysis designed to
 *          be used within the 3D clustering
 *
 */
#ifndef PrincipalComponentsAlg_h
#define PrincipalComponentsAlg_h

// Framework Includes
namespace fhicl {
  class ParameterSet;
}

// LArSoft includes
namespace detinfo {
  class DetectorPropertiesData;
}
namespace geo {
  class Geometry;
}

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

// Eigen
#include <Eigen/Core>

namespace lar_cluster3d {

  /**
   *  @brief  Cluster3D class
   */
  class PrincipalComponentsAlg {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    PrincipalComponentsAlg(fhicl::ParameterSet const& pset);

    /**
     *  @brief Run the Principal Components Analysis
     */
    void PCAAnalysis(const detinfo::DetectorPropertiesData& detProp,
                     const reco::HitPairListPtr& hitPairVector,
                     reco::PrincipalComponents& pca,
                     float doca3DScl = 3.) const;

    void PCAAnalysis_3D(const reco::HitPairListPtr& hitPairList,
                        reco::PrincipalComponents& pca,
                        bool skeletonOnly = false) const;

    void PCAAnalysis_2D(const detinfo::DetectorPropertiesData& detProp,
                        const reco::HitPairListPtr& hitPairVector,
                        reco::PrincipalComponents& pca,
                        bool updateAvePos = false) const;

    void PCAAnalysis_calc3DDocas(const reco::HitPairListPtr& hitPairVector,
                                 const reco::PrincipalComponents& pca) const;

    void PCAAnalysis_calc2DDocas(const reco::Hit2DListPtr& hit2DVector,
                                 const reco::PrincipalComponents& pca) const;

    int PCAAnalysis_reject2DOutliers(const reco::HitPairListPtr& hitPairVector,
                                     reco::PrincipalComponents& pca,
                                     float aveHitDoca) const;

    int PCAAnalysis_reject3DOutliers(const reco::HitPairListPtr& hitPairVector,
                                     const reco::PrincipalComponents& pca,
                                     float aveHitDoca) const;

  private:
    float m_parallel;                ///< means lines are parallel
    const geo::Geometry* m_geometry; // pointer to the Geometry service
  };

} // namespace lar_cluster3d
#endif

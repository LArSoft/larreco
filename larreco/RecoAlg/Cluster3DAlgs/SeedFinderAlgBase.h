/**
 *  @file   SeedFinderAlgBase.h
 *
 *  @brief  This is intended to define an interface to all Seed finder algorithms employed
 *          by the 3D clustering
 *
 */
#ifndef SeedFinderAlgBase_h
#define SeedFinderAlgBase_h

#include "fhiclcpp/fwd.h"
#include "lardataobj/RecoBase/Seed.h"

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

#include <vector>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d {
  typedef std::pair<recob::Seed, reco::HitPairListPtr> SeedHitPairListPair;
  typedef std::vector<SeedHitPairListPair> SeedHitPairListPairVec;

  /**
 *  @brief  SeedFinderAlgBase class
 */
  class SeedFinderAlgBase {
  public:
    virtual ~SeedFinderAlgBase() = default;
    /**
     *  @brief Define the interface to take an input list of 3D hits and return seed candidates
     *         so hits are ordered along the axis
     */
    virtual bool findTrackSeeds(reco::HitPairListPtr& hitPairListPtr,
                                reco::PrincipalComponents& inputPCA,
                                SeedHitPairListPairVec& seedHitPairVec) const = 0;

  protected:
    /**
     *  @brief Define a comparator which will sort hits by arc length along a PCA axis
     */
    struct Sort3DHitsByArcLen3D {
      bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
      {
        return left->getArclenToPoca() < right->getArclenToPoca();
      }
    };

    /**
     *  @brief Define a comparator which will sort hits by the absolute value of arc length
     *         so hits are ordered closed to PCA origin to furthest
     */
    struct Sort3DHitsByAbsArcLen3D {
      bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
      {
        return fabs(left->getArclenToPoca()) < fabs(right->getArclenToPoca());
      }
    };

  private:
  };

} // namespace lar_cluster3d
#endif

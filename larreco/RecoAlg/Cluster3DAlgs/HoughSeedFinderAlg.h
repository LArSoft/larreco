/**
 *  @file   HoughSeedFinderAlg.h
 *
 *  @brief  This is an algorithm for finding recob::Seed objects in 3D clusters
 *
 */
#ifndef HoughSeedFinderAlg_h
#define HoughSeedFinderAlg_h

// Framework includes
namespace fhicl {
  class ParameterSet;
}

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/SeedFinderAlgBase.h"
namespace geo {
  class Geometry;
}

// ROOT includes
#include "TCanvas.h"
class TFrame;
class TVector3;
class TVirtualPad;

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d {

  /**
 *  @brief  HoughSeedFinderAlg class
 */
  class HoughSeedFinderAlg : public SeedFinderAlgBase {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    HoughSeedFinderAlg(fhicl::ParameterSet const& pset);

    /**
     *  @brief Given the list of hits this will search for candidate Seed objects and return them
     */
    bool findTrackSeeds(reco::HitPairListPtr& hitPairListPtr,
                        reco::PrincipalComponents& inputPCA,
                        SeedHitPairListPairVec& seedHitPairVec) const override;

    /**
     *  @brief Given the list of hits this will return the sets of hits which belong on the same line
     */
    bool findTrackHits(reco::HitPairListPtr& hitPairListPtr,
                       reco::PrincipalComponents& inputPCA,
                       reco::HitPairListPtrList& hitPairListPtrList) const;

  private:
    /**
     *  @brief Using Principal Components Axis, look for gaps in a list of 3D hits
     *
     *  @param inputHitList - input list of 3D hits to check
     *  @param pca          - Principal Components Axis to use
     *  @param hitListList  - output list of hit lists which are gap free
     */
    void findHitGaps(reco::HitPairListPtr& inputHitList, reco::HitPairListPtr& outputList) const;

    /**
     *  @brief Forward declaration of some of the objects necessary for hough transform
     */

    class AccumulatorBin;
    class SortHoughClusterList;
    struct SortBinIndexList;

    // Basic structure for holding our accumlator bins (to avoid a full array)
    // structure will be rho bin for first key, theta bin for second key
    typedef std::pair<int, int> BinIndex;
    typedef std::map<BinIndex, AccumulatorBin> RhoThetaAccumulatorBinMap;
    typedef std::list<BinIndex> HoughCluster;
    typedef std::list<HoughCluster> HoughClusterList;

    void HoughRegionQuery(BinIndex& curBin,
                          RhoThetaAccumulatorBinMap& rhoThetaAccumulatorBinMap,
                          HoughCluster& neighborPts,
                          size_t threshold) const;

    void expandHoughCluster(BinIndex& curBin,
                            HoughCluster& neighborPts,
                            HoughCluster& houghCluster,
                            RhoThetaAccumulatorBinMap& rhoThetaAccumulatorBinMap,
                            size_t threshold) const;

    void findHoughClusters(const reco::HitPairListPtr& inputHits,
                           reco::PrincipalComponents& pca,
                           int& nLoops,
                           RhoThetaAccumulatorBinMap& rhoThetaMap,
                           HoughClusterList& clusterList) const;

    /**
     *  @brief Given a list of candidate "seed" 3D hits, build the seed and get associated unique 2D hits
     */
    bool buildSeed(reco::HitPairListPtr& seed3DHits, SeedHitPairListPair& seedHitPair) const;

    void LineFit2DHits(std::set<const reco::ClusterHit2D*>& hitList,
                       double XOrigin,
                       TVector3& Pos,
                       TVector3& Dir,
                       double& ChiDOF) const;

    size_t m_minimum3DHits;   ///<
    int m_thetaBins;          ///<
    int m_rhoBins;            ///<
    size_t m_hiThresholdMin;  ///<
    double m_hiThresholdFrac; ///<
    double m_loThresholdFrac; ///<
    size_t m_numSeed2DHits;   ///<
    double m_numAveDocas;     ///<
    int m_numSkippedHits;     ///<
    int m_maxLoopsPerCluster; ///<
    double m_maximumGap;      ///<

    geo::Geometry const* m_geometry; // pointer to the Geometry service
    PrincipalComponentsAlg m_pcaAlg; // For running Principal Components Analysis

    bool m_displayHist;
    mutable std::vector<std::unique_ptr<TCanvas>> m_Canvases; ///< Graphical trace canvases.
    mutable std::vector<TVirtualPad*> m_Pads;                 ///< View pads in current canvas.
  };

} // namespace lar_cluster3d
#endif

/**
 *  @file   HoughSeedFinderAlg.h
 * 
 *  @brief  This is an algorithm for finding recob::Seed objects in 3D clusters
 * 
 */
#ifndef HoughSeedFinderAlg_h
#define HoughSeedFinderAlg_h

#include "RecoAlg/Cluster3DAlgs/SeedFinderAlgBase.h"
#include "RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Seed.h"
#include "RecoObjects/Cluster3D.h"

// ROOT includes
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{

/**
 *  @brief  HoughSeedFinderAlg class
 */
class HoughSeedFinderAlg : virtual public SeedFinderAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    HoughSeedFinderAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~HoughSeedFinderAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const &pset);

    /**
     *  @brief Given the list of hits this will search for candidate Seed objects and return them
     */
    virtual bool findTrackSeeds(reco::HitPairListPtr&      hitPairListPtr,
                                reco::PrincipalComponents& inputPCA,
                                SeedHitPairListPairVec&    seedHitPairVec);

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

    class  AccumulatorBin;
    class  SortHoughClusterList;
    struct SortBinIndexList;
    
    // Basic structure for holding our accumlator bins (to avoid a full array)
    // structure will be rho bin for first key, theta bin for second key
    typedef std::pair<int, int>                BinIndex;
    typedef std::map<BinIndex, AccumulatorBin> RhoThetaAccumulatorBinMap;
    typedef std::list<BinIndex>                HoughCluster;
    typedef std::list<HoughCluster >           HoughClusterList;
    
    void HoughRegionQuery(BinIndex& curBin, RhoThetaAccumulatorBinMap& rhoThetaAccumulatorBinMap, HoughCluster& neighborPts, size_t threshold) const;
    
    void expandHoughCluster(BinIndex&                  curBin,
                            HoughCluster&              neighborPts,
                            HoughCluster&              houghCluster,
                            RhoThetaAccumulatorBinMap& rhoThetaAccumulatorBinMap,
                            size_t                     threshold) const;
    
    void findHoughClusters(const reco::HitPairListPtr& inputHits,
                           reco::PrincipalComponents&  pca,
                           int&                        nLoops,
                           RhoThetaAccumulatorBinMap&  rhoThetaMap,
                           HoughClusterList&           clusterList);
    
    /**
     *  @brief Given a list of candidate "seed" 3D hits, build the seed and get associated unique 2D hits
     */
    bool buildSeed(reco::HitPairListPtr& seed3DHits, SeedHitPairListPair& seedHitPair) const;

    
    size_t                                 m_minimum3DHits;      ///<
    int                                    m_thetaBins;          ///<
    int                                    m_rhoBins;            ///<
    size_t                                 m_hiThresholdMin;     ///<
    double                                 m_hiThresholdFrac;    ///<
    double                                 m_loThresholdFrac;    ///<
    size_t                                 m_numSeed2DHits;      ///<
    double                                 m_numAveDocas;        ///<
    int                                    m_numSkippedHits;     ///<
    int                                    m_maxLoopsPerCluster; ///<
    double                                 m_maximumGap;         ///<

    geo::Geometry*                         m_geometry;           // pointer to the Geometry service
    util::DetectorProperties*              m_detector;           // Pointer to the detector properties
    
    PrincipalComponentsAlg                 m_pcaAlg;             // For running Principal Components Analysis
    
    bool                                   m_displayHist;
    std::vector<std::unique_ptr<TCanvas> > m_Canvases;           ///< Graphical trace canvases.
    std::vector<TVirtualPad*>              m_Pads;               ///< View pads in current canvas.
};

} // namespace lar_cluster3d
#endif

/**
 *  @file   MinSpanTreeAlg.h
 * 
 *  @brief  This algorithm will create and then cluster 3D hits using DBScan
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef MinSpanTreeAlg_h
#define MinSpanTreeAlg_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Hit3DBuilderAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"

// std includes
#include <vector>
#include <list>
#include <set>
#include <map>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief  MinSpanTreeAlg class definiton
 */
class MinSpanTreeAlg
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    MinSpanTreeAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~MinSpanTreeAlg();

    void reconfigure(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param hitPairClusterMap     A map of hits that have been clustered
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void Cluster3DHits(reco::HitPairList&           hitPairList,
                       reco::ClusterParametersList& clusterParametersList);
    
    /**
     *  @brief Given the results of running DBScan, format the clusters so that they can be
     *         easily transferred back to the larsoft world
     *
     *  @param hitPairClusterMap      map between view and a list of 3D hits
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     *
     *                                The last two parameters are passed through to the FillClusterParams method
     */
    void BuildClusterInfo(reco::ClusterParametersList& clusterParametersList) const;
    
    /**
     *  @brief A generic routine to actually fill the clusterParams
     *
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     */
    using Hit2DToClusterMap = std::unordered_map<const reco::ClusterHit2D*,std::unordered_map<reco::ClusterParameters*,std::set<const reco::ClusterHit3D*>>>;
    
    void FillClusterParams(reco::ClusterParameters& clusterParams, Hit2DToClusterMap& hit2DToClusterMap, double minUniqueFrac = 0., double maxLostFrac=1.) const;

    /**
     *  @brief enumerate the possible values for time checking if monitoring timing
     */
    enum TimeValues {BUILDTHREEDHITS  = 0,
                     BUILDHITTOHITMAP = 1,
                     RUNDBSCAN        = 2,
                     BUILDCLUSTERINFO = 3,
                     PATHFINDING      = 4,
                     NUMTIMEVALUES
    };
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    double getTimeToExecute(TimeValues index) const {return m_timeVector.at(index);}
    
private:
    
    class KdTreeNode;
    
    using KdTreeNodeVec  = std::vector<KdTreeNode>;
    using KdTreeNodeList = std::list<KdTreeNode>;
    using Hit3DVec       = std::vector<const reco::ClusterHit3D*>;
    
    /**
     *  @brief Given an input set of ClusterHit3D objects, build a kd tree structure
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param kdTreeVec             Container for the nodes
     */
    KdTreeNode& BuildKdTree(Hit3DVec::iterator, Hit3DVec::iterator, KdTreeNodeList&, int depth=0) const;
    
    using CandPair    = std::pair<double,const reco::ClusterHit3D*>;
    using CandPairVec = std::vector<CandPair>;
    
    size_t FindNearestNeighbors(const reco::ClusterHit3D*, const KdTreeNode&, CandPairVec&, double&) const;
    bool   FindEntry(const reco::ClusterHit3D*, const KdTreeNode&, CandPairVec&, double&, bool&, int) const;
    bool   FindEntryBrute(const reco::ClusterHit3D*, const KdTreeNode&, int) const;
    
    /**
     *  @brief The bigger question: are two pairs of hits consistent?
     */
    bool consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const;
    
    /**
     *  @brief Driver for Prim's algorithm
     */
    void RunPrimsAlgorithm(reco::HitPairList&, KdTreeNode&, reco::ClusterParametersList&) const;
    
    /**
     *  @brief Prune the obvious ambiguous hits
     */
    void PruneAmbiguousHits(reco::ClusterParameters&, Hit2DToClusterMap&) const;
    
    /**
     *  @brief Algorithm to find the best path through the given cluster
     */
    void FindBestPathInCluster(reco::ClusterParameters&) const;
    
    /**
     *  @brief a depth first search to find longest branches
     */
    reco::HitPairListPtr DepthFirstSearch(const reco::EdgeTuple&, const reco::Hit3DToEdgeMap&, double&) const;
    
    /**
     *  @brief Alternative version of FindBestPathInCluster utilizing an A* algorithm
     */
    void FindBestPathInCluster(reco::ClusterParameters&, KdTreeNode&) const;
    
    /**
     *  @brief Algorithm to find shortest path between two 3D hits
     */
    void AStar(const reco::ClusterHit3D*, const reco::ClusterHit3D*, double alpha, KdTreeNode&, reco::HitPairListPtr&, reco::EdgeList&) const;

    using BestNodeTuple = std::tuple<const reco::ClusterHit3D*,double,double>;
    using BestNodeMap   = std::unordered_map<const reco::ClusterHit3D*,BestNodeTuple>;
    
    void ReconstructBestPath(const reco::ClusterHit3D*, BestNodeMap&, reco::HitPairListPtr&, reco::EdgeList&) const;
    
    double DistanceBetweenNodes(const reco::ClusterHit3D*,const reco::ClusterHit3D*) const;
    
    /**
     *  @brief Find the lowest cost path between two nodes using MST edges
     */
    void LeastCostPath(const reco::EdgeTuple&,
                       const reco::ClusterHit3D*,
                       const reco::Hit3DToEdgeMap&,
                       reco::HitPairListPtr&,
                       reco::EdgeList&,
                       double&) const;
    
    /**
     *  @brief Given an input HitPairList, build out the map of nearest neighbors
     */
    KdTreeNode BuildKdTree(const reco::HitPairList&, KdTreeNodeList&) const;
    
    void CheckHitSorting(reco::ClusterParameters& clusterParams) const;
    
    /**
     *  @brief define data structure for keeping track of channel status
     */
    
    using ChannelStatusVec        = std::vector<size_t>;
    using ChannelStatusByPlaneVec = std::vector<ChannelStatusVec>;
    
    /**
     *  @brief Data members to follow
     */
    
    size_t                               m_minPairPts;
    double                               m_timeAdvanceGap;
    double                               m_numSigmaPeakTime;
    double                               m_pairSigmaPeakTime;
    double                               m_pairMaxDistance;
    size_t                               m_clusterMinHits;
    double                               m_clusterMinUniqueFraction;
    double                               m_clusterMaxLostFraction;
    
    bool                                 m_enableMonitoring;      ///<
    int                                  m_hits;                  ///<
    mutable std::vector<float>           m_timeVector;            ///<
    std::vector<std::vector<float>>      m_wireDir;               ///<
    
    ChannelStatusByPlaneVec              m_channelStatus;
 
    geo::Geometry*                       m_geometry;              //< pointer to the Geometry service
    const detinfo::DetectorProperties*   m_detector;              //< pointer to the detector services
    const lariov::ChannelStatusProvider* m_channelFilter;
    
    PrincipalComponentsAlg               m_pcaAlg;                // For running Principal Components Analysis
};
    
} // namespace lar_cluster3d
#endif

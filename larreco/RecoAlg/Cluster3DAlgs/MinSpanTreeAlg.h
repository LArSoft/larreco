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
 *  @brief define a kd tree node
 */
class KdTreeNode
{
public:
    enum SplitAxis { xPlane,
                     yPlane,
                     zPlane,
                     leaf,
                     null
    };
    
    KdTreeNode(SplitAxis axis, double axisVal, const KdTreeNode& left, const KdTreeNode& right) :
               m_splitAxis(axis),
               m_axisValue(axisVal),
               m_clusterHit3D(0),
               m_leftTree(left),
               m_rightTree(right)
    {}
    
    KdTreeNode(const reco::ClusterHit3D* hit) :
               m_splitAxis(SplitAxis::leaf),
               m_axisValue(0.),
               m_clusterHit3D(hit),
               m_leftTree(*this),
               m_rightTree(*this)
    {}
    
    KdTreeNode() : m_splitAxis(SplitAxis::null),
    m_axisValue(0.),
    m_clusterHit3D(0),
    m_leftTree(*this),
    m_rightTree(*this)
    {}
    
    bool                      isLeafNode()      const {return m_splitAxis == SplitAxis::leaf;}
    bool                      isNullNode()      const {return m_splitAxis == SplitAxis::null;}
    
    SplitAxis                 getSplitAxis()    const {return m_splitAxis;}
    double                    getAxisValue()    const {return m_axisValue;}
    const reco::ClusterHit3D* getClusterHit3D() const {return m_clusterHit3D;}
    const KdTreeNode&         leftTree()        const {return m_leftTree;}
    const KdTreeNode&         rightTree()       const {return m_rightTree;}
    
private:
    
    SplitAxis                 m_splitAxis;
    double                    m_axisValue;
    const reco::ClusterHit3D* m_clusterHit3D;
    const KdTreeNode&         m_leftTree;
    const KdTreeNode&         m_rightTree;
};
    
using KdTreeNodeVec  = std::vector<KdTreeNode>;
using KdTreeNodeList = std::list<KdTreeNode>;
using Hit3DVec       = std::vector<const reco::ClusterHit3D*>;
    
/**
 *  @brief a utility class for keeping track of the state of a hit for DBScan
 */
class MSTScanParams
{
public:
    MSTScanParams() : m_visited(false), m_noise(false), m_inCluster(false), m_count(0) {}
        
    void setVisited()         {m_visited   = true;}
    void setNoise()           {m_noise     = true;}
    void setInCluster()       {m_inCluster = true;}
    void setCount(int count)  {m_count     = count;}
    
    void clearVisited()                   const {m_visited   = false;}
    void incrementCount(size_t count = 1) const {m_count    += count;}
        
    bool   visited()   const {return m_visited;}
    bool   isNoise()   const {return m_noise;}
    bool   inCluster() const {return m_inCluster;}
    size_t getCount()  const {return m_count;}
        
private:
    mutable bool   m_visited;
    bool           m_noise;
    bool           m_inCluster;
    mutable size_t m_count;
};

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
    void ClusterHitsDBScan(reco::HitPairList&           hitPairList,
                           reco::HitPairClusterMap&     hitPairClusterMap,
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
    void BuildClusterInfo(reco::HitPairClusterMap&     hitPairClusterMap,
                          reco::ClusterParametersList& clusterParametersList) const;
    
    /**
     *  @brief A generic routine to actually fill the clusterParams
     *
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     */
    void FillClusterParams(reco::ClusterParameters& clusterParams, double minUniqueFrac = 0., double maxLostFrac=1.) const;

    /**
     *  @brief enumerate the possible values for time checking if monitoring timing
     */
    enum TimeValues {BUILDTHREEDHITS  = 0,
                     BUILDHITTOHITMAP = 1,
                     RUNDBSCAN        = 2,
                     BUILDCLUSTERINFO = 3,
                     NUMTIMEVALUES
    };
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    double getTimeToExecute(TimeValues index) const {return m_timeVector.at(index);}
    
private:
    
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
     *  @brief a depth first search to find longest branches
     */
    using MSTEdgeTuple        = std::tuple<const reco::ClusterHit3D*,const reco::ClusterHit3D*,double>;
    using MSTEdgeList         = std::list<MSTEdgeTuple>;
    using MST3DHitToEdgePair  = std::pair<const reco::ClusterHit3D*, MSTEdgeList>;
    using MST3DHitToEdgeMap   = std::unordered_map<const reco::ClusterHit3D*, MSTEdgeList>;
    using MSTClusterToEdgeMap = std::map<int,MST3DHitToEdgeMap>;
    
    reco::HitPairListPtr DepthFirstSearch(const MSTEdgeTuple&, const MST3DHitToEdgeMap&, double&) const;
    
    /**
     *  @brief The bigger question: are two pairs of hits consistent?
     */
    bool consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const;
    
    typedef std::list<const reco::ClusterHit3D*>                           EpsPairNeighborhoodList;
    typedef std::pair<MSTScanParams, EpsPairNeighborhoodList >             EpsPairNeighborhoodPair;
    typedef std::map<const reco::ClusterHit3D*, EpsPairNeighborhoodPair >  EpsPairNeighborhoodMap;
    typedef std::pair<const reco::ClusterHit3D*, EpsPairNeighborhoodPair > EpsPairNeighborhoodMapPair;
    typedef std::vector<EpsPairNeighborhoodMapPair >                       EpsPairNeighborhoodMapVec;
    
    /**
     *  @brief Driver for Prim's algorithm
     */
    void RunPrimsAlgorithm(reco::HitPairList&, KdTreeNode&, reco::HitPairClusterMap&) const;
    
    /**
     *  @brief Driver for DBScan
     */
    void RunDBScan(EpsPairNeighborhoodMapVec&, reco::HitPairClusterMap&) const;
    
    /**
     *  @brief The primary cluster building section for DBScan
     */
    void expandCluster(EpsPairNeighborhoodMapVec&, EpsPairNeighborhoodMapVec::iterator, reco::HitPairListPtr&, size_t) const;
    
    /**
     *  @brief Given an input HitPairList, build out the map of nearest neighbors
     */
    KdTreeNode BuildKdTree(const reco::HitPairList&, KdTreeNodeList&) const;
    
    /**
     *  @brief Given an input HitPairList, build out the map of nearest neighbors
     */
    size_t BuildNeighborhoodMap(const reco::HitPairList&, EpsPairNeighborhoodMapVec&) const;
    
    /**
     *  @brief define data structure for keeping track of channel status
     */
    
    using ChannelStatusVec       = std::vector<size_t>;
    using ChannelStatusByViewVec = std::vector<ChannelStatusVec>;
    
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
    double                               m_wirePitch[3];
    mutable std::vector<float>           m_timeVector;            ///<
    std::vector<std::vector<double>>     m_wireDir;               ///<
    std::vector<std::vector<double>>     m_wireNormal;            ///<
    
    ChannelStatusByViewVec               m_channelStatus;
 
    geo::Geometry*                       m_geometry;              //< pointer to the Geometry service
    const detinfo::DetectorProperties*   m_detector;              //< pointer to the detector services
    const lariov::ChannelStatusProvider* m_channelFilter;
    
    PrincipalComponentsAlg               m_pcaAlg;                // For running Principal Components Analysis
};
    
} // namespace lar_cluster3d
#endif

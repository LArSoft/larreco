/**
 *  @file   kdTree.h
 * 
 *  @brief  Implements a kdTree for use in clustering
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef kdTree_h
#define kdTree_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "lardata/RecoObjects/Cluster3D.h"

// std includes
#include <vector>
#include <list>
#include <set>
#include <map>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief  kdTree class definiton
 */
class kdTree
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    kdTree(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~kdTree();

    /**
     *  @brief Configure our kdTree...
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    void configure(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief Define some data structures useful for the kd tree
     */
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
    
    using CandPair     = std::pair<double,const reco::ClusterHit3D*>;
    using CandPairList = std::list<CandPair>;
    
    size_t FindNearestNeighbors(const reco::ClusterHit3D*, const KdTreeNode&, CandPairList&, double&) const;
    bool   FindEntry(const reco::ClusterHit3D*, const KdTreeNode&, CandPairList&, double&, bool&, int) const;
    bool   FindEntryBrute(const reco::ClusterHit3D*, const KdTreeNode&, int) const;
    
    /**
     *  @brief Given an input HitPairList, build out the map of nearest neighbors
     */
    KdTreeNode BuildKdTree(const reco::HitPairList&, KdTreeNodeList&) const;
    
    float getTimeToExecute() const {return m_timeToBuild;}
    
private:
    
    /**
     *  @brief The bigger question: are two pairs of hits consistent?
     */
    bool consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const;
    
    double DistanceBetweenNodes(const reco::ClusterHit3D*,const reco::ClusterHit3D*) const;
    
    bool           m_enableMonitoring;      ///<
    mutable float  m_timeToBuild;           ///<
    double         m_pairSigmaPeakTime;

};
    
/**
 *  @brief define a kd tree node
 */
class kdTree::KdTreeNode
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
    
} // namespace lar_cluster3d
#endif

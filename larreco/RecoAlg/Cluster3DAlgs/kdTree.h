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
#include "fhiclcpp/fwd.h"

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

// std includes
#include <list>
#include <utility>
#include <vector>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d {
  /**
 *  @brief  kdTree class definiton
 */
  class kdTree {
  public:
    /**
     *  @brief  Default Constructor
     */
    kdTree()
      : fEnableMonitoring(false)
      , fTimeToBuild(0.)
      , fPairSigmaPeakTime(0.)
      , fRefLeafBestDist(0.)
      , fMaxWireDeltas(0)
    {}

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    kdTree(fhicl::ParameterSet const& pset);

    /**
     *  @brief  Destructor
     */
    ~kdTree();

    /**
     *  @brief Configure our kdTree...
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    void configure(fhicl::ParameterSet const& pset);

    /**
     *  @brief Define some data structures useful for the kd tree
     */
    class KdTreeNode;

    using KdTreeNodeVec = std::vector<KdTreeNode>;
    using KdTreeNodeList = std::list<KdTreeNode>;
    using Hit3DVec = std::vector<const reco::ClusterHit3D*>;

    /**
     *  @brief Given an input set of ClusterHit3D objects, build a kd tree structure
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param kdTreeVec             Container for the nodes
     */
    KdTreeNode& BuildKdTree(Hit3DVec::iterator,
                            Hit3DVec::iterator,
                            KdTreeNodeList&,
                            int depth = 0) const;

    using CandPair = std::pair<double, const reco::ClusterHit3D*>;
    using CandPairList = std::list<CandPair>;

    size_t FindNearestNeighbors(const reco::ClusterHit3D*,
                                const KdTreeNode&,
                                CandPairList&,
                                float&) const;
    bool FindEntry(const reco::ClusterHit3D*, const KdTreeNode&, CandPairList&, float&, bool&, int)
      const;
    bool FindEntryBrute(const reco::ClusterHit3D*, const KdTreeNode&, int) const;

    /**
     *  @brief Given an input HitPairList, build out the map of nearest neighbors
     */
    KdTreeNode BuildKdTree(const reco::HitPairList&, KdTreeNodeList&) const;

    /**
     *  @brief Given an input HitPairListPtr, build out the map of nearest neighbors
     */
    KdTreeNode BuildKdTree(const reco::HitPairListPtr&, KdTreeNodeList&) const;

    float getTimeToExecute() const { return fTimeToBuild; }

  private:
    /**
     *  @brief The bigger question: are two pairs of hits consistent?
     */
    bool consistentPairs(const reco::ClusterHit3D* pair1,
                         const reco::ClusterHit3D* pair2,
                         float& hitSeparation) const;

    float DistanceBetweenNodes(const reco::ClusterHit3D*, const reco::ClusterHit3D*) const;
    float DistanceBetweenNodesYZ(const reco::ClusterHit3D*, const reco::ClusterHit3D*) const;

    bool fEnableMonitoring;     ///<
    mutable float fTimeToBuild; ///<
    float fPairSigmaPeakTime;   ///< Consider hits consistent if "significance" less than this
    float fRefLeafBestDist;     ///< Set neighborhood distance to this when ref leaf found
    int fMaxWireDeltas;         ///< Maximum total number of delta wires
  };

  /**
 *  @brief define a kd tree node
 */
  class kdTree::KdTreeNode {
  public:
    enum SplitAxis { xPlane, yPlane, zPlane, leaf, null };

    KdTreeNode(SplitAxis axis, float axisVal, const KdTreeNode& left, const KdTreeNode& right)
      : m_splitAxis(axis)
      , m_axisValue(axisVal)
      , m_clusterHit3D(0)
      , m_leftTree(left)
      , m_rightTree(right)
    {}

    KdTreeNode(const reco::ClusterHit3D* hit)
      : m_splitAxis(SplitAxis::leaf)
      , m_axisValue(0.)
      , m_clusterHit3D(hit)
      , m_leftTree(*this)
      , m_rightTree(*this)
    {}

    KdTreeNode()
      : m_splitAxis(SplitAxis::null)
      , m_axisValue(0.)
      , m_clusterHit3D(0)
      , m_leftTree(*this)
      , m_rightTree(*this)
    {}

    bool isLeafNode() const { return m_splitAxis == SplitAxis::leaf; }
    bool isNullNode() const { return m_splitAxis == SplitAxis::null; }

    SplitAxis getSplitAxis() const { return m_splitAxis; }
    float getAxisValue() const { return m_axisValue; }
    const reco::ClusterHit3D* getClusterHit3D() const { return m_clusterHit3D; }
    const KdTreeNode& leftTree() const { return m_leftTree; }
    const KdTreeNode& rightTree() const { return m_rightTree; }

  private:
    SplitAxis m_splitAxis;
    float m_axisValue;
    const reco::ClusterHit3D* m_clusterHit3D;
    const KdTreeNode& m_leftTree;
    const KdTreeNode& m_rightTree;
  };

} // namespace lar_cluster3d
#endif

/**
 *  @file   MinSpanTreeAlg.cxx
 *
 *  @brief  Producer module to create 3D clusters from input hits
 *
 */

// Framework Includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterParamsBuilder.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/kdTree.h"

// std includes
#include <iostream>
#include <memory>
#include <unordered_map>

// Eigen includes
#include <Eigen/Core>

#include "TVector3.h"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

  class MinSpanTreeAlg : public IClusterAlg {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit MinSpanTreeAlg(const fhicl::ParameterSet&);

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param hitPairClusterMap     A map of hits that have been clustered
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void Cluster3DHits(reco::HitPairList& hitPairList,
                       reco::ClusterParametersList& clusterParametersList) const override;

    void Cluster3DHits(reco::HitPairListPtr& /* hitPairList */,
                       reco::ClusterParametersList& /* clusterParametersList */) const override
    {}

    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute(TimeValues index) const override { return m_timeVector.at(index); }

  private:
    /**
     *  @brief Driver for Prim's algorithm
     */
    void RunPrimsAlgorithm(reco::HitPairList&,
                           kdTree::KdTreeNode&,
                           reco::ClusterParametersList&) const;

    /**
     *  @brief Prune the obvious ambiguous hits
     */
    void PruneAmbiguousHits(reco::ClusterParameters&, reco::Hit2DToClusterMap&) const;

    /**
     *  @brief Algorithm to find the best path through the given cluster
     */
    void FindBestPathInCluster(reco::ClusterParameters&) const;

    /**
     *  @brief a depth first search to find longest branches
     */
    reco::HitPairListPtr DepthFirstSearch(const reco::EdgeTuple&,
                                          const reco::Hit3DToEdgeMap&,
                                          float&) const;

    /**
     *  @brief Alternative version of FindBestPathInCluster utilizing an A* algorithm
     */
    void FindBestPathInCluster(reco::ClusterParameters&, kdTree::KdTreeNode&) const;

    /**
     *  @brief Algorithm to find shortest path between two 3D hits
     */
    void AStar(const reco::ClusterHit3D*,
               const reco::ClusterHit3D*,
               reco::ClusterParameters&) const;

    using BestNodeTuple = std::tuple<const reco::ClusterHit3D*, float, float>;
    using BestNodeMap = std::unordered_map<const reco::ClusterHit3D*, BestNodeTuple>;

    void ReconstructBestPath(const reco::ClusterHit3D*,
                             BestNodeMap&,
                             reco::HitPairListPtr&,
                             reco::EdgeList&) const;

    float DistanceBetweenNodes(const reco::ClusterHit3D*, const reco::ClusterHit3D*) const;

    /**
     *  @brief Find the lowest cost path between two nodes using MST edges
     */
    void LeastCostPath(const reco::EdgeTuple&,
                       const reco::ClusterHit3D*,
                       reco::ClusterParameters&,
                       float&) const;

    void CheckHitSorting(reco::ClusterParameters& clusterParams) const;

    /**
     *  @brief define data structure for keeping track of channel status
     */
    using ChannelStatusVec = std::vector<size_t>;
    using ChannelStatusByPlaneVec = std::vector<ChannelStatusVec>;

    /**
     *  @brief Data members to follow
     */
    bool m_enableMonitoring;                   ///<
    mutable std::vector<float> m_timeVector;   ///<
    std::vector<std::vector<float>> m_wireDir; ///<

    geo::GeometryCore const* m_geometry; //< pointer to the Geometry provider
    geo::WireReadoutGeom const* m_wireReadoutGeom;

    PrincipalComponentsAlg m_pcaAlg; // For running Principal Components Analysis
    kdTree m_kdTree;                 // For the kdTree

    std::unique_ptr<lar_cluster3d::IClusterParametersBuilder>
      m_clusterBuilder; ///<  Common cluster builder tool
  };

  MinSpanTreeAlg::MinSpanTreeAlg(fhicl::ParameterSet const& pset)
    : m_enableMonitoring{pset.get<bool>("EnableMonitoring", true)}
    , m_geometry{art::ServiceHandle<geo::Geometry const>{}.get()}
    , m_wireReadoutGeom{&art::ServiceHandle<geo::WireReadout const>()->Get()}
    , m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
    , m_kdTree(pset.get<fhicl::ParameterSet>("kdTree"))
  {
    m_timeVector.resize(NUMTIMEVALUES, 0.);

    // Determine the unit directon and normal vectors to the wires
    m_wireDir.resize(3);

    raw::ChannelID_t uChannel(0);
    std::vector<geo::WireID> uWireID = m_wireReadoutGeom->ChannelToWire(uChannel);
    const geo::WireGeo* uWireGeo = m_wireReadoutGeom->WirePtr(uWireID[0]);

    auto const uWireDir = uWireGeo->Direction();
    m_wireDir[0] = {(float)uWireDir.X(), (float)-uWireDir.Z(), (float)uWireDir.Y()};

    raw::ChannelID_t vChannel(2400);
    std::vector<geo::WireID> vWireID = m_wireReadoutGeom->ChannelToWire(vChannel);
    const geo::WireGeo* vWireGeo = m_wireReadoutGeom->WirePtr(vWireID[0]);

    auto const vWireDir = vWireGeo->Direction();
    m_wireDir[1] = {(float)vWireDir.X(), (float)-vWireDir.Z(), (float)vWireDir.Y()};
    m_wireDir[2] = {0., 0., 1.};

    m_clusterBuilder = art::make_tool<lar_cluster3d::IClusterParametersBuilder>(
      pset.get<fhicl::ParameterSet>("ClusterParamsBuilder"));
  }

  void MinSpanTreeAlg::Cluster3DHits(reco::HitPairList& hitPairList,
                                     reco::ClusterParametersList& clusterParametersList) const
  {
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */

    // Zero the time vector
    if (m_enableMonitoring) std::fill(m_timeVector.begin(), m_timeVector.end(), 0.);

    // DBScan is driven of its "epsilon neighborhood". Computing adjacency within DBScan can be time
    // consuming so the idea is the prebuild the adjaceny map and then run DBScan.
    // The following call does this work
    kdTree::KdTreeNodeList kdTreeNodeContainer;
    kdTree::KdTreeNode topNode = m_kdTree.BuildKdTree(hitPairList, kdTreeNodeContainer);

    if (m_enableMonitoring) m_timeVector.at(BUILDHITTOHITMAP) = m_kdTree.getTimeToExecute();

    // Run DBScan to get candidate clusters
    RunPrimsAlgorithm(hitPairList, topNode, clusterParametersList);

    // Initial clustering is done, now trim the list and get output parameters
    cet::cpu_timer theClockBuildClusters;

    // Start clocks if requested
    if (m_enableMonitoring) theClockBuildClusters.start();

    m_clusterBuilder->BuildClusterInfo(clusterParametersList);

    if (m_enableMonitoring) {
      theClockBuildClusters.stop();

      m_timeVector[BUILDCLUSTERINFO] = theClockBuildClusters.accumulated_real_time();
    }

    // Test run the path finding algorithm
    for (auto& clusterParams : clusterParametersList)
      FindBestPathInCluster(clusterParams, topNode);

    mf::LogDebug("MinSpanTreeAlg") << ">>>>> Cluster3DHits done, found "
                                   << clusterParametersList.size() << " clusters" << std::endl;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
  void MinSpanTreeAlg::RunPrimsAlgorithm(reco::HitPairList& hitPairList,
                                         kdTree::KdTreeNode& topNode,
                                         reco::ClusterParametersList& clusterParametersList) const
  {
    // If no hits then no work
    if (hitPairList.empty()) return;

    // Now proceed with building the clusters
    cet::cpu_timer theClockDBScan;

    // Start clocks if requested
    if (m_enableMonitoring) theClockDBScan.start();

    // Initialization
    size_t clusterIdx(0);

    // This will contain our list of edges
    reco::EdgeList curEdgeList;

    // Get the first point
    reco::HitPairList::iterator freeHitItr = hitPairList.begin();
    const reco::ClusterHit3D* lastAddedHit = &(*freeHitItr++);

    lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);

    // Make a cluster...
    clusterParametersList.push_back(reco::ClusterParameters());

    // Get an iterator to the first cluster
    reco::ClusterParametersList::iterator curClusterItr = --clusterParametersList.end();

    // We use pointers here because the objects they point to will change in the loop below
    reco::Hit3DToEdgeMap* curEdgeMap = &(*curClusterItr).getHit3DToEdgeMap();
    reco::HitPairListPtr* curCluster = &(*curClusterItr).getHitPairListPtr();

    // Loop until all hits have been associated to a cluster
    while (1) {
      // and the 3D hit status bits
      lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);

      // Purge the current list to get rid of edges which point to hits already in the cluster
      for (reco::EdgeList::iterator curEdgeItr = curEdgeList.begin();
           curEdgeItr != curEdgeList.end();) {
        if (std::get<1>(*curEdgeItr)->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)
          curEdgeItr = curEdgeList.erase(curEdgeItr);
        else
          curEdgeItr++;
      }

      // Add the lastUsedHit to the current cluster
      curCluster->push_back(lastAddedHit);

      // Set up to find the list of nearest neighbors to the last used hit...
      kdTree::CandPairList CandPairList;
      float bestDistance(1.5); //std::numeric_limits<float>::max());

      // And find them... result will be an unordered list of neigbors
      m_kdTree.FindNearestNeighbors(lastAddedHit, topNode, CandPairList, bestDistance);

      // Copy edges to the current list (but only for hits not already in a cluster)
      for (auto& pair : CandPairList) {
        if (!(pair.second->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)) {
          double edgeWeight = lastAddedHit->getHitChiSquare() * pair.second->getHitChiSquare();

          curEdgeList.push_back(reco::EdgeTuple(lastAddedHit, pair.second, edgeWeight));
        }
      }

      // If the edge list is empty then we have a complete cluster
      if (curEdgeList.empty()) {
        std::cout << "-----------------------------------------------------------------------------"
                     "------------"
                  << std::endl;
        std::cout << "**> Cluster idx: " << clusterIdx++ << " has " << curCluster->size() << " hits"
                  << std::endl;

        // Look for the next "free" hit
        freeHitItr = std::find_if(freeHitItr, hitPairList.end(), [](const auto& hit) {
          return !(hit.getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED);
        });

        // If at end of input list we are done with all hits
        if (freeHitItr == hitPairList.end()) break;

        std::cout << "##################################################################>"
                     "Processing another cluster"
                  << std::endl;

        // Otherwise, get a new cluster and set up
        clusterParametersList.push_back(reco::ClusterParameters());

        curClusterItr = --clusterParametersList.end();

        curEdgeMap = &(*curClusterItr).getHit3DToEdgeMap();
        curCluster = &(*curClusterItr).getHitPairListPtr();
        lastAddedHit = &(*freeHitItr++);
      }
      // Otherwise we are still processing the current cluster
      else {
        // Sort the list of edges by distance
        curEdgeList.sort([](const auto& left, const auto& right) {
          return std::get<2>(left) < std::get<2>(right);
        });

        // Populate the map with the edges...
        reco::EdgeTuple& curEdge = curEdgeList.front();

        (*curEdgeMap)[std::get<0>(curEdge)].push_back(curEdge);
        (*curEdgeMap)[std::get<1>(curEdge)].push_back(
          reco::EdgeTuple(std::get<1>(curEdge), std::get<0>(curEdge), std::get<2>(curEdge)));

        // Update the last hit to be added to the collection
        lastAddedHit = std::get<1>(curEdge);
      }
    }

    if (m_enableMonitoring) {
      theClockDBScan.stop();

      m_timeVector[RUNDBSCAN] = theClockDBScan.accumulated_real_time();
    }
  }

  void MinSpanTreeAlg::FindBestPathInCluster(reco::ClusterParameters& curCluster) const
  {
    reco::HitPairListPtr longestCluster;
    float bestQuality(0.);
    float aveNumEdges(0.);
    size_t maxNumEdges(0);
    size_t nIsolatedHits(0);

    // Now proceed with building the clusters
    cet::cpu_timer theClockPathFinding;

    // Start clocks if requested
    if (m_enableMonitoring) theClockPathFinding.start();

    reco::HitPairListPtr& hitPairList = curCluster.getHitPairListPtr();
    reco::Hit3DToEdgeMap& curEdgeMap = curCluster.getHit3DToEdgeMap();
    reco::EdgeList& bestEdgeList = curCluster.getBestEdgeList();

    // Do some spelunking...
    for (const auto& hit : hitPairList) {
      if (!curEdgeMap[hit].empty() && curEdgeMap[hit].size() == 1) {
        float quality(0.);

        reco::HitPairListPtr tempList =
          DepthFirstSearch(curEdgeMap[hit].front(), curEdgeMap, quality);

        tempList.push_front(std::get<0>(curEdgeMap[hit].front()));

        if (quality > bestQuality) {
          longestCluster = tempList;
          bestQuality = quality;
        }

        nIsolatedHits++;
      }

      aveNumEdges += float(curEdgeMap[hit].size());
      maxNumEdges = std::max(maxNumEdges, curEdgeMap[hit].size());
    }

    aveNumEdges /= float(hitPairList.size());
    std::cout << "----> # isolated hits: " << nIsolatedHits
              << ", longest branch: " << longestCluster.size()
              << ", cluster size: " << hitPairList.size() << ", ave # edges: " << aveNumEdges
              << ", max: " << maxNumEdges << std::endl;

    if (!longestCluster.empty()) {
      hitPairList = longestCluster;
      for (const auto& hit : hitPairList) {
        for (const auto& edge : curEdgeMap[hit])
          bestEdgeList.emplace_back(edge);
      }

      std::cout << "        ====> new cluster size: " << hitPairList.size() << std::endl;
    }

    if (m_enableMonitoring) {
      theClockPathFinding.stop();

      m_timeVector[PATHFINDING] += theClockPathFinding.accumulated_real_time();
    }
  }

  void MinSpanTreeAlg::FindBestPathInCluster(reco::ClusterParameters& clusterParams,
                                             kdTree::KdTreeNode& /* topNode */) const
  {
    // Set up for timing the function
    cet::cpu_timer theClockPathFinding;

    // Start clocks if requested
    if (m_enableMonitoring) theClockPathFinding.start();

    // Trial A* here
    if (clusterParams.getHitPairListPtr().size() > 2) {
      // Get references to what we need....
      reco::HitPairListPtr& curCluster = clusterParams.getHitPairListPtr();
      reco::Hit3DToEdgeMap& curEdgeMap = clusterParams.getHit3DToEdgeMap();

      // Do a quick PCA to determine our parameter "alpha"
      reco::PrincipalComponents pca;
      m_pcaAlg.PCAAnalysis_3D(curCluster, pca);

      // The chances of a failure are remote, still we should check
      if (pca.getSvdOK()) {
        float pcaLen = 3.0 * sqrt(pca.getEigenValues()[2]);
        float pcaWidth = 3.0 * sqrt(pca.getEigenValues()[1]);
        float pcaHeight = 3.0 * sqrt(pca.getEigenValues()[0]);
        const Eigen::Vector3f& pcaCenter = pca.getAvePosition();
        float alpha = std::min(float(1.), std::max(float(0.001), pcaWidth / pcaLen));

        // Create a temporary container for the isolated points
        reco::ProjectedPointList isolatedPointList;

        // Go through and find the isolated points, for those get the projection to the plane of maximum spread
        for (const auto& hit3D : curCluster) {
          // the definition of an isolated hit is that it only has one associated edge
          if (!curEdgeMap[hit3D].empty() && curEdgeMap[hit3D].size() == 1) {
            Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
                                        hit3D->getPosition()[1] - pcaCenter(1),
                                        hit3D->getPosition()[2] - pcaCenter(2));
            Eigen::Vector3f pcaToHit = pca.getEigenVectors() * pcaToHitVec;

            // This sets x,y where x is the longer spread, y the shorter
            isolatedPointList.emplace_back(pcaToHit(2), pcaToHit(1), hit3D);
          }
        }

        std::cout << "************* Finding best path with A* in cluster *****************"
                  << std::endl;
        std::cout << "**> There are " << curCluster.size() << " hits, " << isolatedPointList.size()
                  << " isolated hits, the alpha parameter is " << alpha << std::endl;
        std::cout << "**> PCA len: " << pcaLen << ", wid: " << pcaWidth << ", height: " << pcaHeight
                  << ", ratio: " << pcaHeight / pcaWidth << std::endl;

        // If no isolated points then nothing to do...
        if (isolatedPointList.size() > 1) {
          // Sort the point vec by increasing x, if same then by increasing y.
          isolatedPointList.sort([](const auto& left, const auto& right) {
            return (std::abs(std::get<0>(left) - std::get<0>(right)) >
                    std::numeric_limits<float>::epsilon()) ?
                     std::get<0>(left) < std::get<0>(right) :
                     std::get<1>(left) < std::get<1>(right);
          });

          // Ok, get the two most distance points...
          const reco::ClusterHit3D* startHit = std::get<2>(isolatedPointList.front());
          const reco::ClusterHit3D* stopHit = std::get<2>(isolatedPointList.back());

          std::cout << "**> Sorted " << isolatedPointList.size()
                    << " hits, longest distance: " << DistanceBetweenNodes(startHit, stopHit)
                    << std::endl;

          float cost(std::numeric_limits<float>::max());

          LeastCostPath(curEdgeMap[startHit].front(), stopHit, clusterParams, cost);

          clusterParams.getBestHitPairListPtr().push_front(startHit);

          std::cout << "**> Best path has " << clusterParams.getBestHitPairListPtr().size()
                    << " hits, " << clusterParams.getBestEdgeList().size() << " edges" << std::endl;
        }
      }
      else {
        std::cout << "++++++>>> PCA failure! # hits: " << curCluster.size() << std::endl;
      }
    }

    if (m_enableMonitoring) {
      theClockPathFinding.stop();

      m_timeVector[PATHFINDING] += theClockPathFinding.accumulated_real_time();
    }
  }

  void MinSpanTreeAlg::AStar(const reco::ClusterHit3D* startNode,
                             const reco::ClusterHit3D* goalNode,
                             reco::ClusterParameters& clusterParams) const
  {
    // Recover the list of hits and edges
    reco::HitPairListPtr& pathNodeList = clusterParams.getBestHitPairListPtr();
    reco::EdgeList& bestEdgeList = clusterParams.getBestEdgeList();
    reco::Hit3DToEdgeMap& curEdgeMap = clusterParams.getHit3DToEdgeMap();

    // Find the shortest path from start to goal using an A* algorithm
    // Keep track of the nodes which have already been evaluated
    reco::HitPairListPtr closedList;

    // Keep track of the nodes that have been "discovered" but yet to be evaluated
    reco::HitPairListPtr openList = {startNode};

    // Create a map which, for each node, will tell us the node it can be most efficiencly reached from.
    BestNodeMap bestNodeMap;

    bestNodeMap[startNode] =
      BestNodeTuple(startNode, 0., DistanceBetweenNodes(startNode, goalNode));

    while (!openList.empty()) {
      // The list is not empty so by def we will return something
      reco::HitPairListPtr::iterator currentNodeItr = openList.begin();

      // If the list contains more than one element then we need to find the one with the smallest total estimated cost to the end
      if (openList.size() > 1)
        currentNodeItr = std::min_element(
          openList.begin(), openList.end(), [bestNodeMap](const auto& next, const auto& best) {
            return std::get<2>(bestNodeMap.at(next)) < std::get<2>(bestNodeMap.at(best));
          });

      // Easier to deal directly with the pointer to the node
      const reco::ClusterHit3D* currentNode = *currentNodeItr;

      // Check to see if we have reached the goal and need to evaluate the path
      if (currentNode == goalNode) {
        // The path reconstruction will
        ReconstructBestPath(goalNode, bestNodeMap, pathNodeList, bestEdgeList);

        break;
      }

      // Otherwise need to keep evaluating
      else {
        openList.erase(currentNodeItr);
        currentNode->setStatusBit(reco::ClusterHit3D::PATHCHECKED);

        // Get tuple values for the current node
        const BestNodeTuple& currentNodeTuple = bestNodeMap.at(currentNode);
        float currentNodeScore = std::get<1>(currentNodeTuple);

        // Recover the edges associated to the current point
        const reco::EdgeList& curEdgeList = curEdgeMap[currentNode];

        for (const auto& curEdge : curEdgeList) {
          const reco::ClusterHit3D* candHit3D = std::get<1>(curEdge);

          if (candHit3D->getStatusBits() & reco::ClusterHit3D::PATHCHECKED) continue;

          float tentative_gScore = currentNodeScore + std::get<2>(curEdge);

          // Have we seen the candidate node before?
          BestNodeMap::iterator candNodeItr = bestNodeMap.find(candHit3D);

          if (candNodeItr == bestNodeMap.end()) { openList.push_back(candHit3D); }
          else if (tentative_gScore > std::get<1>(candNodeItr->second))
            continue;

          // Make a guess at score to get to target...
          float guessToTarget = DistanceBetweenNodes(candHit3D, goalNode) / 0.3;

          bestNodeMap[candHit3D] =
            BestNodeTuple(currentNode, tentative_gScore, tentative_gScore + guessToTarget);
        }
      }
    }
  }

  void MinSpanTreeAlg::ReconstructBestPath(const reco::ClusterHit3D* goalNode,
                                           BestNodeMap& bestNodeMap,
                                           reco::HitPairListPtr& pathNodeList,
                                           reco::EdgeList& bestEdgeList) const
  {
    while (std::get<0>(bestNodeMap.at(goalNode)) != goalNode) {
      const reco::ClusterHit3D* nextNode = std::get<0>(bestNodeMap[goalNode]);
      reco::EdgeTuple bestEdge =
        reco::EdgeTuple(goalNode, nextNode, DistanceBetweenNodes(goalNode, nextNode));

      pathNodeList.push_front(goalNode);
      bestEdgeList.push_front(bestEdge);

      goalNode = nextNode;
    }

    pathNodeList.push_front(goalNode);
  }

  void MinSpanTreeAlg::LeastCostPath(const reco::EdgeTuple& curEdge,
                                     const reco::ClusterHit3D* goalNode,
                                     reco::ClusterParameters& clusterParams,
                                     float& showMeTheMoney) const
  {
    // Recover the mapping between hits and edges
    reco::Hit3DToEdgeMap& curEdgeMap = clusterParams.getHit3DToEdgeMap();

    reco::Hit3DToEdgeMap::const_iterator edgeListItr = curEdgeMap.find(std::get<1>(curEdge));

    showMeTheMoney = std::numeric_limits<float>::max();

    if (edgeListItr != curEdgeMap.end() && !edgeListItr->second.empty()) {
      reco::HitPairListPtr& bestNodeList = clusterParams.getBestHitPairListPtr();
      reco::EdgeList& bestEdgeList = clusterParams.getBestEdgeList();

      for (const auto& edge : edgeListItr->second) {
        // skip the self reference
        if (std::get<1>(edge) == std::get<0>(curEdge)) continue;

        // Have we found the droid we are looking for?
        if (std::get<1>(edge) == goalNode) {
          bestNodeList.push_back(goalNode);
          bestEdgeList.push_back(edge);
          showMeTheMoney = std::get<2>(edge);
          break;
        }

        // Keep searching, it is out there somewhere...
        float currentCost(0.);

        LeastCostPath(edge, goalNode, clusterParams, currentCost);

        if (currentCost < std::numeric_limits<float>::max()) {
          showMeTheMoney = std::get<2>(edge) + currentCost;
          break;
        }
      }
    }

    if (showMeTheMoney < std::numeric_limits<float>::max()) {
      clusterParams.getBestHitPairListPtr().push_front(std::get<1>(curEdge));
      clusterParams.getBestEdgeList().push_front(curEdge);
    }
  }

  float MinSpanTreeAlg::DistanceBetweenNodes(const reco::ClusterHit3D* node1,
                                             const reco::ClusterHit3D* node2) const
  {
    const Eigen::Vector3f& node1Pos = node1->getPosition();
    const Eigen::Vector3f& node2Pos = node2->getPosition();
    float deltaNode[] = {
      node1Pos[0] - node2Pos[0], node1Pos[1] - node2Pos[1], node1Pos[2] - node2Pos[2]};

    // Standard euclidean distance
    return std::sqrt(deltaNode[0] * deltaNode[0] + deltaNode[1] * deltaNode[1] +
                     deltaNode[2] * deltaNode[2]);
  }

  reco::HitPairListPtr MinSpanTreeAlg::DepthFirstSearch(const reco::EdgeTuple& curEdge,
                                                        const reco::Hit3DToEdgeMap& hitToEdgeMap,
                                                        float& bestTreeQuality) const
  {
    reco::HitPairListPtr hitPairListPtr;
    float bestQuality(0.);
    float curEdgeWeight = std::max(0.3, std::get<2>(curEdge));
    float curEdgeProj(1. / curEdgeWeight);

    reco::Hit3DToEdgeMap::const_iterator edgeListItr = hitToEdgeMap.find(std::get<1>(curEdge));

    if (edgeListItr != hitToEdgeMap.end()) {
      // The input edge weight has quality factors applied, recalculate just the position difference
      const Eigen::Vector3f& firstHitPos = std::get<0>(curEdge)->getPosition();
      const Eigen::Vector3f& secondHitPos = std::get<1>(curEdge)->getPosition();
      float curEdgeVec[] = {secondHitPos[0] - firstHitPos[0],
                            secondHitPos[1] - firstHitPos[1],
                            secondHitPos[2] - firstHitPos[2]};
      float curEdgeMag = std::sqrt(curEdgeVec[0] * curEdgeVec[0] + curEdgeVec[1] * curEdgeVec[1] +
                                   curEdgeVec[2] * curEdgeVec[2]);

      curEdgeMag = std::max(float(0.1), curEdgeMag);

      for (const auto& edge : edgeListItr->second) {
        // skip the self reference
        if (std::get<1>(edge) == std::get<0>(curEdge)) continue;

        float quality(0.);

        reco::HitPairListPtr tempList = DepthFirstSearch(edge, hitToEdgeMap, quality);

        if (quality > bestQuality) {
          hitPairListPtr = tempList;
          bestQuality = quality;
          curEdgeProj = 1. / curEdgeMag;
        }
      }
    }

    hitPairListPtr.push_front(std::get<1>(curEdge));

    bestTreeQuality += bestQuality + curEdgeProj;

    return hitPairListPtr;
  }

  void MinSpanTreeAlg::PruneAmbiguousHits(reco::ClusterParameters& clusterParams,
                                          reco::Hit2DToClusterMap& hit2DToClusterMap) const
  {

    // Recover the HitPairListPtr from the input clusterParams (which will be the
    // only thing that has been provided)
    reco::HitPairListPtr& hitPairVector = clusterParams.getHitPairListPtr();

    size_t nStartedWith(hitPairVector.size());
    size_t nRejectedHits(0);

    reco::HitPairListPtr goodHits;

    // Loop through the hits and try to week out the clearly ambiguous ones
    for (const auto& hit3D : hitPairVector) {
      // Loop to try to remove ambiguous hits
      size_t n2DHitsIn3DHit(0);
      size_t nThisClusterOnly(0);
      size_t nOtherCluster(0);

      const std::set<const reco::ClusterHit3D*>* otherClusterHits = 0;

      for (const auto& hit2D : hit3D->getHits()) {
        if (!hit2D) continue;

        n2DHitsIn3DHit++;

        if (hit2DToClusterMap[hit2D].size() < 2)
          nThisClusterOnly = hit2DToClusterMap[hit2D][&clusterParams].size();
        else {
          for (const auto& clusterHitMap : hit2DToClusterMap[hit2D]) {
            if (clusterHitMap.first == &clusterParams) continue;

            if (clusterHitMap.second.size() > nOtherCluster) {
              nOtherCluster = clusterHitMap.second.size();
              otherClusterHits = &clusterHitMap.second;
            }
          }
        }
      }

      if (n2DHitsIn3DHit < 3 && nThisClusterOnly > 1 && nOtherCluster > 0) {
        bool skip3DHit(false);

        for (const auto& otherHit3D : *otherClusterHits) {
          size_t nOther2DHits(0);

          for (const auto& otherHit2D : otherHit3D->getHits()) {
            if (!otherHit2D) continue;

            nOther2DHits++;
          }

          if (nOther2DHits > 2) {
            skip3DHit = true;
            nRejectedHits++;
            break;
          }
        }

        if (skip3DHit) continue;
      }

      goodHits.emplace_back(hit3D);
    }

    std::cout << "###>> Input " << nStartedWith << " hits, rejected: " << nRejectedHits
              << std::endl;

    hitPairVector.resize(goodHits.size());
    std::copy(goodHits.begin(), goodHits.end(), hitPairVector.begin());
  }

  struct HitPairClusterOrder {
    bool operator()(const reco::ClusterParametersList::iterator& left,
                    const reco::ClusterParametersList::iterator& right)
    {
      // Watch out for the case where two clusters can have the same number of hits!
      return (*left).getHitPairListPtr().size() > (*right).getHitPairListPtr().size();
    }
  };

  class SetCheckHitOrder {
  public:
    SetCheckHitOrder(const std::vector<size_t>& plane) : m_plane(plane) {}

    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right) const
    {
      // Check if primary view's hit is on the same wire
      if (left->getWireIDs()[m_plane[0]] == right->getWireIDs()[m_plane[0]]) {
        // Same wire but not same hit, order by primary hit time
        if (left->getHits()[m_plane[0]] && right->getHits()[m_plane[0]] &&
            left->getHits()[m_plane[0]] != right->getHits()[m_plane[0]]) {
          return left->getHits()[m_plane[0]]->getHit()->PeakTime() <
                 right->getHits()[m_plane[0]]->getHit()->PeakTime();
        }

        // Primary view is same hit, look at next view's wire
        if (left->getWireIDs()[m_plane[1]] == right->getWireIDs()[m_plane[1]]) {
          // Same wire but not same hit, order by secondary hit time
          if (left->getHits()[m_plane[1]] && right->getHits()[m_plane[1]] &&
              left->getHits()[m_plane[1]] != right->getHits()[m_plane[1]]) {
            return left->getHits()[m_plane[1]]->getHit()->PeakTime() <
                   right->getHits()[m_plane[1]]->getHit()->PeakTime();
          }

          // All that is left is the final view... and this can't be the same hit... (else it is the same 3D hit)
          return left->getWireIDs()[m_plane[2]] < right->getWireIDs()[m_plane[2]];
        }

        return left->getWireIDs()[m_plane[1]] < right->getWireIDs()[m_plane[1]];
      }

      // Order by primary view's wire number
      return left->getWireIDs()[m_plane[0]] < right->getWireIDs()[m_plane[0]];
    }

  private:
    const std::vector<size_t>& m_plane;
  };

  void MinSpanTreeAlg::CheckHitSorting(reco::ClusterParameters& clusterParams) const
  {
    reco::HitPairListPtr& curCluster = clusterParams.getHitPairListPtr();

    // Trial A* here
    if (curCluster.size() > 2) {
      // Do a quick PCA to determine our parameter "alpha"
      reco::PrincipalComponents pca;
      m_pcaAlg.PCAAnalysis_3D(curCluster, pca);

      if (pca.getSvdOK()) {
        const Eigen::Vector3f& pcaAxis = pca.getEigenVectors().row(2);

        std::vector<size_t> closestPlane = {0, 0, 0};
        std::vector<float> bestAngle = {0., 0., 0.};

        for (size_t plane = 0; plane < 3; plane++) {
          const std::vector<float>& wireDir = m_wireDir[plane];

          float dotProd =
            std::fabs(pcaAxis[0] * wireDir[0] + pcaAxis[1] * wireDir[1] + pcaAxis[2] * wireDir[2]);

          if (dotProd > bestAngle[0]) {
            bestAngle[2] = bestAngle[1];
            closestPlane[2] = closestPlane[1];
            bestAngle[1] = bestAngle[0];
            closestPlane[1] = closestPlane[0];
            closestPlane[0] = plane;
            bestAngle[0] = dotProd;
          }
          else if (dotProd > bestAngle[1]) {
            bestAngle[2] = bestAngle[1];
            closestPlane[2] = closestPlane[1];
            closestPlane[1] = plane;
            bestAngle[1] = dotProd;
          }
          else {
            closestPlane[2] = plane;
            bestAngle[2] = dotProd;
          }
        }

        // Get a copy of our 3D hits
        reco::HitPairListPtr localHitList = curCluster;

        // Sort the hits
        localHitList.sort(SetCheckHitOrder(closestPlane));

        // Ok, let's print it all and take a look
        std::cout << "*****************************************************************************"
                     "***************"
                  << std::endl;
        std::cout << "**>>>>> longest axis: " << closestPlane[0] << ", best angle: " << bestAngle[0]
                  << std::endl;
        std::cout << "**>>>>> second  axis: " << closestPlane[1] << ", best angle: " << bestAngle[1]
                  << std::endl;
        std::cout << " " << std::endl;

        reco::HitPairListPtr::iterator firstHitItr = localHitList.begin();
        reco::HitPairListPtr::iterator lastHitItr = localHitList.begin();

        size_t bestPlane = closestPlane[0];

        reco::HitPairListPtr testList;

        while (firstHitItr != localHitList.end()) {
          const reco::ClusterHit3D* currentHit = *firstHitItr;

          // Search for the last matching best view hit
          while (lastHitItr != localHitList.end()) {
            // If a different wire on the best view then we're certainly done
            if (currentHit->getWireIDs()[bestPlane] != (*lastHitItr)->getWireIDs()[bestPlane])
              break;

            // More subtle test to see if same wire but different hit (being careful of case of no hit)
            if (currentHit->getHits()[bestPlane] && (*lastHitItr)->getHits()[bestPlane] &&
                currentHit->getHits()[bestPlane] != (*lastHitItr)->getHits()[bestPlane])
              break;

            // Yet event more subtle test...
            if ((!(currentHit->getHits()[bestPlane]) && (*lastHitItr)->getHits()[bestPlane]) ||
                (currentHit->getHits()[bestPlane] && !((*lastHitItr)->getHits()[bestPlane])))
              break;

            // Not there yet...
            lastHitItr++;
          }

          while (firstHitItr != lastHitItr) {
            currentHit = *++firstHitItr;
          }

          firstHitItr = lastHitItr;
        }
        curCluster = testList;
      }
    }
  }

  DEFINE_ART_CLASS_TOOL(MinSpanTreeAlg)
} // namespace lar_cluster3d

/**
 *  @file   MSTPathFinder_tool.cc
 *
 *  @brief  art Tool for comparing clusters and merging those that are consistent
 *
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "cetlib/cpu_timer.h"
#include "cetlib/search_path.h"

#include "larreco/RecoAlg/Cluster3DAlgs/ConvexHull/ConvexHull.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterModAlg.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterParamsBuilder.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/kdTree.h"

// Eigen
#include <Eigen/Dense>

// Root histograms
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

// std includes
#include <functional>
#include <iostream>
#include <memory>
#include <numeric> // std::accumulate
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

  class MSTPathFinder : virtual public IClusterModAlg {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit MSTPathFinder(const fhicl::ParameterSet&);

    /**
     *  @brief  Destructor
     */
    ~MSTPathFinder();

    void configure(fhicl::ParameterSet const& pset) override;

    /**
     *  @brief Interface for initializing histograms if they are desired
     *         Note that the idea is to put hisgtograms in a subfolder
     *
     *  @param TFileDirectory - the folder to store the hists in
     */
    void initializeHistograms(art::TFileDirectory&) override;

    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void ModifyClusters(reco::ClusterParametersList&) const override;

    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute() const override
    {
      return std::accumulate(fTimeVector.begin(), fTimeVector.end(), 0.);
    }

  private:
    /**
     *  @brief enumerate the possible values for time checking if monitoring timing
     */
    enum TimeValues {
      BUILDTHREEDHITS = 0,
      BUILDHITTOHITMAP = 1,
      RUNDBSCAN = 2,
      BUILDCLUSTERINFO = 3,
      PATHFINDING = 4,
      NUMTIMEVALUES
    };

    /**
     *  @brief Driver for Prim's algorithm
     */
    void RunPrimsAlgorithm(const reco::HitPairListPtr&,
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
               float alpha,
               kdTree::KdTreeNode&,
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

    /**
     *  @brief Add ability to build the convex hull (these needs to be split out! )
     */
    using MinMaxPoints = std::pair<reco::ProjectedPoint, reco::ProjectedPoint>;
    using MinMaxPointPair = std::pair<MinMaxPoints, MinMaxPoints>;

    void buildConvexHull(reco::ClusterParameters&, reco::HitPairListPtr&, int level = 0) const;

    float findConvexHullEndPoints(const reco::EdgeList&,
                                  const reco::ClusterHit3D*,
                                  const reco::ClusterHit3D*) const;

    /**
     *  @brief Data members to follow
     */
    bool fEnableMonitoring;     ///<
    size_t fMinTinyClusterSize; ///< Minimum size for a "tiny" cluster
    float fConvexHullKinkAngle; ///< Angle to declare a kink in convex hull calc
    float fConvexHullMinSep;    ///< Min hit separation to conisder in convex hull

    mutable std::vector<float> fTimeVector; ///<

    geo::Geometry const* fGeometry; //< pointer to the Geometry service

    PrincipalComponentsAlg fPCAAlg; // For running Principal Components Analysis
    kdTree fkdTree;                 // For the kdTree

    std::unique_ptr<lar_cluster3d::IClusterParametersBuilder>
      fClusterBuilder; ///<  Common cluster builder tool
  };

  MSTPathFinder::MSTPathFinder(fhicl::ParameterSet const& pset)
    : fPCAAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
    , fkdTree(pset.get<fhicl::ParameterSet>("kdTree"))
  {
    this->configure(pset);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  MSTPathFinder::~MSTPathFinder() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void MSTPathFinder::configure(fhicl::ParameterSet const& pset)
  {
    fEnableMonitoring = pset.get<bool>("EnableMonitoring", true);
    fMinTinyClusterSize = pset.get<size_t>("MinTinyClusterSize", 40);
    fConvexHullKinkAngle = pset.get<float>("ConvexHullKinkAgle", 0.95);
    fConvexHullMinSep = pset.get<float>("ConvexHullMinSep", 0.65);

    art::ServiceHandle<geo::Geometry const> geometry;

    fGeometry = &*geometry;

    fTimeVector.resize(NUMTIMEVALUES, 0.);

    fClusterBuilder = art::make_tool<lar_cluster3d::IClusterParametersBuilder>(
      pset.get<fhicl::ParameterSet>("ClusterParamsBuilder"));

    return;
  }

  void MSTPathFinder::initializeHistograms(art::TFileDirectory& histDir)
  {
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
    //    fFillHistograms = true;
    //
    //    std::string dirName = "ConvexHullPath";
    //
    //    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    //
    //    // Divide into two sets of hists... those for the overall cluster and
    //    // those for the subclusters
    //    fTopNum3DHits = dir.make<TH1F>("TopNum3DHits",  "Number 3D Hits",  200,    0.,   200.);

    return;
  }

  void MSTPathFinder::ModifyClusters(reco::ClusterParametersList& clusterParametersList) const
  {
    /**
     *  @brief Top level interface tool for performing deghosting and primary path finding
     *  using a minimum spanning tree approach. This is a shell tool, it actually uses the
     *  Minimum Spanning Tree clusering tool...
     */

    // Initial clustering is done, now trim the list and get output parameters
    cet::cpu_timer theClockBuildClusters;

    // Start clocks if requested
    if (fEnableMonitoring) theClockBuildClusters.start();

    // Ok, the idea here is to loop over the input clusters and the process one at a time and then use the MST algorithm
    // to deghost and try to find the best path.
    for (auto& clusterParams : clusterParametersList) {
      // It turns out that computing the convex hull surrounding the points in the 2D projection onto the
      // plane of largest spread in the PCA is a good way to break up the cluster... and we do it here since
      // we (currently) want this to be part of the standard output
      buildConvexHull(clusterParams, clusterParams.getHitPairListPtr());

      // Make sure our cluster has enough hits...
      if (clusterParams.getHitPairListPtr().size() > fMinTinyClusterSize) {
        // DBScan is driven of its "epsilon neighborhood". Computing adjacency within DBScan can be time
        // consuming so the idea is the prebuild the adjaceny map and then run DBScan.
        // The following call does this work
        kdTree::KdTreeNodeList kdTreeNodeContainer;
        kdTree::KdTreeNode topNode =
          fkdTree.BuildKdTree(clusterParams.getHitPairListPtr(), kdTreeNodeContainer);

        if (fEnableMonitoring) fTimeVector.at(BUILDHITTOHITMAP) = fkdTree.getTimeToExecute();

        // We are making subclusters
        reco::ClusterParametersList& daughterParametersList = clusterParams.daughterList();

        // Run DBScan to get candidate clusters
        RunPrimsAlgorithm(clusterParams.getHitPairListPtr(), topNode, daughterParametersList);

        // Initial clustering is done, now trim the list and get output parameters
        cet::cpu_timer theClockBuildClusters;

        // Start clocks if requested
        if (fEnableMonitoring) theClockBuildClusters.start();

        fClusterBuilder->BuildClusterInfo(daughterParametersList);

        if (fEnableMonitoring) {
          theClockBuildClusters.stop();

          fTimeVector[BUILDCLUSTERINFO] = theClockBuildClusters.accumulated_real_time();
        }

        // Test run the path finding algorithm
        for (auto& daughterParams : daughterParametersList)
          FindBestPathInCluster(daughterParams, topNode);
      }
    }

    if (fEnableMonitoring) {
      theClockBuildClusters.stop();

      fTimeVector[BUILDCLUSTERINFO] = theClockBuildClusters.accumulated_real_time();
    }

    mf::LogDebug("MSTPathFinder") << ">>>>> Cluster Path finding done" << std::endl;

    return;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
  void MSTPathFinder::RunPrimsAlgorithm(const reco::HitPairListPtr& hitPairList,
                                        kdTree::KdTreeNode& topNode,
                                        reco::ClusterParametersList& clusterParametersList) const
  {
    // If no hits then no work
    if (hitPairList.empty()) return;

    // Now proceed with building the clusters
    cet::cpu_timer theClockDBScan;

    // Start clocks if requested
    if (fEnableMonitoring) theClockDBScan.start();

    // Initialization
    size_t clusterIdx(0);

    // This will contain our list of edges
    reco::EdgeList curEdgeList;

    // Get the first point
    reco::HitPairListPtr::const_iterator freeHitItr = hitPairList.begin();
    const reco::ClusterHit3D* lastAddedHit = *freeHitItr++;

    lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);

    // Make a cluster...
    clusterParametersList.push_back(reco::ClusterParameters());

    // Get an iterator to the first cluster
    reco::ClusterParameters* curCluster = &clusterParametersList.back();

    // We use pointers here because the objects they point to will change in the loop below
    reco::Hit3DToEdgeMap* curEdgeMap = &curCluster->getHit3DToEdgeMap();
    reco::HitPairListPtr* curClusterHitList = &curCluster->getHitPairListPtr();

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
      curClusterHitList->push_back(lastAddedHit);

      // Set up to find the list of nearest neighbors to the last used hit...
      kdTree::CandPairList CandPairList;
      float bestDistance(1.5); //std::numeric_limits<float>::max());

      // And find them... result will be an unordered list of neigbors
      fkdTree.FindNearestNeighbors(lastAddedHit, topNode, CandPairList, bestDistance);

      // Copy edges to the current list (but only for hits not already in a cluster)
      //        for(auto& pair : CandPairList)
      //            if (!(pair.second->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)) curEdgeList.push_back(reco::EdgeTuple(lastAddedHit,pair.second,pair.first));
      for (auto& pair : CandPairList) {
        if (!(pair.second->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)) {
          double edgeWeight =
            pair.first * lastAddedHit->getHitChiSquare() * pair.second->getHitChiSquare();

          curEdgeList.push_back(reco::EdgeTuple(lastAddedHit, pair.second, edgeWeight));
        }
      }

      // If the edge list is empty then we have a complete cluster
      if (curEdgeList.empty()) {
        std::cout << "-----------------------------------------------------------------------------"
                     "------------"
                  << std::endl;
        std::cout << "**> Cluster idx: " << clusterIdx++ << " has " << curClusterHitList->size()
                  << " hits" << std::endl;

        // Look for the next "free" hit
        freeHitItr = std::find_if(freeHitItr, hitPairList.end(), [](const auto& hit) {
          return !(hit->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED);
        });

        // If at end of input list we are done with all hits
        if (freeHitItr == hitPairList.end()) break;

        std::cout << "##################################################################>"
                     "Processing another cluster"
                  << std::endl;

        // Otherwise, get a new cluster and set up
        clusterParametersList.push_back(reco::ClusterParameters());

        curCluster = &clusterParametersList.back();

        curEdgeMap = &curCluster->getHit3DToEdgeMap();
        curClusterHitList = &curCluster->getHitPairListPtr();
        lastAddedHit = *freeHitItr++;
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

    if (fEnableMonitoring) {
      theClockDBScan.stop();

      fTimeVector[RUNDBSCAN] = theClockDBScan.accumulated_real_time();
    }

    return;
  }

  void MSTPathFinder::FindBestPathInCluster(reco::ClusterParameters& curCluster) const
  {
    reco::HitPairListPtr longestCluster;
    float bestQuality(0.);
    float aveNumEdges(0.);
    size_t maxNumEdges(0);
    size_t nIsolatedHits(0);

    // Now proceed with building the clusters
    cet::cpu_timer theClockPathFinding;

    // Start clocks if requested
    if (fEnableMonitoring) theClockPathFinding.start();

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

    if (fEnableMonitoring) {
      theClockPathFinding.stop();

      fTimeVector[PATHFINDING] += theClockPathFinding.accumulated_real_time();
    }

    return;
  }

  void MSTPathFinder::FindBestPathInCluster(reco::ClusterParameters& clusterParams,
                                            kdTree::KdTreeNode& topNode) const
  {
    // Set up for timing the function
    cet::cpu_timer theClockPathFinding;

    // Start clocks if requested
    if (fEnableMonitoring) theClockPathFinding.start();

    // Trial A* here
    if (clusterParams.getHitPairListPtr().size() > 2) {
      // Do a quick PCA to determine our parameter "alpha"
      fPCAAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());

      // Recover the new fullPCA
      reco::PrincipalComponents& pca = clusterParams.getFullPCA();

      // The chances of a failure are remote, still we should check
      if (pca.getSvdOK()) {
        // Get references to what we need....
        reco::HitPairListPtr& curCluster = clusterParams.getHitPairListPtr();
        reco::Hit3DToEdgeMap& curEdgeMap = clusterParams.getHit3DToEdgeMap();
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

          // Call the AStar function to try to find the best path...
          //                AStar(startHit,stopHit,alpha,topNode,clusterParams);

          float cost(std::numeric_limits<float>::max());

          LeastCostPath(curEdgeMap[startHit].front(), stopHit, clusterParams, cost);

          clusterParams.getBestHitPairListPtr().push_front(startHit);

          std::cout << "**> Best path has " << clusterParams.getBestHitPairListPtr().size()
                    << " hits, " << clusterParams.getBestEdgeList().size() << " edges" << std::endl;
        }

        // Recalculate the PCA based on the hits comprisig the path
        fPCAAlg.PCAAnalysis_3D(clusterParams.getBestHitPairListPtr(), pca);

        // And now compute the convex hull
        buildConvexHull(clusterParams, clusterParams.getBestHitPairListPtr());
      }
      else {
        std::cout << "++++++>>> PCA failure! # hits: " << clusterParams.getHitPairListPtr().size()
                  << std::endl;
      }
    }

    if (fEnableMonitoring) {
      theClockPathFinding.stop();

      fTimeVector[PATHFINDING] += theClockPathFinding.accumulated_real_time();
    }

    return;
  }

  void MSTPathFinder::AStar(const reco::ClusterHit3D* startNode,
                            const reco::ClusterHit3D* goalNode,
                            float alpha,
                            kdTree::KdTreeNode& topNode,
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

    alpha = 1.; //std::max(0.5,alpha);

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

    return;
  }

  void MSTPathFinder::ReconstructBestPath(const reco::ClusterHit3D* goalNode,
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

    return;
  }

  void MSTPathFinder::LeastCostPath(const reco::EdgeTuple& curEdge,
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

    return;
  }

  float MSTPathFinder::DistanceBetweenNodes(const reco::ClusterHit3D* node1,
                                            const reco::ClusterHit3D* node2) const
  {
    const Eigen::Vector3f& node1Pos = node1->getPosition();
    const Eigen::Vector3f& node2Pos = node2->getPosition();
    float deltaNode[] = {
      node1Pos[0] - node2Pos[0], node1Pos[1] - node2Pos[1], node1Pos[2] - node2Pos[2]};

    // Standard euclidean distance
    return std::sqrt(deltaNode[0] * deltaNode[0] + deltaNode[1] * deltaNode[1] +
                     deltaNode[2] * deltaNode[2]);

    // Manhatten distance
    //return std::fabs(deltaNode[0]) + std::fabs(deltaNode[1]) + std::fabs(deltaNode[2]);
    /*
     // Chebyshev distance
     // We look for maximum distance by wires...

     // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
     int wireDeltas[] = {0,0,0};

     for(size_t idx = 0; idx < 3; idx++)
     wireDeltas[idx] = std::abs(int(node1->getWireIDs()[idx].Wire) - int(node2->getWireIDs()[idx].Wire));

     // put wire deltas in order...
     std::sort(wireDeltas, wireDeltas + 3);

     return std::sqrt(deltaNode[0]*deltaNode[0] + 0.09 * float(wireDeltas[2]*wireDeltas[2]));
     */
  }

  reco::HitPairListPtr MSTPathFinder::DepthFirstSearch(const reco::EdgeTuple& curEdge,
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

  void MSTPathFinder::PruneAmbiguousHits(reco::ClusterParameters& clusterParams,
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

      //        reco::ClusterParameters* otherCluster;
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
              //                        otherCluster     = clusterHitMap.first;
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

    return;
  }

  void MSTPathFinder::buildConvexHull(reco::ClusterParameters& clusterParameters,
                                      reco::HitPairListPtr& hitPairListPtr,
                                      int level) const
  {
    // set an indention string
    std::string minuses(level / 2, '-');
    std::string indent(level / 2, ' ');

    indent += minuses;

    // The plan is to build the enclosing 2D polygon around the points in the PCA plane of most spread for this cluster
    // To do so we need to start by building a list of 2D projections onto the plane of most spread...
    reco::PrincipalComponents& pca = clusterParameters.getFullPCA();

    // Recover the parameters from the Principal Components Analysis that we need to project and accumulate
    const Eigen::Vector3f& pcaCenter = pca.getAvePosition();
    reco::ConvexHull& convexHull = clusterParameters.getConvexHull();
    reco::ProjectedPointList& pointList = convexHull.getProjectedPointList();

    // Loop through hits and do projection to plane
    for (const auto& hit3D : hitPairListPtr) {
      Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
                                  hit3D->getPosition()[1] - pcaCenter(1),
                                  hit3D->getPosition()[2] - pcaCenter(2));
      Eigen::Vector3f pcaToHit = pca.getEigenVectors() * pcaToHitVec;

      pointList.emplace_back(pcaToHit(1), pcaToHit(2), hit3D);
    }

    // Sort the point vec by increasing x, then increase y
    pointList.sort([](const auto& left, const auto& right) {
      return (std::abs(std::get<0>(left) - std::get<0>(right)) >
              std::numeric_limits<float>::epsilon()) ?
               std::get<0>(left) < std::get<0>(right) :
               std::get<1>(left) < std::get<1>(right);
    });

    // containers for finding the "best" hull...
    std::vector<ConvexHull> convexHullVec;
    std::vector<reco::ProjectedPointList> rejectedListVec;
    bool increaseDepth(pointList.size() > 3);
    float lastArea(std::numeric_limits<float>::max());

    while (increaseDepth) {
      // Get another convexHull container
      convexHullVec.push_back(ConvexHull(pointList, fConvexHullKinkAngle, fConvexHullMinSep));
      rejectedListVec.push_back(reco::ProjectedPointList());

      const ConvexHull& convexHull = convexHullVec.back();
      reco::ProjectedPointList& rejectedList = rejectedListVec.back();
      const reco::ProjectedPointList& convexHullPoints = convexHull.getConvexHull();

      increaseDepth = false;

      if (convexHull.getConvexHullArea() > 0.) {
        if (convexHullVec.size() < 2 || convexHull.getConvexHullArea() < 0.8 * lastArea) {
          for (auto& point : convexHullPoints) {
            pointList.remove(point);
            rejectedList.emplace_back(point);
          }
          lastArea = convexHull.getConvexHullArea();
          //                increaseDepth = true;
        }
      }
    }

    // do we have a valid convexHull?
    while (!convexHullVec.empty() && convexHullVec.back().getConvexHullArea() < 0.5) {
      convexHullVec.pop_back();
      rejectedListVec.pop_back();
    }

    // If we found the convex hull then build edges around the region
    if (!convexHullVec.empty()) {
      size_t nRejectedTotal(0);
      reco::HitPairListPtr locHitPairListPtr = hitPairListPtr;

      for (const auto& rejectedList : rejectedListVec) {
        nRejectedTotal += rejectedList.size();

        for (const auto& rejectedPoint : rejectedList) {
          if (convexHullVec.back().findNearestDistance(rejectedPoint) > 0.5)
            locHitPairListPtr.remove(std::get<2>(rejectedPoint));
        }
      }

      // Now add "edges" to the cluster to describe the convex hull (for the display)
      reco::ProjectedPointList& convexHullPointList = convexHull.getConvexHullPointList();
      reco::Hit3DToEdgeMap& edgeMap = convexHull.getConvexHullEdgeMap();
      reco::EdgeList& edgeList = convexHull.getConvexHullEdgeList();

      reco::ProjectedPoint lastPoint = convexHullVec.back().getConvexHull().front();

      for (auto& curPoint : convexHullVec.back().getConvexHull()) {
        if (curPoint == lastPoint) continue;

        const reco::ClusterHit3D* lastPoint3D = std::get<2>(lastPoint);
        const reco::ClusterHit3D* curPoint3D = std::get<2>(curPoint);

        float distBetweenPoints = (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) *
                                    (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) +
                                  (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) *
                                    (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) +
                                  (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]) *
                                    (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]);

        distBetweenPoints = std::sqrt(distBetweenPoints);

        reco::EdgeTuple edge(lastPoint3D, curPoint3D, distBetweenPoints);

        convexHullPointList.push_back(curPoint);
        edgeMap[lastPoint3D].push_back(edge);
        edgeMap[curPoint3D].push_back(edge);
        edgeList.emplace_back(edge);

        lastPoint = curPoint;
      }

      // Store the "extreme" points
      const ConvexHull::PointList& extremePoints = convexHullVec.back().getExtremePoints();
      reco::ProjectedPointList& extremePointList = convexHull.getConvexHullExtremePoints();

      for (const auto& point : extremePoints)
        extremePointList.push_back(point);

      // Store the "kink" points
      const reco::ConvexHullKinkTupleList& kinkPoints = convexHullVec.back().getKinkPoints();
      reco::ConvexHullKinkTupleList& kinkPointList = convexHull.getConvexHullKinkPoints();

      for (const auto& kink : kinkPoints)
        kinkPointList.push_back(kink);
    }

    return;
  }

  float MSTPathFinder::findConvexHullEndPoints(const reco::EdgeList& convexHull,
                                               const reco::ClusterHit3D* first3D,
                                               const reco::ClusterHit3D* last3D) const
  {
    float largestDistance(0.);

    // Note that edges are vectors and that the convex hull edge list will be ordered
    // The idea is that the maximum distance from a given edge is to the edge just before the edge that "turns back" towards the current edge
    // meaning that the dot product of the two edges becomes negative.
    reco::EdgeList::const_iterator firstEdgeItr = convexHull.begin();

    while (firstEdgeItr != convexHull.end()) {
      reco::EdgeList::const_iterator nextEdgeItr = firstEdgeItr;

      //        Eigen::Vector2f firstEdgeVec(std::get<3>(*firstEdgeItr),std::get<);
      //        Eigen::Vector2f lastPrimaryVec(lastPCA.getEigenVectors()[0][0],lastPCA.getEigenVectors()[0][1],lastPCA.getEigenVectors()[0][2]);
      //        float           cosToLast = newPrimaryVec.dot(lastPrimaryVec);

      while (++nextEdgeItr != convexHull.end()) {}
    }

    return largestDistance;
  }

  DEFINE_ART_CLASS_TOOL(MSTPathFinder)
} // namespace lar_cluster3d

/**
 *  @file   MinSpanTreeAlg.cxx
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/kdTree.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterParamsBuilder.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <Eigen/Dense>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {
    
class MinSpanTreeAlg : virtual public IClusterAlg
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit MinSpanTreeAlg(const fhicl::ParameterSet&);
    
    /**
     *  @brief  Destructor
     */
    ~MinSpanTreeAlg();
    
    void configure(fhicl::ParameterSet const &pset) override;
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param hitPairClusterMap     A map of hits that have been clustered
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void Cluster3DHits(reco::HitPairList&           hitPairList,
                       reco::ClusterParametersList& clusterParametersList) const override;
    
    void Cluster3DHits(reco::HitPairListPtr&        hitPairList,
                       reco::ClusterParametersList& clusterParametersList) const override {return;}

    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute(TimeValues index) const override {return m_timeVector.at(index);}
    
private:
    
    /**
     *  @brief Driver for Prim's algorithm
     */
    void RunPrimsAlgorithm(reco::HitPairList&, kdTree::KdTreeNode&, reco::ClusterParametersList&) const;
    
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
    reco::HitPairListPtr DepthFirstSearch(const reco::EdgeTuple&, const reco::Hit3DToEdgeMap&, float&) const;
    
    /**
     *  @brief Alternative version of FindBestPathInCluster utilizing an A* algorithm
     */
    void FindBestPathInCluster(reco::ClusterParameters&, kdTree::KdTreeNode&) const;
    
    /**
     *  @brief Algorithm to find shortest path between two 3D hits
     */
    void AStar(const reco::ClusterHit3D*, const reco::ClusterHit3D*, float alpha, kdTree::KdTreeNode&, reco::HitPairListPtr&, reco::EdgeList&) const;
    
    using BestNodeTuple = std::tuple<const reco::ClusterHit3D*,float,float>;
    using BestNodeMap   = std::unordered_map<const reco::ClusterHit3D*,BestNodeTuple>;
    
    void ReconstructBestPath(const reco::ClusterHit3D*, BestNodeMap&, reco::HitPairListPtr&, reco::EdgeList&) const;
    
    float DistanceBetweenNodes(const reco::ClusterHit3D*,const reco::ClusterHit3D*) const;
    
    /**
     *  @brief Find the lowest cost path between two nodes using MST edges
     */
    void LeastCostPath(const reco::EdgeTuple&,
                       const reco::ClusterHit3D*,
                       const reco::Hit3DToEdgeMap&,
                       reco::HitPairListPtr&,
                       reco::EdgeList&,
                       float&) const;
    
    void CheckHitSorting(reco::ClusterParameters& clusterParams) const;
    
    /**
     *  @brief define data structure for keeping track of channel status
     */
    using ChannelStatusVec        = std::vector<size_t>;
    using ChannelStatusByPlaneVec = std::vector<ChannelStatusVec>;
    
    /**
     *  @brief Data members to follow
     */
    bool                                                      m_enableMonitoring;      ///<
    mutable std::vector<float>                                m_timeVector;            ///<
    std::vector<std::vector<float>>                           m_wireDir;               ///<
                         
    geo::Geometry const*                                      m_geometry;              //< pointer to the Geometry service
                         
    PrincipalComponentsAlg                                    m_pcaAlg;                // For running Principal Components Analysis
    kdTree                                                    m_kdTree;                // For the kdTree

    std::unique_ptr<lar_cluster3d::IClusterParametersBuilder> m_clusterBuilder;        ///<  Common cluster builder tool
};

MinSpanTreeAlg::MinSpanTreeAlg(fhicl::ParameterSet const &pset) :
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg")),
    m_kdTree(pset.get<fhicl::ParameterSet>("kdTree"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

MinSpanTreeAlg::~MinSpanTreeAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MinSpanTreeAlg::configure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring         = pset.get<bool>  ("EnableMonitoring",  true  );
    
    art::ServiceHandle<geo::Geometry const> geometry;
    
    m_geometry = &*geometry;
    
    m_timeVector.resize(NUMTIMEVALUES, 0.);
    
    // Determine the unit directon and normal vectors to the wires
    m_wireDir.resize(3);
    
    raw::ChannelID_t uChannel(0);
    std::vector<geo::WireID> uWireID  = m_geometry->ChannelToWire(uChannel);
    const geo::WireGeo*      uWireGeo = m_geometry->WirePtr(uWireID[0]);
    
    TVector3 uWireDir = uWireGeo->Direction();
    
    m_wireDir[0].resize(3);
    m_wireDir[0][0] =  uWireDir[0];
    m_wireDir[0][1] = -uWireDir[2];
    m_wireDir[0][2] =  uWireDir[1];
    
    raw::ChannelID_t vChannel(2400);
    std::vector<geo::WireID> vWireID = m_geometry->ChannelToWire(vChannel);
    const geo::WireGeo* vWireGeo = m_geometry->WirePtr(vWireID[0]);
    
    TVector3 vWireDir = vWireGeo->Direction();
    
    m_wireDir[1].resize(3);
    m_wireDir[1][0] =  vWireDir[0];
    m_wireDir[1][1] = -vWireDir[2];
    m_wireDir[1][2] =  vWireDir[1];
    
    m_wireDir[2].resize(3);
    m_wireDir[2][0] = 0.;
    m_wireDir[2][1] = 0.;
    m_wireDir[2][2] = 1.;

    m_clusterBuilder = art::make_tool<lar_cluster3d::IClusterParametersBuilder>(pset.get<fhicl::ParameterSet>("ClusterParamsBuilder"));
    
    return;
}
    
void MinSpanTreeAlg::Cluster3DHits(reco::HitPairList&           hitPairList,
                                   reco::ClusterParametersList& clusterParametersList) const
{
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */
    
    // Zero the time vector
    if (m_enableMonitoring) std::fill(m_timeVector.begin(),m_timeVector.end(),0.);
    
    // DBScan is driven of its "epsilon neighborhood". Computing adjacency within DBScan can be time
    // consuming so the idea is the prebuild the adjaceny map and then run DBScan.
    // The following call does this work
    kdTree::KdTreeNodeList kdTreeNodeContainer;
    kdTree::KdTreeNode     topNode = m_kdTree.BuildKdTree(hitPairList, kdTreeNodeContainer);
    
    if (m_enableMonitoring) m_timeVector.at(BUILDHITTOHITMAP) = m_kdTree.getTimeToExecute();
    
    // Run DBScan to get candidate clusters
    RunPrimsAlgorithm(hitPairList, topNode, clusterParametersList);
    
    // Initial clustering is done, now trim the list and get output parameters
    cet::cpu_timer theClockBuildClusters;
    
    // Start clocks if requested
    if (m_enableMonitoring) theClockBuildClusters.start();
    
    m_clusterBuilder->BuildClusterInfo(clusterParametersList);
    
    if (m_enableMonitoring)
    {
        theClockBuildClusters.stop();
        
        m_timeVector[BUILDCLUSTERINFO] = theClockBuildClusters.accumulated_real_time();
    }
    
    // Test run the path finding algorithm
    for (auto& clusterParams : clusterParametersList) FindBestPathInCluster(clusterParams, topNode);
    
    mf::LogDebug("Cluster3D") << ">>>>> Cluster3DHits done, found " << clusterParametersList.size() << " clusters" << std::endl;
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
void MinSpanTreeAlg::RunPrimsAlgorithm(reco::HitPairList&           hitPairList,
                                       kdTree::KdTreeNode&          topNode,
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
    reco::HitPairList::iterator freeHitItr   = hitPairList.begin();
    const reco::ClusterHit3D*   lastAddedHit = &(*freeHitItr++);
    
    lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);
    
    // Make a cluster...
    clusterParametersList.push_back(reco::ClusterParameters());
    
    // Get an iterator to the first cluster
    reco::ClusterParametersList::iterator curClusterItr = --clusterParametersList.end();
    
    // We use pointers here because the objects they point to will change in the loop below
    reco::Hit3DToEdgeMap* curEdgeMap = &(*curClusterItr).getHit3DToEdgeMap();
    reco::HitPairListPtr* curCluster = &(*curClusterItr).getHitPairListPtr();
    
    // Loop until all hits have been associated to a cluster
    while(1)
    {
        // and the 3D hit status bits
        lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);
        
        // Purge the current list to get rid of edges which point to hits already in the cluster
        for(reco::EdgeList::iterator curEdgeItr = curEdgeList.begin(); curEdgeItr != curEdgeList.end();)
        {
            if (std::get<1>(*curEdgeItr)->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)
                curEdgeItr = curEdgeList.erase(curEdgeItr);
            else curEdgeItr++;
        }
        
        // Add the lastUsedHit to the current cluster
        curCluster->push_back(lastAddedHit);
        
        // Set up to find the list of nearest neighbors to the last used hit...
        kdTree::CandPairList CandPairList;
        float                bestDistance(1.5); //std::numeric_limits<float>::max());

        // And find them... result will be an unordered list of neigbors
        m_kdTree.FindNearestNeighbors(lastAddedHit, topNode, CandPairList, bestDistance);
        
        // Copy edges to the current list (but only for hits not already in a cluster)
        for(auto& pair : CandPairList)
            if (!(pair.second->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)) curEdgeList.push_back(reco::EdgeTuple(lastAddedHit,pair.second,pair.first));
        
        // If the edge list is empty then we have a complete cluster
        if (curEdgeList.empty())
        {
            std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
            std::cout << "**> Cluster idx: " << clusterIdx++ << " has " << curCluster->size() << " hits" << std::endl;

            // Look for the next "free" hit
            freeHitItr = std::find_if(freeHitItr,hitPairList.end(),[](const auto& hit){return !(hit.getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED);});
            
            // If at end of input list we are done with all hits
            if (freeHitItr == hitPairList.end()) break;
            
            std::cout << "##################################################################>Processing another cluster" << std::endl;
            
            // Otherwise, get a new cluster and set up
            clusterParametersList.push_back(reco::ClusterParameters());
            
            curClusterItr = --clusterParametersList.end();
            
            curEdgeMap   = &(*curClusterItr).getHit3DToEdgeMap();
            curCluster   = &(*curClusterItr).getHitPairListPtr();
            lastAddedHit = &(*freeHitItr++);
        }
        // Otherwise we are still processing the current cluster
        else
        {
            // Sort the list of edges by distance
            curEdgeList.sort([](const auto& left,const auto& right){return std::get<2>(left) < std::get<2>(right);});
            
            // Populate the map with the edges...
            reco::EdgeTuple& curEdge = curEdgeList.front();
            
            (*curEdgeMap)[std::get<0>(curEdge)].push_back(curEdge);
            (*curEdgeMap)[std::get<1>(curEdge)].push_back(reco::EdgeTuple(std::get<1>(curEdge),std::get<0>(curEdge),std::get<2>(curEdge)));
            
            // Update the last hit to be added to the collection
            lastAddedHit = std::get<1>(curEdge);
        }
    }

    if (m_enableMonitoring)
    {
        theClockDBScan.stop();
        
        m_timeVector[RUNDBSCAN] = theClockDBScan.accumulated_real_time();
    }
    
    return;
}

void MinSpanTreeAlg::FindBestPathInCluster(reco::ClusterParameters& curCluster) const
{
    reco::HitPairListPtr longestCluster;
    float                bestQuality(0.);
    float                aveNumEdges(0.);
    size_t               maxNumEdges(0);
    size_t               nIsolatedHits(0);
    
    // Now proceed with building the clusters
    cet::cpu_timer theClockPathFinding;
    
    // Start clocks if requested
    if (m_enableMonitoring) theClockPathFinding.start();
    
    reco::HitPairListPtr& hitPairList  = curCluster.getHitPairListPtr();
    reco::Hit3DToEdgeMap& curEdgeMap   = curCluster.getHit3DToEdgeMap();
    reco::EdgeList&       bestEdgeList = curCluster.getBestEdgeList();
    
    // Do some spelunking...
    for(const auto& hit : hitPairList)
    {
        if (!curEdgeMap[hit].empty() && curEdgeMap[hit].size() == 1)
        {
            float quality(0.);
            
            reco::HitPairListPtr tempList = DepthFirstSearch(curEdgeMap[hit].front(), curEdgeMap, quality);
            
            tempList.push_front(std::get<0>(curEdgeMap[hit].front()));
            
            if (quality > bestQuality)
            {
                longestCluster = tempList;
                bestQuality    = quality;
            }
            
            nIsolatedHits++;
        }
        
        aveNumEdges += float(curEdgeMap[hit].size());
        maxNumEdges  = std::max(maxNumEdges,curEdgeMap[hit].size());
    }
    
    aveNumEdges /= float(hitPairList.size());
    std::cout << "----> # isolated hits: " << nIsolatedHits << ", longest branch: " << longestCluster.size() << ", cluster size: " << hitPairList.size() << ", ave # edges: " << aveNumEdges << ", max: " << maxNumEdges << std::endl;
    
    if (!longestCluster.empty())
    {
        hitPairList = longestCluster;
        for(const auto& hit : hitPairList)
        {
            for(const auto& edge : curEdgeMap[hit]) bestEdgeList.emplace_back(edge);
        }
        
        std::cout << "        ====> new cluster size: " << hitPairList.size() << std::endl;
    }
    
    if (m_enableMonitoring)
    {
        theClockPathFinding.stop();
        
        m_timeVector[PATHFINDING] += theClockPathFinding.accumulated_real_time();
    }
    
    return;
}
    
void MinSpanTreeAlg::FindBestPathInCluster(reco::ClusterParameters& clusterParams, kdTree::KdTreeNode& topNode) const
{
    // Set up for timing the function
    cet::cpu_timer theClockPathFinding;
    
    // Start clocks if requested
    if (m_enableMonitoring) theClockPathFinding.start();
    
    reco::HitPairListPtr& curCluster = clusterParams.getHitPairListPtr();
    reco::Hit3DToEdgeMap& curEdgeMap = clusterParams.getHit3DToEdgeMap();
    
    // Trial A* here
    if (curCluster.size() > 2)
    {
        // Do a quick PCA to determine our parameter "alpha"
        reco::PrincipalComponents pca;
        m_pcaAlg.PCAAnalysis_3D(curCluster, pca);
        
        if (pca.getSvdOK())
        {
            const Eigen::Vector3f&    pcaAxis   = pca.getEigenVectors().row(2);
            float                     pcaLen    = 1.5*sqrt(pca.getEigenValues()[2]);
            float                     pcaWidth  = 1.5*sqrt(pca.getEigenValues()[1]);
            float                     pcaHeight = 1.5*sqrt(pca.getEigenValues()[0]);
            const Eigen::Vector3f&    pcaPos    = pca.getAvePosition();
            float                     alpha     = std::min(float(1.),std::max(float(0.001),pcaWidth/pcaLen));
        
            // The first task is to find the list of hits which are "isolated"
            reco::HitPairListPtr isolatedHitList;
            for(const auto& hit : curCluster)
                if (!curEdgeMap[hit].empty() && curEdgeMap[hit].size() == 1)
                    isolatedHitList.emplace_back(std::get<0>(curEdgeMap[hit].front()));
        
            std::cout << "************* Finding best path with A* in cluster *****************" << std::endl;
            std::cout << "**> There are " << curCluster.size() << " hits, " << isolatedHitList.size() << " isolated hits, the alpha parameter is " << alpha << std::endl;
            std::cout << "**> PCA len: " << pcaLen << ", wid: " << pcaWidth << ", height: " << pcaHeight << ", ratio: " << pcaHeight/pcaWidth << std::endl;
        
            // Goal is to now find separated pairs of isolated hits
            reco::EdgeList edgeList;
            reco::HitPairListPtr::iterator firstItr = isolatedHitList.begin();
            while(firstItr != isolatedHitList.end())
            {
                const reco::ClusterHit3D* firstHit = *firstItr++;
            
                const Eigen::Vector3f& firstPos      = firstHit->getPosition();
                float                  delta1stPca[] = {firstPos[0]-pcaPos[0],firstPos[1]-pcaPos[1],firstPos[2]-pcaPos[2]};
                float                  firstPcaProj  = std::fabs(delta1stPca[0]*pcaAxis[0] + delta1stPca[1]*pcaAxis[1] + delta1stPca[2]*pcaAxis[2]);
            
                if (firstPcaProj < 0.75 * pcaLen) continue;
            
                for(reco::HitPairListPtr::iterator secondItr = firstItr; secondItr != isolatedHitList.end(); secondItr++)
                {
                    const Eigen::Vector3f& secondPos     = (*secondItr)->getPosition();
                    float                  delta2ndPca[] = {secondPos[0]-pcaPos[0],secondPos[1]-pcaPos[1],secondPos[2]-pcaPos[2]};
                    float                  secondPcaProj = std::fabs(delta2ndPca[0]*pcaAxis[0] + delta2ndPca[1]*pcaAxis[1] + delta2ndPca[2]*pcaAxis[2]);
                
                    if (secondPcaProj < 0.75 * pcaLen) continue;
                
                    float deltaPos[] = {secondPos[0]-firstPos[0],secondPos[1]-firstPos[1],secondPos[2]-firstPos[2]};
                    float projection = std::fabs(deltaPos[0]*pcaAxis[0]+deltaPos[1]*pcaAxis[1]+deltaPos[2]*pcaAxis[2]);
                
                    edgeList.emplace_back(firstHit,*secondItr,projection);
                    //edgeList.emplace_back(reco::EdgeTuple(firstHit,*secondItr,DistanceBetweenNodes(firstHit,*secondItr)));
                }
            }
            
            if (edgeList.empty())
            {
                if (isolatedHitList.size() > 20)
                {
                    std::cout << "!!!! What happened???? " << std::endl;
                }
            }
            
            if (!edgeList.empty())
            {
                edgeList.sort([](const auto& left,const auto& right){return std::get<2>(left) > std::get<2>(right);});
        
                // Ok, trial the algorithm by simply looking for the path between the hits at the front of the list
                reco::EdgeTuple& bestEdge = edgeList.front();
        
                std::cout << "**> Sorted " << edgeList.size() << " edges, longest distance: " << DistanceBetweenNodes(std::get<0>(bestEdge),std::get<1>(bestEdge)) << std::endl;
        
                reco::HitPairListPtr& bestHitPairListPtr = clusterParams.getBestHitPairListPtr();
                reco::EdgeList&       bestEdgeList       = clusterParams.getBestEdgeList();
        
                AStar(std::get<0>(bestEdge),std::get<1>(bestEdge),alpha,topNode,bestHitPairListPtr,bestEdgeList);

//                float cost(std::numeric_limits<float>::max());
                
//                LeastCostPath(curEdgeMap[std::get<0>(bestEdge)].front(),std::get<1>(bestEdge),curEdgeMap,bestPathHitList,bestEdgeList,cost);
                
//                bestPathHitList.push_front(std::get<0>(bestEdge));
        
                std::cout << "**> Best path has " << bestHitPairListPtr.size() << " hits, " << bestEdgeList.size() << " edges" << std::endl;
            }
        }
        else
        {
            std::cout << "++++++>>> PCA failure! # hits: " << curCluster.size() << std::endl;
        }
    }
    
    if (m_enableMonitoring)
    {
        theClockPathFinding.stop();
        
        m_timeVector[PATHFINDING] += theClockPathFinding.accumulated_real_time();
    }
    
    return;
}
    
void MinSpanTreeAlg::AStar(const reco::ClusterHit3D* startNode,
                           const reco::ClusterHit3D* goalNode,
                           float                     alpha,
                           kdTree::KdTreeNode&       topNode,
                           reco::HitPairListPtr&     pathNodeList,
                           reco::EdgeList&           bestEdgeList) const
{
    // Find the shortest path from start to goal using an A* algorithm
    // Keep track of the nodes which have already been evaluated
    reco::HitPairListPtr closedList;
    
    // Keep track of the nodes that have been "discovered" but yet to be evaluated
    reco::HitPairListPtr openList = {startNode};
    
    // Create a map which, for each node, will tell us the node it can be most efficiencly reached from.
    BestNodeMap bestNodeMap;
    
    bestNodeMap[startNode] = BestNodeTuple(startNode,0.,DistanceBetweenNodes(startNode,goalNode));

    alpha = 1.; //std::max(0.5,alpha);
    
    while(!openList.empty())
    {
        // The list is not empty so by def we will return something
        reco::HitPairListPtr::iterator currentNodeItr = openList.begin();
        
        // If the list contains more than one element then we need to find the one with the smallest total estimated cost to the end
        if (openList.size() > 1)
            currentNodeItr = std::min_element(openList.begin(),openList.end(),[bestNodeMap](const auto& next, const auto& best){return std::get<2>(bestNodeMap.at(next)) < std::get<2>(bestNodeMap.at(best));});

        // Easier to deal directly with the pointer to the node
        const reco::ClusterHit3D* currentNode = *currentNodeItr;
        
        // Check to see if we have reached the goal and need to evaluate the path
        if (currentNode == goalNode)
        {
            // The path reconstruction will
            ReconstructBestPath(goalNode, bestNodeMap, pathNodeList, bestEdgeList);
            
//            std::cout << "**> Reconstructed best path... ended with " << openList.size() << " hits in openList" << std::endl;
            
            break;
        }
        
        // Otherwise need to keep evaluating
        else
        {
            openList.erase(currentNodeItr);
//            closedList.push_front(currentNode);
            currentNode->setStatusBit(reco::ClusterHit3D::PATHCHECKED);
            
            // Set up to find the list of nearest neighbors to the last used hit...
            kdTree::CandPairList CandPairList;
            float                bestDistance(std::numeric_limits<float>::max());
            
            // And find them... result will be an unordered list of neigbors
            m_kdTree.FindNearestNeighbors(currentNode, topNode, CandPairList, bestDistance);
            
//            std::cout << "**> found " << CandPairList.size() << " nearest neigbhors, bestDistance: " << bestDistance;
//            size_t nAdded(0);
            
            // Get tuple values for the current node
            const BestNodeTuple& currentNodeTuple = bestNodeMap.at(currentNode);
            float                currentNodeScore = std::get<1>(currentNodeTuple);
            
            for(auto& candPair : CandPairList)
            {
                // Ignore those nodes we're already aware of
                //if (std::find(closedList.begin(),closedList.end(),candPair.second) != closedList.end()) continue;
                if (candPair.second->getStatusBits() & reco::ClusterHit3D::PATHCHECKED) continue;
                
                float tentative_gScore = currentNodeScore + candPair.first;
                
                // Have we seen the candidate node before?
                BestNodeMap::iterator candNodeItr = bestNodeMap.find(candPair.second);
                
                if (candNodeItr == bestNodeMap.end())
                {
                    openList.push_back(candPair.second);
//                    nAdded++;
                }
                else if (tentative_gScore > std::get<1>(candNodeItr->second)) continue;

                // Experiment with modification to cost estimate
                const Eigen::Vector3f& currentNodePos  = currentNode->getPosition();
                const Eigen::Vector3f& nextNodePos     = candPair.second->getPosition();
                float                  curNextDelta[]  = {nextNodePos[0]-currentNodePos[0], nextNodePos[1]-currentNodePos[1], nextNodePos[2]-currentNodePos[2]};
                
                const Eigen::Vector3f& goalNodePos     = goalNode->getPosition();
                float                  goalNextDelta[] = {goalNodePos[0]-nextNodePos[0], goalNodePos[1]-nextNodePos[1], goalNodePos[2]-nextNodePos[2]};
                
                float                  curNextMag      = std::sqrt(curNextDelta[0]*curNextDelta[0]   + curNextDelta[1]*curNextDelta[1]   + curNextDelta[2]*curNextDelta[2]);
                float                  goalNextMag     = std::sqrt(goalNextDelta[0]*goalNextDelta[0] + goalNextDelta[1]*goalNextDelta[1] + goalNextDelta[2]*goalNextDelta[2]);
                
                float                  cosTheta        = (curNextDelta[0]*goalNextDelta[0] + curNextDelta[1]*goalNextDelta[1] + curNextDelta[2]*goalNextDelta[2]);
                
                if (cosTheta > 0. || cosTheta < 0.) cosTheta /= (curNextMag * goalNextMag);
                
//                alpha = candPair.second->getMinOverlapFraction();
                cosTheta = 1.;
                
                float hWeight = alpha*goalNextMag/std::max(0.01,0.5*(1.+cosTheta));

                // update our records
                bestNodeMap[candPair.second] = BestNodeTuple(currentNode,tentative_gScore, tentative_gScore + hWeight);
                //bestNodeMap[candPair.second] = BestNodeTuple(currentNode, tentative_gScore, tentative_gScore + alpha*DistanceBetweenNodes(candPair.second,goalNode));
            }
            
//            std::cout << ", added: " << nAdded << ", openList size: " << openList.size() << std::endl;
        }
    }
    
    return;
}

void MinSpanTreeAlg::ReconstructBestPath(const reco::ClusterHit3D* goalNode,
                                         BestNodeMap&              bestNodeMap,
                                         reco::HitPairListPtr&     pathNodeList,
                                         reco::EdgeList&           bestEdgeList) const
{
    while(std::get<0>(bestNodeMap.at(goalNode)) != goalNode)
    {
        const reco::ClusterHit3D* nextNode = std::get<0>(bestNodeMap[goalNode]);
        reco::EdgeTuple           bestEdge = reco::EdgeTuple(goalNode,nextNode,DistanceBetweenNodes(goalNode,nextNode));
        
        pathNodeList.push_front(goalNode);
        bestEdgeList.push_front(bestEdge);
        
        goalNode = nextNode;
    }
    
    pathNodeList.push_front(goalNode);
    
    return;
}
    
void MinSpanTreeAlg::LeastCostPath(const reco::EdgeTuple&      curEdge,
                                   const reco::ClusterHit3D*   goalNode,
                                   const reco::Hit3DToEdgeMap& hitToEdgeMap,
                                   reco::HitPairListPtr&       bestNodeList,
                                   reco::EdgeList&             bestEdgeList,
                                   float&                      showMeTheMoney) const
{
    reco::Hit3DToEdgeMap::const_iterator edgeListItr = hitToEdgeMap.find(std::get<1>(curEdge));
    
    showMeTheMoney = std::numeric_limits<float>::max();
    
    if (edgeListItr != hitToEdgeMap.end())
    {
        for(const auto& edge : edgeListItr->second)
        {
            // skip the self reference
            if (std::get<1>(edge) == std::get<0>(curEdge)) continue;
            
            // Have we found the droid we are looking for?
            if (std::get<1>(edge) == goalNode)
            {
                bestNodeList.push_back(goalNode);
                bestEdgeList.push_back(edge);
                showMeTheMoney = std::get<2>(edge);
                break;
            }
            
            // Keep searching, it is out there somewhere...
            float currentCost(0.);
            
            LeastCostPath(edge,goalNode,hitToEdgeMap,bestNodeList,bestEdgeList,currentCost);
            
            if (currentCost < std::numeric_limits<float>::max())
            {
                showMeTheMoney = std::get<2>(edge) + currentCost;
                break;
            }
        }
    }
    
    if (showMeTheMoney < std::numeric_limits<float>::max())
    {
        bestNodeList.push_front(std::get<1>(curEdge));
        bestEdgeList.push_front(curEdge);
    }
    
    return;
}
    
float MinSpanTreeAlg::DistanceBetweenNodes(const reco::ClusterHit3D* node1,const reco::ClusterHit3D* node2) const
{
    const Eigen::Vector3f& node1Pos    = node1->getPosition();
    const Eigen::Vector3f& node2Pos    = node2->getPosition();
    float                  deltaNode[] = {node1Pos[0]-node2Pos[0], node1Pos[1]-node2Pos[1], node1Pos[2]-node2Pos[2]};
    
    // Standard euclidean distance
    return std::sqrt(deltaNode[0]*deltaNode[0]+deltaNode[1]*deltaNode[1]+deltaNode[2]*deltaNode[2]);
    
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

reco::HitPairListPtr MinSpanTreeAlg::DepthFirstSearch(const reco::EdgeTuple&      curEdge,
                                                      const reco::Hit3DToEdgeMap& hitToEdgeMap,
                                                      float&                      bestTreeQuality) const
{
    reco::HitPairListPtr hitPairListPtr;
    float               bestQuality(0.);
    float               curEdgeWeight = std::max(0.3,std::get<2>(curEdge));
    float               curEdgeProj(1./curEdgeWeight);
    
    reco::Hit3DToEdgeMap::const_iterator edgeListItr = hitToEdgeMap.find(std::get<1>(curEdge));
    
    if (edgeListItr != hitToEdgeMap.end())
    {
        // The input edge weight has quality factors applied, recalculate just the position difference
        const Eigen::Vector3f& firstHitPos  = std::get<0>(curEdge)->getPosition();
        const Eigen::Vector3f& secondHitPos = std::get<1>(curEdge)->getPosition();
        float                  curEdgeVec[] = {secondHitPos[0]-firstHitPos[0],secondHitPos[1]-firstHitPos[1],secondHitPos[2]-firstHitPos[2]};
        float                  curEdgeMag   = std::sqrt(curEdgeVec[0]*curEdgeVec[0]+curEdgeVec[1]*curEdgeVec[1]+curEdgeVec[2]*curEdgeVec[2]);
        
        curEdgeMag = std::max(float(0.1),curEdgeMag);
        
        for(const auto& edge : edgeListItr->second)
        {
            // skip the self reference
            if (std::get<1>(edge) == std::get<0>(curEdge)) continue;
            
            float quality(0.);
            
            reco::HitPairListPtr tempList = DepthFirstSearch(edge,hitToEdgeMap,quality);
            
            if (quality > bestQuality)
            {
                hitPairListPtr = tempList;
                bestQuality    = quality;
                curEdgeProj    = 1./curEdgeMag;
            }
        }
    }
    
    hitPairListPtr.push_front(std::get<1>(curEdge));
    
    bestTreeQuality += bestQuality + curEdgeProj;
    
    return hitPairListPtr;
}
    
void MinSpanTreeAlg::PruneAmbiguousHits(reco::ClusterParameters& clusterParams, reco::Hit2DToClusterMap& hit2DToClusterMap) const
{
    
    // Recover the HitPairListPtr from the input clusterParams (which will be the
    // only thing that has been provided)
    reco::HitPairListPtr& hitPairVector = clusterParams.getHitPairListPtr();
    
    size_t nStartedWith(hitPairVector.size());
    size_t nRejectedHits(0);
    
    reco::HitPairListPtr goodHits;
    
    // Loop through the hits and try to week out the clearly ambiguous ones
    for(const auto& hit3D : hitPairVector)
    {
        // Loop to try to remove ambiguous hits
        size_t n2DHitsIn3DHit(0);
        size_t nThisClusterOnly(0);
        size_t nOtherCluster(0);
        
        //        reco::ClusterParameters* otherCluster;
        const std::set<const reco::ClusterHit3D*>* otherClusterHits = 0;
        
        for(const auto& hit2D : hit3D->getHits())
        {
            if (!hit2D) continue;
            
            n2DHitsIn3DHit++;
            
            if (hit2DToClusterMap[hit2D].size() < 2) nThisClusterOnly = hit2DToClusterMap[hit2D][&clusterParams].size();
            else
            {
                for(const auto& clusterHitMap : hit2DToClusterMap[hit2D])
                {
                    if (clusterHitMap.first == &clusterParams) continue;
                    
                    if (clusterHitMap.second.size() > nOtherCluster)
                    {
                        nOtherCluster    = clusterHitMap.second.size();
                        //                        otherCluster     = clusterHitMap.first;
                        otherClusterHits = &clusterHitMap.second;
                    }
                }
            }
        }
        
        if (n2DHitsIn3DHit < 3 && nThisClusterOnly > 1 && nOtherCluster > 0)
        {
            bool skip3DHit(false);
            
            for(const auto& otherHit3D : *otherClusterHits)
            {
                size_t nOther2DHits(0);
                
                for(const auto& otherHit2D : otherHit3D->getHits())
                {
                    if (!otherHit2D) continue;
                    
                    nOther2DHits++;
                }
                
                if (nOther2DHits > 2)
                {
                    skip3DHit = true;
                    nRejectedHits++;
                    break;
                }
            }
            
            if (skip3DHit) continue;
            
        }
        
        goodHits.emplace_back(hit3D);
    }
    
    std::cout << "###>> Input " << nStartedWith << " hits, rejected: " << nRejectedHits << std::endl;
    
    hitPairVector.resize(goodHits.size());
    std::copy(goodHits.begin(),goodHits.end(),hitPairVector.begin());

    return;
}
    
struct HitPairClusterOrder
{
    bool operator()(const reco::ClusterParametersList::iterator& left, const reco::ClusterParametersList::iterator& right)
    {
        // Watch out for the case where two clusters can have the same number of hits!
        return (*left).getHitPairListPtr().size() > (*right).getHitPairListPtr().size();
    }
};
    
class SetCheckHitOrder
{
public:
    SetCheckHitOrder(const std::vector<size_t>& plane) : m_plane(plane) {}
    
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right) const
    {
        // Check if primary view's hit is on the same wire
        if (left->getWireIDs()[m_plane[0]] == right->getWireIDs()[m_plane[0]])
        {
            // Same wire but not same hit, order by primary hit time
            if (left->getHits()[m_plane[0]] && right->getHits()[m_plane[0]] && left->getHits()[m_plane[0]] != right->getHits()[m_plane[0]])
            {
                return left->getHits()[m_plane[0]]->getHit()->PeakTime() < right->getHits()[m_plane[0]]->getHit()->PeakTime();
            }
            
            // Primary view is same hit, look at next view's wire
            if (left->getWireIDs()[m_plane[1]] == right->getWireIDs()[m_plane[1]])
            {
                // Same wire but not same hit, order by secondary hit time
                if (left->getHits()[m_plane[1]] && right->getHits()[m_plane[1]] && left->getHits()[m_plane[1]] != right->getHits()[m_plane[1]])
                {
                    return left->getHits()[m_plane[1]]->getHit()->PeakTime() < right->getHits()[m_plane[1]]->getHit()->PeakTime();
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
    if (curCluster.size() > 2)
    {
        // Do a quick PCA to determine our parameter "alpha"
        reco::PrincipalComponents pca;
        m_pcaAlg.PCAAnalysis_3D(curCluster, pca);
        
        if (pca.getSvdOK())
        {
            const Eigen::Vector3f& pcaAxis  = pca.getEigenVectors().row(2);
            
            std::vector<size_t> closestPlane = {0, 0, 0 };
            std::vector<float>  bestAngle    = {0.,0.,0.};
            
            for(size_t plane = 0; plane < 3; plane++)
            {
                const std::vector<float>& wireDir = m_wireDir[plane];
                
                float dotProd = std::fabs(pcaAxis[0]*wireDir[0] + pcaAxis[1]*wireDir[1] + pcaAxis[2]*wireDir[2]);
                
                if (dotProd > bestAngle[0])
                {
                    bestAngle[2]    = bestAngle[1];
                    closestPlane[2] = closestPlane[1];
                    bestAngle[1]    = bestAngle[0];
                    closestPlane[1] = closestPlane[0];
                    closestPlane[0] = plane;
                    bestAngle[0]    = dotProd;
                }
                else if (dotProd > bestAngle[1])
                {
                    bestAngle[2]    = bestAngle[1];
                    closestPlane[2] = closestPlane[1];
                    closestPlane[1] = plane;
                    bestAngle[1]    = dotProd;
                }
                else
                {
                    closestPlane[2] = plane;
                    bestAngle[2]    = dotProd;
                }
            }
            
            // Get a copy of our 3D hits
            reco::HitPairListPtr localHitList = curCluster;
            
            // Sort the hits
            localHitList.sort(SetCheckHitOrder(closestPlane));
            
            // Ok, let's print it all and take a look
            std::cout << "********************************************************************************************" << std::endl;
            std::cout << "**>>>>> longest axis: " << closestPlane[0] << ", best angle: " << bestAngle[0] << std::endl;
            std::cout << "**>>>>> second  axis: " << closestPlane[1] << ", best angle: " << bestAngle[1] << std::endl;
            std::cout << " " << std::endl;
            
            reco::HitPairListPtr::iterator firstHitItr = localHitList.begin();
            reco::HitPairListPtr::iterator lastHitItr  = localHitList.begin();
            
            size_t bestPlane = closestPlane[0];
            
            reco::HitPairListPtr testList;
            
            while(firstHitItr != localHitList.end())
            {
                const reco::ClusterHit3D* currentHit = *firstHitItr;
                
                // Search for the last matching best view hit
                while(lastHitItr != localHitList.end())
                {
                    // If a different wire on the best view then we're certainly done
                    if (currentHit->getWireIDs()[bestPlane] != (*lastHitItr)->getWireIDs()[bestPlane]) break;
                    
                    // More subtle test to see if same wire but different hit (being careful of case of no hit)
                    if (currentHit->getHits()[bestPlane] && (*lastHitItr)->getHits()[bestPlane] && currentHit->getHits()[bestPlane] != (*lastHitItr)->getHits()[bestPlane]) break;
                    
                    // Yet event more subtle test...
                    if ((!(currentHit->getHits()[bestPlane]) && (*lastHitItr)->getHits()[bestPlane]) || (currentHit->getHits()[bestPlane] && !((*lastHitItr)->getHits()[bestPlane]))) break;
                    
                    // Not there yet...
                    lastHitItr++;
                }

                // How many hits in this chain?
//                size_t numHits(std::distance(firstHitItr,lastHitItr));
//                float  minOverlapFraction(0.);
                
//                if (numHits > 1)
//                {
//                    reco::HitPairListPtr::iterator bestMinOverlapItr = std::max_element(firstHitItr,lastHitItr,[](const auto& left, const auto& right){return left->getMinOverlapFraction() < right->getMinOverlapFraction();});
//
//                    minOverlapFraction = std::min(0.999*(*bestMinOverlapItr)->getMinOverlapFraction(),0.90);
//                }
                
                while(firstHitItr != lastHitItr)
                {
//                    if (currentHit->getMinOverlapFraction() > minOverlapFraction) testList.push_back(currentHit); //currentHit->setStatusBit(reco::ClusterHit3D::SKELETONHIT);
                    
                    currentHit = *++firstHitItr;
                }
                
                firstHitItr = lastHitItr;
            }
/*
            for(const auto& hit : localHitList)
            {
                std::cout << "- wires: ";
                
                for(size_t idx = 0; idx < 3; idx++)
                {
                    float viewTime = -1.;
                
                    if (hit->getHits()[closestView[idx]]) viewTime = hit->getHits()[closestView[idx]]->getTimeTicks();
                
                    std::cout << closestView[idx] << ":" << hit->getWireIDs()[closestView[idx]].Wire << " - " << viewTime << ", ";
                }
                
                bool isSkeleton = hit->getStatusBits() & reco::ClusterHit3D::SKELETONHIT;
                
                std::cout << "ave time: " << hit->getAvePeakTime() << ", min/max overlap: " << hit->getMinOverlapFraction() << "/" << hit->getMaxOverlapFraction() << ", tagged: " << isSkeleton << std::endl;
                
                if (isSkeleton) testList.push_back(hit);
            }
*/
            curCluster = testList;
        }
    }
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(MinSpanTreeAlg)
} // namespace lar_cluster3d

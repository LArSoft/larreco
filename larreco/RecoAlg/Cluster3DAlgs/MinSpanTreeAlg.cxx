/**
 *  @file   MinSpanTreeAlg.cxx
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/MinSpanTreeAlg.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

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
    
/**
 *  @brief define a kd tree node
 */
class MinSpanTreeAlg::KdTreeNode
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

MinSpanTreeAlg::MinSpanTreeAlg(fhicl::ParameterSet const &pset) :
    m_channelFilter(&art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider()),
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

MinSpanTreeAlg::~MinSpanTreeAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MinSpanTreeAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring         = pset.get<bool>  ("EnableMonitoring",  true  );
    m_minPairPts               = pset.get<size_t>("MinPairPts",                2     );
    m_timeAdvanceGap           = pset.get<double>("TimeAdvanceGap",           50.    );
    m_numSigmaPeakTime         = pset.get<double>("NumSigmaPeakTime",          3.    );
    m_pairSigmaPeakTime        = pset.get<double>("PairSigmaPeakTime",         3.    );
    m_pairMaxDistance          = pset.get<double>("PairMaxDistance",           0.95  );
    m_clusterMinHits           = pset.get<size_t>("ClusterMinHits",            3     );
    m_clusterMinUniqueFraction = pset.get<double>("ClusterMinUniqueFraction",  0.5   );
    m_clusterMaxLostFraction   = pset.get<double>("ClusterMaxLostFraction",    0.5   );
    
    art::ServiceHandle<geo::Geometry> geometry;
    
    m_geometry = &*geometry;
    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    m_timeVector.resize(NUMTIMEVALUES, 0.);
    
    // Determine the unit directon and normal vectors to the wires
    m_wireDir.resize(3);
    
    raw::ChannelID_t uChannel(0);
    std::vector<geo::WireID> uWireID = m_geometry->ChannelToWire(uChannel);
    const geo::WireGeo* uWireGeo = m_geometry->WirePtr(uWireID[0]);
    
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
    
    return;
}
    
void MinSpanTreeAlg::Cluster3DHits(reco::HitPairList&           hitPairList,
                                   reco::ClusterParametersList& clusterParametersList)
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
    
    KdTreeNodeList kdTreeNodeContainer;
    KdTreeNode     topNode = BuildKdTree(hitPairList, kdTreeNodeContainer);
    
    // Run DBScan to get candidate clusters
    RunPrimsAlgorithm(hitPairList, topNode, clusterParametersList);
    
    // Initial clustering is done, now trim the list and get output parameters
    BuildClusterInfo(clusterParametersList);
    
    // Test run the path finding algorithm
    for (auto& clusterParams : clusterParametersList) FindBestPathInCluster(clusterParams, topNode);
    
    mf::LogDebug("Cluster3D") << ">>>>> Cluster3DHits done, found " << clusterParametersList.size() << " clusters" << std::endl;
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
void MinSpanTreeAlg::RunPrimsAlgorithm(reco::HitPairList&           hitPairList,
                                       KdTreeNode&                  topNode,
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
    const reco::ClusterHit3D*   lastAddedHit = (*freeHitItr++).get();
    
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
        reco::EdgeList::iterator curEdgeItr = curEdgeList.begin();
        while(curEdgeItr != curEdgeList.end())
        {
            if (std::get<1>(*curEdgeItr)->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)
                curEdgeItr = curEdgeList.erase(curEdgeItr);
            else curEdgeItr++;
        }
        
        // Add the lastUsedHit to the current cluster
        curCluster->push_back(lastAddedHit);
        
        // Set up to find the list of nearest neighbors to the last used hit...
        CandPairVec candPairVec;
        double      bestDistance(std::numeric_limits<double>::max());

        // And find them... result will be an unordered list of neigbors
        FindNearestNeighbors(lastAddedHit, topNode, candPairVec, bestDistance);
        
        // Copy edges to the current list (but only for hits not already in a cluster)
        for(auto& pair : candPairVec)
            if (!(pair.second->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)) curEdgeList.push_back(reco::EdgeTuple(lastAddedHit,pair.second,pair.first));
        
        // If the edge list is empty then we have a complete cluster
        if (curEdgeList.empty())
        {
            std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
            std::cout << "**> Cluster idx: " << clusterIdx++ << " has " << curCluster->size() << " hits" << std::endl;

            // Look for the next "free" hit
            freeHitItr = std::find_if(freeHitItr,hitPairList.end(),[](const auto& hit){return !(hit->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED);});
            
            // If at end of input list we are done with all hits
            if (freeHitItr == hitPairList.end()) break;
            
            std::cout << "##################################################################>Processing another cluster" << std::endl;
            
            // Otherwise, get a new cluster and set up
            clusterParametersList.push_back(reco::ClusterParameters());
            
            curClusterItr = --clusterParametersList.end();
            
            curEdgeMap   = &(*curClusterItr).getHit3DToEdgeMap();
            curCluster   = &(*curClusterItr).getHitPairListPtr();
            lastAddedHit = (*freeHitItr++).get();
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
            lastAddedHit = std::get<1>(curEdgeList.front());
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
    double               bestQuality(0.);
    double               aveNumEdges(0.);
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
            double quality(0.);
            
            reco::HitPairListPtr tempList = DepthFirstSearch(curEdgeMap[hit].front(), curEdgeMap, quality);
            
            tempList.push_front(std::get<0>(curEdgeMap[hit].front()));
            
            if (quality > bestQuality)
            {
                longestCluster = tempList;
                bestQuality    = quality;
            }
            
            nIsolatedHits++;
        }
        
        aveNumEdges += double(curEdgeMap[hit].size());
        maxNumEdges  = std::max(maxNumEdges,curEdgeMap[hit].size());
    }
    
    aveNumEdges /= double(hitPairList.size());
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
    
void MinSpanTreeAlg::FindBestPathInCluster(reco::ClusterParameters& clusterParams, KdTreeNode& topNode) const
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
            const std::vector<double>& pcaAxis  = pca.getEigenVectors()[0];
            double                     pcaLen   = 1.5*sqrt(pca.getEigenValues()[0]);
            double                     pcaWidth = 1.5*sqrt(pca.getEigenValues()[1]);
            const double*              pcaPos   = pca.getAvePosition();
            double                     alpha    = std::min(1.,std::max(0.001,pcaWidth/pcaLen));
        
            // The first task is to find the list of hits which are "isolated"
            reco::HitPairListPtr isolatedHitList;
            for(const auto& hit : curCluster)
                if (!curEdgeMap[hit].empty() && curEdgeMap[hit].size() == 1)
                    isolatedHitList.emplace_back(std::get<0>(curEdgeMap[hit].front()));
        
            std::cout << "************* Finding best path with A* in cluster *****************" << std::endl;
            std::cout << "**> There are " << curCluster.size() << " hits, " << isolatedHitList.size() << " isolated hits, the alpha parameter is " << alpha << std::endl;
        
            // Goal is to now find separated pairs of isolated hits
            reco::EdgeList edgeList;
            reco::HitPairListPtr::iterator firstItr = isolatedHitList.begin();
            while(firstItr != isolatedHitList.end())
            {
                const reco::ClusterHit3D* firstHit = *firstItr++;
            
                const double* firstPos      = firstHit->getPosition();
                double        delta1stPca[] = {firstPos[0]-pcaPos[0],firstPos[1]-pcaPos[1],firstPos[2]-pcaPos[2]};
                double        firstPcaProj  = std::fabs(delta1stPca[0]*pcaAxis[0] + delta1stPca[1]*pcaAxis[1] + delta1stPca[2]*pcaAxis[2]);
            
                if (firstPcaProj < 0.75 * pcaLen) continue;
            
                for(reco::HitPairListPtr::iterator secondItr = firstItr; secondItr != isolatedHitList.end(); secondItr++)
                {
                    const double* secondPos     = (*secondItr)->getPosition();
                    double        delta2ndPca[] = {secondPos[0]-pcaPos[0],secondPos[1]-pcaPos[1],secondPos[2]-pcaPos[2]};
                    double        secondPcaProj = std::fabs(delta2ndPca[0]*pcaAxis[0] + delta2ndPca[1]*pcaAxis[1] + delta2ndPca[2]*pcaAxis[2]);
                
                    if (secondPcaProj < 0.75 * pcaLen) continue;
                
                    double deltaPos[] = {secondPos[0]-firstPos[0],secondPos[1]-firstPos[1],secondPos[2]-firstPos[2]};
                    double projection = std::fabs(deltaPos[0]*pcaAxis[0]+deltaPos[1]*pcaAxis[1]+deltaPos[2]*pcaAxis[2]);
                
                    edgeList.emplace_back(reco::EdgeTuple(firstHit,*secondItr,projection));
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

//                double cost(std::numeric_limits<double>::max());
                
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
                           double                    alpha,
                           KdTreeNode&               topNode,
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
            CandPairVec candPairVec;
            double      bestDistance(std::numeric_limits<double>::max());
            
            // And find them... result will be an unordered list of neigbors
            FindNearestNeighbors(currentNode, topNode, candPairVec, bestDistance);
            
//            std::cout << "**> found " << candPairVec.size() << " nearest neigbhors, bestDistance: " << bestDistance;
//            size_t nAdded(0);
            
            // Get tuple values for the current node
            const BestNodeTuple& currentNodeTuple = bestNodeMap.at(currentNode);
            double               currentNodeScore = std::get<1>(currentNodeTuple);
            
            for(auto& candPair : candPairVec)
            {
                // Ignore those nodes we're already aware of
                //if (std::find(closedList.begin(),closedList.end(),candPair.second) != closedList.end()) continue;
                if (candPair.second->getStatusBits() & reco::ClusterHit3D::PATHCHECKED) continue;
                
                double tentative_gScore = currentNodeScore + candPair.first;
                
                // Have we seen the candidate node before?
                BestNodeMap::iterator candNodeItr = bestNodeMap.find(candPair.second);
                
                if (candNodeItr == bestNodeMap.end())
                {
                    openList.push_back(candPair.second);
//                    nAdded++;
                }
                else if (tentative_gScore > std::get<1>(candNodeItr->second)) continue;

                // Experiment with modification to cost estimate
                const double* currentNodePos  = currentNode->getPosition();
                const double* nextNodePos     = candPair.second->getPosition();
                double        curNextDelta[]  = {nextNodePos[0]-currentNodePos[0], nextNodePos[1]-currentNodePos[1], nextNodePos[2]-currentNodePos[2]};
                
                const double* goalNodePos     = goalNode->getPosition();
                double        goalNextDelta[] = {goalNodePos[0]-nextNodePos[0], goalNodePos[1]-nextNodePos[1], goalNodePos[2]-nextNodePos[2]};
                
                double        curNextMag      = std::sqrt(curNextDelta[0]*curNextDelta[0]   + curNextDelta[1]*curNextDelta[1]   + curNextDelta[2]*curNextDelta[2]);
                double        goalNextMag     = std::sqrt(goalNextDelta[0]*goalNextDelta[0] + goalNextDelta[1]*goalNextDelta[1] + goalNextDelta[2]*goalNextDelta[2]);
                
                double        cosTheta        = (curNextDelta[0]*goalNextDelta[0] + curNextDelta[1]*goalNextDelta[1] + curNextDelta[2]*goalNextDelta[2]);
                
                if (cosTheta > 0. || cosTheta < 0.) cosTheta /= (curNextMag * goalNextMag);
                
                double        hWeight         = alpha*goalNextMag/std::max(0.01,0.5*(1.+cosTheta));

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
                                   double&                     showMeTheMoney) const
{
    reco::Hit3DToEdgeMap::const_iterator edgeListItr = hitToEdgeMap.find(std::get<1>(curEdge));
    
    showMeTheMoney = std::numeric_limits<double>::max();
    
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
            double currentCost(0.);
            
            LeastCostPath(edge,goalNode,hitToEdgeMap,bestNodeList,bestEdgeList,currentCost);
            
            if (currentCost < std::numeric_limits<double>::max())
            {
                showMeTheMoney = std::get<2>(edge) + currentCost;
                break;
            }
        }
    }
    
    if (showMeTheMoney < std::numeric_limits<double>::max())
    {
        bestNodeList.push_front(std::get<1>(curEdge));
        bestEdgeList.push_front(curEdge);
    }
    
    return;
}
    
double MinSpanTreeAlg::DistanceBetweenNodes(const reco::ClusterHit3D* node1,const reco::ClusterHit3D* node2) const
{
    const double* node1Pos    = node1->getPosition();
    const double* node2Pos    = node2->getPosition();
    double        deltaNode[] = {node1Pos[0]-node2Pos[0], node1Pos[1]-node2Pos[1], node1Pos[2]-node2Pos[2]};
    
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
    
    return std::sqrt(deltaNode[0]*deltaNode[0] + 0.09 * double(wireDeltas[2]*wireDeltas[2]));
 */
}

reco::HitPairListPtr MinSpanTreeAlg::DepthFirstSearch(const reco::EdgeTuple&      curEdge,
                                                      const reco::Hit3DToEdgeMap& hitToEdgeMap,
                                                      double&                     bestTreeQuality) const
{
    reco::HitPairListPtr hitPairListPtr;
    double               bestQuality(0.);
    double               curEdgeWeight = std::max(0.3,std::get<2>(curEdge));
    double               curEdgeProj(1./curEdgeWeight);
    
    reco::Hit3DToEdgeMap::const_iterator edgeListItr = hitToEdgeMap.find(std::get<1>(curEdge));
    
    if (edgeListItr != hitToEdgeMap.end())
    {
        // The input edge weight has quality factors applied, recalculate just the position difference
        const double* firstHitPos  = std::get<0>(curEdge)->getPosition();
        const double* secondHitPos = std::get<1>(curEdge)->getPosition();
        double        curEdgeVec[] = {secondHitPos[0]-firstHitPos[0],secondHitPos[1]-firstHitPos[1],secondHitPos[2]-firstHitPos[2]};
        double        curEdgeMag   = std::sqrt(curEdgeVec[0]*curEdgeVec[0]+curEdgeVec[1]*curEdgeVec[1]+curEdgeVec[2]*curEdgeVec[2]);
        
        curEdgeMag = std::max(0.1,curEdgeMag);
        
        for(const auto& edge : edgeListItr->second)
        {
            // skip the self reference
            if (std::get<1>(edge) == std::get<0>(curEdge)) continue;
            
            double quality(0.);
            
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
    
//------------------------------------------------------------------------------------------------------------------------------------------
MinSpanTreeAlg::KdTreeNode MinSpanTreeAlg::BuildKdTree(const reco::HitPairList& hitPairList,
                                                       KdTreeNodeList&          kdTreeNodeContainer) const
{
    
    // The first task is to build the kd tree
    cet::cpu_timer theClockBuildNeighborhood;
    
    if (m_enableMonitoring) theClockBuildNeighborhood.start();
    
    // The input is a list and we need to copy to a vector so we can sort ranges
    Hit3DVec hit3DVec;
    
    hit3DVec.reserve(hitPairList.size());

    for(const auto& hitPtr : hitPairList) hit3DVec.emplace_back(hitPtr.get());
    
    KdTreeNode topNode = BuildKdTree(hit3DVec.begin(), hit3DVec.end(), kdTreeNodeContainer);
    
    if (m_enableMonitoring)
    {
        theClockBuildNeighborhood.stop();
        m_timeVector[BUILDHITTOHITMAP] = theClockBuildNeighborhood.accumulated_real_time();
    }
    
    return topNode;
}
    
MinSpanTreeAlg::KdTreeNode& MinSpanTreeAlg::BuildKdTree(Hit3DVec::iterator first,
                                                        Hit3DVec::iterator last,
                                                        KdTreeNodeList&    kdTreeNodeContainer,
                                                        int                depth) const
{
    // Ok, so if the input list is more than one element then we have work to do... but if less then handle end condition
    if (std::distance(first,last) < 2)
    {
        if (first != last) kdTreeNodeContainer.emplace_back(KdTreeNode(*first));
        else               kdTreeNodeContainer.emplace_back(KdTreeNode());
//        if (first == last) std::cout << "********************************************* BAD NODE ***************************" << std::endl;
    }
    // Otherwise we need to keep splitting...
    else
    {
        // First task is to find "d" with the largest range. We need to find the min/max for the four dimensions
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxXPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[0] < right->getPosition()[0];});
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxYPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[1] < right->getPosition()[1];});
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxZPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[2] < right->getPosition()[2];});
        
        std::vector<double> rangeVec(3,0.);
        
        rangeVec[0] = (*minMaxXPair.second)->getPosition()[0] - (*minMaxXPair.first)->getPosition()[0];
        rangeVec[1] = (*minMaxYPair.second)->getPosition()[1] - (*minMaxYPair.first)->getPosition()[1];
        rangeVec[2] = (*minMaxZPair.second)->getPosition()[2] - (*minMaxZPair.first)->getPosition()[2];
        
        std::vector<double>::iterator maxRangeItr = std::max_element(rangeVec.begin(),rangeVec.end());
        
        size_t maxRangeIdx = std::distance(rangeVec.begin(),maxRangeItr);
        
        // Sort the list so we can do the split
        std::sort(first,last,[maxRangeIdx](const auto& left, const auto& right){return left->getPosition()[maxRangeIdx] < right->getPosition()[maxRangeIdx];});
        
        size_t             middleElem = std::distance(first,last) / 2;
        Hit3DVec::iterator middleItr  = first;
        
        std::advance(middleItr, middleElem);
        
        // Take care of the special case where the value of the median may be repeated so we actually want to make sure we point at the first occurence
        if (std::distance(first,middleItr) > 1)
        {
            while(middleItr != first+1)
            {
                if (!((*(middleItr-1))->getPosition()[maxRangeIdx] < (*middleItr)->getPosition()[maxRangeIdx])) middleItr--;
                else break;
            }
        }
        
        KdTreeNode::SplitAxis axis[]    = {KdTreeNode::xPlane,KdTreeNode::yPlane,KdTreeNode::zPlane};
        double                axisVal   = 0.5*((*middleItr)->getPosition()[maxRangeIdx] + (*(middleItr-1))->getPosition()[maxRangeIdx]);
        KdTreeNode&           leftNode  = BuildKdTree(first,     middleItr, kdTreeNodeContainer, depth+1);
        KdTreeNode&           rightNode = BuildKdTree(middleItr, last,      kdTreeNodeContainer, depth+1);
    
        kdTreeNodeContainer.push_back(KdTreeNode(axis[maxRangeIdx],axisVal,leftNode,rightNode));
    }
    
    return kdTreeNodeContainer.back();
}
    
size_t MinSpanTreeAlg::FindNearestNeighbors(const reco::ClusterHit3D* refHit, const KdTreeNode& node, CandPairVec& candPairVec, double& bestDist) const
{
    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        double hitSeparation(0.);
        int    wireDeltas[] = {0,0,0};
        
//        std::cout << "###>> nearest neighbor, refHit wires: " << refHit->getWireIDs()[0].Wire << "/" << refHit->getWireIDs()[1].Wire << "/" << refHit->getWireIDs()[2].Wire << ", compare to: " << node.getClusterHit3D()->getWireIDs()[0].Wire << "/" << node.getClusterHit3D()->getWireIDs()[1].Wire << "/" << node.getClusterHit3D()->getWireIDs()[2].Wire << std::endl;
        
        // Is this the droid we are looking for?
        if (refHit == node.getClusterHit3D()) bestDist = 0.5;  // This distance will grab neighbors with delta wire # = 1 in all three planes
        // This is the tight constraint on the hits
        else if (bestDist < std::numeric_limits<double>::max() && consistentPairs(refHit, node.getClusterHit3D(), hitSeparation, wireDeltas))
        {
            candPairVec.emplace_back(CandPair(hitSeparation,node.getClusterHit3D()));
            
            //bestDist = std::max(0.35,std::min(bestDist,hitSeparation));  // This insures we will always consider neighbors with wire # changing in 2 planes
            //bestDist = std::max(0.47,std::min(bestDist,hitSeparation));  // This insures we will always consider neighbors with wire # changing in 2 planes
            bestDist = std::max(0.85,std::min(bestDist,hitSeparation));  // This insures we will always consider neighbors with wire # changing in 2 planes
            
//            std::cout << "###>> nearest neighbor, refHit wires: " << refHit->getWireIDs()[0].Wire << "/" << refHit->getWireIDs()[1].Wire << "/" << refHit->getWireIDs()[2].Wire << ", compare to: " << node.getClusterHit3D()->getWireIDs()[0].Wire << "/" << node.getClusterHit3D()->getWireIDs()[1].Wire << "/" << node.getClusterHit3D()->getWireIDs()[2].Wire << std::endl;
            
//            std::cout << "  ~~~> cand " << candPairVec.size() << ", wire delta u: " << wireDeltas[0] << ", v: " << wireDeltas[1] << ", w: " << wireDeltas[2] << ", sep: " << hitSeparation << ", bestDist: " << bestDist << std::endl;
        }
    }
    // Otherwise we need to keep searching
    else
    {
        double refPosition = refHit->getPosition()[node.getSplitAxis()];
        
        if (refPosition < node.getAxisValue())
        {
            FindNearestNeighbors(refHit, node.leftTree(), candPairVec, bestDist);
            
            if (refPosition + bestDist > node.getAxisValue()) FindNearestNeighbors(refHit, node.rightTree(), candPairVec, bestDist);
        }
        else
        {
            FindNearestNeighbors(refHit, node.rightTree(), candPairVec, bestDist);
            
            if (refPosition - bestDist < node.getAxisValue()) FindNearestNeighbors(refHit, node.leftTree(), candPairVec, bestDist);
        }
    }
    
    return candPairVec.size();
}
    
bool MinSpanTreeAlg::FindEntry(const reco::ClusterHit3D* refHit, const KdTreeNode& node, CandPairVec& candPairVec, double& bestDist, bool& selfNotFound, int depth) const
{
    bool foundEntry(false);
    
    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        double hitSeparation(0.);
        int    wireDeltas[] = {0,0,0};
        
        // Is this the droid we are looking for?
        if (refHit == node.getClusterHit3D()) selfNotFound = false;
        
        // This is the tight constraint on the hits
        if (consistentPairs(refHit, node.getClusterHit3D(), hitSeparation, wireDeltas))
        {
            candPairVec.emplace_back(CandPair(hitSeparation,node.getClusterHit3D()));
            
            if (bestDist < std::numeric_limits<double>::max()) bestDist = std::max(bestDist,hitSeparation);
            else                                               bestDist = std::max(0.5,hitSeparation);
        }
        
        foundEntry = !selfNotFound;
    }
    // Otherwise we need to keep searching
    else
    {
        double refPosition = refHit->getPosition()[node.getSplitAxis()];
        
        if (refPosition < node.getAxisValue())
        {
            foundEntry = FindEntry(refHit, node.leftTree(),  candPairVec, bestDist, selfNotFound, depth+1);
            
            if (!foundEntry && refPosition + bestDist > node.getAxisValue()) foundEntry = FindEntry(refHit, node.rightTree(),  candPairVec, bestDist, selfNotFound, depth+1);
        }
        else
        {
            foundEntry = FindEntry(refHit, node.rightTree(),  candPairVec, bestDist, selfNotFound, depth+1);
            
            if (!foundEntry && refPosition - bestDist < node.getAxisValue()) foundEntry = FindEntry(refHit, node.leftTree(),  candPairVec, bestDist, selfNotFound, depth+1);
        }
    }
    
    return foundEntry;
}
    
bool MinSpanTreeAlg::FindEntryBrute(const reco::ClusterHit3D* refHit, const KdTreeNode& node, int depth) const
{
    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        // This is the tight constraint on the hits
        if (refHit == node.getClusterHit3D()) return true;
    }
    // Otherwise we need to keep searching
    else
    {
        if (FindEntryBrute(refHit, node.leftTree(),  depth+1)) return true;
        if (FindEntryBrute(refHit, node.rightTree(), depth+1)) return true;
    }
    
    return false;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
bool MinSpanTreeAlg::consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const
{
    // Strategy: We consider comparing "hit pairs" which may consist of 2 or 3 actual hits.
    //           Also, if only pairs, they can be U-V, U-W or V-W so we can't assume which views we have
    //           So do a simplified comparison:
    //           1) compare the pair times and require "overlap" (in the sense of hit pair production)
    //           2) look at distance between pairs in each of the wire directions
    
    double pair1PeakTime = pair1->getAvePeakTime();
    double pair1Width    = m_pairSigmaPeakTime * pair1->getSigmaPeakTime();
    double pair2PeakTime = pair2->getAvePeakTime();
    double pair2Width    = m_pairSigmaPeakTime * pair2->getSigmaPeakTime();
    
    double maxUpper      = std::min(pair1PeakTime+pair1Width,pair2PeakTime+pair2Width);
    double minLower      = std::max(pair1PeakTime-pair1Width,pair2PeakTime-pair2Width);
    double pairOverlap   = maxUpper - minLower;
    
    const std::vector<geo::WireID>& pair1Wires = pair1->getWireIDs();
    const std::vector<geo::WireID>& pair2Wires = pair2->getWireIDs();
    
    if ((pair1Wires[0].Wire >  778 && pair1Wires[0].Wire <  784) && (pair2Wires[0].Wire >  778 && pair2Wires[0].Wire <  784) &&
        (pair1Wires[1].Wire >  900 && pair1Wires[1].Wire <  907) && (pair2Wires[1].Wire >  900 && pair2Wires[1].Wire <  907) &&
        (pair1Wires[2].Wire > 1008 && pair1Wires[2].Wire < 1014) && (pair2Wires[2].Wire > 1008 && pair2Wires[2].Wire < 1014))
    {
        std::cout << "++++ Checking pairs: " << pair1Wires[0].Wire << "/" << pair1Wires[1].Wire << "/" << pair1Wires[2].Wire << ", "  << pair2Wires[0].Wire << "/" << pair2Wires[1].Wire << "/" << pair2Wires[2].Wire << std::endl;
        std::cout << "     pair1 time/width: " << pair1PeakTime << ", " << pair1Width << ", pair2 time/width: " << pair2PeakTime << ", " << pair2Width << ", overlap: " << pairOverlap << std::endl;
    }
    
    // Loose constraint to weed out the obviously bad combinations
    if (pairOverlap > 0.1)
    {
        hitSeparation = DistanceBetweenNodes(pair1,pair2);
        
//        size_t hitCount(0);
        
        // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
        for(size_t idx = 0; idx < 3; idx++)
        {
            wireDeltas[idx] = std::abs(int(pair1->getWireIDs()[idx].Wire) - int(pair2->getWireIDs()[idx].Wire));
            
//            if (pair1->getHits()[idx]) hitCount++;
//            if (pair2->getHits()[idx]) hitCount++;
        }
        
        // put wire deltas in order...
        std::sort(wireDeltas, wireDeltas + 3);
        
        // Requirement to be considered a nearest neighbor
        if (wireDeltas[0] < 2 && wireDeltas[1] < 2 && wireDeltas[2] < 3)
        {
//            double overlapFraction = 0.5 * pairOverlap / std::min(pair1Width,pair2Width);

//            hitSeparation /= overlapFraction;
            
            // Scale the hit separation by the number of missing wires
//            for(size_t idx = 0; idx < 6-hitCount; idx++) hitSeparation *= 2.0; //1.1;
            
//            if (wireDeltas[0] == 0) hitSeparation *= 2.0;
            
            hitSeparation = std::max(0.0001,hitSeparation);
            
            return true;
        }
    }
    
    return false;
}
    
void MinSpanTreeAlg::PruneAmbiguousHits(reco::ClusterParameters& clusterParams, Hit2DToClusterMap& hit2DToClusterMap) const
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
        const std::set<const reco::ClusterHit3D*>* otherClusterHits;
        
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
//    SetCheckHitOrder()                                : m_view(std::ve)    {}
    SetCheckHitOrder(const std::vector<size_t>& view) : m_view(view) {}
    
    bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right) const
    {
        // Check if primary view's hit is on the same wire
        if (left->getWireIDs()[m_view[0]] == right->getWireIDs()[m_view[0]])
        {
            // Same wire but not same hit, order by primary hit time
            if (left->getHits()[m_view[0]] && right->getHits()[m_view[0]] && left->getHits()[m_view[0]] != right->getHits()[m_view[0]])
            {
                return left->getHits()[m_view[0]]->getHit().PeakTime() < right->getHits()[m_view[0]]->getHit().PeakTime();
            }
            
            // Primary view is same hit, look at next view's wire
            if (left->getWireIDs()[m_view[1]] == right->getWireIDs()[m_view[1]])
            {
                // Same wire but not same hit, order by secondary hit time
                if (left->getHits()[m_view[1]] && right->getHits()[m_view[1]] && left->getHits()[m_view[1]] != right->getHits()[m_view[1]])
                {
                    return left->getHits()[m_view[1]]->getHit().PeakTime() < right->getHits()[m_view[1]]->getHit().PeakTime();
                }
            
                // All that is left is the final view... and this can't be the same hit... (else it is the same 3D hit)
                return left->getWireIDs()[m_view[2]] < right->getWireIDs()[m_view[2]];
            }
            
            return left->getWireIDs()[m_view[1]] < right->getWireIDs()[m_view[1]];
        }

        // Order by primary view's wire number
        return left->getWireIDs()[m_view[0]] < right->getWireIDs()[m_view[0]];
    }
    
private:
    const std::vector<size_t>& m_view;
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
            const std::vector<double>& pcaAxis  = pca.getEigenVectors()[0];
            
            std::vector<size_t> closestView = {0, 0, 0 };
            std::vector<float>  bestAngle   = {0.,0.,0.};
            
            for(size_t view = 0; view < 3; view++)
            {
                const std::vector<float>& wireDir = m_wireDir[view];
                
                float dotProd = std::fabs(pcaAxis[0]*wireDir[0] + pcaAxis[1]*wireDir[1] + pcaAxis[2]*wireDir[2]);
                
                if (dotProd > bestAngle[0])
                {
                    bestAngle[2]   = bestAngle[1];
                    closestView[2] = closestView[1];
                    bestAngle[1]   = bestAngle[0];
                    closestView[1] = closestView[0];
                    closestView[0] = view;
                    bestAngle[0]   = dotProd;
                }
                else if (dotProd > bestAngle[1])
                {
                    bestAngle[2]   = bestAngle[1];
                    closestView[2] = closestView[1];
                    closestView[1] = view;
                    bestAngle[1]   = dotProd;
                }
                else
                {
                    closestView[2] = view;
                    bestAngle[2]   = dotProd;
                }
            }
            
            // Get a copy of our 3D hits
            reco::HitPairListPtr localHitList = curCluster;
            
            // Sort the hits
            localHitList.sort(SetCheckHitOrder(closestView));
            
            // Ok, let's print it all and take a look
            std::cout << "********************************************************************************************" << std::endl;
            std::cout << "**>>>>> longest axis: " << closestView[0] << ", best angle: " << bestAngle[0] << std::endl;
            std::cout << "**>>>>> second  axis: " << closestView[1] << ", best angle: " << bestAngle[1] << std::endl;
            std::cout << " " << std::endl;
            
            reco::HitPairListPtr::iterator firstHitItr = localHitList.begin();
            reco::HitPairListPtr::iterator lastHitItr  = localHitList.begin();
            
            size_t bestView = closestView[0];
            
            reco::HitPairListPtr testList;
            
            while(firstHitItr != localHitList.end())
            {
                const reco::ClusterHit3D* currentHit = *firstHitItr;
                
                // Search for the last matching best view hit
                while(lastHitItr != localHitList.end())
                {
                    // If a different wire on the best view then we're certainly done
                    if (currentHit->getWireIDs()[bestView] != (*lastHitItr)->getWireIDs()[bestView]) break;
                    
                    // More subtle test to see if same wire but different hit (being careful of case of no hit)
                    if (currentHit->getHits()[bestView] && (*lastHitItr)->getHits()[bestView] && currentHit->getHits()[bestView] != (*lastHitItr)->getHits()[bestView]) break;
                    
                    // Yet event more subtle test...
                    if ((!(currentHit->getHits()[bestView]) && (*lastHitItr)->getHits()[bestView]) || (currentHit->getHits()[bestView] && !((*lastHitItr)->getHits()[bestView]))) break;
                    
                    // Not there yet...
                    lastHitItr++;
                }
                
                // How many hits in this chain?
                size_t numHits(std::distance(firstHitItr,lastHitItr));
                double minOverlapFraction(0.);
                
                if (numHits > 1)
                {
                    reco::HitPairListPtr::iterator bestMinOverlapItr = std::max_element(firstHitItr,lastHitItr,[](const auto& left, const auto& right){return left->getMinOverlapFraction() < right->getMinOverlapFraction();});
                    
                    minOverlapFraction = std::min(0.999*(*bestMinOverlapItr)->getMinOverlapFraction(),0.90);
                }
                
                while(firstHitItr != lastHitItr)
                {
                    if (currentHit->getMinOverlapFraction() > minOverlapFraction) testList.push_back(currentHit); //currentHit->setStatusBit(reco::ClusterHit3D::SKELETONHIT);
                    
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
    
void MinSpanTreeAlg::BuildClusterInfo(reco::ClusterParametersList& clusterParametersList) const
{
    /**
     *  @brief Given a list of a list of candidate cluster hits, build these out into the intermediate
     *         3D cluster objects to pass to the final stage
     *
     *         Note that this routine will also reject unworthy clusters, in particular those that share too
     *         many hits with other clusters. The criteria is that a larger cluster (more hits) will be superior
     *         to a smaller one, if the smaller one shares too many hits with the larger it is zapped.
     *         *** THIS IS AN AREA FOR CONTINUED STUDY ***
     */
    cet::cpu_timer theClockBuildClusterInfo;
    
    if (m_enableMonitoring) theClockBuildClusterInfo.start();
    
    // This is a remote possibility but why not check?
    if (!clusterParametersList.empty())
    {
        // We want to order our clusters on by largest (most number hits) to smallest. So, we'll loop through the clusters,
        // weeding out the unwanted ones and keep track of things in a set of "good" clusters which we'll order
        // by cluster size.
        clusterParametersList.sort();
        
        // The smallest clusters are now at the end, drop those off the back that are less than the mininum necessary
        while(!clusterParametersList.empty() && clusterParametersList.back().getHitPairListPtr().size() < m_clusterMinHits) clusterParametersList.pop_back();
        
        // The next step is to build out a mapping of all 2D hits to clusters
        // Keep track of where the hits get distributed...
        Hit2DToClusterMap hit2DToClusterMap;
        
        reco::ClusterParametersList::iterator clusterItr = clusterParametersList.begin();
        
        for(auto& clusterParams : clusterParametersList)
        {
            for(const auto& hit3D : clusterParams.getHitPairListPtr())
            {
                for(const auto& hit2D : hit3D->getHits())
                {
                    if (!hit2D) continue;
                    
                    hit2DToClusterMap[hit2D][&clusterParams].insert(hit3D);
                }
            }
        }
        
        // Ok, spin through again to remove ambiguous hits
//        for(auto& clusterParams : clusterParametersList) PruneAmbiguousHits(clusterParams,hit2DToClusterMap);
        
        // What remains is an order set of clusters, largest first
        // Now go through and obtain cluster parameters
        clusterItr = clusterParametersList.begin();
        
        while(clusterItr != clusterParametersList.end())
        {
            // Dereference for ease...
            reco::ClusterParameters& clusterParams = *clusterItr;
            
            // Do the actual work of filling the parameters
            FillClusterParams(clusterParams, hit2DToClusterMap, m_clusterMinUniqueFraction, m_clusterMaxLostFraction);
            
            // If this cluster is rejected then the parameters will be empty
            if (clusterParams.getClusterParams().empty() || !clusterParams.getFullPCA().getSvdOK())
            {
                clusterItr = clusterParametersList.erase(clusterItr);
            }
            else clusterItr++;
        }
    }
    
    if (m_enableMonitoring)
    {
        theClockBuildClusterInfo.stop();
        
        m_timeVector[BUILDCLUSTERINFO] = theClockBuildClusterInfo.accumulated_real_time();
    }
    
    return;
}

void MinSpanTreeAlg::FillClusterParams(reco::ClusterParameters& clusterParams, Hit2DToClusterMap& hit2DToClusterMap, double minUniqueFrac, double maxLostFrac) const
{
    /**
     *  @brief Given a list of hits fill out the remaining parameters for this cluster and evaluate the
     *         candidate's worthiness to achieve stardom in the event display
     */
    
    // Recover the HitPairListPtr from the input clusterParams (which will be the
    // only thing that has been provided)
    reco::HitPairListPtr& hitPairVector = clusterParams.getHitPairListPtr();
    
    // To be sure, we should clear the other data members
    clusterParams.getClusterParams().clear();
    clusterParams.getFullPCA() = reco::PrincipalComponents();
    
    // A test of the emergency broadcast system...
//    FindBestPathInCluster(clusterParams);
    CheckHitSorting(clusterParams);
    
    // See if we can avoid duplicates by temporarily transferring to a set
    //std::set<const reco::ClusterHit2D*> hitSet;
    std::vector<const reco::ClusterHit2D*> hitSet;
    
    size_t nTotalHits[]  = {0,0,0};
    size_t nUniqueHits[] = {0,0,0};
    size_t nLostHits[]   = {0,0,0};
    size_t nMultShared2DHits(0);
    size_t nAllHitsShared(0);
    
    std::map<reco::ClusterParameters*,int> clusterHitCountMap;
    
    // Create a list to hold 3D hits which are already in use (criteria below)
    reco::HitPairListPtr usedHitPairList;
    
    // First loop through the 3D hits
    // The goal of this loop is to build a set of unique hits from the hit pairs (which may contain many
    // ambiguous duplicate combinations).
    // The secondary goal is to remove 3D hits marked by hit arbitration to be tossed
    for(const auto& hit3D : hitPairVector)
    {
        size_t nMultClusters(0);
        size_t nHits2D(0);
        size_t nHitsUsed[] = {0,0,0};
        
        // loop over the hits in this 3D Cluster hit
        for(const auto& hit2D : hit3D->getHits())
        {
            if (!hit2D) continue;
            size_t view = hit2D->getHit().View();

            if (hit2D->getStatusBits() & reco::ClusterHit2D::USED) nHitsUsed[view]++;
            else                                                   nUniqueHits[view]++;
            
            // Is this 2D hit shared?
            if (hit2DToClusterMap[hit2D].size() > 1) nMultClusters++;
            for(auto& clusterCntPair : hit2DToClusterMap[hit2D]) clusterHitCountMap[clusterCntPair.first]++;
            
            nTotalHits[view]++;
            nHits2D++;
        }
        
        size_t nHitsAlreadyUsed = std::accumulate(nHitsUsed,nHitsUsed+3,0);
        
        if (nMultClusters > 1)        nMultShared2DHits++;
        if (nMultClusters == nHits2D) nAllHitsShared++;
        
        for(size_t idx=0;idx<3;idx++)
        {
            if (nHitsAlreadyUsed < nHits2D)
            {
                //if (hit3D->getHits()[idx]) hitSet.insert(hit3D->getHits()[idx]);
                if (hit3D->getHits()[idx]) hitSet.push_back(hit3D->getHits()[idx]);
            }
            else nLostHits[idx] += nHitsUsed[idx];
        }
        
        if (nHitsAlreadyUsed == nHits2D) usedHitPairList.emplace_back(hit3D);
    }
    
    int numTotal      = std::accumulate(nTotalHits,nTotalHits+3,0);
    int numUniqueHits = std::accumulate(nUniqueHits,nUniqueHits+3,0);
    int numLostHits   = std::accumulate(nLostHits,nLostHits+3,0);
    
    std::cout << "*********************************************************************" << std::endl;
    std::cout << "**--> cluster: " << &clusterParams << " has " << hitPairVector.size() << " 3D hits, " << numTotal << " 2D hits, match: " << clusterHitCountMap[&clusterParams] << ", shared: " << nMultShared2DHits << ", all: " << nAllHitsShared << std::endl;
    
    for(const auto& clusterCnt : clusterHitCountMap)
    {
        if (clusterCnt.first == &clusterParams) continue;
        std::cout << "      --> cluster " << clusterCnt.first << ", # hits: " << clusterCnt.second << std::endl;
    }
    
    // If we have something left then at this point we make one more check
    // This check is intended to weed out clusters made from isolated groups of ambiguous hits which
    // really belong to a larger cluster
    if (numUniqueHits > 3 && nMultShared2DHits < hitPairVector.size())
    {
        // Look at reject to accept ratio
        //double rejectToAccept = double(numRejected) / double(numAccepted);
        double acceptRatio = double(numUniqueHits) / double(numTotal);
        double lostRatio   = double(numLostHits)   / double(numTotal);
        
        // Also consider the number of hits shared on a given view...
        std::vector<double> uniqueHitVec(3,0.);
        
        for(size_t idx = 0; idx < 3; idx++) uniqueHitVec[idx] = double(nUniqueHits[idx]) / std::max(double(nTotalHits[idx]),1.);
        
        std::sort(uniqueHitVec.begin(),uniqueHitVec.end());
        
        //        double midHitRatio = uniqueHitVec[1];
        
        //        std::cout << "--> # 3D Hits: " << hitPairVector.size() << ", nTot: " << numTotal << ", unique: " << numUniqueHits << ", lost: " << numLostHits << ", accept: " << acceptRatio << ", lost: " << lostRatio << ", mid: " << midHitRatio << ", rats: " << uniqueHitVec[0] << "/" << uniqueHitVec[1] << "/" << uniqueHitVec[2] << std::endl;
        
        acceptRatio = 0.;
        lostRatio   = 0.;
        if(uniqueHitVec[1] > 0.1 && uniqueHitVec[2] > 0.5) acceptRatio = 1.;
        
        // Arbitrary rejection criteria... need to understand
        // Anyway, if we get past this we're making a cluster
        //if (rejectToAccept < rejectFraction)
        if (acceptRatio > minUniqueFrac && lostRatio < maxLostFrac)  // lostRatio cut was 1. - off
        {
            // Add the "good" hits to our cluster parameters
            for(const auto& hit2D : hitSet)
            {
                hit2D->setStatusBit(reco::ClusterHit2D::USED);
                clusterParams.UpdateParameters(hit2D);
            }
            
            size_t nViewsWithHits    = (clusterParams.getClusterParams()[geo::kU].m_hitVector.size() > 0 ? 1 : 0)
            + (clusterParams.getClusterParams()[geo::kV].m_hitVector.size() > 0 ? 1 : 0)
            + (clusterParams.getClusterParams()[geo::kW].m_hitVector.size() > 0 ? 1 : 0);
            size_t nViewsWithMinHits = (clusterParams.getClusterParams()[geo::kU].m_hitVector.size() > 2 ? 1 : 0)
            + (clusterParams.getClusterParams()[geo::kV].m_hitVector.size() > 2 ? 1 : 0)
            + (clusterParams.getClusterParams()[geo::kW].m_hitVector.size() > 2 ? 1 : 0);
            //            // Final selection cut, need at least 3 hits each view
            //            if (nViewsWithHits == 3 && nViewsWithMinHits > 1)
            // Final selection cut, need at least 3 hits each view for at least 2 views
            if (nViewsWithHits > 1 && nViewsWithMinHits > 1)
            {
                // First task is to remove the hits already in use
                if (!usedHitPairList.empty())
                {
                    hitPairVector.sort();
                    usedHitPairList.sort();
                    
                    reco::HitPairListPtr::iterator newListEnd =
                    std::set_difference(hitPairVector.begin(),   hitPairVector.end(),
                                        usedHitPairList.begin(), usedHitPairList.end(),
                                        hitPairVector.begin() );
                    
                    hitPairVector.erase(newListEnd, hitPairVector.end());
                }
                
                // First stage of feature extraction runs here
                m_pcaAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());
                
                // Must have a valid pca
                if (clusterParams.getFullPCA().getSvdOK())
                {
                    // If any hits were thrown away, see if we can rescue them
                    if (!usedHitPairList.empty())
                    {
                        double maxDoca = 2. * sqrt(clusterParams.getFullPCA().getEigenValues()[1]);
                        
                        if (maxDoca < 5.)
                        {
                            size_t curHitVectorSize = hitPairVector.size();
                            
                            m_pcaAlg.PCAAnalysis_calc3DDocas(usedHitPairList, clusterParams.getFullPCA());
                            
                            for(const auto& hit3D : usedHitPairList)
                                if (hit3D->getDocaToAxis() < maxDoca) hitPairVector.push_back(hit3D);
                            
                            if (hitPairVector.size() > curHitVectorSize)
                                m_pcaAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());
                        }
                    }
                    
                    // Set the skeleton PCA to make sure it has some value
                    clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();
                }
            }
        }
    }
    
    return;
}
    

} // namespace lar_cluster3d

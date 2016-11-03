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

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

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
    
    // Determine the unit directon and normal vectors to the wires
    m_wireDir.resize(3);
    m_wireNormal.resize(3);
    
    raw::ChannelID_t uChannel(0);
    std::vector<geo::WireID> uWireID = m_geometry->ChannelToWire(uChannel);
    const geo::WireGeo* uWireGeo = m_geometry->WirePtr(uWireID[0]);
    
    TVector3 uWireDir = uWireGeo->Direction();
    
    m_wireDir[0].resize(3);
    m_wireDir[0][0] = uWireDir[0];
    m_wireDir[0][1] = uWireDir[1];
    m_wireDir[0][2] = uWireDir[2];
    
    m_wireNormal[0].resize(3);
    m_wireNormal[0][0] = 0.;
    m_wireNormal[0][1] = -uWireDir[2];
    m_wireNormal[0][2] =  uWireDir[1];
    
    raw::ChannelID_t vChannel(2400);
    std::vector<geo::WireID> vWireID = m_geometry->ChannelToWire(vChannel);
    const geo::WireGeo* vWireGeo = m_geometry->WirePtr(vWireID[0]);
    
    TVector3 vWireDir = vWireGeo->Direction();
    
    m_wireDir[1].resize(3);
    m_wireDir[1][0] = vWireDir[0];
    m_wireDir[1][1] = vWireDir[1];
    m_wireDir[1][2] = vWireDir[2];
    
    m_wireNormal[1].resize(3);
    m_wireNormal[1][0] = 0.;
    m_wireNormal[1][1] = -vWireDir[2];
    m_wireNormal[1][2] =  vWireDir[1];
    
    m_wireDir[2].resize(3);
    m_wireDir[2][0] = 0.;
    m_wireDir[2][1] = 1.;
    m_wireDir[2][2] = 0.;
    
    m_wireNormal[2].resize(3);
    m_wireNormal[2][0] = 0.;
    m_wireNormal[2][1] = 0.;
    m_wireNormal[2][2] = 1.;
    
    m_wirePitch[0] = m_geometry->WirePitch(0);
    m_wirePitch[1] = m_geometry->WirePitch(1);
    m_wirePitch[2] = m_geometry->WirePitch(2);
    
    m_timeVector.resize(NUMTIMEVALUES, 0.);
}
    
void MinSpanTreeAlg::ClusterHitsDBScan(reco::HitPairList&           hitPairList,
                                       reco::HitPairClusterMap&     hitPairClusterMap,
                                       reco::ClusterParametersList& clusterParametersList)
{
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */
    m_timeVector.resize(NUMTIMEVALUES, 0.);
    
    // The container of pairs and those in each pair's epsilon neighborhood
    EpsPairNeighborhoodMapVec epsPairNeighborhoodMapVec;
    
    epsPairNeighborhoodMapVec.resize(hitPairList.size(), EpsPairNeighborhoodMapPair(0,EpsPairNeighborhoodPair()));  //<-- initialize too!
    
    // DBScan is driven of its "epsilon neighborhood". Computing adjacency within DBScan can be time
    // consuming so the idea is the prebuild the adjaceny map and then run DBScan.
    // The following call does this work
    //BuildNeighborhoodMap(hitPairList, epsPairNeighborhoodMapVec);
    
    KdTreeNodeList kdTreeNodeContainer;
    KdTreeNode     topNode = BuildKdTree(hitPairList, kdTreeNodeContainer);
    
    // Run DBScan to get candidate clusters
    //RunDBScan(epsPairNeighborhoodMapVec, hitPairClusterMap);
    RunPrimsAlgorithm(hitPairList, topNode, hitPairClusterMap);
    
    // Initial clustering is done, now trim the list and get output parameters
    BuildClusterInfo(hitPairClusterMap, clusterParametersList);
    
    mf::LogDebug("Cluster3D") << ">>>>> DBScan done, found " << hitPairClusterMap.size() << " clusters" << std::endl;
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
void MinSpanTreeAlg::RunPrimsAlgorithm(reco::HitPairList& hitPairList, KdTreeNode& topNode, reco::HitPairClusterMap& hitPairClusterMap) const
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
    MSTEdgeList         curEdgeList;
    MSTClusterToEdgeMap hitToEdgeMap;
    
    // Get the first point
    reco::HitPairList::iterator freeHitItr   = hitPairList.begin();
    const reco::ClusterHit3D*   lastAddedHit = (*freeHitItr++).get();
    
    lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);
    
    // Get the first cluster
    MST3DHitToEdgeMap*    curEdgeMap = &hitToEdgeMap[clusterIdx];
    reco::HitPairListPtr* curCluster = &hitPairClusterMap[clusterIdx++];
    
    // Loop until all hits have been associated to a cluster
    while(1)
    {
        // Add the lastUsedHit to the current cluster
        curCluster->push_back(lastAddedHit);

        // Set up to find the list of nearest neighbors to the last used hit...
        CandPairVec candPairVec;
        double      bestDistance(std::numeric_limits<double>::max());

        // And find them... result will be an unordered list of neigbors
        FindNearestNeighbors(lastAddedHit, topNode, candPairVec, bestDistance);
        
        // Copy edges to the current list (but only for hits not already in a cluster)
        for(auto& pair : candPairVec) curEdgeList.push_back(MSTEdgeTuple(lastAddedHit,pair.second,pair.first));
//            if (!(pair.second->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)) curEdgeList.push_back(MSTEdgeTuple(lastUsedHit,pair.second,pair.first));
        
        // If the edge list is empty then we have a complete cluster
        if (curEdgeList.empty())
        {
//            std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
//            std::cout << "**> Cluster idx: " << clusterIdx-1 << " has " << curCluster->size() << " hits, " << hitPairClusterMap[clusterIdx-1].size() << std::endl;
/*
            reco::HitPairListPtr longestCluster;
            double               deepest(0.);
            double               aveNumEdges(0.);
            size_t               maxNumEdges(0);
            size_t               nIsolatedHits(0);
            
            // Do some spelunking...
            for(const auto& hit : *curCluster)
            {
                if (!(*curEdgeMap)[hit].empty() && (*curEdgeMap)[hit].size() < 3) //== 1)
                {
                    double depth(0.);
                    
                    reco::HitPairListPtr tempList = DepthFirstSearch((*curEdgeMap)[hit].front(), *curEdgeMap, depth);
                    
                    tempList.push_front(std::get<0>((*curEdgeMap)[hit].front()));
                    
                    if (depth > deepest)
                    {
                        longestCluster = tempList;
                        deepest        = depth;
                    }
                    
                    nIsolatedHits++;
                }
                
                aveNumEdges += double((*curEdgeMap)[hit].size());
                maxNumEdges  = std::max(maxNumEdges,(*curEdgeMap)[hit].size());
            }
            
            aveNumEdges /= double(curCluster->size());
//            std::cout << "----> # isolated hits: " << nIsolatedHits << ", longest branch: " << longestCluster.size() << ", cluster size: " << curCluster->size() << ", ave # edges: " << aveNumEdges << ", max: " << maxNumEdges << std::endl;
            
            if (!longestCluster.empty())
            {
                *curCluster = longestCluster;
                
//                std::cout << "        ====> new cluster size: " << curCluster->size() << ", " << hitPairClusterMap[clusterIdx-1].size() << std::endl;
            }
*/            
            // Look for the next "free" hit
            freeHitItr = std::find_if(freeHitItr,hitPairList.end(),[](const auto& hit){return !(hit->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED);});
            
            // If at end of input list we are done with all hits
            if (freeHitItr == hitPairList.end()) break;
            
//            std::cout << "##################################################################>Processing another cluster" << std::endl;
            
            // Otherwise, get a new cluster and set up
            curEdgeMap   = &hitToEdgeMap[clusterIdx];
            curCluster   = &hitPairClusterMap[clusterIdx++];
            lastAddedHit = (*freeHitItr++).get();
            
            lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);
        }
        // Otherwise we are still processing the current cluster
        else
        {
            // Sort the list of edges by distance
            curEdgeList.sort([](const auto& left,const auto& right){return std::get<2>(left) < std::get<2>(right);});
            
            // Populate the map with the edges...
            (*curEdgeMap)[std::get<0>(curEdgeList.front())].push_back(curEdgeList.front());
            (*curEdgeMap)[std::get<1>(curEdgeList.front())].push_back(curEdgeList.front());

//            std::cout << "====> Add hit to cluster, hit in cluster: " << std::get<0>(curEdgeList.front())->getWireIDs()[0].Wire << "/" << std::get<0>(curEdgeList.front())->getWireIDs()[1].Wire << "/" << std::get<0>(curEdgeList.front())->getWireIDs()[2].Wire << ", # links: " << (*curEdgeMap)[std::get<0>(curEdgeList.front())].size() << ", added: " << std::get<1>(curEdgeList.front())->getWireIDs()[0].Wire << "/" << std::get<1>(curEdgeList.front())->getWireIDs()[1].Wire << "/" << std::get<1>(curEdgeList.front())->getWireIDs()[2].Wire << ", # links: " << (*curEdgeMap)[std::get<1>(curEdgeList.front())].size() << ", edge len: " << std::get<2>(curEdgeList.front()) << std::endl;
//            std::cout << "      last added hit in cluster: " << lastAddedHit->getWireIDs()[0].Wire << "/" << lastAddedHit->getWireIDs()[1].Wire << "/" << lastAddedHit->getWireIDs()[2].Wire << ", # links: " << (*curEdgeMap)[lastAddedHit].size() << ", current edge list size: " << curEdgeList.size() << std::endl;
//            std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
            
            // Update the last hit to be added to the collection
            lastAddedHit = std::get<1>(curEdgeList.front());
            
            lastAddedHit->setStatusBit(reco::ClusterHit3D::CLUSTERATTACHED);
            
            // Purge the current list to get rid of edges which point to hits already in the cluster
            MSTEdgeList::iterator curEdgeItr = curEdgeList.begin();
            while(curEdgeItr != curEdgeList.end())
            {
                if (std::get<1>(*curEdgeItr)->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED)
                {
                    curEdgeItr = curEdgeList.erase(curEdgeItr);
                }
                else curEdgeItr++;
            }
        }
    }
    
    if (m_enableMonitoring)
    {
        theClockDBScan.stop();
        
        m_timeVector[RUNDBSCAN] = theClockDBScan.accumulated_real_time();
    }
    
    return;
}
    
reco::HitPairListPtr MinSpanTreeAlg::DepthFirstSearch(const MSTEdgeTuple&      curEdge,
                                                      const MST3DHitToEdgeMap& hitToEdgeMap,
                                                      double&                  depth) const
{
    reco::HitPairListPtr hitPairListPtr;
    double               deepest(0.);
    
    MST3DHitToEdgeMap::const_iterator edgeListItr = hitToEdgeMap.find(std::get<1>(curEdge));
    
    if (edgeListItr != hitToEdgeMap.end())
    {
//        const double* firstHitPos  = std::get<0>(curEdge)->getPosition();
//        const double* secondHitPos = std::get<1>(curEdge)->getPosition();
//        double        curEdgeVec[] = {secondHitPos[0]-firstHitPos[0],secondHitPos[1]-firstHitPos[1],secondHitPos[2]-firstHitPos[2]};
//        double        curEdgeMag   = std::sqrt(curEdgeVec[0]*curEdgeVec[0]+curEdgeVec[1]*curEdgeVec[1]+curEdgeVec[2]*curEdgeVec[2]);
        
//        std::cout << "***> edge comparison, # edges to compare: " << edgeListItr->second.size() << ", curEdge wires: " << std::get<1>(curEdge)->getWireIDs()[0].Wire << ", " << std::get<1>(curEdge)->getWireIDs()[1].Wire << ", " << std::get<1>(curEdge)->getWireIDs()[2].Wire << ", mag: " << curEdgeMag << ", sep: " << std::get<2>(curEdge) << std::endl;
        
        for(const auto& edge : edgeListItr->second)
        {
            // skip the self reference
            if (edge == curEdge)
            {
//                std::cout << "    ***> matching edges, skipping" << std::endl;
                continue;
            }

            // Don't consider edges which fold back on themselves
//            const double* nextHitPos     = std::get<1>(edge)->getPosition();
//            double        firstEdgeVec[] = {nextHitPos[0]-firstHitPos[0], nextHitPos[1]-firstHitPos[1], nextHitPos[2]-firstHitPos[2]};
//            double        nextEdgeVec[]  = {nextHitPos[0]-secondHitPos[0],nextHitPos[1]-secondHitPos[1],nextHitPos[2]-secondHitPos[2]};
//            double        firstEdgeMag   = std::sqrt(firstEdgeVec[0]*firstEdgeVec[0]+firstEdgeVec[1]*firstEdgeVec[1]+firstEdgeVec[2]*firstEdgeVec[2]);
//            double        nextEdgeMag    = std::sqrt(nextEdgeVec[0]*nextEdgeVec[0]+nextEdgeVec[1]*nextEdgeVec[1]+nextEdgeVec[2]*nextEdgeVec[2]);
//            double        curDotFirst    = (curEdgeVec[0]*firstEdgeVec[0]+curEdgeVec[1]*firstEdgeVec[1]+curEdgeVec[2]*firstEdgeVec[2]) / (firstEdgeMag);  // Projection of current onto direction of next
//            double        curDotNext     = (curEdgeVec[0]*nextEdgeVec[0]+curEdgeVec[1]*nextEdgeVec[1]+curEdgeVec[2]*nextEdgeVec[2]) / (nextEdgeMag*curEdgeMag);
            
//            std::cout << "     --> projection: " << curDotFirst << ", dot: " << curDotNext << ", wires: " << std::get<1>(edge)->getWireIDs()[0].Wire << ", " << std::get<1>(edge)->getWireIDs()[1].Wire << ", " << std::get<1>(edge)->getWireIDs()[2].Wire << ", mag: " << nextEdgeMag << ", sep: " << std::get<2>(edge) << std::endl;
            
//            if (curDotFirst < 0.) continue;
            
            double tempDepth(0.);
            
            reco::HitPairListPtr tempList = DepthFirstSearch(edge,hitToEdgeMap,tempDepth);
            
            if (tempDepth > deepest)
            {
                hitPairListPtr = tempList;
                deepest        = tempDepth;
            }
        }
    }
    
    hitPairListPtr.push_front(std::get<1>(curEdge));
    
    depth += deepest + 1./std::get<2>(curEdge);
    
    return hitPairListPtr;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
void MinSpanTreeAlg::RunDBScan(EpsPairNeighborhoodMapVec& epsPairNeighborhoodMapVec, reco::HitPairClusterMap& hitPairClusterMap) const
{
    cet::cpu_timer theClockDBScan;
    
    // Start clocks if requested
    if (m_enableMonitoring) theClockDBScan.start();
    
    // With the neighborhood built we can now form "clusters" with DBScan
    int pairClusterIdx(0);
    
    // Clear the cluster list just for something to do here...
    hitPairClusterMap.clear();
    
    // Ok, here we go!
    // We can simply iterate over the map we have just built to loop through the hits "simply"
    for(EpsPairNeighborhoodMapVec::iterator epsPairVecItr  = epsPairNeighborhoodMapVec.begin();
        epsPairVecItr != epsPairNeighborhoodMapVec.end();
        epsPairVecItr++)
    {
        // Skip the null entries (they were filtered out)
        if (!epsPairVecItr->first) continue;
        
        // If this hit has been "visited" already then skip
        if (epsPairVecItr->second.first.visited()) continue;
        
        // We are now visiting it so mark it as so
        epsPairVecItr->second.first.setVisited();
        
        // Check that density is sufficient
        if (epsPairVecItr->second.first.getCount() < m_minPairPts)
        {
            epsPairVecItr->second.first.setNoise();
        }
        else
        {
            // "Create" a new cluster and get a reference to it
            reco::HitPairListPtr& curCluster = hitPairClusterMap[pairClusterIdx++];
            
            // expand the cluster
            expandCluster(epsPairNeighborhoodMapVec, epsPairVecItr, curCluster, m_minPairPts);
        }
    }
    
    if (m_enableMonitoring)
    {
        theClockDBScan.stop();
        
        m_timeVector[RUNDBSCAN] = theClockDBScan.accumulated_real_time();
    }
    
    return;
}
    
void MinSpanTreeAlg::expandCluster(EpsPairNeighborhoodMapVec&          epsNeighborhoodMapVec,
                                   EpsPairNeighborhoodMapVec::iterator epsVecItr,
                                   reco::HitPairListPtr&               curCluster,
                                   size_t                              minPts) const
{
    // This is the main inside loop for the DBScan based clustering algorithm
    //
    // Add the current hit to the current cluster
    epsVecItr->second.first.setInCluster();
    curCluster.push_back(epsVecItr->first);
    
    // Get the list of points in this hit's epsilon neighborhood
    // Note this is a copy so we can modify locally
    EpsPairNeighborhoodList epsNeighborhoodList = epsVecItr->second.second;
    
    while(!epsNeighborhoodList.empty())
    {
        // Dereference the point so we can see in the debugger...
        const reco::ClusterHit3D* neighborPtr = epsNeighborhoodList.front();
        
        // Use that to look up the iterator in the general neighborhood map
        EpsPairNeighborhoodMapVec::iterator curPtEpsVecItr = epsNeighborhoodMapVec.begin();
        
        std::advance(curPtEpsVecItr, neighborPtr->getID());
        
        // If we've not been here before then take action...
        if (!curPtEpsVecItr->second.first.visited())
        {
            curPtEpsVecItr->second.first.setVisited();
            
            // If this epsilon neighborhood of this point is large enough then add its points to our list
            if (curPtEpsVecItr->second.first.getCount() >= minPts)
            {
                // Plan is to loop through the hits in this point's neighborhood and add them to our list provided
                // they have not already been added, or are part of a cluster, etc.
                // So, get the list of points in the neighborhood
                EpsPairNeighborhoodList& locList = curPtEpsVecItr->second.second;
                
                // And loop through them...
                for(EpsPairNeighborhoodList::iterator hitItr = locList.begin(); hitItr != locList.end(); hitItr++)
                {
                    epsNeighborhoodList.push_back(*hitItr);
                }
            }
        }
        
        // If the point is not yet in a cluster then we now add
        if (!curPtEpsVecItr->second.first.inCluster())
        {
            curPtEpsVecItr->second.first.setInCluster();
            curCluster.push_back(curPtEpsVecItr->first);
        }
        
        epsNeighborhoodList.pop_front();
    }
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
KdTreeNode MinSpanTreeAlg::BuildKdTree(const reco::HitPairList& hitPairList, KdTreeNodeList& kdTreeNodeContainer) const
{
    
    // The first task is to build the kd tree
    cet::cpu_timer theClockBuildNeighborhood;
    
    if (m_enableMonitoring) theClockBuildNeighborhood.start();
    
    // Start by constructing a kd tree to locate the 3D hits
    Hit3DVec hit3DVec;
    
    hit3DVec.reserve(hitPairList.size());
    
    Hit3DVec::iterator first = hit3DVec.begin();
    Hit3DVec::iterator last  = hit3DVec.end();
    
    for(const auto& hitPtr : hitPairList) hit3DVec.emplace_back(hitPtr.get());
    
    KdTreeNode topNode = BuildKdTree(hit3DVec.begin(), hit3DVec.end(), kdTreeNodeContainer);
    
    if (m_enableMonitoring)
    {
        theClockBuildNeighborhood.stop();
        m_timeVector[BUILDHITTOHITMAP] = theClockBuildNeighborhood.accumulated_real_time();
    }
    
    return topNode;
}
    
KdTreeNode& MinSpanTreeAlg::BuildKdTree(Hit3DVec::iterator first,
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
        if (refHit == node.getClusterHit3D()) bestDist = 1.0;
        // This is the tight constraint on the hits
        else if (bestDist < std::numeric_limits<double>::max() && !(node.getClusterHit3D()->getStatusBits() & reco::ClusterHit3D::CLUSTERATTACHED) && consistentPairs(refHit, node.getClusterHit3D(), hitSeparation, wireDeltas))
        {
            candPairVec.emplace_back(CandPair(hitSeparation,node.getClusterHit3D()));
            
            bestDist = std::max(0.65,std::min(bestDist,hitSeparation));
            
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
    
size_t MinSpanTreeAlg::BuildNeighborhoodMap(const reco::HitPairList&   hitPairList,
                                            EpsPairNeighborhoodMapVec& epsPairNeighborhoodMapVec) const
{
    /**
     *  @brief build out the epsilon neighborhood map to be used by DBScan
     */
    
    cet::cpu_timer theClockBuildNeighborhood;
    size_t         consistentPairsCnt(0);
    size_t         pairsChecked(0);
    
    if (m_enableMonitoring) theClockBuildNeighborhood.start();

    // Start by constructing a kd tree to locate the 3D hits
    KdTreeNodeList kdTreeNodeContainer;
    Hit3DVec       hit3DVec;
    
    hit3DVec.reserve(hitPairList.size());
    
    for(const auto& hitPtr : hitPairList) hit3DVec.emplace_back(hitPtr.get());
    
    KdTreeNode topNode = BuildKdTree(hit3DVec.begin(), hit3DVec.end(), kdTreeNodeContainer);

    // Make a container to hold candidate consistent pairs
    CandPairVec candPairVec;
    
    //**********************************************************************************
    // Given the list of pairs of hits which are consistent with each other, build out the
    // epsilon neighbor maps
    // Now we go through the pair list and basically repeat the same exercise as above
    // The following assumes that the HitPairList is ordered
    // a) in increasing Z for hits which are not on the "same W wire",
    // b) in increasing U (Y) for hits on the same W wire
    
    for (reco::HitPairList::const_iterator pairItrO = hitPairList.begin(); pairItrO != hitPairList.end(); pairItrO++)
    {
        const reco::ClusterHit3D* hitPairO   = (*pairItrO).get();
        const size_t              hitPairOID = hitPairO->getID();
        
        // Need to initialize the "first" part of the vector pseudo map
        epsPairNeighborhoodMapVec[hitPairOID].first = hitPairO;
        
        // Get reference to the list for this hit so we don't look it up inside the loop
        EpsPairNeighborhoodPair& hitPairOPair(epsPairNeighborhoodMapVec[hitPairOID].second);
        
        // Clear the container
        candPairVec.clear();
        
        double bestDistance(std::numeric_limits<double>::max());

        FindNearestNeighbors(hitPairO, topNode, candPairVec, bestDistance);
        
        if (candPairVec.empty()) continue;
        
        // Sort by the maximum delta wires
        std::sort(candPairVec.begin(),candPairVec.end(),[](const CandPair& left, const CandPair& right){return left.first < right.first;});
        
        size_t minNeighbors = std::min(candPairVec.size()-1,3*m_minPairPts);  // within range of our minimum points...
        double maxDistance  = 1.2 * candPairVec.at(minNeighbors).first;
        
        for(const auto& candPair : candPairVec)
        {
            if (candPair.first > maxDistance) break;
            
            const reco::ClusterHit3D* hitPairI(candPair.second);
            
            hitPairOPair.first.incrementCount();
            hitPairOPair.second.emplace_back(hitPairI);
            
            epsPairNeighborhoodMapVec[hitPairI->getID()].second.first.incrementCount();
            epsPairNeighborhoodMapVec[hitPairI->getID()].second.second.emplace_back(hitPairO);
            
            consistentPairsCnt++;
        }
    }
    
    if (m_enableMonitoring)
    {
        theClockBuildNeighborhood.stop();
        m_timeVector[BUILDHITTOHITMAP] = theClockBuildNeighborhood.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << "Consistent pairs: " << consistentPairsCnt << " of " << pairsChecked << " checked." << std::endl;
    
    return consistentPairsCnt;
}
    
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
    
    // Loose constraint to weed out the obviously bad combinations
    if (pairOverlap > 0.1)
    {
        double pairDistVec[] = {pair1->getPosition()[0] - pair2->getPosition()[0], pair1->getPosition()[1] - pair2->getPosition()[1], pair1->getPosition()[2] - pair2->getPosition()[2]};
        
        hitSeparation = std::sqrt(pairDistVec[0]*pairDistVec[0] + pairDistVec[1]*pairDistVec[1] + pairDistVec[2]*pairDistVec[2]);
        
        size_t hitCount(0);
        
        // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
        for(size_t idx = 0; idx < 3; idx++)
        {
            wireDeltas[idx] = std::abs(int(pair1->getWireIDs()[idx].Wire) - int(pair2->getWireIDs()[idx].Wire));
            
            if (pair1->getHits()[idx]) hitCount++;
            if (pair2->getHits()[idx]) hitCount++;
        }
        
        // put wire deltas in order...
        std::sort(wireDeltas, wireDeltas + 3);
        
        if (wireDeltas[0] < 2 && wireDeltas[1] < 2 && wireDeltas[2] < 3)
        {
            double overlapFraction = 0.5 * pairOverlap / std::min(pair1Width,pair2Width);

            for(size_t idx = 0; idx < 6-hitCount; idx++) hitSeparation *= 1.1;

            hitSeparation /= overlapFraction;
            
            hitSeparation = std::max(0.0001,hitSeparation);
            
//            if (hitSeparation < 0.001) std::cout << "*********> small separation, pairs: " << pair1 << ", " << pair2 << std::endl;
            
            return true;
        }
    }
    
    return false;
}
    
struct HitPairClusterOrder
{
    bool operator()(const reco::HitPairClusterMap::iterator& left, const reco::HitPairClusterMap::iterator& right)
    {
        // Watch out for the case where two clusters can have the same number of hits!
        if (left->second.size() == right->second.size())
            return left->first < right->first;
        
        return left->second.size() > right->second.size();
    }
};
    
void MinSpanTreeAlg::BuildClusterInfo(reco::HitPairClusterMap& hitPairClusterMap, reco::ClusterParametersList& clusterParametersList) const
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
    if (!hitPairClusterMap.empty())
    {
        // We want to order our clusters on by largest (most number hits) to smallest. So, we'll loop through the clusters,
        // weeding out the unwanted ones and keep track of things in a set of "good" clusters which we'll order
        // by cluster size.
        std::set<reco::HitPairClusterMap::iterator, HitPairClusterOrder> hitPairClusterSet;
        
        // Loop over the "Clusters" in our map where this loop serves a double purpose
        // In the first we are weeding out clusters which fall below what we think is the minimum number of hits
        // More importantly, we are transferring the cluster ownership to a set so it can order the clusters by
        // number of associated hits
        for (reco::HitPairClusterMap::iterator mapItr = hitPairClusterMap.begin(); mapItr != hitPairClusterMap.end(); mapItr++)
        {
            // Weed out the little people
            if (mapItr->second.size() < m_clusterMinHits) continue;
            
            // Add to the set
            hitPairClusterSet.insert(mapItr);
        }
        
        // What remains is an order set of clusters, largest first
        // Now go through and obtain cluster parameters
        for(std::set<reco::HitPairClusterMap::iterator, HitPairClusterOrder>::iterator setItr = hitPairClusterSet.begin(); setItr != hitPairClusterSet.end(); setItr++)
        {
            // Recover original map iterator
            reco::HitPairClusterMap::iterator hitPairClusterMapItr = *setItr;
            
            // Create a new cluster params object in the vector
            clusterParametersList.push_back(reco::ClusterParameters(hitPairClusterMapItr));
            
            // Can we get a reference to what we just created?
            reco::ClusterParameters& clusterParams = clusterParametersList.back();
            
            // Do the actual work of filling the parameters
            FillClusterParams(clusterParams, m_clusterMinUniqueFraction, m_clusterMaxLostFraction);
            
            // If this cluster is rejected then the parameters will be empty
            if (clusterParams.m_clusterParams.empty() || !clusterParams.m_fullPCA.getSvdOK())
            {
                clusterParametersList.pop_back();
            }
        }
    }
    
    if (m_enableMonitoring)
    {
        theClockBuildClusterInfo.stop();
        
        m_timeVector[BUILDCLUSTERINFO] = theClockBuildClusterInfo.accumulated_real_time();
    }
    
    return;
}

void MinSpanTreeAlg::FillClusterParams(reco::ClusterParameters& clusterParams, double minUniqueFrac, double maxLostFrac) const
{
    /**
     *  @brief Given a list of hits fill out the remaining parameters for this cluster and evaluate the
     *         candidate's worthiness to achieve stardom in the event display
     */
    
    // Recover the HitPairListPtr from the input clusterParams (which will be the
    // only thing that has been provided)
    reco::HitPairListPtr& hitPairVector = clusterParams.m_hitPairListPtr;
    
    // To be sure, we should clear the other data members
    clusterParams.m_clusterParams.clear();
    clusterParams.m_fullPCA = reco::PrincipalComponents();
    
    // See if we can avoid duplicates by temporarily transferring to a set
    //std::set<const reco::ClusterHit2D*> hitSet;
    std::vector<const reco::ClusterHit2D*> hitSet;
    
    size_t nTotalHits[]  = {0,0,0};
    size_t nUniqueHits[] = {0,0,0};
    size_t nLostHits[]   = {0,0,0};
    
    // Create a list to hold 3D hits which are already in use (criteria below)
    reco::HitPairListPtr usedHitPairList;
    
    // First loop through the 3D hits
    // The goal of this loop is to build a set of unique hits from the hit pairs (which may contain many
    // ambiguous duplicate combinations).
    // The secondary goal is to remove 3D hits marked by hit arbitration to be tossed
    for(const auto& hit3D : hitPairVector)
    {
        size_t nHits2D(0);
        size_t nHitsUsed[] = {0,0,0};
        
        // loop over the hits in this 3D Cluster hit
        for(const auto& hit2D : hit3D->getHits())
        {
            if (!hit2D) continue;
            size_t view = hit2D->getHit().View();
            if (hit2D->getStatusBits() & reco::ClusterHit2D::USED) nHitsUsed[view]++;
            else                                                   nUniqueHits[view]++;
            nTotalHits[view]++;
            nHits2D++;
        }
        
        size_t nHitsAlreadyUsed = std::accumulate(nHitsUsed,nHitsUsed+3,0);
        
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
    
    // If we have something left then at this point we make one more check
    // This check is intended to weed out clusters made from isolated groups of ambiguous hits which
    // really belong to a larger cluster
    if (numUniqueHits > 3)
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
            
            size_t nViewsWithHits    = (clusterParams.m_clusterParams[geo::kU].m_hitVector.size() > 0 ? 1 : 0)
            + (clusterParams.m_clusterParams[geo::kV].m_hitVector.size() > 0 ? 1 : 0)
            + (clusterParams.m_clusterParams[geo::kW].m_hitVector.size() > 0 ? 1 : 0);
            size_t nViewsWithMinHits = (clusterParams.m_clusterParams[geo::kU].m_hitVector.size() > 2 ? 1 : 0)
            + (clusterParams.m_clusterParams[geo::kV].m_hitVector.size() > 2 ? 1 : 0)
            + (clusterParams.m_clusterParams[geo::kW].m_hitVector.size() > 2 ? 1 : 0);
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
                m_pcaAlg.PCAAnalysis_3D(clusterParams.m_hitPairListPtr, clusterParams.m_fullPCA);
                
                // Must have a valid pca
                if (clusterParams.m_fullPCA.getSvdOK())
                {
                    // If any hits were thrown away, see if we can rescue them
                    if (!usedHitPairList.empty())
                    {
                        double maxDoca = 2. * sqrt(clusterParams.m_fullPCA.getEigenValues()[1]);
                        
                        if (maxDoca < 5.)
                        {
                            size_t curHitVectorSize = hitPairVector.size();
                            
                            m_pcaAlg.PCAAnalysis_calc3DDocas(usedHitPairList, clusterParams.m_fullPCA);
                            
                            for(const auto& hit3D : usedHitPairList)
                                if (hit3D->getDocaToAxis() < maxDoca) hitPairVector.push_back(hit3D);
                            
                            if (hitPairVector.size() > curHitVectorSize)
                                m_pcaAlg.PCAAnalysis_3D(clusterParams.m_hitPairListPtr, clusterParams.m_fullPCA);
                        }
                    }
                    
                    // Set the skeleton PCA to make sure it has some value
                    clusterParams.m_skeletonPCA = clusterParams.m_fullPCA;
                }
            }
        }
    }
    
    return;
}
    

} // namespace lar_cluster3d

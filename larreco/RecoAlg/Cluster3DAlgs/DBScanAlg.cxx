/**
 *  @file   Cluster3D_module.cc
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/DBScanAlg.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

DBScanAlg::DBScanAlg(fhicl::ParameterSet const &pset) :
    m_channelFilter(&art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider()),
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DBScanAlg::~DBScanAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DBScanAlg::reconfigure(fhicl::ParameterSet const &pset)
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
    
void DBScanAlg::expandCluster(EpsPairNeighborhoodMapVec&          epsNeighborhoodMapVec,
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
    
void DBScanAlg::ClusterHitsDBScan(reco::HitPairList&           hitPairList,
                                  reco::HitPairClusterMap&     hitPairClusterMap,
                                  reco::ClusterParametersList& clusterParametersList)
{
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */
    cet::cpu_timer theClockDBScan;
    
    m_timeVector.resize(NUMTIMEVALUES, 0.);
    
    // The container of pairs and those in each pair's epsilon neighborhood
    EpsPairNeighborhoodMapVec epsPairNeighborhoodMapVec;
    
    epsPairNeighborhoodMapVec.resize(hitPairList.size(), EpsPairNeighborhoodMapPair(0,EpsPairNeighborhoodPair()));  //<-- initialize too!
    
    // DBScan is driven of its "epsilon neighborhood". Computing adjacency within DBScan can be time
    // consuming so the idea is the prebuild the adjaceny map and then run DBScan.
    // The following call does this work
    BuildNeighborhoodMap(hitPairList, epsPairNeighborhoodMapVec);
    
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
    
    // Initial clustering is done, now trim the list and get output parameters
    BuildClusterInfo(hitPairClusterMap, clusterParametersList);
    
    mf::LogDebug("Cluster3D") << ">>>>> DBScan done, found " << hitPairClusterMap.size() << " clusters" << std::endl;
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
size_t DBScanAlg::BuildNeighborhoodMap(reco::HitPairList& hitPairList, EpsPairNeighborhoodMapVec& epsPairNeighborhoodMapVec) const
{
    /**
     *  @brief build out the epsilon neighborhood map to be used by DBScan
     */
    
    cet::cpu_timer theClockBuildNeighborhood;
    size_t         consistentPairsCnt(0);
    size_t         pairsChecked(0);
    double         hitSeparation;
    
    if (m_enableMonitoring) theClockBuildNeighborhood.start();
    
    // Make a container to hold candidate consistent pairs
    using CandPair = std::pair<double,const reco::ClusterHit3D*>;
    std::vector<CandPair> candPairVec;
    
    int wireDeltas[] = {0,0,0};
    
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

        reco::HitPairList::const_iterator pairItrI = pairItrO;
        
        while (++pairItrI != hitPairList.end())
        {
            const reco::ClusterHit3D* hitPairI = (*pairItrI).get();
            
            // Keep count...
            pairsChecked++;
            
            // This is the tight constraint on the hits
            if (consistentPairs(hitPairO, hitPairI, hitSeparation, wireDeltas))
            {
                candPairVec.emplace_back(CandPair(hitSeparation,hitPairI));
            }
            
            // Is there a loop termination condition? Hits are sorted in pulse start time order
            if (hitPairO->getAvePeakTime() + 3.*hitPairO->getSigmaPeakTime() < hitPairI->getAvePeakTime() - 3.*hitPairI->getSigmaPeakTime()) break;
        }
        
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
    
bool DBScanAlg::consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const
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
        
        if ((wireDeltas[0] < 2 && wireDeltas[1] < 2) || (wireDeltas[0] < 2 && wireDeltas[1] < 3 && wireDeltas[2] < 8))
        {
            double overlapFraction = 0.5 * pairOverlap / std::min(pair1Width,pair2Width);
            
            for(size_t idx = 0; idx < 6-hitCount; idx++) hitSeparation *= 1.1;
            
            hitSeparation /= overlapFraction;
            
            return true;
        }
    }
    
    return false;
}
    
bool DBScanAlg::consistentPairsOrig(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const
{
    // Strategy: We consider comparing "hit pairs" which may consist of 2 or 3 actual hits.
    //           Also, if only pairs, they can be U-V, U-W or V-W so we can't assume which views we have
    //           So do a simplified comparison:
    //           1) compare the pair times and require "overlap" (in the sense of hit pair production)
    //           2) look at distance between pairs in each of the wire directions
    
    double pair1PeakTime = pair1->getAvePeakTime();
    double pair1Width    = 3.*pair1->getSigmaPeakTime();
    double pair2PeakTime = pair2->getAvePeakTime();
    double pair2Width    = 3.*pair2->getSigmaPeakTime();
    
    double maxUpper      = std::min(pair1PeakTime+pair1Width,pair2PeakTime+pair2Width);
    double minLower      = std::max(pair1PeakTime-pair1Width,pair2PeakTime-pair2Width);
    double pairOverlap   = maxUpper - minLower;
    
    // Loose constraint to weed out the obviously bad combinations
    if (pairOverlap > 0.1)
    {
        int numWithOverlap(0);
        int nIdenticalHits(0);
        
        const int maxDeltaWires(2);
        const int maxIdenticalHits(2);
        const int minNumWithOverlap(2);
        
        // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
        for(size_t idx = 0; idx < 3; idx++)
        {
            wireDeltas[idx] = std::abs(int(pair1->getWireIDs()[idx].Wire) - int(pair2->getWireIDs()[idx].Wire));
            
            if (wireDeltas[idx] > maxDeltaWires) break;
            
            if (!pair1->getHits().at(idx) || !pair2->getHits().at(idx)) continue;
            
            const recob::Hit& hit1 = pair1->getHits().at(idx)->getHit();
            const recob::Hit& hit2 = pair2->getHits().at(idx)->getHit();
            
            if (&hit1 == &hit2) nIdenticalHits++;
            
            if (nIdenticalHits > maxIdenticalHits) break;
            
            double hit1PeakTime = hit1.PeakTime();
            double hit1Width    = m_pairSigmaPeakTime * hit1.RMS();
            double hit2PeakTime = hit2.PeakTime();
            double hit2Width    = m_pairSigmaPeakTime * hit2.RMS();
            
            // Allow the hit separation to grow a bit if there are missing wires in between
            if (wireDeltas[idx] > 0)
            {
                double scaleFactor = std::min(1.5 + 0.5 * double(wireDeltas[idx]),3.);
                
                hit1Width = scaleFactor * hit1.RMS();
                hit2Width = scaleFactor * hit2.RMS();
            }
            
            double hitMaxUpper = std::min(hit1PeakTime+hit1Width,hit2PeakTime+hit2Width);
            double hitMinLower = std::max(hit1PeakTime-hit1Width,hit2PeakTime-hit2Width);
            double hitOverlap  = hitMaxUpper - hitMinLower;
            
            if (hitOverlap > 0.1) numWithOverlap++;
        }
        
        // Too many identical hits?
        if (nIdenticalHits > maxIdenticalHits) return false;
        
        // check to make sure the timing is consistent
        if (numWithOverlap < minNumWithOverlap) return false;
        
        std::sort(wireDeltas, wireDeltas + 3);
        
        // At least one view is required to have neighboring wires
        if (wireDeltas[0] > 1 || wireDeltas[1] > 1 || wireDeltas[2] > maxDeltaWires) return false;
        
        return true;
    }
    
    return false;
}
    
bool DBScanAlg::consistentPairsTest(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const
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
        size_t nWiresDiff(0);
        size_t sumDeltaWire(0);
        size_t hitCount(0);
        
        // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
        for(size_t idx = 0; idx < 3; idx++)
        {
            wireDeltas[idx] = std::abs(int(pair1->getWireIDs()[idx].Wire) - int(pair2->getWireIDs()[idx].Wire));
            
            sumDeltaWire += wireDeltas[idx];

            if (wireDeltas[idx] > 0)   nWiresDiff++;
            if (pair1->getHits()[idx]) hitCount++;
            if (pair2->getHits()[idx]) hitCount++;
        }
        
        // put wire deltas in order...
        std::sort(wireDeltas, wireDeltas + 3);
        
//        if ((nWiresDiff == 1 && sumDeltaWire == 1) || (nWiresDiff == 2 && sumDeltaWire > 1 && sumDeltaWire < 4) || (nWiresDiff == 3 && sumDeltaWire == 3))
        if (wireDeltas[0] < 2 && wireDeltas[1] < 2)
        {
            hitSeparation = 0.3 * double(sumDeltaWire);
            
            hitSeparation = double(7 - hitCount) * hitSeparation;
            
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
    
void DBScanAlg::BuildClusterInfo(reco::HitPairClusterMap& hitPairClusterMap, reco::ClusterParametersList& clusterParametersList) const
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

void DBScanAlg::FillClusterParams(reco::ClusterParameters& clusterParams, double minUniqueFrac, double maxLostFrac) const
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
    std::set<const reco::ClusterHit2D*> hitSet;
    
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
                if (hit3D->getHits()[idx]) hitSet.insert(hit3D->getHits()[idx]);
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

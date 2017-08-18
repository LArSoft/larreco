/**
 *  @file   Cluster3D_module.cc
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/ClusterParamsBuilder.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larcore/Geometry/Geometry.h"
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
/**
 *  @brief a utility class for keeping track of the state of a hit for DBScan
 */
class DBScanParams
{
public:
    DBScanParams() : m_visited(false), m_noise(false), m_inCluster(false), m_count(0) {}
    
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
 *  @brief  DBScanAlg class definiton
 */
class DBScanAlg : virtual public IClusterAlg
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit DBScanAlg(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief  Destructor
     */
    ~DBScanAlg();
    
    void configure(const fhicl::ParameterSet&) override;
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void Cluster3DHits(reco::HitPairList&           hitPairList,
                       reco::ClusterParametersList& clusterParametersList) const override;
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    double getTimeToExecute(IClusterAlg::TimeValues index) const override {return m_timeVector.at(index);}
    
private:
    
    /**
     *  @brief The bigger question: are two pairs of hits consistent?
     */
    bool consistentPairsOrig(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const;
    bool consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const;
    
    bool consistentPairsTest(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const;
    
    typedef std::list<const reco::ClusterHit3D*>                           EpsPairNeighborhoodList;
    typedef std::pair<DBScanParams, EpsPairNeighborhoodList >              EpsPairNeighborhoodPair;
    typedef std::map<const reco::ClusterHit3D*, EpsPairNeighborhoodPair >  EpsPairNeighborhoodMap;
    typedef std::pair<const reco::ClusterHit3D*, EpsPairNeighborhoodPair > EpsPairNeighborhoodMapPair;
    typedef std::vector<EpsPairNeighborhoodMapPair >                       EpsPairNeighborhoodMapVec;
    
    /**
     *  @brief the main routine for DBScan
     */
    void expandCluster(EpsPairNeighborhoodMapVec&          nMap,
                       EpsPairNeighborhoodMapVec::iterator vecItr,
                       reco::HitPairListPtr&               cluster,
                       size_t                              minPts) const;
    
    /**
     *  @brief Given an input HitPairList, build out the map of nearest neighbors
     */
    size_t BuildNeighborhoodMap(reco::HitPairList& hitPairList, EpsPairNeighborhoodMapVec& epsPairNeighborhoodMapVec) const;
    
    /**
     *  @brief define data structure for keeping track of channel status
     */
    
    using ChannelStatusVec        = std::vector<size_t>;
    using ChannelStatusByPlaneVec = std::vector<ChannelStatusVec>;
    
    /**
     *  @brief Data members to follow
     */
    size_t                               m_minPairPts;
    double                               m_pairSigmaPeakTime;
    
    bool                                 m_enableMonitoring;      ///<
    double                               m_wirePitch[3];
    mutable std::vector<float>           m_timeVector;            ///<
    std::vector<std::vector<double>>     m_wireDir;               ///<
    std::vector<std::vector<double>>     m_wireNormal;            ///<
    
    geo::Geometry*                       m_geometry;              //< pointer to the Geometry service
    
    ClusterParamsBuilder                 m_clusterBuilder;        // Common cluster builder tool
};

DBScanAlg::DBScanAlg(fhicl::ParameterSet const &pset) :
    m_clusterBuilder(pset.get<fhicl::ParameterSet>("ClusterParamsBuilder"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DBScanAlg::~DBScanAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DBScanAlg::configure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring         = pset.get<bool>  ("EnableMonitoring",  true  );
    m_minPairPts               = pset.get<size_t>("MinPairPts",                2     );
    m_pairSigmaPeakTime        = pset.get<double>("PairSigmaPeakTime",         3.    );
    
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
    
void DBScanAlg::Cluster3DHits(reco::HitPairList&           hitPairList,
                              reco::ClusterParametersList& clusterParametersList) const
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
            clusterParametersList.push_back(reco::ClusterParameters());
            
            reco::HitPairListPtr& curCluster = clusterParametersList.back().getHitPairListPtr();
            
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
    cet::cpu_timer theClockBuildClusters;
    
    // Start clocks if requested
    if (m_enableMonitoring) theClockBuildClusters.start();
    
    m_clusterBuilder.BuildClusterInfo(clusterParametersList);
    
    if (m_enableMonitoring)
    {
        theClockBuildClusters.stop();
        
        m_timeVector[BUILDCLUSTERINFO] = theClockBuildClusters.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> DBScan done, found " << clusterParametersList.size() << " clusters" << std::endl;
    
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
    
DEFINE_ART_CLASS_TOOL(DBScanAlg)
} // namespace lar_cluster3d

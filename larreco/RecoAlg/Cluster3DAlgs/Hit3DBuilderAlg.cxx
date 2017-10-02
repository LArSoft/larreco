/**
 *  @file   Cluster3D_module.cc
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/Hit3DBuilderAlg.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

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

Hit3DBuilderAlg::Hit3DBuilderAlg(fhicl::ParameterSet const &pset) :
    m_channelFilter(&art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider())
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

Hit3DBuilderAlg::~Hit3DBuilderAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Hit3DBuilderAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring = pset.get<bool> ("EnableMonitoring",    true);
    m_numSigmaPeakTime = pset.get<float>("NumSigmaPeakTime",    3.  );
    m_hitWidthSclFctr  = pset.get<float>("HitWidthScaleFactor", 6.  );
    m_deltaPeakTimeSig = pset.get<float>("DeltaPeakTimeSig",    1.7 );
    
    art::ServiceHandle<geo::Geometry> geometry;
    
    m_geometry = &*geometry;
    
    m_wirePitch[0] = m_geometry->WirePitch(0);
    m_wirePitch[1] = m_geometry->WirePitch(1);
    m_wirePitch[2] = m_geometry->WirePitch(2);
    
    m_buildTime = 0.;
}
    
void Hit3DBuilderAlg::BuildChannelStatusVec(PlaneToWireToHitSetMap& planeToWireToHitSetMap)
{
    // This is called each event, clear out the previous version and start over
    if (!m_channelStatus.empty()) m_channelStatus.clear();
    
    m_channelStatus.resize(m_geometry->Nplanes());

    // Loop through views/planes to set the wire length vectors
    for(size_t idx = 0; idx < m_channelStatus.size(); idx++)
    {
        m_channelStatus.at(idx) = ChannelStatusVec(m_geometry->Nwires(idx), 5);
    }
    
    // Loop through the channels and mark those that are "bad"
    for(size_t channel = 0; channel < m_geometry->Nchannels(); channel++)
    {
        if( !m_channelFilter->IsGood(channel))
        {
            std::vector<geo::WireID>                wireIDVec = m_geometry->ChannelToWire(channel);
            geo::WireID                             wireID    = wireIDVec[0];
            lariov::ChannelStatusProvider::Status_t chanStat  = m_channelFilter->Status(channel);
            
            m_channelStatus[wireID.Plane][wireID.Wire] = chanStat;
        }
    }
    
    // add quiet wires in U plane for microboone (this will done "correctly" in near term)
//    PlaneToWireToHitSetMap::iterator plane0HitItr = planeToWireToHitSetMap.find(geo::PlaneID(0,0,0));
    
//    if (plane0HitItr != planeToWireToHitSetMap.end())
//    {
////        WireToHitSetMap& wireToHitSetMap = uPlaneHitItr->second;
    
//        for(size_t idx = 2016; idx < 2096; idx++)  m_channelStatus[0][idx] = 3;
//        for(size_t idx = 2192; idx < 2304; idx++)  m_channelStatus[0][idx] = 3;
//        for(size_t idx = 2352; idx < 2384; idx++)  m_channelStatus[0][idx] = 3;
//        //for(size_t idx = 2016; idx < 2096; idx++) if (wireToHitSetMap.find(idx) == wireToHitSetMap.end()) m_channelStatus[0][idx] = 3;
//        //for(size_t idx = 2192; idx < 2304; idx++) if (wireToHitSetMap.find(idx) == wireToHitSetMap.end()) m_channelStatus[0][idx] = 3;
//        //for(size_t idx = 2352; idx < 2384; idx++) if (wireToHitSetMap.find(idx) == wireToHitSetMap.end()) m_channelStatus[0][idx] = 3;
////      for(size_t idx = 2016; idx < 2384; idx++) m_channelStatus[0][idx] = 3;
//    }
    
    return;
}

    
bool SetPeakHitPairIteratorOrder(const reco::HitPairList::iterator& left, const reco::HitPairList::iterator& right)
{
    return (*left)->getAvePeakTime() < (*right)->getAvePeakTime();
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
    
void Hit3DBuilderAlg::BuildHit3D(PlaneToHitVectorMap& planeToHitVectorMap, PlaneToWireToHitSetMap& planeToWireToHitSetMap, reco::HitPairList& hitPairList)
{
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */
    cet::cpu_timer theClockMakeHits;
    
    if (m_enableMonitoring) theClockMakeHits.start();
    
    // The first task is to take the lists of input 2D hits (a map of view to sorted lists of 2D hits)
    // and then to build a list of 3D hits to be used in downstream processing
    BuildChannelStatusVec(planeToWireToHitSetMap);
    
    size_t numHitPairs = BuildHitPairMap(planeToHitVectorMap, hitPairList);
    
    if (m_enableMonitoring)
    {
        theClockMakeHits.stop();
        
        m_buildTime = theClockMakeHits.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> 3D hit building done, found " << numHitPairs << " 3D Hits" << std::endl;
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
bool SetPositionOrder(const std::unique_ptr<reco::ClusterHit3D>& left, const std::unique_ptr<reco::ClusterHit3D>& right)
{
    // The positions in the Y-Z plane are quantized so take advantage of that for ordering
    // First check that we are in the same "bin" in the z direction
    if (left->getHits().back()->getHit().WireID().Wire == right->getHits().back()->getHit().WireID().Wire) // These hits are "on the same w wire"
    {
        // We can use the U wires as a proxy for ordering hits in increasing Y
        // where we remember that as Y increases the U wire number decreases
        if (left->getHits().front()->getHit().WireID().Wire == right->getHits().front()->getHit().WireID().Wire)
        {
            // In the time direction we are not quantized at a large enough level to check
            return left->getX() < right->getX();
        }
        
        return left->getHits().front()->getHit().WireID().Wire > right->getHits().front()->getHit().WireID().Wire;
    }
    
    return left->getHits().back()->getHit().WireID().Wire < right->getHits().back()->getHit().WireID().Wire;
}
    
bool SetHitStartTimeOrder(const reco::ClusterHit2D* left, const reco::ClusterHit2D* right)
{
    // Sort by "modified start time" of pulse
    return left->getHit().PeakTime() - left->getHit().RMS() < right->getHit().PeakTime() - right->getHit().RMS();
}
    
using  HitVectorItrPair = std::pair<HitVector::iterator,HitVector::iterator>;
    
class SetStartTimeOrder
{
public:
    SetStartTimeOrder()             : m_numRMS(1.)     {}
    SetStartTimeOrder(float numRMS) : m_numRMS(numRMS) {}
    
    bool operator()(const HitVectorItrPair& left, const HitVectorItrPair& right) const
    {
        // Protect against possible issue?
        if (left.first != left.second && right.first != right.second)
        {
            // Sort by "modified start time" of pulse
            return (*left.first)->getTimeTicks() - m_numRMS*(*left.first)->getHit().RMS() < (*right.first)->getTimeTicks() - m_numRMS*(*right.first)->getHit().RMS();
        }
    
        return left.first != left.second;
    }
    
private:
    float m_numRMS;
};
    
bool SetPairStartTimeOrder(const std::unique_ptr<reco::ClusterHit3D>& left, const std::unique_ptr<reco::ClusterHit3D>& right)
{
    // Sort by "modified start time" of pulse
    return left->getAvePeakTime() - left->getSigmaPeakTime() < right->getAvePeakTime() - right->getSigmaPeakTime();
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
size_t Hit3DBuilderAlg::BuildHitPairMap(PlaneToHitVectorMap& planeToHitVectorMap, reco::HitPairList& hitPairList) const
{
    /**
     *  @brief Given input 2D hits, build out the lists of possible 3D hits
     *
     *         The current strategy: ideally all 3D hits would be comprised of a triplet of 2D hits, one from each view
     *         However, we have concern that, in particular, the v-plane may have some inefficiency which we have to be
     *         be prepared to deal with. The idea, then, is to first make the association of hits in the U and W planes
     *         and then look for the match in the V plane. In the event we don't find the match in the V plane then we
     *         will evaluate the situation and in some instances keep the U-W pairs in order to keep efficiency high.
     */
    size_t totalNumHits(0);
    size_t hitPairCntr(0);
    
    size_t nTriplets(0);
    size_t nDeadChanHits(0);
    
    // Set up to loop over cryostats and tpcs...
    for(size_t cryoIdx = 0; cryoIdx < m_geometry->Ncryostats(); cryoIdx++)
    {
        for(size_t tpcIdx = 0; tpcIdx < m_geometry->NTPC(); tpcIdx++)
        {
            PlaneToHitVectorMap::iterator mapItr0 = planeToHitVectorMap.find(geo::PlaneID(cryoIdx,tpcIdx,0));
            PlaneToHitVectorMap::iterator mapItr1 = planeToHitVectorMap.find(geo::PlaneID(cryoIdx,tpcIdx,1));
            PlaneToHitVectorMap::iterator mapItr2 = planeToHitVectorMap.find(geo::PlaneID(cryoIdx,tpcIdx,2));
    
            size_t nPlanesWithHits = (mapItr0 != planeToHitVectorMap.end() ? 1 : 0)
                                   + (mapItr1 != planeToHitVectorMap.end() ? 1 : 0)
                                   + (mapItr2 != planeToHitVectorMap.end() ? 1 : 0);
    
            if (nPlanesWithHits < 2) continue;
    
            HitVector& hitVector0 = mapItr0->second;
            HitVector& hitVector1 = mapItr1->second;
            HitVector& hitVector2 = mapItr2->second;
    
            // We are going to resort the hits into "start time" order...
            std::sort(hitVector0.begin(), hitVector0.end(), SetHitStartTimeOrder);
            std::sort(hitVector1.begin(), hitVector1.end(), SetHitStartTimeOrder);
            std::sort(hitVector2.begin(), hitVector2.end(), SetHitStartTimeOrder);
    
            PlaneHitVectorItrPairVec hitItrVec
                     = {HitVectorItrPair(hitVector0.begin(),hitVector0.end()),
                        HitVectorItrPair(hitVector1.begin(),hitVector1.end()),
                        HitVectorItrPair(hitVector2.begin(),hitVector2.end())
                        };
    
            totalNumHits += BuildHitPairMapByTPC(hitItrVec, hitPairList);
        }
    }
    
    // Return the hit pair list but sorted by z and y positions (faster traversal in next steps)
    hitPairList.sort(SetPairStartTimeOrder);
    
    // Where are we?
    mf::LogDebug("Cluster3D") << "Total number hits: " << totalNumHits << std::endl;
    mf::LogDebug("Cluster3D") << "Created a total of " << hitPairList.size() << " hit pairs, counted: " << hitPairCntr << std::endl;
    mf::LogDebug("Cluster3D") << "-- Triplets: " << nTriplets << ", dead channel pairs: " << nDeadChanHits << std::endl;
    
    return hitPairList.size();
}
    
size_t Hit3DBuilderAlg::BuildHitPairMapByTPC(PlaneHitVectorItrPairVec& hitItrVec, reco::HitPairList& hitPairList) const
{
    /**
     *  @brief Given input 2D hits, build out the lists of possible 3D hits
     *
     *         The current strategy: ideally all 3D hits would be comprised of a triplet of 2D hits, one from each view
     *         However, we have concern that, in particular, the v-plane may have some inefficiency which we have to be
     *         be prepared to deal with. The idea, then, is to first make the association of hits in the U and W planes
     *         and then look for the match in the V plane. In the event we don't find the match in the V plane then we
     *         will evaluate the situation and in some instances keep the U-W pairs in order to keep efficiency high.
     */
    
    // Define functions to set start/end iterators in the loop below
    auto SetStartIterator = [](HitVector::iterator startItr, HitVector::iterator endItr, float rms, float startTime)
    {
        while(startItr != endItr)
        {
            float numRMS(rms);
            if ((*startItr)->getTimeTicks() + numRMS * (*startItr)->getHit().RMS() < startTime) startItr++;
            else break;
        }
        return startItr;
    };
    
    auto SetEndIterator = [](HitVector::iterator firstItr, HitVector::iterator endItr, float rms, float endTime)
    {
        while(firstItr != endItr)
        {
            float numRMS(rms);
            if ((*firstItr)->getTimeTicks() - numRMS * (*firstItr)->getHit().RMS() < endTime) firstItr++;
            else break;
        }
        return firstItr;
    };
    
    size_t nTriplets(0);
    size_t nDeadChanHits(0);
    
    //*********************************************************************************
    // Basically, we try to loop until done...
    while(1)
    {
        // Sort so that the earliest hit time will be the first element, etc.
        std::sort(hitItrVec.begin(),hitItrVec.end(),SetStartTimeOrder(m_numSigmaPeakTime));
        
        // This loop iteration's golden hit
        const reco::ClusterHit2D* goldenHit = *hitItrVec[0].first;
        
        // The range of history... (for this hit)
        float goldenTimeStart = goldenHit->getTimeTicks() - m_numSigmaPeakTime * goldenHit->getHit().RMS() - 0.1;
        float goldenTimeEnd   = goldenHit->getTimeTicks() + m_numSigmaPeakTime * goldenHit->getHit().RMS() + 0.1;
        
        // Set iterators to insure we'll be in the overlap ranges
        HitVector::iterator hitItr1Start = SetStartIterator(hitItrVec[1].first, hitItrVec[1].second, m_numSigmaPeakTime, goldenTimeStart);
        HitVector::iterator hitItr1End   = SetEndIterator( hitItr1Start,        hitItrVec[1].second, m_numSigmaPeakTime, goldenTimeEnd);
        HitVector::iterator hitItr2Start = SetStartIterator(hitItrVec[2].first, hitItrVec[2].second, m_numSigmaPeakTime, goldenTimeStart);
        HitVector::iterator hitItr2End   = SetEndIterator( hitItr2Start,        hitItrVec[2].second, m_numSigmaPeakTime, goldenTimeEnd);
        
        // Since we'll use these many times in the internal loops, pre make the pairs for the second set of hits
        size_t             curHitListSize(hitPairList.size());
        HitMatchPairVecMap pair12Map;
        HitMatchPairVecMap pair13Map;
        
        size_t n12Pairs = findGoodHitPairs(goldenHit, hitItr1Start, hitItr1End, hitPairList, pair12Map);
        size_t n13Pairs = findGoodHitPairs(goldenHit, hitItr2Start, hitItr2End, hitPairList, pair13Map);
        
        nDeadChanHits  += hitPairList.size() - curHitListSize;
        curHitListSize  = hitPairList.size();
        
        if (n12Pairs > n13Pairs) findGoodTriplets(pair12Map, pair13Map, hitPairList);
        else                     findGoodTriplets(pair13Map, pair12Map, hitPairList);
        
        nTriplets += hitPairList.size() - curHitListSize;
        
        hitItrVec[0].first++;
        
        int nPlanesWithHits(0);
        
        for(auto& pair : hitItrVec)
            if (pair.first != pair.second) nPlanesWithHits++;
        
        if (nPlanesWithHits < 2) break;
    }
    
    return hitPairList.size();
}
    
size_t Hit3DBuilderAlg::findGoodHitPairs(const reco::ClusterHit2D* goldenHit,
                                         HitVector::iterator&      startItr,
                                         HitVector::iterator&      endItr,
                                         reco::HitPairList&        hitPairList,
                                         HitMatchPairVecMap&       hitMatchMap) const
{
    // temporary container for dead channel hits
    std::vector<reco::ClusterHit3D> tempDeadChanVec;
    reco::ClusterHit3D              deadChanPair;

    // Temporary container for pairs
    HitMatchPairVec tempPairVec;
    
    while(startItr != endItr)
    {
        reco::ClusterHit2D* hit = *startItr++;
        reco::ClusterHit3D  pair;
        
        makeHitPair(pair, goldenHit, hit, m_hitWidthSclFctr);
        
        if (!(pair.getAvePeakTime() > 0.)) continue;
        
        tempPairVec.emplace_back(HitMatchPair(hit,pair));
    }

    // Can we try to weed out extra hit pairs and keep only the "best"?
    if (tempPairVec.size() > 1)
    {
        // Sort by "significance" of agreement
        std::sort(tempPairVec.begin(),tempPairVec.end(),[](HitMatchPair& left, HitMatchPair& right){return left.second.getDeltaPeakTime()/left.second.getSigmaPeakTime() < right.second.getDeltaPeakTime()/right.second.getSigmaPeakTime();});
    }
    
    for(auto& pair : tempPairVec)
    {
        const reco::ClusterHit2D* hit   = pair.first;
        reco::ClusterHit3D&       pair2 = pair.second;
    
        geo::WireID wireID = hit->getHit().WireID();
    
        hitMatchMap[wireID.Wire].emplace_back(HitMatchPair(hit,pair2));
    }
    
    return tempPairVec.size();
}
    
void Hit3DBuilderAlg::findGoodTriplets(HitMatchPairVecMap& pair12Map, HitMatchPairVecMap& pair13Map, reco::HitPairList& hitPairList) const
{
    // Build triplets from the two lists of hit pairs
    if (!pair12Map.empty())
    {
        // temporary container for dead channel hits
        std::vector<reco::ClusterHit3D> tempDeadChanVec;
        reco::ClusterHit3D              deadChanPair;
        
        // Keep track of which third plane hits have been used
        std::map<const reco::ClusterHit3D*,bool> usedPairMap;
        
        // Initial population of this map with the pair13Map hits
        for(const auto& pair13 : pair13Map)
        {
            for(const auto& hit2Dhit3DPair : pair13.second) usedPairMap[&hit2Dhit3DPair.second] = false;
        }
        
        // The outer loop is over all hit pairs made from the first two plane combinations
        for(const auto& pair12 : pair12Map)
        {
            if (pair12.second.empty()) continue;

            // "Discover" the missing view (and we can't rely on assuming there are hits in the pair13Map at this point)
            size_t missPlane = 0;
            
            if      (!pair12.second.front().second.getHits()[1]) missPlane = 1;
            else if (!pair12.second.front().second.getHits()[2]) missPlane = 2;
            
            geo::WireID missingPlaneID(0,0,missPlane,0);
            
            // This loop is over hit pairs that share the same first two plane wires but may have different
            // hit times on those wires
            for(const auto& hit2Dhit3DPair : pair12.second)
            {
                const reco::ClusterHit3D& pair1  = hit2Dhit3DPair.second;
                
                // Get the wire ID for the nearest wire to the position of this hit
                geo::WireID wireID = NearestWireID(pair1.getPosition(), missingPlaneID);
                
                // populate the map with initial value
                usedPairMap[&pair1] = false;
                
                // Now look up the hit pairs on the wire which matches the current hit pair
                std::map<size_t,HitMatchPairVec>::iterator thirdPlaneHitMapItr = pair13Map.find(wireID.Wire);
                
                // Watch for the interesting case of round off error... which means we look at the next wire
                if (thirdPlaneHitMapItr == pair13Map.end()) thirdPlaneHitMapItr = pair13Map.find(wireID.Wire+1);
                
                // Loop over third plane hits and try to form a triplet
                if (thirdPlaneHitMapItr != pair13Map.end())
                {
                    for(const auto& thirdPlaneHitItr : thirdPlaneHitMapItr->second)
                    {
                        const reco::ClusterHit2D* hit2  = thirdPlaneHitItr.first;
                        const reco::ClusterHit3D& pair2 = thirdPlaneHitItr.second;
                        
                        // If success try for the triplet
                        reco::ClusterHit3D triplet;
                        
                        if (makeHitTriplet(triplet, pair1, hit2))
                        {
                            triplet.setID(hitPairList.size());
                            hitPairList.emplace_back(std::unique_ptr<reco::ClusterHit3D>(new reco::ClusterHit3D(triplet)));
                            usedPairMap[&pair1] = true;
                            usedPairMap[&pair2] = true;
                        }
                    }
                }
            }
        }

        // One more loop through the other pairs to check for sick channels
        for(const auto& pairMapPair : usedPairMap)
        {
            if (pairMapPair.second) continue;

            const reco::ClusterHit3D* pair = pairMapPair.first;

            // Here we look to see if we failed to make a triplet because the partner wire was dead/noisy/sick
            if (makeDeadChannelPair(deadChanPair, *pair, 4, 0, 0.)) tempDeadChanVec.emplace_back(deadChanPair);
        }

        // Handle the dead wire triplets
        if(!tempDeadChanVec.empty())
        {
            // If we have many then see if we can trim down a bit by keeping those with time significance
            if (tempDeadChanVec.size() > 1)
            {
                // Sort by "significance" of agreement
                std::sort(tempDeadChanVec.begin(),tempDeadChanVec.end(),[](const auto& left, const auto& right){return left.getDeltaPeakTime()/left.getSigmaPeakTime() < right.getDeltaPeakTime()/right.getSigmaPeakTime();});
                
                // What is the range of "significance" from first to last?
                float firstSig = tempDeadChanVec.front().getDeltaPeakTime() / tempDeadChanVec.front().getSigmaPeakTime();
                float lastSig  = tempDeadChanVec.back().getDeltaPeakTime()  / tempDeadChanVec.back().getSigmaPeakTime();
                float sigRange = lastSig - firstSig;
                
                if (lastSig > 0.5 * m_deltaPeakTimeSig && sigRange > 0.5)
                {
                    // Declare a maximum of 1.5 * the average of the first and last pairs...
                    float maxSignificance = std::max(0.75 * (firstSig + lastSig),1.0);
                    
                    std::vector<reco::ClusterHit3D>::iterator firstBadElem = std::find_if(tempDeadChanVec.begin(),tempDeadChanVec.end(),[&maxSignificance](const auto& pair){return pair.getDeltaPeakTime()/pair.getSigmaPeakTime() > maxSignificance;});
                    
                    // But only keep the best 10?
                    if (std::distance(tempDeadChanVec.begin(),firstBadElem) > 20) firstBadElem = tempDeadChanVec.begin() + 20;
                    // Keep at least one hit...
                    else if (firstBadElem == tempDeadChanVec.begin()) firstBadElem++;
                    
                    tempDeadChanVec.resize(std::distance(tempDeadChanVec.begin(),firstBadElem));
                }
            }

            for(auto& pair : tempDeadChanVec)
            {
                pair.setID(hitPairList.size());
                hitPairList.emplace_back(std::unique_ptr<reco::ClusterHit3D>(new reco::ClusterHit3D(pair)));
            }
        }
    }
    
    return;
}

bool Hit3DBuilderAlg::makeHitPair(reco::ClusterHit3D&       hitPair,
                                  const reco::ClusterHit2D* hit1,
                                  const reco::ClusterHit2D* hit2,
                                  float                     hitWidthSclFctr,
                                  size_t                    hitPairCntr) const
{
    // Assume failure
    bool result(false);
    
    // We assume in this routine that we are looking at hits in different views
    // The first mission is to check that the wires intersect
    const geo::WireID& hit1WireID = hit1->getHit().WireID();
    const geo::WireID& hit2WireID = hit2->getHit().WireID();
    
    geo::WireIDIntersection widIntersect;
    
    if (m_geometry->WireIDsIntersect(hit1WireID, hit2WireID, widIntersect))
    {
        // Wires intersect so now we can check the timing
        float hit1Peak  = hit1->getTimeTicks();
        float hit1Sigma = hit1->getHit().RMS();
        
        float hit2Peak  = hit2->getTimeTicks();
        float hit2Sigma = hit2->getHit().RMS();

        float hit1Width = hitWidthSclFctr * hit1Sigma;
        float hit2Width = hitWidthSclFctr * hit2Sigma;
        
        // Coarse check hit times are "in range"
        if (fabs(hit1Peak - hit2Peak) <= (hit1Width + hit2Width))
        {
            // Check to see that hit peak times are consistent with each other
            float hit1SigSq     = hit1Sigma * hit1Sigma;
            float hit2SigSq     = hit2Sigma * hit2Sigma;
            float avePeakTime   = (hit1Peak / hit1SigSq + hit2Peak / hit2SigSq) * hit1SigSq * hit2SigSq / (hit1SigSq + hit2SigSq);
            float deltaPeakTime = std::fabs(hit1Peak - hit2Peak);
            float sigmaPeakTime = std::sqrt(hit1SigSq + hit2SigSq);
            
            // delta peak time consistency check here
            if (deltaPeakTime < m_deltaPeakTimeSig * sigmaPeakTime)    // 2 sigma consistency? (do this way to avoid divide)
            {
                float totalCharge   = hit1->getHit().Integral() + hit2->getHit().Integral();
            
                float xPositionHit1(hit1->getXPosition());
                float xPositionHit2(hit2->getXPosition());
                float xPosition = (xPositionHit1 / hit1SigSq + xPositionHit2 / hit2SigSq) * hit1SigSq * hit2SigSq / (hit1SigSq + hit2SigSq);
            
                float position[] = {xPosition, float(widIntersect.y), float(widIntersect.z)};
            
                // If to here then we need to sort out the hit pair code telling what views are used
                unsigned statusBits = 1 << hit1->getHit().WireID().Plane | 1 << hit2->getHit().WireID().Plane;
                
                // handle status bits for the 2D hits
                if (hit1->getStatusBits() & reco::ClusterHit2D::USEDINPAIR) hit1->setStatusBit(reco::ClusterHit2D::SHAREDINPAIR);
                if (hit2->getStatusBits() & reco::ClusterHit2D::USEDINPAIR) hit2->setStatusBit(reco::ClusterHit2D::SHAREDINPAIR);
                
                hit1->setStatusBit(reco::ClusterHit2D::USEDINPAIR);
                hit2->setStatusBit(reco::ClusterHit2D::USEDINPAIR);
                
                std::vector<const reco::ClusterHit2D*> hitVector = {0,0,0};
                
                hitVector.at(hit1->getHit().WireID().Plane) = hit1;
                hitVector.at(hit2->getHit().WireID().Plane) = hit2;
                
                // And get the wire IDs
                std::vector<geo::WireID> wireIDVec = {geo::WireID(0,0,geo::kU,0), geo::WireID(0,0,geo::kV,0), geo::WireID(0,0,geo::kW,0)};
                
                wireIDVec[hit1->getHit().WireID().Plane] = hit1->getHit().WireID();
                wireIDVec[hit2->getHit().WireID().Plane] = hit2->getHit().WireID();
                
                // For compiling at the moment
                std::vector<float> hitDelTSigVec = {0.,0.,0.};
                
                hitDelTSigVec.at(hit1->getHit().WireID().Plane) = deltaPeakTime / sigmaPeakTime;
                hitDelTSigVec.at(hit2->getHit().WireID().Plane) = deltaPeakTime / sigmaPeakTime;
                
                // Create the 3D cluster hit
                hitPair = reco::ClusterHit3D(hitPairCntr,
                                             statusBits,
                                             position,
                                             totalCharge,
                                             avePeakTime,
                                             deltaPeakTime,
                                             0.5 * sigmaPeakTime,
                                             0.,
                                             0.,
                                             hitDelTSigVec,
                                             wireIDVec,
                                             hitVector);
                
                result = true;
            }
        }
    }
    // Send it back
    return result;
}

    
bool Hit3DBuilderAlg::makeHitTriplet(reco::ClusterHit3D&       hitTriplet,
                                     const reco::ClusterHit3D& pair,
                                     const reco::ClusterHit2D* hit) const
{
    // Assume failure
    bool result(false);
    
    // Recover all the hits involved
    const reco::ClusterHit2D* hit0(pair.getHits()[0]);
    const reco::ClusterHit2D* hit1(pair.getHits()[1]);
    
    if      (!hit0) hit0 = pair.getHits()[2];
    else if (!hit1) hit1 = pair.getHits()[2];
    
    // Let's do a quick consistency check on the input hits to make sure we are in range...
    // Require the W hit to be "in range" with the UV Pair
    if (fabs(hit->getTimeTicks() - pair.getAvePeakTime()) < m_hitWidthSclFctr * (pair.getSigmaPeakTime() + hit->getHit().RMS()))
    {
        // Timing in range, now check that the input hit wire "intersects" with the input pair's wires
        geo::WireID wireID   = NearestWireID(pair.getPosition(), hit->getHit().WireID());
        
        // There is an interesting round off issue that we need to watch for...
        if (wireID.Wire == hit->getHit().WireID().Wire || wireID.Wire + 1 == hit->getHit().WireID().Wire)
        {
            // Use the existing code to see the U and W hits are willing to pair with the V hit
            reco::ClusterHit3D pair0h;
            reco::ClusterHit3D pair1h;
            
            // If good pairs made here then we can try to make a triplet
            if (makeHitPair(pair0h, hit0, hit, m_hitWidthSclFctr) && makeHitPair(pair1h, hit1, hit, m_hitWidthSclFctr))
            {
                std::vector<const reco::ClusterHit3D*> pairVec;
            
                pairVec.resize(3);
            
                pairVec[hit->getHit().WireID().Plane]  = &pair;
                pairVec[hit0->getHit().WireID().Plane] = &pair1h;
                pairVec[hit1->getHit().WireID().Plane] = &pair0h;
            
                float deltaZ_w  = pairVec[2]->getPosition()[2] - 0.5 * (pairVec[0]->getPosition()[2] + pairVec[1]->getPosition()[2]);
                float deltaY_uv = pairVec[1]->getPosition()[1] - pairVec[0]->getPosition()[1];

                // The intersection of wires on 3 planes is actually an equilateral triangle... Each pair will have its position at one of the
                // corners, the difference in distance along the z axis will be 1/2 wire spacing, the difference along the y axis is
                // 1/2 wire space / cos(pitch)
                if (std::fabs(std::fabs(deltaZ_w) - 0.5 * m_wirePitch[2]) < .05 && std::fabs(std::fabs(deltaY_uv) - 0.5774 * m_wirePitch[2]) < 0.05)
                {
                    // Weighted average, delta and sigmas
                    float hitSigma      = hit->getHit().RMS();
                    float hit0Sigma     = hit0->getHit().RMS();
                    float hit1Sigma     = hit1->getHit().RMS();
                    float hitWidWeight  = 1. / (hitSigma  * hitSigma);
                    float hit0WidWeight = 1. / (hit0Sigma * hit0Sigma);
                    float hit1WidWeight = 1. / (hit1Sigma * hit1Sigma);
                    float denominator   = 1. / (hitWidWeight + hit0WidWeight + hit1WidWeight);
                    float avePeakTime   = (hit->getTimeTicks() * hitWidWeight + hit0->getTimeTicks() * hit0WidWeight + hit1->getTimeTicks() * hit1WidWeight) * denominator;
                    
                    // The x position is a weighted sum but the y-z position is simply the average
                    float xPosition   = (hit->getXPosition() * hitWidWeight + hit0->getXPosition() * hit0WidWeight + hit1->getXPosition() * hit1WidWeight) * denominator;
                    float position[]  = { xPosition,
                                          float((pair.getPosition()[1] + pair0h.getPosition()[1] + pair1h.getPosition()[1]) / 3.),
                                          float((pair.getPosition()[2] + pair0h.getPosition()[2] + pair1h.getPosition()[2]) / 3.)};
                    float totalCharge = pair.getTotalCharge() + pair0h.getTotalCharge() + pair1h.getTotalCharge();
                    
                    std::vector<const reco::ClusterHit2D*> hitVector(3);
            
                    // Make sure we have the hits
                    hitVector[hit0->getHit().WireID().Plane] = hit0;
                    hitVector[hit1->getHit().WireID().Plane] = hit1;
                    hitVector[hit->getHit().WireID().Plane]  = hit;

                    // And get the wire IDs
                    std::vector<geo::WireID> wireIDVec = {geo::WireID(0,0,geo::kU,0), geo::WireID(0,0,geo::kV,0), geo::WireID(0,0,geo::kW,0)};
                    
                    for(const auto& hit : hitVector)
                    {
                        wireIDVec[hit->getHit().WireID().Plane] = hit->getHit().WireID();
                        
                        if (hit->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);
                        
                        hit->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);
                    }
                    
                    unsigned int statusBits(0x7);
                    
                    // For compiling at the moment
                    std::vector<float> hitDelTSigVec = {0.,0.,0.};
                    
                    float hitPairDeltaT   = std::fabs(hit->getTimeTicks()-pair.getAvePeakTime());
                    float hitPairSig      = std::sqrt(hitSigma*hitSigma + pair.getSigmaPeakTime()*pair.getSigmaPeakTime());
                    float hit0Pair1DeltaT = std::fabs(hit0->getTimeTicks()-pair1h.getAvePeakTime());
                    float hit0Pair1Sig    = std::sqrt(hit0Sigma*hit0Sigma + pair1h.getSigmaPeakTime()*pair1h.getSigmaPeakTime());
                    float hit1Pair0DeltaT = std::fabs(hit1->getTimeTicks()-pair0h.getAvePeakTime());
                    float hit1Pair0Sig    = std::sqrt(hit1Sigma*hit1Sigma + pair0h.getSigmaPeakTime()*pair0h.getSigmaPeakTime());
                    
                    // Want deltaPeakTime and sigmaPeakTime to be the worst of the lot...
                    float deltaPeakTime = hitPairDeltaT;
                    float sigmaPeakTime = hitPairSig;
                    
                    if (deltaPeakTime/sigmaPeakTime < hit0Pair1DeltaT/hit0Pair1Sig)
                    {
                        deltaPeakTime = hit0Pair1DeltaT;
                        sigmaPeakTime = hit0Pair1Sig;
                    }
                    
                    if (deltaPeakTime/sigmaPeakTime < hit1Pair0DeltaT/hit1Pair0Sig)
                    {
                        deltaPeakTime = hit1Pair0DeltaT;
                        sigmaPeakTime = hit1Pair0Sig;
                    }
                    
                    hitDelTSigVec.at(hit->getHit().WireID().Plane ) = hitPairDeltaT   / hitPairSig;
                    hitDelTSigVec.at(hit0->getHit().WireID().Plane) = hit0Pair1DeltaT / hit0Pair1Sig;
                    hitDelTSigVec.at(hit1->getHit().WireID().Plane) = hit1Pair0DeltaT / hit1Pair0Sig;
                    
                    // Create the 3D cluster hit
                    hitTriplet = reco::ClusterHit3D(0,
                                                    statusBits,
                                                    position,
                                                    totalCharge,
                                                    avePeakTime,
                                                    deltaPeakTime,
                                                    sigmaPeakTime,
                                                    0.,
                                                    0.,
                                                    hitDelTSigVec,
                                                    wireIDVec,
                                                    hitVector);
                    
                    result = true;
                }
//                else std::cout << "3D Build --> rejecting triplet by position, delta z: " << deltaZ_w << ", delta y: " << deltaY_uv << std::endl;
            }
        }
    }
    
    // return success/fail
    return result;
}
    
bool Hit3DBuilderAlg::makeDeadChannelPair(reco::ClusterHit3D&       pairOut,
                                          const reco::ClusterHit3D& pair,
                                          size_t                    maxChanStatus,
                                          size_t                    minChanStatus,
                                          float                     minOverlap) const
{
    // Assume failure (most common result)
    bool result(false);
    
    const reco::ClusterHit2D* hit0 = pair.getHits().at(0);
    const reco::ClusterHit2D* hit1 = pair.getHits().at(1);
    
    size_t missPlane(2);
    
    // u plane hit is missing
    if (!hit0)
    {
        hit0      = pair.getHits().at(2);
        missPlane = 0;
    }
    // v plane hit is missing
    else if (!hit1)
    {
        hit1      = pair.getHits().at(2);
        missPlane = 1;
    }
    
    // Which plane is missing?
    geo::WireID wireID0 = hit0->getHit().WireID();
    geo::WireID wireID1 = hit1->getHit().WireID();
    
    // Ok, recover the wireID expected in the third plane...
    geo::WireID wireIn(0,0,missPlane,0);
    geo::WireID wireID = NearestWireID(pair.getPosition(), wireIn);
    
    // There can be a round off issue so check the next wire as well
    bool wireStatus    = m_channelStatus[wireID.Plane][wireID.Wire]   < maxChanStatus && m_channelStatus[wireID.Plane][wireID.Wire]   >= minChanStatus;
    bool wireOneStatus = m_channelStatus[wireID.Plane][wireID.Wire+1] < maxChanStatus && m_channelStatus[wireID.Plane][wireID.Wire+1] >= minChanStatus;
    
    // Make sure they are of at least the minimum status
    if(wireStatus || wireOneStatus)
    {
        // Sort out which is the wire we're dealing with
        if (!wireStatus) wireID.Wire += 1;
    
        // Want to refine position since we "know" the missing wire
        geo::WireIDIntersection widIntersect0;
    
        if (m_geometry->WireIDsIntersect(wireID0, wireID, widIntersect0))
        {
            geo::WireIDIntersection widIntersect1;
        
            if (m_geometry->WireIDsIntersect(wireID1, wireID, widIntersect1))
            {
                float newPosition[] = {pair.getPosition()[0],pair.getPosition()[1],pair.getPosition()[2]};
            
                newPosition[1] = (newPosition[1] + widIntersect0.y + widIntersect1.y) / 3.;
                newPosition[2] = (newPosition[2] + widIntersect0.z + widIntersect1.z) / 3.;
            
                pairOut = pair;
                pairOut.setWireID(wireID);
                pairOut.setPosition(newPosition);
                
                if (hit0->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit0->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);
                if (hit1->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit1->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);
                
                hit0->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);
                hit1->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);
            
                result  = true;
            }
        }
    }
    
    return result;
}
    
const reco::ClusterHit2D* Hit3DBuilderAlg::FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float pairDeltaTimeLimits) const
{
    static const float minCharge(0.);
    
    const reco::ClusterHit2D* bestVHit(0);
    
    float pairAvePeakTime(pair.getAvePeakTime());
    
    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit().Integral() < minCharge) continue;
        
        float hitVPeakTime(hit2D->getTimeTicks());
        float deltaPeakTime(pairAvePeakTime-hitVPeakTime);
        
        if (deltaPeakTime >  pairDeltaTimeLimits) continue;
        
        if (deltaPeakTime < -pairDeltaTimeLimits) break;

        pairDeltaTimeLimits = fabs(deltaPeakTime);
        bestVHit            = hit2D;
    }

    return bestVHit;
}
    
int Hit3DBuilderAlg::FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float range) const
{
    static const float minCharge(0.);
    
    int    numberInRange(0);
    float pairAvePeakTime(pair.getAvePeakTime());
    
    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit().Integral() < minCharge) continue;
        
        float hitVPeakTime(hit2D->getTimeTicks());
        float deltaPeakTime(pairAvePeakTime-hitVPeakTime);
        
        if (deltaPeakTime >  range) continue;
        
        if (deltaPeakTime < -range) break;
        
        numberInRange++;
    }
    
    return numberInRange;
}

geo::WireID Hit3DBuilderAlg::NearestWireID(const float* position, const geo::WireID& wireIDIn) const
{
    geo::WireID wireID(wireIDIn,0);
    
    // Embed the call to the geometry's services nearest wire id method in a try-catch block
    try
    {
        wireID = m_geometry->NearestWireID(position, wireIDIn.Plane);
    }
    catch(std::exception& exc)
    {
        // This can happen, almost always because the coordinates are **just** out of range
        mf::LogWarning("Cluster3D") << "Exception caught finding nearest wire, position - " << exc.what() << std::endl;
        
        // Assume extremum for wire number depending on z coordinate
        if (position[2] < 0.5 * m_geometry->DetLength()) wireID.Wire = 0;
        else                                             wireID.Wire = m_geometry->Nwires(wireIDIn.Plane) - 1;
    }
    
    return wireID;
}
    

} // namespace lar_cluster3d

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
    m_enableMonitoring       = pset.get<bool>  ("EnableMonitoring",  true  );
    m_timeAdvanceGap         = pset.get<double>("TimeAdvanceGap",   50.    );
    m_numSigmaPeakTime       = pset.get<double>("NumSigmaPeakTime",  3.    );
    m_pairOverlapSmall       = pset.get<double>("PairOverlapSmall",  0.5   ); //0.6   );
    m_pairOverlapLarge       = pset.get<double>("PairOverlapLarge",  0.1   ); //0.2   );
    
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
    SetStartTimeOrder()              : m_numRMS(1.)     {}
    SetStartTimeOrder(double numRMS) : m_numRMS(numRMS) {}
    
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
    double m_numRMS;
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
    auto SetStartIterator = [](HitVector::iterator startItr, HitVector::iterator endItr, double rms, double startTime)
    {
        while(startItr != endItr)
        {
            double numRMS(rms);
            if ((*startItr)->getTimeTicks() + numRMS * (*startItr)->getHit().RMS() < startTime) startItr++;
            else break;
        }
        return startItr;
    };
    
    auto SetEndIterator = [](HitVector::iterator firstItr, HitVector::iterator endItr, double rms, double endTime)
    {
        while(firstItr != endItr)
        {
            double numRMS(rms);
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
        
        // The end of time... (for this hit)
        double numRMS = m_numSigmaPeakTime;
        
        // The idea here is that if we have a single gaussian representing an extra wide pulse then
        // we might consider opening the range event a bit more...
        //        if (goldenHit->getHit().DegreesOfFreedom() == 1) numRMS *= 1.5;
        
        double goldenTimeStart = goldenHit->getTimeTicks() - numRMS * goldenHit->getHit().RMS() - 0.1;
        double goldenTimeEnd   = goldenHit->getTimeTicks() + numRMS * goldenHit->getHit().RMS() + 0.1;
        
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
        
        makeHitPair(pair, goldenHit, hit, m_numSigmaPeakTime);
        
        if (!(pair.getAvePeakTime() > 0.)) continue;
        
        tempPairVec.emplace_back(HitMatchPair(hit,pair));
    }
    
    if (tempPairVec.size() > 100)
    {
        std::sort(tempPairVec.begin(),tempPairVec.end(),[](HitMatchPair& left, HitMatchPair& right){return left.second.getMaxOverlapFraction() > right.second.getMaxOverlapFraction();});
    
        //double minOverlap = std::min(0.75 * tempPairVec.front().second.getMaxOverlapFraction(),0.5);
        double minOverlap = 0.75 * tempPairVec.front().second.getMaxOverlapFraction();
        
        HitMatchPairVec::iterator firstBadElem = std::find_if(tempPairVec.begin(),tempPairVec.end(),[&minOverlap](HitMatchPair& pair){return pair.second.getMaxOverlapFraction() < minOverlap;});
        
        tempPairVec.resize(std::distance(tempPairVec.begin(),firstBadElem));
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
            
//            // Get the wire ID for the nearest wire to the position of this hit
//            geo::WireID wireID = NearestWireID(pair12.second.front().second.getPosition(), missingPlaneID);
            
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
            if (makeDeadChannelPair(deadChanPair, *pair, 4, 0))
                tempDeadChanVec.emplace_back(deadChanPair);
        }

        // Handle the dead wire triplets
        if(!tempDeadChanVec.empty())
        {
            // If we have many then see if we can trim down a bit by keeping those with the best overlap
            if (tempDeadChanVec.size() > 20)
            {
                std::sort(tempDeadChanVec.begin(),tempDeadChanVec.end(),[](const reco::ClusterHit3D& left, const reco::ClusterHit3D& right){return left.getMaxOverlapFraction() > right.getMaxOverlapFraction();});
                
                double minOverlap = 0.5 * tempDeadChanVec.front().getMaxOverlapFraction();
                
                std::vector<reco::ClusterHit3D>::iterator firstBadElem = std::find_if(tempDeadChanVec.begin(),tempDeadChanVec.end(),[&minOverlap](const reco::ClusterHit3D& pair){return pair.getMaxOverlapFraction() < minOverlap;});
                
                if (firstBadElem != tempDeadChanVec.end()) tempDeadChanVec.resize(std::distance(tempDeadChanVec.begin(),firstBadElem));
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
                                  double                    hitWidthSclFctr,
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
        // Basically, require that the hit times "overlap"
        // Check here is that they are inconsistent
        double hit1Peak  = hit1->getTimeTicks();
        double hit1Sigma = hit1->getHit().RMS();
        
        double hit2Peak  = hit2->getTimeTicks();
        double hit2Sigma = hit2->getHit().RMS();

        double hit1Width = hitWidthSclFctr * hit1Sigma;
        double hit2Width = hitWidthSclFctr * hit2Sigma;
        
//        if (hit1->getHit().DegreesOfFreedom() == 1) hit1Width *= 2.;
//        if (hit2->getHit().DegreesOfFreedom() == 1) hit2Width *= 2.;
        
        // Check hit times are consistent
        if (fabs(hit1Peak - hit2Peak) <= (hit1Width + hit2Width))
        {
            double maxUpper             = std::min(hit1Peak+hit1Width,hit2Peak+hit2Width);
            double minLower             = std::max(hit1Peak-hit1Width,hit2Peak-hit2Width);
            double overlap              = maxUpper - minLower;
            double overlapFractionSmall = 0.5 * overlap / std::min(hit1Width,hit2Width);  // essentially fraction of small pulse contained in larger
            double overlapFractionLarge = 0.5 * overlap / std::max(hit1Width,hit2Width);  // fraction of larger pulse overlapped by smaller
            double minOverlapLarge      = m_pairOverlapLarge;
            
            if (hit1->getHit().DegreesOfFreedom() == 1 || hit2->getHit().DegreesOfFreedom() == 1) minOverlapLarge = 0.;
            
            // require at least 10% overlap of pulses
            if (overlapFractionSmall > m_pairOverlapSmall && overlapFractionLarge > minOverlapLarge)
            {
                double hit1WidSq     = hit1Width * hit1Width;
                double hit2WidSq     = hit2Width * hit2Width;
                double avePeakTime   = (hit1Peak / hit1WidSq + hit2Peak / hit2WidSq) * hit1WidSq * hit2WidSq / (hit1WidSq + hit2WidSq);
                double deltaPeakTime = fabs(hit1Peak - hit2Peak);
                double sigmaPeakTime = std::sqrt(hit1Sigma*hit1Sigma + hit2Sigma*hit2Sigma);
                double totalCharge   = hit1->getHit().Integral() + hit2->getHit().Integral();
            
                double xPositionHit1(hit1->getXPosition());
                double xPositionHit2(hit2->getXPosition());
                double xPosition = (xPositionHit1 / hit1WidSq + xPositionHit2 / hit2WidSq) * hit1WidSq * hit2WidSq / (hit1WidSq + hit2WidSq);
            
                double position[] = {xPosition, widIntersect.y, widIntersect.z};
            
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
                
                // Create the 3D cluster hit
                hitPair = reco::ClusterHit3D(hitPairCntr,
                                             statusBits,
                                             position,
                                             totalCharge,
                                             avePeakTime,
                                             deltaPeakTime,
                                             sigmaPeakTime,
                                             overlapFractionSmall,
                                             overlapFractionSmall,
                                             0.,
                                             0.,
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
    double pairSigma = m_numSigmaPeakTime * pair.getSigmaPeakTime();
    double hitSigma  = m_numSigmaPeakTime * hit->getHit().RMS();
    
    // Require the W hit to be "in range" with the UV Pair
    if (fabs(hit->getTimeTicks() - pair.getAvePeakTime()) < pairSigma + hitSigma)
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
            if (makeHitPair(pair0h, hit0, hit, m_numSigmaPeakTime) && makeHitPair(pair1h, hit1, hit, m_numSigmaPeakTime))
            {
                std::vector<const reco::ClusterHit3D*> pairVec;
            
                pairVec.resize(3);
            
                pairVec[hit->getHit().WireID().Plane]  = &pair;
                pairVec[hit0->getHit().WireID().Plane] = &pair1h;
                pairVec[hit1->getHit().WireID().Plane] = &pair0h;
            
                double deltaZ_w  = pairVec[2]->getPosition()[2] - 0.5 * (pairVec[0]->getPosition()[2] + pairVec[1]->getPosition()[2]);
                double deltaY_uv = pairVec[1]->getPosition()[1] - pairVec[0]->getPosition()[1];

                // The intersection of wires on 3 planes is actually an equilateral triangle... Each pair will have its position at one of the
                // corners, the difference in distance along the z axis will be 1/2 wire spacing, the difference along the y axis is
                // 1/2 wire space / cos(pitch)
                if (std::fabs(std::fabs(deltaZ_w) - 0.5 * m_wirePitch[2]) < .05 && std::fabs(std::fabs(deltaY_uv) - 0.5774 * m_wirePitch[2]) < 0.05)
                {
                    // Weighted average, delta and sigmas
                    double hitSigma      = hit->getHit().RMS();
                    double hit0Sigma     = hit0->getHit().RMS();
                    double hit1Sigma     = hit1->getHit().RMS();
                    double hitWidWeight  = 1. / (hitSigma  * hitSigma);
                    double hit0WidWeight = 1. / (hit0Sigma * hit0Sigma);
                    double hit1WidWeight = 1. / (hit1Sigma * hit1Sigma);
                    double denominator   = 1. / (hitWidWeight + hit0WidWeight + hit1WidWeight);
                    double avePeakTime   = (hit->getTimeTicks() * hitWidWeight + hit0->getTimeTicks() * hit0WidWeight + hit1->getTimeTicks() * hit1WidWeight) * denominator;
                    
                    // The x position is a weighted sum but the y-z position is simply the average
                    double xPosition  = (hit->getXPosition() * hitWidWeight + hit0->getXPosition() * hit0WidWeight + hit1->getXPosition() * hit1WidWeight) * denominator;
                    double position[] = { xPosition,
                                         (pair.getPosition()[1] + pair0h.getPosition()[1] + pair1h.getPosition()[1]) / 3.,
                                         (pair.getPosition()[2] + pair0h.getPosition()[2] + pair1h.getPosition()[2]) / 3.};
                    
                    double deltaPeakTime = std::max(pair.getDeltaPeakTime(), std::max(pair0h.getDeltaPeakTime(), pair1h.getDeltaPeakTime()));
                    double sigmaPeakTime = std::sqrt(hitSigma*hitSigma + hit0Sigma*hit0Sigma + hit1Sigma*hit1Sigma);
                    double totalCharge   = pair.getTotalCharge() + pair0h.getTotalCharge() + pair1h.getTotalCharge();
        
                    // Overlap fraction... hmmm....
                    double maxOverlapFraction = std::max(pair.getMaxOverlapFraction(), std::max(pair0h.getMaxOverlapFraction(), pair1h.getMaxOverlapFraction()));
                    double minOverlapFraction = std::min(pair.getMaxOverlapFraction(), std::min(pair0h.getMaxOverlapFraction(), pair1h.getMaxOverlapFraction()));
                    
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
                    
                    // Create the 3D cluster hit
                    hitTriplet = reco::ClusterHit3D(0,
                                                    7,
                                                    position,
                                                    totalCharge,
                                                    avePeakTime,
                                                    deltaPeakTime,
                                                    sigmaPeakTime,
                                                    maxOverlapFraction,
                                                    minOverlapFraction,
                                                    0.,
                                                    0.,
                                                    wireIDVec,
                                                    hitVector);
                    
                    result = true;
                }
                else std::cout << "3D Build --> rejecting triplet by position, delta z: " << deltaZ_w << ", delta y: " << deltaY_uv << std::endl;
            }
        }
    }
    
    // return success/fail
    return result;
}
    
bool Hit3DBuilderAlg::makeDeadChannelPair(reco::ClusterHit3D& pairOut, const reco::ClusterHit3D& pair, size_t maxChanStatus, size_t minChanStatus, double minOverlap) const
{
    // Assume failure (most common result)
    bool result(false);
    
    // Check minimum overlap
    if (pair.getMaxOverlapFraction() < minOverlap) return result;
    
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
    bool wireOneStatus = wireStatus; //m_channelStatus[wireID.Plane][wireID.Wire+1] < maxChanStatus && m_channelStatus[wireID.Plane][wireID.Wire+1] >= minChanStatus;
    
    // Make sure they are of at least the minimum status
    if(wireStatus || wireOneStatus)
    {
        // Sort out which is the wire we're dealing with
//        if (!wireStatus) wireID.Wire += 1;
    
        // Want to refine position since we "know" the missing wire
        geo::WireIDIntersection widIntersect0;
    
        if (m_geometry->WireIDsIntersect(wireID0, wireID, widIntersect0))
        {
            geo::WireIDIntersection widIntersect1;
        
            if (m_geometry->WireIDsIntersect(wireID1, wireID, widIntersect1))
            {
                double newPosition[] = {pair.getPosition()[0],pair.getPosition()[1],pair.getPosition()[2]};
            
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
    
const reco::ClusterHit2D* Hit3DBuilderAlg::FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, double pairDeltaTimeLimits) const
{
    static const double minCharge(0.);
    
    const reco::ClusterHit2D* bestVHit(0);
    
    double pairAvePeakTime(pair.getAvePeakTime());
    
    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit().Integral() < minCharge) continue;
        
        double hitVPeakTime(hit2D->getTimeTicks());
        double deltaPeakTime(pairAvePeakTime-hitVPeakTime);
        
        if (deltaPeakTime >  pairDeltaTimeLimits) continue;
        
        if (deltaPeakTime < -pairDeltaTimeLimits) break;

        pairDeltaTimeLimits = fabs(deltaPeakTime);
        bestVHit            = hit2D;
    }

    return bestVHit;
}
    
int Hit3DBuilderAlg::FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, double range) const
{
    static const double minCharge(0.);
    
    int    numberInRange(0);
    double pairAvePeakTime(pair.getAvePeakTime());
    
    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit().Integral() < minCharge) continue;
        
        double hitVPeakTime(hit2D->getTimeTicks());
        double deltaPeakTime(pairAvePeakTime-hitVPeakTime);
        
        if (deltaPeakTime >  range) continue;
        
        if (deltaPeakTime < -range) break;
        
        numberInRange++;
    }
    
    return numberInRange;
}

geo::WireID Hit3DBuilderAlg::NearestWireID(const double* position, const geo::WireID& wireIDIn) const
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

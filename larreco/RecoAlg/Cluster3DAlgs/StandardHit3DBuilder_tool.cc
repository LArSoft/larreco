/**
 *  @file   StandardHit3DBuilder_tool.cc
 * 
 *  @brief  This tool provides "standard" 3D hits built (by this tool) from 2D hits
 * 
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "canvas/Utilities/InputTag.h"

#include "larreco/RecoAlg/Cluster3DAlgs/IHit3DBuilder.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

// Eigen
#include <Eigen/Dense>

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {
    
/**
 *   @brief What follows are several highly useful typedefs which we
 *          want to expose to the outside world
 */

// forward declaration to define an ordering function for our hit set
struct Hit2DSetCompare
{
    bool operator() (const reco::ClusterHit2D*, const reco::ClusterHit2D*) const;
};

using HitVector                   = std::vector<reco::ClusterHit2D*>;
using PlaneToHitVectorMap         = std::map<geo::PlaneID, HitVector>;
using TPCToPlaneToHitVectorMap    = std::map<geo::TPCID, PlaneToHitVectorMap>;
using Hit2DVector                 = std::vector<reco::ClusterHit2D>;
using Hit2DSet                    = std::set<const reco::ClusterHit2D*, Hit2DSetCompare>;
using WireToHitSetMap             = std::map<unsigned int, Hit2DSet>;
using PlaneToWireToHitSetMap      = std::map<geo::PlaneID, WireToHitSetMap>;
using TPCToPlaneToWireToHitSetMap = std::map<geo::TPCID, PlaneToWireToHitSetMap>;
using HitVectorMap                = std::map<size_t, HitVector>;

using HitPairVector               = std::vector<std::unique_ptr<reco::ClusterHit3D>>;

/**
 *  @brief  StandardHit3DBuilder class definiton
 */
class StandardHit3DBuilder : virtual public IHit3DBuilder
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit StandardHit3DBuilder(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief  Destructor
     */
    ~StandardHit3DBuilder();
    
    void configure(const fhicl::ParameterSet&) override;
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void Hit3DBuilder(const art::Event &evt, reco::HitPairList& hitPairList, RecobHitToPtrMap&) const override;
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute(IHit3DBuilder::TimeValues index) const override {return m_timeVector.at(index);}
    
private:

    /**
     *  @brief  Extract the ART hits and the ART hit-particle relationships
     *
     *  @param  evt                   the ART event
     *  @param  hit2DVector           A container for the internal Cluster3D 2D hit objects
     *  @param  PlaneToHitVectorMap   A map between view and the internal Cluster3D 2D hit objects
     *  @param  viewToWireToHitSetMap This maps 2D hits to wires and stores by view
     *  @param  hitToPtrMap           This maps our Cluster2D hits back to art Ptr's to reco Hits
     */
    void CollectArtHits(const art::Event& evt,
                        RecobHitToPtrMap& hitToPtrMap) const;
    
    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    void BuildHit3D(reco::HitPairList& hitPairList) const;

    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    size_t BuildHitPairMap(PlaneToHitVectorMap& planeToHitVectorMap, reco::HitPairList& hitPairList) const;
    
    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    using PlaneHitVectorItrPairVec = std::vector<std::pair<HitVector::iterator,HitVector::iterator>>;
    
    size_t BuildHitPairMapByTPC(PlaneHitVectorItrPairVec& planeHitVectorItrPairVec, reco::HitPairList& hitPairList) const;
    
    /**
     *  @brief This builds a list of candidate hit pairs from lists of hits on two planes
     */
    using HitMatchPair       = std::pair<const reco::ClusterHit2D*,reco::ClusterHit3D>;
    using HitMatchPairVec    = std::vector<HitMatchPair>;
    using HitMatchPairVecMap = std::map<geo::WireID,HitMatchPairVec>;
    
    int findGoodHitPairs(const reco::ClusterHit2D*, HitVector::iterator&, HitVector::iterator&, HitMatchPairVecMap&) const;
    
    /**
     *  @brief This algorithm takes lists of hit pairs and finds good triplets
     */
    void findGoodTriplets(HitMatchPairVecMap&, HitMatchPairVecMap&, reco::HitPairList&, bool = false) const;
    
    /**
     *  @brief Make a HitPair object by checking two hits
     */
    bool makeHitPair(reco::ClusterHit3D&       pairOut,
                     const reco::ClusterHit2D* hit1,
                     const reco::ClusterHit2D* hit2,
                     float                     hitWidthSclFctr = 1.,
                     size_t                    hitPairCntr = 0) const;
    
    /**
     *  @brief Make a 3D HitPair object by checking two hits
     */
    bool makeHitTriplet(reco::ClusterHit3D&       pairOut,
                        const reco::ClusterHit3D& pairIn,
                        const reco::ClusterHit2D* hit2) const;
    
    /**
     *  @brief Make a 3D HitPair object from a valid pair and a dead channel in the missing plane
     */
    bool makeDeadChannelPair(reco::ClusterHit3D& pairOut, const reco::ClusterHit3D& pair, size_t maxStatus = 4, size_t minStatus = 0, float minOverlap=0.2) const;
    
    /**
     *  @brief A utility routine for finding a 2D hit closest in time to the given pair
     */
    const reco::ClusterHit2D* FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float pairDeltaTimeLimits) const;
    
    /**
     *  @brief A utility routine for returning the number of 2D hits from the list in a given range
     */
    int FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float range) const;
    
    /**
     *  @brief Jacket the calls to finding the nearest wire in order to intercept the exceptions if out of range
     */
    geo::WireID NearestWireID(const float* position, const geo::WireID& wireID) const;
    
    /**
     *  @brief Create the internal channel status vector (assume will eventually be event-by-event)
     */
    void BuildChannelStatusVec(PlaneToWireToHitSetMap& planeToWiretoHitSetMap) const;
    
    /**
     *  @brief define data structure for keeping track of channel status
     */
    using ChannelStatusVec        = std::vector<size_t>;
    using ChannelStatusByPlaneVec = std::vector<ChannelStatusVec>;
    
    /**
     *  @brief Data members to follow
     */
    art::InputTag                        m_hitFinderTag;
    float                                m_numSigmaPeakTime;
    float                                m_hitWidthSclFctr;
    float                                m_deltaPeakTimeSig;
    
    bool                                 m_enableMonitoring;      ///<
    float                                m_wirePitch[3];
    mutable std::vector<float>           m_timeVector;            ///<
    
    float                                m_zPosOffset;
    
    // Get instances of the primary data structures needed
    mutable Hit2DVector                  m_clusterHit2DMasterVec;
    mutable PlaneToHitVectorMap          m_planeToHitVectorMap;
    mutable PlaneToWireToHitSetMap       m_planeToWireToHitSetMap;
    

    mutable ChannelStatusByPlaneVec      m_channelStatus;
    mutable size_t                       m_numBadChannels;
    
    geo::Geometry*                       m_geometry;              //< pointer to the Geometry service
    const detinfo::DetectorProperties*   m_detector;              //< Pointer to the detector properties
    const lariov::ChannelStatusProvider* m_channelFilter;
};

StandardHit3DBuilder::StandardHit3DBuilder(fhicl::ParameterSet const &pset) :
    m_channelFilter(&art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider())

{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StandardHit3DBuilder::~StandardHit3DBuilder()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void StandardHit3DBuilder::configure(fhicl::ParameterSet const &pset)
{
    m_hitFinderTag     = pset.get<art::InputTag>("HitFinderTag");
    m_enableMonitoring = pset.get<bool>         ("EnableMonitoring",    true);
    m_numSigmaPeakTime = pset.get<float>        ("NumSigmaPeakTime",    3.  );
    m_hitWidthSclFctr  = pset.get<float>        ("HitWidthScaleFactor", 6.  );
    m_deltaPeakTimeSig = pset.get<float>        ("DeltaPeakTimeSig",    1.7 );
    m_zPosOffset       = pset.get<float>        ("ZPosOffset",          0.0 );
    
    art::ServiceHandle<geo::Geometry> geometry;
    
    m_geometry = &*geometry;
    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();

    m_wirePitch[0] = m_geometry->WirePitch(0);
    m_wirePitch[1] = m_geometry->WirePitch(1);
    m_wirePitch[2] = m_geometry->WirePitch(2);
}
    
void StandardHit3DBuilder::BuildChannelStatusVec(PlaneToWireToHitSetMap& planeToWireToHitSetMap) const
{
    // This is called each event, clear out the previous version and start over
    if (!m_channelStatus.empty()) m_channelStatus.clear();

    m_numBadChannels = 0;
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
            m_numBadChannels++;
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
    
void StandardHit3DBuilder::Hit3DBuilder(const art::Event& evt, reco::HitPairList& hitPairList, RecobHitToPtrMap& clusterHitToArtPtrMap) const
{
    // Clear the internal data structures
    m_clusterHit2DMasterVec.clear();
    m_planeToHitVectorMap.clear();
    m_planeToWireToHitSetMap.clear();

    m_timeVector.resize(NUMTIMEVALUES, 0.);

    // Recover the 2D hits and then organize them into data structures which will be used in the
    // DBscan algorithm for building the 3D clusters
    this->CollectArtHits(evt, clusterHitToArtPtrMap);
    
//    if (m_enableMonitoring) theClockArtHits.stop();
    
    // If there are no hits in our view/wire data structure then do not proceed with the full analysis
    if (!m_planeToWireToHitSetMap.empty())
    {
        // Call the algorithm that builds 3D hits
        this->BuildHit3D(hitPairList);
    }
    
    return;
}

void StandardHit3DBuilder::BuildHit3D(reco::HitPairList& hitPairList) const
{
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */
    cet::cpu_timer theClockMakeHits;
    
    if (m_enableMonitoring) theClockMakeHits.start();
    
    // The first task is to take the lists of input 2D hits (a map of view to sorted lists of 2D hits)
    // and then to build a list of 3D hits to be used in downstream processing
    BuildChannelStatusVec(m_planeToWireToHitSetMap);
    
    size_t numHitPairs = BuildHitPairMap(m_planeToHitVectorMap, hitPairList);
    
    if (m_enableMonitoring)
    {
        theClockMakeHits.stop();
        
        m_timeVector[BUILDTHREEDHITS] = theClockMakeHits.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> 3D hit building done, found " << numHitPairs << " 3D Hits" << std::endl;
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
namespace {
//bool SetHitStartTimeOrder(const reco::ClusterHit2D* left, const reco::ClusterHit2D* right)
//{
//    // Sort by "modified start time" of pulse
//    return left->getHit().PeakTime() - left->getHit().RMS() < right->getHit().PeakTime() - right->getHit().RMS();
//}
    
class SetHitEarliestTimeOrder
{
public:
    SetHitEarliestTimeOrder()             : m_numRMS(1.)     {}
    SetHitEarliestTimeOrder(float numRMS) : m_numRMS(numRMS) {}
    
    bool operator()(const reco::ClusterHit2D* left, const reco::ClusterHit2D* right) const
    {
        return left->getTimeTicks() - m_numRMS * left->getHit().RMS() < right->getTimeTicks() - m_numRMS * right->getHit().RMS();
    }
    
public:
    float m_numRMS;
};

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

size_t StandardHit3DBuilder::BuildHitPairMap(PlaneToHitVectorMap& planeToHitVectorMap, reco::HitPairList& hitPairList) const
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
            
            size_t nPlanesWithHits = (mapItr0 != planeToHitVectorMap.end() && !mapItr0->second.empty() ? 1 : 0)
                                   + (mapItr1 != planeToHitVectorMap.end() && !mapItr1->second.empty() ? 1 : 0)
                                   + (mapItr2 != planeToHitVectorMap.end() && !mapItr2->second.empty() ? 1 : 0);
            
            if (nPlanesWithHits < 2) continue;
            
            HitVector& hitVector0 = mapItr0->second;
            HitVector& hitVector1 = mapItr1->second;
            HitVector& hitVector2 = mapItr2->second;
            
            // We are going to resort the hits into "start time" order...
            std::sort(hitVector0.begin(), hitVector0.end(), SetHitEarliestTimeOrder(m_numSigmaPeakTime)); //SetHitStartTimeOrder);
            std::sort(hitVector1.begin(), hitVector1.end(), SetHitEarliestTimeOrder(m_numSigmaPeakTime)); //SetHitStartTimeOrder);
            std::sort(hitVector2.begin(), hitVector2.end(), SetHitEarliestTimeOrder(m_numSigmaPeakTime)); //SetHitStartTimeOrder);
            
            PlaneHitVectorItrPairVec hitItrVec = {HitVectorItrPair(hitVector0.begin(),hitVector0.end()),
                                                  HitVectorItrPair(hitVector1.begin(),hitVector1.end()),
                                                  HitVectorItrPair(hitVector2.begin(),hitVector2.end())};
            
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

size_t StandardHit3DBuilder::BuildHitPairMapByTPC(PlaneHitVectorItrPairVec& hitItrVec, reco::HitPairList& hitPairList) const
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
        float goldenTimeStart = goldenHit->getTimeTicks() - m_numSigmaPeakTime * goldenHit->getHit().RMS() - std::numeric_limits<float>::epsilon();
        float goldenTimeEnd   = goldenHit->getTimeTicks() + m_numSigmaPeakTime * goldenHit->getHit().RMS() + std::numeric_limits<float>::epsilon();
        
        // Set iterators to insure we'll be in the overlap ranges
        HitVector::iterator hitItr1Start = SetStartIterator(hitItrVec[1].first, hitItrVec[1].second, m_numSigmaPeakTime, goldenTimeStart);
        HitVector::iterator hitItr1End   = SetEndIterator( hitItr1Start,        hitItrVec[1].second, m_numSigmaPeakTime, goldenTimeEnd);
        HitVector::iterator hitItr2Start = SetStartIterator(hitItrVec[2].first, hitItrVec[2].second, m_numSigmaPeakTime, goldenTimeStart);
        HitVector::iterator hitItr2End   = SetEndIterator( hitItr2Start,        hitItrVec[2].second, m_numSigmaPeakTime, goldenTimeEnd);
        
        // Since we'll use these many times in the internal loops, pre make the pairs for the second set of hits
        size_t             curHitListSize(hitPairList.size());
        HitMatchPairVecMap pair12Map;
        HitMatchPairVecMap pair13Map;
        
        size_t n12Pairs = findGoodHitPairs(goldenHit, hitItr1Start, hitItr1End, pair12Map);
        size_t n13Pairs = findGoodHitPairs(goldenHit, hitItr2Start, hitItr2End, pair13Map);

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

int StandardHit3DBuilder::findGoodHitPairs(const reco::ClusterHit2D* goldenHit,
                                           HitVector::iterator&      startItr,
                                           HitVector::iterator&      endItr,
                                           HitMatchPairVecMap&       hitMatchMap) const
{
    int numPairs(0);
    
    // Loop through the input secon hits and make pairs
    while(startItr != endItr)
    {
        reco::ClusterHit2D* hit = *startItr++;
        reco::ClusterHit3D  pair;
        
        // pair returned with a negative ave time is signal of failure
        if (!makeHitPair(pair, goldenHit, hit, m_hitWidthSclFctr)) continue;
        
        geo::WireID wireID = hit->getHit().WireID();
        
        hitMatchMap[wireID].emplace_back(hit,pair);
        
        numPairs++;
    }
    
    return numPairs;
}

void StandardHit3DBuilder::findGoodTriplets(HitMatchPairVecMap& pair12Map, HitMatchPairVecMap& pair13Map, reco::HitPairList& hitPairList, bool tagged) const
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
            
            // Use the planeID for the first hit
            geo::WireID missingPlaneID = pair12.first;
            
            // "Discover" the missing view (and we can't rely on assuming there are hits in the pair13Map at this point)
            size_t missPlane = 0;
            
            if      (!pair12.second.front().second.getHits()[1]) missPlane = 1;
            else if (!pair12.second.front().second.getHits()[2]) missPlane = 2;
            
            missingPlaneID.Plane = missPlane;
            
            // This loop is over hit pairs that share the same first two plane wires but may have different
            // hit times on those wires
            for(const auto& hit2Dhit3DPair : pair12.second)
            {
                const reco::ClusterHit3D& pair1  = hit2Dhit3DPair.second;

                // Get the wire ID for the nearest wire to the position of this hit
                geo::WireID wireID = NearestWireID(pair1.getPosition(), missingPlaneID);
                
                // populate the map with initial value
                usedPairMap[&pair1] = false;
                
                // For TPC's with 60 degree wire pitch the position returned for the the pair will
                // lie between two wires in the missing plane. The call to nearestWireID should return
                // the lower of the pair.
                // So we really want to do a loop here so we can consider both wire combinations
                for(int loopIdx = 0; loopIdx < 2; loopIdx++)
                {
                    // Now look up the hit pairs on the wire which matches the current hit pair
                    HitMatchPairVecMap::iterator thirdPlaneHitMapItr = pair13Map.find(wireID);
               
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
                                hitPairList.emplace_back(new reco::ClusterHit3D(triplet));
                                usedPairMap[&pair1] = true;
                                usedPairMap[&pair2] = true;
                            }
                        }
                    }
                    
                    // Now bump the wire id to the next wire and do this again
                    wireID.Wire += 1;
                }
            }
        }
        
        // One more loop through the other pairs to check for sick channels
        if (m_numBadChannels > 0)
        {
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
                    hitPairList.emplace_back(new reco::ClusterHit3D(pair));
                }
            }
        }
    }
    
    return;
}

bool StandardHit3DBuilder::makeHitPair(reco::ClusterHit3D&       hitPair,
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
        
        // ad hoc correction for most bad fits...
        if (hit1Sigma > 2. * hit1->getHit().PeakAmplitude()) hit1Sigma = 2. * hit1->getHit().PeakAmplitude();
        if (hit2Sigma > 2. * hit2->getHit().PeakAmplitude()) hit2Sigma = 2. * hit2->getHit().PeakAmplitude();

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
                float hitChiSquare  = std::pow((hit1Peak - avePeakTime),2) / hit1SigSq
                                    + std::pow((hit2Peak - avePeakTime),2) / hit2SigSq;
                
                float xPositionHit1(hit1->getXPosition());
                float xPositionHit2(hit2->getXPosition());
                float xPosition = (xPositionHit1 / hit1SigSq + xPositionHit2 / hit2SigSq) * hit1SigSq * hit2SigSq / (hit1SigSq + hit2SigSq);
                
                float position[] = {xPosition, float(widIntersect.y), float(widIntersect.z)-m_zPosOffset};
                
                // If to here then we need to sort out the hit pair code telling what views are used
                unsigned statusBits = 1 << hit1->getHit().WireID().Plane | 1 << hit2->getHit().WireID().Plane;
                
                // handle status bits for the 2D hits
                if (hit1->getStatusBits() & reco::ClusterHit2D::USEDINPAIR) hit1->setStatusBit(reco::ClusterHit2D::SHAREDINPAIR);
                if (hit2->getStatusBits() & reco::ClusterHit2D::USEDINPAIR) hit2->setStatusBit(reco::ClusterHit2D::SHAREDINPAIR);
                
                hit1->setStatusBit(reco::ClusterHit2D::USEDINPAIR);
                hit2->setStatusBit(reco::ClusterHit2D::USEDINPAIR);
                
                reco::ClusterHit2DVec hitVector;
                
                hitVector.resize(3, NULL);
                
                hitVector.at(hit1->getHit().WireID().Plane) = hit1;
                hitVector.at(hit2->getHit().WireID().Plane) = hit2;
                
                unsigned int cryostatIdx = hit1->getHit().WireID().Cryostat;
                unsigned int tpcIdx      = hit1->getHit().WireID().TPC;
                
                // Initialize the wireIdVec
                std::vector<geo::WireID> wireIDVec = {geo::WireID(cryostatIdx,tpcIdx,0,0),
                                                      geo::WireID(cryostatIdx,tpcIdx,1,0),
                                                      geo::WireID(cryostatIdx,tpcIdx,2,0)};
                
                wireIDVec[hit1->getHit().WireID().Plane] = hit1->getHit().WireID();
                wireIDVec[hit2->getHit().WireID().Plane] = hit2->getHit().WireID();
                
                // For compiling at the moment
                std::vector<float> hitDelTSigVec = {0.,0.,0.};
                
                hitDelTSigVec.at(hit1->getHit().WireID().Plane) = deltaPeakTime / sigmaPeakTime;
                hitDelTSigVec.at(hit2->getHit().WireID().Plane) = deltaPeakTime / sigmaPeakTime;
                
                // Create the 3D cluster hit
                hitPair.initialize(hitPairCntr,
                                   statusBits,
                                   position,
                                   totalCharge,
                                   avePeakTime,
                                   deltaPeakTime,
                                   sigmaPeakTime,
                                   hitChiSquare,
                                   0.,
                                   0.,
                                   hitVector,
                                   hitDelTSigVec,
                                   wireIDVec);
                
                result = true;
            }
        }
    }

    // Send it back
    return result;
}


bool StandardHit3DBuilder::makeHitTriplet(reco::ClusterHit3D&       hitTriplet,
                                          const reco::ClusterHit3D& pair,
                                          const reco::ClusterHit2D* hit) const
{
    // Assume failure
    bool result(false);
    
    static const float rmsToSig(1.0); //0.75); //0.57735027);
    
    // We are going to force the wire pitch here, some time in the future we need to fix
    static const double wirePitch(0.3);

    // Recover hit info
    float hitTimeTicks = hit->getTimeTicks();
    float hitSigma     = hit->getHit().RMS();
    
    // Special case check
    if (hitSigma > 2. * hit->getHit().PeakAmplitude()) hitSigma = 2. * hit->getHit().PeakAmplitude();

    // Let's do a quick consistency check on the input hits to make sure we are in range...
    // Require the W hit to be "in range" with the UV Pair
    if (fabs(hitTimeTicks - pair.getAvePeakTime()) < m_hitWidthSclFctr * (pair.getSigmaPeakTime() + hitSigma))
    {
        // Timing in range, now check that the input hit wire "intersects" with the input pair's wires
        geo::WireID wireID = NearestWireID(pair.getPosition(), hit->getHit().WireID());

        // There is an interesting round off issue that we need to watch for...
        if (wireID.Wire == hit->getHit().WireID().Wire || wireID.Wire + 1 == hit->getHit().WireID().Wire)
        {
            // Use the existing code to see the U and W hits are willing to pair with the V hit
            reco::ClusterHit3D pair0h;
            reco::ClusterHit3D pair1h;
            
            // Recover all the hits involved
            const reco::ClusterHit2DVec& pairHitVec = pair.getHits();
            const reco::ClusterHit2D*    hit0       = pairHitVec.at(0);
            const reco::ClusterHit2D*    hit1       = pairHitVec.at(1);
            
            if      (!hit0) hit0 = pairHitVec.at(2);
            else if (!hit1) hit1 = pairHitVec.at(2);
            
            // If good pairs made here then we can try to make a triplet
            if (makeHitPair(pair0h, hit0, hit, m_hitWidthSclFctr) && makeHitPair(pair1h, hit1, hit, m_hitWidthSclFctr))
            {
                // We want to make sure the 3 sets of pair are really consistent
                // For TPC's with a 60 degree pitch the wire "intersection" will be an equilateral triangle
                // with equal sides of length  wire pitch / sin(60 degrees) / 2
                // So we make sure all three achieve this
                Eigen::Vector2f pair0hYZVec(pair0h.getPosition()[1],pair0h.getPosition()[2]);
                Eigen::Vector2f pair1hYZVec(pair1h.getPosition()[1],pair1h.getPosition()[2]);
                Eigen::Vector2f pairYZVec(pair.getPosition()[1],pair.getPosition()[2]);
                
                std::vector<float> sideVec = {(pair0hYZVec - pair1hYZVec).norm(),(pair1hYZVec - pairYZVec).norm(),(pairYZVec   - pair0hYZVec).norm()};

                // The three sides will not be identically equal because of numeric issues. It is really sufficient to simply
                // check that the longest side is less than the wire pitch
                if (*std::max_element(sideVec.begin(),sideVec.end()) < wirePitch)
                {
                    // Get a copy of the input hit vector (note the order is by plane - by definition)
                    reco::ClusterHit2DVec hitVector = pair.getHits();
                    
                    // include the new hit
                    hitVector.at(hit->getHit().WireID().Plane) = hit;
                    
                    // Set up to get average peak time, hitChiSquare, etc.
                    unsigned int statusBits(0x7);
                    float        totalCharge(0.);
                    float        avePeakTime(0.);
                    float        weightSum(0.);
                    float        xPosition(0.);
                    
                    // And get the wire IDs
                    std::vector<geo::WireID> wireIDVec = {geo::WireID(0,0,geo::kU,0), geo::WireID(0,0,geo::kV,0), geo::WireID(0,0,geo::kW,0)};
                    
                    // First loop through the hits to get WireIDs and calculate the averages
                    for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
                    {
                        const reco::ClusterHit2D* hit2D = hitVector.at(planeIdx);
                        
                        wireIDVec.at(planeIdx) = hit2D->getHit().WireID();
                        
                        if (hit2D->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit2D->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);
                        
                        hit2D->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);
                        
                        float hitRMS   = rmsToSig * hit2D->getHit().RMS();
                        float weight   = 1. / (hitRMS * hitRMS);
                        float peakTime = hit2D->getTimeTicks();
                        
                        avePeakTime += peakTime * weight;
                        xPosition   += hit2D->getXPosition() * weight;
                        weightSum   += weight;
                        totalCharge += hit2D->getHit().Integral();
                    }
                    
                    avePeakTime /= weightSum;
                    xPosition   /= weightSum;
                    
                    float position[]  = { xPosition,
                                          float((pairYZVec[0] + pair0hYZVec[0] + pair1hYZVec[0]) / 3.),
                                          float((pairYZVec[1] + pair0hYZVec[1] + pair1hYZVec[1]) / 3.)};

                    // Armed with the average peak time, now get hitChiSquare and the sig vec
                    float              hitChiSquare(0.);
                    float              sigmaPeakTime(std::sqrt(1./weightSum));
                    std::vector<float> hitDelTSigVec;
                    
                    for(const auto& hit2D : hitVector)
                    {
                        float hitRMS    = rmsToSig * hit2D->getHit().RMS();
                        float combRMS   = std::sqrt(hitRMS*hitRMS - sigmaPeakTime*sigmaPeakTime);
                        float peakTime  = hit2D->getTimeTicks();
                        float deltaTime = peakTime - avePeakTime;
                        float hitSig    = deltaTime / combRMS; //hitRMS;
                        
                        hitChiSquare += hitSig * hitSig;
                        
                        hitDelTSigVec.emplace_back(std::fabs(hitSig));
                    }
                    
                    // Usurping "deltaPeakTime" to be the maximum pull
                    float deltaPeakTime = *std::max_element(hitDelTSigVec.begin(),hitDelTSigVec.end());
                    
                    // Create the 3D cluster hit
                    hitTriplet.initialize(0,
                                          statusBits,
                                          position,
                                          totalCharge,
                                          avePeakTime,
                                          deltaPeakTime,
                                          sigmaPeakTime,
                                          hitChiSquare,
                                          0.,
                                          0.,
                                          hitVector,
                                          hitDelTSigVec,
                                          wireIDVec);
                    
                    result = true;
                }
            }
        }
    }
    
    // return success/fail
    return result;
}

bool StandardHit3DBuilder::makeDeadChannelPair(reco::ClusterHit3D&       pairOut,
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
    geo::WireID wireIn(wireID0.Cryostat,wireID0.TPC,missPlane,0);
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
                newPosition[2] = (newPosition[2] + widIntersect0.z + widIntersect1.z - 2. * m_zPosOffset) / 3.;
                
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

const reco::ClusterHit2D* StandardHit3DBuilder::FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float pairDeltaTimeLimits) const
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

int StandardHit3DBuilder::FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float range) const
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

geo::WireID StandardHit3DBuilder::NearestWireID(const float* position, const geo::WireID& wireIDIn) const
{
    geo::WireID wireID = wireIDIn;
    
    // Embed the call to the geometry's services nearest wire id method in a try-catch block
    try
    {
        // Switch from NearestWireID to this method to avoid the roundoff error issues...
        double distanceToWire = m_geometry->Plane(wireIDIn).WireCoordinate(position);
        
        wireID.Wire = int(distanceToWire);
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

//------------------------------------------------------------------------------------------------------------------------------------------
bool SetHitTimeOrder(const reco::ClusterHit2D* left, const reco::ClusterHit2D* right)
{
    // Sort by "modified start time" of pulse
    return left->getHit().PeakTime() < right->getHit().PeakTime();
}

bool Hit2DSetCompare::operator() (const reco::ClusterHit2D* left, const reco::ClusterHit2D* right) const
{
    return left->getHit().PeakTime() < right->getHit().PeakTime();
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
void StandardHit3DBuilder::CollectArtHits(const art::Event& evt,
                                          RecobHitToPtrMap& hitToPtrMap) const
{
    /**
     *  @brief Recover the 2D hits from art and fill out the local data structures for the 3D clustering
     */
    art::Handle< std::vector<recob::Hit> > recobHitHandle;
    evt.getByLabel(m_hitFinderTag, recobHitHandle);
    
    if (!recobHitHandle.isValid()) return;
    
    cet::cpu_timer theClockMakeHits;
    
    if (m_enableMonitoring) theClockMakeHits.start();
    
    // We'll want to correct the hit times for the plane offsets
    // (note this is already taken care of when converting to position)
    std::map<geo::PlaneID, double> planeOffsetMap;
    
    // Reserve memory for the hit vector
    m_clusterHit2DMasterVec.reserve(recobHitHandle->size());
    
    // Initialize the plane to hit vector map
    for(size_t cryoIdx = 0; cryoIdx < m_geometry->Ncryostats(); cryoIdx++)
    {
        for(size_t tpcIdx = 0; tpcIdx < m_geometry->NTPC(); tpcIdx++)
        {
            m_planeToHitVectorMap[geo::PlaneID(cryoIdx,tpcIdx,0)] = HitVector();
            m_planeToHitVectorMap[geo::PlaneID(cryoIdx,tpcIdx,1)] = HitVector();
            m_planeToHitVectorMap[geo::PlaneID(cryoIdx,tpcIdx,2)] = HitVector();
            
            // What we want here are the relative offsets between the planes
            // Note that plane 0 is assumed the "first" plane and is the reference
            planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,0)] = 0.;
            planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,1)] = m_detector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,1))
                                                           - m_detector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0));
            planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,2)] = m_detector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,2))
                                                           - m_detector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0));
            
            std::cout << "***> plane 0 offset: " << planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,0)] << ", plane 1: " << planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,1)] << ", plane 2: " << planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,2)] << std::endl;
            std::cout << "     Det prop plane 0: " << m_detector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0)) << ", plane 1: "  << m_detector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,1)) << ", plane 2: " << m_detector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,2)) << ", Trig: " << m_detector->TriggerOffset() << std::endl;
        }
    }

    // Cycle through the recob hits to build ClusterHit2D objects and insert
    // them into the map
    for (size_t cIdx = 0; cIdx < recobHitHandle->size(); cIdx++)
    {
        art::Ptr<recob::Hit> recobHit(recobHitHandle, cIdx);
        
        // Skip junk hits
//        if (recobHit->DegreesOfFreedom() > 1 && recobHit->Multiplicity() > 1 && (recobHit->RMS() < 3.8 || recobHit->PeakAmplitude() < 10)) continue;
        
        const geo::WireID& hitWireID(recobHit->WireID());
        
        double hitPeakTime(recobHit->PeakTime() - planeOffsetMap[recobHit->WireID().planeID()]);
        double xPosition(m_detector->ConvertTicksToX(recobHit->PeakTime(), hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat));
        
        m_clusterHit2DMasterVec.emplace_back(0, 0., 0., xPosition, hitPeakTime, *recobHit);
        
        m_planeToHitVectorMap[recobHit->WireID().planeID()].push_back(&m_clusterHit2DMasterVec.back());
        m_planeToWireToHitSetMap[recobHit->WireID().planeID()][recobHit->WireID().Wire].insert(&m_clusterHit2DMasterVec.back());
        
        const recob::Hit* recobHitPtr = recobHit.get();
        hitToPtrMap[recobHitPtr]      = recobHit;
    }
    
    // Make a loop through to sort the recover hits in time order
    for(auto& hitVectorMap : m_planeToHitVectorMap)
        std::sort(hitVectorMap.second.begin(), hitVectorMap.second.end(), SetHitTimeOrder);
    
    if (m_enableMonitoring)
    {
        theClockMakeHits.stop();
        
        m_timeVector[COLLECTARTHITS] = theClockMakeHits.accumulated_real_time();
    }

    mf::LogDebug("Cluster3D") << ">>>>> Number of ART hits: " << m_clusterHit2DMasterVec.size() << std::endl;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
DEFINE_ART_CLASS_TOOL(StandardHit3DBuilder)
} // namespace lar_cluster3d

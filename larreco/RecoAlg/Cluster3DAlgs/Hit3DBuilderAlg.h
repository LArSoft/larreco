/**
 *  @file   Hit3DBuilderAlg.h
 * 
 *  @brief  This algorithm will create and then cluster 3D hits using DBScan
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef Hit3DBuilderAlg_h
#define Hit3DBuilderAlg_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

// std includes
#include <vector>
#include <list>
#include <set>
#include <map>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
    
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
 *  @brief  Hit3DBuilderAlg class definiton
 */
class Hit3DBuilderAlg
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    Hit3DBuilderAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~Hit3DBuilderAlg();

    void reconfigure(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    void BuildHit3D(PlaneToHitVectorMap& planeToHitVectorMap, PlaneToWireToHitSetMap& planeToWireToHitSetMap, reco::HitPairList& hitPairList);
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute() const {return m_buildTime;}
    
private:
    
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
    using HitMatchPairVecMap = std::map<size_t,HitMatchPairVec>;
    
    size_t findGoodHitPairs(const reco::ClusterHit2D*, HitVector::iterator&, HitVector::iterator&, reco::HitPairList&, HitMatchPairVecMap&) const;
    
    /**
     *  @brief This algorithm takes lists of hit pairs and finds good triplets
     */
    void findGoodTriplets(HitMatchPairVecMap&, HitMatchPairVecMap&, reco::HitPairList&) const;
    
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
    void BuildChannelStatusVec(PlaneToWireToHitSetMap& planeToWiretoHitSetMap);
    
    /**
     *  @brief define data structure for keeping track of channel status
     */
    
    using ChannelStatusVec        = std::vector<size_t>;
    using ChannelStatusByPlaneVec = std::vector<ChannelStatusVec>;
    
    /**
     *  @brief Data members to follow
     */

    float                                m_numSigmaPeakTime;
    float                                m_hitWidthSclFctr;
    float                                m_deltaPeakTimeSig;
    
    bool                                 m_enableMonitoring;      ///<
    float                                m_wirePitch[3];
    float                                m_buildTime;             ///<
    
    ChannelStatusByPlaneVec              m_channelStatus;
 
    geo::Geometry*                       m_geometry;              //< pointer to the Geometry service
    const lariov::ChannelStatusProvider* m_channelFilter;
};
    
} // namespace lar_cluster3d
#endif

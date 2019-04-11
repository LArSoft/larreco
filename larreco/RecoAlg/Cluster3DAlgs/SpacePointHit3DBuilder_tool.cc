/**
 *  @file   SpacePointHit3DBuilder_tool.cc
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
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/ArtDataHelper/HitCreator.h"

// std includes
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

/**
 *  @brief  SpacePointHit3DBuilder class definiton
 */
class SpacePointHit3DBuilder : virtual public IHit3DBuilder
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit SpacePointHit3DBuilder(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief  Destructor
     */
    ~SpacePointHit3DBuilder();
    
    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::EDProducer*) override;

    void configure(const fhicl::ParameterSet&) override;
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void Hit3DBuilder(art::EDProducer&, art::Event &evt, reco::HitPairList& hitPairList, RecobHitToPtrMap&) override;
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute(IHit3DBuilder::TimeValues index) const override {return fTimeVector.at(index);}
    
private:

    using Hit2DVector                 = std::vector<reco::ClusterHit2D>;

    /**
     *  @brief Data members to follow
     */
    art::InputTag                        fSpacePointProducerLabel;
    art::InputTag                        fHitProducerLabel;
    bool                                 fDoWireAssns;
    bool                                 fDoRawDigitAssns;

    bool                                 fEnableMonitoring;       ///<
    mutable std::vector<float>           fTimeVector;             ///<
    
    // Get instances of the primary data structures needed
    mutable Hit2DVector                  m_clusterHit2DMasterVec;
    
    const geo::Geometry*                 fGeometry;               //< pointer to the Geometry service
    const detinfo::DetectorProperties*   fDetector;               //< Pointer to the detector properties
};

SpacePointHit3DBuilder::SpacePointHit3DBuilder(fhicl::ParameterSet const &pset)
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

SpacePointHit3DBuilder::~SpacePointHit3DBuilder()
{
}
    
void SpacePointHit3DBuilder::produces(art::EDProducer* producer)
{
    producer->produces< std::vector<recob::Hit>>();
    
    if (fDoWireAssns)     producer->produces< art::Assns<recob::Wire,   recob::Hit>>();
    if (fDoRawDigitAssns) producer->produces< art::Assns<raw::RawDigit, recob::Hit>>();
    
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void SpacePointHit3DBuilder::configure(fhicl::ParameterSet const &pset)
{
    fSpacePointProducerLabel = pset.get<art::InputTag>("SpacePointProducerLabel"  );
    fHitProducerLabel        = pset.get<art::InputTag>("HitProducerLabel"         );
    fDoWireAssns             = pset.get<bool         >("DoWireAssns",         true);
    fDoRawDigitAssns         = pset.get<bool         >("DoRawDigitAssns",     true);
    fEnableMonitoring        = pset.get<bool>         ("EnableMonitoring",    true);
    
    art::ServiceHandle<geo::Geometry const> geometry;
    
    fGeometry = &*geometry;
    fDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
}
    

void SpacePointHit3DBuilder::Hit3DBuilder(art::EDProducer& prod, art::Event& evt, reco::HitPairList& hitPairList, RecobHitToPtrMap& recobHitToArtPtrMap)
{
    /**
     *  @brief Recover the 2D hits from art and fill out the local data structures for the 3D clustering
     */
    
    fTimeVector.resize(NUMTIMEVALUES, 0.);
    
    cet::cpu_timer theClockMakeHits;
    
    if (fEnableMonitoring) theClockMakeHits.start();

    // Start by recovering the associations between space points and hits
    art::Handle< art::Assns<recob::Hit, recob::SpacePoint> > hitSpacePointAssnsHandle;
    evt.getByLabel(fSpacePointProducerLabel, hitSpacePointAssnsHandle);
    
    if (!hitSpacePointAssnsHandle.isValid()) return;
    
    // Get a hit refiner for the output hit collection
    recob::HitRefinerAssociator hitRefiner(prod, evt, fHitProducerLabel, true, true);
    
    // We need to spin through the associations first to build a map between the SpacePoints and
    // the WireID associated to the collection plane... where for APA style TPCs this will be
    // unambiguous (note that this is not true for ICARUS!)
    using SpacePointToWireIDMap = std::unordered_map<const recob::SpacePoint*,geo::WireID>;
    
    SpacePointToWireIDMap spacePointToWireIDMap;
    
    for(const auto& assnPair : *hitSpacePointAssnsHandle)
    {
        if (assnPair.first->SignalType() == geo::kCollection)
            spacePointToWireIDMap[assnPair.second.get()] = assnPair.first->WireID();
    }

    // First step is to loop through and get a mapping between space points and associated hits
    // and, importantly, a list of unique hits (and mapping between art ptr and hit)
    using OldHitToNewHitMap   = std::map<const recob::Hit*,const recob::Hit*>;
    using SpacePointHitVecMap = std::map<const recob::SpacePoint*, std::vector<const recob::Hit*>>;
    
    OldHitToNewHitMap   oldHitToNewHitMap;
    SpacePointHitVecMap spacePointHitVecMap;
    
    // We need a container for our new hits...
    std::unique_ptr<std::vector<recob::Hit>> newHitVecPtr(new std::vector<recob::Hit>);
    
    // reserve a chunk of memory... cannot be more hits than 3 x # spacer points...
    newHitVecPtr->reserve(3 * hitSpacePointAssnsHandle->size());
    
    // Use this handy art utility to make art::Ptr objects to the new recob::Hits for use in the output phase
    art::PtrMaker<recob::Hit> ptrMaker(evt);

    for(auto& assnPair : *hitSpacePointAssnsHandle)
    {
        const art::Ptr<recob::SpacePoint> spacePoint = assnPair.second;
        const art::Ptr<recob::Hit>&       recobHit   = assnPair.first;
        
        // If we have seen this hit before then no need to create new hit
        if (oldHitToNewHitMap.find(recobHit.get()) == oldHitToNewHitMap.end())
        {
            // Recover the reference WireID from our previous map
            geo::WireID refWireID  = spacePointToWireIDMap[spacePoint.get()];
            geo::WireID thisWireID = recobHit.get()->WireID();
            
            // Recover the list of possible WireIDs from the geometry service
            const std::vector<geo::WireID>& wireIDs = fGeometry->ChannelToWire(recobHit.get()->Channel());
            
            // Loop to find match
            for(const auto& wireID : wireIDs)
            {
                if (wireID.TPC != refWireID.TPC || wireID.Cryostat != refWireID.Cryostat) continue;
                thisWireID = wireID;
                break;
            }
            
            // Create and save the new recob::Hit with the correct WireID
            newHitVecPtr->emplace_back(recob::HitCreator(*recobHit.get(),thisWireID).copy());
            
            // Recover a pointer to it...
            recob::Hit* newHit = &newHitVecPtr->back();
            
            spacePointHitVecMap[spacePoint.get()].push_back(newHit);
            
            recobHitToArtPtrMap[newHit]       = ptrMaker(newHitVecPtr->size()-1);
            oldHitToNewHitMap[recobHit.get()] = newHit;
        }
        else spacePointHitVecMap[spacePoint.get()].push_back(oldHitToNewHitMap[recobHit.get()]);
    }
    
    // We'll want to correct the hit times for the plane offsets
    // (note this is already taken care of when converting to position)
    std::map<geo::PlaneID, double> planeOffsetMap;
    
    // Initialize the plane to hit vector map
    for(size_t cryoIdx = 0; cryoIdx < fGeometry->Ncryostats(); cryoIdx++)
    {
        for(size_t tpcIdx = 0; tpcIdx < fGeometry->NTPC(); tpcIdx++)
        {
            // What we want here are the relative offsets between the planes
            // Note that plane 0 is assumed the "first" plane and is the reference
            planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,0)] = 0.;
            planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,1)] = fDetector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,1))
                                                           - fDetector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0));
            planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,2)] = fDetector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,2))
                                                           - fDetector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0));
            
            std::cout << "***> plane 0 offset: " << planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,0)] << ", plane 1: " << planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,1)] << ", plane 2: " << planeOffsetMap[geo::PlaneID(cryoIdx,tpcIdx,2)] << std::endl;
            std::cout << "     Det prop plane 0: " << fDetector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0)) << ", plane 1: "  << fDetector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,1)) << ", plane 2: " << fDetector->GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,2)) << ", Trig: " << fDetector->TriggerOffset() << std::endl;
        }
    }

    // We need temporary mapping from recob::Hit's to our 2D hits
    using RecobHitTo2DHitMap = std::map<const recob::Hit*,const reco::ClusterHit2D*>;
    
    RecobHitTo2DHitMap recobHitTo2DHitMap;
    
    // Set the size of the container for our hits
    m_clusterHit2DMasterVec.clear();
    m_clusterHit2DMasterVec.reserve(oldHitToNewHitMap.size());

    // Now go throught the list of unique hits and create the 2D hits we'll use
    for(auto& hitPair : oldHitToNewHitMap)
    {
        const recob::Hit* recobHit = hitPair.second;
        const geo::WireID& hitWireID(recobHit->WireID());
        
        double hitPeakTime(recobHit->PeakTime() - planeOffsetMap.at(hitWireID.planeID())); //planeOffsetMap[hitWireID.planeID()]);
        double xPosition(fDetector->ConvertTicksToX(recobHit->PeakTime(), hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat));
        
        m_clusterHit2DMasterVec.emplace_back(0, 0., 0., xPosition, hitPeakTime, hitWireID, recobHit);
        
        recobHitTo2DHitMap[recobHit]  = &m_clusterHit2DMasterVec.back();
    }
    
    // Now we can go through the space points and build our 3D hits
    for(auto& pointPair : spacePointHitVecMap)
    {
        const recob::SpacePoint*              spacePoint  = pointPair.first;
        const std::vector<const recob::Hit*>& recobHitVec = pointPair.second;
        
        if (recobHitVec.size() != 3)
        {
            std::cout << "************>>>>>> do not have 3 hits associated to space point! " << recobHitVec.size() << " ***************" << std::endl;
            continue;
        }
        
        std::vector<const reco::ClusterHit2D*> hit2DVec(recobHitVec.size());
        
        for(const auto& recobHit : recobHitVec)
        {
            const reco::ClusterHit2D* hit2D = recobHitTo2DHitMap.at(recobHit);
            
            hit2DVec[hit2D->WireID().Plane] = hit2D;
        }
        
        // Weighted average, delta, sigmas, chisquare, kitchen sink, refrigerator for beer, etc.
        float avePeakTime(0.);
        float weightSum(0.);

        for(const auto& hit2D : hit2DVec)
        {
            float hitSigma = hit2D->getHit()->RMS();
            float weight   = 1. / (hitSigma * hitSigma);

            avePeakTime   += weight * hit2D->getTimeTicks();
            weightSum     += weight;
        }
        
        avePeakTime   /= weightSum;

        // Armed with the average peak time, now get hitChiSquare and the sig vec
        float              hitChiSquare(0.);
        float              sigmaPeakTime(std::sqrt(1./weightSum));
        
        for(const auto& hit2D : hit2DVec)
        {
            float hitRMS    = hit2D->getHit()->RMS();
            float combRMS   = std::sqrt(hitRMS*hitRMS - sigmaPeakTime*sigmaPeakTime);
            float peakTime  = hit2D->getTimeTicks();
            float deltaTime = peakTime - avePeakTime;
            float hitSig    = deltaTime / combRMS;
            
            hitChiSquare += hitSig * hitSig;
        }

        // The x position is a weighted sum but the y-z position is simply the average
        Eigen::Vector3f position(float(spacePoint->XYZ()[0]), float(spacePoint->XYZ()[1]), float(spacePoint->XYZ()[2]));
        float totalCharge = hit2DVec[0]->getHit()->Integral() + hit2DVec[1]->getHit()->Integral() + hit2DVec[2]->getHit()->Integral();
            
        reco::ClusterHit2DVec hitVector;
        
        hitVector.resize(3,NULL);
            
        // Make sure we have the hits
        hitVector.at(hit2DVec[0]->WireID().Plane) = hit2DVec[0];
        hitVector.at(hit2DVec[1]->WireID().Plane) = hit2DVec[1];
        hitVector.at(hit2DVec[2]->WireID().Plane) = hit2DVec[2];
            
        // And get the wire IDs
        std::vector<geo::WireID> wireIDVec = {geo::WireID(0,0,geo::kU,0), geo::WireID(0,0,geo::kV,0), geo::WireID(0,0,geo::kW,0)};
            
        for(const auto& hit : hitVector)
        {
            wireIDVec[hit->WireID().Plane] = hit->WireID();
                
            if (hit->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);
                
            hit->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);
        }
            
        unsigned int statusBits(0x7);
            
        // For compiling at the moment
        std::vector<float> hitDelTSigVec = {0.,0.,0.};

        hitDelTSigVec[0] = std::fabs(hitVector[0]->getTimeTicks() - 0.5 * (hitVector[1]->getTimeTicks() + hitVector[2]->getTimeTicks()));
        hitDelTSigVec[1] = std::fabs(hitVector[1]->getTimeTicks() - 0.5 * (hitVector[2]->getTimeTicks() + hitVector[0]->getTimeTicks()));
        hitDelTSigVec[2] = std::fabs(hitVector[2]->getTimeTicks() - 0.5 * (hitVector[0]->getTimeTicks() + hitVector[1]->getTimeTicks()));
        
        float deltaPeakTime = *std::min_element(hitDelTSigVec.begin(),hitDelTSigVec.end());

        // Create the 3D cluster hit
        hitPairList.emplace_back(0,
                                 statusBits,
                                 position,
                                 totalCharge,
                                 avePeakTime,
                                 deltaPeakTime,
                                 sigmaPeakTime,
                                 hitChiSquare,
                                 0.,
                                 0.,
                                 0.,
                                 hitVector,
                                 hitDelTSigVec,
                                 wireIDVec);
    }
    
    // Now we give the new hits to the refinery
    // Note that one advantage of using this utility is that it handles the
    // Hit/Wire and Hit/RawDigit associations all behind the scenes for us
    hitRefiner.use_hits(std::move(newHitVecPtr));

    // Output the new hit collection to the event
    hitRefiner.put_into();

    if (fEnableMonitoring)
    {
        theClockMakeHits.stop();
        
        fTimeVector[BUILDTHREEDHITS] = theClockMakeHits.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> 3D hit building done, found " << hitPairList.size() << " 3D Hits" << std::endl;

    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
    
DEFINE_ART_CLASS_TOOL(SpacePointHit3DBuilder)
} // namespace lar_cluster3d

/**
 *  @file   SpacePointHit3DBuilder_tool.cc
 *
 *  @brief  This tool provides "standard" 3D hits built (by this tool) from 2D hits
 *
 */

// Framework Includes
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IHit3DBuilder.h"

// std includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

// Ack!
#include "TTree.h"

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
    virtual void produces(art::ProducesCollector&) override;

    void configure(const fhicl::ParameterSet&) override;

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void Hit3DBuilder(art::Event &evt, reco::HitPairList& hitPairList, RecobHitToPtrMap&) override;

    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute(IHit3DBuilder::TimeValues index) const override {return fTimeVector.at(index);}

private:

    /**
     *  @brief clear the tuple vectors before processing next event
     */
    void clear();

    /**
     * @brief Perform charge integration between limits
     */
    float chargeIntegral(float,float,float,float,int,int) const;

    using Hit2DVector                 = std::vector<reco::ClusterHit2D>;

    /**
     *  @brief Data members to follow
     */
    art::InputTag                        fSpacePointProducerLabel;
    art::InputTag                        fHitProducerLabel;
    bool                                 fDoWireAssns;
    bool                                 fDoRawDigitAssns;
    float                                m_maxHit3DChiSquare;     ///< Provide ability to select hits based on "chi square"
    bool                                 m_outputHistograms;      ///< Take the time to create and fill some histograms for diagnostics

    bool                                 fEnableMonitoring;       ///<
    mutable std::vector<float>           fTimeVector;             ///<

    // Define some basic histograms
    TTree*                               m_tupleTree;             ///< output analysis tree

    mutable std::vector<float>           m_deltaTimeVec;
    mutable std::vector<float>           m_chiSquare3DVec;
    mutable std::vector<float>           m_maxPullVec;
    mutable std::vector<float>           m_overlapFractionVec;
    mutable std::vector<float>           m_overlapRangeVec;
    mutable std::vector<float>           m_maxSideVecVec;
    mutable std::vector<float>           m_pairWireDistVec;
    mutable std::vector<float>           m_smallChargeDiffVec;
    mutable std::vector<int>             m_smallIndexVec;
    mutable std::vector<float>           m_qualityMetricVec;
    mutable std::vector<float>           m_spacePointChargeVec;
    mutable std::vector<float>           m_hitAsymmetryVec;

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

void SpacePointHit3DBuilder::produces(art::ProducesCollector& collector)
{
    collector.produces< std::vector<recob::Hit>>();

    if (fDoWireAssns)     collector.produces< art::Assns<recob::Wire,   recob::Hit>>();
    if (fDoRawDigitAssns) collector.produces< art::Assns<raw::RawDigit, recob::Hit>>();

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
    m_maxHit3DChiSquare      = pset.get<float        >("MaxHitChiSquare",     6.0 );
    m_outputHistograms       = pset.get<bool         >("OutputHistograms",   false );

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;

    if (m_outputHistograms)
    {
        m_tupleTree = tfs->make<TTree>("Hit3DBuilderTree", "Tree by StandardHit3DBuilder");

        clear();

        m_tupleTree->Branch("DeltaTime2D",     "std::vector<float>", &m_deltaTimeVec);
        m_tupleTree->Branch("ChiSquare3D",     "std::vector<float>", &m_chiSquare3DVec);
        m_tupleTree->Branch("MaxPullValue",    "std::vector<float>", &m_maxPullVec);
        m_tupleTree->Branch("OverlapFraction", "std::vector<float>", &m_overlapFractionVec);
        m_tupleTree->Branch("OverlapRange",    "std::vector<float>", &m_overlapRangeVec);
        m_tupleTree->Branch("MaxSideVec",      "std::vector<float>", &m_maxSideVecVec);
        m_tupleTree->Branch("PairWireDistVec", "std::vector<float>", &m_pairWireDistVec);
        m_tupleTree->Branch("SmallChargeDiff", "std::vector<float>", &m_smallChargeDiffVec);
        m_tupleTree->Branch("SmallChargeIdx",  "std::vector<int>",   &m_smallIndexVec);
        m_tupleTree->Branch("QualityMetric",   "std::vector<float>", &m_qualityMetricVec);
        m_tupleTree->Branch("SPCharge",        "std::vector<float>", &m_spacePointChargeVec);
        m_tupleTree->Branch("HitAsymmetry",    "std::vector<float>", &m_hitAsymmetryVec);
    }

    art::ServiceHandle<geo::Geometry const> geometry;

    fGeometry = &*geometry;
    fDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

void SpacePointHit3DBuilder::clear()
{
    m_deltaTimeVec.clear();
    m_chiSquare3DVec.clear();
    m_maxPullVec.clear();
    m_overlapFractionVec.clear();
    m_overlapRangeVec.clear();
    m_maxSideVecVec.clear();
    m_pairWireDistVec.clear();
    m_smallChargeDiffVec.clear();
    m_smallIndexVec.clear();
    m_qualityMetricVec.clear();
    m_spacePointChargeVec.clear();
    m_hitAsymmetryVec.clear();

    return;
}

void SpacePointHit3DBuilder::Hit3DBuilder(art::Event& evt, reco::HitPairList& hitPairList, RecobHitToPtrMap& recobHitToArtPtrMap)
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
    recob::HitRefinerAssociator hitRefiner(evt, fHitProducerLabel, fDoWireAssns, fDoRawDigitAssns);

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

        reco::ClusterHit2DVec hitVector(recobHitVec.size());

        for(const auto& recobHit : recobHitVec)
        {
            const reco::ClusterHit2D* hit2D = recobHitTo2DHitMap.at(recobHit);

            hitVector[hit2D->WireID().Plane] = hit2D;
        }

        // Set up to get average peak time, hitChiSquare, etc.
        unsigned int statusBits(0x7);
        float        avePeakTime(0.);
        float        weightSum(0.);

        // And get the wire IDs
        std::vector<geo::WireID> wireIDVec = {geo::WireID(), geo::WireID(), geo::WireID()};

        // First loop through the hits to get WireIDs and calculate the averages
        for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
        {
            const reco::ClusterHit2D* hit2D = hitVector[planeIdx];

            wireIDVec[planeIdx] = hit2D->WireID();

            if (hit2D->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit2D->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);

            hit2D->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);

            float hitRMS   = hit2D->getHit()->RMS();
            float weight   = 1. / (hitRMS * hitRMS);
            float peakTime = hit2D->getTimeTicks();

            avePeakTime += peakTime * weight;
            weightSum   += weight;
        }

        avePeakTime /= weightSum;

        // Armed with the average peak time, now get hitChiSquare and the sig vec
        float              hitChiSquare(0.);
        float              sigmaPeakTime(std::sqrt(1./weightSum));
        std::vector<float> hitDelTSigVec;

        for(const auto& hit2D : hitVector)
        {
            float hitRMS    = hit2D->getHit()->RMS();
            float combRMS   = std::sqrt(hitRMS*hitRMS - sigmaPeakTime*sigmaPeakTime);
            float peakTime  = hit2D->getTimeTicks();
            float deltaTime = peakTime - avePeakTime;
            float hitSig    = deltaTime / combRMS;

            hitChiSquare += hitSig * hitSig;

            hitDelTSigVec.emplace_back(std::fabs(hitSig));
        }

        if (m_outputHistograms) m_chiSquare3DVec.push_back(hitChiSquare);

        // Need to determine the hit overlap ranges
        int lowMinIndex(std::numeric_limits<int>::max());
        int lowMaxIndex(std::numeric_limits<int>::min());
        int hiMinIndex(std::numeric_limits<int>::max());
        int hiMaxIndex(std::numeric_limits<int>::min());

        // This loop through hits to find min/max values for the common overlap region
        for(const auto& hit2D : hitVector)
        {
            int   hitStart = hit2D->getHit()->PeakTime() - 2. * hit2D->getHit()->RMS() - 0.5;
            int   hitStop  = hit2D->getHit()->PeakTime() + 2. * hit2D->getHit()->RMS() + 0.5;

            lowMinIndex = std::min(hitStart,    lowMinIndex);
            lowMaxIndex = std::max(hitStart,    lowMaxIndex);
            hiMinIndex  = std::min(hitStop + 1, hiMinIndex);
            hiMaxIndex  = std::max(hitStop + 1, hiMaxIndex);
        }

        // Keep only "good" hits...
        if (hitChiSquare < m_maxHit3DChiSquare && hiMinIndex > lowMaxIndex)
        {
            // One more pass through hits to get charge
            std::vector<float> chargeVec;

            for(const auto& hit2D : hitVector)
                chargeVec.push_back(chargeIntegral(hit2D->getHit()->PeakTime(),hit2D->getHit()->PeakAmplitude(),hit2D->getHit()->RMS(),1.,lowMaxIndex,hiMinIndex));

            float totalCharge     = std::accumulate(chargeVec.begin(),chargeVec.end(),0.) / float(chargeVec.size());
            float overlapRange    = float(hiMinIndex - lowMaxIndex);
            float overlapFraction = overlapRange / float(hiMaxIndex - lowMinIndex);

            // Set up to compute the charge asymmetry
            std::vector<float> smallestChargeDiffVec;
            std::vector<float> chargeAveVec;
            float              smallestDiff(std::numeric_limits<float>::max());
            size_t             chargeIndex(0);

            for(size_t idx = 0; idx < 3; idx++)
            {
                size_t leftIdx  = (idx + 2) % 3;
                size_t rightIdx = (idx + 1) % 3;

                smallestChargeDiffVec.push_back(std::abs(chargeVec[leftIdx] - chargeVec[rightIdx]));
                chargeAveVec.push_back(float(0.5 * (chargeVec[leftIdx] + chargeVec[rightIdx])));

                if (smallestChargeDiffVec.back() < smallestDiff)
                {
                    smallestDiff = smallestChargeDiffVec.back();
                    chargeIndex  = idx;
                }

                // Take opportunity to look at peak time diff
                if (m_outputHistograms)
                {
                    float deltaPeakTime = hitVector[leftIdx]->getTimeTicks() - hitVector[rightIdx]->getTimeTicks();

                    m_deltaTimeVec.push_back(deltaPeakTime);
                }
            }

            float chargeAsymmetry = (chargeAveVec[chargeIndex] - chargeVec[chargeIndex]) / (chargeAveVec[chargeIndex] + chargeVec[chargeIndex]);

            // If this is true there has to be a negative charge that snuck in somehow
            if (chargeAsymmetry < -1. || chargeAsymmetry > 1.)
            {
                const geo::WireID& hitWireID = hitVector[chargeIndex]->WireID();

                std::cout << "============> Charge asymmetry out of range: " << chargeAsymmetry << " <============" << std::endl;
                std::cout << "     hit C: " << hitWireID.Cryostat << ", TPC: " << hitWireID.TPC << ", Plane: " << hitWireID.Plane << ", Wire: " << hitWireID.Wire << std::endl;
                std::cout << "     charge: " << chargeVec[0] << ", " << chargeVec[1] << ", " << chargeVec[2] << std::endl;
                std::cout << "     index: " << chargeIndex << ", smallest diff: " << smallestDiff << std::endl;
                continue;
            }

            // Usurping "deltaPeakTime" to be the maximum pull
            float deltaPeakTime = *std::max_element(hitDelTSigVec.begin(),hitDelTSigVec.end());

            if (m_outputHistograms)
            {
                m_smallChargeDiffVec.push_back(smallestDiff);
                m_smallIndexVec.push_back(chargeIndex);
                m_maxPullVec.push_back(deltaPeakTime);
                m_qualityMetricVec.push_back(hitChiSquare);
                m_spacePointChargeVec.push_back(totalCharge);
                m_overlapFractionVec.push_back(overlapFraction);
                m_overlapRangeVec.push_back(overlapRange);
                m_hitAsymmetryVec.push_back(chargeAsymmetry);
            }

            Eigen::Vector3f position(float(spacePoint->XYZ()[0]), float(spacePoint->XYZ()[1]), float(spacePoint->XYZ()[2]));

            // Create the 3D cluster hit
            hitPairList.emplace_back(0,
                                     statusBits,
                                     position,
                                     totalCharge,
                                     avePeakTime,
                                     deltaPeakTime,
                                     sigmaPeakTime,
                                     hitChiSquare,
                                     overlapFraction,
                                     chargeAsymmetry,
                                     0.,
                                     0.,
                                     hitVector,
                                     hitDelTSigVec,
                                     wireIDVec);
        }
    }

    // Now we give the new hits to the refinery
    // Note that one advantage of using this utility is that it handles the
    // Hit/Wire and Hit/RawDigit associations all behind the scenes for us
    hitRefiner.use_hits(std::move(newHitVecPtr));

    // Output the new hit collection to the event
    hitRefiner.put_into();

    // Handle tree output too
    m_tupleTree->Fill();

    clear();

    if (fEnableMonitoring)
    {
        theClockMakeHits.stop();

        fTimeVector[BUILDTHREEDHITS] = theClockMakeHits.accumulated_real_time();
    }

    mf::LogDebug("Cluster3D") << ">>>>> 3D hit building done, found " << hitPairList.size() << " 3D Hits" << std::endl;

    return;
}

float SpacePointHit3DBuilder::chargeIntegral(float peakMean,
                                             float peakAmp,
                                             float peakSigma,
                                             float areaNorm,
                                             int   low,
                                             int   hi) const
{
    float integral(0);

    for(int sigPos = low; sigPos < hi; sigPos++)
    {
        float arg = (float(sigPos) - peakMean + 0.5) / peakSigma;
        integral += peakAmp * std::exp(-0.5 * arg * arg);
    }

    return integral;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DEFINE_ART_CLASS_TOOL(SpacePointHit3DBuilder)
} // namespace lar_cluster3d

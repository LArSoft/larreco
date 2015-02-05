/**
 *  @file   Cluster3D_module.cc
 *
 *          Class:       Cluster3D
 *          Module Type: Producer
 * 
 *  @brief  Producer module to create 3D clusters from input recob::Hit objects
 *
 *          This producer module will drive the 3D association of recob::Hit objects
 *          to form 3D clusters. This information will be output as:
 *          1) a PFParticle to anchor all the other objects (as associations)
 *          2) three recob::Cluster objects representing 2D hit clusterswith associations 
 *             to the 2D hits comprising each of them
 *          3) One or more recob::PCAxis objects representing the Principal Components
 *             Analysis output of the space points associated to the 3D objects
 *          4) recob::SpacePoints representing the accepted 3D points for each PFParticle
 *          5) recob::Seed objects and associated seed hits representing candidate straight
 *             line segments in the space point collection. 
 *
 *          The module has two main sections
 *          1) Find the 3D clusters of 3D hits
 *          2) Get the output objects for each of the 3D clusters
 *          See the code below for more detail on these steps
 *
 *          Note that the general 3D cluster suite of algorithms make extensive use of a set of data objects 
 *          which contain volatile data members. At the end of the routine these are used to make the output
 *          LArSoft data products described above. See LarData/RecoObjects/Cluster3D.h
 *   
 *          Configuration parameters:
 *          HitFinderModuleLabel:         the producer module responsible for making the recob:Hits to use
 *          EnableMonitoring:             if true then basic monitoring of the module performed
 *          DBScanAlg:                    Parameter block required by the 3D clustering algorithm
 *          PrincipalComponentsAlg:       Parameter block required by the Principal Components Analysis Algorithm
 *          SkeletonAlg:                  Parameter block required by the 3D skeletonization algorithm
 *          SeedFinderAlg:                Parameter block required by the Hough Seed Finder algorithm
 *          PCASeedFinderAlg:             Parameter block required by the PCA Seed Finder algorithm
 *          ParrallelHitsAlg:             Parameter block required by the parallel hits algorithm
 *
 *          The current producer module does not try to analyze or break apart PFParticles
 *          so, for example, all tracks emanating from a common vertex will be associated
 *          to a single PFParticle
 *
 *  @author usher@slac.stanford.edu 
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCTruth.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/PCAxis.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Seed.h"
#include "RecoObjects/Cluster3D.h"
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "RecoAlg/Cluster3DAlgs/HoughSeedFinderAlg.h"
#include "RecoAlg/Cluster3DAlgs/PCASeedFinderAlg.h"
#include "RecoAlg/Cluster3DAlgs/ParallelHitsSeedFinderAlg.h"
#include "RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "RecoAlg/Cluster3DAlgs/SkeletonAlg.h"
#include "RecoAlg/Cluster3DAlgs/DBScanAlg.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterRecoUtil/OverriddenClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"
#include "ClusterFinder/ClusterCreator.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
//------------------------------------------------------------------------------------------------------------------------------------------
// Start with typedefs and definitions of some utility classes

typedef std::map<const recob::Hit*, art::Ptr<recob::Hit> > RecobHitToPtrMap;
typedef std::vector<art::Ptr<recob::Hit> >                 RecobHitVector;
    
typedef std::pair<reco::PrincipalComponents, reco::HitPairClusterMap::iterator> PCAHitPairClusterMapPair;
    
/**
 *  @brief A utility class used in construction of 3D clusters
 */
class RecobClusterParameters
{
public:
    RecobClusterParameters() : m_startTime(999999.),
        m_sigmaStartTime(1.),
        m_endTime(0.),
        m_sigmaEndTime(1.),
        m_startWire(9999999),
        m_endWire(0),
        m_view(geo::kUnknown)
    {
        m_hitVector.clear();
    }
        
    void UpdateParameters(const reco::ClusterHit2D* hit);
        
    double         m_startTime;
    double         m_sigmaStartTime;
    double         m_endTime;
    double         m_sigmaEndTime;
    unsigned int   m_startWire;
    unsigned int   m_endWire;
    geo::View_t    m_view;
    HitVectorConst m_hitVector;
};
    
typedef std::map<geo::View_t, RecobClusterParameters> ViewToClusterParamsMap;

/**
 *  @brief Class wrapping the above and containing volatile information to characterize the cluster
 */
class ClusterParameters
{
public:
    ClusterParameters(reco::HitPairClusterMap::iterator& mapItr) : m_hitPairClusterMapItr(mapItr)
    {
        m_clusterParams.clear();
    }
    
    void UpdateParameters(const reco::ClusterHit2D* hit)
    {
        m_clusterParams[hit->getHit().View()].UpdateParameters(hit);
    }
    
    ViewToClusterParamsMap            m_clusterParams;
    reco::HitPairClusterMap::iterator m_hitPairClusterMapItr;
    reco::PrincipalComponents         m_fullPCA;
    reco::PrincipalComponents         m_skeletonPCA;
};
    
typedef std::list<ClusterParameters> ClusterParametersList;

//------------------------------------------------------------------------------------------------------------------------------------------
// Definition of the producer module here
    
/**
 *  @brief  Definition of the Cluster3D class
 */
class Cluster3D : public art::EDProducer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset - reference to the parameters used by this module and its algorithms
     */
    Cluster3D(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~Cluster3D();

    /**
     *  @brief declare the standard art functions that we'll implement in this producer module
     */
    void beginJob();
    void endJob();
    void produce(art::Event &evt);
    void reconfigure(fhicl::ParameterSet const &pset);

private:

    /**
     *  @brief  Event Preparation 
     * 
     *  @param  evt  the ART event 
     */
    void PrepareEvent(const art::Event &evt);  

    /**
     *  @brief  Extract the ART hits and the ART hit-particle relationships
     * 
     *  @param  evt                   the ART event
     *  @param  hit2DVector           A container for the internal Cluster3D 2D hit objects
     *  @param  viewToHitVectorMap    A map between view and the internal Cluster3D 2D hit objects
     *  @param  viewToWireToHitSetMap This maps 2D hits to wires and stores by view
     */
    void CollectArtHits(art::Event &evt, Hit2DVector& hit2DVector, ViewToHitVectorMap& viewToHitVector, ViewToWireToHitSetMap& viewToWireToHitSetMap);

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();
    
    /**
     *  @brief Given the results of running DBScan, format the clusters so that they can be 
     *         easily transferred back to the larsoft world
     *
     *  @param hitPairClusterMap      map between view and a list of 3D hits
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     *  @param doca3DScl              Gives a maximum 3D doca to the PCA axis cut
     *
     *                                The last two parameters are passed through to the FillClusterParams method
     */
    void BuildClusterInfo(reco::HitPairClusterMap& hitPairClusterMap, ClusterParametersList& clusterParametersList, double rejectFraction = 0.5, double doca3DScl = 5.);
    
    /**
     *  @brief A generic routine to actually fill the clusterParams
     *
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     *  @param doca3DScl              Gives a maximum 3D doca to the PCA axis cut
     */
    void FillClusterParams(ClusterParameters& clusterParams, double rejectFraction = 0.5, double doca3DScl = 3.);
    
    /**
     *  @brief An interface to the seed finding algorithm
     *
     *  @param evt          the ART event
     *  @param cluster      structure of information representing a single cluster
     *  @param seedVec      the output vector of candidate seeds
     *  @param seedHitAssns the associations between the seeds and the 2D hits making them
     */
    void findTrackSeeds(art::Event&                         evt,
                        ClusterParameters&                  cluster,
                        std::vector<recob::Seed>&           seedVec,
                        art::Assns<recob::Seed,recob::Hit>& seedHitAssns);
    /**
     *  @brief Produces the art output from all the work done in this producer module
     *
     *  @param evt                   the ART event
     *  @param hitPairList           List of all 3D Hits in internal Cluster3D format
     *  @param clusterParametersList Data structure containing the cluster information to output
     */
    void ProduceArtClusters(art::Event &evt, HitPairList& hitPairList, ClusterParametersList& clusterParametersList);

    /**
     *   Algorithm parameters
     */
    bool                      m_enableMonitoring;      ///<
    std::string               m_hitfinderModuleLabel;  ///<

    /**
     *   Tree variables for output
     */
    TTree*                    m_pRecoTree;             ///<
    int                       m_run;                   ///<
    int                       m_event;                 ///<
    int                       m_hits;                  ///< Keeps track of the number of hits seen
    float                     m_totalTime;             ///< Keeps track of total execution time
    float                     m_artHitsTime;           ///< Keeps track of time to recover hits
    float                     m_makeHitsTime;          ///< Keeps track of time to build 3D hits
    float                     m_buildNeighborhoodTime; ///< Keeps track of time to build epsilon neighborhood
    float                     m_dbscanTime;            ///< Keeps track of time to run DBScan
    float                     m_finishTime;            ///< Keeps track of time to run output module
    
    /** 
     *   Other useful variables
     */
    RecobHitToPtrMap          m_hitToPtrMap;           ///<  Mapping for recovering art::Ptr to hits

    geo::Geometry*            m_geometry;              ///<  pointer to the Geometry service
    util::DetectorProperties* m_detector;              ///<  Pointer to the detector properties
    
    DBScanAlg                 m_dbScanAlg;             ///<  Algorithm to cluster hits
    PrincipalComponentsAlg    m_pcaAlg;                ///<  Principal Components algorithm
    SkeletonAlg               m_skeletonAlg;           ///<  Skeleton point finder
    HoughSeedFinderAlg        m_seedFinderAlg;         ///<  Seed finder
    PCASeedFinderAlg          m_pcaSeedFinderAlg;      ///<  Use PCA axis to find seeds
    ParallelHitsSeedFinderAlg m_parallelHitsAlg;       ///<  Deal with parallel hits clusters
};

DEFINE_ART_MODULE(Cluster3D)

} // namespace lar_cluster3d

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

Cluster3D::Cluster3D(fhicl::ParameterSet const &pset) :
    m_dbScanAlg(pset.get<fhicl::ParameterSet>("DBScanAlg")),
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg")),
    m_skeletonAlg(pset.get<fhicl::ParameterSet>("SkeletonAlg")),
    m_seedFinderAlg(pset.get<fhicl::ParameterSet>("SeedFinderAlg")),
    m_pcaSeedFinderAlg(pset.get<fhicl::ParameterSet>("PCASeedFinderAlg")),
    m_parallelHitsAlg(pset.get<fhicl::ParameterSet>("ParallelHitsAlg"))
{
    this->reconfigure(pset);

    produces< std::vector<recob::PCAxis> >();
    produces< std::vector<recob::PFParticle> >();
    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::SpacePoint> >();
    produces< std::vector<recob::Seed> >();
    produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();
    produces< art::Assns<recob::Seed,       recob::Hit> >();
    produces< art::Assns<recob::Cluster,    recob::Hit> >();
    produces< art::Assns<recob::SpacePoint, recob::Hit> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

Cluster3D::~Cluster3D()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::reconfigure(fhicl::ParameterSet const &pset)
{
    m_hitfinderModuleLabel = pset.get<std::string>("HitFinderModuleLabel","gaushit");
    m_enableMonitoring     = pset.get<bool>("EnableMonitoring",false);
    
    m_dbScanAlg.reconfigure(pset.get<fhicl::ParameterSet>("DBScanAlg"));
    m_pcaAlg.reconfigure(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));
    m_skeletonAlg.reconfigure(pset.get<fhicl::ParameterSet>("SkeletonAlg"));
    m_seedFinderAlg.reconfigure(pset.get<fhicl::ParameterSet>("SeedFinderAlg"));
    m_pcaSeedFinderAlg.reconfigure(pset.get<fhicl::ParameterSet>("PCASeedFinderAlg"));
    m_parallelHitsAlg.reconfigure(pset.get<fhicl::ParameterSet>("ParallelHitsAlg"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::beginJob()
{
    /**
     *  @brief beginJob will be tasked with initializing monitoring, in necessary, but also to init the 
     *         geometry and detector services (and this probably needs to go in a "beginEvent" method?)
     */
    if (m_enableMonitoring)
        this->InitializeMonitoring();
    
    art::ServiceHandle<geo::Geometry>            geometry;
    art::ServiceHandle<util::DetectorProperties> detectorProperties;
    
    m_geometry = &*geometry;
    m_detector = &*detectorProperties;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::produce(art::Event &evt)
{
    /**
     *  @brief Producer method for reovering the 2D hits and driving the 3D reconstruciton
     */
    mf::LogInfo("Cluster3D") << " *** Cluster3D::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] Starting Now! *** " << std::endl;

    // Set up for monitoring the timing... at some point this should be removed in favor of
    // external profilers
    cet::cpu_timer theClockTotal;
    cet::cpu_timer theClockArtHits;
    cet::cpu_timer theClockFinish;
    
    if (m_enableMonitoring)
    {
        theClockTotal.start();
        theClockArtHits.start();
    }
    
    // This really only does anything if we are monitoring since it clears our tree variables
    this->PrepareEvent(evt);

    // Get instances of the primary data structures needed
    Hit2DVector                  clusterHit2DMasterVec;
    ViewToHitVectorMap           viewToHitVectorMap;
    ViewToWireToHitSetMap        viewToWireToHitSetMap;
    reco::HitPairClusterMap      hitPairClusterMap;
    ClusterParametersList        clusterParametersList;
    std::auto_ptr< HitPairList > hitPairList(new HitPairList); // Potentially lots of hits, use heap instead of stack
    
    // Paranoia... clear clear cleeeaaar
    m_hitToPtrMap.clear();
    
    // Recover the 2D hits and then organize them into data structures which will be used in the
    // DBscan algorithm for building the 3D clusters
    this->CollectArtHits(evt, clusterHit2DMasterVec, viewToHitVectorMap, viewToWireToHitSetMap);
    
    if (m_enableMonitoring) theClockArtHits.stop();
    
    // If there are no hits in our view/wire data structure then do not proceed with the full analysis
    if (!viewToWireToHitSetMap.empty())
    {
        // Call the main workhorse algorithm for building the local version of candidate 3D clusters
        m_dbScanAlg.ClusterHitsDBScan(viewToHitVectorMap, viewToWireToHitSetMap, *hitPairList, hitPairClusterMap);
        
        // Given the work above, process and build the list of 3D clusters to output
        BuildClusterInfo(hitPairClusterMap, clusterParametersList);
    }
    
    if(m_enableMonitoring) theClockFinish.start();

    // Call the module that does the end processing (of which there is quite a bit of work!)
    // This goes here to insure that something is always written to the data store
    ProduceArtClusters(evt, *hitPairList, clusterParametersList);
    
    if (m_enableMonitoring) theClockFinish.stop();
    
    // We are done with this map
    m_hitToPtrMap.clear();
    
    // If monitoring then deal with the fallout
    if (m_enableMonitoring)
    {
        theClockTotal.stop();

        m_run                   = evt.run();
        m_event                 = evt.id().event();
        m_totalTime             = theClockTotal.accumulated_real_time();
        m_artHitsTime           = theClockArtHits.accumulated_real_time();
        m_makeHitsTime          = m_dbScanAlg.getTimeToExecute(DBScanAlg::BUILDTHREEDHITS);
        m_buildNeighborhoodTime = m_dbScanAlg.getTimeToExecute(DBScanAlg::BUILDHITTOHITMAP);
        m_dbscanTime            = m_dbScanAlg.getTimeToExecute(DBScanAlg::RUNDBSCAN);
        m_finishTime            = theClockFinish.accumulated_real_time();
        m_hits                  = static_cast<int>(clusterHit2DMasterVec.size());
        m_pRecoTree->Fill();
        
        mf::LogDebug("Cluster3D") << "*** Cluster3D total time: " << m_totalTime << ", art: " << m_artHitsTime << ", make: " << m_makeHitsTime
        << ", build: " << m_buildNeighborhoodTime << ", dbscan: " << m_dbscanTime << ", finish: " << m_finishTime << std::endl;
    }
    
    // Will we ever get here? ;-)
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run",                  &m_run,                   "run/I");
    m_pRecoTree->Branch("event",                &m_event,                 "event/I");
    m_pRecoTree->Branch("hits",                 &m_hits,                  "hits/I");
    m_pRecoTree->Branch("totalTime",            &m_totalTime,             "time/F");
    m_pRecoTree->Branch("artHitsTime",          &m_artHitsTime,           "time/F");
    m_pRecoTree->Branch("makeHitsTime",         &m_makeHitsTime,          "time/F");
    m_pRecoTree->Branch("buildneigborhoodTime", &m_buildNeighborhoodTime, "time/F");
    m_pRecoTree->Branch("dbscanTime",           &m_dbscanTime,            "time/F");
    m_pRecoTree->Branch("finishTime",           &m_finishTime,            "time/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::PrepareEvent(const art::Event &evt)
{
    m_run                   = evt.run();
    m_event                 = evt.id().event();
    m_hits                  = 0;
    m_totalTime             = 0.f;
    m_artHitsTime           = 0.f;
    m_makeHitsTime          = 0.f;
    m_buildNeighborhoodTime = 0.f;
    m_dbscanTime            = 0.f;
    m_finishTime            = 0.f;
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
void Cluster3D::CollectArtHits(art::Event&            evt,
                               Hit2DVector&           hitVector,
                               ViewToHitVectorMap&    viewToHitVectorMap,
                               ViewToWireToHitSetMap& viewToWireToHitSetMap)
{
    /**
     *  @brief Recover the 2D hits from art and fill out the local data structures for the 3D clustering
     */
    art::Handle< std::vector<recob::Hit> > recobHitHandle;
    evt.getByLabel(m_hitfinderModuleLabel, recobHitHandle);
    
    if (!recobHitHandle.isValid()) return;
    
    // We'll need the offsets for each plane
    std::map<geo::View_t, double> viewOffsetMap;
    
    viewOffsetMap[geo::kU] = m_detector->GetXTicksOffset(geo::kU, 0, 0)-m_detector->TriggerOffset();
    viewOffsetMap[geo::kV] = m_detector->GetXTicksOffset(geo::kV, 0, 0)-m_detector->TriggerOffset();
    viewOffsetMap[geo::kW] = m_detector->GetXTicksOffset(geo::kW, 0, 0)-m_detector->TriggerOffset();
    
    // Reserve memory for the hit vector
    hitVector.reserve(recobHitHandle->size());
    
    // Cycle through the recob hits to build ClusterHit2D objects and insert
    // them into the map
    for (size_t cIdx = 0; cIdx < recobHitHandle->size(); cIdx++)
    {
        art::Ptr<recob::Hit> recobHit(recobHitHandle, cIdx);
        
        const geo::WireID& hitWireID(recobHit->WireID());
        
        double hitPeakTime(recobHit->PeakTime() - viewOffsetMap[recobHit->View()]);
        double xPosition(m_detector->ConvertTicksToX(recobHit->PeakTime(), hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat));
        
        hitVector.emplace_back(reco::ClusterHit2D(0, 0., 0., xPosition, hitPeakTime, *recobHit));

        viewToHitVectorMap[recobHit->View()].push_back(&hitVector.back());
        viewToWireToHitSetMap[recobHit->View()][recobHit->WireID().Wire].insert(&hitVector.back());
        
        const recob::Hit* recobHitPtr = recobHit.get();
        m_hitToPtrMap[recobHitPtr]    = recobHit;
    }
    
    // Make a loop through to sort the recover hits in time order
    std::sort(viewToHitVectorMap[geo::kU].begin(), viewToHitVectorMap[geo::kU].end(), SetHitTimeOrder);
    std::sort(viewToHitVectorMap[geo::kV].begin(), viewToHitVectorMap[geo::kV].end(), SetHitTimeOrder);
    std::sort(viewToHitVectorMap[geo::kW].begin(), viewToHitVectorMap[geo::kW].end(), SetHitTimeOrder);

    mf::LogDebug("Cluster3D") << ">>>>> Number of ART hits: " << recobHitHandle->size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void RecobClusterParameters::UpdateParameters(const reco::ClusterHit2D* clusterHit)
{
    /**
     *  @brief a utility routine for building 3D clusters to keep basic info up to date
     *         (a candidate for a better way to do this)
     */
    const recob::Hit& hit = clusterHit->getHit();
    
    // Need to keep track of stuff so we can form cluster
    if (hit.WireID().Wire < m_startWire)
    {
        m_startWire      = hit.WireID().Wire;
        m_startTime      = hit.PeakTimeMinusRMS();
        m_sigmaStartTime = hit.SigmaPeakTime();
    }
    
    if (hit.WireID().Wire > m_endWire)
    {
        m_endWire      = hit.WireID().Wire;
        m_endTime      = hit.PeakTimePlusRMS();
        m_sigmaEndTime = hit.SigmaPeakTime();
    }
    
    m_view         = hit.View();
    
    m_hitVector.push_back(clusterHit);
    
    return;
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
    
void Cluster3D::BuildClusterInfo(reco::HitPairClusterMap& hitPairClusterMap, ClusterParametersList& clusterParametersList, double rejectFraction, double doca3DScl)
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
    
    // This is a remote possibility but why not check?
    if (!hitPairClusterMap.empty())
    {
        size_t minHitsPerCluster(3);
        
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
            if (mapItr->second.size() < minHitsPerCluster) continue;
            
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
            clusterParametersList.push_back(ClusterParameters(hitPairClusterMapItr));
            
            // Can we get a reference to what we just created?
            ClusterParameters& clusterParams = clusterParametersList.back();
            
            // Do the actual work of filling the parameters
            FillClusterParams(clusterParams, rejectFraction, doca3DScl);
            
            // If this cluster is rejected then the parameters will be empty
            if (clusterParams.m_clusterParams.empty() || !clusterParams.m_fullPCA.getSvdOK())
            {
                clusterParametersList.pop_back();
            }
        }
    }
    
    return;
}
    
void Cluster3D::FillClusterParams(ClusterParameters& clusterParams, double rejectFraction, double doca3DScl)
{
    /**
     *  @brief Given a list of hits fill out the remaining parameters for this cluster and evaluate the
     *         candidate's worthiness to achieve stardom in the event display
     */
    
    // Recover the HitPairListPtr from the input clusterParams (which will be the
    // only thing that has been provided)
    reco::HitPairListPtr& hitPairVector = clusterParams.m_hitPairClusterMapItr->second;
    
    // To be sure, we should clear the other data members
    clusterParams.m_clusterParams.clear();
    clusterParams.m_fullPCA = reco::PrincipalComponents();
    
    // See if we can avoid duplicates by temporarily transferring to a set
    std::set<const reco::ClusterHit2D*> hitSet;
    
    int numTotal(0);
    int numUniqueHits(0);
    int numSharedHits(0);
    int numLostHits(0);
    
    // First loop through the 3D hits
    // The goal of this loop is to build a set of unique hits from the hit pairs (which may contain many
    // ambiguous duplicate combinations).
    // The secondary goal is to remove 3D hits marked by hit arbitration to be tossed
    for(const auto& hit3D : hitPairVector)
    {
        int nHitsAlreadyUsed(0);
        
        numTotal += hit3D->getHits().size();
        
        // loop over the hits in this 3D Cluster hit
        for(const auto& hit2D : hit3D->getHits())
        {
            if (hit2D->getStatusBits() & 0x1) nHitsAlreadyUsed++;
        }
        
        numSharedHits += nHitsAlreadyUsed;
        
        if (nHitsAlreadyUsed < 2) //int(hit3D->getHits().size()))
        {
            numUniqueHits += int(hit3D->getHits().size()) - nHitsAlreadyUsed;
            for(const auto& hit2D : hit3D->getHits()) hitSet.insert(hit2D);
        }
        else
        {
            numLostHits += hit3D->getHits().size();
            hit3D->setStatusBit(0x80000000);
        }
    }
    
    // If we have something left then at this point we make one more check
    // This check is intended to weed out clusters made from isolated groups of ambiguous hits which
    // really belong to a larger cluster
    if (numUniqueHits > 3)
    {
        // Look at reject to accept ratio
        //double rejectToAccept = double(numRejected) / double(numAccepted);
        double acceptRatio = double(numUniqueHits) / double(numTotal);
        double lostRatio   = double(numLostHits)   / double(numTotal);
        
        // Arbitrary rejection criteria... need to understand
        // Anyway, if we get past this we're making a cluster
        //if (rejectToAccept < rejectFraction)
        if (acceptRatio > rejectFraction && lostRatio < 1.)
        {
            // Why do we need to explicitly define this?
            unsigned int usedBit(0x1);
            
            // Add the "good" hits to our cluster parameters
            for(const auto& hit2D : hitSet)
            {
                hit2D->setStatusBit(usedBit);
                clusterParams.UpdateParameters(hit2D);
            }
            
            // Final selection cut, need at least 3 hits each view
            if (clusterParams.m_clusterParams[geo::kU].m_hitVector.size() > 3 &&
                clusterParams.m_clusterParams[geo::kV].m_hitVector.size() > 3 &&
                clusterParams.m_clusterParams[geo::kW].m_hitVector.size() > 3)
            {
                // First stage of feature extraction runs here
                //PCAAnalysis(clusterParams, doca3DScl);
                m_pcaAlg.PCAAnalysis_3D(clusterParams.m_hitPairClusterMapItr->second, clusterParams.m_fullPCA);
                
                // Set the skeleton PCA to make sure it has some value
                clusterParams.m_skeletonPCA = clusterParams.m_fullPCA;
            }
        }
    }

    return;
}
    
void Cluster3D::findTrackSeeds(art::Event&                         evt,
                               ClusterParameters&                  cluster,
                               std::vector<recob::Seed>&           seedVec,
                               art::Assns<recob::Seed,recob::Hit>& seedHitAssns)
{
    /**
     *  @brief This method provides an interface to various algorithms for finding candiate 
     *         recob::Seed objects and, as well, their candidate related seed hits
     */
    
    // Make sure we are using the right pca
    reco::PrincipalComponents& fullPCA        = cluster.m_fullPCA;
    reco::PrincipalComponents& skeletonPCA    = cluster.m_skeletonPCA;
    reco::HitPairListPtr       hitPairListPtr = cluster.m_hitPairClusterMapItr->second;
    reco::HitPairListPtr       skeletonListPtr;

    // We want to work with the "skeleton" hits so first step is to call the algorithm to
    // recover only these hits from the entire input collection
    m_skeletonAlg.GetSkeletonHits(hitPairListPtr, skeletonListPtr);
    
    // Skeleton hits are nice but we can do better if we then make a pass through to "average"
    // the skeleton hits position in the Y-Z plane
    m_skeletonAlg.AverageSkeletonPositions(skeletonListPtr);
    
    SeedHitPairListPairVec seedHitPairVec;
    
    // Some combination of the elements below will be used to determine which seed finding algorithm
    // to pursue below
    double eigenVal0 = 3. * sqrt(skeletonPCA.getEigenValues()[0]);
    double eigenVal1 = 3. * sqrt(skeletonPCA.getEigenValues()[1]);
    double eigenVal2 = 3. * sqrt(skeletonPCA.getEigenValues()[2]);
    double transRMS  = sqrt(std::pow(eigenVal1,2) + std::pow(eigenVal2,2));
    
    bool   foundGoodSeed(false);

    // Choose a method for finding the seeds based on the PCA that was run...
    // Currently we have an ad hoc if-else block which I hope will be improved soon!
    if (fabs(fullPCA.getEigenVectors()[2][0]) > 0.999 && 3. * sqrt(fullPCA.getEigenValues()[1]) > 20.)
    {
        // In this case we have a track moving relatively parallel to the wire plane with lots of
        // ambiguous 3D hits. Your best bet here is to use the "parallel hits" algorithm to get the
        // best axis and seeds
        // This algorithm does not fail (foundGoodSeed will always return true)
        foundGoodSeed = m_parallelHitsAlg.findTrackSeeds(skeletonListPtr, skeletonPCA, seedHitPairVec);
    }
    else if (eigenVal0 > 40. && transRMS < 5.)
    {
        // If the input cluster is relatively "straight" then chances are it is a single straight track,
        // probably a CR muon, and we can simply use the PCA to determine the seed
        // This algorithm will check both "ends" of the input hits and if the angles become inconsistent
        // then it will "fail"
        foundGoodSeed = m_pcaSeedFinderAlg.findTrackSeeds(skeletonListPtr, skeletonPCA, seedHitPairVec);
    }
    
    // In the event the above two methods failed then we hit it with the real seed finder
    if (!foundGoodSeed)
    {
        // If here then we have a complicated 3D cluster and we'll use the hough transform algorithm to
        // return a list of candidate seeds and seed hits
        m_seedFinderAlg.findTrackSeeds(skeletonListPtr, skeletonPCA, seedHitPairVec);
    }

    // Go through the returned lists and build out the art friendly seeds and hits
    for(const auto& seedHitPair : seedHitPairVec)
    {
        seedVec.push_back(seedHitPair.first);
        
        // We use a set here because our 3D hits can share 2D hits
        // The set will make sure we get unique combinations of 2D hits
        std::set<art::Ptr<recob::Hit> > seedHitSet;
        
        for(const auto& hit3D : seedHitPair.second)
        {
            for(const auto& hit2D : hit3D->getHits())
            {
                const recob::Hit* recobHit = &hit2D->getHit();
                
                seedHitSet.insert(m_hitToPtrMap[recobHit]);
            }
        }
        
        RecobHitVector seedHitVec;
        
        for(const auto& hit2D : seedHitSet) seedHitVec.push_back(hit2D);
        
        util::CreateAssn(*this, evt, seedVec, seedHitVec, seedHitAssns);
    }
    
    return;
}

struct ClusterOrder
{
    bool operator()(HitClusterMap::const_iterator left, HitClusterMap::const_iterator right)
    {
        return left->second.size() > right->second.size();
    }
};
    
void Cluster3D::ProduceArtClusters(art::Event &evt, HitPairList& hitPairVector, ClusterParametersList& clusterParametersList)
{
    /**
     *  @brief The workhorse to take the candidate 3D clusters and produce all of the necessary art output
     */
    
    mf::LogDebug("Cluster3D") << " *** Cluster3D::ProduceArtClusters() *** " << std::endl;
    
    std::unique_ptr< std::vector<recob::PCAxis> >     artPCAxisVector( new std::vector<recob::PCAxis>         );
    std::unique_ptr< std::vector<recob::PFParticle> > artPFParticleVector( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::Cluster> >    artClusterVector( new std::vector<recob::Cluster>       );
    std::unique_ptr< std::vector<recob::SpacePoint> > artSpacePointVector( new std::vector<recob::SpacePoint> );
    std::unique_ptr< std::vector<recob::Seed> >       artSeedVector( new std::vector<recob::Seed>             );
    
    std::unique_ptr< art::Assns<recob::Cluster,    recob::Hit> >         artClusterAssociations(    new art::Assns<recob::Cluster,    recob::Hit>          );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis> >      artPFPartAxisAssociations( new art::Assns<recob::PFParticle, recob::PCAxis>       );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >     artPFPartClusAssociations( new art::Assns<recob::PFParticle, recob::Cluster>      );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> >  artPFPartSPAssociations(   new art::Assns<recob::PFParticle, recob::SpacePoint>   );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed> >        artPFPartSeedAssociations( new art::Assns<recob::PFParticle, recob::Seed>         );
    std::unique_ptr< art::Assns<recob::Seed,       recob::Hit> >         artSeedHitAssociations(    new art::Assns<recob::Seed,       recob::Hit>         );
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> >         artSPHitAssociations(      new art::Assns<recob::SpacePoint, recob::Hit>          );
    
    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here, except that we override selected items
    // (so, thanks to metaprogramming, we finally have wrappers of wrappers);
    // configuration would happen here, but we are using the default
    // configuration for that algorithm
    using OverriddenClusterParamsAlg_t
      = cluster::OverriddenClusterParamsAlg<cluster::StandardClusterParamsAlg>;
    cluster::ClusterParamsImportWrapper<OverriddenClusterParamsAlg_t>
      ClusterParamAlgo;
    
    // Create id for space points
    int    spacePointID(0);
    int    pcaAxisID(0);
    
    // Make sure there is something to do here!
    if (!clusterParametersList.empty())
    {
        // Grab an instance of the Geometry as it is needed for the association of hits
        art::ServiceHandle<util::LArProperties> theLarProperties;

        // indices for the clusters created
        int    clusterIdx(0);
        size_t pfParticleIdx(0);
        
        // Keep track of the position of the latest added 2D clusters
        size_t clusterStart(0);
        size_t clusterEnd(0);
        
        // Create id for space points
        int    spacePointID(0);
        
        // This is the loop over candidate 3D clusters
        for(auto& clusterParameters : clusterParametersList)
        {
            // It should be straightforward at this point to transfer information from our vector of clusters
            // to the larsoft objects... of course we still have some work to do first, in particular to
            // find the candidate seeds and their seed hits
            
            // We keep track of 2 PCA axes, the first is the "full" PCA run over all the 3D hits in the
            // candidate cluster. The second will be that derived from just using the "skeleton" hits.
            // Make a copy of the full PCA to keep that, then get a reference for the skeleton PCA
            reco::PrincipalComponents& fullPCA     = clusterParameters.m_fullPCA;
            reco::PrincipalComponents& skeletonPCA = clusterParameters.m_skeletonPCA;
            
            // The chances of getting here and this condition not being true are probably zero... but check anyway
            if (!fullPCA.getSvdOK())
            {
                mf::LogDebug("Cluster3D") << "--> no feature extraction done on this cluster!!" << std::endl;
                continue;
            }
            
            // As tracks become more parallel to the wire plane the number of "ambiguous" 3D hits can increase
            // rapidly. Now that we have more information we can go back through these hits and do a better job
            // selecting "the right ones". Here we call the "medial skeleton" algorithm which uses a modification
            // of a standard medial skeleton procedure to get the 3D hits we want
            // But note that even this is hopeless in the worst case and, in fact, it can be a time waster
            // So bypass when you recognize that condition
            if (!(fabs(fullPCA.getEigenVectors()[2][0]) > 0.999 && 3. * sqrt(fullPCA.getEigenValues()[1]) > 20.))
            {
                int nSkeletonPoints = m_skeletonAlg.FindMedialSkeleton(clusterParameters.m_hitPairClusterMapItr->second);
            
                // If enough skeleton points then rerun pca with only those
                if (nSkeletonPoints > 10)
                {
                    // Now rerun the principal components axis on just those points
                    m_pcaAlg.PCAAnalysis_3D(clusterParameters.m_hitPairClusterMapItr->second, skeletonPCA, true);
                
                    // If there was a failure (can that happen?) then restore the full PCA
                    if (!skeletonPCA.getSvdOK()) skeletonPCA = fullPCA;
                }
            }
            
            // Start loop over views to build out the hit lists and the 2D cluster objects
            for(ViewToClusterParamsMap::const_iterator viewItr = clusterParameters.m_clusterParams.begin(); viewItr != clusterParameters.m_clusterParams.end(); viewItr++)
            {
                const RecobClusterParameters& clusParams = viewItr->second;
                
                // We love looping. In this case, our list of hits is comprised of "ClusterHits" and we need to get a RecobHitVector instead...
                RecobHitVector recobHits;
                
                for(HitVectorConst::const_iterator hitItr = clusParams.m_hitVector.begin(); hitItr != clusParams.m_hitVector.end(); hitItr++)
                {
                    art::Ptr<recob::Hit> hitPtr = m_hitToPtrMap[&(*hitItr)->getHit()];
                    recobHits.push_back(hitPtr);
                }
                
                // And sorting! Sorting is good for the mind, soul and body
                // ooopsss... don't do this else event display will look funky
//                std::sort(recobHits.begin(), recobHits.end());
                
                // Get the tdc/wire slope... from the unit vector...
                double dTdW(0.);
                double startWire(clusParams.m_startWire);
                double endWire(clusParams.m_endWire);
                double startTime(clusParams.m_startTime);
                double endTime(clusParams.m_endTime);
                double wirePitch(m_geometry->WirePitch(clusParams.m_view));
                double midWire(0);

                // Get wire number corresponding to the current position
                try
                {
                    midWire = 1.*m_geometry->NearestWire(skeletonPCA.getAvePosition(), clusParams.m_view);
                } catch (cet::exception& e)
                {
                    midWire  = skeletonPCA.getAvePosition()[2];
                    midWire /= wirePitch;
                }
                
                // Now sort out the slope in this view
                // This follows the code in the display pacakage
                double thetaWire     = m_geometry->Plane(clusParams.m_view).Wire(0).ThetaZ();
                double driftVelocity = theLarProperties->DriftVelocity();
                double timeTick      = m_detector->SamplingRate()*1.e-3;
                
                //rotate coord system CCW around x-axis by pi-thetawire
                //   new yprime direction is perpendicular to the wire direction
                //   in the same plane as the wires and in the direction of
                //   increasing wire number
                //use yprime-component of dir cos in rotated coord sys to get
                //   dTdW (number of time ticks per unit of wire pitch)
                double rotAng = 3.1416-thetaWire;
                double yPrime = std::cos(rotAng)*skeletonPCA.getEigenVectors()[0][1]+std::sin(rotAng)*skeletonPCA.getEigenVectors()[0][2];
                
                if (fabs(yPrime) < 0.0000001) yPrime = 0.0000001;
                
                dTdW = skeletonPCA.getEigenVectors()[0][0]*wirePitch/(driftVelocity*timeTick*yPrime);
                
                // Finally, use the time corresponding to this position
                double midTime = m_detector->ConvertXToTicks(skeletonPCA.getAvePosition()[0], clusParams.m_view, 0, 0);
                
                // plane ID is not a part of clusParams... get the one from the first hit
                geo::PlaneID plane; // invalid by default
                if (!recobHits.empty())
                  plane = recobHits.front()->WireID().planeID();
                
                // Ok, now adjust everything to draw a nice long line through our cluster
                double deltaWires = 0.5 * (endWire - startWire);
                
                startWire = midWire - deltaWires;
                endWire   = midWire + deltaWires;
                startTime = midTime - dTdW * deltaWires;
                endTime   = midTime + dTdW * deltaWires;
                
                // feed the algorithm with all the cluster hits
                ClusterParamAlgo.ImportHits(recobHits);
                
                // override the end angles: instead of using the standard
                // algorithm, we use a precomputed value
                ClusterParamAlgo.OverrideParameter
                  (OverriddenClusterParamsAlg_t::cpStartAngle, std::tan(dTdW));
                ClusterParamAlgo.OverrideParameter
                  (OverriddenClusterParamsAlg_t::cpEndAngle, std::tan(dTdW));
                
                // create the recob::Cluster directly in the vector
                cluster::ClusterCreator artCluster(
                  ClusterParamAlgo,                     // algo
                  startWire,                            // start_wire
                  0.,                                   // sigma_start_wire
                  startTime,                            // start_tick
                  clusParams.m_sigmaStartTime,          // sigma_start_tick
                  endWire,                              // end_wire
                  0.,                                   // sigma_end_wire,
                  endTime,                              // end_tick
                  clusParams.m_sigmaEndTime,            // sigma_end_tick
                  clusterIdx++,                         // ID
                  clusParams.m_view,                    // view
                  plane,                                // plane
                  recob::Cluster::Sentry                // sentry
                  );
                
                artClusterVector->emplace_back(artCluster.move());
                             
                util::CreateAssn(*this, evt, *artClusterVector, recobHits, *artClusterAssociations);
                clusterEnd++;
            }
            
            // Deal with converting the Hit Pairs to art
            // Recover the hit pairs and start looping! Love to loop!
            reco::HitPairListPtr& clusHitPairVector = clusterParameters.m_hitPairClusterMapItr->second;
            
            // Last, let's try to get seeds for tracking..
            // Keep track of how many we have so far
            size_t numSeedsStart = artSeedVector->size();
            
            // Call the magical algorith to do the dirty work
            findTrackSeeds(evt, clusterParameters, *artSeedVector, *artSeedHitAssociations);
            
            // Right now error matrix is uniform...
            double spError[] = {1., 0., 1., 0., 0., 1.};
            
            // Keep track of current start for space points
            int spacePointStart(spacePointID);
            
            // Copy these hits to the vector to be stored with the event
            for (auto& hitPair : clusHitPairVector)
            {
                // Don't make space point if this hit was "rejected"
                if (hitPair->getStatusBits() & 0x80000000) continue;
                
                double chisq = 1.;    // secret handshake...
                
                if      ((hitPair->getStatusBits() & 0x30000000) == 0x10000000) chisq = -1.;  // pure skeleton point
                else if ((hitPair->getStatusBits() & 0x30000000) == 0x20000000) chisq = -2.;  // pure edge point
                else if ((hitPair->getStatusBits() & 0x30000000) == 0x30000000) chisq = -3.;  // both skeleton and edge

                if      ((hitPair->getStatusBits() & 0x40000000) == 0x40000000) chisq = -4.;  // Seed point
                
                if (hitPair->getHits().size() < 3) chisq = -10.;
                
                // Mark this hit pair as in use
                hitPair->setStatusBit(0x800000);
                double spacePointPos[] = {hitPair->getPosition()[0],hitPair->getPosition()[1],hitPair->getPosition()[2]};
                artSpacePointVector->push_back(recob::SpacePoint(spacePointPos, spError, chisq, spacePointID++));
            }
            
            // Empty daughter vector for now
            std::vector<size_t> nullVector;
            
            // Create the PFParticle to tie the pieces together
            size_t parentID(recob::PFParticle::kPFParticlePrimary);
            
            recob::PFParticle pfParticle(13, pfParticleIdx++, parentID, nullVector);
            artPFParticleVector->push_back(pfParticle);
            
            // Look at making the PCAxis and associations - for both the skeleton (the first) and the full
            recob::PCAxis skelPcAxis(skeletonPCA.getSvdOK(),
                                     skeletonPCA.getNumHitsUsed(),
                                     skeletonPCA.getEigenValues(),
                                     skeletonPCA.getEigenVectors(),
                                     skeletonPCA.getAvePosition(),
                                     skeletonPCA.getAveHitDoca(),
                                     pcaAxisID++);
            
            artPCAxisVector->push_back(skelPcAxis);
            
            recob::PCAxis fullPcAxis(fullPCA.getSvdOK(),
                                     fullPCA.getNumHitsUsed(),
                                     fullPCA.getEigenValues(),
                                     fullPCA.getEigenVectors(),
                                     fullPCA.getAvePosition(),
                                     fullPCA.getAveHitDoca(),
                                     pcaAxisID++);
            
            artPCAxisVector->push_back(fullPcAxis);
            
            util::CreateAssn(*this, evt, *artPFParticleVector, *artPCAxisVector, *artPFPartAxisAssociations, artPCAxisVector->size()-2, artPCAxisVector->size());
            
            // Create associations to the PFParticle
            util::CreateAssn(*this, evt, *artPFParticleVector, *artSeedVector, *artPFPartSeedAssociations, numSeedsStart, artSeedVector->size());
            
            // Make associations to the 2D cluster objects
            util::CreateAssn(*this, evt, *artPFParticleVector, *artClusterVector, *artPFPartClusAssociations, clusterStart, clusterEnd);
            
            // Make associations to the SpacePoints
            util::CreateAssn(*this, evt, *artPFParticleVector, *artSpacePointVector, *artPFPartSPAssociations, spacePointStart, spacePointID);
            
            // Update the start/end indices
            clusterStart = clusterEnd;
        }
    }
    
    // Right now error matrix is uniform...
    double spError[] = {1., 0., 1., 0., 0., 1.};
    
    // Run through the HitPairVector and add any unused hit pairs to the list
    for(auto& hitPair : hitPairVector)
    {
        unsigned statusBits = hitPair->getStatusBits();
        
        if (!(statusBits & 0x800000))
        {
            double spacePointPos[] = {hitPair->getPosition()[0],hitPair->getPosition()[1],hitPair->getPosition()[2]};
            artSpacePointVector->push_back(recob::SpacePoint(spacePointPos, spError, 1., spacePointID++));
        }
    }
    
    // Finaly done, now output everything to art
    evt.put(std::move(artPCAxisVector));
    evt.put(std::move(artPFParticleVector));
    evt.put(std::move(artClusterVector));
    evt.put(std::move(artSpacePointVector));
    evt.put(std::move(artSeedVector));
    evt.put(std::move(artPFPartAxisAssociations));
    evt.put(std::move(artPFPartClusAssociations));
    evt.put(std::move(artClusterAssociations));
    evt.put(std::move(artPFPartSPAssociations));
    evt.put(std::move(artPFPartSeedAssociations));
    evt.put(std::move(artSeedHitAssociations));
    evt.put(std::move(artSPHitAssociations));
    
    return;
}


} // namespace lar_cluster3d

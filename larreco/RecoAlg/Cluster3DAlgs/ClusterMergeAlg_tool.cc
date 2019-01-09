/**
 *  @file   ClusterMergeAlg.cxx
 * 
 *  @brief  Algorithm for comparing clusters and merging those that are consistent
 * 
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/IClusterModAlg.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// Root includes
#include "TH1F.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {
    
class ClusterMergeAlg : virtual public IClusterModAlg
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit ClusterMergeAlg(const fhicl::ParameterSet&);
    
    /**
     *  @brief  Destructor
     */
    ~ClusterMergeAlg();
    
    void configure(fhicl::ParameterSet const &pset) override;
    
    /**
     *  @brief Interface for initializing histograms if they are desired
     *         Note that the idea is to put hisgtograms in a subfolder
     *
     *  @param TFileDirectory - the folder to store the hists in
     */
    void initializeHistograms(art::TFileDirectory&) override;

    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void ModifyClusters(reco::ClusterParametersList&) const override;
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute() const override {return fTimeToProcess;}
    
private:
    
    bool linearClusters(reco::ClusterParameters&, reco::ClusterParameters&) const;
    
    bool linearClusters2(const reco::PrincipalComponents&, const reco::PrincipalComponents&) const;
    
    bool mergeClusters(reco::ClusterParameters&, reco::ClusterParameters&) const;
    
    float closestApproach(const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, Eigen::Vector3f&, Eigen::Vector3f&, Eigen::Vector3f&) const;
    
    const reco::ClusterHit3D* findClosestHit3D(const Eigen::Vector3f&, const Eigen::Vector3f&, const reco::HitPairListPtr&) const;
    
    /**
     *  @brief Data members to follow
     */
    bool                                 fEnableMonitoring;       ///< If true then turn on monitoring (e.g. timing)
    float                                fAxisAngleScaleFactor;   ///< Addition scaling on the angle selection cut
    float                                fMinTransEigenVal;       ///< Set a mininum allowed value for the transverse eigen values
    float                                fNumTransEigens;         ///< "Num sigma" for cutting on transverse eigen values
    float                                fMaxDOCASeparation;      ///< Minimum DOCA between PCA primary axes to keep from being silly
    float                                fMinEigenToProcess;      ///< Don't look anymore at clusters below this size
    bool                                 fOutputHistograms;       ///< Take the time to create and fill some histograms for diagnostics
    
    
    mutable float                        fTimeToProcess;          ///< Keep track of how long it took to run this algorithm
    
    std::vector<TH1F*>                   fFirstEigenValueHists;   ///< First Cluster eigen value
    std::vector<TH1F*>                   fNextEigenValueHists;    ///< Next Cluster eigen value
    TH1F*                                fNumMergedClusters;      ///< How many clusters were merged?
    TH1F*                                f1stEndPtLenRatHist;     ///< Ratio of distance from first to next PCA to first eigen value
    TH1F*                                f1stPocaProj1Hist;       ///< First endpoint trans 1 doca to next
    TH1F*                                f1stPocaProj2Hist;       ///< First endpoint trans 2 doca to next
    
    TH1F*                                f1st2ndPocaLenHist;      ///< Length of POCA firs axis to 2nd PCA center
    TH1F*                                f2ndPProjEigenHist;      ///< Poca vec projection of next eigen vals
    TH1F*                                f2ndProjEigenRatHist;    ///< Ratio of projected to first eigen values
    TH1F*                                f2ndVProjEigenHist;      ///< First center/ Next center vec projection of next eigen vals
    TH1F*                                f2ndVProjEigenRatHist;   ///< Ratio of projected to first eigen values
    TH1F*                                f1st2ndVecLenHist;       ///< First center/ Next center vector length
    TH1F*                                f1stDProjEigenHist;      ///< First/Next projected first PCA eigen vals
    TH1F*                                f1stProjEigenRatHist;    ///< Ratio of projected to first eigen values
    TH1F*                                f1st2ndPocaRatHist;      ///< ratio of POCA to projected eigen
    TH1F*                                f1st2ndVecRatHist;       ///< ratio of first/next distance to projected eigen
    TH1F*                                f1st2ndVecGapHist;       ///< First center/ Next center vector "gap"

    TH1F*                                fAxesDocaHist;           ///< Doca between two axes
    TH1F*                                f1stDocaArcLRatHist;     ///< Arclen ratio to first POCA
    TH1F*                                f1stDocaArcLRatPHist;    ///< Arclen ratio to first POCA
    TH1F*                                f2ndDocaArcLRatHist;     ///< Arclen ratio to second POCA
    TH1F*                                f2ndDocaArcLRatPHist;    ///< Arclen ratio to second POCA
    
    TH1F*                                fAxisDocaFirstProj1Hist; ///< DOCA projection on 1st axis
    TH1F*                                fAxisDocaFirstProj2Hist; ///< DOCA projection on 2nd axis
    TH1F*                                fAxisDocaNextProj1Hist;  ///< DOCA projection on 1st axis
    TH1F*                                fAxisDocaNextProj2Hist;  ///< DOCA projection on 2nd axis

    TH1F*                                fG1stEndPtLenRatHist;    ///< Ratio of distance from first to next PCA to first eigen value
    TH1F*                                fG1stPocaProj1RatHist;   ///< First endpoint trans 1 doca to next
    TH1F*                                fG1stPocaProj2RatHist;   ///< First endpoint trans 2 doca to next
    TH1F*                                fGAxesDocaHist;          ///< Doca between two axes
    TH1F*                                fG1stDocaArcLRatHist;    ///< Arclen ratio to first POCA
    TH1F*                                fG2ndDocaArcLRatHist;    ///< Arclen ratio to second POCA

    geo::Geometry*                       fGeometry;               //< pointer to the Geometry service
    
    PrincipalComponentsAlg               fPCAAlg;                 // For running Principal Components Analysis
};

ClusterMergeAlg::ClusterMergeAlg(fhicl::ParameterSet const &pset) :
    fPCAAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterMergeAlg::~ClusterMergeAlg()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergeAlg::configure(fhicl::ParameterSet const &pset)
{
    fEnableMonitoring     = pset.get<bool> ("EnableMonitoring",     true  );
    fAxisAngleScaleFactor = pset.get<float>("AxisAngleScaleFactor", 5.    );
    fMinTransEigenVal     = pset.get<float>("MinTransEigenVal",     0.09  );
    fNumTransEigens       = pset.get<float>("NumTransEigenVals",    10.   );
    fMaxDOCASeparation    = pset.get<float>("MaxDOCASeparation",    30.   );
    fMinEigenToProcess    = pset.get<float>("MinEigenToProcess",    2.0   );
    fOutputHistograms     = pset.get<bool> ("OutputHistograms",     false );
    
    art::ServiceHandle<geo::Geometry> geometry;
    
    fGeometry = &*geometry;
    
    fTimeToProcess = 0.;

    // If asked, define some histograms
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;
        
        // Make a directory for these histograms
        art::TFileDirectory dir = tfs->mkdir("MergeClusters");
        
        fFirstEigenValueHists.resize(3,nullptr);
        fNextEigenValueHists.resize(3,nullptr);

        for(size_t idx : {0, 1, 2})
        {
            fFirstEigenValueHists[idx] = dir.make<TH1F>(Form("FEigen1st%1zu",idx),"Eigen Val", 200, 0., 1024./std::max(float(16*idx),float(1)));
            fNextEigenValueHists[idx]  = dir.make<TH1F>(Form("FEigen2nd%1zu",idx),"Eigen Val", 200, 0., 1024./std::max(float(16*idx),float(1)));
        }
        
        fNumMergedClusters      = dir.make<TH1F>("NumMergedClus",      "Number Merged",           200,    0., 200.);
        
        f1stEndPtLenRatHist     = dir.make<TH1F>("FirstNextLenRat",    "First/Next Length ratio", 200,    0.,  20.);
        f1stPocaProj1Hist       = dir.make<TH1F>("FirstProj1Rat",      "First EP Doca 1 ratio",   200,  -50.,  50.);
        f1stPocaProj2Hist       = dir.make<TH1F>("FirstProj2Rat",      "First EP Doca 1 ratio",   200,  -50.,  50.);
        
        f1st2ndPocaLenHist      = dir.make<TH1F>("First2ndPocaLen",    "First/Next POCA len",     200,    0., 100.);
        f2ndPProjEigenHist      = dir.make<TH1F>("ScndProjEigen",      "POCA Projected Eigen",    200,    0., 100.);
        f2ndProjEigenRatHist    = dir.make<TH1F>("ScndProjEigenRat",   "2nd Proj Eigen Rat",      200,    0.,   2.);
        f2ndVProjEigenHist      = dir.make<TH1F>("ScndVProjEigen",     "POCA Projected Eigen",    200,    0., 100.);
        f2ndVProjEigenRatHist   = dir.make<TH1F>("ScndVProjEigenRat",  "2nd Proj Eigen Rat",      200,    0.,   2.);
        f1st2ndVecLenHist       = dir.make<TH1F>("First2ndVecLen",     "First/Next Vec len",      200,    0., 400.);
        f1stDProjEigenHist      = dir.make<TH1F>("FirstProjEigen",     "First/Next Proj Eigen",   200,    0., 400.);
        f1stProjEigenRatHist    = dir.make<TH1F>("FirstProjEigenRat",  "First Proj Eigen Rat",    200,    0.,   2.);
        f1st2ndPocaRatHist      = dir.make<TH1F>("First2ndPocaProj",   "First/Next Poca/Eigen",   200,    0.,  50.);
        f1st2ndVecRatHist       = dir.make<TH1F>("First2ndVecProj",    "First/Next Proj/Eigen",   200,    0.,  50.);
        f1st2ndVecGapHist       = dir.make<TH1F>("First2ndVecGap",     "First/Next Vec len",      200, -400., 400.);
 
        fAxesDocaHist           = dir.make<TH1F>("AxesDoca",           "Axes doca",               200,    0., 200.);
        f1stDocaArcLRatHist     = dir.make<TH1F>("Doca1stArcLRat",     "First arclen ratio",      200,  -25.,  25.);
        f1stDocaArcLRatPHist    = dir.make<TH1F>("Doca1stArcLRatP",    "First arclen P ratio",    200,  -25.,  25.);
        f2ndDocaArcLRatHist     = dir.make<TH1F>("Doca2ndArcLRat",     "First arclen ratio",      200,  -25.,  25.);
        f2ndDocaArcLRatPHist    = dir.make<TH1F>("Doca2ndArcLRatP",    "First arclen P ratio",    200,  -25.,  25.);
        
        fAxisDocaFirstProj1Hist = dir.make<TH1F>("AxisDocaFProj1",     "Axis Doca First Proj 1",  200,    0., 100.);
        fAxisDocaFirstProj2Hist = dir.make<TH1F>("AxisDocaFProj2",     "Axis Doca First Proj 1",  200,    0., 100.);
        fAxisDocaNextProj1Hist  = dir.make<TH1F>("AxisDocaNProj1",     "Axis Doca First Proj 1",  200,    0., 100.);
        fAxisDocaNextProj2Hist  = dir.make<TH1F>("AxisDocaNProj2",     "Axis Doca First Proj 1",  200,    0., 100.);

        fG1stEndPtLenRatHist    = dir.make<TH1F>("GFirstNextLenRat",   "First/Next Length ratio", 200,    0., 20.);
        fG1stPocaProj1RatHist   = dir.make<TH1F>("GFirstEPDoca1Rat",   "First EP Doca 1 ratio",   200,   -1.,  1.);
        fG1stPocaProj2RatHist   = dir.make<TH1F>("GFirstEPDoca2Rat",   "First EP Doca 1 ratio",   200,   -1.,  1.);
 
        fGAxesDocaHist          = dir.make<TH1F>("GAxesDoca",          "Axes doca",               250,    0., 500.);
        fG1stDocaArcLRatHist    = dir.make<TH1F>("GDoca1stArcLRat",    "First arclen ratio",      200,  -25., 25.);
        fG2ndDocaArcLRatHist    = dir.make<TH1F>("GDoca2ndArcLRat",    "First arclen ratio",      200,  -25., 25.);
    }
    
    return;
}
    
void ClusterMergeAlg::initializeHistograms(art::TFileDirectory&)
{
    return;
}

void ClusterMergeAlg::ModifyClusters(reco::ClusterParametersList& clusterParametersList) const
{
    /**
     *  @brief Top level interface for algorithm to consider pairs of clusters from the input
     *         list and determine if they are consistent with each other and, therefore, should
     *         be merged. This is done by looking at the PCA for each cluster and looking at the
     *         projection of the primary axis along the vector connecting their centers.
     */
    
    // Initial clustering is done, now trim the list and get output parameters
    cet::cpu_timer theClockBuildClusters;
    
    // Start clocks if requested
    if (fEnableMonitoring) theClockBuildClusters.start();
    
    // Resort by the largest first eigen value (so the "longest" cluster)
    clusterParametersList.sort([](auto& left, auto& right){return left.getFullPCA().getEigenValues()[0] > right.getFullPCA().getEigenValues()[0];});
    
    // The idea is to continually loop through all clusters until we get to the point where we are no longer doing any merging.
    // If two clusters are merged then we need to recycle through again because it may be that the new cluster can match when
    // it might not have in the past
    size_t                                lastClusterListCount = clusterParametersList.size() + 1;
    reco::ClusterParametersList::iterator lastFirstClusterItr  = clusterParametersList.begin();
    
    int numMergedClusters = 0;
    
    int outerLoopCount(0);

    while(clusterParametersList.size() != lastClusterListCount)
    {
        // Update the last count
        lastClusterListCount = clusterParametersList.size();
        
        // Keep track of the first cluster iterator each pass through
        reco::ClusterParametersList::iterator firstClusterItr = lastFirstClusterItr++;
        
        std::cout << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "        Outer loop # " << outerLoopCount++ << std::endl;
        
        int clusterIdx(0);

        // Loop through the clusters
        while(firstClusterItr != clusterParametersList.end())
        {
            reco::ClusterParameters&              firstClusterParams = *firstClusterItr;
            reco::ClusterParametersList::iterator nextClusterItr     = firstClusterItr;
            
            // Take a brief interlude to do some histogramming
            if (fOutputHistograms)
            {
                const reco::PrincipalComponents& firstPCA = firstClusterParams.getFullPCA();
                
                Eigen::Vector3f firstEigenVals(         2.0 * sqrt(firstPCA.getEigenValues()[0]),
                                               std::min(2.0 * sqrt(std::max(firstPCA.getEigenValues()[1],fMinTransEigenVal)), 100.),
                                               std::min(2.0 * sqrt(std::max(firstPCA.getEigenValues()[2],fMinTransEigenVal)), 50.));
                
                for(size_t idx = 0; idx < 3; idx++) fFirstEigenValueHists[idx]->Fill(firstEigenVals[idx],1.);
            }

            // Once you get down to the smallest clusters if they haven't already been absorbed there is no need to check them
            if (firstClusterParams.getFullPCA().getEigenValues()[0] < fMinEigenToProcess) break;
            
            std::cout << "************** First Cluster ptr: " << &firstClusterParams << ", internal count: " << clusterIdx++ << " **************" << std::endl;
            std::cout << "   PCA eigen: " << 1.5*std::sqrt(firstClusterParams.getFullPCA().getEigenValues()[0]) << ", " << 1.5*std::sqrt(firstClusterParams.getFullPCA().getEigenValues()[1]) << ", " << 1.5*std::sqrt(firstClusterParams.getFullPCA().getEigenValues()[2])  << std::endl;
            
            // want the next one...
            nextClusterItr++;
            
            while(nextClusterItr != clusterParametersList.end())
            {
                reco::ClusterParameters& nextClusterParams = *nextClusterItr;
                
                // On any given loop through here it **might** be that the first cluster has been modified. So can't cache
                // the parameters, need the curret ones
                if (linearClusters(firstClusterParams,nextClusterParams))
                {
                    if (mergeClusters(firstClusterParams, nextClusterParams))
                    {
                        Eigen::Vector3f firstPosToNextPosVec  = nextClusterParams.getFullPCA().getAvePosition() - firstClusterParams.getFullPCA().getAvePosition();
                        std::cout << "   ---> Merged cluster ptr: " << &nextClusterParams << ", dist: " << firstPosToNextPosVec.norm() << ", eigen: " << 1.5*std::sqrt(nextClusterParams.getFullPCA().getEigenValues()[0]) << ", " << 1.5*std::sqrt(nextClusterParams.getFullPCA().getEigenValues()[1]) << ", " <<  1.5*std::sqrt(nextClusterParams.getFullPCA().getEigenValues()[2]) << std::endl;
                        std::cout << "   ~~~> New PCA eigen: " << 1.5*std::sqrt(firstClusterParams.getFullPCA().getEigenValues()[0]) << ", " << 1.5*std::sqrt(firstClusterParams.getFullPCA().getEigenValues()[1]) << ", " << 1.5*std::sqrt(firstClusterParams.getFullPCA().getEigenValues()[2])  << std::endl;

                        // Now remove the "next" cluster
                        nextClusterItr = clusterParametersList.erase(nextClusterItr);
                        
                        // Restart loop?
                        nextClusterItr = firstClusterItr;
                        
                        numMergedClusters++;
                    }
                }
                
                nextClusterItr++;
            }
            
            firstClusterItr++;
        }
    }
    
    std::cout << "***>> Merged " << numMergedClusters << std::endl;
    
    if (fOutputHistograms) fNumMergedClusters->Fill(numMergedClusters, 1.);
    
    
    if (fEnableMonitoring)
    {
        theClockBuildClusters.stop();
        
        fTimeToProcess = theClockBuildClusters.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> Merge clusters done, found " << clusterParametersList.size() << " clusters" << std::endl;
    
    return;
}
    
bool ClusterMergeAlg::linearClusters(reco::ClusterParameters& firstCluster, reco::ClusterParameters& nextCluster) const
{
    // Assume failure
    bool consistent(false);
    
    // The goal here is to compare the two input PCA's and determine if they are effectively colinear and
    // within reasonable range to consider merging them. Note that a key assumption is that the first input
    // PCA is from the "bigger" cluster, the second is "smaller" and may have a less reliable PCA.
    
    // Dereference the PCA's
    const reco::PrincipalComponents& firstPCA = firstCluster.getFullPCA();
    const reco::PrincipalComponents& nextPCA  = nextCluster.getFullPCA();

    // Recover the positions of the centers of the two clusters
    const Eigen::Vector3f& firstCenter = firstPCA.getAvePosition();
    const Eigen::Vector3f& nextCenter  = nextPCA.getAvePosition();
    
    // And form a vector between the two centers
    Eigen::Vector3f firstPosToNextPosVec  = nextCenter - firstCenter;
    Eigen::Vector3f firstPosToNextPosUnit = firstPosToNextPosVec.normalized();
    
    // Now get the first PCA's primary axis and since we'll use them get all of them at once...
    Eigen::Vector3f firstAxis0(firstPCA.getEigenVectors().row(0));
    Eigen::Vector3f firstAxis1(firstPCA.getEigenVectors().row(1));
    Eigen::Vector3f firstAxis2(firstPCA.getEigenVectors().row(2));
    
    // Will want the distance of closest approach of the next cluser's center to the primary, start by finding arc length
    float arcLenToNextDoca = firstPosToNextPosVec.dot(firstAxis0);
    
    // Adopt the convention that the cluster axis is in same direction as vector from first to next centers
    // And preserve the overall orientation of the PCA by flipping all if we flip the first one
    if (arcLenToNextDoca < 0.)
    {
        firstAxis0       = -firstAxis0;
        firstAxis1       = -firstAxis1;
        firstAxis2       = -firstAxis2;
        arcLenToNextDoca = -arcLenToNextDoca;
    }
    
    // Recover the eigen values of the first and second PCAs for selection cuts
    Eigen::Vector3f firstEigenVals(         1.5 * sqrt(         firstPCA.getEigenValues()[0]                          ),
                                   std::min(1.5 * sqrt(std::max(firstPCA.getEigenValues()[1],fMinTransEigenVal)), 100.),
                                   std::min(1.5 * sqrt(std::max(firstPCA.getEigenValues()[2],fMinTransEigenVal)),  50.));
    
    // Develop metric for max allowed angle from the eigen values
    float rMaxFirst    = std::sqrt(firstEigenVals[1] * firstEigenVals[1] + firstEigenVals[2] * firstEigenVals[2]);
    float cosMaxFirst  = std::max(firstEigenVals[0] / std::sqrt(firstEigenVals[0] * firstEigenVals[0] + fNumTransEigens * rMaxFirst * fNumTransEigens * rMaxFirst),float(0.8));
    float cosFirstAxis = firstAxis0.dot(firstPosToNextPosUnit);
    
    // Get a measure of the projection of the three eigen values along the axis from the first to next PCA centers
    float firstToNextProjEigen = firstEigenVals[0] * std::abs(firstPosToNextPosUnit.dot(firstAxis0))
                               + firstEigenVals[1] * std::abs(firstPosToNextPosUnit.dot(firstAxis1))
                               + firstEigenVals[2] * std::abs(firstPosToNextPosUnit.dot(firstAxis2));

    // This makes first selection:
    if (cosFirstAxis > cosMaxFirst && arcLenToNextDoca < 3. * firstToNextProjEigen)
    {
        // Recover the axes for the next PCA and make sure pointing convention is observed
        Eigen::Vector3f nextAxis0(nextPCA.getEigenVectors().row(0));
        Eigen::Vector3f nextAxis1(nextPCA.getEigenVectors().row(1));
        Eigen::Vector3f nextAxis2(nextPCA.getEigenVectors().row(2));
        
        // Get the arc length from next to first centers - note convention for PCA orientation
        float arcLenToFirstDoca = -firstPosToNextPosVec.dot(nextAxis0);
        
        // Convention on axis direction applied again
        if (arcLenToFirstDoca > 0.)
        {
            nextAxis0         = -nextAxis0;
            nextAxis1         = -nextAxis1;
            nextAxis2         = -nextAxis2;
            arcLenToFirstDoca = -arcLenToFirstDoca;
        }
        
        // Get the eigen values again
        Eigen::Vector3f nextEigenVals(         1.5 * sqrt(nextPCA.getEigenValues()[0]),
                                      std::min(1.5 * sqrt(std::max(nextPCA.getEigenValues()[1],fMinTransEigenVal)), 64.),
                                      std::min(1.5 * sqrt(std::max(nextPCA.getEigenValues()[2],fMinTransEigenVal)), 32.));
        
        // Develop metric for max allowed angle from the eigen values
        float rMaxNext    = std::sqrt(nextEigenVals[1] * nextEigenVals[1] + nextEigenVals[2] * nextEigenVals[2]);
        float cosMaxNext  = std::max(nextEigenVals[0] / std::sqrt(nextEigenVals[0] * nextEigenVals[0] + fNumTransEigens * rMaxNext * fNumTransEigens * rMaxNext),float(0.8));
        float cosNextAxis = nextAxis0.dot(firstPosToNextPosUnit);
        
        // To really get the gap we need to look at the projected distance between the last hit in the first cluster and the first hit
        // in the next cluster...
        const reco::ClusterHit3D* closestNextHit3D      = findClosestHit3D(firstCenter, firstPosToNextPosUnit, nextCluster.getHitPairListPtr());
        Eigen::Vector3f           nextPosToFirstPosUnit = -firstPosToNextPosUnit;
        const reco::ClusterHit3D* closestFirstHit3D     = findClosestHit3D(nextCenter, nextPosToFirstPosUnit, firstCluster.getHitPairListPtr());
        Eigen::Vector3f           firstToNextGapVec     = closestNextHit3D->getPosition() - closestFirstHit3D->getPosition();

        // Guesstimate the gap between the two PCAs
        float firstToNextGap = firstToNextGapVec.dot(firstPosToNextPosUnit);
        
        // Determine the projected "eigen distance"
        // First for the next cluster...
        float nextToFirstProjEigen = nextEigenVals[0] * std::abs(nextPosToFirstPosUnit.dot(nextAxis0))
                                   + nextEigenVals[1] * std::abs(nextPosToFirstPosUnit.dot(nextAxis1))
                                   + nextEigenVals[2] * std::abs(nextPosToFirstPosUnit.dot(nextAxis2));
        
        std::cout << "  --> dCenters: " << firstPosToNextPosVec.norm() << ", gap: " << firstToNextGap << ", next eig: " << nextEigenVals[0] << ", proj: " << nextToFirstProjEigen << ", cos: " << cosNextAxis << ", first eig: " << firstEigenVals[0] << ", proj: " << firstToNextProjEigen << ", cos: " << cosFirstAxis << std::endl;

        // This makes first selection:
        if (cosNextAxis > cosMaxNext && firstToNextGap < nextToFirstProjEigen)
        {
            // We need to watch out for the case of colinearity which can screw up a calculation of the distance of closest approach
            // Basically, we look to see if the angle between the two axes is less than 5 degrees
            if (firstAxis0.dot(nextAxis0) > 0.9961)
            {
                // Check the transverse distance
                Eigen::Vector3f pocaVec = nextCenter - (firstCenter + arcLenToNextDoca * firstAxis0);
                
//                std::cout << "**>> cos: " << firstAxis0.dot(nextAxis0) << ", dist: " << firstPosToNextPosVec.norm() << ", 1st eigen: " << firstEigenVals[0] << ", 2nd: " << nextEigenVals[0] << ", poca len: " << pocaVec.norm() << ", max: " << rMaxFirst << std::endl;
                
                if (pocaVec.norm() < rMaxFirst) consistent = true;
            }
            
            // Might be nearly colinear but outside of the above off
            if (!consistent)
            {
                // Closest approach calculaiton results vectors
                Eigen::Vector3f firstPoca;
                Eigen::Vector3f nextPoca;
                Eigen::Vector3f firstToNextVec;
                
                // Recover the doca of the two axes and their points of closest approach along each axis
                float lineDoca = closestApproach(firstCenter, firstAxis0, nextCenter, nextAxis0, firstPoca, nextPoca, firstToNextVec);
                
                // The below returned the signed arc lengths to their respective pocas
                float arcLenToFirstPoca = (firstPoca - firstCenter).dot(firstAxis0);
                float arcLenToNextPoca  = (nextPoca  - nextCenter ).dot(nextAxis0);
                
                if (fOutputHistograms)
                {
                    fAxesDocaHist->Fill(lineDoca, 1.);
                    f1stDocaArcLRatHist->Fill(arcLenToFirstPoca/firstEigenVals[0],1.);
                    f2ndDocaArcLRatHist->Fill(arcLenToNextPoca/nextEigenVals[0],1.);
                }
                
                if (lineDoca < fMaxDOCASeparation &&
                    arcLenToFirstPoca > 0. && arcLenToFirstPoca          < 2. * firstEigenVals[0] &&
                    arcLenToNextPoca  < 0. && std::abs(arcLenToNextPoca) < 5. * nextEigenVals[0]    )
                {
                    consistent = true;
                }
            }
        }
    }

    return consistent;
}

bool ClusterMergeAlg::linearClusters2(const reco::PrincipalComponents& firstPCA, const reco::PrincipalComponents& nextPCA) const
{
    // Assume failure
    bool consistent(false);
    
    // The goal here is to compare the two input PCA's and determine if they are effectively colinear and
    // within reasonable range to consider merging them. Note that a key assumption is that the first input
    // PCA is from the "bigger" cluster, the second is "smaller" and may have a less reliable PCA.
    
    // Recover the positions of the centers of the two clusters
    const Eigen::Vector3f& firstCenter = firstPCA.getAvePosition();
    const Eigen::Vector3f& nextCenter  = nextPCA.getAvePosition();
    
    // And form a vector between the two centers
    Eigen::Vector3f firstPosToNextPosVec  = nextCenter - firstCenter;
    Eigen::Vector3f firstPosToNextPosUnit = firstPosToNextPosVec.normalized();
    
    // Now get the first PCA's primary axis and since we'll use them get all of them at once...
    Eigen::Vector3f firstAxis0(firstPCA.getEigenVectors().row(0));
    Eigen::Vector3f firstAxis1(firstPCA.getEigenVectors().row(1));
    Eigen::Vector3f firstAxis2(firstPCA.getEigenVectors().row(2));
    
    // Will want the distance of closest approach of the next cluser's center to the primary, start by finding arc length
    float arcLenToNextDoca = firstPosToNextPosVec.dot(firstAxis0);
    
    // Adopt the convention that the cluster axis is in same direction as vector from first to next centers
    // And preserve the overall orientation of the PCA by flipping all if we flip the first one
    if (arcLenToNextDoca < 0.)
    {
        firstAxis0       = -firstAxis0;
        firstAxis1       = -firstAxis1;
        firstAxis2       = -firstAxis2;
        arcLenToNextDoca = -arcLenToNextDoca;
    }
    
    // Recover the eigen values of the first and second PCAs for selection cuts
    Eigen::Vector3f firstEigenVals(         1.5 * sqrt(firstPCA.getEigenValues()[0]),
                                   std::min(1.5 * sqrt(std::max(firstPCA.getEigenValues()[1],fMinTransEigenVal)), 64.),
                                   std::min(1.5 * sqrt(std::max(firstPCA.getEigenValues()[2],fMinTransEigenVal)), 32.));
    
    // First task is to develop something that determines that the center of the next PCA is essentially colinear with the first one.
    // In an ideal world we would imagine the tolerance for this being automagically generated... for example by using the transverse
    // eigen values...
    // First up is to determine the distance the next PCAs center is from the first axis
    Eigen::Vector3f firstPocaToNextCenter = nextCenter - (firstCenter + arcLenToNextDoca * firstAxis0);
    
    // This should be perpendicular to the first axis, check:
    float cosPocaToFirst = firstAxis0.dot(firstPocaToNextCenter.normalized());
    
    if (std::abs(cosPocaToFirst) > 2.e-4)
    {
        std::cout << "POCA vector is not perpendicular to the first axis: " << cosPocaToFirst << std::endl;
    }
    
    // Get the projections of this vector along both transverse axes
    float firstPocaProj1 = firstPocaToNextCenter.dot(firstAxis1);
    float firstPocaProj2 = firstPocaToNextCenter.dot(firstAxis2);
    
    // We are going to then compare these with the eigen values, which will define an acceptance cone along the direction
    // of the first axis. If the projectsion fall in range then we are good to go to the next step
    // If the next PCA is some distance from the first then we want to scale the the eigenvalues accordingly
    float pocaScaleFactor = std::max(arcLenToNextDoca / firstEigenVals[0],float(1.));
    
    // Do some histogramming if asked
    if (fOutputHistograms)
    {
        f1stEndPtLenRatHist->Fill(std::min(arcLenToNextDoca/firstEigenVals[0],float(19.99)),1.);
        f1stPocaProj1Hist->Fill(firstPocaProj1/(pocaScaleFactor * firstEigenVals[1]),1.);
        f1stPocaProj2Hist->Fill(firstPocaProj2/(pocaScaleFactor * firstEigenVals[2]),1.);
    }
    
    // Ok, check colinearity
    if (std::abs(firstPocaProj1) < fNumTransEigens * pocaScaleFactor * firstEigenVals[1] &&
        std::abs(firstPocaProj2) < fNumTransEigens * pocaScaleFactor * firstEigenVals[2]   )
    {
        // Recover the axes for the next PCA and make sure pointing convention is observed
        Eigen::Vector3f nextAxis0(nextPCA.getEigenVectors().row(0));
        Eigen::Vector3f nextAxis1(nextPCA.getEigenVectors().row(1));
        Eigen::Vector3f nextAxis2(nextPCA.getEigenVectors().row(2));
        
        // Get the arc length from next to first centers - note convention for PCA orientation
        float arcLenToFirstDoca = -firstPosToNextPosVec.dot(nextAxis0);
        
        // Convention on axis direction applied again
        if (arcLenToFirstDoca > 0.)
        {
            nextAxis0         = -nextAxis0;
            nextAxis1         = -nextAxis1;
            nextAxis2         = -nextAxis2;
            arcLenToFirstDoca = -arcLenToFirstDoca;
        }
        
        // Handle the special case that the cluster we are comparing to is "embedded" in the current cluster
        if (arcLenToNextDoca < firstEigenVals[0] && std::abs(firstPocaProj1) < firstEigenVals[1] && std::abs(firstPocaProj2) < 2. * firstEigenVals[2])
        {
            // Either the primary or the secondary PCA axis will be more or less aligned with the first primary
            // (the secondary alignment is a special case for isochronous clusters)
            if (firstAxis0.dot(nextAxis0) > 0.95 || std::abs(firstAxis1.dot(nextAxis1)) > 0.95) consistent = true;
        }
        
        // If no consistency don't give up, more detailed checks may be necessary
        if (!consistent)
        {
            // First check if these two clusters can possibly be consistent where we'll look at the projected distance
            // vs the projected distance from the PCA center to its "shell" defined by the eigen values
            
            // Get the eigen values again
            Eigen::Vector3f nextEigenVals(         1.5 * sqrt(nextPCA.getEigenValues()[0]),
                                          std::min(1.5 * sqrt(std::max(nextPCA.getEigenValues()[1],fMinTransEigenVal)), 64.),
                                          std::min(1.5 * sqrt(std::max(nextPCA.getEigenValues()[2],fMinTransEigenVal)), 32.));
            
            // Convert the poca vector to a unit vector but remember its magnitude
            float           firstPocaLen              = firstPocaToNextCenter.norm();
            Eigen::Vector3f firstPocaToNextCenterUnit = firstPocaToNextCenter.normalized();

            // Armed with this info we can get projected eigen value distances...
            float nextPocaProjEigen = nextEigenVals[0] * std::abs(firstPocaToNextCenterUnit.dot(nextAxis0))
                                    + nextEigenVals[1] * std::abs(firstPocaToNextCenterUnit.dot(nextAxis1))
                                    + nextEigenVals[2] * std::abs(firstPocaToNextCenterUnit.dot(nextAxis2));
            
            float nextToFirstProjEigen = nextEigenVals[0] * std::abs(firstPosToNextPosUnit.dot(nextAxis0))
                                       + nextEigenVals[1] * std::abs(firstPosToNextPosUnit.dot(nextAxis1))
                                       + nextEigenVals[2] * std::abs(firstPosToNextPosUnit.dot(nextAxis2));
            
            // Also look at the same for the first PCA to the next
            float firstToNextProjEigen = firstEigenVals[0] * std::abs(firstPosToNextPosUnit.dot(firstAxis0))
                                       + firstEigenVals[1] * std::abs(firstPosToNextPosUnit.dot(firstAxis1))
                                       + firstEigenVals[2] * std::abs(firstPosToNextPosUnit.dot(firstAxis2));
            
            float firstPosToNextPosLen = firstPosToNextPosVec.norm();
            float firstPosToNextPosGap = firstPosToNextPosLen - firstToNextProjEigen - nextToFirstProjEigen;
            
            if (fOutputHistograms)
            {
                for(size_t idx = 0; idx < 3; idx++)
                {
                    fNextEigenValueHists[idx]->Fill(firstEigenVals[idx],1.);
                }
                
                f1st2ndPocaLenHist->Fill(firstPocaLen, 1.);
                f2ndPProjEigenHist->Fill(nextPocaProjEigen, 1.);
                f2ndVProjEigenHist->Fill(nextToFirstProjEigen, 1.);
                f2ndProjEigenRatHist->Fill(nextPocaProjEigen/nextEigenVals[0], 1.);
                f2ndVProjEigenRatHist->Fill(nextToFirstProjEigen/nextEigenVals[0], 1.);
                f1st2ndVecLenHist->Fill(firstPosToNextPosLen, 1.);
                f1st2ndVecGapHist->Fill(firstPosToNextPosGap, 1.);
                f1stDProjEigenHist->Fill(firstToNextProjEigen, 1.);
                f1stProjEigenRatHist->Fill(firstToNextProjEigen/firstEigenVals[0], 1.);
                f1st2ndPocaRatHist->Fill(firstPocaLen/nextPocaProjEigen, 1.);
                f1st2ndVecRatHist->Fill(firstPosToNextPosLen/firstToNextProjEigen, 1.);
            }
            
            // Require both distance between clusters and projected distance to first axis to be "in range"
//            if (firstPosToNextPosLen < 4. * firstToNextProjEigen && firstPocaLen < 20. * nextPocaProjEigen)
            if (firstToNextProjEigen/firstEigenVals[0] > 0.9 && nextToFirstProjEigen/nextEigenVals[0] > 0.95)
            {
                // Strategy now will be to find the points of closest approach and compare these to the with of the PCAs
                // If this good then look at the gap along the vector from the first to second PCA center.
                // Require the gap to be less than the smaller PCA primary eigenvalue
                // Now we recover the doca of the two primary axes of the input clusters
                
                // Closest approach calculaiton results vectors
                Eigen::Vector3f firstPoca;
                Eigen::Vector3f nextPoca;
                Eigen::Vector3f firstToNextVec;
                
                // Recover the doca of the two axes and their points of closest approach along each axis
                float lineDoca = closestApproach(firstCenter, firstAxis0, nextCenter, nextAxis0, firstPoca, nextPoca, firstToNextVec);
                
                // The below returned the signed arc lengths to their respective pocas
                float arcLenToFirstPoca = (firstPoca - firstCenter).dot(firstAxis0);
                float arcLenToNextPoca  = (nextPoca  - nextCenter ).dot(nextAxis0);
                
                if (fOutputHistograms)
                {
                    fAxesDocaHist->Fill(lineDoca, 1.);
                    f1stDocaArcLRatHist->Fill(arcLenToFirstPoca/firstEigenVals[0],1.);
                    f1stDocaArcLRatPHist->Fill(arcLenToFirstPoca/firstToNextProjEigen,1.);
                    f2ndDocaArcLRatHist->Fill(arcLenToNextPoca/nextEigenVals[0],1.);
                    f2ndDocaArcLRatPHist->Fill(arcLenToNextPoca/nextPocaProjEigen,1.);
                }
                
                if (lineDoca < fMaxDOCASeparation &&
                    arcLenToFirstPoca > 0. && arcLenToFirstPoca          <  5. * firstEigenVals[0] &&
                    arcLenToNextPoca  < 0. && std::abs(arcLenToNextPoca) < 10. * nextEigenVals[0])
                {
                    // Compute docas to both PCA axes
                    float firstToNextPocaVecProj1 = std::fabs(firstAxis1.dot(firstToNextVec));
                    float firstToNextPocaVecProj2 = std::fabs(firstAxis2.dot(firstToNextVec));
                    float nextToFirstPocaVecProj1 = std::fabs(nextAxis1.dot(firstToNextVec));
                    float nextToFirstPocaVecProj2 = std::fabs(nextAxis2.dot(firstToNextVec));
                    
                    // Get angular expansion factors based on the points of closest approach of the two axes
                    float firstAngleScaleFactor = std::min(fAxisAngleScaleFactor, std::max(std::abs(arcLenToFirstPoca) / firstEigenVals(0),float(1.)));
                    float nextAngleScaleFactor  = std::min(fAxisAngleScaleFactor, std::max(std::abs(arcLenToNextPoca)  / nextEigenVals(0), float(1.)));
                    
                    if (fOutputHistograms)
                    {
                        fAxisDocaFirstProj1Hist->Fill(firstToNextPocaVecProj1, 1.);
                        fAxisDocaFirstProj2Hist->Fill(firstToNextPocaVecProj2, 1.);
                        fAxisDocaNextProj1Hist->Fill(nextToFirstPocaVecProj1, 1.);
                        fAxisDocaNextProj2Hist->Fill(nextToFirstPocaVecProj2, 1.);
                    }
                    
                    // Check these are consistent with the eigen values
                    if (firstToNextPocaVecProj1 < firstAngleScaleFactor*firstEigenVals[1] && firstToNextPocaVecProj2 < firstAngleScaleFactor*firstEigenVals[2] &&
                        nextToFirstPocaVecProj1 < nextAngleScaleFactor*nextEigenVals[1]   && nextToFirstPocaVecProj2 < nextAngleScaleFactor*nextEigenVals[2])
                    {
                        // So far so good, now look at the gap...
                        // Last up try to compute the gap between the two clusters. This will be taken as
                        // the projections of the eigenvectors of length given by their eigenvalues along
                        // the vector between the two centers
                        float projectionFirstClus  = firstEigenVals(0) * firstPosToNextPosUnit.dot(firstAxis0);
                        float projectionNextClus   = nextEigenVals(0)  * firstPosToNextPosUnit.dot(nextAxis0);
                        
                        float gap = firstPosToNextPosLen - projectionFirstClus - projectionNextClus;
                        
                        if (gap < 20000. * nextEigenVals(0)) consistent = true;
                        
                        if (fOutputHistograms)
                        {
                            fG1stEndPtLenRatHist->Fill(std::min(arcLenToNextDoca/firstEigenVals[0],float(19.99)),1.);
                            fG1stPocaProj1RatHist->Fill(std::abs(firstPocaProj1)/(pocaScaleFactor * firstEigenVals[1]),1.);
                            fG1stPocaProj2RatHist->Fill(std::abs(firstPocaProj2)/(pocaScaleFactor * firstEigenVals[2]),1.);
                            
                            fGAxesDocaHist->Fill(lineDoca, 1.);
                            fG1stDocaArcLRatHist->Fill(arcLenToFirstPoca/firstEigenVals[0],1.);
                            fG2ndDocaArcLRatHist->Fill(arcLenToNextPoca/nextEigenVals[0],1.);
                        }
                    }
                }
            }
        }
    }

    
    return consistent;
}

bool ClusterMergeAlg::mergeClusters(reco::ClusterParameters& firstClusterParams, reco::ClusterParameters& nextClusterParams) const
{
    bool merged(false);
    
    // Merge the next cluster into the first one
    // Do this by making a local copy of the input cluster parameters
    reco::ClusterParameters clusterParams = firstClusterParams;
    
    // Get the hits
    reco::HitPairListPtr& hitPairListPtr = clusterParams.getHitPairListPtr();
    
    for(const auto* hit : nextClusterParams.getHitPairListPtr())
    {
        hitPairListPtr.push_back(hit);
        
        for(const auto* hit2D : hit->getHits())
            if (hit2D) clusterParams.UpdateParameters(hit2D);
    }
    
    // Recalculate the PCA
    fPCAAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());

//    if (newEigen[0] > 0.8 * origEigen[0] && newEigen[1] * origEigen[0] < 5. * newEigen[0] * origEigen[1])
    {
        // Must have a valid pca
        if (clusterParams.getFullPCA().getSvdOK())
        {
            // Finish out the merging here
            reco::Hit3DToEdgeMap& firstEdgeMap = clusterParams.getHit3DToEdgeMap();
            reco::Hit3DToEdgeMap& nextEdgeMap  = nextClusterParams.getHit3DToEdgeMap();
            
            for(const auto& pair : nextEdgeMap) firstEdgeMap[pair.first] = pair.second;
            
            // Set the skeleton PCA to make sure it has some value
            clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();
            
            // Copy back to the input
            firstClusterParams = clusterParams;
            
            merged = true;
        }
    }
    
    return merged;
}
    
float ClusterMergeAlg::closestApproach(const Eigen::Vector3f& P0, const Eigen::Vector3f& u0,
                                       const Eigen::Vector3f& P1, const Eigen::Vector3f& u1,
                                       Eigen::Vector3f&       poca0,
                                       Eigen::Vector3f&       poca1,
                                       Eigen::Vector3f&       firstNextVec) const
{
    // Technique is to compute the arclength to each point of closest approach
    Eigen::Vector3f w0 = P0 - P1;
    float a(1.);
    float b(u0.dot(u1));
    float c(1.);
    float d(u0.dot(w0));
    float e(u1.dot(w0));
    float den(a * c - b * b);
    
    if (std::abs(den) > std::numeric_limits<float>::epsilon())
    {
        float arcLen0 = (b * e - c * d) / den;
        float arcLen1 = (a * e - b * d) / den;
        
        poca0 = P0 + arcLen0 * u0;
        poca1 = P1 + arcLen1 * u1;
    }
    else
    {
        poca0 = P0;
        poca1 = P1;
    }
    
    firstNextVec = poca1 - poca0;
    
    return firstNextVec.norm();
}
    
const reco::ClusterHit3D* ClusterMergeAlg::findClosestHit3D(const Eigen::Vector3f& refPoint, const Eigen::Vector3f& refVector, const reco::HitPairListPtr& hitList) const
{
    const reco::ClusterHit3D* nearestHit3D(hitList.front());
    float                     closest(std::numeric_limits<float>::max());
    
    for(const auto& hit3D : hitList)
    {
        Eigen::Vector3f refToHitVec = hit3D->getPosition() - refPoint;
        float           arcLenToHit = refToHitVec.dot(refVector);
        
        if (arcLenToHit < closest)
        {
            nearestHit3D = hit3D;
            closest      = arcLenToHit;
        }
    }
    
    return nearestHit3D;
}

    
DEFINE_ART_CLASS_TOOL(ClusterMergeAlg)
} // namespace lar_cluster3d

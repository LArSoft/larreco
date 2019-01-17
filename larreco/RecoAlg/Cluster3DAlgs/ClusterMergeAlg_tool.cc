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

    bool mergeClusters(reco::ClusterParameters&, reco::ClusterParameters&) const;
    
    float closestApproach(const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, Eigen::Vector3f&, Eigen::Vector3f&, Eigen::Vector3f&) const;
    
    const reco::ClusterHit3D* findClosestHit3D(const Eigen::Vector3f&, const Eigen::Vector3f&, const reco::HitPairListPtr&) const;
    
    const reco::ClusterHit3D* findFurthestHit3D(const Eigen::Vector3f&, const Eigen::Vector3f&, const reco::HitPairListPtr&) const;

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
    
    TH1F*                                f1stTo2ndPosLenHist;     ///< Distance between cluster centers
    
    TH1F*                                fRMaxFirstHist;          ///< radius of "cylinder" first cluster
    TH1F*                                fCosMaxFirstHist;        ///< Cos angle beteeen cylinder axis and edge
    TH1F*                                fCosFirstAxisHist;       ///< Cos angle between vector between centers and first axis
    TH1F*                                f1stTo2ndProjEigenHist;  ///< arc length along vector between centers from first center to cylinder
    
    TH1F*                                fRMaxNextHist;           ///< radius of "cylinder" next cluster
    TH1F*                                fCosMaxNextHist;         ///< Cos angle beteeen cylinder axis and edge
    TH1F*                                fCosNextAxisHist;        ///< Cos angle between vector between centers and next axis
    TH1F*                                f2ndTo1stProjEigenHist;  ///< arc length along vector between centers from next center to cylinder
    
    TH1F*                                fGapBetweenClusHist;     ///< Gap between clusters
    TH1F*                                fGapRatToLenHist;        ///< Ratio of gap to distance between centers

    TH1F*                                fAxesDocaHist;           ///< Closest distance between two primary axes
    TH1F*                                f1stDocaArcLRatHist;     ///< arc length to POCA for DOCA first cluster
    TH1F*                                f2ndDocaArcLRatHist;     ///< arc length to POCA for DOCA second cluster
    
    TH1F*                                f1stTo2ndPosLenRatHist;  ///< Ratio of distance between centers and first proj eigen
    TH1F*                                fGapRatHist;             ///< Ratio of gap and next proj eigen

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
    fNumTransEigens       = pset.get<float>("NumTransEigenVals",    5.    );
    fMaxDOCASeparation    = pset.get<float>("MaxDOCASeparation",    10.   );
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
        
        fNumMergedClusters      = dir.make<TH1F>("NumMergedClus",      "Number Merged",             200,    0., 200.);
        
        f1stTo2ndPosLenHist     = dir.make<TH1F>("1stTo2ndPosLen",     "Distance between Clusters", 250,    0., 500.);
        
        fRMaxFirstHist          = dir.make<TH1F>("rMaxFirst",          "Radius of First Cluster",   200,    0., 100.);
        fCosMaxFirstHist        = dir.make<TH1F>("CosMaxFirst",        "Cos Angle First Cyl/Axis",  200,    0., 1.  );
        fCosFirstAxisHist       = dir.make<TH1F>("CosFirstAxis",       "Cos Angle First Next/Axis", 200,    0., 1.  );
        f1stTo2ndProjEigenHist  = dir.make<TH1F>("1stTo2ndProjEigen",  "Projected Distance First",  200,    0., 200.);
        
        fRMaxNextHist           = dir.make<TH1F>("rMaxNext",           "Radius of Next Cluster",    200,    0., 100.);
        fCosMaxNextHist         = dir.make<TH1F>("CosMaxNext",         "Cos Angle Next Cyl/Axis",   200,    0., 1.  );
        fCosNextAxisHist        = dir.make<TH1F>("CosNextAxis",        "Cos Angle Next Next/Axis",  200,    0., 1.  );
        f2ndTo1stProjEigenHist  = dir.make<TH1F>("2ndTo1stProjEigen",  "Projected Distance Next",   200,    0., 200.);
        
        fGapBetweenClusHist     = dir.make<TH1F>("ClusterGap",         "Gap Between Clusters",      400, -200., 200.);
        fGapRatToLenHist        = dir.make<TH1F>("GapRatToLen",        "Ratio Gap to Distance",     400,  -20.,  20.);

        fAxesDocaHist           = dir.make<TH1F>("AxesDocaHist",       "DOCA",                      200,    0., 200.);
        f1stDocaArcLRatHist     = dir.make<TH1F>("ALenPOCA1Hist",      "Arc Len to POCA 1",         400, -200., 200.);
        f2ndDocaArcLRatHist     = dir.make<TH1F>("ALenPOCA2Hist",      "Arc Len to POCA 2",         400, -200., 200.);
        
        f1stTo2ndPosLenRatHist  = dir.make<TH1F>("1stTo2ndPosLenRat",  "Ratio clus dist to eigen",  200.,   0.,  20.);
        fGapRatHist             = dir.make<TH1F>("GapRat",             "Ratio Gap to Next Eigen",   400,  -20.,  20.);
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
    clusterParametersList.sort([](auto& left, auto& right){return left.getFullPCA().getEigenValues()[2] > right.getFullPCA().getEigenValues()[2];});
    
    // The idea is to continually loop through all clusters until we get to the point where we are no longer doing any merging.
    // If two clusters are merged then we need to recycle through again because it may be that the new cluster can match when
    // it might not have in the past
    size_t                                lastClusterListCount = clusterParametersList.size() + 1;
    reco::ClusterParametersList::iterator lastFirstClusterItr  = clusterParametersList.begin();
    
    int numMergedClusters(0);
    int nOutsideLoops(0);

    while(clusterParametersList.size() != lastClusterListCount)
    {
        // Update the last count
        lastClusterListCount = clusterParametersList.size();
        
        // Keep track of the first cluster iterator each pass through
        reco::ClusterParametersList::iterator firstClusterItr = lastFirstClusterItr++;

        // Loop through the clusters
        while(firstClusterItr != clusterParametersList.end())
        {
            reco::ClusterParameters&              firstClusterParams = *firstClusterItr;
            reco::ClusterParametersList::iterator nextClusterItr     = firstClusterItr;
            
            // Take a brief interlude to do some histogramming
            if (fOutputHistograms)
            {
                const reco::PrincipalComponents& firstPCA = firstClusterParams.getFullPCA();
                
                Eigen::Vector3f firstEigenVals(std::min(1.5 * std::sqrt(std::max(firstPCA.getEigenValues()[0],fMinTransEigenVal)),  50.),
                                               std::min(1.5 * std::sqrt(std::max(firstPCA.getEigenValues()[1],fMinTransEigenVal)), 100.),
                                                        1.5 * std::sqrt(         firstPCA.getEigenValues()[2]));
                
                for(size_t idx = 0; idx < 3; idx++) fFirstEigenValueHists[idx]->Fill(firstEigenVals[idx],1.);
            }

            // Once you get down to the smallest clusters if they haven't already been absorbed there is no need to check them
            if (firstClusterParams.getFullPCA().getEigenValues()[2] < fMinEigenToProcess) break;
            
            while(++nextClusterItr != clusterParametersList.end())
            {
                reco::ClusterParameters& nextClusterParams = *nextClusterItr;
                
                // On any given loop through here it **might** be that the first cluster has been modified. So can't cache
                // the parameters, need the curret ones
                if (linearClusters(firstClusterParams,nextClusterParams))
                {
                    if (mergeClusters(firstClusterParams, nextClusterParams))
                    {
                        // Now remove the "next" cluster
                        nextClusterItr = clusterParametersList.erase(nextClusterItr);
                        
                        // Our new cluster may be larger than those before it so we should try to move backwards to make sure to
                        // give now smaller clusters the chance to get merged as well
                        reco::ClusterParametersList::iterator biggestItr = firstClusterItr;
                        
                        // Step backwards until we hit the beginning of the list
                        while(biggestItr-- != clusterParametersList.begin())
                        {
                            // The list has already been sorted largest to smallest so if we find a cluster that is "larger" then
                            // we have hit the limit
                            if ((*biggestItr).getFullPCA().getEigenValues()[2] > (*firstClusterItr).getFullPCA().getEigenValues()[2])
                            {
                                // Want to insert after the cluster with the larger primary axes.. but there is
                                // no point in doing anything if the new cluster is already in the right spot
                                if (++biggestItr != firstClusterItr) clusterParametersList.splice(biggestItr, clusterParametersList, firstClusterItr);
                                
                                break;
                            }
                        }
                        
                        // Restart loop?
                        nextClusterItr = firstClusterItr;
                        
                        numMergedClusters++;
                    }
                }
            }
            
            firstClusterItr++;
        }
        
        nOutsideLoops++;
    }
    
    std::cout << "==> # merged: " << numMergedClusters << ", in " << nOutsideLoops << " outside loops " << std::endl;
    
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
    float           firstPosToNextPosLen  = firstPosToNextPosVec.norm();
    
    // Now get the first PCA's primary axis and since we'll use them get all of them at once...
    Eigen::Vector3f firstAxis0(firstPCA.getEigenVectors().row(0));
    Eigen::Vector3f firstAxis1(firstPCA.getEigenVectors().row(1));
    Eigen::Vector3f firstAxis2(firstPCA.getEigenVectors().row(2));
    
    // Will want the distance of closest approach of the next cluser's center to the primary, start by finding arc length
    float arcLenToNextDoca = firstPosToNextPosVec.dot(firstAxis2);
    
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
    Eigen::Vector3f firstEigenVals(std::min(1.5 * sqrt(std::max(firstPCA.getEigenValues()[0],fMinTransEigenVal)),  50.),
                                   std::min(1.5 * sqrt(std::max(firstPCA.getEigenValues()[1],fMinTransEigenVal)), 100.),
                                            1.5 * sqrt(         firstPCA.getEigenValues()[2]));
    
    // We treat the PCA as defining a cyclinder of radius given by the quadrature sum of the two transverse eigenvalues
    float rMaxFirst    = fNumTransEigens * std::sqrt(firstEigenVals[0] * firstEigenVals[0] + firstEigenVals[1] * firstEigenVals[1]);
    float cosMaxFirst  = std::max(firstEigenVals[2] / std::sqrt(firstEigenVals[2] * firstEigenVals[2] + rMaxFirst * rMaxFirst),float(0.8));
    float cosFirstAxis = firstAxis2.dot(firstPosToNextPosUnit);
    
    // Now calculate a measure of the length inside the cylider along the vector between the cluster centers
    float firstToNextProjEigen = firstEigenVals[2];
    
    // There are two cases to consider, the first is that the vector exits out the side of the cylinder
    if (cosFirstAxis < cosMaxFirst)
    {
        // In this case we need the sign of the angle between the vector and the cluster primary axis
        float sinFirstAxis = std::sqrt(1. - cosFirstAxis * cosFirstAxis);
        
        // Here the length will be given by the radius, not the primary eigen value
        firstToNextProjEigen = rMaxFirst / sinFirstAxis;
    }
    else firstToNextProjEigen /= cosFirstAxis;
    
    // A brief interlude to fill a few histograms
    if (fOutputHistograms)
    {
        f1stTo2ndPosLenHist->Fill(firstPosToNextPosLen, 1.);
        
        fRMaxFirstHist->Fill(rMaxFirst, 1.);
        fCosMaxFirstHist->Fill(cosMaxFirst, 1.);
        fCosFirstAxisHist->Fill(cosFirstAxis, 1.);
        f1stTo2ndProjEigenHist->Fill(firstToNextProjEigen, 1.);
    }

    // This makes first selection, it should eliminate most of the junk cases:
    if (firstPosToNextPosLen < 3. * firstToNextProjEigen)
    {
        // Recover the axes for the next PCA and make sure pointing convention is observed
        Eigen::Vector3f nextAxis0(nextPCA.getEigenVectors().row(0));
        Eigen::Vector3f nextAxis1(nextPCA.getEigenVectors().row(1));
        Eigen::Vector3f nextAxis2(nextPCA.getEigenVectors().row(2));
        
        // Get the arc length from next to first centers - note convention for PCA orientation
        float arcLenToFirstDoca = -firstPosToNextPosVec.dot(nextAxis2);
        
        // Convention on axis direction applied again
        if (arcLenToFirstDoca > 0.)
        {
            nextAxis0         = -nextAxis0;
            nextAxis1         = -nextAxis1;
            nextAxis2         = -nextAxis2;
            arcLenToFirstDoca = -arcLenToFirstDoca;
        }
        
        // Get the eigen values again
        Eigen::Vector3f nextEigenVals(std::min(1.5 * sqrt(std::max(nextPCA.getEigenValues()[0],fMinTransEigenVal)),  50.),
                                      std::min(1.5 * sqrt(std::max(nextPCA.getEigenValues()[1],fMinTransEigenVal)), 100.),
                                               1.5 * sqrt(         nextPCA.getEigenValues()[2]));
        
        // Repeat the calculation of the length of the vector through the cluster "cylinder"...
        float rMaxNext    = fNumTransEigens * std::sqrt(nextEigenVals[0] * nextEigenVals[0] + nextEigenVals[1] * nextEigenVals[1]);
        float cosMaxNext  = std::max(nextEigenVals[2] / std::sqrt(nextEigenVals[2] * nextEigenVals[2] + rMaxNext * rMaxNext),float(0.8));
        float cosNextAxis = nextAxis2.dot(firstPosToNextPosUnit);
        
        // Now calculate a measure of the length inside the cylider along the vector between the cluster centers
        float nextToFirstProjEigen = nextEigenVals[2];
        
        // There are two cases to consider, the first is that the vector exits out the side of the cylinder
        if (cosNextAxis < cosMaxNext)
        {
            // In this case we need the sign of the angle between the vector and the cluster primary axis
            float sinNextAxis = std::sqrt(1. - cosNextAxis * cosNextAxis);
            
            // Here the length will be given by the radius, not the primary eigen value
            nextToFirstProjEigen = rMaxNext / sinNextAxis;
        }
        else nextToFirstProjEigen /= cosNextAxis;

        // Compute the "gap" between the clusters
        float gapBetweenClusters = firstPosToNextPosLen - firstToNextProjEigen - nextToFirstProjEigen;

        // Allow a generous gap but significantly derate as angle to next cluster increases
        //float nextToFirstScaleFactor = 6. * firstAxis2.dot(nextAxis2);
        float nextToFirstScaleFactor = 3. * firstAxis2.dot(nextAxis2);
        
        nextToFirstScaleFactor = 6.;
        
        // Another brief interlude to fill some histograms
        if (fOutputHistograms)
        {
            for(size_t idx = 0; idx < 3; idx++) fNextEigenValueHists[idx]->Fill(nextEigenVals[idx],1.);
            
            fRMaxNextHist->Fill(rMaxNext, 1.);
            fCosMaxNextHist->Fill(cosMaxNext, 1.);
            fCosNextAxisHist->Fill(cosNextAxis, 1.);
            f2ndTo1stProjEigenHist->Fill(nextToFirstProjEigen, 1.);
            
            fGapBetweenClusHist->Fill(gapBetweenClusters, 1.);
            fGapRatToLenHist->Fill(gapBetweenClusters/firstPosToNextPosLen, 1.);
        }

        // If the gap is negative then the clusters "overlap" (one might be embedded or the PCA axes extend past the space points)
        // We simply check that they are not too far apart
        if (gapBetweenClusters < nextToFirstScaleFactor * nextToFirstProjEigen)
        {
            // Look for the case that we have an embedded cluster
            // An embedded cluster will have its center very close to the center of the surrounding cluster. In most cases its primary axis will be close
            // to that of the surrounding cluster. However, a case to be sensitve to is a small cluster that gets its axis rotated ~90 degrees... so try to
            // stay away from angle cuts...
            // So, we are looking for clusters that lie within the length of the first cluster, with a center that is close to the first cluster's primary
            // axis and which should be fully "contained" in the first cluster.
            
            Eigen::Vector3f nextFirstAxisPOCA = firstCenter + arcLenToNextDoca * firstAxis2;
            Eigen::Vector3f nextDocaVec       = nextCenter - nextFirstAxisPOCA;
            Eigen::Vector3f nextDocaVecUnit   = nextDocaVec.normalized();
            
            // Get a measure of the projection of the three eigen values along the axis from the next center to the first primary axis
            float nextDocaVecProjEigen = nextEigenVals[0] * std::abs(nextDocaVecUnit.dot(nextAxis0))
                                       + nextEigenVals[1] * std::abs(nextDocaVecUnit.dot(nextAxis1))
                                       + nextEigenVals[2] * std::abs(nextDocaVecUnit.dot(nextAxis2));
            
            float nextDocaScaleFactor  = 5.; //6. * firstAxis2.dot(nextAxis2);
            
            // Now we start the search for matching clusters that are essentially side by side
            // These clusters should be in the forward direction (hence an angle cut), arc length along the primary axis to the center should place it
            // "outside" the first cluster and looking back from the next cluster we can't be too far away....
            //else if (cosFirstAxis > cosMaxFirst && std::abs(arcLenToFirstDoca) > 0.8 * nextEigenVals[2] && firstPosToNextPosLen < 10. * nextToFirstProjEigen)
            if (cosFirstAxis > cosMaxFirst && nextDocaVec.norm() < nextDocaScaleFactor * nextDocaVecProjEigen)
            {
                // Now check the distance of closest approach betweent the two vectors
                // Closest approach calculaiton results vectors
                Eigen::Vector3f firstPoca;
                Eigen::Vector3f nextPoca;
                Eigen::Vector3f firstToNextVec;
            
                // Recover the doca of the two axes and their points of closest approach along each axis
                float lineDoca = closestApproach(firstCenter, firstAxis2, nextCenter, nextAxis2, firstPoca, nextPoca, firstToNextVec);
                
                // Scaling factor to increase doca distance as distance grows
                float rMaxScaleFactor = 2. * std::min(firstPosToNextPosLen / firstToNextProjEigen, float(3.));
            
                if (fOutputHistograms)
                {
                    // The below returned the signed arc lengths to their respective pocas
                    float arcLenToFirstPoca = (firstPoca - firstCenter).dot(firstAxis2);
                    float arcLenToNextPoca  = (nextPoca  - nextCenter ).dot(nextAxis2);
                    
                    fAxesDocaHist->Fill(lineDoca, 1.);
                    f1stDocaArcLRatHist->Fill(arcLenToFirstPoca/firstEigenVals[2],1.);
                    f2ndDocaArcLRatHist->Fill(arcLenToNextPoca/nextEigenVals[2],1.);
                }
 
                if (lineDoca < rMaxScaleFactor * rMaxFirst)
                {
                    consistent = true;
                    
                    if (fOutputHistograms)
                    {
                        f1stTo2ndPosLenRatHist->Fill(firstPosToNextPosLen/firstToNextProjEigen, 1.);
                        fGapRatHist->Fill(gapBetweenClusters/nextToFirstProjEigen, 1.);
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
    
    // Make sure lines are not colinear
    if (std::abs(den) > 10. * std::numeric_limits<float>::epsilon())
    {
        float arcLen0 = (b * e - c * d) / den;
        float arcLen1 = (a * e - b * d) / den;
        
        poca0 = P0 + arcLen0 * u0;
        poca1 = P1 + arcLen1 * u1;
    }
    // Otherwise, take the poca to be the distance to the line at the second point
    else
    {
        float arcLenToNextPoint = w0.dot(u0);
        
        poca0 = P0 + arcLenToNextPoint * u0;
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
    
const reco::ClusterHit3D* ClusterMergeAlg::findFurthestHit3D(const Eigen::Vector3f& refPoint, const Eigen::Vector3f& refVector, const reco::HitPairListPtr& hitList) const
{
    const reco::ClusterHit3D* nearestHit3D(hitList.front());
    float                     furthest(-std::numeric_limits<float>::max());
    
    for(const auto& hit3D : hitList)
    {
        Eigen::Vector3f refToHitVec = hit3D->getPosition() - refPoint;
        float           arcLenToHit = refToHitVec.dot(refVector);
        
        if (arcLenToHit > furthest)
        {
            nearestHit3D = hit3D;
            furthest     = arcLenToHit;
        }
    }
    
    return nearestHit3D;
}

    
DEFINE_ART_CLASS_TOOL(ClusterMergeAlg)
} // namespace lar_cluster3d

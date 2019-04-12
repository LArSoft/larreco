/**
 *  @file   ConvexHullPathFinder_tool.cc
 *
 *  @brief  art Tool for comparing clusters and merging those that are consistent
 *
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/IClusterModAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/ConvexHull/ConvexHull.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"

// Eigen
#include <Eigen/Dense>

// Root histograms
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

class ConvexHullPathFinder : virtual public IClusterModAlg
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit ConvexHullPathFinder(const fhicl::ParameterSet&);

    /**
     *  @brief  Destructor
     */
    ~ConvexHullPathFinder();

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

    /**
     *  @brief Use PCA to try to find path in cluster
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    reco::ClusterParametersList::iterator subDivideCluster(reco::ClusterParameters&              cluster,
                                                           reco::PrincipalComponents&            lastPCA,
                                                           reco::ClusterParametersList::iterator positionItr,
                                                           reco::ClusterParametersList&          outputClusterList,
                                                           int                                   level = 0) const;

    using HitOrderTuple     = std::tuple<float,float,reco::ProjectedPoint>;
    using HitOrderTupleList = std::list<HitOrderTuple>;

    bool makeCandidateCluster(Eigen::Vector3f&,
                              reco::ClusterParameters&,
                              reco::HitPairListPtr::iterator,
                              reco::HitPairListPtr::iterator,
                              int) const;

    bool makeCandidateCluster(Eigen::Vector3f&,
                              reco::ClusterParameters&,
                              HitOrderTupleList&,
                              int) const;

    bool completeCandidateCluster(Eigen::Vector3f&, reco::ClusterParameters&, int) const;

    bool breakClusterByKinks(     reco::ClusterParameters&, reco::ClusterParametersList&, int level) const;
    bool breakClusterByKinksTrial(reco::ClusterParameters&, reco::ClusterParametersList&, int level) const;
    bool breakClusterByMaxDefect( reco::ClusterParameters&, reco::ClusterParametersList&, int level) const;
    bool breakClusterInHalf(      reco::ClusterParameters&, reco::ClusterParametersList&, int level) const;
    bool breakClusterAtBigGap(    reco::ClusterParameters&, reco::ClusterParametersList&, int level) const;

    float closestApproach(const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          Eigen::Vector3f&,
                          Eigen::Vector3f&) const;

    using MinMaxPoints    = std::pair<reco::ProjectedPoint,reco::ProjectedPoint>;
    using MinMaxPointPair = std::pair<MinMaxPoints,MinMaxPoints>;

    void buildConvexHull(reco::ClusterParameters& clusterParameters, int level = 0)     const;

    float findConvexHullEndPoints(const reco::EdgeList&, const reco::ClusterHit3D*, const reco::ClusterHit3D*) const;

    using KinkTuple    = std::tuple<int, reco::ConvexHullKinkTuple, HitOrderTupleList, HitOrderTupleList>;
    using KinkTupleVec = std::vector<KinkTuple>;

    void orderHitsAlongEdge(const reco::ProjectedPointList&, const reco::ProjectedPoint&, const Eigen::Vector2f&, HitOrderTupleList&) const;

    void pruneHitOrderTupleLists(HitOrderTupleList&, HitOrderTupleList&) const;

    void fillConvexHullHists(reco::ClusterParameters&, bool) const;

    /**
     *  @brief FHICL parameters
     */
    bool                                        fEnableMonitoring;      ///<
    size_t                                      fMinTinyClusterSize;    ///< Minimum size for a "tiny" cluster
    float                                       fMinGapSize;            ///< Minimum gap size to break at gaps
    float                                       fMinEigen0To1Ratio;     ///< Minimum ratio of eigen 0 to 1 to continue breaking
    float                                       fConvexHullKinkAngle;   ///< Angle to declare a kink in convex hull calc
    float                                       fConvexHullMinSep;      ///< Min hit separation to conisder in convex hull
    mutable float                               fTimeToProcess;         ///<

    /**
     *  @brief Histogram definitions
     */
    bool                                        fFillHistograms;

    TH1F*                                       fTopNum3DHits;
    TH1F*                                       fTopNumEdges;
    TH1F*                                       fTopEigen21Ratio;
    TH1F*                                       fTopEigen20Ratio;
    TH1F*                                       fTopEigen10Ratio;
    TH1F*                                       fTopPrimaryLength;
    TH1F*                                       fTopExtremeSep;
    TH1F*                                       fTopConvexCosEdge;
    TH1F*                                       fTopConvexEdgeLen;

    TH1F*                                       fSubNum3DHits;
    TH1F*                                       fSubNumEdges;
    TH1F*                                       fSubEigen21Ratio;
    TH1F*                                       fSubEigen20Ratio;
    TH1F*                                       fSubEigen10Ratio;
    TH1F*                                       fSubPrimaryLength;
    TH1F*                                       fSubCosToPrevPCA;
    TH1F*                                       fSubCosExtToPCA;
    TH1F*                                       fSubConvexCosEdge;
    TH1F*                                       fSubConvexEdgeLen;
    TH1F*                                       fSubMaxDefect;
    TH1F*                                       fSubUsedDefect;

    /**
     *  @brief Tools
     */
    std::unique_ptr<lar_cluster3d::IClusterAlg> fClusterAlg;            ///<  Algorithm to do 3D space point clustering
    PrincipalComponentsAlg                      fPCAAlg;                // For running Principal Components Analysis
};

ConvexHullPathFinder::ConvexHullPathFinder(fhicl::ParameterSet const &pset) :
    fPCAAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ConvexHullPathFinder::~ConvexHullPathFinder()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConvexHullPathFinder::configure(fhicl::ParameterSet const &pset)
{
    fEnableMonitoring     = pset.get<bool>  ("EnableMonitoring",   true );
    fMinTinyClusterSize   = pset.get<size_t>("MinTinyClusterSize", 40   );
    fMinGapSize           = pset.get<float >("MinClusterGapSize",   2.0 );
    fMinEigen0To1Ratio    = pset.get<float >("MinEigen0To1Ratio",  10.0 );
    fConvexHullKinkAngle  = pset.get<float >("ConvexHullKinkAgle",  0.95);
    fConvexHullMinSep     = pset.get<float >("ConvexHullMinSep",    0.65);
    fClusterAlg           = art::make_tool<lar_cluster3d::IClusterAlg>(pset.get<fhicl::ParameterSet>("ClusterAlg"));

    fTimeToProcess = 0.;

    return;
}

void ConvexHullPathFinder::initializeHistograms(art::TFileDirectory& histDir)
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
    fFillHistograms = true;

    std::string dirName = "ConvexHullPath";

    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());

    // Divide into two sets of hists... those for the overall cluster and
    // those for the subclusters
    fTopNum3DHits     = dir.make<TH1F>("TopNum3DHits",  "Number 3D Hits",  200,    0.,   200.);
    fTopNumEdges      = dir.make<TH1F>("TopNumEdges",   "Number Edges",    200,    0.,   200.);
    fTopEigen21Ratio  = dir.make<TH1F>("TopEigen21Rat", "Eigen 2/1 Ratio", 100,    0.,     1.);
    fTopEigen20Ratio  = dir.make<TH1F>("TopEigen20Rat", "Eigen 2/0 Ratio", 100,    0.,     1.);
    fTopEigen10Ratio  = dir.make<TH1F>("TopEigen10Rat", "Eigen 1/0 Ratio", 100,    0.,     1.);
    fTopPrimaryLength = dir.make<TH1F>("TopPrimaryLen", "Primary Length",  200,    0.,   200.);
    fTopExtremeSep    = dir.make<TH1F>("TopExtremeSep", "Extreme Dist",    200,    0.,   200.);
    fTopConvexCosEdge = dir.make<TH1F>("TopConvexCos",  "CH Edge Cos",     100,   -1.,     1.);
    fTopConvexEdgeLen = dir.make<TH1F>("TopConvexEdge", "CH Edge Len",     200,    0.,    50.);

    fSubNum3DHits     = dir.make<TH1F>("SubNum3DHits",  "Number 3D Hits",  200,    0.,   200.);
    fSubNumEdges      = dir.make<TH1F>("SubNumEdges",   "Number Edges",    200,    0.,   200.);
    fSubEigen21Ratio  = dir.make<TH1F>("SubEigen21Rat", "Eigen 2/1 Ratio", 100,    0.,     1.);
    fSubEigen20Ratio  = dir.make<TH1F>("SubEigen20Rat", "Eigen 2/0 Ratio", 100,    0.,     1.);
    fSubEigen10Ratio  = dir.make<TH1F>("SubEigen10Rat", "Eigen 1/0 Ratio", 100,    0.,     1.);
    fSubPrimaryLength = dir.make<TH1F>("SubPrimaryLen", "Primary Length",  200,    0.,   200.);
    fSubCosToPrevPCA  = dir.make<TH1F>("SubCosToPrev",  "Cos(theta)",      101,    0.,     1.01);
    fSubCosExtToPCA   = dir.make<TH1F>("SubCosExtPCA",  "Cos(theta)",      102,   -1.01,   1.01);
    fSubConvexCosEdge = dir.make<TH1F>("SubConvexCos",  "CH Edge Cos",     100,   -1.,     1.);
    fSubConvexEdgeLen = dir.make<TH1F>("SubConvexEdge", "CH Edge Len",     200,    0.,    50.);
    fSubMaxDefect     = dir.make<TH1F>("SubMaxDefect",  "Max Defect",      100,    0.,    50.);
    fSubUsedDefect    = dir.make<TH1F>("SubUsedDefect", "Used Defect",     100,    0.,    50.);

    return;
}

void ConvexHullPathFinder::ModifyClusters(reco::ClusterParametersList& clusterParametersList) const
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

    // This is the loop over candidate 3D clusters
    // Note that it might be that the list of candidate clusters is modified by splitting
    // So we use the following construct to make sure we get all of them
    reco::ClusterParametersList::iterator clusterParametersListItr = clusterParametersList.begin();

    while(clusterParametersListItr != clusterParametersList.end())
    {
        // Dereference to get the cluster paramters
        reco::ClusterParameters& clusterParameters = *clusterParametersListItr;

        // It turns out that computing the convex hull surrounding the points in the 2D projection onto the
        // plane of largest spread in the PCA is a good way to break up the cluster... and we do it here since
        // we (currently) want this to be part of the standard output
        buildConvexHull(clusterParameters);

        // Make sure our cluster has enough hits...
        if (clusterParameters.getHitPairListPtr().size() > fMinTinyClusterSize)
        {
            // Get an interim cluster list
            reco::ClusterParametersList reclusteredParameters;

            // Call the main workhorse algorithm for building the local version of candidate 3D clusters
            //******** Remind me why we need to call this at this point when the same hits will be used? ********
            //fClusterAlg->Cluster3DHits(clusterParameters.getHitPairListPtr(), reclusteredParameters);
            reclusteredParameters.push_back(clusterParameters);

            // Only process non-empty results
            if (!reclusteredParameters.empty())
            {
                // Loop over the reclustered set
                for (auto& cluster : reclusteredParameters)
                {
                    // It turns out that computing the convex hull surrounding the points in the 2D projection onto the
                    // plane of largest spread in the PCA is a good way to break up the cluster... and we do it here since
                    // we (currently) want this to be part of the standard output
                    buildConvexHull(cluster, 2);

                    // Break our cluster into smaller elements...
                    subDivideCluster(cluster, cluster.getFullPCA(), cluster.daughterList().end(), cluster.daughterList(), 0);

                    // Add the daughters to the cluster
                    clusterParameters.daughterList().insert(clusterParameters.daughterList().end(),cluster);

                    // If filling histograms we do the main cluster here
                    if (fFillHistograms)
                    {
                        reco::PrincipalComponents& fullPCA        = cluster.getFullPCA();
                        std::vector<double>        eigenValVec    = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
                                                                     3. * std::sqrt(fullPCA.getEigenValues()[1]),
                                                                     3. * std::sqrt(fullPCA.getEigenValues()[2])};
                        double                     eigen2To1Ratio = eigenValVec[0] / eigenValVec[1];
                        double                     eigen1To0Ratio = eigenValVec[1] / eigenValVec[2];
                        double                     eigen2To0Ratio = eigenValVec[2] / eigenValVec[2];
                        int                        num3DHits      = cluster.getHitPairListPtr().size();
                        int                        numEdges       = cluster.getConvexHull().getConvexHullEdgeList().size();

                        fTopNum3DHits->Fill(std::min(num3DHits,199), 1.);
                        fTopNumEdges->Fill(std::min(numEdges,199),   1.);
                        fTopEigen21Ratio->Fill(eigen2To1Ratio, 1.);
                        fTopEigen20Ratio->Fill(eigen2To0Ratio, 1.);
                        fTopEigen10Ratio->Fill(eigen1To0Ratio, 1.);
                        fTopPrimaryLength->Fill(std::min(eigenValVec[2],199.), 1.);
//                        fTopExtremeSep->Fill(std::min(edgeLen,199.), 1.);
                        fillConvexHullHists(clusterParameters, true);
                    }
                }
            }
        }

        // Go to next cluster parameters object
        clusterParametersListItr++;
    }

    if (fEnableMonitoring)
    {
        theClockBuildClusters.stop();

        fTimeToProcess = theClockBuildClusters.accumulated_real_time();
    }

    mf::LogDebug("Cluster3D") << ">>>>> Cluster Path finding done" << std::endl;

    return;
}

reco::ClusterParametersList::iterator ConvexHullPathFinder::subDivideCluster(reco::ClusterParameters&              clusterToBreak,
                                                                             reco::PrincipalComponents&            lastPCA,
                                                                             reco::ClusterParametersList::iterator positionItr,
                                                                             reco::ClusterParametersList&          outputClusterList,
                                                                             int                                   level) const
{
    // This is a recursive routine to divide an input cluster, according to the maximum defect point of
    // the convex hull until we reach the point of no further improvement.
    // The assumption is that the input cluster is fully formed already, this routine then simply
    // divides, if successful division into two pieces it then stores the results

    // No point in doing anything if we don't have enough space points
    if (clusterToBreak.getHitPairListPtr().size() > size_t(2 * fMinTinyClusterSize))
    {
        // set an indention string
//        std::string pluses(level/2, '+');
//        std::string indent(level/2, ' ');
//
//        indent += pluses;

        // We want to find the convex hull vertices that lie furthest from the line to/from the extreme points
        // To find these we:
        // 1) recover the extreme points
        // 2) form the vector between them
        // 3) loop through the vertices and keep track of distance to this vector
        // 4) Sort the resulting list by furthest points and select the one we want
        reco::ConvexHull& convexHull = clusterToBreak.getConvexHull();

        reco::ProjectedPointList::const_iterator extremePointListItr = convexHull.getConvexHullExtremePoints().begin();

        const reco::ClusterHit3D*  firstEdgeHit  = std::get<2>(*extremePointListItr++);
        const reco::ClusterHit3D*  secondEdgeHit = std::get<2>(*extremePointListItr);
        Eigen::Vector3f            edgeVec(secondEdgeHit->getPosition()[0] - firstEdgeHit->getPosition()[0],
                                           secondEdgeHit->getPosition()[1] - firstEdgeHit->getPosition()[1],
                                           secondEdgeHit->getPosition()[2] - firstEdgeHit->getPosition()[2]);
//        double                     edgeLen       = edgeVec.norm();

        // normalize it
        edgeVec.normalize();

        // Recover the PCA for the input cluster
        reco::PrincipalComponents& fullPCA     = clusterToBreak.getFullPCA();
        Eigen::Vector3f            fullPrimaryVec(fullPCA.getEigenVectors().row(2));

        // Calculate the doca to the PCA primary axis for each 3D hit
        // Importantly, this also gives us the arclength along that axis to the hit
//        fPCAAlg.PCAAnalysis_calc3DDocas(clusHitPairVector, fullPCA);

        // Sort the hits along the PCA
//        clusHitPairVector.sort([](const auto& left, const auto& right){return left->getArclenToPoca() < right->getArclenToPoca();});

        // Get a temporary  container to hol
        reco::ClusterParametersList tempClusterParametersList;

        // Try breaking clusters by finding the "maximum defect" point.
        // If this fails the fallback in the event of still large clusters is to split in half
        // If starting with the top level cluster then we first try to break at the kinks
        if (level == 0)
        {
            //        if (!breakClusterByMaxDefect(clusterToBreak, clusHitPairVector, tempClusterParametersList, level))
            if (!breakClusterByKinks(clusterToBreak, tempClusterParametersList, level))
            {
                // Look to see if we can break at a gap
                if (!breakClusterAtBigGap(clusterToBreak, tempClusterParametersList, level))
                {
                    // It might be that we have a large deviation in the convex hull...
                    if (!breakClusterByMaxDefect(clusterToBreak, tempClusterParametersList, level))
                    {
                        std::vector<double> eigenValVec = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
                                                           3. * std::sqrt(fullPCA.getEigenValues()[1]),
                                                           3. * std::sqrt(fullPCA.getEigenValues()[2])};

                        // Well, we don't want "flippers" so make sure the edge has some teeth to it
                        //if (edgeLen > 10.) breakClusterInHalf(clusterToBreak, clusHitPairVector, tempClusterParametersList, level);
                        if (eigenValVec[2] > fMinEigen0To1Ratio * eigenValVec[1]) breakClusterInHalf(clusterToBreak, tempClusterParametersList, level);
                    }
                }
            }

        }
        // Otherwise, change the order
        else
        {
            //        if (!breakClusterByMaxDefect(clusterToBreak, clusHitPairVector, tempClusterParametersList, level))
            if (!breakClusterAtBigGap(clusterToBreak, tempClusterParametersList, level))
            {
                // Look to see if we can break at a gap
                if (!breakClusterByKinks(clusterToBreak, tempClusterParametersList, level))
                {
                    // It might be that we have a large deviation in the convex hull...
                    if (!breakClusterByMaxDefect(clusterToBreak, tempClusterParametersList, level))
                    {
                        std::vector<double> eigenValVec = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
                                                           3. * std::sqrt(fullPCA.getEigenValues()[1]),
                                                           3. * std::sqrt(fullPCA.getEigenValues()[2])};

                        // Well, we don't want "flippers" so make sure the edge has some teeth to it
                        //if (edgeLen > 10.) breakClusterInHalf(clusterToBreak, clusHitPairVector, tempClusterParametersList, level);
                        if (eigenValVec[2] > fMinEigen0To1Ratio * eigenValVec[1]) breakClusterInHalf(clusterToBreak, tempClusterParametersList, level);
                    }
                }
            }

        }

        // Can only end with no candidate clusters or two so don't
        for(auto& clusterParams : tempClusterParametersList)
        {
            size_t curOutputClusterListSize = outputClusterList.size();

            positionItr = subDivideCluster(clusterParams, fullPCA, positionItr, outputClusterList, level+4);

            // If the cluster we sent in was successfully broken then the position iterator will be shifted
            // This means we don't want to restore the current cluster here
            if (curOutputClusterListSize < outputClusterList.size()) continue;

            // The current cluster was not further subdivided so we store its info here
            // I think this is where we fill out the rest of the parameters?
            // Start by adding the 2D hits...
            // See if we can avoid duplicates by temporarily transferring to a set
            std::set<const reco::ClusterHit2D*> hitSet;

            // Loop through 3D hits to get a set of unique 2D hits
            for(const auto& hit3D : clusterParams.getHitPairListPtr())
            {
                for(const auto& hit2D : hit3D->getHits())
                {
                    if (hit2D) hitSet.insert(hit2D);
                }
            }

            // Now add these to the new cluster
            for(const auto& hit2D : hitSet)
            {
                hit2D->setStatusBit(reco::ClusterHit2D::USED);
                clusterParams.UpdateParameters(hit2D);
            }

            positionItr = outputClusterList.insert(positionItr,clusterParams);

            // Are we filling histograms
            if (fFillHistograms)
            {
                std::vector<double> eigenValVec = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
                                                   3. * std::sqrt(fullPCA.getEigenValues()[1]),
                                                   3. * std::sqrt(fullPCA.getEigenValues()[2])};

                // Recover the new fullPCA
                reco::PrincipalComponents& newFullPCA = clusterParams.getFullPCA();

                Eigen::Vector3f newPrimaryVec(fullPCA.getEigenVectors().row(2));
                Eigen::Vector3f lastPrimaryVec(newFullPCA.getEigenVectors().row(2));

                int             num3DHits      = clusterParams.getHitPairListPtr().size();
                int             numEdges       = clusterParams.getConvexHull().getConvexHullEdgeList().size();
                float           cosToLast      = newPrimaryVec.dot(lastPrimaryVec);
                double          eigen2To1Ratio = eigenValVec[0] / eigenValVec[1];
                double          eigen1To0Ratio = eigenValVec[1] / eigenValVec[2];
                double          eigen2To0Ratio = eigenValVec[0] / eigenValVec[2];

                fSubNum3DHits->Fill(std::min(num3DHits,199), 1.);
                fSubNumEdges->Fill(std::min(numEdges,199),   1.);
                fSubEigen21Ratio->Fill(eigen2To1Ratio, 1.);
                fSubEigen20Ratio->Fill(eigen2To0Ratio, 1.);
                fSubEigen10Ratio->Fill(eigen1To0Ratio, 1.);
                fSubCosToPrevPCA->Fill(cosToLast, 1.);
                fSubPrimaryLength->Fill(std::min(eigenValVec[2],199.), 1.);
                fSubCosExtToPCA->Fill(fullPrimaryVec.dot(edgeVec), 1.);
            }

            // The above points to the element, want the next element
            positionItr++;
        }
    }

    return positionItr;
}

bool ConvexHullPathFinder::makeCandidateCluster(Eigen::Vector3f&               primaryPCA,
                                                reco::ClusterParameters&       candCluster,
                                                reco::HitPairListPtr::iterator firstHitItr,
                                                reco::HitPairListPtr::iterator lastHitItr,
                                                int                            level) const
{
    std::string indent(level/2, ' ');

    reco::HitPairListPtr& hitPairListPtr = candCluster.getHitPairListPtr();

    // size the container...
    hitPairListPtr.resize(std::distance(firstHitItr,lastHitItr));

    // and copy the hits into the container
    std::copy(firstHitItr,lastHitItr,hitPairListPtr.begin());

    // Will we want to store this cluster?
    bool keepThisCluster(false);

    // Must have a valid pca
    if (completeCandidateCluster(primaryPCA, candCluster, level))
    {
        // Recover the new fullPCA
        reco::PrincipalComponents& newFullPCA = candCluster.getFullPCA();

        // Need to check if the PCA direction has been reversed
        Eigen::Vector3f newPrimaryVec(newFullPCA.getEigenVectors().row(2));

        std::vector<double> eigenValVec      = {3. * std::sqrt(newFullPCA.getEigenValues()[0]),
                                                3. * std::sqrt(newFullPCA.getEigenValues()[1]),
                                                3. * std::sqrt(newFullPCA.getEigenValues()[2])};
        double              cosNewToLast     = std::abs(primaryPCA.dot(newPrimaryVec));
        double              eigen2To1Ratio   = eigenValVec[0] / eigenValVec[1];
        double              eigen1To0Ratio   = eigenValVec[1] / eigenValVec[2];

        // Create a rough cut intended to tell us when we've reached the land of diminishing returns
//        if (candCluster.getBestEdgeList().size() > 4 && cosNewToLast > 0.25 && eigen2To1Ratio < 0.9 && eigen2To0Ratio > 0.001)
        if (candCluster.getConvexHull().getConvexHullEdgeList().size() > 4 && cosNewToLast > 0.25 && eigen2To1Ratio > 0.01 && eigen2To1Ratio < 0.99 && eigen1To0Ratio < 0.5)
        {
            keepThisCluster = true;
        }
    }

    return keepThisCluster;
}

bool ConvexHullPathFinder::makeCandidateCluster(Eigen::Vector3f&         primaryPCA,
                                                reco::ClusterParameters& candCluster,
                                                HitOrderTupleList&       orderedList,
                                                int                      level) const
{
    std::string indent(level/2, ' ');

    reco::HitPairListPtr& hitPairListPtr = candCluster.getHitPairListPtr();

    // Fill the list the old fashioned way...
    for(const auto& tupleVal : orderedList) hitPairListPtr.emplace_back(std::get<2>(std::get<2>(tupleVal)));

    // Will we want to store this cluster?
    bool keepThisCluster(false);

    // Must have a valid pca
    if (completeCandidateCluster(primaryPCA, candCluster, level))
    {
        // Recover the new fullPCA
        reco::PrincipalComponents& newFullPCA = candCluster.getFullPCA();

        // Need to check if the PCA direction has been reversed
        Eigen::Vector3f newPrimaryVec(newFullPCA.getEigenVectors().row(2));

        std::vector<double> eigenValVec      = {3. * std::sqrt(newFullPCA.getEigenValues()[0]),
                                                3. * std::sqrt(newFullPCA.getEigenValues()[1]),
                                                3. * std::sqrt(newFullPCA.getEigenValues()[2])};
        double              cosNewToLast     = std::abs(primaryPCA.dot(newPrimaryVec));
        double              eigen2To1Ratio   = eigenValVec[0] / eigenValVec[1];
        double              eigen1To0Ratio   = eigenValVec[1] / eigenValVec[2];

        // Create a rough cut intended to tell us when we've reached the land of diminishing returns
        //        if (candCluster.getBestEdgeList().size() > 4 && cosNewToLast > 0.25 && eigen2To1Ratio < 0.9 && eigen2To0Ratio > 0.001)
        if (candCluster.getBestEdgeList().size() > 4 && cosNewToLast > 0.25 && eigen2To1Ratio > 0.01 && eigen2To1Ratio < 0.99 && eigen1To0Ratio < 0.5)
        {
            keepThisCluster = true;
        }
    }

    return keepThisCluster;
}

bool ConvexHullPathFinder::completeCandidateCluster(Eigen::Vector3f& primaryPCA, reco::ClusterParameters& candCluster, int level) const
{
    // First stage of feature extraction runs here
    fPCAAlg.PCAAnalysis_3D(candCluster.getHitPairListPtr(), candCluster.getFullPCA());

    // Recover the new fullPCA
    reco::PrincipalComponents& newFullPCA = candCluster.getFullPCA();

    // Will we want to store this cluster?
    bool keepThisCluster(false);

    // Must have a valid pca
    if (newFullPCA.getSvdOK())
    {
        // Need to check if the PCA direction has been reversed
        Eigen::Vector3f newPrimaryVec(newFullPCA.getEigenVectors().row(2));

        // If the PCA's are opposite the flip the axes
        if (primaryPCA.dot(newPrimaryVec) < 0.)
        {
            for(size_t vecIdx = 0; vecIdx < 3; vecIdx++) newFullPCA.flipAxis(vecIdx);
        }

        // Set the skeleton PCA to make sure it has some value
        candCluster.getSkeletonPCA() = candCluster.getFullPCA();

        // Be sure to compute the oonvex hull surrounding the now broken cluster
        buildConvexHull(candCluster, level+2);

        keepThisCluster = true;
    }

    return keepThisCluster;
}

bool ConvexHullPathFinder::breakClusterByKinks(reco::ClusterParameters& clusterToBreak, reco::ClusterParametersList& outputClusterList, int level) const
{
    // Set up container to keep track of edges
    using HitKinkTuple    = std::tuple<int, reco::HitPairListPtr::iterator>;
    using HitKinkTupleVec = std::vector<HitKinkTuple>;

    // Recover our hits
    reco::HitPairListPtr& hitList = clusterToBreak.getHitPairListPtr();

    // Set up container to keep track of edges
    HitKinkTupleVec kinkTupleVec;

    reco::ConvexHull&              convexHull    = clusterToBreak.getConvexHull();
    reco::ConvexHullKinkTupleList& kinkPointList = convexHull.getConvexHullKinkPoints();

    for(auto& kink : kinkPointList)
    {
        const reco::ClusterHit3D* hit3D = std::get<2>(std::get<0>(kink));

        reco::HitPairListPtr::iterator kinkItr = std::find(hitList.begin(),hitList.end(),hit3D);

        if (kinkItr == hitList.end()) continue;

        int numStartToKink = std::distance(hitList.begin(),kinkItr);
        int numKinkToEnd   = std::distance(kinkItr, hitList.end());
        int minNumHits     = std::min(numStartToKink,numKinkToEnd);

        if (minNumHits > int(fMinTinyClusterSize)) kinkTupleVec.emplace_back(minNumHits,kinkItr);
    }

    // No work if the list is empty
    if (!kinkTupleVec.empty())
    {
        std::sort(kinkTupleVec.begin(),kinkTupleVec.end(),[](const auto& left,const auto& right){return std::get<0>(left) > std::get<0>(right);});

        // Recover the kink point
        reco::HitPairListPtr::iterator kinkItr = std::get<1>(kinkTupleVec.front());

        // Set up to split the input cluster
        outputClusterList.push_back(reco::ClusterParameters());

        reco::ClusterParameters& clusterParams1  = outputClusterList.back();

        reco::PrincipalComponents& fullPCA(clusterToBreak.getFullPCA());
        Eigen::Vector3f            fullPrimaryVec(fullPCA.getEigenVectors().row(2));

        if (makeCandidateCluster(fullPrimaryVec, clusterParams1, hitList.begin(), kinkItr, level))
        {
            outputClusterList.push_back(reco::ClusterParameters());

            reco::ClusterParameters& clusterParams2  = outputClusterList.back();

            makeCandidateCluster(fullPrimaryVec, clusterParams2, kinkItr, hitList.end(), level);

            if (fFillHistograms)
            {
                fillConvexHullHists(clusterParams1, false);
                fillConvexHullHists(clusterParams2, false);
            }
        }

        // If we did not make 2 clusters then be sure to clear the output list
        if (outputClusterList.size() != 2) outputClusterList.clear();
    }

    return !outputClusterList.empty();
}

bool ConvexHullPathFinder::breakClusterByKinksTrial(reco::ClusterParameters& clusterToBreak, reco::ClusterParametersList& outputClusterList, int level) const
{
    // Set up container to keep track of edges
    KinkTupleVec kinkTupleVec;

    reco::ConvexHull&              convexHull    = clusterToBreak.getConvexHull();
    reco::ProjectedPointList&      pointList     = convexHull.getProjectedPointList();
    reco::ConvexHullKinkTupleList& kinkPointList = convexHull.getConvexHullKinkPoints();

    for(auto& kink : kinkPointList)
    {
        // Make an instance of the vec value to avoid copying if we can...
        kinkTupleVec.push_back(KinkTuple());

        KinkTuple& kinkTuple = kinkTupleVec.back();

        std::get<1>(kinkTuple) = kink;

        // Recover vectors, want them pointing away from intersection point
        Eigen::Vector2f    firstEdge  = -std::get<1>(kink);
        HitOrderTupleList& firstList  = std::get<2>(kinkTuple);
        HitOrderTupleList& secondList = std::get<3>(kinkTuple);

        orderHitsAlongEdge(pointList, std::get<0>(kink), firstEdge, firstList);

        if (firstList.size() > fMinTinyClusterSize)
        {
            Eigen::Vector2f secondEdge = std::get<2>(kink);

            orderHitsAlongEdge(pointList, std::get<0>(kink), secondEdge, secondList);

            if (secondList.size() > fMinTinyClusterSize)
                std::get<0>(kinkTuple) = std::min(firstList.size(),secondList.size());
        }

        // Special handling...
        if (firstList.size() + secondList.size() > pointList.size())
        {
            if (firstList.size() > secondList.size()) pruneHitOrderTupleLists(firstList,secondList);
            else                                      pruneHitOrderTupleLists(secondList,firstList);

            std::get<0>(kinkTuple) = std::min(firstList.size(),secondList.size());
        }

        if (std::get<0>(kinkTuple) < int(fMinTinyClusterSize)) kinkTupleVec.pop_back();
    }

    // No work if the list is empty
    if (!kinkTupleVec.empty())
    {
        // If more than one then want the kink with the most elements both sizes
        std::sort(kinkTupleVec.begin(),kinkTupleVec.end(),[](const auto& left,const auto& right){return std::get<0>(left) > std::get<0>(right);});

        // Recover the kink point
        KinkTuple& kinkTuple = kinkTupleVec.front();

        // Set up to split the input cluster
        outputClusterList.push_back(reco::ClusterParameters());

        reco::ClusterParameters& clusterParams1  = outputClusterList.back();

        reco::PrincipalComponents& fullPCA(clusterToBreak.getFullPCA());
        Eigen::Vector3f            fullPrimaryVec(fullPCA.getEigenVectors().row(2));

        if (makeCandidateCluster(fullPrimaryVec, clusterParams1, std::get<2>(kinkTuple), level))
        {
            outputClusterList.push_back(reco::ClusterParameters());

            reco::ClusterParameters& clusterParams2  = outputClusterList.back();

            makeCandidateCluster(fullPrimaryVec, clusterParams2, std::get<3>(kinkTuple), level);

            if (fFillHistograms)
            {
                fillConvexHullHists(clusterParams1, false);
                fillConvexHullHists(clusterParams2, false);
            }
        }

        // If we did not make 2 clusters then be sure to clear the output list
        if (outputClusterList.size() != 2) outputClusterList.clear();
    }

    return !outputClusterList.empty();
}

void ConvexHullPathFinder::orderHitsAlongEdge(const reco::ProjectedPointList& hitList,
                                              const reco::ProjectedPoint& point,
                                              const Eigen::Vector2f& edge,
                                              HitOrderTupleList& orderedList) const
{
    // Use the input kink point as the start point of the edge
    Eigen::Vector2f kinkPos(std::get<0>(point),std::get<1>(point));

    // Loop over the input hits
    for (const auto& hit : hitList)
    {
        // Now we need to calculate the doca and poca...
        // Start by getting this hits position
        Eigen::Vector2f hitPos(std::get<0>(hit),std::get<1>(hit));

        // Form a TVector from this to the cluster average position
        Eigen::Vector2f hitToKinkVec = hitPos - kinkPos;

        // With this we can get the arclength to the doca point
        float arcLenToPoca = hitToKinkVec.dot(edge);

        // Require the hit to not be past the kink point
        if (arcLenToPoca < 0.) continue;

        // Get the coordinates along the axis for this point
        Eigen::Vector2f pocaPos = kinkPos + arcLenToPoca * edge;

        // Now get doca and poca
        Eigen::Vector2f pocaPosToHitPos = hitPos - pocaPos;
        float           pocaToAxis      = pocaPosToHitPos.norm();

        std::cout << "-- arcLenToPoca: " << arcLenToPoca << ", doca: " << pocaToAxis << std::endl;

        orderedList.emplace_back(arcLenToPoca,pocaToAxis,hit);
    }

    // Sort the list in order of ascending distance from the kink point
    orderedList.sort([](const auto& left,const auto& right){return std::get<0>(left) < std::get<0>(right);});

    return;
}

void ConvexHullPathFinder::pruneHitOrderTupleLists(HitOrderTupleList& shortList, HitOrderTupleList& longList) const
{
    // Assume the first list is the short one, so we loop through the elements of that list..
    HitOrderTupleList::iterator shortItr = shortList.begin();

    while(shortItr != shortList.end())
    {
        // Recover the search key
        const reco::ClusterHit3D* hit3D = std::get<2>(std::get<2>(*shortItr));

        // Ok, find the corresponding point in the other list...
        HitOrderTupleList::iterator longItr = std::find_if(longList.begin(),longList.end(),[&hit3D](const auto& elem){return hit3D == std::get<2>(std::get<2>(elem));});

        if (longItr != longList.end())
        {
            if (std::get<1>(*longItr) < std::get<1>(*shortItr))
            {
                shortItr = shortList.erase(shortItr);
            }
            else
            {
                longItr = longList.erase(longItr);
                shortItr++;
            }
        }
        else shortItr++;
    }


    return;
}

bool ConvexHullPathFinder::breakClusterByMaxDefect(reco::ClusterParameters& clusterToBreak, reco::ClusterParametersList& outputClusterList, int level) const
{
    // Set up container to keep track of edges
    using DistEdgeTuple    = std::tuple<float, const reco::EdgeTuple*>;
    using DistEdgeTupleVec = std::vector<DistEdgeTuple>;

    DistEdgeTupleVec distEdgeTupleVec;

    reco::ProjectedPointList::const_iterator extremePointListItr = clusterToBreak.getConvexHull().getConvexHullExtremePoints().begin();

    const reco::ClusterHit3D*  firstEdgeHit  = std::get<2>(*extremePointListItr++);
    const reco::ClusterHit3D*  secondEdgeHit = std::get<2>(*extremePointListItr);
    Eigen::Vector3f            edgeVec(secondEdgeHit->getPosition()[0] - firstEdgeHit->getPosition()[0],
                                       secondEdgeHit->getPosition()[1] - firstEdgeHit->getPosition()[1],
                                       secondEdgeHit->getPosition()[2] - firstEdgeHit->getPosition()[2]);
    double                     edgeLen       = edgeVec.norm();

    // normalize it
    edgeVec.normalize();

    // Now loop through all the edges and search for the furthers point
    for(const auto& edge : clusterToBreak.getBestEdgeList())
    {
        const reco::ClusterHit3D* nextEdgeHit = std::get<0>(edge);  // recover the first point

        // Create vector to this point from the longest edge
        Eigen::Vector3f hitToEdgeVec(nextEdgeHit->getPosition()[0] - firstEdgeHit->getPosition()[0],
                                     nextEdgeHit->getPosition()[1] - firstEdgeHit->getPosition()[1],
                                     nextEdgeHit->getPosition()[2] - firstEdgeHit->getPosition()[2]);

        // Get projection
        float hitProjection = hitToEdgeVec.dot(edgeVec);

        // Require that the point is really "opposite" the longest edge
        if (hitProjection > 0. && hitProjection < edgeLen)
        {
            Eigen::Vector3f distToHitVec = hitToEdgeVec - hitProjection * edgeVec;
            float           distToHit    = distToHitVec.norm();

            distEdgeTupleVec.emplace_back(distToHit,&edge);
        }
    }

    std::sort(distEdgeTupleVec.begin(),distEdgeTupleVec.end(),[](const auto& left,const auto& right){return std::get<0>(left) > std::get<0>(right);});

    reco::PrincipalComponents& fullPCA(clusterToBreak.getFullPCA());
    Eigen::Vector3f            fullPrimaryVec(fullPCA.getEigenVectors().row(2));

    // Recover our hits
    reco::HitPairListPtr& hitList = clusterToBreak.getHitPairListPtr();

    // Calculate the doca to the PCA primary axis for each 3D hit
    // Importantly, this also gives us the arclength along that axis to the hit
    fPCAAlg.PCAAnalysis_calc3DDocas(hitList, fullPCA);

    // Sort the hits along the PCA
    hitList.sort([](const auto& left, const auto& right){return left->getArclenToPoca() < right->getArclenToPoca();});

    // Get a temporary  container to hol
    float usedDefectDist(0.);

    for(const auto& distEdgeTuple : distEdgeTupleVec)
    {
        const reco::EdgeTuple&    edgeTuple = *std::get<1>(distEdgeTuple);
        const reco::ClusterHit3D* edgeHit   = std::get<0>(edgeTuple);

        usedDefectDist = std::get<0>(distEdgeTuple);

        // Now find the hit identified above as furthest away
        reco::HitPairListPtr::iterator vertexItr = std::find(hitList.begin(),hitList.end(),edgeHit);

        // Make sure enough hits either side, otherwise we just keep the current cluster
        if (vertexItr == hitList.end() || std::distance(hitList.begin(),vertexItr) < int(fMinTinyClusterSize) || std::distance(vertexItr,hitList.end()) < int(fMinTinyClusterSize)) continue;

        outputClusterList.push_back(reco::ClusterParameters());

        reco::ClusterParameters& clusterParams1  = outputClusterList.back();

        if (makeCandidateCluster(fullPrimaryVec, clusterParams1, hitList.begin(), vertexItr, level))
        {
            outputClusterList.push_back(reco::ClusterParameters());

            reco::ClusterParameters& clusterParams2  = outputClusterList.back();

            if (makeCandidateCluster(fullPrimaryVec, clusterParams2, vertexItr, hitList.end(), level))
            {
                if (fFillHistograms)
                {
                    fSubMaxDefect->Fill(std::get<0>(distEdgeTupleVec.front()), 1.);
                    fSubUsedDefect->Fill(usedDefectDist, 1.);
                }
                break;
            }
        }

        // If here then we could not make two valid clusters and so we try again
        outputClusterList.clear();
    }

   return !outputClusterList.empty();
}

bool ConvexHullPathFinder::breakClusterInHalf(reco::ClusterParameters& clusterToBreak, reco::ClusterParametersList& outputClusterList, int level) const
{
    reco::PrincipalComponents& fullPCA(clusterToBreak.getFullPCA());
    Eigen::Vector3f            fullPrimaryVec(fullPCA.getEigenVectors().row(2));

    // Recover our hits
    reco::HitPairListPtr& hitList = clusterToBreak.getHitPairListPtr();

    // Calculate the doca to the PCA primary axis for each 3D hit
    // Importantly, this also gives us the arclength along that axis to the hit
    fPCAAlg.PCAAnalysis_calc3DDocas(hitList, fullPCA);

    // Sort the hits along the PCA
    hitList.sort([](const auto& left, const auto& right){return left->getArclenToPoca() < right->getArclenToPoca();});

    reco::HitPairListPtr::iterator vertexItr = hitList.begin();

    std::advance(vertexItr, hitList.size()/2);

    outputClusterList.push_back(reco::ClusterParameters());

    reco::ClusterParameters& clusterParams1  = outputClusterList.back();

    if (makeCandidateCluster(fullPrimaryVec, clusterParams1, hitList.begin(), vertexItr, level))
    {
        outputClusterList.push_back(reco::ClusterParameters());

        reco::ClusterParameters& clusterParams2  = outputClusterList.back();

        makeCandidateCluster(fullPrimaryVec, clusterParams2, vertexItr, hitList.end(), level);
    }

    if (outputClusterList.size() != 2) outputClusterList.clear();

    return !outputClusterList.empty();
}

bool ConvexHullPathFinder::breakClusterAtBigGap(reco::ClusterParameters& clusterToBreak, reco::ClusterParametersList& outputClusterList, int level) const
{
    // Idea here is to scan the input hit list (assumed ordered along the current PCA) and look for "large" gaps
    // Here a gap is determined when the hits were ordered by their distance along the primary PCA to their doca to it.
    reco::PrincipalComponents& fullPCA(clusterToBreak.getFullPCA());

    // Recover our hits
    reco::HitPairListPtr& hitList = clusterToBreak.getHitPairListPtr();

    // Calculate the doca to the PCA primary axis for each 3D hit
    // Importantly, this also gives us the arclength along that axis to the hit
    fPCAAlg.PCAAnalysis_calc3DDocas(hitList, fullPCA);

    // Sort the hits along the PCA
    hitList.sort([](const auto& left, const auto& right){return left->getArclenToPoca() < right->getArclenToPoca();});

    // Loop through the input hit list and keep track of first hit of largest gap
    reco::HitPairListPtr::iterator bigGapHitItr = hitList.begin();
    float                                biggestGap   = 0.;

    reco::HitPairListPtr::iterator lastHitItr = hitList.begin();

    for(reco::HitPairListPtr::iterator hitItr = hitList.begin(); hitItr != hitList.end(); hitItr++)
    {
        float currentGap = std::abs((*hitItr)->getArclenToPoca() - (*lastHitItr)->getArclenToPoca());

        if (currentGap > biggestGap)
        {
            bigGapHitItr = hitItr;      // Note that this is an iterator and will be the "end" going from begin, and "begin" for second half
            biggestGap   = currentGap;
        }

        lastHitItr = hitItr;
    }

    // Require some minimum gap size...
    if (biggestGap > fMinGapSize)
    {
        outputClusterList.push_back(reco::ClusterParameters());

        reco::ClusterParameters& clusterParams1  = outputClusterList.back();

        reco::PrincipalComponents& fullPCA(clusterToBreak.getFullPCA());
        Eigen::Vector3f            fullPrimaryVec(fullPCA.getEigenVectors().row(2));

        if (makeCandidateCluster(fullPrimaryVec, clusterParams1, hitList.begin(), bigGapHitItr, level))
        {
            outputClusterList.push_back(reco::ClusterParameters());

            reco::ClusterParameters& clusterParams2  = outputClusterList.back();

            makeCandidateCluster(fullPrimaryVec, clusterParams2, bigGapHitItr, hitList.end(), level);
        }

        if (outputClusterList.size() != 2) outputClusterList.clear();
    }

    return !outputClusterList.empty();
}


void ConvexHullPathFinder::buildConvexHull(reco::ClusterParameters& clusterParameters, int level) const
{
    // set an indention string
    std::string minuses(level/2, '-');
    std::string indent(level/2, ' ');

    indent += minuses;

    // The plan is to build the enclosing 2D polygon around the points in the PCA plane of most spread for this cluster
    // To do so we need to start by building a list of 2D projections onto the plane of most spread...
    reco::PrincipalComponents& pca = clusterParameters.getFullPCA();

    // Recover the parameters from the Principal Components Analysis that we need to project and accumulate
    const Eigen::Vector3f&    pcaCenter  = pca.getAvePosition();
    reco::ConvexHull&         convexHull = clusterParameters.getConvexHull();
    reco::ProjectedPointList& pointList  = convexHull.getProjectedPointList();

    // Loop through hits and do projection to plane
    for(const auto& hit3D : clusterParameters.getHitPairListPtr())
    {
        Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
                                    hit3D->getPosition()[1] - pcaCenter(1),
                                    hit3D->getPosition()[2] - pcaCenter(2));
        Eigen::Vector3f pcaToHit = pca.getEigenVectors() * pcaToHitVec;

        pointList.emplace_back(dcel2d::Point(pcaToHit(1),pcaToHit(2),hit3D));
    }

    // Sort the point vec by increasing x, then increase y
    pointList.sort([](const auto& left, const auto& right){return (std::abs(std::get<0>(left) - std::get<0>(right)) > std::numeric_limits<float>::epsilon()) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);});

    // containers for finding the "best" hull...
    std::vector<ConvexHull>               convexHullVec;
    std::vector<reco::ProjectedPointList> rejectedListVec;
    bool                                  increaseDepth(pointList.size() > 3);
    float                                 lastArea(std::numeric_limits<float>::max());

    while(increaseDepth)
    {
        // Get another convexHull container
        convexHullVec.push_back(ConvexHull(pointList, fConvexHullKinkAngle, fConvexHullMinSep));
        rejectedListVec.push_back(reco::ProjectedPointList());

        const ConvexHull&               convexHull       = convexHullVec.back();
        reco::ProjectedPointList&       rejectedList     = rejectedListVec.back();
        const reco::ProjectedPointList& convexHullPoints = convexHull.getConvexHull();

        increaseDepth = false;

        if (convexHull.getConvexHullArea() > 0.)
        {
            if (convexHullVec.size() < 2 || convexHull.getConvexHullArea() < 0.8 * lastArea)
            {
                for(auto& point : convexHullPoints)
                {
                    pointList.remove(point);
                    rejectedList.emplace_back(point);
                }
                lastArea      = convexHull.getConvexHullArea();
//                increaseDepth = true;
            }
        }
    }

    // do we have a valid convexHull?
    while(!convexHullVec.empty() && convexHullVec.back().getConvexHullArea() < 0.5)
    {
        convexHullVec.pop_back();
        rejectedListVec.pop_back();
    }

    // If we found the convex hull then build edges around the region
    if (!convexHullVec.empty())
    {
        size_t               nRejectedTotal(0);
        reco::HitPairListPtr hitPairListPtr = clusterParameters.getHitPairListPtr();

        for(const auto& rejectedList : rejectedListVec)
        {
            nRejectedTotal += rejectedList.size();

            for(const auto& rejectedPoint : rejectedList)
            {
                if (convexHullVec.back().findNearestDistance(rejectedPoint) > 0.5)
                    hitPairListPtr.remove(std::get<2>(rejectedPoint));
            }
        }

        // Now add "edges" to the cluster to describe the convex hull (for the display)
        reco::ProjectedPointList& convexHullPointList = convexHull.getConvexHullPointList();
        reco::Hit3DToEdgeMap&     edgeMap             = convexHull.getConvexHullEdgeMap();
        reco::EdgeList&           edgeList            = convexHull.getConvexHullEdgeList();

        reco::ProjectedPoint lastPoint = convexHullVec.back().getConvexHull().front();

        for(auto& curPoint : convexHullVec.back().getConvexHull())
        {
            if (curPoint == lastPoint) continue;

            const reco::ClusterHit3D* lastPoint3D = std::get<2>(lastPoint);
            const reco::ClusterHit3D* curPoint3D  = std::get<2>(curPoint);

            float distBetweenPoints = (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) * (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0])
                                    + (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) * (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1])
                                    + (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]) * (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]);

            distBetweenPoints = std::sqrt(distBetweenPoints);

            reco::EdgeTuple edge(lastPoint3D,curPoint3D,distBetweenPoints);

            convexHullPointList.push_back(curPoint);
            edgeMap[lastPoint3D].push_back(edge);
            edgeMap[curPoint3D].push_back(edge);
            edgeList.emplace_back(edge);

            lastPoint = curPoint;
        }

        // Store the "extreme" points
        const ConvexHull::PointList& extremePoints    = convexHullVec.back().getExtremePoints();
        reco::ProjectedPointList&    extremePointList = convexHull.getConvexHullExtremePoints();

        for(const auto& point : extremePoints) extremePointList.push_back(point);

        // Store the "kink" points
        const reco::ConvexHullKinkTupleList& kinkPoints    = convexHullVec.back().getKinkPoints();
        reco::ConvexHullKinkTupleList&       kinkPointList = convexHull.getConvexHullKinkPoints();

        for(const auto& kink : kinkPoints) kinkPointList.push_back(kink);
    }

    return;
}

void ConvexHullPathFinder::fillConvexHullHists(reco::ClusterParameters& clusterParameters, bool top) const
{
    reco::ProjectedPointList& convexHullPoints = clusterParameters.getConvexHull().getConvexHullPointList();

    if (convexHullPoints.size() > 2)
    {
        reco::ProjectedPointList::iterator pointItr = convexHullPoints.begin();

        // Advance to the second to last element
        std::advance(pointItr, convexHullPoints.size() - 2);

        reco::ProjectedPoint lastPoint = *pointItr++;

        // Reset pointer to the first element
        pointItr = convexHullPoints.begin();

        reco::ProjectedPoint curPoint = *pointItr++;
        Eigen::Vector2f      lastEdge(std::get<0>(curPoint) - std::get<0>(lastPoint), std::get<1>(curPoint) - std::get<1>(lastPoint));

        lastEdge.normalize();

        while(pointItr != convexHullPoints.end())
        {
            reco::ProjectedPoint& nextPoint = *pointItr++;

            Eigen::Vector2f nextEdge(std::get<0>(nextPoint) - std::get<0>(curPoint), std::get<1>(nextPoint) - std::get<1>(curPoint));
            float           nextEdgeLen = nextEdge.norm();

            nextEdge.normalize();

            float cosLastNextEdge = lastEdge.dot(nextEdge);

            if (top)
            {
                fTopConvexCosEdge->Fill(cosLastNextEdge, 1.);
                fTopConvexEdgeLen->Fill(std::min(nextEdgeLen,float(49.9)), 1.);
            }
            else
            {
                fSubConvexCosEdge->Fill(cosLastNextEdge, 1.);
                fSubConvexEdgeLen->Fill(std::min(nextEdgeLen,float(49.9)), 1.);
            }

            if (nextEdgeLen > fConvexHullMinSep) lastEdge = nextEdge;

            curPoint = nextPoint;
        }
    }

    return;
}

float ConvexHullPathFinder::closestApproach(const Eigen::Vector3f& P0,
                                            const Eigen::Vector3f& u0,
                                            const Eigen::Vector3f& P1,
                                            const Eigen::Vector3f& u1,
                                            Eigen::Vector3f&       poca0,
                                            Eigen::Vector3f&       poca1) const
{
    // Technique is to compute the arclength to each point of closest approach
    Eigen::Vector3f w0 = P0 - P1;
    float a(1.);
    float b(u0.dot(u1));
    float c(1.);
    float d(u0.dot(w0));
    float e(u1.dot(w0));
    float den(a * c - b * b);

    float arcLen0 = (b * e - c * d) / den;
    float arcLen1 = (a * e - b * d) / den;

    poca0 = P0 + arcLen0 * u0;
    poca1 = P1 + arcLen1 * u1;

    return (poca0 - poca1).norm();
}

float ConvexHullPathFinder::findConvexHullEndPoints(const reco::EdgeList& convexHull, const reco::ClusterHit3D* first3D, const reco::ClusterHit3D* last3D) const
{
    float largestDistance(0.);

    // Note that edges are vectors and that the convex hull edge list will be ordered
    // The idea is that the maximum distance from a given edge is to the edge just before the edge that "turns back" towards the current edge
    // meaning that the dot product of the two edges becomes negative.
    reco::EdgeList::const_iterator firstEdgeItr = convexHull.begin();

    while(firstEdgeItr != convexHull.end())
    {
        reco::EdgeList::const_iterator nextEdgeItr = firstEdgeItr;

//        Eigen::Vector2f firstEdgeVec(std::get<3>(*firstEdgeItr),std::get<);
//        Eigen::Vector2f lastPrimaryVec(lastPCA.getEigenVectors()[0][0],lastPCA.getEigenVectors()[0][1],lastPCA.getEigenVectors()[0][2]);
//        float           cosToLast = newPrimaryVec.dot(lastPrimaryVec);

        while(++nextEdgeItr != convexHull.end())
        {

        }
    }

    return largestDistance;
}

DEFINE_ART_CLASS_TOOL(ConvexHullPathFinder)
} // namespace lar_cluster3d

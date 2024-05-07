/**
 *  @file   VoronoiPathFinder_tool.cc
 *
 *  @brief  art Tool for comparing clusters and merging those that are consistent
 *
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileDirectory.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/RecoAlg/Cluster3DAlgs/ConvexHull/ConvexHull.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterModAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/Voronoi.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"

// Eigen
#include <Eigen/Core>

// Root histograms
#include "TH1F.h"

// std includes
#include <iostream>
#include <memory>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

  class VoronoiPathFinder : virtual public IClusterModAlg {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit VoronoiPathFinder(const fhicl::ParameterSet&);

    /**
     *  @brief  Destructor
     */
    ~VoronoiPathFinder();

    void configure(fhicl::ParameterSet const& pset) override;

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
    float getTimeToExecute() const override { return fTimeToProcess; }

  private:
    /**
     *  @brief Use PCA to try to find path in cluster
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    reco::ClusterParametersList::iterator breakIntoTinyBits(
      reco::ClusterParameters& cluster,
      reco::PrincipalComponents& lastPCA,
      reco::ClusterParametersList::iterator positionItr,
      reco::ClusterParametersList& outputClusterList,
      int level = 0) const;

    /**
     *  @brief Use PCA to try to find path in cluster
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    reco::ClusterParametersList::iterator subDivideCluster(
      reco::ClusterParameters& cluster,
      reco::ClusterParametersList::iterator positionItr,
      reco::ClusterParametersList& outputClusterList,
      int level = 0) const;

    bool makeCandidateCluster(Eigen::Vector3f&,
                              reco::ClusterParameters&,
                              reco::HitPairListPtr::iterator,
                              reco::HitPairListPtr::iterator,
                              int) const;

    float closestApproach(const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          Eigen::Vector3f&,
                          Eigen::Vector3f&) const;

    using Point = std::tuple<float, float, const reco::ClusterHit3D*>;
    using PointList = std::list<Point>;
    using MinMaxPoints = std::pair<Point, Point>;
    using MinMaxPointPair = std::pair<MinMaxPoints, MinMaxPoints>;

    void buildConvexHull(reco::ClusterParameters& clusterParameters, int level = 0) const;
    void buildVoronoiDiagram(reco::ClusterParameters& clusterParameters) const;

    /**
     *  @brief FHICL parameters
     */
    bool fEnableMonitoring;       ///<
    size_t fMinTinyClusterSize;   ///< Minimum size for a "tiny" cluster
    mutable float fTimeToProcess; ///<

    /**
     *  @brief Histogram definitions
     */
    bool fFillHistograms;

    TH1F* fTopNum3DHits;
    TH1F* fTopNumEdges;
    TH1F* fTopEigen21Ratio;
    TH1F* fTopEigen20Ratio;
    TH1F* fTopEigen10Ratio;
    TH1F* fTopPrimaryLength;

    TH1F* fSubNum3DHits;
    TH1F* fSubNumEdges;
    TH1F* fSubEigen21Ratio;
    TH1F* fSubEigen20Ratio;
    TH1F* fSubEigen10Ratio;
    TH1F* fSubPrimaryLength;
    TH1F* fSubCosToPrevPCA;
    TH1F* fSubCosExtToPCA;
    TH1F* fSubMaxDefect;
    TH1F* fSubUsedDefect;

    /**
     *  @brief Tools
     */
    std::unique_ptr<lar_cluster3d::IClusterAlg>
      fClusterAlg;                  ///<  Algorithm to do 3D space point clustering
    PrincipalComponentsAlg fPCAAlg; // For running Principal Components Analysis
  };

  VoronoiPathFinder::VoronoiPathFinder(fhicl::ParameterSet const& pset)
    : fPCAAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
  {
    this->configure(pset);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  VoronoiPathFinder::~VoronoiPathFinder() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void VoronoiPathFinder::configure(fhicl::ParameterSet const& pset)
  {
    fEnableMonitoring = pset.get<bool>("EnableMonitoring", true);
    fMinTinyClusterSize = pset.get<size_t>("MinTinyClusterSize", 40);
    fClusterAlg =
      art::make_tool<lar_cluster3d::IClusterAlg>(pset.get<fhicl::ParameterSet>("ClusterAlg"));

    fTimeToProcess = 0.;

    return;
  }

  void VoronoiPathFinder::initializeHistograms(art::TFileDirectory& histDir)
  {
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
    fFillHistograms = true;

    std::string dirName = "VoronoiPath";

    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());

    // Divide into two sets of hists... those for the overall cluster and
    // those for the subclusters
    fTopNum3DHits = dir.make<TH1F>("TopNum3DHits", "Number 3D Hits", 200, 0., 200.);
    fTopNumEdges = dir.make<TH1F>("TopNumEdges", "Number Edges", 200, 0., 200.);
    fTopEigen21Ratio = dir.make<TH1F>("TopEigen21Rat", "Eigen 2/1 Ratio", 100, 0., 1.);
    fTopEigen20Ratio = dir.make<TH1F>("TopEigen20Rat", "Eigen 2/0 Ratio", 100, 0., 1.);
    fTopEigen10Ratio = dir.make<TH1F>("TopEigen10Rat", "Eigen 1/0 Ratio", 100, 0., 1.);
    fTopPrimaryLength = dir.make<TH1F>("TopPrimaryLen", "Primary Length", 200, 0., 200.);

    fSubNum3DHits = dir.make<TH1F>("SubNum3DHits", "Number 3D Hits", 200, 0., 200.);
    fSubNumEdges = dir.make<TH1F>("SubNumEdges", "Number Edges", 200, 0., 200.);
    fSubEigen21Ratio = dir.make<TH1F>("SubEigen21Rat", "Eigen 2/1 Ratio", 100, 0., 1.);
    fSubEigen20Ratio = dir.make<TH1F>("SubEigen20Rat", "Eigen 2/0 Ratio", 100, 0., 1.);
    fSubEigen10Ratio = dir.make<TH1F>("SubEigen10Rat", "Eigen 1/0 Ratio", 100, 0., 1.);
    fSubPrimaryLength = dir.make<TH1F>("SubPrimaryLen", "Primary Length", 200, 0., 200.);
    fSubCosToPrevPCA = dir.make<TH1F>("SubCosToPrev", "Cos(theta)", 101, 0., 1.01);
    fSubCosExtToPCA = dir.make<TH1F>("SubCosExtPCA", "Cos(theta)", 102, -1.01, 1.01);
    fSubMaxDefect = dir.make<TH1F>("SubMaxDefect", "Max Defect", 100, 0., 50.);
    fSubUsedDefect = dir.make<TH1F>("SubUsedDefect", "Used Defect", 100, 0., 50.);

    return;
  }

  void VoronoiPathFinder::ModifyClusters(reco::ClusterParametersList& clusterParametersList) const
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

    int countClusters(0);

    // This is the loop over candidate 3D clusters
    // Note that it might be that the list of candidate clusters is modified by splitting
    // So we use the following construct to make sure we get all of them
    reco::ClusterParametersList::iterator clusterParametersListItr = clusterParametersList.begin();

    while (clusterParametersListItr != clusterParametersList.end()) {
      // Dereference to get the cluster paramters
      reco::ClusterParameters& clusterParameters = *clusterParametersListItr;

      std::cout << "**> Looking at Cluster " << countClusters++
                << ", # hits: " << clusterParameters.getHitPairListPtr().size() << std::endl;

      // It turns out that computing the convex hull surrounding the points in the 2D projection onto the
      // plane of largest spread in the PCA is a good way to break up the cluster... and we do it here since
      // we (currently) want this to be part of the standard output
      buildVoronoiDiagram(clusterParameters);

      // Make sure our cluster has enough hits...
      if (clusterParameters.getHitPairListPtr().size() > fMinTinyClusterSize) {
        // Get an interim cluster list
        reco::ClusterParametersList reclusteredParameters;

        // Call the main workhorse algorithm for building the local version of candidate 3D clusters
        //******** Remind me why we need to call this at this point when the same hits will be used? ********
        //fClusterAlg->Cluster3DHits(clusterParameters.getHitPairListPtr(), reclusteredParameters);
        reclusteredParameters.push_back(clusterParameters);

        std::cout << ">>>>>>>>>>> Reclustered to " << reclusteredParameters.size()
                  << " Clusters <<<<<<<<<<<<<<<" << std::endl;

        // Only process non-empty results
        if (!reclusteredParameters.empty()) {
          // Loop over the reclustered set
          for (auto& cluster : reclusteredParameters) {
            std::cout << "****> Calling breakIntoTinyBits with "
                      << cluster.getHitPairListPtr().size() << " hits" << std::endl;

            // It turns out that computing the convex hull surrounding the points in the 2D projection onto the
            // plane of largest spread in the PCA is a good way to break up the cluster... and we do it here since
            // we (currently) want this to be part of the standard output
            buildConvexHull(cluster, 2);

            // Break our cluster into smaller elements...
            subDivideCluster(cluster, cluster.daughterList().end(), cluster.daughterList(), 4);

            std::cout << "****> Broke Cluster with " << cluster.getHitPairListPtr().size()
                      << " into " << cluster.daughterList().size() << " sub clusters";
            for (auto& clus : cluster.daughterList())
              std::cout << ", " << clus.getHitPairListPtr().size();
            std::cout << std::endl;

            // Add the daughters to the cluster
            clusterParameters.daughterList().insert(clusterParameters.daughterList().end(),
                                                    cluster);

            // If filling histograms we do the main cluster here
            if (fFillHistograms) {
              reco::PrincipalComponents& fullPCA = cluster.getFullPCA();
              std::vector<double> eigenValVec = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
                                                 3. * std::sqrt(fullPCA.getEigenValues()[1]),
                                                 3. * std::sqrt(fullPCA.getEigenValues()[2])};
              double eigen2To1Ratio = eigenValVec[0] / eigenValVec[1];
              double eigen1To0Ratio = eigenValVec[1] / eigenValVec[2];
              double eigen2To0Ratio = eigenValVec[0] / eigenValVec[2];
              int num3DHits = cluster.getHitPairListPtr().size();
              int numEdges = cluster.getBestEdgeList().size();

              fTopNum3DHits->Fill(std::min(num3DHits, 199), 1.);
              fTopNumEdges->Fill(std::min(numEdges, 199), 1.);
              fTopEigen21Ratio->Fill(eigen2To1Ratio, 1.);
              fTopEigen20Ratio->Fill(eigen2To0Ratio, 1.);
              fTopEigen10Ratio->Fill(eigen1To0Ratio, 1.);
              fTopPrimaryLength->Fill(std::min(eigenValVec[0], 199.), 1.);
            }
          }
        }
      }

      // Go to next cluster parameters object
      clusterParametersListItr++;
    }

    if (fEnableMonitoring) {
      theClockBuildClusters.stop();

      fTimeToProcess = theClockBuildClusters.accumulated_real_time();
    }

    mf::LogDebug("Cluster3D") << ">>>>> Cluster Path finding done" << std::endl;

    return;
  }

  reco::ClusterParametersList::iterator VoronoiPathFinder::breakIntoTinyBits(
    reco::ClusterParameters& clusterToBreak,
    reco::PrincipalComponents& lastPCA,
    reco::ClusterParametersList::iterator positionItr,
    reco::ClusterParametersList& outputClusterList,
    int level) const
  {
    // This needs to be a recursive routine...
    // Idea is to take the input cluster and order 3D hits by arclength along PCA primary axis
    // If the cluster is above the minimum number of hits then we divide into two and call ourself
    // with the two halves. This means we form a new cluster with hits and PCA and then call ourself
    // If the cluster is below the minimum then we can't break any more, simply add this cluster to
    // the new output list.

    // set an indention string
    std::string pluses(level / 2, '+');
    std::string indent(level / 2, ' ');

    indent += pluses;

    reco::ClusterParametersList::iterator inputPositionItr = positionItr;

    // To make best use of this we'll also want the PCA for this cluster... so...
    // Recover the prime ingredients
    reco::PrincipalComponents& fullPCA = clusterToBreak.getFullPCA();
    std::vector<double> eigenValVec = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
                                       3. * std::sqrt(fullPCA.getEigenValues()[1]),
                                       3. * std::sqrt(fullPCA.getEigenValues()[2])};
    Eigen::Vector3f fullPrimaryVec(fullPCA.getEigenVectors().row(2));
    Eigen::Vector3f lastPrimaryVec(lastPCA.getEigenVectors().row(2));

    double cosNewToLast = std::abs(fullPrimaryVec.dot(lastPrimaryVec));
    double eigen2To1Ratio = eigenValVec[2] / eigenValVec[1];
    double eigen1To0Ratio = eigenValVec[1] / eigenValVec[2];
    double eigen2To0Ratio = eigenValVec[0] / eigenValVec[2];
    double eigen2And1Ave = 0.5 * (eigenValVec[1] + eigenValVec[0]);
    double eigenAveTo0Ratio = eigen2And1Ave / eigenValVec[2];

    bool storeCurrentCluster(true);
    int minimumClusterSize(fMinTinyClusterSize);

    std::cout << indent << ">>> breakIntoTinyBits with "
              << clusterToBreak.getHitPairListPtr().size() << " input hits, "
              << clusterToBreak.getBestEdgeList().size() << " edges, rat21: " << eigen2To1Ratio
              << ", rat20: " << eigen2To0Ratio << ", rat10: " << eigen1To0Ratio
              << ", ave0: " << eigenAveTo0Ratio << std::endl;
    std::cout << indent << "   --> eigen 0/1/2: " << eigenValVec[0] << "/" << eigenValVec[1] << "/"
              << eigenValVec[2] << ", cos: " << cosNewToLast << std::endl;

    // Create a rough cut intended to tell us when we've reached the land of diminishing returns
    if (clusterToBreak.getBestEdgeList().size() > 5 && cosNewToLast > 0.25 &&
        eigen2To1Ratio < 0.9 && eigen2To0Ratio > 0.001 &&
        clusterToBreak.getHitPairListPtr().size() > size_t(2 * minimumClusterSize)) {
      // We want to find the convex hull vertices that lie furthest from the line to/from the extreme points
      // To find these we:
      // 1) recover the extreme points
      // 2) form the vector between them
      // 3) loop through the vertices and keep track of distance to this vector
      // 4) Sort the resulting list by furthest points and select the one we want
      reco::ProjectedPointList::const_iterator extremePointListItr =
        clusterToBreak.getConvexHull().getConvexHullExtremePoints().begin();

      const reco::ClusterHit3D* firstEdgeHit = std::get<2>(*extremePointListItr++);
      const reco::ClusterHit3D* secondEdgeHit = std::get<2>(*extremePointListItr);
      Eigen::Vector3f edgeVec(secondEdgeHit->getPosition()[0] - firstEdgeHit->getPosition()[0],
                              secondEdgeHit->getPosition()[1] - firstEdgeHit->getPosition()[1],
                              secondEdgeHit->getPosition()[2] - firstEdgeHit->getPosition()[2]);
      double edgeLen = edgeVec.norm();

      // normalize it
      edgeVec.normalize();

      // Recover the list of 3D hits associated to this cluster
      reco::HitPairListPtr& clusHitPairVector = clusterToBreak.getHitPairListPtr();

      // Calculate the doca to the PCA primary axis for each 3D hit
      // Importantly, this also gives us the arclength along that axis to the hit
      fPCAAlg.PCAAnalysis_calc3DDocas(clusHitPairVector, fullPCA);

      // Sort the hits along the PCA
      clusHitPairVector.sort([](const auto& left, const auto& right) {
        return left->getArclenToPoca() < right->getArclenToPoca();
      });

      // Set up container to keep track of edges
      using DistEdgeTuple = std::tuple<float, const reco::EdgeTuple*>;
      using DistEdgeTupleVec = std::vector<DistEdgeTuple>;

      DistEdgeTupleVec distEdgeTupleVec;

      // Now loop through all the edges and search for the furthers point
      for (const auto& edge : clusterToBreak.getBestEdgeList()) {
        const reco::ClusterHit3D* nextEdgeHit = std::get<0>(edge); // recover the first point

        // Create vector to this point from the longest edge
        Eigen::Vector3f hitToEdgeVec(nextEdgeHit->getPosition()[0] - firstEdgeHit->getPosition()[0],
                                     nextEdgeHit->getPosition()[1] - firstEdgeHit->getPosition()[1],
                                     nextEdgeHit->getPosition()[2] -
                                       firstEdgeHit->getPosition()[2]);

        // Get projection
        float hitProjection = hitToEdgeVec.dot(edgeVec);

        // Require that the point is really "opposite" the longest edge
        if (hitProjection > 0. && hitProjection < edgeLen) {
          Eigen::Vector3f distToHitVec = hitToEdgeVec - hitProjection * edgeVec;
          float distToHit = distToHitVec.norm();

          distEdgeTupleVec.emplace_back(distToHit, &edge);
        }
      }

      std::sort(
        distEdgeTupleVec.begin(), distEdgeTupleVec.end(), [](const auto& left, const auto& right) {
          return std::get<0>(left) > std::get<0>(right);
        });

      for (const auto& distEdgeTuple : distEdgeTupleVec) {
        const reco::EdgeTuple& edgeTuple = *std::get<1>(distEdgeTuple);
        const reco::ClusterHit3D* edgeHit = std::get<0>(edgeTuple);

        // Now find the hit identified above as furthest away
        reco::HitPairListPtr::iterator vertexItr =
          std::find(clusHitPairVector.begin(), clusHitPairVector.end(), edgeHit);

        // Make sure enough hits either side, otherwise we just keep the current cluster
        if (vertexItr == clusHitPairVector.end() ||
            std::distance(clusHitPairVector.begin(), vertexItr) < minimumClusterSize ||
            std::distance(vertexItr, clusHitPairVector.end()) < minimumClusterSize)
          continue;

        // Now we create a list of pairs of iterators to the start and end of each subcluster
        using Hit3DItrPair =
          std::pair<reco::HitPairListPtr::iterator, reco::HitPairListPtr::iterator>;
        using VertexPairList = std::list<Hit3DItrPair>;

        VertexPairList vertexPairList;

        vertexPairList.emplace_back(Hit3DItrPair(clusHitPairVector.begin(), vertexItr));
        vertexPairList.emplace_back(Hit3DItrPair(vertexItr, clusHitPairVector.end()));

        storeCurrentCluster = false;

        // Ok, now loop through our pairs
        for (auto& hit3DItrPair : vertexPairList) {
          reco::ClusterParameters clusterParams;
          reco::HitPairListPtr& hitPairListPtr = clusterParams.getHitPairListPtr();

          std::cout << indent << "+>    -- building new cluster, size: "
                    << std::distance(hit3DItrPair.first, hit3DItrPair.second) << std::endl;

          // size the container...
          hitPairListPtr.resize(std::distance(hit3DItrPair.first, hit3DItrPair.second));

          // and copy the hits into the container
          std::copy(hit3DItrPair.first, hit3DItrPair.second, hitPairListPtr.begin());

          // First stage of feature extraction runs here
          fPCAAlg.PCAAnalysis_3D(hitPairListPtr, clusterParams.getFullPCA());

          // Recover the new fullPCA
          reco::PrincipalComponents& newFullPCA = clusterParams.getFullPCA();

          // Must have a valid pca
          if (newFullPCA.getSvdOK()) {
            std::cout << indent << "+>    -- >> cluster has a valid Full PCA" << std::endl;

            // If the PCA's are opposite the flip the axes
            if (fullPrimaryVec.dot(newFullPCA.getEigenVectors().row(2)) < 0.) {
              for (size_t vecIdx = 0; vecIdx < 3; vecIdx++)
                newFullPCA.flipAxis(vecIdx);
            }

            // Set the skeleton PCA to make sure it has some value
            clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();

            // Be sure to compute the oonvex hull surrounding the now broken cluster
            buildConvexHull(clusterParams, level + 2);

            positionItr =
              breakIntoTinyBits(clusterParams, fullPCA, positionItr, outputClusterList, level + 4);
          }
        }

        // If successful in breaking the cluster then we are done, otherwise try to loop
        // again taking the next furthest hit
        break;
      }
    }

    // First question, are we done breaking?
    if (storeCurrentCluster) {
      // I think this is where we fill out the rest of the parameters?
      // Start by adding the 2D hits...
      // See if we can avoid duplicates by temporarily transferring to a set
      std::set<const reco::ClusterHit2D*> hitSet;

      // Loop through 3D hits to get a set of unique 2D hits
      for (const auto& hit3D : clusterToBreak.getHitPairListPtr()) {
        for (const auto& hit2D : hit3D->getHits()) {
          if (hit2D) hitSet.insert(hit2D);
        }
      }

      // Now add these to the new cluster
      for (const auto& hit2D : hitSet) {
        hit2D->setStatusBit(reco::ClusterHit2D::USED);
        clusterToBreak.UpdateParameters(hit2D);
      }

      std::cout << indent << "*********>>> storing new subcluster of size "
                << clusterToBreak.getHitPairListPtr().size() << std::endl;

      positionItr = outputClusterList.insert(positionItr, clusterToBreak);

      // Are we filling histograms
      if (fFillHistograms) {
        int num3DHits = clusterToBreak.getHitPairListPtr().size();
        int numEdges = clusterToBreak.getBestEdgeList().size();
        Eigen::Vector3f newPrimaryVec(fullPCA.getEigenVectors().row(2));
        Eigen::Vector3f lastPrimaryVec(lastPCA.getEigenVectors().row(2));
        float cosToLast = newPrimaryVec.dot(lastPrimaryVec);

        fSubNum3DHits->Fill(std::min(num3DHits, 199), 1.);
        fSubNumEdges->Fill(std::min(numEdges, 199), 1.);
        fSubEigen21Ratio->Fill(eigen2To1Ratio, 1.);
        fSubEigen20Ratio->Fill(eigen2To0Ratio, 1.);
        fSubEigen10Ratio->Fill(eigen1To0Ratio, 1.);
        fSubCosToPrevPCA->Fill(cosToLast, 1.);
        fSubPrimaryLength->Fill(std::min(eigenValVec[2], 199.), 1.);
      }

      // The above points to the element, want the next element
      positionItr++;
    }
    else if (inputPositionItr != positionItr) {
      std::cout << indent << "***** DID NOT STORE A CLUSTER *****" << std::endl;
    }

    return positionItr;
  }

  reco::ClusterParametersList::iterator VoronoiPathFinder::subDivideCluster(
    reco::ClusterParameters& clusterToBreak,
    reco::ClusterParametersList::iterator positionItr,
    reco::ClusterParametersList& outputClusterList,
    int level) const
  {
    // This is a recursive routine to divide an input cluster, according to the maximum defect point of
    // the convex hull until we reach the point of no further improvement.
    // The assumption is that the input cluster is fully formed already, this routine then simply
    // divides, if successful division into two pieces it then stores the results

    // No point in doing anything if we don't have enough space points
    if (clusterToBreak.getHitPairListPtr().size() > size_t(2 * fMinTinyClusterSize)) {
      // set an indention string
      std::string pluses(level / 2, '+');
      std::string indent(level / 2, ' ');

      indent += pluses;

      // To make best use of this we'll also want the PCA for this cluster... so...
      // Recover the prime ingredients
      reco::PrincipalComponents& fullPCA = clusterToBreak.getFullPCA();
      std::vector<double> eigenValVec = {3. * std::sqrt(fullPCA.getEigenValues()[0]),
                                         3. * std::sqrt(fullPCA.getEigenValues()[1]),
                                         3. * std::sqrt(fullPCA.getEigenValues()[2])};
      Eigen::Vector3f fullPrimaryVec(fullPCA.getEigenVectors().row(2));

      // We want to find the convex hull vertices that lie furthest from the line to/from the extreme points
      // To find these we:
      // 1) recover the extreme points
      // 2) form the vector between them
      // 3) loop through the vertices and keep track of distance to this vector
      // 4) Sort the resulting list by furthest points and select the one we want
      reco::ProjectedPointList::const_iterator extremePointListItr =
        clusterToBreak.getConvexHull().getConvexHullExtremePoints().begin();

      const reco::ClusterHit3D* firstEdgeHit = std::get<2>(*extremePointListItr++);
      const reco::ClusterHit3D* secondEdgeHit = std::get<2>(*extremePointListItr);
      Eigen::Vector3f edgeVec(secondEdgeHit->getPosition()[0] - firstEdgeHit->getPosition()[0],
                              secondEdgeHit->getPosition()[1] - firstEdgeHit->getPosition()[1],
                              secondEdgeHit->getPosition()[2] - firstEdgeHit->getPosition()[2]);
      double edgeLen = edgeVec.norm();

      // normalize it
      edgeVec.normalize();

      // Recover the list of 3D hits associated to this cluster
      reco::HitPairListPtr& clusHitPairVector = clusterToBreak.getHitPairListPtr();

      // Calculate the doca to the PCA primary axis for each 3D hit
      // Importantly, this also gives us the arclength along that axis to the hit
      fPCAAlg.PCAAnalysis_calc3DDocas(clusHitPairVector, fullPCA);

      // Sort the hits along the PCA
      clusHitPairVector.sort([](const auto& left, const auto& right) {
        return left->getArclenToPoca() < right->getArclenToPoca();
      });

      // Set up container to keep track of edges
      using DistEdgeTuple = std::tuple<float, const reco::EdgeTuple*>;
      using DistEdgeTupleVec = std::vector<DistEdgeTuple>;

      DistEdgeTupleVec distEdgeTupleVec;

      // Now loop through all the edges and search for the furthers point
      for (const auto& edge : clusterToBreak.getBestEdgeList()) {
        const reco::ClusterHit3D* nextEdgeHit = std::get<0>(edge); // recover the first point

        // Create vector to this point from the longest edge
        Eigen::Vector3f hitToEdgeVec(nextEdgeHit->getPosition()[0] - firstEdgeHit->getPosition()[0],
                                     nextEdgeHit->getPosition()[1] - firstEdgeHit->getPosition()[1],
                                     nextEdgeHit->getPosition()[2] -
                                       firstEdgeHit->getPosition()[2]);

        // Get projection
        float hitProjection = hitToEdgeVec.dot(edgeVec);

        // Require that the point is really "opposite" the longest edge
        if (hitProjection > 0. && hitProjection < edgeLen) {
          Eigen::Vector3f distToHitVec = hitToEdgeVec - hitProjection * edgeVec;
          float distToHit = distToHitVec.norm();

          distEdgeTupleVec.emplace_back(distToHit, &edge);
        }
      }

      std::sort(
        distEdgeTupleVec.begin(), distEdgeTupleVec.end(), [](const auto& left, const auto& right) {
          return std::get<0>(left) > std::get<0>(right);
        });

      // Get a temporary  container to hol
      reco::ClusterParametersList tempClusterParametersList;
      float usedDefectDist(0.);

      for (const auto& distEdgeTuple : distEdgeTupleVec) {
        const reco::EdgeTuple& edgeTuple = *std::get<1>(distEdgeTuple);
        const reco::ClusterHit3D* edgeHit = std::get<0>(edgeTuple);

        usedDefectDist = std::get<0>(distEdgeTuple);

        // Now find the hit identified above as furthest away
        reco::HitPairListPtr::iterator vertexItr =
          std::find(clusHitPairVector.begin(), clusHitPairVector.end(), edgeHit);

        // Make sure enough hits either side, otherwise we just keep the current cluster
        if (vertexItr == clusHitPairVector.end() ||
            std::distance(clusHitPairVector.begin(), vertexItr) < int(fMinTinyClusterSize) ||
            std::distance(vertexItr, clusHitPairVector.end()) < int(fMinTinyClusterSize))
          continue;

        tempClusterParametersList.push_back(reco::ClusterParameters());

        reco::ClusterParameters& clusterParams1 = tempClusterParametersList.back();

        if (makeCandidateCluster(
              fullPrimaryVec, clusterParams1, clusHitPairVector.begin(), vertexItr, level)) {
          tempClusterParametersList.push_back(reco::ClusterParameters());

          reco::ClusterParameters& clusterParams2 = tempClusterParametersList.back();

          if (makeCandidateCluster(
                fullPrimaryVec, clusterParams2, vertexItr, clusHitPairVector.end(), level))
            break;
        }

        // If here then we could not make two valid clusters and so we try again
        tempClusterParametersList.clear();
      }

      // Fallback in the event of still large clusters but not defect points
      if (tempClusterParametersList.empty()) {
        std::cout << indent << "===> no cluster cands, edgeLen: " << edgeLen
                  << ", # hits: " << clusHitPairVector.size()
                  << ", max defect: " << std::get<0>(distEdgeTupleVec.front()) << std::endl;

        usedDefectDist = 0.;

        if (edgeLen > 20.) {
          reco::HitPairListPtr::iterator vertexItr = clusHitPairVector.begin();

          std::advance(vertexItr, clusHitPairVector.size() / 2);

          tempClusterParametersList.push_back(reco::ClusterParameters());

          reco::ClusterParameters& clusterParams1 = tempClusterParametersList.back();

          if (makeCandidateCluster(
                fullPrimaryVec, clusterParams1, clusHitPairVector.begin(), vertexItr, level)) {
            tempClusterParametersList.push_back(reco::ClusterParameters());

            reco::ClusterParameters& clusterParams2 = tempClusterParametersList.back();

            if (!makeCandidateCluster(
                  fullPrimaryVec, clusterParams2, vertexItr, clusHitPairVector.end(), level))
              tempClusterParametersList.clear();
          }
        }
      }

      // Can only end with no candidate clusters or two so don't
      for (auto& clusterParams : tempClusterParametersList) {
        size_t curOutputClusterListSize = outputClusterList.size();

        positionItr = subDivideCluster(clusterParams, positionItr, outputClusterList, level + 4);

        std::cout << indent << "Output cluster list prev: " << curOutputClusterListSize
                  << ", now: " << outputClusterList.size() << std::endl;

        // If the cluster we sent in was successfully broken then the position iterator will be shifted
        // This means we don't want to restore the current cluster here
        if (curOutputClusterListSize < outputClusterList.size()) continue;

        // I think this is where we fill out the rest of the parameters?
        // Start by adding the 2D hits...
        // See if we can avoid duplicates by temporarily transferring to a set
        std::set<const reco::ClusterHit2D*> hitSet;

        // Loop through 3D hits to get a set of unique 2D hits
        for (const auto& hit3D : clusterParams.getHitPairListPtr()) {
          for (const auto& hit2D : hit3D->getHits()) {
            if (hit2D) hitSet.insert(hit2D);
          }
        }

        // Now add these to the new cluster
        for (const auto& hit2D : hitSet) {
          hit2D->setStatusBit(reco::ClusterHit2D::USED);
          clusterParams.UpdateParameters(hit2D);
        }

        std::cout << indent << "*********>>> storing new subcluster of size "
                  << clusterParams.getHitPairListPtr().size() << std::endl;

        positionItr = outputClusterList.insert(positionItr, clusterParams);

        // Are we filling histograms
        if (fFillHistograms) {
          // Recover the new fullPCA
          reco::PrincipalComponents& newFullPCA = clusterParams.getFullPCA();

          Eigen::Vector3f newPrimaryVec(fullPCA.getEigenVectors().row(2));
          Eigen::Vector3f lastPrimaryVec(newFullPCA.getEigenVectors().row(2));

          int num3DHits = clusterParams.getHitPairListPtr().size();
          int numEdges = clusterParams.getBestEdgeList().size();
          float cosToLast = newPrimaryVec.dot(lastPrimaryVec);
          double eigen2To1Ratio = eigenValVec[0] / eigenValVec[1];
          double eigen1To0Ratio = eigenValVec[1] / eigenValVec[2];
          double eigen2To0Ratio = eigenValVec[2] / eigenValVec[2];

          fSubNum3DHits->Fill(std::min(num3DHits, 199), 1.);
          fSubNumEdges->Fill(std::min(numEdges, 199), 1.);
          fSubEigen21Ratio->Fill(eigen2To1Ratio, 1.);
          fSubEigen20Ratio->Fill(eigen2To0Ratio, 1.);
          fSubEigen10Ratio->Fill(eigen1To0Ratio, 1.);
          fSubCosToPrevPCA->Fill(cosToLast, 1.);
          fSubPrimaryLength->Fill(std::min(eigenValVec[2], 199.), 1.);
          fSubCosExtToPCA->Fill(fullPrimaryVec.dot(edgeVec), 1.);
          fSubMaxDefect->Fill(std::get<0>(distEdgeTupleVec.front()), 1.);
          fSubUsedDefect->Fill(usedDefectDist, 1.);
        }

        // The above points to the element, want the next element
        positionItr++;
      }
    }

    return positionItr;
  }

  bool VoronoiPathFinder::makeCandidateCluster(Eigen::Vector3f& primaryPCA,
                                               reco::ClusterParameters& candCluster,
                                               reco::HitPairListPtr::iterator firstHitItr,
                                               reco::HitPairListPtr::iterator lastHitItr,
                                               int level) const
  {
    std::string indent(level / 2, ' ');

    reco::HitPairListPtr& hitPairListPtr = candCluster.getHitPairListPtr();

    std::cout << indent
              << "+>    -- building new cluster, size: " << std::distance(firstHitItr, lastHitItr)
              << std::endl;

    // size the container...
    hitPairListPtr.resize(std::distance(firstHitItr, lastHitItr));

    // and copy the hits into the container
    std::copy(firstHitItr, lastHitItr, hitPairListPtr.begin());

    // First stage of feature extraction runs here
    fPCAAlg.PCAAnalysis_3D(hitPairListPtr, candCluster.getFullPCA());

    // Recover the new fullPCA
    reco::PrincipalComponents& newFullPCA = candCluster.getFullPCA();

    // Will we want to store this cluster?
    bool keepThisCluster(false);

    // Must have a valid pca
    if (newFullPCA.getSvdOK()) {
      std::cout << indent << "+>    -- >> cluster has a valid Full PCA" << std::endl;

      // Need to check if the PCA direction has been reversed
      Eigen::Vector3f newPrimaryVec(newFullPCA.getEigenVectors().row(2));

      // If the PCA's are opposite the flip the axes
      if (primaryPCA.dot(newPrimaryVec) < 0.) {
        for (size_t vecIdx = 0; vecIdx < 3; vecIdx++)
          newFullPCA.flipAxis(vecIdx);

        newPrimaryVec = Eigen::Vector3f(newFullPCA.getEigenVectors().row(2));
      }

      // Set the skeleton PCA to make sure it has some value
      candCluster.getSkeletonPCA() = candCluster.getFullPCA();

      // Be sure to compute the oonvex hull surrounding the now broken cluster
      buildConvexHull(candCluster, level + 2);

      std::vector<double> eigenValVec = {3. * std::sqrt(newFullPCA.getEigenValues()[0]),
                                         3. * std::sqrt(newFullPCA.getEigenValues()[1]),
                                         3. * std::sqrt(newFullPCA.getEigenValues()[2])};
      double cosNewToLast = std::abs(primaryPCA.dot(newPrimaryVec));
      double eigen2To1Ratio = eigenValVec[2] / eigenValVec[1];
      double eigen1To0Ratio = eigenValVec[1] / eigenValVec[2];
      double eigen2To0Ratio = eigenValVec[0] / eigenValVec[2];
      double eigen2And1Ave = 0.5 * (eigenValVec[1] + eigenValVec[0]);
      double eigenAveTo0Ratio = eigen2And1Ave / eigenValVec[2];

      std::cout << indent << ">>> subDivideClusters with " << candCluster.getHitPairListPtr().size()
                << " input hits, " << candCluster.getBestEdgeList().size()
                << " edges, rat21: " << eigen2To1Ratio << ", rat20: " << eigen2To0Ratio
                << ", rat10: " << eigen1To0Ratio << ", ave0: " << eigenAveTo0Ratio << std::endl;
      std::cout << indent << "   --> eigen 0/1/2: " << eigenValVec[0] << "/" << eigenValVec[1]
                << "/" << eigenValVec[2] << ", cos: " << cosNewToLast << std::endl;

      // Create a rough cut intended to tell us when we've reached the land of diminishing returns
      //        if (candCluster.getBestEdgeList().size() > 4 && cosNewToLast > 0.25 && eigen2To1Ratio < 0.9 && eigen2To0Ratio > 0.001)
      if (candCluster.getBestEdgeList().size() > 4 && cosNewToLast > 0.25 &&
          eigen2To1Ratio > 0.01 && eigen2To1Ratio < 0.99 && eigen1To0Ratio < 0.5) {
        keepThisCluster = true;
      }
    }

    return keepThisCluster;
  }

  void VoronoiPathFinder::buildConvexHull(reco::ClusterParameters& clusterParameters,
                                          int level) const
  {
    // set an indention string
    std::string minuses(level / 2, '-');
    std::string indent(level / 2, ' ');

    indent += minuses;

    // The plan is to build the enclosing 2D polygon around the points in the PCA plane of most spread for this cluster
    // To do so we need to start by building a list of 2D projections onto the plane of most spread...
    reco::PrincipalComponents& pca = clusterParameters.getFullPCA();

    // Recover the parameters from the Principal Components Analysis that we need to project and accumulate
    Eigen::Vector3f pcaCenter(
      pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);

    //dcel2d::PointList pointList;
    using Point = std::tuple<float, float, const reco::ClusterHit3D*>;
    using PointList = std::list<Point>;

    reco::ConvexHull& convexHull = clusterParameters.getConvexHull();
    PointList pointList = convexHull.getProjectedPointList();

    // Loop through hits and do projection to plane
    for (const auto& hit3D : clusterParameters.getHitPairListPtr()) {
      Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
                                  hit3D->getPosition()[1] - pcaCenter(1),
                                  hit3D->getPosition()[2] - pcaCenter(2));
      Eigen::Vector3f pcaToHit = pca.getEigenVectors() * pcaToHitVec;

      pointList.emplace_back(dcel2d::Point(pcaToHit(0), pcaToHit(1), hit3D));
    }

    // Sort the point vec by increasing x, then increase y
    pointList.sort([](const auto& left, const auto& right) {
      return (std::abs(std::get<0>(left) - std::get<0>(right)) >
              std::numeric_limits<float>::epsilon()) ?
               std::get<0>(left) < std::get<0>(right) :
               std::get<1>(left) < std::get<1>(right);
    });

    // containers for finding the "best" hull...
    std::vector<ConvexHull> convexHullVec;
    std::vector<PointList> rejectedListVec;
    bool increaseDepth(pointList.size() > 5);
    float lastArea(std::numeric_limits<float>::max());

    while (increaseDepth) {
      // Get another convexHull container
      convexHullVec.push_back(ConvexHull(pointList));
      rejectedListVec.push_back(PointList());

      const ConvexHull& convexHull = convexHullVec.back();
      PointList& rejectedList = rejectedListVec.back();
      const PointList& convexHullPoints = convexHull.getConvexHull();

      increaseDepth = false;

      if (convexHull.getConvexHullArea() > 0.) {
        std::cout << indent << "-> built convex hull, 3D hits: " << pointList.size() << " with "
                  << convexHullPoints.size() << " vertices"
                  << ", area: " << convexHull.getConvexHullArea() << std::endl;
        std::cout << indent << "-> -Points:";
        for (const auto& point : convexHullPoints)
          std::cout << " (" << std::get<0>(point) << "," << std::get<1>(point) << ")";
        std::cout << std::endl;

        if (convexHullVec.size() < 2 || convexHull.getConvexHullArea() < 0.8 * lastArea) {
          for (auto& point : convexHullPoints) {
            pointList.remove(point);
            rejectedList.emplace_back(point);
          }
          lastArea = convexHull.getConvexHullArea();
          //                increaseDepth = true;
        }
      }
    }

    // do we have a valid convexHull?
    while (!convexHullVec.empty() && convexHullVec.back().getConvexHullArea() < 0.5) {
      convexHullVec.pop_back();
      rejectedListVec.pop_back();
    }

    // If we found the convex hull then build edges around the region
    if (!convexHullVec.empty()) {
      size_t nRejectedTotal(0);
      reco::HitPairListPtr hitPairListPtr = clusterParameters.getHitPairListPtr();

      for (const auto& rejectedList : rejectedListVec) {
        nRejectedTotal += rejectedList.size();

        for (const auto& rejectedPoint : rejectedList) {
          std::cout << indent << "-> -- Point is "
                    << convexHullVec.back().findNearestDistance(rejectedPoint)
                    << " from nearest edge" << std::endl;

          if (convexHullVec.back().findNearestDistance(rejectedPoint) > 0.5)
            hitPairListPtr.remove(std::get<2>(rejectedPoint));
        }
      }

      std::cout << indent << "-> Removed " << nRejectedTotal << " leaving " << pointList.size()
                << "/" << hitPairListPtr.size() << " points" << std::endl;

      // Now add "edges" to the cluster to describe the convex hull (for the display)
      reco::Hit3DToEdgeMap& edgeMap = convexHull.getConvexHullEdgeMap();
      reco::EdgeList& bestEdgeList = convexHull.getConvexHullEdgeList();

      Point lastPoint = convexHullVec.back().getConvexHull().front();

      for (auto& curPoint : convexHullVec.back().getConvexHull()) {
        if (curPoint == lastPoint) continue;

        const reco::ClusterHit3D* lastPoint3D = std::get<2>(lastPoint);
        const reco::ClusterHit3D* curPoint3D = std::get<2>(curPoint);

        float distBetweenPoints = (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) *
                                    (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) +
                                  (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) *
                                    (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) +
                                  (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]) *
                                    (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]);

        distBetweenPoints = std::sqrt(distBetweenPoints);

        reco::EdgeTuple edge(lastPoint3D, curPoint3D, distBetweenPoints);

        edgeMap[lastPoint3D].push_back(edge);
        edgeMap[curPoint3D].push_back(edge);
        bestEdgeList.emplace_back(edge);

        lastPoint = curPoint;
      }

      // Store the "extreme" points
      const ConvexHull::PointList& extremePoints = convexHullVec.back().getExtremePoints();
      reco::ProjectedPointList& extremePointList = convexHull.getConvexHullExtremePoints();

      for (const auto& point : extremePoints)
        extremePointList.push_back(point);
    }

    return;
  }

  void VoronoiPathFinder::buildVoronoiDiagram(reco::ClusterParameters& clusterParameters) const
  {
    // The plan is to build the enclosing 2D polygon around the points in the PCA plane of most spread for this cluster
    // To do so we need to start by building a list of 2D projections onto the plane of most spread...
    reco::PrincipalComponents& pca = clusterParameters.getFullPCA();

    // Recover the parameters from the Principal Components Analysis that we need to project and accumulate
    Eigen::Vector3f pcaCenter(
      pca.getAvePosition()[0], pca.getAvePosition()[1], pca.getAvePosition()[2]);

    dcel2d::PointList pointList;

    // Loop through hits and do projection to plane
    for (const auto& hit3D : clusterParameters.getHitPairListPtr()) {
      Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
                                  hit3D->getPosition()[1] - pcaCenter(1),
                                  hit3D->getPosition()[2] - pcaCenter(2));
      Eigen::Vector3f pcaToHit = pca.getEigenVectors() * pcaToHitVec;

      pointList.emplace_back(dcel2d::Point(pcaToHit(1), pcaToHit(2), hit3D));
    }

    // Sort the point vec by increasing x, then increase y
    pointList.sort([](const auto& left, const auto& right) {
      return (std::abs(std::get<0>(left) - std::get<0>(right)) >
              std::numeric_limits<float>::epsilon()) ?
               std::get<0>(left) < std::get<0>(right) :
               std::get<1>(left) < std::get<1>(right);
    });

    std::cout << "  ==> Build V diagram, sorted point list contains " << pointList.size() << " hits"
              << std::endl;

    // Set up the voronoi diagram builder
    voronoi2d::VoronoiDiagram voronoiDiagram(clusterParameters.getHalfEdgeList(),
                                             clusterParameters.getVertexList(),
                                             clusterParameters.getFaceList());

    // And make the diagram
    voronoiDiagram.buildVoronoiDiagram(pointList);

    // Recover the voronoi diagram vertex list and the container to store them in
    dcel2d::VertexList& vertexList = clusterParameters.getVertexList();

    // Now get the inverse of the rotation matrix so we can get the vertex positions,
    // which lie in the plane of the two largest PCA axes, in the standard coordinate system
    Eigen::Matrix3f rotationMatrixInv = pca.getEigenVectors().inverse();

    // Translate and fill
    for (auto& vertex : vertexList) {
      Eigen::Vector3f coords = rotationMatrixInv * vertex.getCoords();

      coords += pcaCenter;

      vertex.setCoords(coords);
    }

    // Now do the Convex Hull
    // Now add "edges" to the cluster to describe the convex hull (for the display)
    reco::Hit3DToEdgeMap& edgeMap = clusterParameters.getHit3DToEdgeMap();
    reco::EdgeList& bestEdgeList = clusterParameters.getBestEdgeList();

    //    const dcel2d::PointList& edgePoints = voronoiDiagram.getConvexHull();
    //    PointList                localList;
    //
    //    for(const auto& edgePoint : edgePoints) localList.emplace_back(std::get<0>(edgePoint),std::get<1>(edgePoint),std::get<2>(edgePoint));
    //
    //    // Sort the point vec by increasing x, then increase y
    //    localList.sort([](const auto& left, const auto& right){return (std::abs(std::get<0>(left) - std::get<0>(right)) > std::numeric_limits<float>::epsilon()) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);});
    //
    //    // Why do I need to do this?
    //    ConvexHull convexHull(localList);
    //
    //    Point lastPoint = convexHull.getConvexHull().front();
    dcel2d::Point lastPoint = voronoiDiagram.getConvexHull().front();

    //    std::cout << "@@@@>> Build convex hull, voronoi handed " << edgePoints.size() << " points, convexHull cut to " << convexHull.getConvexHull().size() << std::endl;

    //    for(auto& curPoint : convexHull.getConvexHull())
    for (auto& curPoint : voronoiDiagram.getConvexHull()) {
      if (curPoint == lastPoint) continue;

      const reco::ClusterHit3D* lastPoint3D = std::get<2>(lastPoint);
      const reco::ClusterHit3D* curPoint3D = std::get<2>(curPoint);

      float distBetweenPoints = (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) *
                                  (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) +
                                (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) *
                                  (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) +
                                (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]) *
                                  (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]);

      distBetweenPoints = std::sqrt(distBetweenPoints);

      reco::EdgeTuple edge(lastPoint3D, curPoint3D, distBetweenPoints);

      edgeMap[lastPoint3D].push_back(edge);
      edgeMap[curPoint3D].push_back(edge);
      bestEdgeList.emplace_back(edge);

      lastPoint = curPoint;
    }

    std::cout << "****> vertexList containted " << vertexList.size() << " vertices for "
              << clusterParameters.getHitPairListPtr().size() << " hits" << std::endl;

    return;
  }

  float VoronoiPathFinder::closestApproach(const Eigen::Vector3f& P0,
                                           const Eigen::Vector3f& u0,
                                           const Eigen::Vector3f& P1,
                                           const Eigen::Vector3f& u1,
                                           Eigen::Vector3f& poca0,
                                           Eigen::Vector3f& poca1) const
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

  DEFINE_ART_CLASS_TOOL(VoronoiPathFinder)
} // namespace lar_cluster3d

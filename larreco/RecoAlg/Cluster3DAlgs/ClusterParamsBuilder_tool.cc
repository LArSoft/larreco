/**
 *  @file   ClusterParamsBuilder.cc
 *
 *  @brief  A tool to select good clusters and fill their output parameters
 *
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterParamsBuilder.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"

// std includes
#include <iostream>
#include <memory>
#include <unordered_set>

namespace lar_cluster3d {

  /**
 *  @brief  ClusterParamsBuilder class definiton
 */
  class ClusterParamsBuilder : virtual public IClusterParametersBuilder {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit ClusterParamsBuilder(fhicl::ParameterSet const& pset);

    /**
     *  @brief  Destructor
     */
    virtual ~ClusterParamsBuilder();

    void configure(const fhicl::ParameterSet&) override;

    /**
     *  @brief Given the results of running DBScan, format the clusters so that they can be
     *         easily transferred back to the larsoft world
     *
     *  @param hitPairClusterMap      map between view and a list of 3D hits
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     *
     *                                The last two parameters are passed through to the FillClusterParams method
     */
    void BuildClusterInfo(reco::ClusterParametersList& clusterParametersList) const override;

    /**
     *  @brief Fill the cluster parameters (expose to outside world for case of splitting/merging clusters)
     *
     *  @param ClusterParameters   The cluster parameters container to be modified
     *  @param Hit2DToClusterMap   Map to keep track of 2D hit to cluster association
     *  @param double              minimum fraction of unique hits
     *  @param double              maximum fraction of "lost" hits
     */

    void FillClusterParams(reco::ClusterParameters&,
                           reco::Hit2DToClusterMap&,
                           double minUniqueFrac = 0.,
                           double maxLostFrac = 1.) const override;

  private:
    /**
     *  @brief Is a cluster "good" and worth keeping?
     *
     *  @param ClusterParameters   The cluster parameters of cluster to be checked
     *  @param Hit2DToClusterMap   Map to keep track of 2D hit to cluster association
     */

    bool keepThisCluster(reco::ClusterParameters&, const reco::Hit2DToClusterMap&) const;

    void storeThisCluster(reco::ClusterParameters&, reco::Hit2DToClusterMap&) const;

    void removeUsedHitsFromMap(reco::ClusterParameters&,
                               reco::HitPairListPtr&,
                               reco::Hit2DToClusterMap&) const;

    /**
     *  @brief Data members to follow
     */
    size_t m_clusterMinHits;
    double m_clusterMinUniqueFraction;
    double m_clusterMaxLostFraction;

    PrincipalComponentsAlg m_pcaAlg; // For running Principal Components Analysis
  };

  //------------------------------------------------------------------------------------------------------------------------------------------
  // implementation follows

  ClusterParamsBuilder::ClusterParamsBuilder(fhicl::ParameterSet const& pset)
    : m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
  {
    this->configure(pset);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  ClusterParamsBuilder::~ClusterParamsBuilder() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void ClusterParamsBuilder::configure(fhicl::ParameterSet const& pset)
  {
    m_clusterMinHits = pset.get<size_t>("ClusterMinHits", 3);
    m_clusterMinUniqueFraction = pset.get<double>("ClusterMinUniqueFraction", 0.5);
    m_clusterMaxLostFraction = pset.get<double>("ClusterMaxLostFraction", 0.5);

    return;
  }

  void ClusterParamsBuilder::BuildClusterInfo(
    reco::ClusterParametersList& clusterParametersList) const
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
    if (!clusterParametersList.empty()) {
      // We want to order our clusters on by largest (most number hits) to smallest. So, we'll loop through the clusters,
      // weeding out the unwanted ones and keep track of things in a set of "good" clusters which we'll order
      // by cluster size.
      clusterParametersList.sort();

      // The smallest clusters are now at the end, drop those off the back that are less than the mininum necessary
      while (!clusterParametersList.empty() &&
             clusterParametersList.back().getHitPairListPtr().size() < m_clusterMinHits)
        clusterParametersList.pop_back();

      // The next step is to build out a mapping of all 2D hits to clusters
      // Keep track of where the hits get distributed...
      reco::Hit2DToClusterMap hit2DToClusterMap;

      reco::ClusterParametersList::iterator clusterItr = clusterParametersList.begin();

      //        for(auto& clusterParams : clusterParametersList)
      //        {
      //            for(const auto& hit3D : clusterParams.getHitPairListPtr())
      //            {
      //                for(const auto& hit2D : hit3D->getHits())
      //                {
      //                    if (!hit2D) continue;
      //
      //                    hit2DToClusterMap[hit2D][&clusterParams].insert(hit3D);
      //                }
      //            }
      //        }

      // Ok, spin through again to remove ambiguous hits
      //        for(auto& clusterParams : clusterParametersList) PruneAmbiguousHits(clusterParams,hit2DToClusterMap);

      // What remains is an order set of clusters, largest first
      // Now go through and obtain cluster parameters
      //        clusterItr = clusterParametersList.begin();

      //        while(clusterItr != clusterParametersList.end())
      //        {
      //            // Dereference for ease...
      //            reco::ClusterParameters& clusterParams = *clusterItr;
      //
      //            // Do the actual work of filling the parameters
      //            FillClusterParams(clusterParams, hit2DToClusterMap, m_clusterMinUniqueFraction, m_clusterMaxLostFraction);
      //
      //            // If this cluster is rejected then the parameters will be empty
      //            if (clusterParams.getClusterParams().empty() || !clusterParams.getFullPCA().getSvdOK())
      //            {
      //                clusterItr = clusterParametersList.erase(clusterItr);
      //            }
      //            else clusterItr++;
      //        }

      while (clusterItr != clusterParametersList.end()) {
        // Dereference for ease...
        reco::ClusterParameters& clusterParams = *clusterItr;

        if (keepThisCluster(clusterParams, hit2DToClusterMap)) {
          storeThisCluster(clusterParams, hit2DToClusterMap);
          clusterItr++;
        }
        else
          clusterItr = clusterParametersList.erase(clusterItr);
      }
    }

    return;
  }

  bool ClusterParamsBuilder::keepThisCluster(reco::ClusterParameters& clusterParams,
                                             const reco::Hit2DToClusterMap& hit2DToClusterMap) const
  {
    // Try to keep simple by looking at the 2D hits associated to the cluster and checking to see how many, by plane, are already
    // in use. Reject clusters where too many hits are shared.

    bool keepThisCluster = false;

    // Define some handy data structures for counting the number of times hits get used and shared
    using HitCountMap = std::unordered_map<const reco::ClusterHit2D*, int>;
    using PlaneHitCountMapVec = std::vector<HitCountMap>;

    PlaneHitCountMapVec totalPlaneHitCountMapVec(3); // counts total number of hits
    PlaneHitCountMapVec sharedPlaneHitCountMapVec(
      3); // this is the number shared with a bigger cluster
    PlaneHitCountMapVec uniquePlaneHitCountMapVec(
      3); // this is the number unique to this cluster (so far)

    // Go through the hits and check usage...
    for (const auto& hit3D : clusterParams.getHitPairListPtr()) {
      for (const auto& hit2D : hit3D->getHits()) {
        if (!hit2D) continue;

        size_t hitPlane = hit2D->WireID().Plane;

        totalPlaneHitCountMapVec[hitPlane][hit2D]++;

        reco::Hit2DToClusterMap::const_iterator hit2DToClusIter = hit2DToClusterMap.find(hit2D);

        if (hit2DToClusIter != hit2DToClusterMap.end()) {
          sharedPlaneHitCountMapVec[hitPlane][hit2D]++;
        }
        else
          uniquePlaneHitCountMapVec[hitPlane][hit2D]++;
      }
    }

    // First try... look at fractions of unique hits each plane
    std::vector<float> uniqueFractionVec(3, 0.);

    for (size_t idx = 0; idx < 3; idx++) {
      if (!totalPlaneHitCountMapVec[idx].empty())
        uniqueFractionVec[idx] = float(uniquePlaneHitCountMapVec[idx].size()) /
                                 float(totalPlaneHitCountMapVec[idx].size());
    }

    float overallFraction = uniqueFractionVec[0] * uniqueFractionVec[1] * uniqueFractionVec[2];
    float maxFraction = *std::max_element(uniqueFractionVec.begin(), uniqueFractionVec.end());

    if (maxFraction > 0.9 || overallFraction > 0.2) keepThisCluster = true;

    return keepThisCluster;
  }

  void ClusterParamsBuilder::storeThisCluster(reco::ClusterParameters& clusterParams,
                                              reco::Hit2DToClusterMap& hit2DToClusterMap) const
  {
    // See if we can avoid duplicates by temporarily transferring to a set
    std::unordered_set<const reco::ClusterHit2D*> hitSet;

    // first task is to mark the hits and update the hit to cluster mapping
    for (const auto& hit3D : clusterParams.getHitPairListPtr()) {
      for (const auto& hit2D : hit3D->getHits()) {
        if (!hit2D) continue;

        hit2DToClusterMap[hit2D][&clusterParams].insert(hit3D);
        hitSet.insert(hit2D);
      }
    }

    // First stage of feature extraction runs here
    m_pcaAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());

    // Must have a valid pca
    if (clusterParams.getFullPCA().getSvdOK()) {
      // Set the skeleton PCA to make sure it has some value
      clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();

      // Add the "good" hits to our cluster parameters
      for (const auto& hit2D : hitSet) {
        hit2D->setStatusBit(reco::ClusterHit2D::USED);
        clusterParams.UpdateParameters(hit2D);
      }
    }

    return;
  }

  void ClusterParamsBuilder::FillClusterParams(reco::ClusterParameters& clusterParams,
                                               reco::Hit2DToClusterMap& hit2DToClusterMap,
                                               double minUniqueFrac,
                                               double maxLostFrac) const
  {
    /**
     *  @brief Given a list of hits fill out the remaining parameters for this cluster and evaluate the
     *         candidate's worthiness to achieve stardom in the event display
     */
    // Recover the HitPairListPtr from the input clusterParams (which will be the
    // only thing that has been provided)
    reco::HitPairListPtr& hitPairVector = clusterParams.getHitPairListPtr();

    // To be sure, we should clear the other data members
    clusterParams.getClusterParams().clear();
    clusterParams.getFullPCA() = reco::PrincipalComponents();

    // A test of the emergency broadcast system...
    //    FindBestPathInCluster(clusterParams);
    //    CheckHitSorting(clusterParams);

    // See if we can avoid duplicates by temporarily transferring to a set
    std::set<const reco::ClusterHit2D*> hitSet;

    // Ultimately we want to keep track of the number of unique 2D hits in this cluster
    // So use a vector (by plane) of sets of hits
    std::vector<size_t> planeHit2DVec;
    std::vector<size_t> planeUniqueHit2DVec;

    planeHit2DVec.resize(3);
    planeUniqueHit2DVec.resize(3);

    // Map from 2D hits to associated 3D hits
    reco::Hit2DToHit3DListMap& hit2DToHit3DListMap = clusterParams.getHit2DToHit3DListMap();

    // The map from 2D to 3D hits will contain unique entries for 2D hits so we can do some quick accounting here
    for (const auto& hitMapPair : hit2DToHit3DListMap) {
      size_t plane = hitMapPair.first->WireID().Plane;

      planeHit2DVec[plane] += hitMapPair.second.size();
      if (!(hitMapPair.first->getStatusBits() & reco::ClusterHit2D::USED))
        planeUniqueHit2DVec[plane] += hitMapPair.second.size();
    }

    // Get totals
    int numTotalHits(0);
    int numUniqueHits(0);

    // Also consider the number of hits shared on a given view...
    std::vector<float> uniqueHitFracVec(3, 0.);
    int nPlanesWithHits(0);
    int nPlanesWithUniqueHits(0);
    size_t minPlane(0);
    size_t minPlaneCnt = planeUniqueHit2DVec[0];

    // Loop through the planes
    for (int idx = 0; idx < 3; idx++) {
      // numerology
      numTotalHits += planeHit2DVec[idx];
      numUniqueHits += planeUniqueHit2DVec[idx];

      if (planeHit2DVec[idx] > 0) nPlanesWithHits++;
      if (planeUniqueHit2DVec[idx] > 0) nPlanesWithUniqueHits++;

      // Compute the fraction of unique hits in this plane
      uniqueHitFracVec[idx] =
        float(planeUniqueHit2DVec[idx]) / std::max(float(planeHit2DVec[idx]), float(1.));

      // Finding the plane with the fewest hits
      if (planeHit2DVec[idx] < minPlaneCnt) {
        minPlaneCnt = planeHit2DVec[idx];
        minPlane = idx;
      }
    }

    // If we have something left then at this point we make one more check
    // This check is intended to weed out clusters made from isolated groups of ambiguous hits which
    // really belong to a larger cluster
    if (numUniqueHits > 0.25 * numTotalHits && nPlanesWithHits > 1 && nPlanesWithUniqueHits > 1) {
      // Sorts lowest to highest
      std::sort(uniqueHitFracVec.begin(), uniqueHitFracVec.end());

      float acceptRatio = 0.;

      if (uniqueHitFracVec[0] * uniqueHitFracVec[1] > 0.25) acceptRatio = 1.;

      float uniqueFraction = uniqueHitFracVec[0] * uniqueHitFracVec[1] * uniqueHitFracVec[2];

      // Arbitrary rejection criteria... need to understand
      if (uniqueFraction > 0.6 && acceptRatio > 0.) {
        // Create a list to hold 3D hits which are already in use (criteria below)
        reco::HitPairListPtr usedHitPairList;

        //            std::cout << "--------> Starting 3D hit removal, # 2D hits/plane: " << planeHit2DVec[0] << "/" << planeHit2DVec[1] << "/" << planeHit2DVec[2] << std::endl;

        // Have survived laugh test, do final processing...
        // In this first loop go through all the 2D hits and identify the 3D hits that are candidates for deletion
        for (auto& pair : hit2DToHit3DListMap) {
          // Check to see if this can happen
          if (pair.second.empty()) {
            std::cout << "<<<<<< no matching 3D hits for reco hit in final hit processing >>>>>>"
                      << std::endl;
            continue;
          }

          // Which plane for this hit?
          size_t hitPlane = pair.first->WireID().Plane;

          // Only reject hits on the planes not the fewest 2D hits and really only do this if more than a couple
          if (hitPlane != minPlane && pair.second.size() > 2) {
            // If this hit is associated to a number of 3D hits then do some arbitration
            // Start by sorting the 3D hits by "significance"
            // --> Really should do this by the significance of adding the hit we are looking at?
            //pair.second.sort([hitPlane](const auto& left, const auto& right){return left->getHitDelTSigVec()[hitPlane] < right->getHitDelTSigVec()[hitPlane];});
            pair.second.sort([](const auto& left, const auto& right) {
              return left->getHitChiSquare() < right->getHitChiSquare();
            });

            ////std::cout << "~~~~> Checking hit removal, # matches: " << pair.second.size() << ", first params: " << pair.second.front()->getHitDelTSigVec()[hitPlane] << ", last params: "<< pair.second.back()->getHitDelTSigVec()[hitPlane];
            //std::cout << "~~~~> Checking hit removal, # matches: " << pair.second.size() << ", first params: " << pair.second.front()->getHitChiSquare() << ", last params: "<< pair.second.back()->getHitChiSquare();

            // From sorted list, determine a rejection value to eliminate bad hits
            //float cutDeltaTSig = std::min(2.0,std::max(0.5, double((pair.second.front()->getHitDelTSigVec()[hitPlane]))));
            float cutDeltaTSig =
              std::min(2.0, std::max(0.5, double(pair.second.front()->getHitChiSquare())));

            //std::cout << ", cutDeltaTSig: " << cutDeltaTSig;

            cutDeltaTSig = 10.;

            // And here go through the process of eliminating it
            //reco::HitPairListPtr::iterator firstBadHitItr = std::find_if(pair.second.begin(),pair.second.end(),[hitPlane,cutDeltaTSig](const auto& hitPtr){return hitPtr->getHitDelTSigVec()[hitPlane] > cutDeltaTSig;});
            reco::HitPairListPtr::iterator firstBadHitItr = std::find_if(
              pair.second.begin(), pair.second.end(), [cutDeltaTSig](const auto& hitPtr) {
                return hitPtr->getHitChiSquare() > cutDeltaTSig;
              });

            // We need to worry about cutting too many hits... use this loop to try to expand the range in a reasonable fashion
            //                    while(std::distance(pair.second.begin(),firstBadHitItr) < int(pair.second.size()/3) && cutDeltaTSig < 0.5)
            //                    {
            //                        float candDeltaTSig = (*firstBadHitItr)->getHitDelTSigVec()[hitPlane];
            //
            //                        if (candDeltaTSig > 2. * cutDeltaTSig) break;
            //
            //                        firstBadHitItr++;
            //                        cutDeltaTSig = candDeltaTSig;
            //                    }

            reco::HitPairListPtr rejectCandList;

            std::copy(firstBadHitItr, pair.second.end(), std::back_inserter(rejectCandList));

            //std::cout << ", bad hits: " << rejectCandList.size() << std::endl;

            // Remove the 3D hits from all the lists
            for (const auto& hit3D : rejectCandList) {
              bool rejectThisHit(true);
              std::vector<std::pair<reco::HitPairListPtr&, reco::HitPairListPtr::iterator>>
                deleteVec;

              for (const auto& hit2D : hit3D->getHits()) {
                // Watch for null hit (dead channels)
                if (!hit2D) continue;

                reco::HitPairListPtr& removeHitList = hit2DToHit3DListMap[hit2D];

                // Don't allow all the 3D hits associated to this 2D hit to be rejected?
                if (removeHitList.size() < 2) {
                  //std::cout << "   ---> remove list too small, size: " << removeHitList.size() << " for hit: " << hit2D << ", pair.first: " << pair.first << std::endl;
                  rejectThisHit = false;
                  break;
                }

                reco::HitPairListPtr::iterator removeItr =
                  std::find(removeHitList.begin(), removeHitList.end(), hit3D);

                if (removeItr != removeHitList.end())
                  deleteVec.emplace_back(removeHitList, removeItr);
                //else std::cout << "======>> Did not find 3D hit to remove from list for 2D hit! <<+++++++++" << std::endl;
              }

              if (rejectThisHit) {
                for (auto& rejectPair : deleteVec)
                  rejectPair.first.erase(rejectPair.second);

                usedHitPairList.push_back(hit3D);
              }
            }
          }

          hitSet.insert(pair.first);
        }

        // Now we go through the list of candidates and delete those which are unworthy of being processed...
        if (!usedHitPairList.empty()) {
          // Loop through the hits watching out for double counting
          const reco::ClusterHit3D* lastHit3D = 0;

          for (const auto& hit3D : usedHitPairList) {
            if (hit3D == lastHit3D) continue;

            reco::HitPairListPtr::iterator hit3DItr =
              std::find(hitPairVector.begin(), hitPairVector.end(), hit3D);

            if (hit3DItr != hitPairVector.end()) {
              // Mark the hit
              hit3D->setStatusBit(reco::ClusterHit3D::REJECTEDHIT);

              // Remove from the cluster's hit container
              hitPairVector.erase(hit3DItr);

              // If the clustering algorithm includes edges then need to get rid of those as well
              if (!clusterParams.getHit3DToEdgeMap().empty()) {
                reco::Hit3DToEdgeMap& edgeMap = clusterParams.getHit3DToEdgeMap();

                edgeMap.erase(edgeMap.find(hit3D));
              }
            }

            lastHit3D = hit3D;
          }
          //                removeUsedHitsFromMap(clusterParams, usedHitPairList, hit2DToClusterMap);
        }

        // First stage of feature extraction runs here
        m_pcaAlg.PCAAnalysis_3D(hitPairVector, clusterParams.getFullPCA());

        // Must have a valid pca
        if (clusterParams.getFullPCA().getSvdOK()) {
          // Set the skeleton PCA to make sure it has some value
          clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();

          // Add the "good" hits to our cluster parameters
          for (const auto& hit2D : hitSet) {
            hit2D->setStatusBit(reco::ClusterHit2D::USED);
            clusterParams.UpdateParameters(hit2D);
          }
        }
      }
    }

    return;
  }

  void ClusterParamsBuilder::removeUsedHitsFromMap(reco::ClusterParameters& clusterParams,
                                                   reco::HitPairListPtr& usedHitPairList,
                                                   reco::Hit2DToClusterMap& hit2DToClusterMap) const
  {
    // Clean up our hit to cluster map
    for (const auto& hit3D : usedHitPairList) {
      for (const auto& hit2D : hit3D->getHits()) {
        if (!hit2D) continue;

        reco::Hit2DToClusterMap::iterator hitToClusMapItr = hit2DToClusterMap.find(hit2D);

        // I am pretty sure this can't happen but let's check anyway...
        if (hitToClusMapItr == hit2DToClusterMap.end()) {
          std::cout << "*********** COULD NOT FIND ENTRY FOR 2D HIT! **************" << std::endl;
          break;
        }

        // Ok, the same hit can be shared in the same cluster so must be careful here
        // First loop over clusters looking for match
        reco::ClusterToHitPairSetMap::iterator clusToHit3DMapItr =
          hitToClusMapItr->second.find(&clusterParams);

        // This also can't happen
        if (clusToHit3DMapItr == hitToClusMapItr->second.end()) {
          std::cout << "************ DUCK! THE SKY HAS FALLEN!! *********" << std::endl;
          break;
        }

        // If this hit is shared by more than one 3D hit then pick the right one
        if (clusToHit3DMapItr->second.size() > 1) {
          reco::HitPairSetPtr::iterator hit3DItr = clusToHit3DMapItr->second.find(hit3D);

          clusToHit3DMapItr->second.erase(hit3DItr);
        }
        else
          hitToClusMapItr->second.erase(clusToHit3DMapItr);
      }
    }

    return;
  }

  DEFINE_ART_CLASS_TOOL(ClusterParamsBuilder)
} // namespace lar_cluster3d

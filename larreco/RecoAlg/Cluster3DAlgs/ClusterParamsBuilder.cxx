/**
 *  @file   Cluster3D_module.cc
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "cetlib/search_path.h"

#include "larreco/RecoAlg/Cluster3DAlgs/ClusterParamsBuilder.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

ClusterParamsBuilder::ClusterParamsBuilder(fhicl::ParameterSet const &pset) :
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterParamsBuilder::~ClusterParamsBuilder()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterParamsBuilder::configure(fhicl::ParameterSet const &pset)
{
    m_clusterMinHits           = pset.get<size_t>("ClusterMinHits",            3     );
    m_clusterMinUniqueFraction = pset.get<double>("ClusterMinUniqueFraction",  0.5   );
    m_clusterMaxLostFraction   = pset.get<double>("ClusterMaxLostFraction",    0.5   );
    
    return;
}
    
void ClusterParamsBuilder::BuildClusterInfo(reco::ClusterParametersList& clusterParametersList) const
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
    if (!clusterParametersList.empty())
    {
        // We want to order our clusters on by largest (most number hits) to smallest. So, we'll loop through the clusters,
        // weeding out the unwanted ones and keep track of things in a set of "good" clusters which we'll order
        // by cluster size.
        clusterParametersList.sort();
        
        // The smallest clusters are now at the end, drop those off the back that are less than the mininum necessary
        while(!clusterParametersList.empty() && clusterParametersList.back().getHitPairListPtr().size() < m_clusterMinHits) clusterParametersList.pop_back();
        
        // The next step is to build out a mapping of all 2D hits to clusters
        // Keep track of where the hits get distributed...
        reco::Hit2DToClusterMap hit2DToClusterMap;
        
        reco::ClusterParametersList::iterator clusterItr = clusterParametersList.begin();
        
        for(auto& clusterParams : clusterParametersList)
        {
            for(const auto& hit3D : clusterParams.getHitPairListPtr())
            {
                for(const auto& hit2D : hit3D->getHits())
                {
                    if (!hit2D) continue;
                    
                    hit2DToClusterMap[hit2D][&clusterParams].insert(hit3D);
                }
            }
        }
        
        // Ok, spin through again to remove ambiguous hits
        //        for(auto& clusterParams : clusterParametersList) PruneAmbiguousHits(clusterParams,hit2DToClusterMap);
        
        // What remains is an order set of clusters, largest first
        // Now go through and obtain cluster parameters
        clusterItr = clusterParametersList.begin();
        
        while(clusterItr != clusterParametersList.end())
        {
            // Dereference for ease...
            reco::ClusterParameters& clusterParams = *clusterItr;
            
            // Do the actual work of filling the parameters
            FillClusterParams(clusterParams, hit2DToClusterMap, m_clusterMinUniqueFraction, m_clusterMaxLostFraction);
            
            // If this cluster is rejected then the parameters will be empty
            if (clusterParams.getClusterParams().empty() || !clusterParams.getFullPCA().getSvdOK())
            {
                clusterItr = clusterParametersList.erase(clusterItr);
            }
            else clusterItr++;
        }
    }
    
    return;
}

void ClusterParamsBuilder::FillClusterParams(reco::ClusterParameters& clusterParams,
                                             reco::Hit2DToClusterMap& hit2DToClusterMap,
                                             double                   minUniqueFrac,
                                             double                   maxLostFrac) const
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
    //std::vector<const reco::ClusterHit2D*> hitSet;
    
    std::vector<int> nTotalHitsVec  = {0,0,0}; // Keep track of the total number of hits per plane
    std::vector<int> nUniqueHitsVec = {0,0,0}; // Keep track of the number of unique hits per plane
    
    int              numLostHits(0);
    int              numUnique3D(0);
    int              numOneShared3D(0);
    int              numTwoShared3D(0);
    int              numThreeShared3D(0);
    int              numAllShared3D(0);
    
    int              otherClusterMaxSize(0);
//    int              otherClusterMaxSharedHits(0);
    int              otherClusterMinSize(std::numeric_limits<int>::max());
//    int              otherClusterMinSharedHits(0);
    
    // Ultimately we want to keep track of the number of unique 2D hits in this cluster
    // So use a vector (by plane) of sets of hits
    std::vector<std::set<const reco::ClusterHit2D*>> planeHit2DSetVec;
    std::vector<std::set<const reco::ClusterHit2D*>> planeUniqueHit2DSetVec;

    planeHit2DSetVec.resize(3);
    planeUniqueHit2DSetVec.resize(3);

    // Keep track of associated clusters
    std::set<const reco::ClusterParameters*> totalOtherClusterSet;
    
    // Map from 2D hits to associated 3D hits
    std::map<const reco::ClusterHit2D*, reco::HitPairListPtr> hit2DToHit3DListMap;
    
    // Create a list to hold 3D hits which are already in use (criteria below)
    reco::HitPairListPtr usedHitPairList;
    
    std::cout << "*** building cluster, starting with " << hitPairVector.size() << " 3D hits" << std::endl;

    // Loop through all the 3D hits in this cluster and determine the total number of hits per plane,
    // the number of unique and lost hits.
    for(const auto& hit3D : hitPairVector)
    {
        int nHits2D(0);
        int nSharedHits(0);
        int nUniqueHits(0);
        int nSharedHits_test(0);
        
        // Keep track of associated clusters
        std::set<const reco::ClusterParameters*> associatedClusterSet;
        
        // First loop through hits to categorize
        for(const auto& hit2D : hit3D->getHits())
        {
            if (!hit2D) continue;
            
            size_t plane = hit2D->getHit().WireID().Plane;

            // Check if this hit is shared on another cluster
            if (hit2DToClusterMap.at(hit2D).size() > 1)
            {
                for(const auto& clusterPair : hit2DToClusterMap.at(hit2D))
                {
                    if (clusterPair.first != &clusterParams)
                    {
                        associatedClusterSet.insert(clusterPair.first);
                        totalOtherClusterSet.insert(clusterPair.first);
                        
                        if (int(clusterPair.first->getHitPairListPtr().size()) > otherClusterMaxSize)
                        {
                            otherClusterMaxSize       = clusterPair.first->getHitPairListPtr().size();
//                            otherClusterMaxSharedHits = clusterPair.second.size();
                        }
                        if (int(clusterPair.first->getHitPairListPtr().size()) < otherClusterMinSize)
                        {
                            otherClusterMinSize       = clusterPair.first->getHitPairListPtr().size();
//                            otherClusterMinSharedHits = clusterPair.second.size();
                        }
                    }
                }
                
                nSharedHits_test++;
            }
            else nUniqueHits++;

            if (hit2D->getStatusBits() & reco::ClusterHit2D::USED) nSharedHits++;
            
            nTotalHitsVec.at(plane)++;
            
            nHits2D++;
            
            // Keep track of the hits on each plane
            planeHit2DSetVec[plane].insert(hit2D);
            if (!(hit2D->getStatusBits() & reco::ClusterHit2D::USED)) planeUniqueHit2DSetVec[plane].insert(hit2D);
        }
        
        if      (nUniqueHits == nHits2D) numUnique3D++;
        else if (nSharedHits_test == 1)  numOneShared3D++;
        else if (nSharedHits_test == 2)  numTwoShared3D++;
        else                             numThreeShared3D++;
        
        if (nSharedHits == nHits2D) numAllShared3D++;
        
        // If more than one shared hit then candidate for removal...
        if (nSharedHits_test > 4 && nSharedHits > 0)
        {
            usedHitPairList.emplace_back(hit3D);
            
            numLostHits += nHits2D;
        }
        else
        {
            // One more intermediate step... make map of 2D to 3D hits for final stage
            for(const auto& hit2D : hit3D->getHits())
            {
                if (!hit2D) continue;
            
                hit2DToHit3DListMap[hit2D].push_back(hit3D);
                
                size_t plane = hit2D->getHit().WireID().Plane;
                
                if (!(hit2D->getStatusBits() & reco::ClusterHit2D::USED)) nUniqueHitsVec.at(plane)++;
            }
        }
    }

    // Get totals
    int numTotalHits(0);
    int numUniqueHits(0);
    
    for(int idx = 0; idx < 3; idx++)
    {
        numTotalHits += planeHit2DSetVec.at(idx).size();
        numUniqueHits += planeUniqueHit2DSetVec.at(idx).size();
    }
    
    // If we have something left then at this point we make one more check
    // This check is intended to weed out clusters made from isolated groups of ambiguous hits which
    // really belong to a larger cluster
    if (numUniqueHits > 0)
    {
        // Look at reject to accept ratio
        //double rejectToAccept = double(numRejected) / double(numAccepted);
        float acceptRatio = float(numUniqueHits) / float(numTotalHits);
        float lostRatio   = float(numLostHits)   / float(numTotalHits);
        
        // Also consider the number of hits shared on a given view...
        std::vector<float> uniqueHitVec(3,0.);
        
        for(size_t idx = 0; idx < 3; idx++) uniqueHitVec[idx] = float(planeUniqueHit2DSetVec.at(idx).size()) / std::max(float(planeHit2DSetVec.at(idx).size()),float(1.));
        
        // Sorts lowest to highest
        std::sort(uniqueHitVec.begin(),uniqueHitVec.end());
        
        acceptRatio = 0.;
        lostRatio   = 0.;
        //if(uniqueHitVec[1] > 0.1 && uniqueHitVec[2] > 0.5) acceptRatio = 1.;
        if(uniqueHitVec[0] * uniqueHitVec[1] > 0.25) acceptRatio = 1.;
        
//        float allSharedFraction = numAllShared3D / float(hitPairVector.size());
//        float uniqueFraction = float(numUniqueHits) / float(numTotalHits);
        float uniqueFraction = uniqueHitVec[0] * uniqueHitVec[1] * uniqueHitVec[2];

//        std::cout << "**--> # 3D Hits: " << hitPairVector.size() << ", nTot: " << numTotalHits << " " << nTotalHitsVec[0] << "/" << nTotalHitsVec[1] << "/" << nTotalHitsVec[2] << ", unique: " << numUniqueHits << " " << nUniqueHitsVec[0] << "/" << nUniqueHitsVec[1] << "/" << nUniqueHitsVec[2] << ", lost: " << usedHitPairList.size() << ", accept: " << acceptRatio << ", rats: " << uniqueHitVec[0] << "/" << uniqueHitVec[1] << "/" << uniqueHitVec[2] << std::endl;
//        std::cout << "      -- nUnique3D: " << numUnique3D << ", 1 shared: " << numOneShared3D << ", 2 shared: " << numTwoShared3D << ", 3 shared: " << numThreeShared3D << ", all shared: " << numAllShared3D << std::endl;
//        std::cout << "      -- total # other clusters: " << totalOtherClusterSet.size() << ", max size: " << otherClusterMaxSize << ", shared: " << otherClusterMaxSharedHits << ", min size: " << otherClusterMinSize << ", shared: " << otherClusterMinSharedHits << std::endl;
        
        // Arbitrary rejection criteria... need to understand
        // Anyway, if we get past this we're making a cluster
        //if (rejectToAccept < rejectFraction)
        //if (numUnique3D+numOneShared3D > 0 && allSharedFraction < 1.0 && acceptRatio > minUniqueFrac && lostRatio < maxLostFrac)  // lostRatio cut was 1. - off
        if (uniqueFraction > 0.6 && acceptRatio > 0. && lostRatio > -1.)
        {
            int nPlanesWithHits    = std::accumulate(nUniqueHitsVec.begin(), nUniqueHitsVec.end(), 0, [](int accum, int total){return total > 0 ? ++accum : accum;});
            int nPlanesWithMinHits = std::accumulate(nUniqueHitsVec.begin(), nUniqueHitsVec.end(), 0, [](int accum, int total){return total > 2 ? ++accum : accum;});
            
//            std::cout << "      --> nPlanesWithHits: " << nPlanesWithHits << ", nPlanesWithMinHits: " << nPlanesWithMinHits << ", unique prod: " << uniqueHitVec[1] * uniqueHitVec[2] << std::endl;

            //            // Final selection cut, need at least 3 hits each view
            //            if (nViewsWithHits == 3 && nViewsWithMinHits > 1)
            // Final selection cut, need at least 3 hits each view for at least 2 views
            //if (nPlanesWithHits > 1 && nPlanesWithMinHits > 1)
            if (nPlanesWithHits > -1 && nPlanesWithMinHits > -1)
            {
                std::cout << "--------> Starting 3D hit removal, # 2D hits/plane: " << planeHit2DSetVec[0].size() << "/" << planeHit2DSetVec[1].size() << "/" << planeHit2DSetVec[2].size() << std::endl;
                // Have survived laugh test, do final processing...
                // In this first loop go through all the 2D hits and identify the 3D hits that are candidates for deletion
                for(auto& pair : hit2DToHit3DListMap)
                {
                    // Check to see if this can happen
                    if (pair.second.empty())
                    {
                        std::cout << "<<<<<< no matching 3D hits for reco hit in final hit processing >>>>>>" << std::endl;
                        continue;
                    }
                    
                    // No point in doing anything if not more than one entry
                    if (pair.second.size() < 2) continue;

                    // If this hit is associated to a number of 3D hits then do some arbitration
                    pair.second.sort([](const auto& left, const auto& right){return left->getDeltaPeakTime()/left->getSigmaPeakTime() < right->getDeltaPeakTime()/right->getSigmaPeakTime();});
                    
                    std::cout << "~~~~> Checking hit removal, # matches: " << pair.second.size() << ", first params: " << pair.second.front()->getDeltaPeakTime() << "/" << pair.second.front()->getSigmaPeakTime() << ", last params: "<< pair.second.back()->getDeltaPeakTime() << "/" << pair.second.back()->getSigmaPeakTime();
                    
                    //float cutDeltaTSig = std::min(1.1,std::max(0.2, 3. * (pair.second.front()->getDeltaPeakTime() / pair.second.front()->getSigmaPeakTime())));
                    float cutDeltaTSig = std::min(1.7,std::max(0.5, 2.0 * (pair.second.front()->getDeltaPeakTime() / pair.second.front()->getSigmaPeakTime())));
                    
                    std::cout << ", cutDeltaTSig: " << cutDeltaTSig;
                    
                    reco::HitPairListPtr::iterator firstBadHitItr = std::find_if(pair.second.begin(),pair.second.end(),[cutDeltaTSig](const auto& hitPtr){return hitPtr->getDeltaPeakTime() > cutDeltaTSig * hitPtr->getSigmaPeakTime();});
                    
                    reco::HitPairListPtr rejectCandList;
                    
                    std::copy(firstBadHitItr,pair.second.end(),std::back_inserter(rejectCandList));
                    
                    std::cout << ", bad hits: " << rejectCandList.size() << std::endl;
                    
                    // Remove the 3D hits from all the lists
                    for(const auto& hit3D : rejectCandList)
                    {
                        bool rejectThisHit(true);
                        std::vector<std::pair<reco::HitPairListPtr&,reco::HitPairListPtr::iterator>> deleteVec;
                        
                        for(const auto& hit2D : hit3D->getHits())
                        {
                            // Watch for null hit (dead channels)
                            if (!hit2D) continue;
                            
                            reco::HitPairListPtr& removeHitList = hit2DToHit3DListMap.at(hit2D);

                            // Don't allow all the 3D hits associated to this 2D hit to be rejected?
                            if (removeHitList.size() > 0) //< 2)
                            {
                                std::cout << "   ---> remove list too small, size: " << removeHitList.size() << " for hit: " << hit2D << ", pair.first: " << pair.first << std::endl;
                                rejectThisHit = false;
                                break;
                            }
                            
                            reco::HitPairListPtr::iterator removeItr = std::find(removeHitList.begin(),removeHitList.end(),hit3D);
                            
                            if (removeItr != removeHitList.end()) deleteVec.emplace_back(removeHitList,removeItr);
                            else std::cout << "======>> Did not find 3D hit to remove from list for 2D hit! <<+++++++++" << std::endl;
                        }
                        
                        if (rejectThisHit)
                        {
                            for(auto& rejectPair : deleteVec) rejectPair.first.erase(rejectPair.second);
                            
                            usedHitPairList.push_back(hit3D);
                        }
                    }

                    hitSet.insert(pair.first);
                }
                
                // Now we go through the list of candidates and delete those which are unworthy of being processed...
                if (!usedHitPairList.empty())
                {
                    // Loop through the hits watching out for double counting
                    const reco::ClusterHit3D* lastHit3D = 0;
                    
                    for(const auto& hit3D : usedHitPairList)
                    {
                        if (hit3D == lastHit3D) continue;
                        
                        // Loop through the 2D hits to cross check that each has more than one 3D hit associated to it
                        bool okToDelete(true);
/*
                        for(const auto& hit2D : hit3D->getHits())
                        {
                            if (!hit2D) continue;
                            
                            if (hit2DToHit3DListMap.at(hit2D).size() < 2)
                            {
                                okToDelete = false;
                                break;
                            }
                        }
*/
                        if (okToDelete)
                        {
                            reco::HitPairListPtr::iterator hit3DItr = std::find(hitPairVector.begin(),hitPairVector.end(),hit3D);
                            
                            if (hit3DItr != hitPairVector.end())
                            {
                                // Mark the hit
                                hit3D->setStatusBit(reco::ClusterHit3D::REJECTEDHIT);
                                
                                // Remove from the cluster's hit container
                                hitPairVector.erase(hit3DItr);
                                
                                // If the clustering algorithm includes edges then need to get rid of those as well
                                if (!clusterParams.getHit3DToEdgeMap().empty())
                                {
                                    reco::Hit3DToEdgeMap& edgeMap = clusterParams.getHit3DToEdgeMap();
                                    
                                    edgeMap.erase(edgeMap.find(hit3D));
                                }
                            }
                        }
                        
                        lastHit3D = hit3D;
                    }
//                    removeUsedHitsFromMap(clusterParams, usedHitPairList, hit2DToClusterMap);
                }

                // First stage of feature extraction runs here
                m_pcaAlg.PCAAnalysis_3D(hitPairVector, clusterParams.getFullPCA());
                
                // Must have a valid pca
                if (clusterParams.getFullPCA().getSvdOK())
                {
                    std::cout << "      ******> Valid Full PCA calculated over " << clusterParams.getHitPairListPtr().size() << " hits" << std::endl;
                    
                    // Set the skeleton PCA to make sure it has some value
                    clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();
                    
                    // Add the "good" hits to our cluster parameters
                    for(const auto& hit2D : hitSet)
                    {
                        hit2D->setStatusBit(reco::ClusterHit2D::USED);
                        clusterParams.UpdateParameters(hit2D);
                    }
                }
            }
        }
    }
    
    return;
}
    
void ClusterParamsBuilder::removeUsedHitsFromMap(reco::ClusterParameters& clusterParams,
                                                 reco::HitPairListPtr&    usedHitPairList,
                                                 reco::Hit2DToClusterMap& hit2DToClusterMap) const
{
    // Clean up our hit to cluster map
    for(const auto& hit3D : usedHitPairList)
    {
        for(const auto& hit2D : hit3D->getHits())
        {
            if (!hit2D) continue;
            
            reco::Hit2DToClusterMap::iterator hitToClusMapItr = hit2DToClusterMap.find(hit2D);
            
            // I am pretty sure this can't happen but let's check anyway...
            if (hitToClusMapItr == hit2DToClusterMap.end())
            {
                std::cout << "*********** COULD NOT FIND ENTRY FOR 2D HIT! **************" << std::endl;
                break;
            }
            
            // Ok, the same hit can be shared in the same cluster so must be careful here
            // First loop over clusters looking for match
            reco::ClusterToHitPairSetMap::iterator clusToHit3DMapItr = hitToClusMapItr->second.find(&clusterParams);
            
            // This also can't happen
            if (clusToHit3DMapItr == hitToClusMapItr->second.end())
            {
                std::cout << "************ DUCK! THE SKY HAS FALLEN!! *********" << std::endl;
                break;
            }
            
            // If this hit is shared by more than one 3D hit then pick the right one
            if (clusToHit3DMapItr->second.size() > 1)
            {
                reco::HitPairSetPtr::iterator hit3DItr = clusToHit3DMapItr->second.find(hit3D);
                
                clusToHit3DMapItr->second.erase(hit3DItr);
            }
            else
                hitToClusMapItr->second.erase(clusToHit3DMapItr);
            
        }
    }

    return;
}
    
} // namespace lar_cluster3d

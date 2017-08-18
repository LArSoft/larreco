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
    //std::set<const reco::ClusterHit2D*> hitSet;
    std::vector<const reco::ClusterHit2D*> hitSet;
    
    size_t nTotalHits[]  = {0,0,0};
    size_t nUniqueHits[] = {0,0,0};
    size_t nLostHits[]   = {0,0,0};
    size_t nMultShared2DHits(0);
    size_t nAllHitsShared(0);
    
    std::map<reco::ClusterParameters*,int> clusterHitCountMap;
    
    // Create a list to hold 3D hits which are already in use (criteria below)
    reco::HitPairListPtr usedHitPairList;
    
    // First loop through the 3D hits
    // The goal of this loop is to build a set of unique hits from the hit pairs (which may contain many
    // ambiguous duplicate combinations).
    // The secondary goal is to remove 3D hits marked by hit arbitration to be tossed
    for(const auto& hit3D : hitPairVector)
    {
        size_t nMultClusters(0);
        size_t nHits2D(0);
        size_t nHitsUsed[] = {0,0,0};
        
        // loop over the hits in this 3D Cluster hit
        for(const auto& hit2D : hit3D->getHits())
        {
            if (!hit2D) continue;
            
            size_t plane = hit2D->getHit().WireID().Plane;
            
            // If 2D hit is marked as already used then record as such...
            if (hit2D->getStatusBits() & reco::ClusterHit2D::USED) nHitsUsed[plane]++;
            else                                                   nUniqueHits[plane]++;
            
            // Is this 2D hit shared? (don't forget to ignore the self reference)
            if (hit2DToClusterMap[hit2D].size() > 1)             nMultClusters++;
            for(auto& clusterCntPair : hit2DToClusterMap[hit2D]) clusterHitCountMap[clusterCntPair.first]++;
            
            nTotalHits[plane]++;
            nHits2D++;
        }
        
        size_t nHitsAlreadyUsed = std::accumulate(nHitsUsed,nHitsUsed+3,0);
        
        if (nMultClusters > 0)        nMultShared2DHits++;
        if (nMultClusters == nHits2D) nAllHitsShared++;
        
        for(size_t idx=0;idx<3;idx++)
        {
            if (nHitsAlreadyUsed < nHits2D)
            {
                //if (hit3D->getHits()[idx]) hitSet.insert(hit3D->getHits()[idx]);
                if (hit3D->getHits()[idx]) hitSet.push_back(hit3D->getHits()[idx]);
            }
            else nLostHits[idx] += nHitsUsed[idx];
        }
        
        if (nHitsAlreadyUsed == nHits2D) usedHitPairList.emplace_back(hit3D);
    }
    
    int numTotal      = std::accumulate(nTotalHits,nTotalHits+3,0);
    int numUniqueHits = std::accumulate(nUniqueHits,nUniqueHits+3,0);
    int numLostHits   = std::accumulate(nLostHits,nLostHits+3,0);
    
    std::cout << "*********************************************************************" << std::endl;
    std::cout << "**--> cluster: " << &clusterParams << " has " << hitPairVector.size() << " 3D hits, " << numTotal << ", numUnique: " << numUniqueHits << " 2D hits, match: " << clusterHitCountMap[&clusterParams] << ", shared: " << nMultShared2DHits << ", all: " << nAllHitsShared << std::endl;
    
    for(const auto& clusterCnt : clusterHitCountMap)
    {
        if (clusterCnt.first == &clusterParams) continue;
        std::cout << "      --> cluster " << clusterCnt.first << ", # hits: " << clusterCnt.second << std::endl;
    }
    
    // If we have something left then at this point we make one more check
    // This check is intended to weed out clusters made from isolated groups of ambiguous hits which
    // really belong to a larger cluster
    if (numUniqueHits > 3 && nMultShared2DHits < hitPairVector.size())
    {
        // Look at reject to accept ratio
        //double rejectToAccept = double(numRejected) / double(numAccepted);
        double acceptRatio = double(numUniqueHits) / double(numTotal);
        double lostRatio   = double(numLostHits)   / double(numTotal);
        
        // Also consider the number of hits shared on a given view...
        std::vector<double> uniqueHitVec(3,0.);
        
        for(size_t idx = 0; idx < 3; idx++) uniqueHitVec[idx] = double(nUniqueHits[idx]) / std::max(double(nTotalHits[idx]),1.);
        
        std::sort(uniqueHitVec.begin(),uniqueHitVec.end());
        
        double midHitRatio = uniqueHitVec[1];
        
        std::cout << "**--> # 3D Hits: " << hitPairVector.size() << ", nTot: " << numTotal << ", unique: " << numUniqueHits << ", lost: " << numLostHits << ", accept: " << acceptRatio << ", lost: " << lostRatio << ", mid: " << midHitRatio << ", rats: " << uniqueHitVec[0] << "/" << uniqueHitVec[1] << "/" << uniqueHitVec[2] << std::endl;
        
        acceptRatio = 0.;
        lostRatio   = 0.;
        if(uniqueHitVec[1] > 0.1 && uniqueHitVec[2] > 0.5) acceptRatio = 1.;
        
        // Arbitrary rejection criteria... need to understand
        // Anyway, if we get past this we're making a cluster
        //if (rejectToAccept < rejectFraction)
        if (acceptRatio > minUniqueFrac && lostRatio < maxLostFrac)  // lostRatio cut was 1. - off
        {
            // Add the "good" hits to our cluster parameters
            for(const auto& hit2D : hitSet)
            {
                hit2D->setStatusBit(reco::ClusterHit2D::USED);
                clusterParams.UpdateParameters(hit2D);
            }
            
            size_t nPlanesWithHits    = (clusterParams.getClusterParams()[0].m_hitVector.size() > 0 ? 1 : 0)
            + (clusterParams.getClusterParams()[1].m_hitVector.size() > 0 ? 1 : 0)
            + (clusterParams.getClusterParams()[2].m_hitVector.size() > 0 ? 1 : 0);
            size_t nPlanesWithMinHits = (clusterParams.getClusterParams()[0].m_hitVector.size() > 2 ? 1 : 0)
            + (clusterParams.getClusterParams()[1].m_hitVector.size() > 2 ? 1 : 0)
            + (clusterParams.getClusterParams()[2].m_hitVector.size() > 2 ? 1 : 0);
            //            // Final selection cut, need at least 3 hits each view
            //            if (nViewsWithHits == 3 && nViewsWithMinHits > 1)
            // Final selection cut, need at least 3 hits each view for at least 2 views
            if (nPlanesWithHits > 1 && nPlanesWithMinHits > 1)
            {
                // First task is to remove the hits already in use
                if (!usedHitPairList.empty())
                {
                    hitPairVector.sort();
                    usedHitPairList.sort();
                    
                    reco::HitPairListPtr::iterator newListEnd =
                    std::set_difference(hitPairVector.begin(),   hitPairVector.end(),
                                        usedHitPairList.begin(), usedHitPairList.end(),
                                        hitPairVector.begin() );
                    
                    hitPairVector.erase(newListEnd, hitPairVector.end());
                }
                
                // First stage of feature extraction runs here
                m_pcaAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());
                
                // Must have a valid pca
                if (clusterParams.getFullPCA().getSvdOK())
                {
                    std::cout << "**--> Valid Full PCA calculated" << std::endl;
                    
                    // If any hits were thrown away, see if we can rescue them
                    if (!usedHitPairList.empty())
                    {
                        double maxDoca = 2. * sqrt(clusterParams.getFullPCA().getEigenValues()[1]);
                        
                        if (maxDoca < 5.)
                        {
                            size_t curHitVectorSize = hitPairVector.size();
                            
                            m_pcaAlg.PCAAnalysis_calc3DDocas(usedHitPairList, clusterParams.getFullPCA());
                            
                            for(const auto& hit3D : usedHitPairList)
                                if (hit3D->getDocaToAxis() < maxDoca) hitPairVector.push_back(hit3D);
                            
                            if (hitPairVector.size() > curHitVectorSize)
                                m_pcaAlg.PCAAnalysis_3D(clusterParams.getHitPairListPtr(), clusterParams.getFullPCA());
                        }
                    }
                    
                    // Set the skeleton PCA to make sure it has some value
                    clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();
                }
            }
        }
    }
    
    return;
}
    
} // namespace lar_cluster3d

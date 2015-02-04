/**
 *  @file   ParallelHitsSeedFinderAlg.cxx
 * 
 *  @brief  Implementation of the Seed Finder Algorithm
 *          The intent of this algorithm is to take an input list of 3D space points and from those
 *          to find candidate track start points and directions
 */

// The main include
#include "RecoAlg/Cluster3DAlgs/ParallelHitsSeedFinderAlg.h"
// Framework Includes

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoObjects/Cluster3D.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

// ROOT includes
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2D.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

ParallelHitsSeedFinderAlg::ParallelHitsSeedFinderAlg(fhicl::ParameterSet const &pset) :
     m_maxNumEdgeHits(300),
     m_gapDistance(5.),
     m_numSeed2DHits(80),
     m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->reconfigure(pset);
    
    art::ServiceHandle<geo::Geometry>            geometry;
    art::ServiceHandle<util::DetectorProperties> detectorProperties;
    
    m_geometry = &*geometry;
    m_detector = &*detectorProperties;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParallelHitsSeedFinderAlg::~ParallelHitsSeedFinderAlg()
{
}
    
void ParallelHitsSeedFinderAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_maxNumEdgeHits = pset.get<size_t>("MaxNumEdgeHits", 300);
    m_gapDistance    = pset.get<double>("GapDistance",     5.);
    m_numSeed2DHits  = pset.get<size_t>("NumSeed2DHits",   80);
    m_pcaAlg.reconfigure(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));
}
    
bool ParallelHitsSeedFinderAlg::findTrackSeeds(reco::HitPairListPtr&      inputHitPairListPtr,
                                               reco::PrincipalComponents& inputPCA,
                                               SeedHitPairListPairVec&    seedHitPairVec)
{
    // This routine can fail...
    bool foundGoodSeeds(false);
    
    // Make sure we are using the right pca
    reco::PrincipalComponents pca = inputPCA;
    
    if (pca.getSvdOK())
    {
        // This routine is typically called when there are LOTS of hits... so we are going to try
        // to reduce the number of operations on the full list as much as possible. However, we
        // really need the input hit list to be sorted by th input PCA so do that now
        m_pcaAlg.PCAAnalysis_calc3DDocas(inputHitPairListPtr, pca);
        
        // Use this info to sort the hits along the principle axis
        // Note that this will sort hits from the "middle" to the "outside"
        inputHitPairListPtr.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
        
        // Now, we like to keep 20% of the total number of hits at the ends
        size_t numHitsToDrop = 4 * inputHitPairListPtr.size() / 5;
        size_t numEdgeHits   = std::min(size_t((inputHitPairListPtr.size() - numHitsToDrop) / 2), m_maxNumEdgeHits);
        
        // Get an iterator for the hits
        reco::HitPairListPtr::iterator edgeHitItr = inputHitPairListPtr.begin();
        
        std::advance(edgeHitItr, numEdgeHits);
        
        // Make a container for copying in the edge hits and size it to hold all the hits
        reco::HitPairListPtr hit3DList;
        
        hit3DList.resize(2*numEdgeHits);
        
        // Copy the low edge hit range into the list
        reco::HitPairListPtr::iterator nextHit3DItr = std::copy(inputHitPairListPtr.begin(), edgeHitItr, hit3DList.begin());
        
        // Now advance the iterator into the main container and copy the rest of the elements
        std::advance(edgeHitItr, inputHitPairListPtr.size() - 2 * numEdgeHits);
        
        std::copy(edgeHitItr, inputHitPairListPtr.end(), nextHit3DItr);
        
        // Now rerun the PCA
        reco::PrincipalComponents edgePCA;
        
        m_pcaAlg.PCAAnalysis_3D(hit3DList, edgePCA, true);
        
        if (edgePCA.getSvdOK())
        {
            // Still looking to point "down"
            if (edgePCA.getEigenVectors()[0][1] > 0.) edgePCA.flipAxis(0);
            
            // Recompute the 3D docas and arc lengths
            m_pcaAlg.PCAAnalysis_calc3DDocas(hit3DList, edgePCA);
            
            // Now sort hits in regular order
            hit3DList.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
            
            // Use a set to count 2D hits by keeping only a single copy
            std::set<const reco::ClusterHit2D*> hit2DSet;
        
            // Look for numSeed2DHits which are continuous
            double lastArcLen = hit3DList.front()->getArclenToPoca();
        
            reco::HitPairListPtr::iterator startItr = hit3DList.begin();
            reco::HitPairListPtr::iterator lastItr  = hit3DList.begin();
        
            while(++lastItr != hit3DList.end())
            {
                const reco::ClusterHit3D* hit3D  = *lastItr;
                double                    arcLen = hit3D->getArclenToPoca();
            
                if (fabs(arcLen-lastArcLen) > m_gapDistance)
                {
                    startItr = lastItr;
                    hit2DSet.clear();
                }
            
                for(const auto& hit2D : hit3D->getHits())
                {
                    hit2DSet.insert(hit2D);
                }
           
                if (hit2DSet.size() > m_numSeed2DHits) break;
            
                lastArcLen = arcLen;
            }
        
            if (hit2DSet.size() > m_numSeed2DHits)
            {
                if (startItr != hit3DList.begin()) hit3DList.erase(hit3DList.begin(), startItr);
                if (lastItr  != hit3DList.end()  ) hit3DList.erase(lastItr,           hit3DList.end());
        
                reco::PrincipalComponents seedPca;
        
                m_pcaAlg.PCAAnalysis_3D(hit3DList, seedPca, true);
        
                if (seedPca.getSvdOK())
                {
                    // Still looking to point "down"
                    if (seedPca.getEigenVectors()[0][1] > 0.) seedPca.flipAxis(0);
                    
                    // We want our seed to be in the middle of our seed points
                    // First step to getting there is to reorder the hits along the
                    // new pca axis
                    m_pcaAlg.PCAAnalysis_calc3DDocas(hit3DList, seedPca);
                    
                    // Use this info to sort the hits along the principle axis - note by absolute value of arc length
                    hit3DList.sort(SeedFinderAlgBase::Sort3DHitsByAbsArcLen3D());
                    
                    // Now translate the seedCenter by the arc len to the first hit
                    double seedDir[3]   = {seedPca.getEigenVectors()[0][0], seedPca.getEigenVectors()[0][1], seedPca.getEigenVectors()[0][2]};
                    double seedStart[3] = {hit3DList.front()->getX(), hit3DList.front()->getY(), hit3DList.front()->getZ()};
                    
                    for(const auto& hit3D : hit3DList) hit3D->setStatusBit(0x40000000);
            
                    seedHitPairVec.emplace_back(std::pair<recob::Seed, reco::HitPairListPtr>(recob::Seed(seedStart, seedDir), hit3DList));
                    
                    inputPCA = edgePCA;
                    
                    foundGoodSeeds = true;
                }
            }
        }
    }
    
    return foundGoodSeeds;
}

} // namespace lar_cluster3d

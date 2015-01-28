/**
 *  @file   PCASeedFinderAlg.cxx
 * 
 *  @brief  Implementation of the Seed Finder Algorithm
 *
 *          The intent of this algorithm is to take an input list of 3D space points and from those
 *          to find candidate track start points and directions
 */

// The main include
#include "RecoAlg/Cluster3DAlgs/PCASeedFinderAlg.h"
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

PCASeedFinderAlg::PCASeedFinderAlg(fhicl::ParameterSet const &pset) :
     m_gapDistance(5.),
     m_numSeed2DHits(80),
     m_minAllowedCosAng(0.7),
     m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->reconfigure(pset);
    
    art::ServiceHandle<geo::Geometry>            geometry;
    art::ServiceHandle<util::DetectorProperties> detectorProperties;
    
    m_geometry = &*geometry;
    m_detector = &*detectorProperties;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PCASeedFinderAlg::~PCASeedFinderAlg()
{
}
    
void PCASeedFinderAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    m_gapDistance      = pset.get<double>("GapDistance",       5.);
    m_numSeed2DHits    = pset.get<size_t>("NumSeed2DHits",     80);
    m_minAllowedCosAng = pset.get<double>("MinAllowedCosAng", 0.7);
    m_pcaAlg.reconfigure(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));
}    
    
bool PCASeedFinderAlg::findTrackSeeds(reco::HitPairListPtr&      inputHitPairListPtr,
                                      reco::PrincipalComponents& inputPCA,
                                      SeedHitPairListPairVec&    seedHitPairVec)
{
    bool foundGoodSeed(false);
    
    // Make sure we are using the right pca
    reco::HitPairListPtr hitPairListPtr = inputHitPairListPtr;
    
    // Make a local copy of the input pca
    reco::PrincipalComponents pca = inputPCA;

    // We also require that there be some spread in the data, otherwise not worth running?
    double eigenVal0 = 3. * sqrt(pca.getEigenValues()[0]);
    double eigenVal1 = 3. * sqrt(pca.getEigenValues()[1]);
    
    if (eigenVal0 > 5. && eigenVal1 > 0.001)
    {
        // Presume CR muons will be "downward going"...
        if (pca.getEigenVectors()[0][1] > 0.) pca.flipAxis(0);
        
        // Use the following to set the 3D doca and arclength for each hit
        m_pcaAlg.PCAAnalysis_calc3DDocas(hitPairListPtr, pca);
        
        // Use this info to sort the hits along the principle axis
        hitPairListPtr.sort(SeedFinderAlgBase::Sort3DHitsByArcLen3D());
        //hitPairListPtr.sort(PcaSort3DHitsByAbsArcLen3D());
        
        // Make a local copy of the 3D hits
        reco::HitPairListPtr hit3DList;
        
        hit3DList.resize(hitPairListPtr.size());
        
        std::copy(hitPairListPtr.begin(), hitPairListPtr.end(), hit3DList.begin());
        
        reco::PrincipalComponents seedPca = pca;
        
        if (getHitsAtEnd(hit3DList, seedPca))
        {
            // We can use the same routine to check hits at the opposite end to make sure
            // we have consistency between both ends of the track.
            // So... follow the same general program as above
            reco::HitPairListPtr checkList;
            
            checkList.resize(hitPairListPtr.size());
            
            std::copy(hitPairListPtr.begin(), hitPairListPtr.end(), checkList.begin());
            
            std::reverse(checkList.begin(), checkList.end());
            
            reco::PrincipalComponents checkPca = pca;
            
            if (getHitsAtEnd(checkList, checkPca))
            {
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
                
                foundGoodSeed = true;
            }
        }
    }
    
    return foundGoodSeed;
}
    
bool PCASeedFinderAlg::getHitsAtEnd(reco::HitPairListPtr& hit3DList, reco::PrincipalComponents& seedPca) const
{
    bool foundGoodSeedHits(false);
    
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
        
        // On input, the seedPca will contain the original values so we can recover the original axis now
        TVector3 planeVec0(seedPca.getEigenVectors()[0][0],seedPca.getEigenVectors()[0][1],seedPca.getEigenVectors()[0][2]);
        
        m_pcaAlg.PCAAnalysis_3D(hit3DList, seedPca, true);
        
        if (seedPca.getSvdOK())
        {
            // Still looking to point "down"
            if (seedPca.getEigenVectors()[0][1] > 0.) seedPca.flipAxis(0);
            
            // Check that the seed PCA we have found is consistent with the input PCA
            TVector3 primarySeedAxis(seedPca.getEigenVectors()[0][0],seedPca.getEigenVectors()[0][1],seedPca.getEigenVectors()[0][2]);
            
            double cosAng = primarySeedAxis.Dot(planeVec0);
            
            // If the proposed seed axis is not relatively aligned with the input PCA then
            // we should not be using this method to return seeds. Check that here
            if (cosAng > m_minAllowedCosAng) foundGoodSeedHits = true;
        }
    }
    
    return foundGoodSeedHits;
}
    

} // namespace lar_cluster3d

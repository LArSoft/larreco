/**
 *  @file   kdTree.cxx
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/kdTree.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <unordered_map>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

kdTree::kdTree(fhicl::ParameterSet const &pset)
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

kdTree::~kdTree()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void kdTree::configure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring  = pset.get<bool>  ("EnableMonitoring",  true);
    m_pairSigmaPeakTime = pset.get<double>("PairSigmaPeakTime", 3.  );

    m_timeToBuild = 0;
    
    return;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
kdTree::KdTreeNode kdTree::BuildKdTree(const reco::HitPairList& hitPairList,
                                       KdTreeNodeList&          kdTreeNodeContainer) const
{
    
    // The first task is to build the kd tree
    cet::cpu_timer theClockBuildNeighborhood;
    
    if (m_enableMonitoring) theClockBuildNeighborhood.start();
    
    // The input is a list and we need to copy to a vector so we can sort ranges
    Hit3DVec hit3DVec;
    
    hit3DVec.reserve(hitPairList.size());

    for(const auto& hitPtr : hitPairList) hit3DVec.emplace_back(hitPtr.get());
    
    KdTreeNode topNode = BuildKdTree(hit3DVec.begin(), hit3DVec.end(), kdTreeNodeContainer);
    
    if (m_enableMonitoring)
    {
        theClockBuildNeighborhood.stop();
        m_timeToBuild = theClockBuildNeighborhood.accumulated_real_time();
    }
    
    return topNode;
}
    
kdTree::KdTreeNode& kdTree::BuildKdTree(Hit3DVec::iterator first,
                                                        Hit3DVec::iterator last,
                                                        KdTreeNodeList&    kdTreeNodeContainer,
                                                        int                depth) const
{
    // Ok, so if the input list is more than one element then we have work to do... but if less then handle end condition
    if (std::distance(first,last) < 2)
    {
        if (first != last) kdTreeNodeContainer.emplace_back(KdTreeNode(*first));
        else               kdTreeNodeContainer.emplace_back(KdTreeNode());
//        if (first == last) std::cout << "********************************************* BAD NODE ***************************" << std::endl;
    }
    // Otherwise we need to keep splitting...
    else
    {
        // First task is to find "d" with the largest range. We need to find the min/max for the four dimensions
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxXPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[0] < right->getPosition()[0];});
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxYPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[1] < right->getPosition()[1];});
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxZPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[2] < right->getPosition()[2];});
        
        std::vector<double> rangeVec(3,0.);
        
        rangeVec[0] = (*minMaxXPair.second)->getPosition()[0] - (*minMaxXPair.first)->getPosition()[0];
        rangeVec[1] = (*minMaxYPair.second)->getPosition()[1] - (*minMaxYPair.first)->getPosition()[1];
        rangeVec[2] = (*minMaxZPair.second)->getPosition()[2] - (*minMaxZPair.first)->getPosition()[2];
        
        std::vector<double>::iterator maxRangeItr = std::max_element(rangeVec.begin(),rangeVec.end());
        
        size_t maxRangeIdx = std::distance(rangeVec.begin(),maxRangeItr);
        
        // Sort the list so we can do the split
        std::sort(first,last,[maxRangeIdx](const auto& left, const auto& right){return left->getPosition()[maxRangeIdx] < right->getPosition()[maxRangeIdx];});
        
        size_t             middleElem = std::distance(first,last) / 2;
        Hit3DVec::iterator middleItr  = first;
        
        std::advance(middleItr, middleElem);
        
        // Take care of the special case where the value of the median may be repeated so we actually want to make sure we point at the first occurence
        if (std::distance(first,middleItr) > 1)
        {
            while(middleItr != first+1)
            {
                if (!((*(middleItr-1))->getPosition()[maxRangeIdx] < (*middleItr)->getPosition()[maxRangeIdx])) middleItr--;
                else break;
            }
        }
        
        KdTreeNode::SplitAxis axis[]    = {KdTreeNode::xPlane,KdTreeNode::yPlane,KdTreeNode::zPlane};
        double                axisVal   = 0.5*((*middleItr)->getPosition()[maxRangeIdx] + (*(middleItr-1))->getPosition()[maxRangeIdx]);
        KdTreeNode&           leftNode  = BuildKdTree(first,     middleItr, kdTreeNodeContainer, depth+1);
        KdTreeNode&           rightNode = BuildKdTree(middleItr, last,      kdTreeNodeContainer, depth+1);
    
        kdTreeNodeContainer.push_back(KdTreeNode(axis[maxRangeIdx],axisVal,leftNode,rightNode));
    }
    
    return kdTreeNodeContainer.back();
}
    
size_t kdTree::FindNearestNeighbors(const reco::ClusterHit3D* refHit, const KdTreeNode& node, CandPairVec& candPairVec, double& bestDist) const
{
    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        double hitSeparation(0.);
        int    wireDeltas[] = {0,0,0};
        
//        std::cout << "###>> nearest neighbor, refHit wires: " << refHit->getWireIDs()[0].Wire << "/" << refHit->getWireIDs()[1].Wire << "/" << refHit->getWireIDs()[2].Wire << ", compare to: " << node.getClusterHit3D()->getWireIDs()[0].Wire << "/" << node.getClusterHit3D()->getWireIDs()[1].Wire << "/" << node.getClusterHit3D()->getWireIDs()[2].Wire << std::endl;
        
        // Is this the droid we are looking for?
        if (refHit == node.getClusterHit3D()) bestDist = 0.5;  // This distance will grab neighbors with delta wire # = 1 in all three planes
        // This is the tight constraint on the hits
        else if (bestDist < std::numeric_limits<double>::max() && consistentPairs(refHit, node.getClusterHit3D(), hitSeparation, wireDeltas))
        {
            candPairVec.emplace_back(CandPair(hitSeparation,node.getClusterHit3D()));
            
            //bestDist = std::max(0.35,std::min(bestDist,hitSeparation));  // This insures we will always consider neighbors with wire # changing in 2 planes
            //bestDist = std::max(0.47,std::min(bestDist,hitSeparation));  // This insures we will always consider neighbors with wire # changing in 2 planes
            bestDist = std::max(0.85,std::min(bestDist,hitSeparation));  // This insures we will always consider neighbors with wire # changing in 2 planes
            
//            std::cout << "###>> nearest neighbor, refHit wires: " << refHit->getWireIDs()[0].Wire << "/" << refHit->getWireIDs()[1].Wire << "/" << refHit->getWireIDs()[2].Wire << ", compare to: " << node.getClusterHit3D()->getWireIDs()[0].Wire << "/" << node.getClusterHit3D()->getWireIDs()[1].Wire << "/" << node.getClusterHit3D()->getWireIDs()[2].Wire << std::endl;
            
//            std::cout << "  ~~~> cand " << candPairVec.size() << ", wire delta u: " << wireDeltas[0] << ", v: " << wireDeltas[1] << ", w: " << wireDeltas[2] << ", sep: " << hitSeparation << ", bestDist: " << bestDist << std::endl;
        }
    }
    // Otherwise we need to keep searching
    else
    {
        double refPosition = refHit->getPosition()[node.getSplitAxis()];
        
        if (refPosition < node.getAxisValue())
        {
            FindNearestNeighbors(refHit, node.leftTree(), candPairVec, bestDist);
            
            if (refPosition + bestDist > node.getAxisValue()) FindNearestNeighbors(refHit, node.rightTree(), candPairVec, bestDist);
        }
        else
        {
            FindNearestNeighbors(refHit, node.rightTree(), candPairVec, bestDist);
            
            if (refPosition - bestDist < node.getAxisValue()) FindNearestNeighbors(refHit, node.leftTree(), candPairVec, bestDist);
        }
    }
    
    return candPairVec.size();
}
    
bool kdTree::FindEntry(const reco::ClusterHit3D* refHit, const KdTreeNode& node, CandPairVec& candPairVec, double& bestDist, bool& selfNotFound, int depth) const
{
    bool foundEntry(false);
    
    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        double hitSeparation(0.);
        int    wireDeltas[] = {0,0,0};
        
        // Is this the droid we are looking for?
        if (refHit == node.getClusterHit3D()) selfNotFound = false;
        
        // This is the tight constraint on the hits
        if (consistentPairs(refHit, node.getClusterHit3D(), hitSeparation, wireDeltas))
        {
            candPairVec.emplace_back(CandPair(hitSeparation,node.getClusterHit3D()));
            
            if (bestDist < std::numeric_limits<double>::max()) bestDist = std::max(bestDist,hitSeparation);
            else                                               bestDist = std::max(0.5,hitSeparation);
        }
        
        foundEntry = !selfNotFound;
    }
    // Otherwise we need to keep searching
    else
    {
        double refPosition = refHit->getPosition()[node.getSplitAxis()];
        
        if (refPosition < node.getAxisValue())
        {
            foundEntry = FindEntry(refHit, node.leftTree(),  candPairVec, bestDist, selfNotFound, depth+1);
            
            if (!foundEntry && refPosition + bestDist > node.getAxisValue()) foundEntry = FindEntry(refHit, node.rightTree(),  candPairVec, bestDist, selfNotFound, depth+1);
        }
        else
        {
            foundEntry = FindEntry(refHit, node.rightTree(),  candPairVec, bestDist, selfNotFound, depth+1);
            
            if (!foundEntry && refPosition - bestDist < node.getAxisValue()) foundEntry = FindEntry(refHit, node.leftTree(),  candPairVec, bestDist, selfNotFound, depth+1);
        }
    }
    
    return foundEntry;
}
    
bool kdTree::FindEntryBrute(const reco::ClusterHit3D* refHit, const KdTreeNode& node, int depth) const
{
    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        // This is the tight constraint on the hits
        if (refHit == node.getClusterHit3D()) return true;
    }
    // Otherwise we need to keep searching
    else
    {
        if (FindEntryBrute(refHit, node.leftTree(),  depth+1)) return true;
        if (FindEntryBrute(refHit, node.rightTree(), depth+1)) return true;
    }
    
    return false;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

bool kdTree::consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, double& hitSeparation, int* wireDeltas) const
{
    // Strategy: We consider comparing "hit pairs" which may consist of 2 or 3 actual hits.
    //           Also, if only pairs, they can be U-V, U-W or V-W so we can't assume which views we have
    //           So do a simplified comparison:
    //           1) compare the pair times and require "overlap" (in the sense of hit pair production)
    //           2) look at distance between pairs in each of the wire directions
    
    double pair1PeakTime = pair1->getAvePeakTime();
    double pair1Width    = m_pairSigmaPeakTime * pair1->getSigmaPeakTime();
    double pair2PeakTime = pair2->getAvePeakTime();
    double pair2Width    = m_pairSigmaPeakTime * pair2->getSigmaPeakTime();
    
    double maxUpper      = std::min(pair1PeakTime+pair1Width,pair2PeakTime+pair2Width);
    double minLower      = std::max(pair1PeakTime-pair1Width,pair2PeakTime-pair2Width);
    double pairOverlap   = maxUpper - minLower;
    
    // Loose constraint to weed out the obviously bad combinations
    if (pairOverlap > 0.0)
    {
        hitSeparation = DistanceBetweenNodes(pair1,pair2);
        
        //        size_t hitCount(0);
        
        // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
        for(size_t idx = 0; idx < 3; idx++)
        {
            wireDeltas[idx] = std::abs(int(pair1->getWireIDs()[idx].Wire) - int(pair2->getWireIDs()[idx].Wire));
            
            //            if (pair1->getHits()[idx]) hitCount++;
            //            if (pair2->getHits()[idx]) hitCount++;
        }
        
        // put wire deltas in order...
        std::sort(wireDeltas, wireDeltas + 3);
        
        // Requirement to be considered a nearest neighbor
        if (wireDeltas[0] < 2 && wireDeltas[1] < 2 && wireDeltas[2] < 3)
        {
            //            double overlapFraction = 0.5 * pairOverlap / std::min(pair1Width,pair2Width);
            
            //            hitSeparation /= overlapFraction;
            
            // Scale the hit separation by the number of missing wires
            //            for(size_t idx = 0; idx < 6-hitCount; idx++) hitSeparation *= 2.0; //1.1;
            
            //            if (wireDeltas[0] == 0) hitSeparation *= 2.0;
            
            hitSeparation = std::max(0.0001,hitSeparation);
            
            return true;
        }
    }
    
    return false;
}
    
double kdTree::DistanceBetweenNodes(const reco::ClusterHit3D* node1,const reco::ClusterHit3D* node2) const
{
    const double* node1Pos    = node1->getPosition();
    const double* node2Pos    = node2->getPosition();
    double        deltaNode[] = {node1Pos[0]-node2Pos[0], node1Pos[1]-node2Pos[1], node1Pos[2]-node2Pos[2]};
    
    // Standard euclidean distance
    return std::sqrt(deltaNode[0]*deltaNode[0]+deltaNode[1]*deltaNode[1]+deltaNode[2]*deltaNode[2]);
    
    // Manhatten distance
    //return std::fabs(deltaNode[0]) + std::fabs(deltaNode[1]) + std::fabs(deltaNode[2]);
    /*
     // Chebyshev distance
     // We look for maximum distance by wires...
     
     // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
     int wireDeltas[] = {0,0,0};
     
     for(size_t idx = 0; idx < 3; idx++)
     wireDeltas[idx] = std::abs(int(node1->getWireIDs()[idx].Wire) - int(node2->getWireIDs()[idx].Wire));
     
     // put wire deltas in order...
     std::sort(wireDeltas, wireDeltas + 3);
     
     return std::sqrt(deltaNode[0]*deltaNode[0] + 0.09 * double(wireDeltas[2]*wireDeltas[2]));
     */
}

} // namespace lar_cluster3d

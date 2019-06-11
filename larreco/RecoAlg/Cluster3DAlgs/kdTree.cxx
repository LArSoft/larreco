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
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// std includes
#include <functional>
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
    fEnableMonitoring  = pset.get<bool> ("EnableMonitoring",  true);
    fPairSigmaPeakTime = pset.get<float>("PairSigmaPeakTime", 3.  );
    fRefLeafBestDist   = pset.get<float>("RefLeafBestDist",   0.5 );
    fMaxWireDeltas     = pset.get<int>  ("MaxWireDeltas",     3   );

    fTimeToBuild = 0;

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------
kdTree::KdTreeNode kdTree::BuildKdTree(const reco::HitPairList& hitPairList,
                                       KdTreeNodeList&          kdTreeNodeContainer) const
{
    // The first task is to build the kd tree
    cet::cpu_timer theClockBuildNeighborhood;

    if (fEnableMonitoring) theClockBuildNeighborhood.start();

    // The input is a list and we need to copy to a vector so we can sort ranges
    Hit3DVec hit3DVec;

    hit3DVec.reserve(hitPairList.size());

    for(const auto& hit : hitPairList) hit3DVec.emplace_back(&hit);

    KdTreeNode topNode = BuildKdTree(hit3DVec.begin(), hit3DVec.end(), kdTreeNodeContainer);

    if (fEnableMonitoring)
    {
        theClockBuildNeighborhood.stop();
        fTimeToBuild = theClockBuildNeighborhood.accumulated_real_time();
    }

    return topNode;
}

//------------------------------------------------------------------------------------------------------------------------------------------
kdTree::KdTreeNode kdTree::BuildKdTree(const reco::HitPairListPtr& hitPairList,
                                       KdTreeNodeList&             kdTreeNodeContainer) const
{

    // The first task is to build the kd tree
    cet::cpu_timer theClockBuildNeighborhood;

    if (fEnableMonitoring) theClockBuildNeighborhood.start();

    // The input is a list and we need to copy to a vector so we can sort ranges
    //Hit3DVec hit3DVec{std::begin(hitPairList),std::end(hitPairList)};
    Hit3DVec hit3DVec;

    hit3DVec.reserve(hitPairList.size());

    for(const auto& hit3D : hitPairList)
    {
        // Make sure all the bits used by the clustering stage have been cleared
        hit3D->clearStatusBits(~(reco::ClusterHit3D::HITINVIEW0 | reco::ClusterHit3D::HITINVIEW1 | reco::ClusterHit3D::HITINVIEW2));
        for(const auto& hit2D : hit3D->getHits())
            if (hit2D) hit2D->clearStatusBits(0xFFFFFFFF);
        hit3DVec.emplace_back(hit3D);
    }

    KdTreeNode topNode = BuildKdTree(hit3DVec.begin(), hit3DVec.end(), kdTreeNodeContainer);

    if (fEnableMonitoring)
    {
        theClockBuildNeighborhood.stop();
        fTimeToBuild = theClockBuildNeighborhood.accumulated_real_time();
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
        if (first != last) kdTreeNodeContainer.emplace_back(*first);
        else               kdTreeNodeContainer.emplace_back(KdTreeNode());
    }
    // Otherwise we need to keep splitting...
    else
    {
        // First task is to find "d" with the largest range. We need to find the min/max for the four dimensions
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxXPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[0] < right->getPosition()[0];});
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxYPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[1] < right->getPosition()[1];});
        std::pair<Hit3DVec::iterator,Hit3DVec::iterator> minMaxZPair = std::minmax_element(first,last,[](const reco::ClusterHit3D* left, const reco::ClusterHit3D* right){return left->getPosition()[2] < right->getPosition()[2];});

        std::vector<float> rangeVec(3,0.);

        rangeVec[0] = (*minMaxXPair.second)->getPosition()[0] - (*minMaxXPair.first)->getPosition()[0];
        rangeVec[1] = (*minMaxYPair.second)->getPosition()[1] - (*minMaxYPair.first)->getPosition()[1];
        rangeVec[2] = (*minMaxZPair.second)->getPosition()[2] - (*minMaxZPair.first)->getPosition()[2];

        std::vector<float>::iterator maxRangeItr = std::max_element(rangeVec.begin(),rangeVec.end());

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
        float                 axisVal   = 0.5*((*middleItr)->getPosition()[maxRangeIdx] + (*(middleItr-1))->getPosition()[maxRangeIdx]);
        KdTreeNode&           leftNode  = BuildKdTree(first,     middleItr, kdTreeNodeContainer, depth+1);
        KdTreeNode&           rightNode = BuildKdTree(middleItr, last,      kdTreeNodeContainer, depth+1);

        kdTreeNodeContainer.push_back(KdTreeNode(axis[maxRangeIdx],axisVal,leftNode,rightNode));
    }

    return kdTreeNodeContainer.back();
}

size_t kdTree::FindNearestNeighbors(const reco::ClusterHit3D* refHit, const KdTreeNode& node, CandPairList& CandPairList, float& bestDist) const
{
    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        // Is this the droid we are looking for?
        if (refHit == node.getClusterHit3D()) bestDist = fRefLeafBestDist;  // This distance will grab neighbors with delta wire # = 1 in all three planes
        // This is the tight constraint on the hits
        else if (consistentPairs(refHit, node.getClusterHit3D(), bestDist))
        {
            CandPairList.emplace_back(bestDist,node.getClusterHit3D());

            bestDist = std::max(fRefLeafBestDist, bestDist);  // This insures we will always consider neighbors with wire # changing in 2 planes
        }
    }
    // Otherwise we need to keep searching
    else
    {
        float refPosition = refHit->getPosition()[node.getSplitAxis()];

        if (refPosition < node.getAxisValue())
        {
            FindNearestNeighbors(refHit, node.leftTree(), CandPairList, bestDist);

            if (refPosition + bestDist > node.getAxisValue()) FindNearestNeighbors(refHit, node.rightTree(), CandPairList, bestDist);
        }
        else
        {
            FindNearestNeighbors(refHit, node.rightTree(), CandPairList, bestDist);

            if (refPosition - bestDist < node.getAxisValue()) FindNearestNeighbors(refHit, node.leftTree(), CandPairList, bestDist);
        }
    }

    return CandPairList.size();
}

bool kdTree::FindEntry(const reco::ClusterHit3D* refHit, const KdTreeNode& node, CandPairList& CandPairList, float& bestDist, bool& selfNotFound, int depth) const
{
    bool foundEntry(false);

    // If at a leaf then time to decide to add hit or not
    if (node.isLeafNode())
    {
        float hitSeparation(0.);

        // Is this the droid we are looking for?
        if (refHit == node.getClusterHit3D()) selfNotFound = false;

        // This is the tight constraint on the hits
        if (consistentPairs(refHit, node.getClusterHit3D(), hitSeparation))
        {
            CandPairList.emplace_back(hitSeparation,node.getClusterHit3D());

            if (bestDist < std::numeric_limits<float>::max()) bestDist = std::max(bestDist,hitSeparation);
            else                                              bestDist = std::max(float(0.5),hitSeparation);
        }

        foundEntry = !selfNotFound;
    }
    // Otherwise we need to keep searching
    else
    {
        float refPosition = refHit->getPosition()[node.getSplitAxis()];

        if (refPosition < node.getAxisValue())
        {
            foundEntry = FindEntry(refHit, node.leftTree(),  CandPairList, bestDist, selfNotFound, depth+1);

            if (!foundEntry && refPosition + bestDist > node.getAxisValue()) foundEntry = FindEntry(refHit, node.rightTree(),  CandPairList, bestDist, selfNotFound, depth+1);
        }
        else
        {
            foundEntry = FindEntry(refHit, node.rightTree(),  CandPairList, bestDist, selfNotFound, depth+1);

            if (!foundEntry && refPosition - bestDist < node.getAxisValue()) foundEntry = FindEntry(refHit, node.leftTree(),  CandPairList, bestDist, selfNotFound, depth+1);
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

bool kdTree::consistentPairs(const reco::ClusterHit3D* pair1, const reco::ClusterHit3D* pair2, float& bestDist) const
{
    // Strategy: We consider comparing "hit pairs" which may consist of 2 or 3 actual hits.
    //           Also, if only pairs, they can be U-V, U-W or V-W so we can't assume which views we have
    //           So do a simplified comparison:
    //           1) compare the pair times and require "overlap" (in the sense of hit pair production)
    //           2) look at distance between pairs in each of the wire directions

    bool consistent(false);

    if (bestDist < std::numeric_limits<float>::max() && pair1->getWireIDs()[0].Cryostat == pair2->getWireIDs()[0].Cryostat && pair1->getWireIDs()[0].TPC == pair2->getWireIDs()[0].TPC)
    {
        // Loose constraint to weed out the obviously bad combinations
        // So this is not strictly correct but is close enough and should save computation time...
        if (std::fabs(pair1->getAvePeakTime() - pair2->getAvePeakTime()) < fPairSigmaPeakTime * (pair1->getSigmaPeakTime() + pair2->getSigmaPeakTime()))
        {
            int wireDeltas[] = {0, 0, 0};

            // Now go through the hits and compare view by view to look for delta wire and tigher constraint on delta t
            for(size_t idx = 0; idx < 3; idx++)
                wireDeltas[idx] = std::abs(int(pair1->getWireIDs()[idx].Wire) - int(pair2->getWireIDs()[idx].Wire));

            // put wire deltas in order...
            std::sort(wireDeltas, wireDeltas + 3);

            bool checkSeparation(false);

            // Because we hve sorted this is all we need to check
            if (wireDeltas[2] < 2) checkSeparation = true;

            // Requirement to be considered a nearest neighbor
            if (checkSeparation)
            {
                float hitSeparation = std::max(float(0.0001),DistanceBetweenNodesYZ(pair1,pair2));

                // Final cut...
                if (hitSeparation < bestDist)
                {
                    bestDist = hitSeparation;
                    consistent = true;
                }
            }
        }
    }

    return consistent;
}


float kdTree::DistanceBetweenNodesYZ(const reco::ClusterHit3D* node1,const reco::ClusterHit3D* node2) const
{
    const Eigen::Vector3f& node1Pos    = node1->getPosition();
    const Eigen::Vector3f& node2Pos    = node2->getPosition();
    float                  deltaNode[] = {node1Pos[0]-node2Pos[0], node1Pos[1]-node2Pos[1], node1Pos[2]-node2Pos[2]};
    float                  yzDist2     = deltaNode[1]*deltaNode[1] + deltaNode[2]*deltaNode[2];

    // Standard euclidean distance
    return std::sqrt(yzDist2);
}

float kdTree::DistanceBetweenNodes(const reco::ClusterHit3D* node1,const reco::ClusterHit3D* node2) const
{
    const Eigen::Vector3f& node1Pos    = node1->getPosition();
    const Eigen::Vector3f& node2Pos    = node2->getPosition();
    float                  deltaNode[] = {node1Pos[0]-node2Pos[0], node1Pos[1]-node2Pos[1], node1Pos[2]-node2Pos[2]};
    float                  yzDist2     = deltaNode[1]*deltaNode[1] + deltaNode[2]*deltaNode[2];
    float                  xDist2      = deltaNode[0]*deltaNode[0];

    // Standard euclidean distance
    return std::sqrt(xDist2 + yzDist2);

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

     return std::sqrt(deltaNode[0]*deltaNode[0] + 0.09 * float(wireDeltas[2]*wireDeltas[2]));
     */
}

} // namespace lar_cluster3d

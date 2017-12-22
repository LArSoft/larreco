/**
 *  @file   BeachLine.cxx
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes

#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/BeachLine.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

BSTNode::BSTNode(IEvent* event, BSTNode* parent, BSTNode* leftChild, BSTNode* rightChild)
{
    m_event       = event;
    m_parent      = parent;
    m_leftChild   = leftChild;
    m_rightChild  = rightChild;
    m_predecessor = NULL;
    m_successor   = NULL;
}

BeachLine::BeachLine(int numSites) : m_root(0), m_numSites(numSites)
{
    m_nodeVec.reserve(4*m_numSites);
}
    
BSTNode* BeachLine::findBestLeaf(IEvent* event, BSTNode* topNode) const
{
    // Assumption: a leaf will have NULL child pointers so the idea is to
    // follow the left or right child pointers until we get to a leaf
    BSTNode* node = topNode;
    
    while(node->getLeftChild() && node->getRightChild())
    {
        // Which path do we follow down the tree?
        if (node->getEvent()->newSiteToLeft(event, node->getPredecessor()->getEvent(), node->getSuccessor()->getEvent()))
                node = findBestLeaf(event, node->getLeftChild());
        else    node = findBestLeaf(event, node->getRightChild());
    }
    
    return node;
}

BSTNode* BeachLine::insertNewLeaf(IEvent* event, BSTNode* node)
{
    // Have we found a null pointer?
    if (node == NULL)
    {
        m_nodeVec.push_back(BSTNode(event));
        node = &m_nodeVec.back();
        return node;
    }
    
    // Right now assume that the input node is the leaf at which we want to insert the new data
    // For the beachline, this means we are inserting a new arc into the beachline by dividing the
    // current arc. So we are going to replace the input leaf with a subtree having three leaves
    // (two breakpoints)...
    // Start by creating a node for the new arc
    m_nodeVec.push_back(BSTNode(event));  // This will be the new site point
    
    BSTNode* newLeaf = &m_nodeVec.back();
    
    m_nodeVec.push_back(*node);  // This will be the new left leaf (the original arc)
    
    BSTNode* leftLeaf = &m_nodeVec.back();
    
    m_nodeVec.push_back(BSTNode(event));  // This will be the breakpoint between the left and new leaves
    
    BSTNode* breakNode = &m_nodeVec.back();
    
    m_nodeVec.push_back(BSTNode(event)); // Finally, this is the breakpoint between new and right leaves
    
    BSTNode* topNode = &m_nodeVec.back();
    
    // Set this to be the king of the local world
    topNode->setParent(node->getParent());

    // Make sure the parent points to this new node
    if (node->getParent())
    {
        if (node->getParent()->getLeftChild() == node) node->getParent()->setLeftChild(topNode);
        else                                           node->getParent()->setRightChild(topNode);
    }

    // Set our children
    topNode->setLeftChild(breakNode);
    topNode->setRightChild(node);

    // Original node is now child of the top node
    node->setParent(topNode);

    // Now set the children of the left breakpoint node
    breakNode->setLeftChild(leftLeaf);
    breakNode->setRightChild(newLeaf);
    
    // Finally set the leaves parents
    leftLeaf->setParent(breakNode);
    newLeaf->setParent(newLeaf);
    
    // Now we wire up traversal chain, going left to right
    leftLeaf->setPredecessor(node->getPredecessor());
    leftLeaf->setSuccessor(breakNode);
    
    breakNode->setPredecessor(leftLeaf);
    breakNode->setSuccessor(newLeaf);
    
    newLeaf->setPredecessor(breakNode);
    newLeaf->setSuccessor(topNode);
    
    node->setPredecessor(topNode);
    
    return topNode;
}

    
} // namespace lar_cluster3d

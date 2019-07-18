/**
 *  @file   BeachLine.cxx
 *
 *  @brief  Producer module to create 3D clusters from input hits
 *
 */

// Framework Includes

#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/BeachLine.h"

// std includes
#include <iostream>
#include <limits>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace voronoi2d {

BSTNode::BSTNode(IEvent* event, BSTNode* parent, BSTNode* leftChild, BSTNode* rightChild)
{
    m_depth       = 0;
    m_event       = event;
    m_parent      = parent;
    m_leftChild   = leftChild;
    m_rightChild  = rightChild;
    m_predecessor = NULL;
    m_successor   = NULL;
    m_associated  = NULL;

    m_event->setBSTNode(this);

    // Reset depth
    setDepth();
}

void BSTNode::setDepth()
{
    if (m_leftChild && m_rightChild)
    {
        int maxDepth = std::max(m_leftChild->getDepth(),m_rightChild->getDepth());

        m_depth = maxDepth + 1;
    }
    else m_depth = 0;

    // If we change depth at this level then need to ripple it up through the tree
    if (m_parent) m_parent->setDepth();

    return;
}

BSTNode* BeachLine::insertNewLeaf(IEvent* event)
{
    // Find the insertion point for the new event
    BSTNode* node = findBestLeaf(event, m_root);

    // Insert it
    node = insertNewLeaf(event, node);

    // Now rebalance starting at this node
    rebalance(node);

    // Check beach line integrity
//    checkBeachLine(event->xPos() - 0.000001);

    return node;
}

BSTNode* BeachLine::findBestLeaf(const IEvent* event, BSTNode* topNode) const
{
    // Assumption: a leaf will have NULL child pointers so the idea is to
    // follow the left or right child pointers until we get to a leaf
    BSTNode* node = topNode;

    // A leaf is found when the node has no children
    if (node && node->getLeftChild() && node->getRightChild())
    {
        // This node represents a breakpoint between two arcs and we can
        // recover these immediately by getting the predecessor and successor leafs
        BSTNode* rightLeaf = node->getSuccessor();
        BSTNode* leftLeaf  = node->getPredecessor();

        // Which path do we follow down the tree?
        if (m_utilities.newSiteToLeft(event, leftLeaf->getEvent(), rightLeaf->getEvent()))
            node = findBestLeaf(event, node->getLeftChild());
        else
            node = findBestLeaf(event, node->getRightChild());
    }

    return node;
}

BSTNode* BeachLine::insertNewLeaf(IEvent* event, BSTNode* node)
{
    // The idea of this function is to insert a new Site Event into the beach line
    // where it is assumed that the input node is the matched arc into which we
    // insert the new site event.
    // The function then returns the new leaf created which represents the new arc

    // Have we found a null pointer?
    if (node == NULL)
    {
        m_nodeVec.push_back(BSTNode(event));
        node = &m_nodeVec.back();
        return node;
    }

    // Check if a circle event had been definied for the node we are splitting,
    // if so then we need to invalidate that circle event
    if (node->getAssociated())
    {
        node->getAssociated()->getEvent()->setInvalid();
        node->getAssociated()->setAssociated(NULL);
        node->setAssociated(NULL);
    }

    // Right now assume that the input node is the leaf at which we want to insert the new data
    // For the beachline, this means we are inserting a new arc into the beachline by dividing the
    // current arc. So we are going to replace the input leaf with a subtree having three leaves
    // (two breakpoints)...
    // Start by creating a node for the new arc
    m_nodeVec.push_back(BSTNode(event));  // This will be the new site point

    BSTNode* newLeaf = &m_nodeVec.back();

    m_nodeVec.push_back(BSTNode(*node));  // This will be the new left leaf (the original arc)

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
    // But if there is no parent then the topnode is the new root
    else m_root = topNode;

    // Set our children
    topNode->setLeftChild(breakNode);
    topNode->setRightChild(node);

    // Original node is now child of the top node
    node->setParent(topNode);

    // If there was an associated circle event to this node, invalidate it
    if (node->getAssociated())
    {
        node->getAssociated()->setAssociated(NULL);
        node->getAssociated()->getEvent()->setInvalid();
        node->setAssociated(NULL);
        leftLeaf->setAssociated(NULL);
    }

    // Now set the parent and children of the left breakpoint node
    breakNode->setParent(topNode);
    breakNode->setLeftChild(leftLeaf);
    breakNode->setRightChild(newLeaf);

    // Finally set the leaves parents
    leftLeaf->setParent(breakNode);
    newLeaf->setParent(breakNode);

    // Now we wire up traversal chain, going left to right
    leftLeaf->setPredecessor(node->getPredecessor());
    leftLeaf->setSuccessor(breakNode);

    if (node->getPredecessor()) node->getPredecessor()->setSuccessor(leftLeaf);

    breakNode->setPredecessor(leftLeaf);
    breakNode->setSuccessor(newLeaf);

    newLeaf->setPredecessor(breakNode);
    newLeaf->setSuccessor(topNode);

    topNode->setPredecessor(newLeaf);
    topNode->setSuccessor(node);
    node->setPredecessor(topNode);

    // Finally, reset the depths
    // By definition the break node will have depth of 1
    breakNode->setDepth();

    // Check the beachline integrity
    checkBeachLine(newLeaf->getEvent()->xPos()-0.00000001);

    return newLeaf;
}

BSTNode* BeachLine::removeLeaf(BSTNode* node)
{
    // The input node is assumed to be a leaf (arc) and is the disappearing arc
    // between a leaf (arc) to the left and one to the right. There are breakpoints
    // between the arcs.
    // One of the intervening breakpoints is the parent of the leaf to remove
    BSTNode* nodeParent = node->getParent();

    // parent of the parent
    BSTNode* grandParent = nodeParent->getParent();

    // Parent is either left or right child, node is left or right child.... my brain is dizzy
    BSTNode* sibling = nodeParent->getRightChild();

    if (node == sibling) sibling = nodeParent->getLeftChild();

    if (nodeParent == grandParent->getLeftChild()) grandParent->setLeftChild(sibling);
    else                                           grandParent->setRightChild(sibling);

    sibling->setParent(grandParent);

    // Now we need to deal with the predecessor/successor chain.
    // This should be straightforward, we are removing the middle arc and the immediate parent
    // breakpoint with the grandparent becoming the new breakpoint. It should be that we simply
    // set the left arc's successor, the right arc's predecessor and make sure the grandparent points
    // to the right objects.
    // Also note that if here there MUST be a left and right arc
    BSTNode* arcLeft  = node->getPredecessor()->getPredecessor();
    BSTNode* arcRight = node->getSuccessor()->getSuccessor();

    // Note as well that any circle events for those arcs are now invalid
    if (arcLeft->getAssociated())
    {
        arcLeft->getAssociated()->getEvent()->setInvalid();
        arcLeft->getAssociated()->setAssociated(NULL);
        arcLeft->setAssociated(NULL);
    }

    if (arcRight->getAssociated())
    {
        arcRight->getAssociated()->getEvent()->setInvalid();
        arcRight->getAssociated()->setAssociated(NULL);
        arcRight->setAssociated(NULL);
    }

    // Basically, we need to connect the left and right arcs to their common break point
    // What breakpoint that is will be determined by weather the arc we are removing is a
    // left or right child
    if (node == nodeParent->getLeftChild())
    {
        // In this case, the right arc's predecessor becomes the node's predecessor
        // The left arc is all set
        arcRight->setPredecessor(node->getPredecessor());
        node->getPredecessor()->setSuccessor(arcRight);
    }
    else
    {
        // Here the left arc's successor is what needs to be changed
        arcLeft->setSuccessor(node->getSuccessor());
        node->getSuccessor()->setPredecessor(arcLeft);
    }

    // zap the pointers for the removed nodes
    node->setParent(NULL);
    nodeParent->setLeftChild(NULL);
    nodeParent->setRightChild(NULL);
    nodeParent->setParent(NULL);
    nodeParent->setSuccessor(NULL);
    nodeParent->setPredecessor(NULL);
    node->setSuccessor(NULL);
    node->setPredecessor(NULL);
    node->setHalfEdge(NULL);
    node->setFace(NULL);

    // Reset the depth
    grandParent->setDepth();

    // Check beach line integrity
//    checkBeachLine(node->getAssociated()->getEvent()->xPos()-0.000001);

    // Rebalance
    rebalance(grandParent);

    // Check beach line integrity
    checkBeachLine(node->getAssociated()->getEvent()->xPos()-0.00000001);

    // Return the new breakpoint
    return arcLeft->getSuccessor();
}

int BeachLine::countNodes() const
{
    int nodeCount(0);

    countNodes(m_root, nodeCount);

    return nodeCount;
}

int BeachLine::countLeaves() const
{
    int leafCount(0);

    countLeaves(m_root, leafCount);

    return leafCount;
}

void BeachLine::countNodes(const BSTNode* node, int& nodeCount) const
{
    if (node)
    {
        if (node->getLeftChild())  countNodes(node->getLeftChild(),  nodeCount);
        if (node->getRightChild()) countNodes(node->getRightChild(), nodeCount);

        if ((node->getLeftChild() != NULL) != (node->getRightChild() != NULL))
        {
            std::cout << "****** Tree has one branch but not the other! *******" << std::endl;
        }

        nodeCount++;
    }

    return;
}

void BeachLine::countLeaves(const BSTNode* node, int& leafCount) const
{
    // If not a leaf then still have children to search
    if (node->getLeftChild() && node->getRightChild())
    {
        countLeaves(node->getLeftChild(),  leafCount);
        countLeaves(node->getRightChild(), leafCount);
    }
    else leafCount += 1;

    return;
}

int BeachLine::traverseBeach() const
{
    int leafCount(0);

    // Starting with the root, dive down until we find a leaf
    BSTNode* node = m_root;

    // Basically, follow the left branch down until we have no more children
    while(node->getLeftChild()) node = node->getLeftChild();

    // Note that construction we should be a leaf, now we traverse across
    // the beach line in both directions to get the leaf count
    if (node)
    {
        leafCount += traverseBeachLeft(node->getPredecessor());
        leafCount += traverseBeachRight(node->getSuccessor());

        // just to be sure...
        if (!node->getLeftChild() && !node->getRightChild()) leafCount++;
    }

    return leafCount;
}

void BeachLine::checkBeachLine(double beachLine) const
{
    // Starting with the root, dive down until we find the leftmost leaf
    BSTNode* node = m_root;

    if (!node) return;

    // Basically, follow the left branch down until we have no more children
    while(node->getLeftChild()) node = node->getLeftChild();

    // Keep track of breakpoints
    double   lastBreakPointY = -std::numeric_limits<double>::max();
    int      nBadCompares(0);
    int      nNodes(0);
    int      nBreakPoints(0);
    int      nLeaves(0);
    BSTNode* lastBreakPoint(NULL);

    const double tolerance(1.e-5);

    // This is the start of the beach line, we now traverse across and and check status
    // of each breakpoint's position
    while(node->getSuccessor())
    {
        // Is this a breakpoint?
        if (node->getLeftChild() && node->getRightChild())
        {
            RootsPair roots;
            double breakPoint = m_utilities.computeBreak(beachLine, node->getPredecessor()->getEvent(), node->getSuccessor()->getEvent(), roots);

            if (breakPoint + tolerance < lastBreakPointY)
            {
                std::cout << "  #####>> Beach line check gets bad breakpoint, last: " << lastBreakPointY << ", new: " << breakPoint << ", roots: " << roots.first << "/" << roots.second << std::endl;
                std::cout << "          Current left arc x,y: " << node->getPredecessor()->getEvent()->xPos() << ", " << node->getPredecessor()->getEvent()->yPos() << ", right arc x,y: " << node->getSuccessor()->getEvent()->xPos() << ", " << node->getSuccessor()->getEvent()->yPos() << ", beachLine: " << beachLine;
                if (node->getPredecessor()->getAssociated()) std::cout << ", left: "  << node->getPredecessor()->getAssociated()->getEvent()->isValid();
                if (node->getSuccessor()->getAssociated())   std::cout << ", right: " << node->getSuccessor()->getAssociated()->getEvent()->isValid();
                std::cout << std::endl;
                if (lastBreakPoint)
                {
                    std::cout << "          Previous left arc x,y: " << lastBreakPoint->getPredecessor()->getEvent()->xPos() << ", " << lastBreakPoint->getPredecessor()->getEvent()->yPos() << ", right arc x,y: " << lastBreakPoint->getSuccessor()->getEvent()->xPos() << ", " << lastBreakPoint->getSuccessor()->getEvent()->yPos();
                    if (lastBreakPoint->getPredecessor()->getAssociated()) std::cout << ", left: "  << lastBreakPoint->getPredecessor()->getAssociated()->getEvent()->isValid();
                    if (lastBreakPoint->getSuccessor()->getAssociated())   std::cout << ", right: " << lastBreakPoint->getSuccessor()->getAssociated()->getEvent()->isValid();
                    std::cout << std::endl;
                }
                nBadCompares++;
            }

            lastBreakPointY = breakPoint;
            lastBreakPoint  = node;
            nBreakPoints++;
        }
        else
        {
            // Confirm that the next leaf in the beachline is the next according to the tree
            if (node->getSuccessor())
            {
                BSTNode* temp = node;
                while(temp->getParent() && temp != temp->getParent()->getLeftChild()) temp = temp->getParent();
                if (temp->getParent()) temp = temp->getParent()->getRightChild();
                while(temp->getLeftChild()) temp = temp->getLeftChild();

                if (node->getSuccessor()->getSuccessor() != temp || node != temp->getPredecessor()->getPredecessor())
                {
                    std::cout << "          --> Successor tree/beach mismatch, leaf # " << nLeaves << ", node: " << node << ", " << node->getEvent()->xPos() << "/" << node->getEvent()->yPos() << ", s: " << node->getSuccessor() << ", ss: " << node->getSuccessor()->getSuccessor() << std::endl;
                    std::cout << "              temp: " << temp << ", " << temp->getEvent()->xPos() << "/" << temp->getEvent()->yPos() << ", p: " << temp->getPredecessor() << ", pp: " << temp->getPredecessor()->getPredecessor() << std::endl;
                }
            }

            if (node->getPredecessor())
            {
                BSTNode* temp = node;
                while(temp->getParent() && temp != temp->getParent()->getRightChild()) temp = temp->getParent();
                if (temp->getParent()) temp = temp->getParent()->getLeftChild();
                while(temp->getRightChild()) temp = temp->getRightChild();

                if (node->getPredecessor()->getPredecessor() != temp || node != temp->getSuccessor()->getSuccessor())
                {
                    std::cout << "          --> Predecessor tree/beach mismatch, leaf # " << nLeaves << ", node: " << node << ", " << node->getEvent()->xPos() << "/" << node->getEvent()->yPos() << ", p: " << node->getPredecessor() << ", pp: " << node->getPredecessor()->getPredecessor() << std::endl;
                    std::cout << "              temp: " << temp << ", " << temp->getEvent()->xPos() << "/" << temp->getEvent()->yPos() << ", s: " << temp->getSuccessor() << ", ss: " << temp->getSuccessor()->getSuccessor() << std::endl;
                }
            }

            nLeaves++;
        }

        nNodes++;

        node = node->getSuccessor();
    }

    if (nBadCompares > 0) std::cout << "=======>> Beach line check resulted in " << nBadCompares << " bad compares of " << nBreakPoints << " break points checked, with " << nLeaves << " leaves" << std::endl;

//    std::cout << "-------------------------------------------------------------------------------------------------------" << std::endl;

    return;
}

int BeachLine::traverseBeachLeft(BSTNode* node) const
{
    int leafCount(0);

    if (node)
    {
        // Keep traversing
        if (node->getPredecessor()) leafCount += traverseBeachLeft(node->getPredecessor());

        // Are we also a leaf?
        if (!node->getLeftChild() && !node->getRightChild()) leafCount++;
    }

    return leafCount;
}

int BeachLine::traverseBeachRight(BSTNode* node) const
{
    int leafCount(0);

    if (node)
    {
        // Keep traversing
        if (node->getSuccessor()) leafCount += traverseBeachRight(node->getSuccessor());

        // Are we also a leaf?
        if (!node->getLeftChild() && !node->getRightChild()) leafCount++;
    }

    return leafCount;
}

int BeachLine::getTreeDepth(const BSTNode* node) const
{
    int depth(0);

    // Node exists and its not a leaf
    if (node && node->getLeftChild() && node->getRightChild())
    {
        depth = std::max(getTreeDepth(node->getLeftChild()),getTreeDepth(node->getRightChild()));

        depth++;
    }
    else if (node && (node->getLeftChild() || node->getRightChild()))
    {
        std::cout << "****** Found a node which only one child: " << node << ", L/R: " << node->getLeftChild() << "/" << node->getRightChild() << std::endl;
    }

    return depth;
}

void BeachLine::rebalance(BSTNode* node)
{
    // The idea is to rebalance starting with the current node and the walking back up the branch
    // until we reach the ultimate parent.
    // First, if at internal node then check depth down either branch
    if (node->getLeftChild() && node->getRightChild())
    {
        int depthLeft  = getTreeDepth(node->getLeftChild());
        int depthRight = getTreeDepth(node->getRightChild());

        if (depthLeft != node->getLeftChild()->getDepth() || depthRight != node->getRightChild()->getDepth())
        {
            std::cout << "       --> node depth: " << getTreeDepth(node) << ", left/right: " << depthLeft << "/" << node->getLeftChild()->getDepth() << ", " << depthRight << "/" << node->getRightChild()->getDepth() << ", parent/depth " << node->getParent() << "/" << getTreeDepth(node->getParent()) << std::endl;
        }

        // If left branch is longer then we rotate with the left child
        if      (depthLeft - depthRight > 2) node = rotateWithLeftChild(node);
        else if (depthRight - depthLeft > 2) node = rotateWithRightChild(node);
    }

    // Ok now rebalance the parent unless we are the root
    if (node->getParent()) rebalance(node->getParent());
    // In which case update the internal root node pointer
    else m_root = node;

    return;
}

BSTNode* BeachLine::rotateWithLeftChild(BSTNode* node)
{
    // Here we rebalance by rotating the root node with its left child
    BSTNode* newTopNode = node->getLeftChild();
    BSTNode* parent     = node->getParent();

    // Check if there is a parent and if so make sure it points at the new node
    if (parent)
    {
        if (parent->getLeftChild() == node) parent->setLeftChild(newTopNode);
        else                                parent->setRightChild(newTopNode);
    }
    // if no parent this the new root
    else m_root = newTopNode;

    // Swap parents (for the root node the parent is null)
    newTopNode->setParent(parent);
    node->setParent(newTopNode);

    // Reset the children
    BSTNode* childToSwitch = newTopNode->getRightChild();

    childToSwitch->setParent(node);
    node->setLeftChild(childToSwitch);
    newTopNode->setRightChild(node);

    // Reset node depth
    node->getLeftChild()->setDepth();

    return newTopNode;
}

BSTNode* BeachLine::rotateWithRightChild(BSTNode* node)
{
    // Here we rebalance by rotating the root node with its left child
    BSTNode* newTopNode = node->getRightChild();
    BSTNode* parent     = node->getParent();

    // Check if there is a parent and if so make sure it points at the new node
    if (parent)
    {
        if (parent->getLeftChild() == node) parent->setLeftChild(newTopNode);
        else                                parent->setRightChild(newTopNode);
    }
    // if no parent this the new root
    else m_root = newTopNode;

    // Swap parents (for the root node the parent is null)
    newTopNode->setParent(parent);
    node->setParent(newTopNode);

    // Reset the children
    BSTNode* childToSwitch = newTopNode->getLeftChild();

    childToSwitch->setParent(node);
    node->setRightChild(childToSwitch);
    newTopNode->setLeftChild(node);

    // Reset node depths
    node->getRightChild()->setDepth();

    return newTopNode;
}


} // namespace lar_cluster3d

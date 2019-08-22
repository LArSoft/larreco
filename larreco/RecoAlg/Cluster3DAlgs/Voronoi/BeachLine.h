/**
 *  @file   BeachLine.h
 *
 *  @brief  Represents the beachline implemented as a self balancing binary search tree
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef BeachLine_h
#define BeachLine_h

#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/IEvent.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/EventUtilities.h"
namespace dcel2d { class Face; class HalfEdge; }

// std includes
#include <list>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace voronoi2d
{
/**
 *  @brief  BSTNode class definiton specifically for use in constructing
 *          Voronoi diagrams. We are trying to follow the prescription
 *          described in "Computational Geometry" by Mark de Berg, et al.
 *
 *          Note that in this implementation the internal nodes of the tree
 *          will describe the breakpoints in the beach line and the leaves of
 *          the tree will describe the arcs (site points).
 */
class BSTNode
{
public:
    /**
     *  @brief  Constructor
     */
    BSTNode() :
        m_depth(0),
        m_event(NULL),
        m_parent(NULL),
        m_leftChild(NULL),
        m_rightChild(NULL),
        m_predecessor(NULL),
        m_successor(NULL),
        m_associated(NULL),
        m_halfEdge(NULL),
        m_face(NULL)
    {}

    BSTNode(IEvent* event) :
        m_depth(0),
        m_event(event),
        m_parent(NULL),
        m_leftChild(NULL),
        m_rightChild(NULL),
        m_predecessor(NULL),
        m_successor(NULL),
        m_associated(NULL),
        m_halfEdge(NULL),
        m_face(NULL)
    {
        if (m_event) m_event->setBSTNode(this);
    }

    BSTNode(IEvent*, BSTNode*, BSTNode*, BSTNode*);

    /**
     *  @brief recover the data members
     */
    int               getDepth()       const {return m_depth;}
    IEvent*           getEvent()       const {return m_event;}
    BSTNode*          getParent()      const {return m_parent;}
    BSTNode*          getLeftChild()   const {return m_leftChild;}
    BSTNode*          getRightChild()  const {return m_rightChild;}
    BSTNode*          getPredecessor() const {return m_predecessor;}
    BSTNode*          getSuccessor()   const {return m_successor;}
    BSTNode*          getAssociated()  const {return m_associated;}

    dcel2d::HalfEdge* getHalfEdge()    const {return m_halfEdge;}
    dcel2d::Face*     getFace()        const {return m_face;}

    /**
     *  @brief Allow setting of the points
     */
    void setParent(BSTNode* node)              {m_parent      = node;}
    void setLeftChild(BSTNode* node)           {m_leftChild   = node;}
    void setRightChild(BSTNode* node)          {m_rightChild  = node;}
    void setPredecessor(BSTNode* node)         {m_predecessor = node;}
    void setSuccessor(BSTNode* node)           {m_successor   = node;}
    void setAssociated(BSTNode* node)          {m_associated  = node;}

    void setHalfEdge(dcel2d::HalfEdge* half)   {m_halfEdge    = half;}
    void setFace(dcel2d::Face* face)           {m_face        = face;}

    void setDepth(int depth)                   {m_depth       = depth;}
    void setDepth();

    /**
     *  @brief Provide override definition for ordering
     */
    bool operator<(const BSTNode&) const;

private:
    int               m_depth;        // Keep track of depth of nodes to this one
    IEvent*           m_event;        // Pointer to the event object
    BSTNode*          m_parent;       // Tree traversal - parent node
    BSTNode*          m_leftChild;    // Tree traversal - left child node
    BSTNode*          m_rightChild;   // Tree traversal - right child node
    BSTNode*          m_predecessor;  // Beachline traversal - predecessor
    BSTNode*          m_successor;    // Beachline traversal - successor
    BSTNode*          m_associated;   // This allows handling of circle events
    dcel2d::HalfEdge* m_halfEdge;     // If a breakpoint then we associate halfedges
    dcel2d::Face*     m_face;         // If a leaf then we associated faces
};

using BSTNodeList = std::list<BSTNode>;

/**
 * @brief This defines the actual beach line. The idea is to implement this as a
 *        self balancing binary search tree.
 */

class BeachLine
{
public:
    BeachLine() : m_root(NULL) {m_nodeVec.clear();}

    bool           isEmpty()                         const {return m_root == NULL;}
    void           setEmpty()                              {m_root = NULL;}
    const BSTNode* getTopNode()                      const {return m_root;}
    BSTNode*       findBestLeaf(const IEvent* event) const {return findBestLeaf(event, m_root);}
    BSTNode*       insertNewLeaf(IEvent*);
    BSTNode*       removeLeaf(BSTNode*);
    int            getHeight()                       const {return getTreeDepth(m_root);}
    int            countNodes()                      const;
    int            countLeaves()                     const;
    int            traverseBeach()                   const;

private:
    BSTNode* insertNewLeaf(IEvent*, BSTNode*);

    BSTNode* findBestLeaf(const IEvent*, BSTNode*) const;

    void countNodes(const BSTNode*, int&) const;
    void countLeaves(const BSTNode*, int&) const;

    int traverseBeachLeft(BSTNode*) const;
    int traverseBeachRight(BSTNode*) const;

    void checkBeachLine(double) const;

    /**
     *  @brief This recovers the depth of longest branch in the tree below input node
     */
    int getTreeDepth(const BSTNode*) const;

    /**
     *  @brief Tree balancing functions
     */
    void     rebalance(BSTNode*);
    BSTNode* rotateWithLeftChild(BSTNode*);
    BSTNode* rotateWithRightChild(BSTNode*);

    BSTNode*       m_root;      // the root of all evil, er, the top node
    BSTNodeList    m_nodeVec;   // Use this to keep track of the nodes

    EventUtilities m_utilities;
};

} // namespace lar_cluster3d
#endif

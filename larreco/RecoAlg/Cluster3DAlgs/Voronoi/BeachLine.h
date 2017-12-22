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

// std includes
#include <vector>
#include <algorithm>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief Define a virtual interface to the beachline "events"
 */
class IEvent
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IEvent() noexcept = default;
    
    /**
     *  @brief Interface for configuring the particular algorithm tool
     */
    
    virtual void setInvalid() = 0;
    
    virtual bool isSite()   const = 0;
    virtual bool isCircle() const = 0;
    virtual bool isValid()  const = 0;
    
    virtual float beachLinePos() const = 0;
    virtual float xPos()         const = 0;

    virtual bool operator<(const IEvent& right) const = 0;
    
    virtual bool newSiteToLeft(const IEvent*, const IEvent*, const IEvent*) const = 0;
};

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
        m_event(NULL), m_parent(NULL), m_leftChild(NULL), m_rightChild(NULL), m_predecessor(NULL), m_successor(NULL)
    {}
    
    BSTNode(IEvent* event) :
    m_event(event), m_parent(NULL), m_leftChild(NULL), m_rightChild(NULL), m_predecessor(NULL), m_successor(NULL)
    {}
    
    BSTNode(IEvent*, BSTNode*, BSTNode*, BSTNode*);

    /**
     *  @brief  Virtual Destructor
     */
    ~BSTNode() {}

    /**
     *  @brief recover the data members
     */
    IEvent*  getEvent()       const {return m_event;}
    BSTNode* getParent()      const {return m_parent;}
    BSTNode* getLeftChild()   const {return m_leftChild;}
    BSTNode* getRightChild()  const {return m_rightChild;}
    BSTNode* getPredecessor() const {return m_predecessor;}
    BSTNode* getSuccessor()   const {return m_successor;}
    
    /**
     *  @brief Allow setting of the points
     */
    void setParent(BSTNode* node)      {m_parent      = node;}
    void setLeftChild(BSTNode* node)   {m_leftChild   = node;}
    void setRightChild(BSTNode* node)  {m_rightChild  = node;}
    void setPredecessor(BSTNode* node) {m_predecessor = node;}
    void setSuccessor(BSTNode* node)   {m_successor   = node;}
    
    /**
     *  @brief Provide override definition for ordering
     */
    bool operator<(const BSTNode&) const;

private:
    IEvent*   m_event;       // Pointer to the event object
    BSTNode* m_parent;       // Tree traversal - parent node
    BSTNode* m_leftChild;    // Tree traversal - left child node
    BSTNode* m_rightChild;   // Tree traversal - right child node
    BSTNode* m_predecessor;  // Beachline traversal - predecessor
    BSTNode* m_successor;    // Beachline traversal - successor
};

/**
 * @brief This defines the actual beach line. The idea is to implement this as a
 *        self balancing binary search tree.
 *
 */
    
class BeachLine
{
public:
    BeachLine(int numSites);
    
    bool isEmpty() const {return m_root == NULL;}
    
    void setEmpty() {m_root = NULL;}
    
    BSTNode* findBestLeaf(IEvent* event) const {return findBestLeaf(event, m_root);}
    
    void insertNewLeaf(IEvent* event) {m_root = insertNewLeaf(event, m_root);}
    
    int getHeight(BSTNode*) const;
    
private:
    BSTNode* insertNewLeaf(IEvent*, BSTNode*);
    
    BSTNode* findBestLeaf(IEvent*, BSTNode*) const;

    BSTNode*             m_root;      // the root of all evil, er, the top node
    int                  m_numSites;  // Number of sites to process
    std::vector<BSTNode> m_nodeVec;   // Use this to keep track of the nodes
};
    
} // namespace lar_cluster3d
#endif

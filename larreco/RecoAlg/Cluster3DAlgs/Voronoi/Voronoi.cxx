/**
 *  @file   VoronoiDiagram.cxx
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// std includes
#include <iostream>
#include <numeric>
#include <functional>

// Framework Includes
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/Voronoi.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/DCEL.h"
#include "larreco/RecoAlg/Cluster3DAlgs/ConvexHull/ConvexHull.h"

// LArSoft includes

// boost includes
#include <boost/range/adaptor/reversed.hpp>
#include <boost/polygon/voronoi.hpp>

// Declare this here for boost
namespace boost
{
    namespace polygon
    {
        template <>
        struct geometry_concept<dcel2d::Point>
        {
            typedef point_concept type;
        };
        
        template <>
        struct point_traits<dcel2d::Point>
        {
            typedef int coordinate_type;
            
            static inline coordinate_type get(const dcel2d::Point& point, orientation_2d orient)
            {
                return (orient == HORIZONTAL) ? std::get<1>(point) : std::get<0>(point);
            }
        };
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace voronoi2d {

VoronoiDiagram::VoronoiDiagram(dcel2d::HalfEdgeList& halfEdgeList, dcel2d::VertexList& vertexList, dcel2d::FaceList& faceList) :
    fHalfEdgeList(halfEdgeList),
    fVertexList(vertexList),
    fFaceList(faceList),
    fXMin(0.),
    fXMax(0.),
    fYMin(0.),
    fYMax(0.),
    fVoronoiDiagramArea(0.)
{
    fHalfEdgeList.clear();
    fVertexList.clear();
    fFaceList.clear();
    fPointList.clear();
    fSiteEventList.clear();
    fCircleEventList.clear();
    fCircleNodeList.clear();
    fConvexHullList.clear();
    
    // And the area
    fVoronoiDiagramArea = Area();
}

//------------------------------------------------------------------------------------------------------------------------------------------

VoronoiDiagram::~VoronoiDiagram()
{
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
bool VoronoiDiagram::isLeft(const dcel2d::Point& p0, const dcel2d::Point& p1, const dcel2d::Point& pCheck) const
{
    // Use the cross product to determine if the check point lies to the left, on or right
    // of the line defined by points p0 and p1
    return crossProduct(p0, p1, pCheck) > 0;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

double VoronoiDiagram::crossProduct(const dcel2d::Point& p0, const dcel2d::Point& p1, const dcel2d::Point& p2) const
{
    // Define a quick 2D cross product here since it will used quite a bit!
    double deltaX  = std::get<0>(p1) - std::get<0>(p0);
    double deltaY  = std::get<1>(p1) - std::get<1>(p0);
    double dCheckX = std::get<0>(p2) - std::get<0>(p0);
    double dCheckY = std::get<1>(p2) - std::get<1>(p0);
    
    return ((deltaX * dCheckY) - (deltaY * dCheckX));
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
double VoronoiDiagram::Area() const
{
    double area(0.);
    
    // Compute the area by taking advantage of
    // 1) the ability to decompose a convex hull into triangles,
    // 2) the ability to use the cross product to calculate the area
    // So, the technique is to pick a point (for technical reasons we use 0,0)
    // and then sum the signed area of triangles formed from this point to two adjecent
    // vertices on the convex hull.
    dcel2d::Point center(0.,0.,NULL);
    dcel2d::Point lastPoint = fPointList.front();
    
    for(const auto& point : fPointList)
    {
        if (point != lastPoint) area += 0.5 * crossProduct(center,lastPoint,point);
        
        lastPoint = point;
    }
    
    return area;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool compareSiteEventPtrs(const IEvent* left, const IEvent* right) {return *left < *right;}

//------------------------------------------------------------------------------------------------------------------------------------------

void VoronoiDiagram::buildVoronoiDiagram(const dcel2d::PointList& pointList)
{
    // Insure all the local data structures have been cleared
    fHalfEdgeList.clear();
    fVertexList.clear();
    fFaceList.clear();
    fPointList.clear();
    fSiteEventList.clear();
    fCircleEventList.clear();
    fCircleNodeList.clear();
    fNumBadCircles = 0;
    
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "==> # input points: " << pointList.size() << std::endl;
    
    // Define the priority queue to contain our events
    EventQueue eventQueue(compareSiteEventPtrs);
    
    // Now populate the event queue with site events
    for(const auto& point : pointList)
    {
        fSiteEventList.emplace_back(point);
        IEvent* iEvent = &fSiteEventList.back();
        eventQueue.push(iEvent);
    }

    // Declare the beachline which will contain the BSTNode objects for site events
    BeachLine beachLine;
    
    // Finally, we need a container for our circle event BSTNodes
    BSTNodeList circleNodeList;
    
    // Now process the queue
    while(!eventQueue.empty())
    {
        IEvent* event = eventQueue.top();
        
        eventQueue.pop();
        
        // If a site or circle event then handle appropriately
        if      (event->isSite())  handleSiteEvents(beachLine, eventQueue, event);
        else if (event->isValid()) handleCircleEvents(beachLine, eventQueue, event);
    }
    
    std::cout << "*******> # input points: " << pointList.size() << ", remaining leaves: " << beachLine.countLeaves() << ", " << beachLine.traverseBeach() << ", # bad circles: " << fNumBadCircles << std::endl;
    std::cout << "           Faces: " << fFaceList.size() << ", Vertices: " << fVertexList.size() << ", # half edges: " << fHalfEdgeList.size() << std::endl;
    
    // Get the bounding box
    findBoundingBox(fVertexList);
    
    std::cout << "           Range min/maxes, x: " << fXMin << ", " << fXMax << ", y: " << fYMin << ", " << fYMax << std::endl;

    // Terminate the infinite edges
    terminateInfiniteEdges(beachLine, std::get<0>(pointList.front()));
    
    // Look for open faces
    int nOpenFaces(0);
    
    std::map<int,int> edgeCountMap;
    std::map<int,int> openCountMap;

    for (const auto& face : fFaceList)
    {
        int                     nEdges(1);
        bool                    closed(false);
        const dcel2d::HalfEdge* startEdge = face.getHalfEdge();
        
        // Count forwards
        for(const dcel2d::HalfEdge* halfEdge = startEdge->getNextHalfEdge(); halfEdge && !closed;)
        {
            if (halfEdge->getFace() != &face)
            {
                std::cout << "  ===> halfEdge does not match face: " << halfEdge << ", face: " << halfEdge->getFace() << ", base: " << &face << std::endl;
            }
            
            if (halfEdge == startEdge)
            {
                closed = true;
                break;
            }
            nEdges++;
            halfEdge = halfEdge->getNextHalfEdge();
        }
        
        // Count backwards. Note if face is closed then nothing to do here
        if (!closed)
        {
            for(const dcel2d::HalfEdge* halfEdge = startEdge->getLastHalfEdge(); halfEdge && !closed;)
            {
                if (halfEdge->getFace() != &face)
                {
                    std::cout << "  ===> halfEdge does not match face: " << halfEdge << ", face: " << halfEdge->getFace() << ", base: " << &face << std::endl;
                }
                
                if (halfEdge == startEdge)
                {
                    closed = true;
                    break;
                }
                nEdges++;
                halfEdge = halfEdge->getLastHalfEdge();
            }
        }
        
        if (!closed)
        {
            nOpenFaces++;
            openCountMap[nEdges]++;
        }
        
        edgeCountMap[nEdges]++;
    }
    
    std::cout << "==> Found " << nOpenFaces << " open faces from total of " << fFaceList.size() << std::endl;
    for(const auto& edgeCount : edgeCountMap) std::cout << "    -> all edges,  # edges: " << edgeCount.first << ", count: " << edgeCount.second << std::endl;
    for(const auto& edgeCount : openCountMap) std::cout << "    -> open edges, # edges: " << edgeCount.first << ", count: " << edgeCount.second << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;

    // Clear internal containers that are no longer useful
    fSiteEventList.clear();
    fCircleEventList.clear();
    fCircleNodeList.clear();
    
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VoronoiDiagram::buildVoronoiDiagramBoost(const dcel2d::PointList& pointList)
{
    // Insure all the local data structures have been cleared
    fHalfEdgeList.clear();
    fVertexList.clear();
    fFaceList.clear();
    fPointList.clear();
    fSiteEventList.clear();
    fCircleEventList.clear();
    fCircleNodeList.clear();
    fNumBadCircles = 0;
    
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "==> # input points: " << pointList.size() << std::endl;

    // Construct out voronoi diagram
    boost::polygon::voronoi_diagram<double> vd;
    boost::polygon::construct_voronoi(pointList.begin(),pointList.end(),&vd);
    
    // Make maps for translating from boost to me (or we can rewrite our code in boost... for now maps)
    using BoostEdgeToEdgeMap     = std::map<const boost::polygon::voronoi_edge<double>*,   dcel2d::HalfEdge*>;
    using BoostVertexToVertexMap = std::map<const boost::polygon::voronoi_vertex<double>*, dcel2d::Vertex*>;
    using BoostCellToFaceMap     = std::map<const boost::polygon::voronoi_cell<double>*,   dcel2d::Face*>;
    
    BoostEdgeToEdgeMap     boostEdgeToEdgeMap;
    BoostVertexToVertexMap boostVertexToVertexMap;
    BoostCellToFaceMap     boostCellToFaceMap;
    
    // Loop over the edges
    for(const auto& edge : vd.edges())
    {
        const boost::polygon::voronoi_edge<double>* twin = edge.twin();
        
        boostTranslation(pointList, &edge, twin, boostEdgeToEdgeMap, boostVertexToVertexMap, boostCellToFaceMap);
        boostTranslation(pointList, twin, &edge, boostEdgeToEdgeMap, boostVertexToVertexMap, boostCellToFaceMap);
    }
    
    //std::cout << "==> Found " << nOpenFaces << " open faces from total of " << fFaceList.size() << std::endl;
    //for(const auto& edgeCount : edgeCountMap) std::cout << "    -> all edges,  # edges: " << edgeCount.first << ", count: " << edgeCount.second << std::endl;
    //for(const auto& edgeCount : openCountMap) std::cout << "    -> open edges, # edges: " << edgeCount.first << ", count: " << edgeCount.second << std::endl;
    std::cout << "==> Boost returns " << fHalfEdgeList.size() << " edges, " << fFaceList.size() << " faces, " << fVertexList.size() << " vertices " << std::endl;
    std::cout << "                  " << vd.edges().size()    << " edges, " << vd.cells().size() << " faces, " << vd.vertices().size() << " vertices" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    std::cout << "******************************************************************************************************************" << std::endl;
    
    // Clear internal containers that are no longer useful
    fSiteEventList.clear();
    fCircleEventList.clear();
    fCircleNodeList.clear();
    
    return;
}
    
void VoronoiDiagram::boostTranslation(const dcel2d::PointList&                    pointList,
                                      const boost::polygon::voronoi_edge<double>* edge,
                                      const boost::polygon::voronoi_edge<double>* twin,
                                      BoostEdgeToEdgeMap&                         boostEdgeToEdgeMap,
                                      BoostVertexToVertexMap&                     boostVertexToVertexMap,
                                      BoostCellToFaceMap&                         boostCellToFaceMap)
{
    dcel2d::HalfEdge* halfEdge = NULL;
    dcel2d::HalfEdge* twinEdge = NULL;
   
    if (boostEdgeToEdgeMap.find(edge) != boostEdgeToEdgeMap.end()) halfEdge = boostEdgeToEdgeMap.at(edge);
    else
    {
        fHalfEdgeList.emplace_back();
    
        halfEdge = &fHalfEdgeList.back();
    
        boostEdgeToEdgeMap[edge] = halfEdge;
    }
    
    if (boostEdgeToEdgeMap.find(twin) != boostEdgeToEdgeMap.end()) twinEdge = boostEdgeToEdgeMap.at(twin);
    else
    {
        fHalfEdgeList.emplace_back();
        
        twinEdge = &fHalfEdgeList.back();
        
        boostEdgeToEdgeMap[twin] = twinEdge;
    }
    
    // Do the primary half edge first
    const boost::polygon::voronoi_vertex<double>* boostVertex = edge->vertex1();
    dcel2d::Vertex*                               vertex      = NULL;
    
    // note we can have a null vertex (infinite edge)
    if (boostVertex)
    {
        if (boostVertexToVertexMap.find(boostVertex) == boostVertexToVertexMap.end())
        {
            dcel2d::Coords coords(boostVertex->y(),boostVertex->x(),0.);
        
            fVertexList.emplace_back(coords, halfEdge);
        
            vertex = &fVertexList.back();
        
            boostVertexToVertexMap[boostVertex] = vertex;
        }
        else vertex = boostVertexToVertexMap.at(boostVertex);
    }
    
    const boost::polygon::voronoi_cell<double>* boostCell = edge->cell();
    dcel2d::Face*                               face      = NULL;
    
    if (boostCellToFaceMap.find(boostCell) == boostCellToFaceMap.end())
    {
        dcel2d::PointList::const_iterator pointItr = pointList.begin();
        int                               pointIdx = boostCell->source_index();
        
        std::advance(pointItr, pointIdx);
        
        const dcel2d::Point& point = *pointItr;
        dcel2d::Coords       coords(std::get<0>(point),std::get<1>(point),0.);
        
        fFaceList.emplace_back(halfEdge,coords,std::get<2>(point));
        
        face = &fFaceList.back();
        
        boostCellToFaceMap[boostCell] = face;
    }
    
    halfEdge->setTargetVertex(vertex);
    halfEdge->setFace(face);
    halfEdge->setTwinHalfEdge(twinEdge);
    
    // For the prev/next half edges we can have two cases, so check:
    if (boostEdgeToEdgeMap.find(edge->next()) != boostEdgeToEdgeMap.end())
    {
        dcel2d::HalfEdge* nextEdge = boostEdgeToEdgeMap.at(edge->next());
        
        halfEdge->setNextHalfEdge(nextEdge);
        nextEdge->setLastHalfEdge(halfEdge);
    }
    
    if (boostEdgeToEdgeMap.find(edge->prev()) != boostEdgeToEdgeMap.end())
    {
        dcel2d::HalfEdge* lastEdge = boostEdgeToEdgeMap.at(edge->prev());
        
        halfEdge->setLastHalfEdge(lastEdge);
        lastEdge->setNextHalfEdge(halfEdge);
    }
    
    return;
}

void VoronoiDiagram::handleSiteEvents(BeachLine&  beachLine,
                                      EventQueue& eventQueue,
                                      IEvent*     siteEvent)
{
    // Insert the new site event into the beach line and recover the leaf for the
    // new arc in the beach line
    // NOTE: invalidation of any possible circle events occurs in the call to this routine
    BSTNode* newLeaf = beachLine.insertNewLeaf(siteEvent);
    
    // Create a new vertex for each site event
    fFaceList.emplace_back(dcel2d::Face(NULL,siteEvent->getCoords(),std::get<2>(siteEvent->getPoint())));
    
    dcel2d::Face* face = &fFaceList.back();
    
    newLeaf->setFace(face);
    
    // If we are the first site added to the BST then we are done
    if (!(newLeaf->getPredecessor() || newLeaf->getSuccessor())) return;
    
    // So, now we deal with creating the edges that will be mapped out by the two breakpoints as they
    // move with the beachline. Note that this is a single edge pair (edge and its twin) between the
    // two breakpoints that have been created.
    fHalfEdgeList.push_back(dcel2d::HalfEdge());
    
    dcel2d::HalfEdge* halfEdge = &fHalfEdgeList.back();
    
    fHalfEdgeList.push_back(dcel2d::HalfEdge());
    
    dcel2d::HalfEdge* halfTwin = &fHalfEdgeList.back();
    
    // Each half edge is the twin of the other
    halfEdge->setTwinHalfEdge(halfTwin);
    halfTwin->setTwinHalfEdge(halfEdge);

    // Point the breakpoint to the new edge
    newLeaf->getPredecessor()->setHalfEdge(halfEdge);
    newLeaf->getSuccessor()->setHalfEdge(halfTwin);
    
    // The second half edge corresponds to our face for the left breakpoint...
    face->setHalfEdge(halfTwin);
    halfTwin->setFace(face);

    // Update for the left leaf
    BSTNode* leftLeaf = newLeaf->getPredecessor()->getPredecessor();
    
    // This needs to be checked since the first face will not have an edge set to it
    if (!leftLeaf->getFace()->getHalfEdge()) leftLeaf->getFace()->setHalfEdge(halfEdge);
    
    halfEdge->setFace(leftLeaf->getFace());
    
    // Try to make circle events either side of this leaf
    makeLeftCircleEvent( eventQueue, newLeaf, siteEvent->xPos());
    makeRightCircleEvent(eventQueue, newLeaf, siteEvent->xPos());

    return;
}
    
void VoronoiDiagram::handleCircleEvents(BeachLine&  beachLine,
                                        EventQueue& eventQueue,
                                        IEvent*     circleEvent)
{
    BSTNode* circleNode = circleEvent->getBSTNode();
    BSTNode* arcNode    = circleNode->getAssociated();

    // Recover the half edges for the breakpoints either side of the arc we're removing
    dcel2d::HalfEdge* leftHalfEdge  = arcNode->getPredecessor()->getHalfEdge();
    dcel2d::HalfEdge* rightHalfEdge = arcNode->getSuccessor()->getHalfEdge();
    
    if (leftHalfEdge->getTwinHalfEdge()->getFace() != rightHalfEdge->getFace())
    {
        std::cout << ">>>> Face mismatch in circle handling, left face: " << leftHalfEdge->getFace() << ", left twin: " << leftHalfEdge->getTwinHalfEdge()->getFace() << ", right: " << rightHalfEdge->getFace() << ", right twin: " << rightHalfEdge->getTwinHalfEdge()->getFace() << ", arc face: " << arcNode->getFace() << std::endl;
    }

    // Create the new vertex point
    // Don't forget we need to swap coordinates when returning to the outside world
    dcel2d::Coords vertexPos(circleEvent->circleCenter());
    
    fVertexList.push_back(dcel2d::Vertex(vertexPos,leftHalfEdge->getTwinHalfEdge()));
    
    dcel2d::Vertex* vertex = &fVertexList.back();
    
    // The edges we obtained "emanate" from the new vertex point,
    // their twins will point to it
    leftHalfEdge->getTwinHalfEdge()->setTargetVertex(vertex);
    rightHalfEdge->getTwinHalfEdge()->setTargetVertex(vertex);

    // Remove the arc which will return the new breakpoint
    // NOTE: invalidation of any possible circle events occurs in the call to this routine
    BSTNode* newBreakPoint = beachLine.removeLeaf(arcNode);
    
    // Now create edges associated with the new break point
    fHalfEdgeList.push_back(dcel2d::HalfEdge());
    
    dcel2d::HalfEdge* halfEdgeOne = &fHalfEdgeList.back();
    
    fHalfEdgeList.push_back(dcel2d::HalfEdge());
    
    dcel2d::HalfEdge* halfEdgeTwo = &fHalfEdgeList.back();
    
    // Associate the first edge with the new break point
    newBreakPoint->setHalfEdge(halfEdgeTwo);
    
    // These are twins
    halfEdgeOne->setTwinHalfEdge(halfEdgeTwo);
    halfEdgeTwo->setTwinHalfEdge(halfEdgeOne);
    
    // The second points to the vertex we just made
    halfEdgeTwo->setTargetVertex(vertex);
    vertex->setHalfEdge(halfEdgeOne);

    // halfEdgeTwo points to the vertex, face to left, so pairs with
    // the leftHalfEdge (leaving the vertex, face to left)
    halfEdgeTwo->setNextHalfEdge(leftHalfEdge);
    leftHalfEdge->setLastHalfEdge(halfEdgeTwo);
    halfEdgeTwo->setFace(leftHalfEdge->getFace());
    
    // halfEdgeOne emanates from the vertex, face to left, so pairs with
    // the twin of the rightHalfEdge (pointing to vertex, face to left)
    halfEdgeOne->setLastHalfEdge(rightHalfEdge->getTwinHalfEdge());
    rightHalfEdge->getTwinHalfEdge()->setNextHalfEdge(halfEdgeOne);
    halfEdgeOne->setFace(rightHalfEdge->getTwinHalfEdge()->getFace());

    // Finally, the twin of the leftHalfEdge points to the vertex (face to left)
    // and will pair with the rightHalfEdge emanating from the vertex
    // In this case they should already be sharing the same face
    rightHalfEdge->setLastHalfEdge(leftHalfEdge->getTwinHalfEdge());
    leftHalfEdge->getTwinHalfEdge()->setNextHalfEdge(rightHalfEdge);
    
    // We'll try to make circle events with the remnant arcs so get the pointers now
    // Note that we want the former left and right arcs to be the middle arcs in this case
    BSTNode* leftArc  = newBreakPoint->getPredecessor();
    BSTNode* rightArc = newBreakPoint->getSuccessor();
    
    // Look for a new circle candidates
    // In this case, we'd like the right arc to be the middle of the right circle,
    // the left arc to be the middle of the left circle, hence the logic below
    makeRightCircleEvent(eventQueue, leftArc, circleEvent->xPos());
    makeLeftCircleEvent(eventQueue, rightArc, circleEvent->xPos());

    return;
}
    
void VoronoiDiagram::makeLeftCircleEvent(EventQueue&  eventQueue,
                                         BSTNode*     leaf,
                                         double       beachLine)
{
    
    // Check status of triplet of site events to the left of this new leaf
    if (leaf->getPredecessor())
    {
        BSTNode* midLeaf = leaf->getPredecessor()->getPredecessor();
        
        if (midLeaf && midLeaf->getPredecessor())
        {
            BSTNode* edgeLeaf = midLeaf->getPredecessor()->getPredecessor();
            
            // edge leaves must be different
            if (leaf->getEvent()->getPoint() == edgeLeaf->getEvent()->getPoint()) return;
            
            if (edgeLeaf)
            {
                IEvent* circleEvent = makeCircleEvent(edgeLeaf, midLeaf, leaf, beachLine);
                
                // Did we succeed in making a circle event?
                if (circleEvent)
                {
                    // Add to the circle node list
                    fCircleNodeList.emplace_back(circleEvent);
                    
                    BSTNode* circleNode = &fCircleNodeList.back();
                    
                    // If there was an associated circle event to this node, invalidate it
                    if (midLeaf->getAssociated())
                    {
                        midLeaf->getAssociated()->setAssociated(NULL);
                        midLeaf->getAssociated()->getEvent()->setInvalid();
                    }
                    
                    // Now reset to point at the new one
                    midLeaf->setAssociated(circleNode);
                    circleNode->setAssociated(midLeaf);

                    eventQueue.push(circleEvent);
                }
            }
        }
    }
    
    return;
}
    
void VoronoiDiagram::makeRightCircleEvent(EventQueue& eventQueue,
                                          BSTNode*    leaf,
                                          double      beachLine)
{
    // Check status of triplet of site events to the left of this new leaf
    if (leaf->getSuccessor())
    {
        BSTNode* midLeaf = leaf->getSuccessor()->getSuccessor();
        
        if (midLeaf && midLeaf->getSuccessor())
        {
            BSTNode* edgeLeaf = midLeaf->getSuccessor()->getSuccessor();
            
            if (leaf->getEvent()->getPoint() == edgeLeaf->getEvent()->getPoint()) return;

            if (edgeLeaf)
            {
                IEvent* circleEvent = makeCircleEvent(leaf, midLeaf, edgeLeaf, beachLine);
                
                // Did we succeed in making a circle event?
                if (circleEvent)
                {
                    // Add to the circle node list
                    fCircleNodeList.emplace_back(circleEvent);
                    
                    BSTNode* circleNode = &fCircleNodeList.back();
                    
                    // If there was an associated circle event to this node, invalidate it
                    if (midLeaf->getAssociated())
                    {
                        midLeaf->getAssociated()->setAssociated(NULL);
                        midLeaf->getAssociated()->getEvent()->setInvalid();
                    }
                    
                    // Now reset to point at the new one
                    midLeaf->setAssociated(circleNode);
                    circleNode->setAssociated(midLeaf);

                    eventQueue.push(circleEvent);
                }
            }
        }
    }
    
    return;
}

IEvent* VoronoiDiagram::makeCircleEvent(BSTNode* arc1, BSTNode* arc2, BSTNode* arc3, double beachLinePos)
{
    // It might be that we don't create a new circle
    IEvent* circle = 0;
    
    // First step is to calculate the center and radius of the circle determined by the three input site events
    dcel2d::Coords center;
    double         radius;
    double         deltaR;
    
    IEvent* p1 = arc1->getEvent();
    IEvent* p2 = arc2->getEvent();
    IEvent* p3 = arc3->getEvent();

    // Compute the center of the circle. Note that this will also automagically check that breakpoints are
    // converging so that this is a circle we want
    if (computeCircleCenter( p1->getCoords(), p2->getCoords(), p3->getCoords(), center,  radius, deltaR))
    {
        double circleBottomX = center[0] - radius;
        
        // Now check if the bottom of this circle lies below the beach line
        if (beachLinePos >= circleBottomX - 10. * deltaR)
        {
            // Making a circle event!
            dcel2d::Point circleBottom(circleBottomX, center[1], NULL);
        
            fCircleEventList.emplace_back(circleBottom, center);
        
            circle = &fCircleEventList.back();
        }
        else if (circleBottomX - beachLinePos < 1.e-4)
            std::cout << "==> Circle close, beachLine: " << beachLinePos << ", circleBottomX: " << circleBottomX << ", deltaR: " << deltaR << ", d: " << circleBottomX - beachLinePos << std::endl;
    }
    
    return circle;
}

bool VoronoiDiagram::computeCircleCenter(const dcel2d::Coords& p1,
                                         const dcel2d::Coords& p2,
                                         const dcel2d::Coords& p3,
                                         dcel2d::Coords&       center,
                                         double&               radius,
                                         double&               delta) const
{
    // The method is to translate the three points to a system where the first point is at the origin. Then we
    // are looking for a circle that passes through the origin and the two remaining (translated) points. In
    // this mode the circle radius is the distance from the origin to the circle center which simplifies the
    // calculation.
    double xCoord = p1[0];
    double yCoord = p1[1];
    
    double x2     = p2[0] - xCoord;
    double y2     = p2[1] - yCoord;
    double x3     = p3[0] - xCoord;
    double y3     = p3[1] - yCoord;
    
    double det    = x2 * y3 - x3 * y2;
    
    // Points are colinear so cannot make a circle
    // Or, if negative then the midpoint is "left" of line between first and third points meaning the
    // circle curvature is the "wrong" way for making a circle event
    if (det <= 0.) // std::numeric_limits<double>::epsilon())
    //if (!(std::abs(det) > 0.)) // std::numeric_limits<double>::epsilon())
    {
        if (det > -std::numeric_limits<double>::epsilon())
            std::cout << "      --->Circle failure, det: " << det << ", mid x: " << p2[0] << ", y: " << p2[1]
                                                                  << ",   l x: " << p1[0] << ", y: " << p1[1]
                                                                  << ",   r x: " << p3[0] << ", y: " << p3[1] << std::endl;
        
        return false;
    }
    
    double p2sqr   = x2 * x2 + y2 * y2;
    double p3sqr   = x3 * x3 + y3 * y3;
    
    double cxpr    = 0.5 * (y3 * p2sqr - y2 * p3sqr) / det;
    double cypr    = 0.5 * (x2 * p3sqr - x3 * p2sqr) / det;
    
    radius    = std::sqrt(cxpr * cxpr + cypr * cypr);
    center[0] = cxpr + xCoord;
    center[1] = cypr + yCoord;
    center[2] = 0.;

    // So... roundoff error can cause degenerate circles to get missed
    // We need to attack that problem by calculating the radius all possible
    // ways and then take the largest...
    dcel2d::Coords p1Rad = p1 - center;
    dcel2d::Coords p2Rad = p2 - center;
    dcel2d::Coords p3Rad = p3 - center;

    std::vector<float> radSqrVec;

    radSqrVec.push_back(p1Rad.norm());
    radSqrVec.push_back(p2Rad.norm());
    radSqrVec.push_back(p3Rad.norm());

    double maxRadius = *std::max_element(radSqrVec.begin(),radSqrVec.end());
    
    delta  = std::max(5.e-7, maxRadius - radius);
    
    if (radius > 1000.)
    {
//        std::cout << "       ***> Radius = " << radius << ", circ x,y: " << center[0] << ", " << center[1] << ", p1: " << xCoord << "," << yCoord << ", p2: " << p2[0] << "," << p2[1] << ", p3: " << p3[0] << "," << p3[1] << ", det: " << det << std::endl;
        fNumBadCircles++;
    }
    
    return true;
}

bool VoronoiDiagram::computeCircleCenter2(const dcel2d::Coords& p1,
                                          const dcel2d::Coords& p2,
                                          const dcel2d::Coords& p3,
                                          dcel2d::Coords&       center,
                                          double&               radius,
                                          double&               delta) const
{
    // Compute the circle center as the intersection of the two perpendicular bisectors of rays between the points
    double slope12 = (p2[1] - p1[1]) / (p2[0] - p1[0]);
    double slope32 = (p3[1] - p2[1]) / (p3[0] - p2[0]);
    
    if (std::abs(slope12 - slope32) <= std::numeric_limits<double>::epsilon())
    {
        std::cout << "       >>>> Matching slopes! points: (" << p1[0] << "," << p1[1] << "), ("<< p2[0] << "," << p2[1] << "), ("<< p3[0] << "," << p3[1] << ")" << std::endl;
        
        return false;
    }
    
    center[0] = (slope12 * slope32 * (p3[1] - p1[1]) + slope12 * (p2[0] + p3[0]) - slope32 * (p1[0] + p2[0])) / (2. * (slope12 - slope32));
    center[1] = 0.5 * (p1[1] + p2[1]) - (center[0] - 0.5 * (p1[0] + p2[0])) / slope12;
    center[2] = 0.;
    radius    = std::sqrt(std::pow((p2[0] - center[0]),2) + std::pow((p2[1] - center[1]),2));

    if (radius > 100.)
    {
        std::cout << "       ***> Rad2 = " << radius << ", circ x,y: " << center[0] << "," << center[1] << ", p1: " << p1[0] << "," << p1[1] << ", p2: " << p2[0] << "," << p2[1] << ", p3: " << p3[0] << ", " << p3[1] << std::endl;
    }
    
    return true;
}

bool VoronoiDiagram::computeCircleCenter3(const dcel2d::Coords& p1,
                                          const dcel2d::Coords& p2,
                                          const dcel2d::Coords& p3,
                                          dcel2d::Coords&       center,
                                          double&               radius,
                                          double&               delta) const
{
    // Yet another bisector method to calculate the circle center...
    double temp = p2[0] * p2[0] + p2[1] * p2[1];
    double p1p2 = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2.0;
    double p2p3 = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2.0;
    double det  = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1]);
    
    if (std::abs(det) < 10. * std::numeric_limits<double>::epsilon())
    {
        std::cout << "       >>>> Determinant zero! points: (" << p1[0] << "," << p1[1] << "), ("<< p2[0] << "," << p2[1] << "), ("<< p3[0] << "," << p3[1] << ")" << std::endl;
        
        return false;
    }

    det = 1. / det;
    
    center[0] = (p1p2 * (p2[1] - p3[1]) - p2p3 * (p1[1] - p2[1])) * det;
    center[1] = (p2p3 * (p1[0] - p2[0])   - p1p2 * (p2[0] - p3[0]))   * det;
    center[2] = 0.;
    
    radius        = std::sqrt(std::pow((p1[0] - center[0]),2) + std::pow((p1[1] - center[1]),2));
    
    if (radius > 100.)
    {
        std::cout << "       ***> Rad3 = " << radius << ", circ x,y: " << center[0] << "," << center[1] << ", p1: " << p1[0] << "," << p1[1] << ", p2: " << p2[0] << "," << p2[1] << ", p3: " << p3[0] << ", " << p3[1] << std::endl;
    }
    
    return true;
}
    
void VoronoiDiagram::terminateInfiniteEdges(BeachLine& beachLine, double beachLinePos)
{
    // Need to complete processing of the beachline, the remaning leaves represent the site points with "infinite"
    // edges which we need to terminate at our bounding box.
    // For now let's just do step one which is to find the convex hull
    const BSTNode* node = beachLine.getTopNode();
    
    // Do the convex hull independently but before terminating edges
    getConvexHull(node);
    
    if (node)
    {
        // Run down the left branches to find the left most leaf
        while(node->getLeftChild()) node = node->getLeftChild();
        
        // Things to keep track of...
        int    nodeCount(0);
        int    nBadBreaks(0);
        double lastBreakPoint = std::numeric_limits<double>::lowest();
        
        // Now we are going to traverse across the beach line and process the remaining leaves
        // This will identify the infinite edges and find the convex hull
        while(node)
        {
            std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
            
            // Do we have a leaf?
            if (!node->getLeftChild() && !node->getRightChild())
            {
                std::cout << "node " << nodeCount << " has leaf: " << node << ", face: " << node->getFace() << ", pos: " << node->getEvent()->getCoords()[0] << "/" << node->getEvent()->getCoords()[1] << std::endl;
            }
            // Otherwise it is a break point
            else
            {
                RootsPair roots;
                double breakPoint  = fUtilities.computeBreak(beachLinePos, node->getPredecessor()->getEvent(), node->getSuccessor()->getEvent(), roots);
                double leftArcVal  = fUtilities.computeArcVal(beachLinePos, breakPoint, node->getPredecessor()->getEvent());
                double rightArcVal = fUtilities.computeArcVal(beachLinePos, breakPoint, node->getSuccessor()->getEvent());
                
                dcel2d::HalfEdge* halfEdge = node->getHalfEdge();
                dcel2d::HalfEdge* halfTwin = halfEdge->getTwinHalfEdge();
                dcel2d::Vertex*   vertex   = halfEdge->getTargetVertex();

                std::cout << "node " << nodeCount << " has break: " << breakPoint << ", beachPos: " << beachLinePos << " (last: " << lastBreakPoint << "), leftArcVal: " << leftArcVal << ", rightArcVal: " << rightArcVal << std::endl;
                if (lastBreakPoint > breakPoint) std::cout << "     ***** lastBreakPoint larger than current break point *****" << std::endl;
                std::cout << "     left arc x/y: " << node->getPredecessor()->getEvent()->xPos() << "/" << node->getPredecessor()->getEvent()->yPos() << ", right arc x/y: "  << node->getSuccessor()->getEvent()->xPos() << "/" << node->getSuccessor()->getEvent()->yPos() << std::endl;
                std::cout << "     halfEdge: " << halfEdge << ", target vtx: " << vertex << ", face: " << halfEdge->getFace() << ", twin: " << halfTwin << ", twin tgt: " << halfTwin->getTargetVertex() << ", twin face: " << halfTwin->getFace() << ", next: " << halfEdge->getNextHalfEdge() << ", last: " << halfEdge->getLastHalfEdge() << std::endl;
                
                const dcel2d::Coords& leftLeafPos  = node->getPredecessor()->getEvent()->getCoords();
                const dcel2d::Coords& rightLeafPos = node->getSuccessor()->getEvent()->getCoords();
                Eigen::Vector3f lrPosDiff          = rightLeafPos - leftLeafPos;
                Eigen::Vector3f lrPosSum           = rightLeafPos + leftLeafPos;
                
                lrPosDiff.normalize();
                lrPosSum *= 0.5;
                
                dcel2d::Coords leafPos = node->getSuccessor()->getEvent()->getCoords();
                dcel2d::Coords leafDir = lrPosDiff;
                
                if (vertex)
                {
                    dcel2d::Coords breakDir(leafDir[1],-leafDir[0],0.);
                    dcel2d::Coords breakPos  = lrPosSum;
                    dcel2d::Coords vertexPos = vertex->getCoords();
                
                    // Get the arclength from the break position to the intersection of the two lines
                    dcel2d::Coords breakToLeafPos = leafPos - vertexPos;
                    double         arcLenToLine   = breakToLeafPos.dot(breakDir);
                
                    // Now get position
                    dcel2d::Coords breakVertexPos = vertexPos + arcLenToLine * breakDir;

                    std::cout << "     halfEdge position: " << breakPos[0] << "/" << breakPos[1] << ", vertex: " << vertexPos[0] << "/" << vertexPos[1] << ", end: " << breakVertexPos[0] << "/" << breakVertexPos[1] << ", dir: " << breakDir[0] << "/" << breakDir[1] << ", arclen: " << arcLenToLine  << std::endl;

                    // At this point we need to create a vertex and then terminate the half edges here...
                    fVertexList.push_back(dcel2d::Vertex(breakVertexPos,halfEdge));
                
                    dcel2d::Vertex* breakVertex = &fVertexList.back();
                
                    halfTwin->setTargetVertex(breakVertex);
                }
                else
                {
                    std::cout << "****** null vertex!!! Skipping to next node... *********" << std::endl;
                }

                if (lastBreakPoint > breakPoint) nBadBreaks++;

                lastBreakPoint = breakPoint;
            }
            
            if (node->getAssociated())
            {
                BSTNode* associated = node->getAssociated();
                IEvent*  iEvent     = associated->getEvent();
                
                std::cout << "     -> associated circle: " << iEvent->isCircle() << ", is valid: " << iEvent->isValid() << std::endl;
            }
            
            node = node->getSuccessor();
            nodeCount++;
        }

        // Now that we have the convex hull, loop over vertices to see if they are inside or outside of the convex hull
        dcel2d::VertexList::iterator curVertexItr = fVertexList.begin();
        
        size_t nVerticesInitial = fVertexList.size();
        
        // Loop over all vertices to begin with
        while(curVertexItr != fVertexList.end())
        {
            // Dereference vertex
            dcel2d::Vertex& vertex = *curVertexItr;
            
            bool outsideHull = !isInsideConvexHull(vertex);

            // Do we need to drop this vertex?
            if (outsideHull)
            {
                curVertexItr = fVertexList.erase(curVertexItr);
            }
            else curVertexItr++;
        }
        
        mergeDegenerateVertices();
        ComputeFaceArea();

        std::cout << "Loop over beachline done, saved " << fConvexHullList.size() << " arcs, encountered " << nBadBreaks << " bad break points" << std::endl;
        std::cout << "-- started with " << nVerticesInitial << " vertices, found " << fVertexList.size() << " inside convex hull" << std::endl;
    }
    
    return;
}

void VoronoiDiagram::getConvexHull(const BSTNode* topNode)
{
    // Assume the input node is the top of the binary search tree and represents
    // the beach line at the end of the sweep algorithm
    // The convex hull will then be the leaf elements along the beach line.
    // Note that this algorithm will give you the convex hull in a clockwise manner
    if (topNode)
    {
        const BSTNode* node = topNode;
        
        // Initialize the center
        fConvexHullCenter = dcel2d::Coords(0.,0.,0.);
        
        // Run down the left branches to find the right most leaf
        // We do it this way to make sure we get a CCW convex hull
        while(node->getRightChild()) node = node->getRightChild();
        
        // Include a check on convexity...
        Eigen::Vector3f prevVec(0.,0.,0.);
        dcel2d::Coords  lastPoint(0.,0.,0.);

        // Now we are going to traverse across the beach line and process the remaining leaves
        // This will identify the convex hull
        while(node)
        {
            if (!node->getLeftChild() && !node->getRightChild())
            {
                // Add point to the convex hull list
                fConvexHullList.emplace_back(node->getEvent()->getPoint());
                
                // Mark the face as on the convex hull
                node->getFace()->setOnConvexHull();
                
                dcel2d::Coords  nextPoint = node->getEvent()->getCoords();
                Eigen::Vector3f curVec    = nextPoint - lastPoint;
                
                curVec.normalize();
                
                double dotProd = prevVec.dot(curVec);
                
                std::cout << "--> lastPoint: " << lastPoint[0] << "/" << lastPoint[1] << ", tan: " << std::atan2(lastPoint[1],lastPoint[0]) << ", curPoint: " << nextPoint[0] << "/" << nextPoint[1] << ", tan: " << std::atan2(nextPoint[1],nextPoint[0]) << ", dot: " << dotProd << std::endl;
                
                prevVec   = curVec;
                lastPoint = nextPoint;
            }
            
            node = node->getPredecessor();
        }

        // Annoyingly, our algorithm does not contain only the convex hull points and so we need to skim out the renegades...
        lar_cluster3d::ConvexHull::PointList localList;
        
        for(const auto& edgePoint : fConvexHullList) localList.emplace_back(std::get<0>(edgePoint),std::get<1>(edgePoint),std::get<2>(edgePoint));
        
        // Sort the point vec by increasing x, then increase y
        localList.sort([](const auto& left, const auto& right){return (std::abs(std::get<0>(left) - std::get<0>(right)) > std::numeric_limits<float>::epsilon()) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);});
        
        // Why do I need to do this?
        lar_cluster3d::ConvexHull convexHull(localList);
        
        // Clear the convex hull list...
//        fConvexHullList.clear();
        
        std::cout << "~~~>> there are " << convexHull.getConvexHull().size() << " convex hull points and " << fConvexHullList.size() << " infinite cells" << std::endl;
        
        // Now rebuild it
        for(const auto& hullPoint : convexHull.getConvexHull())
        {
//            fConvexHullList.emplace_back(std::get<0>(hullPoint),std::get<1>(hullPoint),std::get<2>(hullPoint));

            std::cout << "~~~ Convex hull Point: " << std::get<0>(hullPoint) << ", " << std::get<1>(hullPoint) << std::endl;
        }
    }
    
    return;
}
    
bool VoronoiDiagram::isInsideConvexHull(const dcel2d::Vertex& vertex) const
{
    bool          insideHull(true);
    dcel2d::Point vertexPos(vertex.getCoords()[0],vertex.getCoords()[1],NULL);
    
    // Now check to see if the vertex is inside the convex hull
    dcel2d::PointList::const_iterator hullItr    = fConvexHullList.begin();
    dcel2d::Point                     firstPoint = *hullItr++;
    
    // We assume here the convex hull is stored in a CCW order
    while(hullItr != fConvexHullList.end())
    {
        dcel2d::Point secondPoint = *hullItr++;

        // CCW order means we check to see if this point lies to left of line from first to second point
        // in order for the point to be "inside" the convex hull
        // Note that we are looking to reject points that are outside...
        if (!isLeft(firstPoint,secondPoint,vertexPos))
        {
            insideHull = false;
            break;
        }
        
        firstPoint = secondPoint;
    }

    return insideHull;
}
    
bool VoronoiDiagram::isOutsideConvexHull(const dcel2d::Vertex&             vertex,
                                         dcel2d::PointList::const_iterator firstHullPointItr,
                                         dcel2d::Coords&                   intersection,
                                         double&                           distToConvexHull) const
{
    bool          outsideHull(false);
    dcel2d::Point vertexPos(vertex.getCoords()[0],vertex.getCoords()[1],NULL);
    
    distToConvexHull  = std::numeric_limits<double>::max();
    firstHullPointItr = fConvexHullList.begin();
    
    // Now check to see if the vertex is inside the convex hull
    dcel2d::PointList::const_iterator hullItr       = firstHullPointItr;
    dcel2d::PointList::const_iterator firstPointItr = hullItr++;
    
    while(hullItr != fConvexHullList.end())
    {
        dcel2d::PointList::const_iterator secondPointItr = hullItr++;
        
        // Dereference some stuff
        double xPrevToPoint = (std::get<0>(vertexPos)       - std::get<0>(*firstPointItr));
        double yPrevToPoint = (std::get<1>(vertexPos)       - std::get<1>(*firstPointItr));
        double xPrevToCur   = (std::get<0>(*secondPointItr) - std::get<0>(*firstPointItr));
        double yPrevToCur   = (std::get<1>(*secondPointItr) - std::get<1>(*firstPointItr));
        double edgeLength   = std::sqrt(xPrevToCur * xPrevToCur + yPrevToCur * yPrevToCur);
        
        // Find projection onto convex hull edge
        double projection = ((xPrevToPoint * xPrevToCur) + (yPrevToPoint * yPrevToCur)) / edgeLength;
        
        // DOCA point
        dcel2d::Point docaPoint(std::get<0>(*firstPointItr) + projection * xPrevToCur / edgeLength,
                                std::get<1>(*firstPointItr) + projection * yPrevToCur / edgeLength, 0);
        
        if (projection > edgeLength) docaPoint = *secondPointItr;
        if (projection < 0)          docaPoint = *firstPointItr;
        
        double xDocaDist = std::get<0>(vertexPos) - std::get<0>(docaPoint);
        double yDocaDist = std::get<1>(vertexPos) - std::get<1>(docaPoint);
        double docaDist  = xDocaDist * xDocaDist + yDocaDist * yDocaDist;
        
        if (docaDist < distToConvexHull)
        {
            firstHullPointItr = firstPointItr;
            intersection      = dcel2d::Coords(std::get<0>(docaPoint),std::get<1>(docaPoint),0.);
            distToConvexHull  = docaDist;
        }

        // Check to see if this point is outside the convex hull
        if (isLeft(*firstPointItr,*secondPointItr,vertexPos)) outsideHull = true;
        
        firstPointItr = secondPointItr;
    }
    
    return outsideHull;
}
    
void VoronoiDiagram::mergeDegenerateVertices()
{
    dcel2d::HalfEdgeList::iterator edgeItr = fHalfEdgeList.begin();
    
    while(edgeItr != fHalfEdgeList.end())
    {
        dcel2d::HalfEdge* halfEdge = &(*edgeItr);
        dcel2d::HalfEdge* twinEdge = halfEdge->getTwinHalfEdge();
        
        // Make sure we are not looking at an infinite edge
        if (halfEdge->getTargetVertex() && twinEdge->getTargetVertex())
        {
            dcel2d::Coords vtxPosDiff = halfEdge->getTargetVertex()->getCoords() - twinEdge->getTargetVertex()->getCoords();
            
            if (vtxPosDiff.norm() < 1.e-3)
            {
                std::cout << "***>> found a degenerate vertex! " << halfEdge->getTargetVertex()->getCoords()[0] << "/" << halfEdge->getTargetVertex()->getCoords()[1] << ", d: " << vtxPosDiff.norm() << std::endl;
                edgeItr++;
            }
            else edgeItr++;
        }
        else edgeItr++;
    }
    
    return;
}
    
    
double VoronoiDiagram::ComputeFaceArea()
{
    // Compute the area by taking advantage of
    // 1) the ability to decompose a convex hull into triangles,
    // 2) the ability to use the cross product to calculate the area
    // So the idea is to loop through all the faces and then follow the edges
    // around the face to compute the area of the face.
    // Note that a special case are the "infinite faces" which lie on the
    // convex hull of the event. Skip these for now...
    
    double totalArea(0.);
    int    nNonInfiniteFaces(0);
    double smallestArea(std::numeric_limits<double>::max());
    double largestArea(0.);
    
    std::vector<std::pair<double,const dcel2d::Face*>> areaFaceVec;
    
    areaFaceVec.reserve(fFaceList.size());
    
    for(auto& face : fFaceList)
    {
//        const dcel2d::Coords&   faceCoords = face.getCoords();
        const dcel2d::HalfEdge* halfEdge   = face.getHalfEdge();
        double            faceArea(0.);
        int               numEdges(0);
        bool              doNext     = true;
        
        dcel2d::Coords faceCenter(0.,0.,0.);
        
        while(doNext)
        {
            if (halfEdge->getTargetVertex())
                faceCenter += halfEdge->getTargetVertex()->getCoords();

            numEdges++;
            
            halfEdge = halfEdge->getNextHalfEdge();
            
            if (!halfEdge)
            {
                faceArea = std::numeric_limits<double>::max();
                doNext   = false;
            }
            
            if (halfEdge == face.getHalfEdge()) doNext = false;
        }
        
        faceCenter /= numEdges;
        
        halfEdge = face.getHalfEdge();
        doNext   = true;

        while(doNext)
        {
            const dcel2d::HalfEdge* twinEdge = halfEdge->getTwinHalfEdge();
            
            if (!halfEdge->getTargetVertex() || !twinEdge->getTargetVertex())
            {
                faceArea = std::numeric_limits<double>::max();
                break;
            }
            
            // Recover the two vertex points
            const dcel2d::Coords& p1 = halfEdge->getTargetVertex()->getCoords();
            const dcel2d::Coords& p2 = twinEdge->getTargetVertex()->getCoords();

            // Define a quick 2D cross product here since it will used quite a bit!
            double dp1p0X = p1[0] - faceCenter[0];
            double dp1p0Y = p1[1] - faceCenter[1];
            double dp2p0X = p2[0] - faceCenter[0];
            double dp2p0Y = p2[1] - faceCenter[1];
            
            //faceArea += dp1p0X * dp2p0Y - dp1p0Y * dp2p0X;
            double crossProduct = dp1p0X * dp2p0Y - dp1p0Y * dp2p0X;
            
            faceArea += crossProduct;
            
            if (crossProduct < 0.)
            {
                dcel2d::Coords edgeVec = p1 - p2;
                std::cout << "--- negative cross: " << crossProduct << ", edgeLen: " << edgeVec.norm() << ", x/y: " << edgeVec[0] << "/" << edgeVec[1] << std::endl;
            }
            
//            numEdges++;

            halfEdge = halfEdge->getNextHalfEdge();
            
            if (!halfEdge)
            {
                faceArea = std::numeric_limits<double>::max();
                break;
            }
            
            if (halfEdge == face.getHalfEdge()) doNext = false;
        }
        
        areaFaceVec.emplace_back(faceArea,&face);
        
        if (faceArea < std::numeric_limits<double>::max() && faceArea > 0.)
        {
            nNonInfiniteFaces++;
            totalArea    += faceArea;
            smallestArea  = std::min(faceArea,smallestArea);
            largestArea   = std::max(faceArea,largestArea);
        }
        
        if (faceArea < 1.e-4) std::cout << "---> face area <= 0: " << faceArea << ", with " << numEdges << " edges" << std::endl;
        
        face.setFaceArea(faceArea);
    }
    
    // Calculate a truncated mean...
    std::sort(areaFaceVec.begin(),areaFaceVec.end(),[](const auto& left, const auto& right){return left.first < right.first;});
    
    std::vector<std::pair<double,const dcel2d::Face*>>::iterator firstItr = std::find_if(areaFaceVec.begin(),areaFaceVec.end(),[](const auto& val){return val.first > 0.;});
    std::vector<std::pair<double,const dcel2d::Face*>>::iterator lastItr = std::find_if(areaFaceVec.begin(),areaFaceVec.end(),[](const auto& val){return !(val.first < std::numeric_limits<double>::max());});

    size_t nToKeep = 0.8 * std::distance(firstItr,lastItr);
    
    std::cout << ">>>>> nToKeep: " << nToKeep << ", last non infinite entry: " << std::distance(areaFaceVec.begin(),lastItr) << std::endl;
    
    double totalTruncArea = std::accumulate(firstItr,firstItr+nToKeep,0.,[](auto sum, const auto& val){return sum+val.first;});
    double aveTruncArea   = totalTruncArea / double(nToKeep);
    
    if (nNonInfiniteFaces > 0) std::cout << ">>>> Face area for " << nNonInfiniteFaces << ", ave area: " << totalArea / nNonInfiniteFaces << ", ave trunc area: " << aveTruncArea << ", ratio: " << totalTruncArea / totalArea << ", smallest: " << smallestArea << ", largest: " << largestArea << std::endl;
    else std::cout << ">>>>> there are no non infinite faces" << std::endl;
    
    return totalArea;
}


void VoronoiDiagram::findBoundingBox(const dcel2d::VertexList& vertexList)
{
    // Find extremes in x to start
    std::pair<dcel2d::VertexList::const_iterator,dcel2d::VertexList::const_iterator> minMaxItrX = std::minmax_element(vertexList.begin(),vertexList.end(),[](const auto& left, const auto& right){return left.getCoords()[0] < right.getCoords()[0];});
    
    fXMin = minMaxItrX.first->getCoords()[0];
    fXMax = minMaxItrX.second->getCoords()[0];

    // To get the extremes in y we need to make a pass through the list
    std::pair<dcel2d::VertexList::const_iterator,dcel2d::VertexList::const_iterator> minMaxItrY = std::minmax_element(vertexList.begin(),vertexList.end(),[](const auto& left, const auto& right){return left.getCoords()[1] < right.getCoords()[1];});
    
    fYMin = minMaxItrY.first->getCoords()[1];
    fYMax = minMaxItrY.second->getCoords()[1];

    return;
}

VoronoiDiagram::PointPair VoronoiDiagram::findNearestEdge(const dcel2d::Point& point, double& closestDistance) const
{
    // The idea is to find the nearest edge of the convex hull, defined by
    // two adjacent vertices of the hull, to the input point.
    // As near as I can tell, the best way to do this is to apply brute force...
    // Idea will be to iterate over pairs of points
    dcel2d::PointList::const_iterator curPointItr = fPointList.begin();
    dcel2d::Point                     prevPoint   = *curPointItr++;
    dcel2d::Point                     curPoint    = *curPointItr;
    
    // Set up the winner
    PointPair closestEdge(prevPoint,curPoint);
    
    closestDistance = std::numeric_limits<double>::max();
    
    // curPointItr is meant to point to the second point
    while(curPointItr != fPointList.end())
    {
        if (curPoint != prevPoint)
        {
            // Dereference some stuff
            double xPrevToPoint = (std::get<0>(point)    - std::get<0>(prevPoint));
            double yPrevToPoint = (std::get<1>(point)    - std::get<1>(prevPoint));
            double xPrevToCur   = (std::get<0>(curPoint) - std::get<0>(prevPoint));
            double yPrevToCur   = (std::get<1>(curPoint) - std::get<1>(prevPoint));
            double edgeLength   = std::sqrt(xPrevToCur * xPrevToCur + yPrevToCur * yPrevToCur);

            // Find projection onto convex hull edge
            double projection = ((xPrevToPoint * xPrevToCur) + (yPrevToPoint * yPrevToCur)) / edgeLength;
            
            // DOCA point
            dcel2d::Point docaPoint(std::get<0>(prevPoint) + projection * xPrevToCur / edgeLength,
                                    std::get<1>(prevPoint) + projection * yPrevToCur / edgeLength, 0);
            
            if (projection > edgeLength) docaPoint = curPoint;
            if (projection < 0)          docaPoint = prevPoint;
            
            double xDocaDist = std::get<0>(point) - std::get<0>(docaPoint);
            double yDocaDist = std::get<1>(point) - std::get<1>(docaPoint);
            double docaDist  = xDocaDist * xDocaDist + yDocaDist * yDocaDist;
            
            if (docaDist < closestDistance)
            {
                closestEdge     = PointPair(prevPoint,curPoint);
                closestDistance = docaDist;
            }
        }
        
        prevPoint = curPoint;
        curPoint  = *curPointItr++;
    }
    
    closestDistance = std::sqrt(closestDistance);

    // Convention is convex hull vertices sorted in counter clockwise fashion so if the point
    // lays to the left of the nearest edge then it must be an interior point
    if (isLeft(closestEdge.first,closestEdge.second,point)) closestDistance = -closestDistance;
    
    return closestEdge;
}

double VoronoiDiagram::findNearestDistance(const dcel2d::Point& point) const
{
    double closestDistance;
    
    findNearestEdge(point,closestDistance);
    
    return closestDistance;
}

} // namespace lar_cluster3d

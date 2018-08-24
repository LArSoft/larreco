/**
 *  @file   VoronoiDiagram.h
 * 
 *  @brief  Implements a VoronoiDiagram for use in clustering
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef VoronoiDiagram_h
#define VoronoiDiagram_h

// std includes
#include <list>
#include <algorithm>
#include <queue>

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/SweepEvent.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/BeachLine.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/EventUtilities.h"

#include <boost/polygon/voronoi.hpp>

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
//------------------------------------------------------------------------------------------------------------------------------------------

namespace voronoi2d
{
/**
 *  @brief  VoronoiDiagram class definiton
 */
class VoronoiDiagram
{
public:
    using PointPair       = std::pair<dcel2d::Point,dcel2d::Point>;
    using MinMaxPointPair = std::pair<PointPair,PointPair>;
    
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    VoronoiDiagram(dcel2d::HalfEdgeList&, dcel2d::VertexList&, dcel2d::FaceList&);

    /**
     *  @brief  Destructor
     */
    virtual ~VoronoiDiagram();
    
    /**
     *  @brief Given an input set of 2D points construct a 2D voronoi diagram
     *
     *  @param PointList The input list of 2D points
     */
    void buildVoronoiDiagram(const dcel2d::PointList&);
    
    void buildVoronoiDiagramBoost(const dcel2d::PointList&);
    
    /**
     *  @brief Recover the list of faces
     */
    const dcel2d::FaceList& getFaceList() const {return fFaceList;}
    
    /**
     *  @brief Recover the list of vertices
     */
    const dcel2d::VertexList& getVertexList() const {return fVertexList;}

    /**
     *  @brief recover the area of the convex hull
     */
    double getVoronoiDiagramArea() const {return fVoronoiDiagramArea;}
    
    /**
     *  @brief recover the point list representing the convex hull
     */
    const dcel2d::PointList& getConvexHull() const {return fConvexHullList;}
    
    /**
     *  @brief Given an input Point, find the nearest edge
     */
    PointPair getExtremePoints() const;

    /**
     *  @brief Given an input Point, find the nearest edge
     */
    PointPair findNearestEdge(const dcel2d::Point&, double&) const;
    
    /**
     *  @brief Given an input point, find the distance to the nearest edge/point
     */
    double findNearestDistance(const dcel2d::Point&) const;

private:
    
    using EventQueue = std::priority_queue<IEvent*, std::vector<IEvent*>, bool(*)(const IEvent*,const IEvent*)>;
    
    /**
     *  @brief There are two types of events in the queue, here we handle site events
     */
    void handleSiteEvents(BeachLine&, EventQueue&, IEvent*);
    
    /**
     *  @brief There are two types of events in the queue, here we handle circle events
     */
    void handleCircleEvents(BeachLine&, EventQueue&, IEvent*);
    
    void makeLeftCircleEvent(EventQueue&, BSTNode*, double);
    void makeRightCircleEvent(EventQueue&, BSTNode*, double);

    /**
     *  @brief There are two types of events in the queue, here we handle circle events
     */
    IEvent* makeCircleEvent(BSTNode*, BSTNode*, BSTNode*, double);
    
    bool computeCircleCenter(const dcel2d::Coords&, const dcel2d::Coords&, const dcel2d::Coords&, dcel2d::Coords&, double&, double&) const;

    bool computeCircleCenter2(const dcel2d::Coords&, const dcel2d::Coords&, const dcel2d::Coords&, dcel2d::Coords&, double&, double&) const;
    
    bool computeCircleCenter3(const dcel2d::Coords&, const dcel2d::Coords&, const dcel2d::Coords&, dcel2d::Coords&, double&, double&) const;
    
    /**
     *  @brief this function recovers the convex hull
     */
    void getConvexHull(const BSTNode*);
    
    /**
     *  @brief this aims to process remaining elements in the beachline after the event queue has been depleted
     */
    void terminateInfiniteEdges(BeachLine&, double);
    
    /**
     *  @brief Is a vertex inside the convex hull - meant to be a fast check
     */
    bool isInsideConvexHull(const dcel2d::Vertex&) const;
    
    /**
     *  @brief is vertex outside the convex hull and if so return some useful information
     */
    bool isOutsideConvexHull(const dcel2d::Vertex&, dcel2d::PointList::const_iterator, dcel2d::Coords&, double&) const;

    /**
     *  @brief Find the min/max values in x-y to use as a bounding box
     */
    void findBoundingBox(const dcel2d::VertexList&);
    
    /**
     * @brief Translate boost to dcel
     */
    using BoostEdgeToEdgeMap     = std::map<const boost::polygon::voronoi_edge<double>*,   dcel2d::HalfEdge*>;
    using BoostVertexToVertexMap = std::map<const boost::polygon::voronoi_vertex<double>*, dcel2d::Vertex*>;
    using BoostCellToFaceMap     = std::map<const boost::polygon::voronoi_cell<double>*,   dcel2d::Face*>;

    void boostTranslation(const dcel2d::PointList&,
                          const boost::polygon::voronoi_edge<double>*,
                          const boost::polygon::voronoi_edge<double>*,
                          BoostEdgeToEdgeMap&,
                          BoostVertexToVertexMap&,
                          BoostCellToFaceMap&);

    /**
     *  @brief merge degenerate vertices (found by zero length edges)
     */
    void mergeDegenerateVertices();
    
    /**
     *  @brief Compute the area of the faces
     */
    double ComputeFaceArea();
    
    /**
     *  @brief Gets the cross product of line from p0 to p1 and p0 to p2
     */
    double crossProduct(const dcel2d::Point& p0, const dcel2d::Point& p1, const dcel2d::Point& p2) const;

    /**
     *  @brief Computes the area of the enclosed convext hull
     */
    double Area() const;

    /**
     *  @brief Determines whether a point is to the left, right or on line specifiec by p0 and p1
     */
    bool isLeft(const dcel2d::Point& p0, const dcel2d::Point& p1, const dcel2d::Point& pCheck) const;
    
    dcel2d::HalfEdgeList& fHalfEdgeList;
    dcel2d::VertexList&   fVertexList;
    dcel2d::FaceList&     fFaceList;

    dcel2d::PointList     fPointList;
    SiteEventList         fSiteEventList;       //< Container for site events
    CircleEventList       fCircleEventList;     //< Container for circle events
    BSTNodeList           fCircleNodeList;      //< Container for the circle "nodes"
    
    dcel2d::PointList     fConvexHullList;      //< Points representing the convex hull
    dcel2d::Coords        fConvexHullCenter;    //< Center of the convex hull
    
    double                fXMin;                //< Bounding box min value x
    double                fXMax;                //< Bounding box max value x
    double                fYMin;                //< Bounding box min value y
    double                fYMax;                //< Bounding box max value y
    mutable int           fNumBadCircles;       //< Number bad circles
    double                fVoronoiDiagramArea;
    
    EventUtilities        fUtilities;           //< Handy functions to operate on arcs

};
    
} // namespace lar_cluster3d
#endif

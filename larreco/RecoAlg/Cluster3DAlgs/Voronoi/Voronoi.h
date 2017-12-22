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

// LArSoft includes
#include "lardata/RecoObjects/Cluster3D.h"

// std includes
#include <list>
#include <algorithm>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief  VoronoiDiagram class definiton
 */
class VoronoiDiagram
{
public:
    /**
     *  @brief Definitions used by the VoronoiDiagram algorithm
     */
    using Point           = std::tuple<float,float,const reco::ClusterHit3D*>;
    using PointList       = std::list<Point>;
    using PointPair       = std::pair<Point,Point>;
    using MinMaxPointPair = std::pair<PointPair,PointPair>;
    
    /**
     *  @brief Internal class definition to facilitate construction of diagram
     */
    class Event : public Point
    {
    public:
        Event(const Point& point) : Point(point) {m_valid = true;}
        ~Event() {}
        
        void setInvalid() {m_valid = false;}
        
        bool isSite()   const {return std::get<2>(*this);}
        bool isCircle() const {return !std::get<2>(*this);}
        bool isValid()  const {return m_valid;}
        
        bool operator<(const Event& right) const {return std::get<0>(*this) < std::get<0>(right);}
    private:
        mutable bool m_valid;
    };
    
    /**
     *  @brief internal class to define an arc
     */
    class Arc : public Event
    {
    public:
        Arc(const Event* event) : Event(*event) {}
        ~Arc() {}
        
        // Define the sorting operator for arcs
        bool operator<(const Arc& right) const {return std::get<1>(*this) < std::get<1>(right);}
        float operator()(float x, float l) const;
    private:
        Event* m_circle;
    };
    
    using ArcVec = std::vector<Arc>;
    
    /**
     * @brief this will describe break points
     */
    class BreakPoint
    {
    public:
        BreakPoint(Arc* left, Arc* right) : m_left(left), m_right(right) {}
        ~BreakPoint() {}
        
        Arc* getLeftArc()  const {return m_left;}
        Arc* getRightArc() const {return m_right;}
        
        float*       getBreakPoint(const Event*);
        const float* getBreakPoint() const {return m_breakpoint;}
        
        void setLeftArc(Arc* arc)  {m_left = arc;}
        void setRightArc(Arc* arc) {m_right = arc;}
        
    private:
        Arc*  m_left;
        Arc*  m_right;
        float m_breakpoint[2];
    };
    
    using BeachLine = std::list<BreakPoint>;

    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    VoronoiDiagram(const PointList&);

    /**
     *  @brief  Destructor
     */
    virtual ~VoronoiDiagram();
    
    /**
     *  @brief recover the list of points used to build convex hull
     */
    const PointList& getPointsList() {return fPoints;}
    
    /**
     *  @brief recover the list of convex hull vertices
     */
    const PointList& getVoronoiDiagram() const {return fVoronoiDiagram;}
    
    /**
     *  @brief find the ends of the convex hull (along its x axis)
     */
    const MinMaxPointPair& getMinMaxPointPair() const {return fMinMaxPointPair;}
    
    /**
     *  @brief recover the area of the convex hull
     */
    float getVoronoiDiagramArea() const {return fVoronoiDiagramArea;}
    
    /**
     *  @brief Given an input Point, find the nearest edge
     */
    PointPair findNearestEdge(const Point&, float&) const;
    
    /**
     *  @brief Given an input point, find the distance to the nearest edge/point
     */
    float findNearestDistance(const Point&) const;

private:
    
    /**
     *  @brief Given an input set of 2D points build a convex hull around them
     *
     *  @param PointList           The input list of 2D points to build hull around
     */
    void buildVoronoiDiagram(const PointList&);
    
    /**
     * @brief Find the best matching arc given a site point
     */
    ArcVec::iterator FindBestArc(ArcVec&, const Event*, float*) const;

    /**
     * @brief Find the nearest breakpoint
     */
    BeachLine::iterator FindBestArc(BeachLine&, const Event*, float*) const;

    /**
     *  @brief Gets the cross product of line from p0 to p1 and p0 to p2
     */
    float crossProduct(const Point& p0, const Point& p1, const Point& p2) const;

    /**
     *  @brief Computes the area of the enclosed convext hull
     */
    float Area() const;

    /**
     *  @brief Determines whether a point is to the left, right or on line specifiec by p0 and p1
     */
    bool isLeft(const Point& p0, const Point& p1, const Point& pCheck) const;

    const PointList& fPoints;
    PointList        fVoronoiDiagram;
    MinMaxPointPair  fMinMaxPointPair;
    float            fVoronoiDiagramArea;

};
    
} // namespace lar_cluster3d
#endif

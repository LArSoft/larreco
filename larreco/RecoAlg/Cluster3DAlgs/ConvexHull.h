/**
 *  @file   ConvexHull.h
 * 
 *  @brief  Implements a ConvexHull for use in clustering
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef ConvexHull_h
#define ConvexHull_h

// Algorithm includes
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

// std includes
#include <list>
#include <algorithm>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
/**
 *  @brief  ConvexHull class definiton
 */
class ConvexHull
{
public:
    /**
     *  @brief Definitions used by the ConvexHull algorithm
     */
    using Point           = std::tuple<float,float,const reco::ClusterHit3D*>;
    using PointList       = std::list<Point>;
    using PointPair       = std::pair<Point,Point>;
    using MinMaxPointPair = std::pair<PointPair,PointPair>;
    
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    ConvexHull(const PointList&);

    /**
     *  @brief  Destructor
     */
    virtual ~ConvexHull();
    
    /**
     *  @brief recover the list of points used to build convex hull
     */
    const PointList& getPointsList() {return fPoints;}
    
    /**
     *  @brief recover the list of convex hull vertices
     */
    const PointList& getConvexHull() const {return fConvexHull;}
    
    /**
     *  @brief find the ends of the convex hull (along its x axis)
     */
    const MinMaxPointPair& getMinMaxPointPair() const {return fMinMaxPointPair;}
    
    /**
     *  @brief recover the area of the convex hull
     */
    float getConvexHullArea() const {return fConvexHullArea;}
    
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
    void getConvexHull(const PointList&);

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
    PointList        fConvexHull;
    MinMaxPointPair  fMinMaxPointPair;
    float            fConvexHullArea;

};
    
} // namespace lar_cluster3d
#endif

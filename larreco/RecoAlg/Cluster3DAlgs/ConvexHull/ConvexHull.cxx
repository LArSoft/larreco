/**
 *  @file   ConvexHull.cxx
 *
 *  @brief  Producer module to create 3D clusters from input hits
 *
 */

// Framework Includes

#include "larreco/RecoAlg/Cluster3DAlgs/ConvexHull/ConvexHull.h"

// LArSoft includes

// boost includes
#include <boost/range/adaptor/reversed.hpp>

// Eigen
#include <Eigen/Core>

// std includes
#include <cmath>
#include <limits>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

  ConvexHull::ConvexHull(const PointList& pointListIn, float kinkAngle, float minEdgeDistance)
    : fKinkAngle(kinkAngle)
    , fMinEdgeDistance(minEdgeDistance)
    , fPoints(pointListIn)
    , fConvexHullArea(0.)
  {
    fConvexHull.clear();

    // Build out the convex hull around the input points
    getConvexHull(fPoints);

    // And the area
    fConvexHullArea = Area();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  ConvexHull::~ConvexHull() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool ConvexHull::isLeft(const Point& p0, const Point& p1, const Point& pCheck) const
  {
    // Use the cross product to determine if the check point lies to the left, on or right
    // of the line defined by points p0 and p1
    return crossProduct(p0, p1, pCheck) > 0;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  float ConvexHull::crossProduct(const Point& p0, const Point& p1, const Point& p2) const
  {
    // Define a quick 2D cross product here since it will used quite a bit!
    float deltaX = std::get<0>(p1) - std::get<0>(p0);
    float deltaY = std::get<1>(p1) - std::get<1>(p0);
    float dCheckX = std::get<0>(p2) - std::get<0>(p0);
    float dCheckY = std::get<1>(p2) - std::get<1>(p0);

    return ((deltaX * dCheckY) - (deltaY * dCheckX));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  float ConvexHull::Area() const
  {
    float area(0.);

    // Compute the area by taking advantage of
    // 1) the ability to decompose a convex hull into triangles,
    // 2) the ability to use the cross product to calculate the area
    // So, the technique is to pick a point which we take as the center of the polygon
    // and then sum the signed area of triangles formed from this point to two adjecent
    // vertices on the convex hull.
    float x = 0.;
    float y = 0.;
    float n = float(fConvexHull.size());

    for (const auto& point : fConvexHull) {
      x += std::get<0>(point);
      y += std::get<1>(point);
    }

    Point center(x / n, y / n, 0);

    Point lastPoint = fConvexHull.front();

    for (const auto& point : fConvexHull) {
      if (point != lastPoint) area += 0.5 * crossProduct(center, lastPoint, point);

      lastPoint = point;
    }

    return area;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void ConvexHull::getConvexHull(const PointList& pointList)
  {
    // Start by identifying the min/max points
    Point pMinMin = pointList.front();

    // Find the point with maximum y sharing the same x value...
    PointList::const_iterator pMinMaxItr =
      std::find_if(pointList.begin(), pointList.end(), [pMinMin](const auto& elem) {
        return std::get<0>(elem) != std::get<0>(pMinMin);
      });

    Point pMinMax = *(--pMinMaxItr); // This could be / probably is the same point

    fMinMaxPointPair.first.first = pMinMin;
    fMinMaxPointPair.first.second = pMinMax;

    Point pMaxMax = pointList.back(); // get the last point

    // Find the point with minimum y sharing the same x value...
    PointList::const_reverse_iterator pMaxMinItr =
      std::find_if(pointList.rbegin(), pointList.rend(), [pMaxMax](const auto& elem) {
        return std::get<0>(elem) != std::get<0>(pMaxMax);
      });

    Point pMaxMin = *(--pMaxMinItr); // This could be / probably is the same point

    fMinMaxPointPair.second.first = pMaxMin;
    fMinMaxPointPair.second.second = pMaxMax;

    // Get the lower convex hull
    PointList lowerHullList;

    lowerHullList.push_back(pMinMin);

    // loop over points in the set
    for (auto curPoint : pointList) {
      // First check that we even want to consider this point
      if (isLeft(pMinMin, pMaxMin, curPoint)) continue;

      // In words: treat "lowerHullVec" as a stack. If there is only one entry then
      // push the current point onto the stack. If more than one entry then check if
      // the current point is to the left of the line from the last two points in the
      // stack. If it is then add the current point to the stack, if it is not then
      // pop the last point off the stack and try again.
      while (lowerHullList.size() > 1) {
        Point lastPoint = lowerHullList.back();

        lowerHullList.pop_back();

        Point nextToLastPoint = lowerHullList.back();

        if (isLeft(nextToLastPoint, lastPoint, curPoint)) {
          lowerHullList.push_back(lastPoint);
          break;
        }
      }

      lowerHullList.push_back(curPoint);
    }

    // Now get the upper hull
    PointList upperHullList;

    upperHullList.push_back(pMaxMax);

    for (auto curPoint : boost::adaptors::reverse(pointList)) {
      // First check that we even want to consider this point
      // Remember that we are going "backwards" so still want
      // curPoint to lie to the "left"
      if (isLeft(pMaxMax, pMinMax, curPoint)) continue;

      // Replicate the above but going the other direction...
      while (upperHullList.size() > 1) {
        Point lastPoint = upperHullList.back();

        upperHullList.pop_back();

        Point nextToLastPoint = upperHullList.back();

        if (isLeft(nextToLastPoint, lastPoint, curPoint)) {
          upperHullList.push_back(lastPoint);
          break;
        }
      }

      upperHullList.push_back(curPoint);
    }

    // Now we merge the two lists into the output list
    std::copy(lowerHullList.begin(), lowerHullList.end(), std::back_inserter(fConvexHull));

    if (pMaxMin == pMaxMax) upperHullList.pop_front();

    std::copy(upperHullList.begin(), upperHullList.end(), std::back_inserter(fConvexHull));

    if (pMinMin != pMinMax) fConvexHull.push_back(pMinMin);

    return;
  }

  const ConvexHull::PointList& ConvexHull::getExtremePoints()
  {
    PointList::const_iterator nextPointItr = fConvexHull.begin();
    PointList::const_iterator firstPointItr = nextPointItr++;

    float maxSeparation(0.);

    // Make sure the current list has been cleared
    fExtremePoints.clear();

    // For finding the two farthest points
    PointPair extremePoints(Point(0, 0, NULL), Point(0, 0, NULL));

    while (nextPointItr != fConvexHull.end()) {
      // Get a vector representing the edge from the first to the current next point
      Point firstPoint = *firstPointItr++;
      Point nextPoint = *nextPointItr;
      Eigen::Vector2f firstEdge(std::get<0>(*firstPointItr) - std::get<0>(firstPoint),
                                std::get<1>(*firstPointItr) - std::get<1>(firstPoint));

      // normalize it
      firstEdge.normalize();

      PointList::const_iterator endPointItr = nextPointItr;

      while (++endPointItr != fConvexHull.end()) {
        Point endPoint = *endPointItr;
        Eigen::Vector2f nextEdge(std::get<0>(endPoint) - std::get<0>(nextPoint),
                                 std::get<1>(endPoint) - std::get<1>(nextPoint));

        // normalize it
        nextEdge.normalize();

        // Have we found the turnaround point?
        if (firstEdge.dot(nextEdge) < 0.) {
          Eigen::Vector2f separation(std::get<0>(nextPoint) - std::get<0>(firstPoint),
                                     std::get<1>(nextPoint) - std::get<1>(firstPoint));
          float separationDistance = separation.norm();

          if (separationDistance > maxSeparation) {
            extremePoints.first = firstPoint;
            extremePoints.second = nextPoint;
            maxSeparation = separationDistance;
          }

          // Regardless of thise being the maximum distance we have hit a turnaround point so
          // we need to break out of this loop
          break;
        }

        nextPointItr = endPointItr;
        nextPoint = endPoint;
      }

      // If we have hit the end of the convex hull without finding a turnaround point then we are not
      // going to find one so break out of the main loop
      if (endPointItr == fConvexHull.end()) break;

      // Need to make sure we don't overrun the next point
      if (firstPointItr == nextPointItr) nextPointItr++;
    }

    fExtremePoints.push_back(extremePoints.first);
    fExtremePoints.push_back(extremePoints.second);

    return fExtremePoints;
  }

  const reco::ConvexHullKinkTupleList& ConvexHull::getKinkPoints()
  {
    // Goal here is to isolate the points where we see a large deviation in the contour defined by the
    // convex hull. The complications are numerous, chief is that some deviations can be artificially
    // "rounded" by a series of very small steps around a large corner. So we need to protect against that.

    // Make sure the current list has been cleared
    fKinkPoints.clear();

    // Need a mininum number of points/edges
    if (fConvexHull.size() > 3) {
      // Idea will be to traverse the convex hull and keep track of all points where there is a
      // "kink" which will be defined as a "large" angle between adjacent edges.
      // Recall that construxtion of the convex hull results in the same point at the start and
      // end of the list
      // Getting the initial point requires some contortions because the convex hull point list will
      // contain the same point at both ends of the list (why?)
      PointList::iterator pointItr = fConvexHull.begin();

      // Advance to the second to last element
      std::advance(pointItr, fConvexHull.size() - 2);

      Point lastPoint = *pointItr++;

      // Reset pointer to the first element
      pointItr = fConvexHull.begin();

      Point curPoint = *pointItr++;
      Eigen::Vector2f lastEdge(std::get<0>(curPoint) - std::get<0>(lastPoint),
                               std::get<1>(curPoint) - std::get<1>(lastPoint));

      lastEdge.normalize();

      while (pointItr != fConvexHull.end()) {
        Point& nextPoint = *pointItr++;

        Eigen::Vector2f nextEdge(std::get<0>(nextPoint) - std::get<0>(curPoint),
                                 std::get<1>(nextPoint) - std::get<1>(curPoint));

        if (nextEdge.norm() > fMinEdgeDistance) {
          nextEdge.normalize();

          if (lastEdge.dot(nextEdge) < fKinkAngle)
            fKinkPoints.emplace_back(curPoint, lastEdge, nextEdge);

          lastEdge = nextEdge;
        }

        curPoint = nextPoint;
      }
    }

    return fKinkPoints;
  }

  ConvexHull::PointPair ConvexHull::findNearestEdge(const Point& point,
                                                    float& closestDistance) const
  {
    // The idea is to find the nearest edge of the convex hull, defined by
    // two adjacent vertices of the hull, to the input point.
    // As near as I can tell, the best way to do this is to apply brute force...
    // Idea will be to iterate over pairs of points
    PointList::const_iterator curPointItr = fConvexHull.begin();
    Point prevPoint = *curPointItr++;
    Point curPoint = *curPointItr;

    // Set up the winner
    PointPair closestEdge(prevPoint, curPoint);

    closestDistance = std::numeric_limits<float>::max();

    // curPointItr is meant to point to the second point
    while (curPointItr != fConvexHull.end()) {
      if (curPoint != prevPoint) {
        // Dereference some stuff
        float xPrevToPoint = (std::get<0>(point) - std::get<0>(prevPoint));
        float yPrevToPoint = (std::get<1>(point) - std::get<1>(prevPoint));
        float xPrevToCur = (std::get<0>(curPoint) - std::get<0>(prevPoint));
        float yPrevToCur = (std::get<1>(curPoint) - std::get<1>(prevPoint));
        float edgeLength = std::sqrt(xPrevToCur * xPrevToCur + yPrevToCur * yPrevToCur);

        // Find projection onto convex hull edge
        float projection = ((xPrevToPoint * xPrevToCur) + (yPrevToPoint * yPrevToCur)) / edgeLength;

        // DOCA point
        Point docaPoint(std::get<0>(prevPoint) + projection * xPrevToCur / edgeLength,
                        std::get<1>(prevPoint) + projection * yPrevToCur / edgeLength,
                        0);

        if (projection > edgeLength) docaPoint = curPoint;
        if (projection < 0) docaPoint = prevPoint;

        float xDocaDist = std::get<0>(point) - std::get<0>(docaPoint);
        float yDocaDist = std::get<1>(point) - std::get<1>(docaPoint);
        float docaDist = xDocaDist * xDocaDist + yDocaDist * yDocaDist;

        if (docaDist < closestDistance) {
          closestEdge = PointPair(prevPoint, curPoint);
          closestDistance = docaDist;
        }
      }

      prevPoint = curPoint;
      curPoint = *curPointItr++;
    }

    closestDistance = std::sqrt(closestDistance);

    // Convention is convex hull vertices sorted in counter clockwise fashion so if the point
    // lays to the left of the nearest edge then it must be an interior point
    if (isLeft(closestEdge.first, closestEdge.second, point)) closestDistance = -closestDistance;

    return closestEdge;
  }

  float ConvexHull::findNearestDistance(const Point& point) const
  {
    float closestDistance;

    findNearestEdge(point, closestDistance);

    return closestDistance;
  }

} // namespace lar_cluster3d

/**
 *  @file   VoronoiDiagram.cxx
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes

#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/Voronoi.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/SweepEvent.h"

// LArSoft includes

// boost includes
#include <boost/range/adaptor/reversed.hpp>

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

VoronoiDiagram::VoronoiDiagram(const PointList& pointListIn) : fPoints(pointListIn), fVoronoiDiagramArea(0.)
{
    fVoronoiDiagram.clear();
    
    // Build out the convex hull around the input points
    buildVoronoiDiagram(fPoints);
    
    // And the area
    fVoronoiDiagramArea = Area();
}

//------------------------------------------------------------------------------------------------------------------------------------------

VoronoiDiagram::~VoronoiDiagram()
{
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
bool VoronoiDiagram::isLeft(const Point& p0, const Point& p1, const Point& pCheck) const
{
    // Use the cross product to determine if the check point lies to the left, on or right
    // of the line defined by points p0 and p1
    return crossProduct(p0, p1, pCheck) > 0;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

float VoronoiDiagram::crossProduct(const Point& p0, const Point& p1, const Point& p2) const
{
    // Define a quick 2D cross product here since it will used quite a bit!
    float deltaX  = std::get<0>(p1) - std::get<0>(p0);
    float deltaY  = std::get<1>(p1) - std::get<1>(p0);
    float dCheckX = std::get<0>(p2) - std::get<0>(p0);
    float dCheckY = std::get<1>(p2) - std::get<1>(p0);
    
    return ((deltaX * dCheckY) - (deltaY * dCheckX));
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
float VoronoiDiagram::Area() const
{
    float area(0.);
    
    // Compute the area by taking advantage of
    // 1) the ability to decompose a convex hull into triangles,
    // 2) the ability to use the cross product to calculate the area
    // So, the technique is to pick a point (for technical reasons we use 0,0)
    // and then sum the signed area of triangles formed from this point to two adjecent
    // vertices on the convex hull.
    Point center(0.,0.,0);
    Point lastPoint = fVoronoiDiagram.front();
    
    for(const auto& point : fVoronoiDiagram)
    {
        if (point != lastPoint) area += 0.5 * crossProduct(center,lastPoint,point);
        
        lastPoint = point;
    }
    
    return area;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool compareEventPtrs(const VoronoiDiagram::Event* left, const VoronoiDiagram::Event* right) {return (std::get<0>(*left) != std::get<0>(*right)) ? std::get<0>(*left) > std::get<0>(*right) : std::get<1>(*left) > std::get<1>(*right);}

//------------------------------------------------------------------------------------------------------------------------------------------

void VoronoiDiagram::buildVoronoiDiagram(const PointList& pointList)
{
    // Start by identifying the min/max points
    Point pMinMin = pointList.front();
    
    // Find the point with maximum y sharing the same x value...
    PointList::const_iterator pMinMaxItr = std::find_if(pointList.begin(),pointList.end(),[pMinMin](const auto& elem){return std::get<0>(elem) != std::get<0>(pMinMin);});
    
    Point pMinMax = *(--pMinMaxItr);  // This could be / probably is the same point
    
    fMinMaxPointPair.first.first  = pMinMin;
    fMinMaxPointPair.first.second = pMinMax;
    
    Point pMaxMax = pointList.back();  // get the last point
    
    // Find the point with minimum y sharing the same x value...
    PointList::const_reverse_iterator pMaxMinItr = std::find_if(pointList.rbegin(),pointList.rend(),[pMaxMax](const auto& elem){return std::get<0>(elem) != std::get<0>(pMaxMax);});
    
    Point pMaxMin = *(--pMaxMinItr);  // This could be / probably is the same point
    
    fMinMaxPointPair.second.first  = pMaxMin;
    fMinMaxPointPair.second.second = pMaxMax;
    
    // Build events and store in an stl vector
    // At the same time, build the priority_queue so we only do one loop here
    std::vector<Event> eventVec;
    bool(* fCompareEventPtrsPtr)(const Event* left, const Event* right) = compareEventPtrs;
    std::priority_queue<const Event*, std::vector<const Event*>, bool(*)(const VoronoiDiagram::Event* left, const VoronoiDiagram::Event* right)> eventQueue(fCompareEventPtrsPtr);
    
    // Make sure we don't reallocate the vector when filling
    eventVec.reserve(2*pointList.size());

    for(const auto& point : pointList)
    {
        eventVec.push_back(Event(point));
        eventQueue.push(&eventVec.back());
    }
    
    // Define the beach line and a container for the arcs
    BeachLine beachLine;
    ArcVec    arcVec;
    
    arcVec.reserve(2*pointList.size());
    
    // Seed it
    arcVec.push_back(Arc(eventQueue.top()));
    eventQueue.pop();

    beachLine.push_back(BreakPoint(0,&arcVec.back()));
    beachLine.push_back(BreakPoint(&arcVec.back(),0));
    
    // Advance the sweep line
    while(!eventQueue.empty())
    {
        const Event* event = eventQueue.top();
        
        // If a site event then we are adding a new arc to the beach line
        if (event->isSite())
        {
            float breakPoint[] = {0.,0.};
            
            // Search for the last break point less than the current site point
            BeachLine::iterator breakPointItr = FindBestArc(beachLine, event, breakPoint);
            
            // just checking
            std::cout << "===>Breakpoint: " << breakPoint[1] << ", event: " << std::get<1>(*event) << std::endl;
            
            // The arc to the "left" of this break point is the one to deal with here
            // Basically, we split that one into two and insert the new point, this will
            // necessarily create two new breakpoints
            arcVec.push_back(Arc(event));
            
            Arc* origLeftArc = breakPointItr->getLeftArc();
            Arc* newArc      = &arcVec.back();
            
            breakPointItr->setLeftArc(newArc);
            breakPointItr = beachLine.insert(breakPointItr, BreakPoint(origLeftArc,newArc));
        }
        
        eventQueue.pop();
    }
    
    return;
}
    
VoronoiDiagram::ArcVec::iterator VoronoiDiagram::FindBestArc(ArcVec& arcVec, const Event* site, float* arcPoint) const
{
    ArcVec::iterator bestItr  = arcVec.begin();
    float            bestDist = std::numeric_limits<float>::max();
    
    for(ArcVec::iterator arcItr = arcVec.begin(); arcItr != arcVec.end(); arcItr++)
    {
        float deltaX    = std::get<0>(*site) - std::get<0>(*arcItr);
        float deltaY    = std::get<1>(*site) - std::get<1>(*arcItr);
        float distToArc = 0.5 * (deltaX + deltaY * deltaY / deltaX);
        
        if (distToArc < bestDist)
        {
            bestItr  = arcItr;
            bestDist = distToArc;
        }
    }
    
    arcPoint[0] = std::get<0>(*bestItr);
    arcPoint[1] = std::get<1>(*bestItr);

    return bestItr;
}
    
VoronoiDiagram::BeachLine::iterator VoronoiDiagram::FindBestArc(BeachLine& beachLine, const Event* site, float* breakPoint) const
{
    BeachLine::iterator beachItr = beachLine.begin();
    
    std::cout << "** beachLine size: " << beachLine.size() << ", x,y: " << std::get<0>(*site) << ", " << std::get<1>(*site);
    
    while(beachItr != beachLine.end() && std::get<1>(*site) > (*beachItr).getBreakPoint(site)[1])
    {
        beachItr++;
    }
    
    breakPoint[0] = (*beachItr).getBreakPoint()[0];
    breakPoint[1] = (*beachItr).getBreakPoint()[1];
    
    std::cout << ", breakpoint: " << breakPoint[0] << ", " << breakPoint[1] << std::endl;
    
    return beachItr;
}
    
float VoronoiDiagram::Arc::operator()(float x, float l) const
{
    float posX  = std::get<1>(*this);
    float posY  = std::get<0>(*this);
    float denom = 2. * (posY - l);
    
    return (x * x - 2. * posX * x + posX * posX + posY * posY - l * l) / denom;
}

float* VoronoiDiagram::BreakPoint::getBreakPoint(const Event* site)
{
    // Initialize, this would be return if the right arc is missing
    m_breakpoint[0] = std::get<0>(*site);
    m_breakpoint[1] = std::numeric_limits<float>::max();
    
    // Can only compute a breakpoint if we have two arcs
    if (m_left && m_right)
    {
        // Note that the input site is used to give us the position of the sweep line
        // Start by getting the delta x values (remembering the sweep is in the x direction)
        float ly      = std::get<0>(*site);
        float deltaY1 = std::get<0>(*m_left)  - ly;
        float deltaY2 = std::get<0>(*m_right) - ly;
        
        // if the two are the same then the arcs are side-by-side and intersection is right in the middle
        if (abs(deltaY1 - deltaY2) < std::numeric_limits<float>::epsilon())
        {
            m_breakpoint[1] = 0.5 * (std::get<1>(*m_right) - std::get<1>(*m_left));
            m_breakpoint[0] = (*m_left)(m_breakpoint[1], ly);
        }
        // otherwise, we do the full calculation
        else
        {
            // set up for quadratic equation
            float p1x = std::get<1>(*m_left);
            float p1y = std::get<0>(*m_left);
            float p2x = std::get<1>(*m_right);
            float p2y = std::get<0>(*m_right);
            float a   = deltaY2 - deltaY1;
            float b   = 2. * (p2x * deltaY1 - p1x * deltaY2);
            float c   = deltaY2 * (p1x * p1x + p1y * p1y - ly * ly) - deltaY1 * (p2x * p2x + p2y * p2y - ly * ly);
            
            float radical = b * b - 4. * a * c;
            
            if (radical < 0.)
            {
                std::cout << "This is a problem... radical: " << radical << ", a: " << a << ", b: " << b << ", c: " << c << std::endl;
                radical = 0.;
            }
            
            radical = sqrt(radical);
            
            float xIntersect  = 0.5 * (-b + radical) / a;
            
            if (xIntersect < p1x || xIntersect > p2x) xIntersect = 0.5 * (-b - radical) / a;

            m_breakpoint[0] = (*m_left)(xIntersect, ly);
            m_breakpoint[1] = xIntersect;
        }
    }
    // Watch for missing left arc, need to set break point to negative infinity
    else if (!m_left) m_breakpoint[1] = std::numeric_limits<float>::lowest();
    
    return m_breakpoint;
}

VoronoiDiagram::PointPair VoronoiDiagram::findNearestEdge(const Point& point, float& closestDistance) const
{
    // The idea is to find the nearest edge of the convex hull, defined by
    // two adjacent vertices of the hull, to the input point.
    // As near as I can tell, the best way to do this is to apply brute force...
    // Idea will be to iterate over pairs of points
    PointList::const_iterator curPointItr = fVoronoiDiagram.begin();
    Point                     prevPoint   = *curPointItr++;
    Point                     curPoint    = *curPointItr;
    
    // Set up the winner
    PointPair closestEdge(prevPoint,curPoint);
    
    closestDistance = std::numeric_limits<float>::max();
    
    // curPointItr is meant to point to the second point
    while(curPointItr != fVoronoiDiagram.end())
    {
        if (curPoint != prevPoint)
        {
            // Dereference some stuff
            float xPrevToPoint = (std::get<0>(point)    - std::get<0>(prevPoint));
            float yPrevToPoint = (std::get<1>(point)    - std::get<1>(prevPoint));
            float xPrevToCur   = (std::get<0>(curPoint) - std::get<0>(prevPoint));
            float yPrevToCur   = (std::get<1>(curPoint) - std::get<1>(prevPoint));
            float edgeLength   = std::sqrt(xPrevToCur * xPrevToCur + yPrevToCur * yPrevToCur);

            // Find projection onto convex hull edge
            float projection = ((xPrevToPoint * xPrevToCur) + (yPrevToPoint * yPrevToCur)) / edgeLength;
            
            // DOCA point
            Point docaPoint(std::get<0>(prevPoint) + projection * xPrevToCur / edgeLength,
                            std::get<1>(prevPoint) + projection * yPrevToCur / edgeLength, 0);
            
            if (projection > edgeLength) docaPoint = curPoint;
            if (projection < 0)          docaPoint = prevPoint;
            
            float xDocaDist = std::get<0>(point) - std::get<0>(docaPoint);
            float yDocaDist = std::get<1>(point) - std::get<1>(docaPoint);
            float docaDist  = xDocaDist * xDocaDist + yDocaDist * yDocaDist;
            
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

float VoronoiDiagram::findNearestDistance(const Point& point) const
{
    float     closestDistance;
    
    findNearestEdge(point,closestDistance);
    
    return closestDistance;
}

} // namespace lar_cluster3d

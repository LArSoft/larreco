/**
 *  @file   Event.cxx
 * 
 *  @brief  Producer module to create 3D clusters from input hits
 * 
 */

// Framework Includes

#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/SweepEvent.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {
    
bool SiteEvent::newSiteToLeft(const IEvent* newSite, const IEvent* leftArc, const IEvent* rightArc) const
{
    // Note that the input site is used to give us the position of the sweep line
    // Start by getting the delta x values (remembering the sweep is in the x direction)
    float ly         = newSite->beachLinePos();
    float deltaY1    = leftArc->beachLinePos()  - ly;
    float deltaY2    = rightArc->beachLinePos() - ly;
    float breakPoint = -std::numeric_limits<float>::max();
    
    // if the two are the same then the arcs are side-by-side and intersection is right in the middle
    if (abs(deltaY1 - deltaY2) < std::numeric_limits<float>::epsilon()) breakPoint = 0.5 * (rightArc->xPos() - leftArc->xPos());
    
    // otherwise, we do the full calculation
    else
    {
        // set up for quadratic equation
        float p1x = leftArc->xPos();
        float p1y = leftArc->beachLinePos();
        float p2x = rightArc->xPos();
        float p2y = rightArc->beachLinePos();
        float a   = deltaY2 - deltaY1;
        float b   = 2. * (p2x * deltaY1 - p1x * deltaY2);
        float c   = deltaY2 * (p1x * p1x + p1y * p1y - ly * ly) - deltaY1 * (p2x * p2x + p2y * p2y - ly * ly);
        
        float radical = b * b - 4. * a * c;
        
        if (radical < 0.)
        {
            std::cout << "***********************************************************" << std::endl;
            std::cout << "This is a problem... radical: " << radical << ", a: " << a << ", b: " << b << ", c: " << c << std::endl;
            std::cout << "***********************************************************" << std::endl;
            radical = 0.;
        }
        
        radical = sqrt(radical);
        
        float xIntersect  = 0.5 * (-b + radical) / a;
        
        if (xIntersect < p1x || xIntersect > p2x) xIntersect = 0.5 * (-b - radical) / a;
        
        breakPoint = xIntersect;
    }
    
    return newSite->xPos() < breakPoint;
}

} // namespace lar_cluster3d

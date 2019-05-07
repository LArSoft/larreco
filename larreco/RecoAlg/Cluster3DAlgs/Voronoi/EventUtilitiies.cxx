/**
 *  @file   Event.cxx
 *
 *  @brief  Producer module to create 3D clusters from input hits
 *
 */

// Framework Includes

#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/EventUtilities.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/BeachLine.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"

// std includes
#include <functional>
#include <memory>
#include <queue>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace voronoi2d {

double EventUtilities::computeArcVal(double beachPos, double yPos, const IEvent* arc) const
{
    // Note that if the input arc site point lies on the beach line then the arc position is infinite
    double arcVal   = std::numeric_limits<double>::max();
    double deltaxLx = arc->xPos() - beachPos;

    if (std::abs(deltaxLx) > std::numeric_limits<double>::epsilon())
    {
        double deltayPy = yPos - arc->yPos();
        double sumPxLx  = arc->xPos() + beachPos;

        arcVal = 0.5 * (deltayPy * deltayPy / deltaxLx + sumPxLx);
    }

    return arcVal;
}

double EventUtilities::computeBreak(const double beachLinePos, const IEvent* leftArc, const IEvent* rightArc, RootsPair& roots) const
{
    // Given arcs to the left and right of this node (meaning we are a breakpoint), compute the
    // current coordinates of the breakpoint based on the input beachline position
    double lx         = beachLinePos;
    double deltaX1    = leftArc->xPos()  - lx;
    double deltaX2    = rightArc->xPos() - lx;
    double breakPoint = -std::numeric_limits<double>::max();

    // if the two are the same then the arcs are side-by-side and intersection is right in the middle
    if (std::abs(deltaX1 - deltaX2) < std::numeric_limits<double>::epsilon()) breakPoint = 0.5 * (rightArc->yPos() + leftArc->yPos());

    // otherwise, we do the full calculation
    else
    {
        // set up for quadratic equation
        double p1x     = leftArc->xPos();
        double p1y     = leftArc->yPos();
        double p2x     = rightArc->xPos();
        double p2y     = rightArc->yPos();
        double a       = p2x - p1x;
        double b       = 2. * (p2y * deltaX1 - p1y * deltaX2);
        double c       = deltaX2 * (p1y * p1y + deltaX1 * (p1x + lx)) - deltaX1 * (p2y * p2y + deltaX2 * (p2x + lx));
        double radical = std::max(0.,b * b - 4. * a * c);

        if (radical > 0.) radical = sqrt(radical);

        double rootPos = 0.5 * (-b + radical) / a;
        double rootNeg = 0.5 * (-b - radical) / a;

        roots.first  = std::min(rootPos, rootNeg);
        roots.second = std::max(rootPos, rootNeg);

        // Ah yes, the eternal question... which solution?
        // Generally, we think we want the solution which is "between" the left and the right
        // However, it can be that the right arc is the remnant of an arc whose center lies to the
        // left of the center of the left arc, and vice versa.
        if (p1x < p2x) breakPoint = roots.second;
        else           breakPoint = roots.first;
    }

    return breakPoint;
}

bool EventUtilities::newSiteToLeft(const IEvent* newSite, const IEvent* leftArc, const IEvent* rightArc) const
{
    // Note that the input site is used to give us the position of the sweep line
    // Using the coordinates of the beach line then recover the current breakpoint between the two
    // input arcs
    RootsPair roots;

    double breakPoint = computeBreak(newSite->xPos()-0.000001, leftArc, rightArc, roots);

    return newSite->yPos() < breakPoint;
}

} // namespace lar_cluster3d

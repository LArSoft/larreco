/**
 *  @file   EventUtilities.h
 * 
 *  @brief  Provides some basic functions operating in IEvent class objects
 *
 *  @author usher@slac.stanford.edu
 * 
 */
#ifndef EventUtilities_voronoi2d_h
#define EventUtilities_voronoi2d_h

// Get the beach line definitions
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/IEvent.h"

// std includes
#include <vector>
#include <list>
#include <algorithm>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace voronoi2d
{
using RootsPair = std::pair<double,double>;
    
/**
 *  @brief Internal class definitions to facilitate construction of diagram
 */
class EventUtilities
{
public:
    EventUtilities() {}
    ~EventUtilities() {}

    double computeArcVal(const double, const double, const IEvent*)             const;
    double computeBreak(const double, const IEvent*, const IEvent*, RootsPair&) const;
    bool   newSiteToLeft(const IEvent*, const IEvent*, const IEvent*)           const;
private:
};
    
} // namespace lar_cluster3d
#endif

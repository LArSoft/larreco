/**
 *  @file   DCEL2D.h
 *
 *  @brief  Definitions for a doubly connected edge list
 *          This will define a vertex, half edge and face
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IEvent_voronoi2d_h
#define IEvent_voronoi2d_h

#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/DCEL.h"

// std includes
#include <vector>
#include <algorithm>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace voronoi2d
{
/**
 *  @brief Define a virtual interface to the beachline "events"
 */
class BSTNode;

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

    virtual void                  setInvalid()         const = 0;

    virtual bool                  isSite()             const = 0;
    virtual bool                  isCircle()           const = 0;
    virtual bool                  isValid()            const = 0;

    virtual const dcel2d::Point&  getPoint()           const = 0;
    virtual double                xPos()               const = 0;
    virtual double                yPos()               const = 0;
    virtual const dcel2d::Coords& getCoords()          const = 0;
    virtual const dcel2d::Coords& circleCenter()       const = 0;

    virtual BSTNode*              getBSTNode()         const = 0;
    virtual void                  setBSTNode(BSTNode*)       = 0;

    virtual bool operator<(const IEvent& right) const = 0;
};

} // namespace lar_cluster3d
#endif

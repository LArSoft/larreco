/**
 *  @file   Event.h
 *
 *  @brief  Represents the Event implemented as a self balancing binary search tree
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef Event_h
#define Event_h

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/DCEL.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/IEvent.h"
namespace voronoi2d { class BSTNode; }

// std includes
#include <list>
#include <tuple>

// Eigen includes
#include <Eigen/Core>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace voronoi2d
{
/**
 *  @brief Internal class definitions to facilitate construction of diagram
 */
class SiteEvent : public dcel2d::Point, virtual public IEvent
{
    /**
     *  @brief This defines "Site" events which are generated from the
     *         input points. This implements the "IEvent" interface
     */
public:
    SiteEvent(const dcel2d::Point& point) : dcel2d::Point(point), m_valid(true), m_node(NULL)
    {
        m_coords = dcel2d::Coords(std::get<0>(point),std::get<1>(point),0.);
    }

    void                  setInvalid()              const override {m_valid = false;}
    void                  setBSTNode(BSTNode* node)       override {m_node = node;}

    bool                  isSite()                  const override {return std::get<2>(*this) != NULL;}
    bool                  isCircle()                const override {return std::get<2>(*this) == NULL;}
    bool                  isValid()                 const override {return m_valid;}
    const dcel2d::Point&  getPoint()                const override {return *this;}
    double                xPos()                    const override {return m_coords[0];}
    double                yPos()                    const override {return m_coords[1];}
    const dcel2d::Coords& getCoords()               const override {return m_coords;}
    const dcel2d::Coords& circleCenter()            const override {return m_coords;}
    BSTNode*              getBSTNode()              const override {return m_node;}

    bool operator<(const IEvent& right) const override {return xPos() < right.xPos();}
private:
    dcel2d::Coords m_coords;
    mutable bool   m_valid;
    BSTNode*       m_node;
};

class CircleEvent : public dcel2d::Point, virtual public IEvent
{
    /**
     *  @brief This defines "Circle" events which are generated during the
     *         "sweep" of the beach line and define when an arc on the
     *         beach line will disappear
     */
public:
    CircleEvent(const dcel2d::Point& point, const dcel2d::Coords& center) :
        dcel2d::Point(point),
        m_valid(true),
        m_node(NULL)
    {
        m_circleCenter = center;
    }
    ~CircleEvent() {}

    void                  setInvalid()              const override {m_valid = false;}
    void                  setBSTNode(BSTNode* node)       override {m_node = node;}

    bool                  isSite()                  const override {return std::get<2>(*this) != NULL;}
    bool                  isCircle()                const override {return std::get<2>(*this) == NULL;}
    bool                  isValid()                 const override {return m_valid;}
    const dcel2d::Point&  getPoint()                const override {return *this;}
    double                xPos()                    const override {return std::get<0>(*this);}
    double                yPos()                    const override {return std::get<1>(*this);}
    const dcel2d::Coords& getCoords()               const override {return m_circleCenter;}
    const dcel2d::Coords& circleCenter()            const override {return m_circleCenter;}
    BSTNode*              getBSTNode()              const override {return m_node;}

    bool operator<(const IEvent& right) const override {return xPos() < right.xPos();}
private:
    dcel2d::Coords m_circleCenter;
    mutable bool   m_valid;
    BSTNode*       m_node;
};

using SiteEventList   = std::list<SiteEvent>;
using CircleEventList = std::list<CircleEvent>;

} // namespace lar_cluster3d
#endif

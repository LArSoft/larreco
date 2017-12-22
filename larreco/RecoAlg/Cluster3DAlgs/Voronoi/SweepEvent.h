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

// Get the beach line definitions
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/BeachLine.h"

// LArSoft includes
#include "lardata/RecoObjects/Cluster3D.h"

// std includes
#include <vector>
#include <list>
#include <algorithm>
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
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
class SiteEvent : public Point, virtual public IEvent
{
public:
    SiteEvent(const Point& point) : Point(point) {m_valid = true;}
    ~SiteEvent() {}
    
    void setInvalid() override {m_valid = false;}
    
    bool  isSite()       const override {return std::get<2>(*this);}
    bool  isCircle()     const override {return !std::get<2>(*this);}
    bool  isValid()      const override {return m_valid;}
    float beachLinePos() const override {return std::get<0>(*this);}
    float xPos()         const override {return std::get<1>(*this);}

    bool operator<(const IEvent& right) const override {return beachLinePos() < right.beachLinePos();}
    
    bool newSiteToLeft(const IEvent*, const IEvent*, const IEvent*) const override;
private:
    mutable bool m_valid;
};
    
} // namespace lar_cluster3d
#endif

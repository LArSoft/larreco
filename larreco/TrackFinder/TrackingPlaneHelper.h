#ifndef TRACKINGPLANEHELPER_H
#define TRACKINGPLANEHELPER_H

#include "lardataobj/RecoBase/TrackingPlane.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardata/RecoObjects/SurfWireX.h"
#include "larcore/Geometry/WireGeo.h"
#include "TMath.h"

namespace recob {
  namespace tracking {

    inline Plane makePlane(recob::tracking::Point_t const& pos, recob::tracking::Vector_t const& dir) { return Plane(pos, dir); }

    Plane makePlane(recob::Trajectory::TrajectoryPoint_t const& s) { return Plane(s.position, s.direction()); }

    Plane makePlane(trkf::SurfWireX const& s) { return Plane(Point_t(s.x0(),s.y0(),s.z0()), Vector_t(0,-std::sin(s.phi()),std::cos(s.phi()))); }

    Plane makePlane(geo::WireGeo const& wgeom) { 
      double xyz[3] = {0.};
      wgeom.GetCenter(xyz); 
      double phi = TMath::PiOver2() - wgeom.ThetaZ();
      return Plane(Point_t(0.,xyz[1], xyz[2]), Vector_t(0,-std::sin(phi),std::cos(phi)));
    }

  }
}

#endif

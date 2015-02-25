//////////////////////////////////////////////////////////////////////
///
/// TrackTrajectoryAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm fitting a 3D trajectory through a set of hit pairs
///
////////////////////////////////////////////////////////////////////////
#ifndef TRACKTRAJECTORYALG_H
#define TRACKTRAJECTORYALG_H

#include <math.h>
#include <algorithm>
#include <vector>

#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoAlg/TrackLineFitAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "TVectorD.h"

namespace trkf {

  class TrackTrajectoryAlg {
    public:

    TrackTrajectoryAlg();
    
    virtual ~TrackTrajectoryAlg();
    
    void TrackTrajectory(
      std::array<std::vector<std::pair<double, geo::WireID>>,3>& trkXW,
      std::vector<TVector3>& TrajPos, std::vector<TVector3>& TrajDir,
      std::array<std::vector<double>,3> trkChg = std::array<std::vector<double>,3>());
    
    private:

    art::ServiceHandle<geo::Geometry> geom;

    TrackLineFitAlg fTrackLineFitAlg;
    
  }; // class TrackTrajectoryAlg

} // namespace trkf

#endif // ifndef TRACKTRAJECTORYALG_H

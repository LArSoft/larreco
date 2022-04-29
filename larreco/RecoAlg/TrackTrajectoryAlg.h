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

#include <array>
#include <vector>

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/TrackLineFitAlg.h"
namespace geo { struct WireID; }

#include "TVector3.h"

namespace trkf {

  class TrackTrajectoryAlg {
    public:

    void TrackTrajectory(std::array<std::vector<geo::WireID>,3> trkWID,
                         std::array<std::vector<double>,3> trkX,
                         std::array<std::vector<double>,3> trkXErr,
                         std::vector<TVector3>& TrajPos, std::vector<TVector3>& TrajDir);

    private:

    art::ServiceHandle<geo::Geometry const> geom;

    double minX;
    unsigned short minXPln;
    double maxX;
    unsigned short maxXPln;
    bool prt;

    unsigned short fMaxTrajPoints;  // maximum number of trajectory points
    double fHitWidthFactor;         // scales the number of trajectory points to the hit rms


    TrackLineFitAlg fTrackLineFitAlg;

    void ShortTrackTrajectory(std::array<std::vector<geo::WireID>,3> trkWID,
                              std::array<std::vector<double>,3> trkX,
                              std::array<std::vector<double>,3> trkXErr,
                              std::vector<TVector3>& TrajPos, std::vector<TVector3>& TrajDir);

  }; // class TrackTrajectoryAlg

} // namespace trkf

#endif // ifndef TRACKTRAJECTORYALG_H

//////////////////////////////////////////////////////////////////////
///
/// TrackLineFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 3D line given a number of points in 3 wire planes
///
////////////////////////////////////////////////////////////////////////
#ifndef TRACKLINEFITALG_H
#define TRACKLINEFITALG_H

#include <vector>

#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/WireReadout.h"

namespace geo {
  struct WireID;
}

class TVector3;

namespace trkf {

  class TrackLineFitAlg {
  public:
    void TrkLineFit(std::vector<geo::WireID>& hitWID,
                    std::vector<double>& hitX,
                    std::vector<double>& hitXErr,
                    double XOrigin,
                    TVector3& Pos,
                    TVector3& Dir,
                    float& ChiDOF) const;

  private:
    geo::WireReadoutGeom const* wireReadoutGeom = &art::ServiceHandle<geo::WireReadout>()->Get();

  }; // class TrackLineFitAlg

} // namespace trkf

#endif // ifndef TRACKLINEFITALG_H

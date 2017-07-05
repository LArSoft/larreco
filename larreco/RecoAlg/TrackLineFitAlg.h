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

#include <math.h>
#include <algorithm>
#include <vector>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// ROOT includes
#include "TMath.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TMatrixD.h"

namespace trkf {

  class TrackLineFitAlg {
    public:

    TrackLineFitAlg();
    
    virtual ~TrackLineFitAlg();
    
    void TrkLineFit(std::vector<geo::WireID>& hitWID, std::vector<double>& hitX, std::vector<double>& hitXErr,
                    double XOrigin, TVector3& Pos, TVector3& Dir, float& ChiDOF);
    
    private:

    art::ServiceHandle<geo::Geometry> geom;
    
  }; // class TrackLineFitAlg

} // namespace trkf

#endif // ifndef TRACKLINEFITALG_H

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTracking
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), June 2016
//
// Single track reconstruction toolkit based on Projection Matching Algorithm. Uses cluster collections
// to find single tracks (modularized version of our original code).
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgTracking_h
#define PMAlgTracking_h


// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"

// ROOT & C++
#include <memory>

namespace pma
{
	recob::Track convertFrom(const pma::Track3D& src, unsigned int tidx);
}

#endif

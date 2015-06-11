////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Projection Matching Algorithm
// see RecoAlg/PMAlg/PmaTrack3D.h for more details.
//
//      Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
//      Based on the "Precise 3D track reco..." AHEP (2013) 260820, with all the tricks that we
//      developed later and with the work for the full-event topology optimization that is still
//      under construction.
//
// Progress:
//    May-June 2015:  track finding and validation, quite a conservative iterative merging
//                    of matching clusters and growing tracks, no attempts to consciously
//                    build multi-track structures yet, however:
//                    cosmic tracking works fine since they are sets of independent tracks
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ProjectionMatchingAlg_h
#define ProjectionMatchingAlg_h

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "RecoAlg/PMAlg/Utilities.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"

// ROOT & C++
#include <memory>

namespace pma
{
	class ProjectionMatchingAlg;
}

class pma::ProjectionMatchingAlg
{
public:

	ProjectionMatchingAlg(const fhicl::ParameterSet& pset);
	virtual ~ProjectionMatchingAlg(void);

	void reconfigure(const fhicl::ParameterSet& p);

private:

  // Parameters used in the algorithm

  // Geometry and detector properties
  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

};

#endif

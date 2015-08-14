////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgVertexing
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), August 2015
//
// 3D vertex finding for Projection Matching Algorithm
//
//      Uses collection of pma::Track3D to find vertex candidates, then joins track in these points
//      and reoptimizes full structure of tracks.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgVertexing_h
#define PMAlgVertexing_h

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
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "RecoAlg/PMAlg/PmaTrack3D.h"

// ROOT & C++
#include <memory>

namespace pma
{
	class PMAlgVertexing;
}

class pma::PMAlgVertexing
{
public:

	PMAlgVertexing(const fhicl::ParameterSet& pset);
	virtual ~PMAlgVertexing(void);

	void reconfigure(const fhicl::ParameterSet& pset);

	/// Copy input tracks, find 3D vertices, connect tracks, break them or flip if needed,
	/// reoptimize track structures. Result is available as a collection of new tracks, vertices
	/// should be found as track interconnections (use GetVertices function).
	/// Optionally a collection of already known 3D vertex candidates may be provided.
	/// Return value is the number of vertices.
	size_t run(const std::vector< pma::Track3D* >& trk_input,
	           const std::vector< TVector3 >& vtx_input = std::vector< TVector3 >());

	const std::vector< pma::Track3D* >& getTracks(void) const { return fOutTracks; }

	std::vector< recob::Vertex > getVertices(void) const;

private:

	std::vector< pma::Track3D* > fOutTracks;
	void cleanTracks(void);

	// Parameters used in the algorithm

	double fInputVtxDist2D; // use vtx given at input if dist. [cm] to track in all 2D projections is below this max. value
	double fInputVtxDistY;  // use vtx given at input if dist. [cm] to track in 3D-Y is below this max. value
};

#endif

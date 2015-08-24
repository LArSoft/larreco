////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgVertexing
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), August 2015
//
// 3D vertex finding for Projection Matching Algorithm
//
//      Uses collection of pma::Track3D to find vertex candidates, then joins tracks in these points
//      and reoptimizes full structure of tracks.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgVertexing_h
#define PMAlgVertexing_h

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
	virtual ~PMAlgVertexing(void); // delete last produced tracks

	void reconfigure(const fhicl::ParameterSet& pset);

	/// Copy input tracks, find 3D vertices, connect tracks, break them or flip if needed,
	/// reoptimize track structures. Result is available as a collection of new tracks, vertices
	/// should be found as track interconnections (use GetVertices function).
	size_t run(const std::vector< pma::Track3D* >& trk_input);

	/// Copy input tracks, use provided 3D vertices to connect tracks, break tracks or flip if
	/// needed, reoptimize track structures. Result is available as a collection of new tracks,
	/// input vertices that were actually associated to tracks are copied to the output
	/// collection (use GetVertices function).
	size_t run(const std::vector< pma::Track3D* >& trk_input,
	           const std::vector< TVector3 >& vtx_input);

	const std::vector< pma::Track3D* >& getTracks(void) const { return fOutTracks; }

	const std::vector< TVector3 >& getVertices(void) const { return fOutVertices; }

private:
	bool findOneVtx(void);

	std::vector< TVector3 > fOutVertices;
	std::vector< pma::Track3D* > fOutTracks;
	std::vector< pma::Track3D* > fShortTracks;
	std::vector< pma::Track3D* > fEmTracks;
	void cleanTracks(void);

	// Parameters used in the algorithm

	double fInputVtxDist2D; // use vtx given at input if dist. [cm] to track in all 2D projections is below this max. value
	double fInputVtxDistY;  // use vtx given at input if dist. [cm] to track in 3D-Y is below this max. value
};

#endif

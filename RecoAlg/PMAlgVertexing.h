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

#include "RecoAlg/PMAlg/PmaTrkCandidate.h"
#include "RecoAlg/PMAlg/PmaVtxCandidate.h"

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
	virtual ~PMAlgVertexing(void); // delete last produced tracks (if not passed to output)

	void reconfigure(const fhicl::ParameterSet& pset);

	void reset(void) { cleanTracks(); }

	/// Copy input tracks, find 3D vertices, connect tracks, break them or flip if needed,
	/// reoptimize track structures. Result is returned as a collection of new tracks, that
	/// replaces content of trk_input (old tracks are deleted).
	/// Vertices can be accessed with getVertices function.
	size_t run(pma::trk_candidates& trk_input);

	/// Copy input tracks, use provided 3D vertices to connect tracks, break tracks or flip if
	/// needed, reoptimize track structures. Result is returned as a collection of new tracks,
	/// that replaces content of trk_input (old tracks are deleted).
	/// Input vertices that were actually associated to tracks are copied to the output
	/// collection (use getVertices function).
	size_t run(pma::trk_candidates& trk_input,
	           const std::vector< TVector3 >& vtx_input);

	std::vector< std::pair< TVector3, std::vector< size_t > > >
		getVertices(const pma::trk_candidates& tracks) const;

private:
	std::vector< pma::VtxCandidate > firstPassCandidates(void);
	std::vector< pma::VtxCandidate > secondPassCandidates(void);
	size_t makeVertices(std::vector< pma::VtxCandidate >& candidates);

	/// Get dQ/dx sequence to detect various features.
	std::vector< std::pair<double, double> > getdQdx(const pma::Track3D& trk) const;

	/// Get convolution value.
	double convolute(size_t idx, size_t len, double* adc, double const* shape) const;

	/// Check if colinear in 3D and dQ/dx with no significant step.
	bool isSingleParticle(pma::Track3D* trk1, pma::Track3D* trk2) const;

	/// Find elastic scattering vertices on tracks, merge back tracks that were split
	/// during vertex finding. 3D angle between two tracks and dQ/dx is checked.
	void mergeBrokenTracks(pma::trk_candidates& trk_input) const;

	/// Split track and add vertex and reoptimize when dQ/dx step detected.
	void splitMergedTracks(pma::trk_candidates& trk_input) const;

	pma::trk_candidates fOutTracks;
	pma::trk_candidates fShortTracks;
	pma::trk_candidates fEmTracks;
	void cleanTracks(void);

	void sortTracks(const pma::trk_candidates& trk_input);
	void collectTracks(pma::trk_candidates& result);

	// Parameters used in the algorithm

	double fMinTrackLength;   // min. length of tracks used to find vtx candidates (short tracks attached later)

	// just to remember:
	//double fInputVtxDist2D; // use vtx given at input if dist. [cm] to track in all 2D projections is below this max. value
	//double fInputVtxDistY;  // use vtx given at input if dist. [cm] to track in 3D-Y is below this max. value
};

#endif

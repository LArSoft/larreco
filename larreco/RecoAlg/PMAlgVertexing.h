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

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"
#include "larreco/RecoAlg/PMAlg/PmaVtxCandidate.h"

// ROOT & C++
#include <memory>

namespace pma
{
	class PMAlgVertexing;
}

class pma::PMAlgVertexing
{
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Atom<double> MinTrackLength {
			Name("MinTrackLength"),
			Comment("min. length of tracks used to find vtx candidates (short tracks attached later)")
		};

		fhicl::Atom<bool> FindKinks {
			Name("FindKinks"),
			Comment("detect significant kinks on long tracks")
		};

		fhicl::Atom<double> KinkMinDeg {
			Name("KinkMinDeg"),
			Comment("min. angle [deg] in XY of a kink")
		};

		fhicl::Atom<double> KinkMinStd {
			Name("KinkMinStd"),
			Comment("threshold in no. of stdev of all segment angles needed to tag a kink")
		};
    };

	PMAlgVertexing(const Config& config);
	void reconfigure(const Config& config);

	PMAlgVertexing(const fhicl::ParameterSet& pset) :
		PMAlgVertexing(fhicl::Table<Config>(pset, {})())
	{}

	virtual ~PMAlgVertexing(void); // delete last produced tracks (if not passed to output)

	void reset(void) { cleanTracks(); }

	/// Copy input tracks, find 3D vertices, connect tracks, break them or flip if needed,
	/// reoptimize track structures. Result is returned as a collection of new tracks, that
	/// replaces content of trk_input (old tracks are deleted).
	/// Vertices can be accessed with getVertices function.
	size_t run(pma::TrkCandidateColl & trk_input);

	/// Copy input tracks, use provided 3D vertices to connect tracks, break tracks or flip if
	/// needed, reoptimize track structures. Result is returned as a collection of new tracks,
	/// that replaces content of trk_input (old tracks are deleted).
	/// Input vertices that were actually associated to tracks are copied to the output
	/// collection (use getVertices function).
	size_t run(pma::TrkCandidateColl & trk_input,
	           const std::vector< TVector3 >& vtx_input);

	std::vector< std::pair< TVector3, std::vector< std::pair< size_t, bool > > > >
		getVertices(const pma::TrkCandidateColl& tracks, bool onlyBranching = false) const;

	std::vector< std::pair< TVector3, size_t > > getKinks(const pma::TrkCandidateColl& tracks) const;

private:
	bool has(const std::vector<size_t>& v, size_t idx) const
	{
		for (auto c : v) if (c == idx) return true;
		return false;
	}

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
	void mergeBrokenTracks(pma::TrkCandidateColl& trk_input) const;

	/// Split track and add vertex and reoptimize when dQ/dx step detected.
	void splitMergedTracks(pma::TrkCandidateColl& trk_input) const;

	/// Remove penalty on the angle if kink detected and reopt track.
	void findKinksOnTracks(pma::TrkCandidateColl& trk_input) const;

	pma::TrkCandidateColl fOutTracks;
	pma::TrkCandidateColl fShortTracks;
	pma::TrkCandidateColl fEmTracks;
	void cleanTracks(void);

	void sortTracks(const pma::TrkCandidateColl & trk_input);
	void collectTracks(pma::TrkCandidateColl & result);

	// Parameters used in the algorithm

	double fMinTrackLength;   // min. length of tracks used to find vtx candidates (short tracks attached later)

	bool fFindKinks;          // detect significant kinks on long tracks (need min. 5 nodes to collect angle stats)
    double fKinkMinDeg;       // min. angle [deg] in XY of a kink
	double fKinkMinStd;       // threshold in no. of stdev of all segment angles needed to tag a kink

	// just to remember:
	//double fInputVtxDist2D; // use vtx given at input if dist. [cm] to track in all 2D projections is below this max. value
	//double fInputVtxDistY;  // use vtx given at input if dist. [cm] to track in 3D-Y is below this max. value
};

#endif

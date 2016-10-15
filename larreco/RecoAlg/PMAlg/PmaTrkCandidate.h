/**
 *  @file   PmaTrkCandidate.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Track finding helper for the Projection Matching Algorithm
 *
 *          Candidate for 3D track. Used to test 2D cluster associations, validadion result, MSE value.
 *          See PmaTrack3D.h file for details.
 */

#ifndef TrkCandidate_h
#define TrkCandidate_h

#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"

namespace pma
{
	// these types to be replaced with use of feature proposed in redmine #12602
	typedef std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > view_hitmap;
	typedef std::map< unsigned int, view_hitmap > tpc_view_hitmap;
	typedef std::map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

	class TrkCandidate;
	class TrkCandidateColl;
}

class pma::TrkCandidate
{
public:
	TrkCandidate(void);
	TrkCandidate(pma::Track3D* trk, int key = -1, int tid = -1);

	bool IsValid(void) const { return fTrack; }

	bool IsGood(void) const { return fGood; }
	void SetGood(bool b) { fGood = b; }

	pma::Track3D* Track(void) const { return fTrack; }
	void SetTrack(pma::Track3D* trk);
	void DeleteTrack(void);

	const std::vector< size_t > & Clusters(void) const { return fClusters; }
	std::vector< size_t > & Clusters(void) { return fClusters; }

	/// Get key of an external object (like a source PFParticle) associated to this track candidate.
	int Key(void) const { return fKey; }

	/// Set key of an external object associated to this track candidate.
	void SetKey(int key) { fKey = key; }

	int TreeId(void) const { return fTreeId; }
	void SetTreeId(int id) { fTreeId = id; }

	double Mse(void) const { return fMse; }
	void SetMse(double m) { fMse = m; }

	double Validation(void) const { return fValidation; }
	void SetValidation(double v) { fValidation = v; }

	int Parent(void) const { return fParent; }
	void SetParent(int idx) { fParent = idx; }

	const std::vector< size_t > & Daughters(void) const { return fDaughters; }
	std::vector< size_t > & Daughters(void) { return fDaughters; }

private:
	int fParent;
	std::vector< size_t > fDaughters;

	pma::Track3D* fTrack;
	std::vector< size_t > fClusters;
	int fKey, fTreeId;

	double fMse, fValidation;

	bool fGood;
};

class pma::TrkCandidateColl
{
public:
	size_t size(void) const { return fCandidates.size(); }
	void resize(size_t n) { return fCandidates.resize(n); }
	bool empty(void) const { return fCandidates.empty(); }

	void push_back(const TrkCandidate & trk) { fCandidates.push_back(trk); }
	void erase_at(size_t pos) { fCandidates.erase(fCandidates.begin() + pos); }
	void clear(void) { fCandidates.clear(); }

	TrkCandidate & operator[] (size_t i) { return fCandidates[i]; }
	TrkCandidate const & operator[] (size_t i) const { return fCandidates[i]; }

	TrkCandidate & front(void) { return fCandidates.front(); }
	TrkCandidate const & front(void) const { return fCandidates.front(); }

	TrkCandidate & back(void) { return fCandidates.back(); }
	TrkCandidate const & back(void) const { return fCandidates.back(); }

	std::vector< TrkCandidate > const & tracks(void) const { return fCandidates; }
	std::vector< TrkCandidate > & tracks(void) { return fCandidates; }

    std::vector< TrkCandidate > const & parents(void) const { return fParents; }

	int getCandidateIndex(pma::Track3D const * candidate) const;
	void setParentDaughterConnections(void);

	void setTreeId(int id, size_t trkIdx, bool isRoot = true);
	int setTreeIds(void);

	void flipTreesToCoordinate(size_t coordinate);
	void flipTreesByDQdx();

	pma::Track3D* getTreeCopy(pma::TrkCandidateColl & dst, size_t trkIdx, bool isRoot = true);

private:
	std::vector< TrkCandidate > fCandidates;
	std::vector< TrkCandidate > fParents;
};

#endif


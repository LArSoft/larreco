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

#include "RecoAlg/PMAlg/PmaTrack3D.h"

namespace pma
{
	class TrkCandidate;
	typedef std::vector< TrkCandidate > trk_candidates;

	int getCandidateIndex(pma::trk_candidates const & tracks, pma::Track3D const * candidate);
	void setParentDaughterConnections(pma::trk_candidates& tracks);

	void setTreeId(pma::trk_candidates & tracks, int id, size_t trkIdx, bool isRoot = true);
	int setTreeIds(pma::trk_candidates & tracks);

	pma::Track3D* getTreeCopy(pma::trk_candidates & dst, const pma::trk_candidates & src, size_t trkIdx, bool isRoot = true);
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

#endif


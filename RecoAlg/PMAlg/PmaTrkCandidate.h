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
}

class pma::TrkCandidate
{
public:

	TrkCandidate(void);
	TrkCandidate(pma::Track3D* trk, int key = -1);

	bool IsValid(void) const { return fTrack; }

	bool IsGood(void) const { return fGood; }
	void SetGood(bool b) { fGood = b; }

	pma::Track3D* Track(void) const { return fTrack; }
	void SetTrack(pma::Track3D* trk);
	void DeleteTrack(void);

	const std::vector< size_t > & Clusters(void) const { return fClusters; }
	std::vector< size_t > & Clusters(void) { return fClusters; }

	int Key(void) const { return fKey; }
	void SetKey(int key) { fKey = key; }

	double Mse(void) const { return fMse; }
	void SetMse(double m) { fMse = m; }

	double Validation(void) const { return fValidation; }
	void SetValidation(double v) { fValidation = v; }

private:
	pma::Track3D* fTrack;
	std::vector< size_t > fClusters;
	int fKey;
	double fMse;
	double fValidation;
	bool fGood;
};

#endif


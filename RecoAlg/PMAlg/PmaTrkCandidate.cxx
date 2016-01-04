/**
 *  @file   PmaTrkCandidate.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Track finding helper for the Projection Matching Algorithm
 *
 *          Candidate for 3D track. Used to test 2D cluster associations, validadion result, MSE value.
 *          See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/PmaTrkCandidate.h"

pma::TrkCandidate::TrkCandidate(void) :
	fParent(-1),
	fTrack(0),
	fKey(-1),
	fMse(0), fValidation(0),
	fGood(false)
{
}
// ------------------------------------------------------

pma::TrkCandidate::TrkCandidate(pma::Track3D* trk, int key) :
	fParent(-1),
	fTrack(trk),
	fKey(key),
	fMse(0), fValidation(0),
	fGood(false)
{
}
// ------------------------------------------------------

void pma::TrkCandidate::SetTrack(pma::Track3D* trk)
{
	if (fTrack) delete fTrack;
	fTrack = trk;
}
// ------------------------------------------------------

void pma::TrkCandidate::DeleteTrack(void)
{
	if (fTrack) delete fTrack;
	fTrack = 0;
}
// ------------------------------------------------------

int pma::getCandidateIndex(pma::trk_candidates const & tracks, pma::Track3D const * candidate)
{
	for (size_t t = 0; t < tracks.size(); ++t)
		if (tracks[t].Track() == candidate) return t;
	return -1;
}

void pma::setParentDaughterConnections(pma::trk_candidates& tracks)
{
	for (size_t t = 0; t < tracks.size(); ++t)
	{
		pma::Track3D const * trk = tracks[t].Track();
		pma::Node3D const * firstNode = trk->Nodes().front();
		if (firstNode->Prev())
		{
			pma::Track3D const * parentTrk = static_cast< pma::Segment3D* >(firstNode->Prev())->Parent();
			tracks[t].SetParent(pma::getCandidateIndex(tracks, parentTrk));
		}
		for (auto node : trk->Nodes())
			for (size_t i = 0; i < node->NextCount(); ++i)
		{
			pma::Track3D const * daughterTrk = static_cast< pma::Segment3D* >(node->Next(i))->Parent();
			if (daughterTrk != trk)
			{
				int idx = pma::getCandidateIndex(tracks, daughterTrk);
				if (idx >= 0) tracks[t].Daughters().push_back((size_t)idx);
			}
		}
	}
}
// ------------------------------------------------------


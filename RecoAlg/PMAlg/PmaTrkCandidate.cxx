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
	fTrack(0),
	fMse(0), fValidation(0),
	fGood(false)
{

}

void pma::TrkCandidate::SetTrack(pma::Track3D* trk)
{
	if (fTrack) delete fTrack;
	fTrack = trk;
}

void pma::TrkCandidate::DeleteTrack(void)
{
	if (fTrack) delete fTrack;
	fTrack = 0;
}


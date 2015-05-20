/**
 *  @file   PmaTrack3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 */

#include "RecoAlg/PMAlg/PmaTrack3D.h"

pma::Track3D::Track3D(void)
{
}

bool pma::Track3D::HasRefPoint(TVector3* p) const
{
	for (size_t i = 0; i < fAssignedPoints.size(); i++)
		if (fAssignedPoints[i] == p) return true;
	return false;
}


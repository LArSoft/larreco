/**
 *  @file   PmaTrack3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 *          See PmaTrack3D.h file for details.
 */

#include "RecoAlg/PMAlg/PmaTrack3D.h"

pma::Track3D::Track3D(void)
{
}

pma::Track3D::~Track3D(void)
{
	for (size_t i = 0; i < fHits.size(); i++) delete fHits[i];
	for (size_t i = 0; i < fAssignedPoints.size(); i++) delete fAssignedPoints[i];
}

bool pma::Track3D::push_back(const recob::Hit& hit)
{
	pma::Hit3D* h3d = new pma::Hit3D(hit);
	for (size_t i = 0; i < fHits.size(); i++)
	{
		pma::Hit3D* hit_i = fHits[i];
		if ((h3d->PeakTime() == hit_i->PeakTime()) &&
		    (h3d->Wire() == hit_i->Wire()) &&
		    (h3d->View2D() == hit_i->View2D()) &&
		    (h3d->TPC() == hit_i->TPC())) return false;
	}
	fHits.push_back(h3d);
	return true;
}

bool pma::Track3D::HasRefPoint(TVector3* p) const
{
	for (size_t i = 0; i < fAssignedPoints.size(); i++)
		if (fAssignedPoints[i] == p) return true;
	return false;
}


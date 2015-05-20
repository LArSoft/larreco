/**
 *  @file   PmaTrack3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 *          Based on the "Precise 3D track reco..." AHEP (2013) 260820, with all tricks that we
 *          developed later and with the work for full-event topology optimization that is still
 *          under construction.
 *
 *          Current status: basic functionality for single segment optimization.
 */

#ifndef PmaTrack3D_h
#define PmaTrack3D_h

#include "RecoAlg/PMAlg/PmaHit3D.h"
#include "RecoAlg/PMAlg/PmaNode3D.h"
#include "RecoAlg/PMAlg/PmaSegment3D.h"

namespace pma
{
	class Track3D;
}

class pma::Track3D
{
public:
	Track3D(void);

	pma::Hit3D* & operator [] (size_t index) { return fHits[index]; }
	pma::Hit3D* const & operator [] (size_t index) const { return fHits[index]; }

	size_t size() const { return fHits.size(); }

	bool HasRefPoint(TVector3* p) const;

private:
	std::vector< pma::Hit3D* > fHits;
	std::vector< TVector3* > fAssignedPoints;
};

#endif

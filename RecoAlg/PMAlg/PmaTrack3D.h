/**
 *  @file   PmaTrack3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 */

#ifndef PmaTrack3D_h
#define PmaTrack3D_h

#include "RecoAlg/PMAlg/PmaHit3D.h"

namespace pma
{
	class Track3D;
}

class pma::Track3D
{
public:
	Track3D(void);

private:
	std::vector< pma::Hit3D > fHits;
};

#endif

/**
 *  @file   PmaHit3D.h
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 */

#ifndef PmaHit3D_h
#define PmaHit3D_h


#include "RecoBase/Hit.h"

namespace pma
{
	class Hit3D;

	class Track3D;
}

class pma::Hit3D : public recob::Hit
{
	friend class Track3D;

public:
	Hit3D(void);

private:

}

#endif

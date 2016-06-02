/**
 *  @file   TssHit2D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Hit pos in cm and original recob hit ptr.
 */

#ifndef TssHit2D_h
#define TssHit2D_h

#include "lardataobj/RecoBase/Hit.h"

#include "larcore/Geometry/Geometry.h"

#include <functional>

#include "TVector2.h"

namespace tss
{
	class Hit2D;
	//struct bTrajectory3DOrderLess;
}

class tss::Hit2D
{
	//friend struct bTrajectory3DOrderLess;

public:
	Hit2D(void);
	Hit2D(const art::Ptr< recob::Hit > & src);
	Hit2D(const tss::Hit2D & src);

	art::Ptr< recob::Hit > Hit2DPtr(void) const { return fHit; }

	TVector2 const & Point2D(void) const { return fPoint2D; }

	unsigned int Cryo(void) const { return fHit->WireID().Cryostat; }
	unsigned int TPC(void) const { return fHit->WireID().TPC; }
	unsigned int View(void) const { return fPlane; }
	unsigned int Wire(void) const { return fWire; }
	float PeakTime(void) const { return fHit->PeakTime(); }
	int StartTick(void) const { return fHit->StartTick(); }
	int EndTick(void) const { return fHit->EndTick(); }

	float SummedADC(void) const { return fHit->SummedADC(); }
	float GetAmplitude(void) const { return fHit->PeakAmplitude(); }

private:

	art::Ptr< recob::Hit > fHit;  // source 2D hit

	unsigned int fPlane, fWire;

	TVector2 fPoint2D;       // hit position in 2D wire view, scaled to [cm]

};

#endif


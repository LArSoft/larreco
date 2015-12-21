/**
 *  @file   TssHit2D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Hit pos in cm and original recob hit ptr.
 */

#include "TssHit2D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "lardata/DetectorInfo/DetectorProperties.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"


tss::Hit2D::Hit2D(void) :
	fPlane(0), fWire(0),
	fPoint2D(0, 0)
{
}

tss::Hit2D::Hit2D(const art::Ptr< recob::Hit > & src) :
	fHit(src)
{
	fPlane = src->WireID().Plane;
	fWire = src->WireID().Wire;

	fPoint2D = pma::WireDriftToCm(
		fWire, src->PeakTime(), fPlane,
		src->WireID().TPC, src->WireID().Cryostat);
}

tss::Hit2D::Hit2D(const tss::Hit2D & src) :
	fHit(src.fHit),
	fPlane(src.fPlane), fWire(src.fWire),
	fPoint2D(src.fPoint2D)
{
}


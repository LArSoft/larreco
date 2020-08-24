/**
 *  @file   TssHit2D.cxx
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Hit pos in cm and original recob hit ptr.
 */

#include "TssHit2D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

tss::Hit2D::Hit2D(detinfo::DetectorPropertiesData const& detProp, const art::Ptr<recob::Hit>& src)
  : fHit(src)
{
  fPlane = src->WireID().Plane;
  fWire = src->WireID().Wire;

  fPoint2D = pma::WireDriftToCm(
    detProp, fWire, src->PeakTime(), fPlane, src->WireID().TPC, src->WireID().Cryostat);
}

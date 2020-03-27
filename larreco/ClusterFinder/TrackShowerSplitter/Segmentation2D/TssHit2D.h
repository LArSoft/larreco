/**
 *  @file   TssHit2D.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Hit pos in cm and original recob hit ptr.
 */

#ifndef TssHit2D_h
#define TssHit2D_h

#include "canvas/Persistency/Common/Ptr.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
namespace detinfo {
  class DetectorPropertiesData;
}

#include "TVector2.h"

namespace tss {
  class Hit2D;
}

class tss::Hit2D {
public:
  Hit2D(detinfo::DetectorPropertiesData const& detProp, const art::Ptr<recob::Hit>& src);

  art::Ptr<recob::Hit>
  Hit2DPtr() const
  {
    return fHit;
  }

  TVector2 const&
  Point2D() const
  {
    return fPoint2D;
  }

  unsigned int
  Cryo() const
  {
    return fHit->WireID().Cryostat;
  }
  unsigned int
  TPC() const
  {
    return fHit->WireID().TPC;
  }
  unsigned int
  View() const
  {
    return fPlane;
  }
  unsigned int
  Wire() const
  {
    return fWire;
  }
  float
  PeakTime() const
  {
    return fHit->PeakTime();
  }
  int
  StartTick() const
  {
    return fHit->StartTick();
  }
  int
  EndTick() const
  {
    return fHit->EndTick();
  }

  float
  SummedADC() const
  {
    return fHit->SummedADC();
  }
  float
  GetAmplitude() const
  {
    return fHit->PeakAmplitude();
  }

private:
  art::Ptr<recob::Hit> fHit; // source 2D hit

  unsigned int fPlane, fWire;

  TVector2 fPoint2D; // hit position in 2D wire view, scaled to [cm]
};

#endif

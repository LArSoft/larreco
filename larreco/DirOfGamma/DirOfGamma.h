////////////////////////////////////////////////////////////////////////
// Class:
// Module Type:
// File: DirOfGamma.h
// Authors: dorota.stefan@cern.ch robert.sulej@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef DirOfGamma_h
#define DirOfGamma_h

#include "lardataobj/RecoBase/Hit.h"
namespace detinfo {
  class DetectorPropertiesData;
}

#include "canvas/Persistency/Common/Ptr.h"

#include "TVector2.h"

#include <vector>

namespace ems {
  class Hit2D;
  class Bin2D;
  class EndPoint;
  class DirOfGamma;
  class bDistCentMore2D;
  class bDistCentLess2D;
}

class ems::Hit2D {
public:
  Hit2D(detinfo::DetectorPropertiesData const& detProp, art::Ptr<recob::Hit> src);

  TVector2 const&
  GetPointCm() const
  {
    return fPoint;
  }
  double
  GetCharge() const
  {
    return fCharge;
  }

  art::Ptr<recob::Hit> const&
  GetHitPtr() const
  {
    return fHit;
  }

private:
  double fCharge;

  TVector2 fPoint;

  art::Ptr<recob::Hit> fHit;
};

class ems::Bin2D {
public:
  Bin2D(const TVector2& center);

  void Add(Hit2D* hit);

  void Sort();

  void SortLess();

  double
  GetTotCharge() const
  {
    return fTotCharge;
  }

  unsigned int
  Size() const
  {
    return fSize;
  }

  std::vector<Hit2D*> const&
  GetHits2D() const
  {
    return fHits2D;
  }

  const TVector2&
  GetCenter() const
  {
    return fCenter2D;
  }

  std::vector<art::Ptr<recob::Hit>> GetIniHits(const double radius = 10.0,
                                               const unsigned int nhits = 10) const;

private:
  const TVector2& fCenter2D;
  std::vector<Hit2D*> fHits2D;
  double fTotCharge;
  unsigned int fSize;
};

class ems::EndPoint {
public:
  EndPoint(const Hit2D& center, const std::vector<Hit2D*>& hits, unsigned int nbins);

  TVector2 const&
  GetPosition() const
  {
    return fCenter2D.GetPointCm();
  }

  double GetAsymmetry() const;

  double
  GetMaxCharge() const
  {
    return fMaxCharge;
  }

  Bin2D const&
  MaxChargeBin() const
  {
    return fBins[fMaxChargeIdBin];
  }

  std::vector<Bin2D> const&
  GetBins2D() const
  {
    return fBins;
  }

  art::Ptr<recob::Hit> const&
  GetHit() const
  {
    return fCenter2D.GetHitPtr();
  }

  const std::vector<art::Ptr<recob::Hit>>
  GetIniHits() const
  {
    return MaxChargeBin().GetIniHits();
  }

  size_t const&
  GetPlane() const
  {
    return fPlane;
  }
  size_t const&
  GetTPC() const
  {
    return fTpc;
  }
  size_t const&
  GetCryo() const
  {
    return fCryo;
  }

private:
  Hit2D fCenter2D;
  std::vector<Hit2D*> fPoints2D;
  size_t fNbins;

  double fMaxCharge;
  double fMeanCharge;

  std::vector<Bin2D> fBins;

  size_t fMaxChargeIdBin;

  void FillBins();
  void ComputeMaxCharge();
  void ComputeMeanCharge();

  size_t fPlane;
  size_t fTpc;
  size_t fCryo;
};

class ems::DirOfGamma {
public:
  DirOfGamma(const detinfo::DetectorPropertiesData& detProp,
             const std::vector<art::Ptr<recob::Hit>>& src,
             unsigned int nbins,
             unsigned int idcl);
  ~DirOfGamma()
  {
    for (unsigned int i = 0; i < fPoints2D.size(); ++i)
      delete fPoints2D[i];
  }

  TVector2 const&
  GetBaryCenterCm() const
  {
    return fBaryCenter;
  }

  std::vector<Hit2D*> const&
  GetHits2D() const
  {
    return fPoints2D;
  }

  std::vector<EndPoint> const&
  GetCandidates() const
  {
    return fCandidates;
  }

  void
  SetIdCandidate(size_t id)
  {

    fIsCandidateIDset = true;
    fCandidateID = id;
    fStartHit = fCandidates[id].GetHit();
    fStartPoint = fCandidates[id].GetPosition();
    fIniHits = fCandidates[id].MaxChargeBin().GetIniHits();
  }

  const size_t
  GetIdCandidate()
  {
    return fCandidateID;
  }

  art::Ptr<recob::Hit> const&
  GetFirstHit() const
  {
    return fStartHit;
  }

  std::vector<art::Ptr<recob::Hit>> const&
  GetHits()
  {
    return fHits;
  }

  TVector2 const&
  GetFirstPoint() const
  {
    return fStartPoint;
  }

  std::vector<art::Ptr<recob::Hit>> const&
  GetIniHits() const
  {
    return fIniHits;
  }

  size_t const
  GetIdCl() const
  {
    return fIdCl;
  }

private:
  size_t fNbins;
  size_t fIdCl;
  size_t fCandidateID;

  std::vector<Hit2D*> fPoints2D;
  std::vector<Bin2D> fBins;
  std::vector<EndPoint> fCandidates;

  art::Ptr<recob::Hit> fStartHit;
  TVector2 fStartPoint;
  std::vector<art::Ptr<recob::Hit>> fIniHits;
  std::vector<art::Ptr<recob::Hit>> fHits;

  void FindInitialPart();

  void FillBins();

  bool FindCandidates();
  void ComputeBaryCenter();
  void ComputeMaxDist();
  void ComputeMaxCharge();
  void ComputeFinalValues();

  TVector2 fBaryCenter;

  float fNormDist;
  float fNormCharge;

  bool fIsCandidateIDset;
};

class ems::bDistCentMore2D : public std::binary_function<Hit2D*, Hit2D*, bool> {
public:
  bDistCentMore2D(const TVector2& c) : center(c) {}

  bool
  operator()(Hit2D* p1, Hit2D* p2)
  {
    double dx = p1->GetPointCm().X() - center.X();
    double dy = p1->GetPointCm().Y() - center.Y();
    double b1 = dx * dx + dy * dy;
    dx = p2->GetPointCm().X() - center.X();
    dy = p2->GetPointCm().Y() - center.Y();
    double b2 = dx * dx + dy * dy;

    return b1 > b2;
  }

private:
  TVector2 center;
};

class ems::bDistCentLess2D : public std::binary_function<Hit2D*, Hit2D*, bool> {
public:
  bDistCentLess2D(const TVector2& c) : center(c) {}

  bool
  operator()(Hit2D* p1, Hit2D* p2)
  {
    double dx = p1->GetPointCm().X() - center.X();
    double dy = p1->GetPointCm().Y() - center.Y();
    double b1 = dx * dx + dy * dy;
    dx = p2->GetPointCm().X() - center.X();
    dy = p2->GetPointCm().Y() - center.Y();
    double b2 = dx * dx + dy * dy;

    return b1 < b2;
  }

private:
  TVector2 center;
};

#endif

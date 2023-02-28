/**
 *  @file   PmaVtxCandidate.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Vertex finding helper for the Projection Matching Algorithm
 *
 *          Candidate for 3D vertex. Used to test intersections and join tracks in vertices.
 *          See PmaTrack3D.h file for details.
 */

#ifndef VtxCandidate_h
#define VtxCandidate_h

#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"

namespace pma {
  class VtxCandidate;
}

class pma::VtxCandidate {
public:
  static constexpr double kMaxDistToTrack{4.0}; // max. dist. track to center to create vtx
  static constexpr double kMinDistToNode{2.0};  // min. dist. to node needed to split segment

  VtxCandidate(double segMinLength = 0.5)
    : tracksJoined(false)
    , fSegMinLength(segMinLength)
    , fMse(0.0)
    , fMse2D(0.0)
    , fCenter(0., 0., 0.)
    , fErr(0., 0., 0.)
  {}

  bool Has(pma::Track3D* trk) const;

  bool Has(const VtxCandidate& other) const;

  bool IsAttached(pma::Track3D* trk) const;

  bool IsAttached(const VtxCandidate& other) const;

  bool HasLoops() const;

  bool Add(const pma::TrkCandidate& trk);

  double ComputeMse2D();

  double Test(const VtxCandidate& other) const;

  double MaxAngle(double minLength = 0.0) const;

  size_t Size() const { return fAssigned.size(); }
  size_t Size(double minLength) const;

  bool MergeWith(const VtxCandidate& other);

  double Compute();

  bool JoinTracks(detinfo::DetectorPropertiesData const& detProp,
                  pma::TrkCandidateColl& tracks,
                  pma::TrkCandidateColl& src);

  const TVector3& Center() const { return fCenter; }
  double Mse() const { return fMse; }
  double Mse2D() const { return fMse2D; }

  std::pair<pma::Track3D*, size_t> Track(size_t i) const
  {
    return std::pair<pma::Track3D*, size_t>(fAssigned[i].first.Track(), fAssigned[i].second);
  }

private:
  bool has(const std::vector<int>& v, int id) const
  {
    for (auto c : v)
      if (c == id) return true;
    return false;
  }

  bool tracksJoined;
  double fSegMinLength, fMse, fMse2D;
  std::vector<std::pair<pma::TrkCandidate, size_t>> fAssigned;
  TVector3 fCenter, fErr;
};

#endif

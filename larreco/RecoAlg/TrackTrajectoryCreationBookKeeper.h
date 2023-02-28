#ifndef TRACKTRAJECTORYCREATIONBOOKKEEPER_H
#define TRACKTRAJECTORYCREATIONBOOKKEEPER_H

#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

namespace trkmkr {

  /**
   * @file  larreco/RecoAlg/TrackTrajectoryCreationBookKeeper.h
   * @class trkmkr::TrackTrajectoryCreationBookKeeper
   *
   * @brief Helper class to aid the creation of a recob::TrackTrajectory, keeping data vectors in sync.
   *
   * Helper class to aid the creation of a recob::TrackTrajectory, keeping data vectors (Point_t, Vector_t, PointFlags_t, Hit) in sync.
   * Elements of those vectors are added sequentially using the addPoint functions.
   * Once all points have been added a call to the function finalizeTrackTrajectory, builds the track moving the content of the vectors.
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  using namespace recob;
  using Point_t = tracking::Point_t;
  using Vector_t = tracking::Vector_t;
  using PointFlags_t = TrackTrajectory::PointFlags_t;

  class TrackTrajectoryCreationBookKeeper {
  public:
    /// Constructor: needs reference to output hit vector, and hasMomenta bool (true if Vector_t are momenta, false if they are directions).
    TrackTrajectoryCreationBookKeeper(std::vector<art::Ptr<Hit>>& outhits, bool hasMomenta)
      : hasMomenta_(hasMomenta), hits(&outhits)
    {
      hits->clear();
    }
    //
    //@{
    /// Avoid copies of this object
    TrackTrajectoryCreationBookKeeper(const TrackTrajectoryCreationBookKeeper&) = delete;
    TrackTrajectoryCreationBookKeeper(TrackTrajectoryCreationBookKeeper&&) = delete;
    TrackTrajectoryCreationBookKeeper& operator=(const TrackTrajectoryCreationBookKeeper&) = delete;
    TrackTrajectoryCreationBookKeeper& operator=(TrackTrajectoryCreationBookKeeper&&) = delete;
    //@}
    //
    //@{
    /// Add a single point; different version of the functions are provided using const references or rvalue references.
    void addPoint(const Point_t& point,
                  const Vector_t& vect,
                  art::Ptr<Hit> hit,
                  const PointFlags_t& flag)
    {
      positions.push_back(point);
      momenta.push_back(vect);
      hits->push_back(hit);
      flags.push_back(flag);
    }
    void addPoint(Point_t&& point, Vector_t&& vect, art::Ptr<Hit> hit, PointFlags_t&& flag)
    {
      positions.push_back(std::move(point));
      momenta.push_back(std::move(vect));
      hits->push_back(hit);
      flags.push_back(std::move(flag));
    }
    //@}
    //
    /// Get the finalized recob::TrackTrajectory object; internal data vectors are moved so no more points should be added.
    TrackTrajectory finalizeTrackTrajectory()
    {
      return TrackTrajectory(
        std::move(positions), std::move(momenta), std::move(flags), hasMomenta_);
    }
    //
  private:
    bool hasMomenta_;
    std::vector<art::Ptr<Hit>>* hits;
    std::vector<Point_t> positions;
    std::vector<Vector_t> momenta;
    std::vector<PointFlags_t> flags;
    //
  };

}
#endif

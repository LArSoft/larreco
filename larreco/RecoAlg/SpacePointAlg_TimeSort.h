#ifndef SPACEPOINTALG_TIMESORT_H
#define SPACEPOINTALG_TIMESORT_H
/*!
 * Title:   SpacePointAlg_TimeSort class
 * Author:  wketchum@lanl.gov
 * Inputs:  std::vector<recob::Hit> (one for each plane)
 * Outputs: std::vector<recob::SpacePoint>
 *
 * Description:
 * This space point algorithm tries to improve speed by
 * (1) keeping hit collections distinct among planes;
 * (2) sorting hit collections by time; and,
 * (3) having a lookup table for (y,z) coordinate positions.
 * The original use case for this code was with the TTHitFinder,
 * which produces an incredibly large number of hits per plane,
 * making some sorted space point alg more attractive.
 *
 * This code is totally microboone specific, btw.
 */

// LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
namespace detinfo {
  class DetectorPropertiesData;
}

// Frmaework includes
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/fwd.h"

//boost includes
#include "boost/multi_array.hpp"

#include <memory>
#include <vector>

namespace sppt {

  class SpacePointAlg_TimeSort {
  public:
    explicit SpacePointAlg_TimeSort(fhicl::ParameterSet const& pset);

    void setTimeOffsets(detinfo::DetectorPropertiesData const& detProp);
    void fillCoordinatesArrays();

    void createSpacePoints(
      detinfo::DetectorPropertiesData const& detProp,
      std::vector<art::Ptr<recob::Hit>>& hitVec_U,
      std::vector<art::Ptr<recob::Hit>>& hitVec_V,
      std::vector<art::Ptr<recob::Hit>>& hitVec_Y,
      std::unique_ptr<std::vector<recob::SpacePoint>>& spptCollection,
      std::unique_ptr<std::vector<std::vector<art::Ptr<recob::Hit>>>>& spptAssociatedHits);

  private:
    float fTimeDiffMax; /// Maximum allowed time difference
    float fYDiffMax;    /// Maximum allowed y-coordinate difference
    float fZDiffMax;    /// Maximum allowed z-coordinate difference

    bool TIME_OFFSET_SET{false};
    bool COORDINATES_FILLED{false};

    double TIME_OFFSET_U;
    double TIME_OFFSET_V;
    double TIME_OFFSET_Y;
    double TICKS_TO_X;

    boost::multi_array<double, 2> coordinates_UV_y;
    boost::multi_array<double, 2> coordinates_UV_z;
    boost::multi_array<double, 2> coordinates_UY_y;
    boost::multi_array<double, 2> coordinates_UY_z;

    void sortHitsByTime(std::vector<art::Ptr<recob::Hit>>& hits_handle) const;

  }; //class SpacePointAlg_TimeSort

} //end sppt namespace

#endif

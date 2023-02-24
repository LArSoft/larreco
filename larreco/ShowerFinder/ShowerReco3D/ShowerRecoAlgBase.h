/**
 * \file ShowerRecoAlgBase.h
 *
 * \ingroup ShowerReco3D
 *
 * \brief Class def header for a class ShowerRecoAlgBase
 *
 * @author kazuhiro
 */

/** \addtogroup ShowerReco3D

    @{*/
#ifndef RECOTOOL_SHOWERRECOALGBASE_H
#define RECOTOOL_SHOWERRECOALGBASE_H

#include <vector>

#include "larcorealg/Geometry/fwd.h"
#include "lardata/Utilities/PxUtils.h"
#include "lardataobj/RecoBase/Shower.h"
namespace calo {
  class CalorimetryAlg;
}
namespace cluster {
  class ClusterParamsAlg;
}
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace showerreco {

  struct ShowerCluster_t {
    util::PxPoint start_point;
    util::PxPoint end_point;
    double angle_2d;
    unsigned short plane_id;
    std::vector<::util::PxHit> hit_vector;
  };

  /**
     \class ShowerRecoAlgBase
     User defined class ShowerRecoAlgBase ... these comments are used to generate
     doxygen documentation!
  */
  class ShowerRecoAlgBase {
  public:
    virtual ~ShowerRecoAlgBase() = default;

    /// Function to reset algorithm, to be called @ beginning of each event
    virtual void Reset();

    /// Setter for a matched combination of clusters
    virtual void AppendInputClusters(const std::vector<cluster::ClusterParamsAlg>& cpan_v);

    /// Execute reconstruction
    std::vector<recob::Shower> Reconstruct(geo::GeometryCore const& geom,
                                           geo::WireReadoutGeom const& wireReadoutGeom,
                                           detinfo::DetectorClocksData const& clockData,
                                           detinfo::DetectorPropertiesData const& detProp);

  protected:
    /// Function to reorganize input cluster information
    virtual void ProcessInputClusters() {}

    /// Function to reconstruct one shower
    virtual ::recob::Shower RecoOneShower(
      geo::GeometryCore const& geom,
      geo::WireReadoutGeom const& wireReadoutGeom,
      detinfo::DetectorClocksData const& clockData,
      detinfo::DetectorPropertiesData const& detProp,
      const std::vector<showerreco::ShowerCluster_t>& clusters) = 0;

  protected:
    /// Input clusters
    std::vector<std::vector<showerreco::ShowerCluster_t>> fInputClusters;
  };
}

#endif
/** @} */ // end of doxygen group

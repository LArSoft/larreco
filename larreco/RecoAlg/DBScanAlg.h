/////////////////////////////////////////////////////////////////
//  \fileDBScanAlg.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
#ifndef DBSCANALG_H
#define DBSCANALG_H

#include <set>
#include <stdint.h>
#include <vector>

#include "RStarTree/RStarTree.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // for WireID
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace fhicl {
  class ParameterSet;
}

namespace recob {
  class Hit;
}

// RStarTree related infrastructure
//
// Our core objects have a physical extent (i.e. there are not
// points), but a R*-tree should be able to deal with that.
using RTree = RStarTree<uint32_t, 2, 32, 64>; // payload is just an index
using BoundingBox = RTree::BoundingBox;

struct dbsPoint {
  double x, y;
  double dx, dy;
  dbsPoint(double X = 0.0, double Y = 0.0, double dX = 0.0, double dY = 0.0)
    : x(X), y(Y), dx(dX), dy(dY){};
  BoundingBox bounds() const;
  void Expand(double DX, double DY)
  {
    dx += DX;
    dy += DY;
  };
};

namespace cluster {

  //---------------------------------------------------------------
  class DBScanAlg {
  public:
    explicit DBScanAlg(fhicl::ParameterSet const& pset);

    void InitScan(
      const detinfo::DetectorClocksData& clockData,
      const detinfo::DetectorPropertiesData& detProp,
      const std::vector<art::Ptr<recob::Hit>>& allhits,
      std::set<uint32_t> badChannels,
      const std::vector<geo::WireID>& wireids = std::vector<geo::WireID>()); //wireids is optional
    double getSimilarity(const std::vector<double> v1, const std::vector<double> v2);
    std::vector<unsigned int> findNeighbors(unsigned int pid, double threshold, double threshold2);
    void computeSimilarity();
    void run_cluster();
    double getSimilarity2(const std::vector<double> v1, const std::vector<double> v2);
    void computeSimilarity2();
    double getWidthFactor(const std::vector<double> v1, const std::vector<double> v2);
    void computeWidthFactor();

    std::vector<std::vector<unsigned int>> fclusters; ///< collection of something
    std::vector<std::vector<double>> fps;            ///< the collection of points we are working on
    std::vector<unsigned int> fpointId_to_clusterId; ///< mapping point_id -> clusterId
    std::vector<std::vector<double>> fsim;           ///<
    std::vector<std::vector<double>> fsim2;          ///<
    std::vector<std::vector<double>> fsim3;          ///<
    double fMaxWidth;

    RTree fRTree;
    std::vector<dbsPoint> fRect;

  private:
    // eps radius
    // Two points are neighbors if the distance
    // between them does not exceed threshold value.
    double fEps;
    double fEps2;
    //minimum number of points
    unsigned int fMinPts;
    // Which clustering to run
    unsigned int fClusterMethod;  ///< Which clustering method to use
    unsigned int fDistanceMetric; ///< Which distance metric to use

    // noise vector
    std::vector<bool> fnoise;
    std::vector<bool> fvisited;
    std::vector<double> fWirePitch;    ///< the pitch of the wires in each plane
    std::set<uint32_t> fBadChannels;   ///< set of bad channels in this detector
    std::vector<uint32_t> fBadWireSum; ///< running total of bad channels. Used for fast intervening
                                       ///< dead wire counting ala
                                       ///< fBadChannelSum[m]-fBadChannelSum[n].

    // Three differnt version of the clustering code
    void run_dbscan_cluster();
    void run_FN_cluster();
    void run_FN_naive_cluster();

    // Helper routined for run_dbscan_cluster() names and
    // responsibilities taken directly from the paper
    bool ExpandCluster(unsigned int point /* to be added */,
                       unsigned int clusterID /* which is being expanded */);
    std::set<unsigned int> RegionQuery(unsigned int point);
    // Helper for the accelerated run_FN_cluster()
    std::vector<unsigned int> RegionQuery_vector(unsigned int point);

  }; // class DBScanAlg
} // namespace

#endif // ifndef DBSCANALG_H

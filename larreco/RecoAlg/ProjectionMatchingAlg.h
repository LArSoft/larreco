////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Projection Matching Algorithm
// see RecoAlg/PMAlg/PmaTrack3D.h for more details.
//
//      Build 3D segments and whole tracks by matching an object 2D projections to hits, simultaneously
//      in multiple wire planes. Based on the algorithm first presented in "Precise 3D track reco..."
//      AHEP (2013) 260820, with all the tricks that we developed later and with the work for the full-event
//      topology optimization that is still under construction now (2015).
//
//      The algorithm class provides functionality to build a track from selected hits. These
//      can be detailed tracks or just simple segments (if the number of nodes to add is set to 0).
//      The parameters of optimization algorithm, fixed nodes and 3D reference points can be configured here.
//      Please, check the track making module to find a way of selecting appropriate clusteres:
//        PMAlgTrackMaker_module.cc
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ProjectionMatchingAlg_h
#define ProjectionMatchingAlg_h

// Framework includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"
namespace fhicl {
  class ParameterSet;
}

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
namespace detinfo {
  class DetectorPropertiesData;
}
namespace geo {
  class GeometryCore;
  class TPCGeo;
}
namespace img {
  class DataProviderAlg;
} 
namespace lariov {
  class ChannelStatusProvider;
}

// ROOT & C++
#include <map>
#include <memory>
#include <vector>
class TH1F;

namespace pma {
  class ProjectionMatchingAlg;
}

class pma::ProjectionMatchingAlg {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<double> OptimizationEps{
      Name("OptimizationEps"),
      Comment("relative change of the obj.fn which stops optimization after adding a node")};

    fhicl::Atom<double> FineTuningEps{
      Name("FineTuningEps"),
      Comment("relative change of the obj.fn which stops fine-tuning of optimized track")};

    fhicl::Atom<double> TrkValidationDist2D{
      Name("TrkValidationDist2D"),
      Comment("max. distance [cm] used in the track validation in the third plane")};

    fhicl::Atom<double> HitTestingDist2D{
      Name("HitTestingDist2D"),
      Comment("max. distance [cm] used in testing compatibility of hits with the track")};

    fhicl::Atom<double> MinTwoViewFraction{
      Name("MinTwoViewFraction"),
      Comment("min. fraction of track length covered with hits from many 2D views intertwinted "
              "with each other")};

    fhicl::Atom<double> NodeMargin3D{
      Name("NodeMargin3D"),
      Comment("margin in [cm] around TPC for allowed track node positions")};

    fhicl::Atom<double> HitWeightU{Name("HitWeightU"), Comment("weights used for hits in U plane")};

    fhicl::Atom<double> HitWeightV{Name("HitWeightV"), Comment("weights used for hits in V plane")};

    fhicl::Atom<double> HitWeightZ{Name("HitWeightZ"), Comment("weights used for hits in Z plane")};
  };

  ProjectionMatchingAlg(const Config& config);

  ProjectionMatchingAlg(const fhicl::ParameterSet& pset)
    : ProjectionMatchingAlg(fhicl::Table<Config>(pset, {})())
  {}

  /// Calculate the fraction of the track that is close to non-empty pixel
  /// (above thr value) in the ADC image of the testView (a view that was not
  /// used to build the track).
  double validate_on_adc(const detinfo::DetectorPropertiesData& detProp,
                         const lariov::ChannelStatusProvider& channelStatus,
                         const pma::Track3D& trk,
                         const img::DataProviderAlg& adcImage,
                         float thr) const;

  /// Calculate the fraction of the track that is closer than
  /// fTrkValidationDist2D to any hit from hits in the testView (a view that was
  /// not used to build the track). Creates also histograms of values in pixels
  /// for the passing and rejected points on the track, so the threshold value
  /// for the ADC-based calibration can be estimated.
  double validate_on_adc_test(const detinfo::DetectorPropertiesData& detProp,
                              const lariov::ChannelStatusProvider& channelStatus,
                              const pma::Track3D& trk,
                              const img::DataProviderAlg& adcImage,
                              const std::vector<art::Ptr<recob::Hit>>& hits,
                              TH1F* histoPassing,
                              TH1F* histoRejected) const;

  /// Calculate the fraction of the track that is closer than
  /// fTrkValidationDist2D to any hit from hits in their plane (a plane that was
  /// not used to build the track). Hits should be preselected, so all belong to
  /// the same plane.
  double validate(const detinfo::DetectorPropertiesData& detProp,
                  const lariov::ChannelStatusProvider& channelStatus,
                  const pma::Track3D& trk,
                  const std::vector<art::Ptr<recob::Hit>>& hits) const;

  /// Calculate the fraction of the 3D segment that is closer than
  /// fTrkValidationDist2D to any hit from hits in the testPlane of TPC/Cryo.
  /// Hits from the testPlane are preselected by this function among all
  /// provided (so a bit slower than fn above).
  double validate(const detinfo::DetectorPropertiesData& detProp,
                  const lariov::ChannelStatusProvider& channelStatus,
                  const TVector3& p0,
                  const TVector3& p1,
                  const std::vector<art::Ptr<recob::Hit>>& hits,
                  unsigned int testView,
                  unsigned int tpc,
                  unsigned int cryo) const;

  /// Calculate the fraction of trajectory seen by two 2D projections at least; even a
  /// prfect track starts/stops with the hit from one 2D view, then hits from other views
  /// come, which results with the fraction value high, but always < 1.0; wrong cluster
  /// matchings or incomplete tracks give significantly lower values.
  double twoViewFraction(pma::Track3D& trk) const;

  /// Count the number of hits that are closer than eps * fHitTestingDist2D to the track 2D projection.
  unsigned int
  testHits(detinfo::DetectorPropertiesData const& detProp,
           const pma::Track3D& trk,
           const std::vector<art::Ptr<recob::Hit>>& hits,
           double eps = 1.0) const
  {
    return trk.TestHits(detProp, hits, eps * fHitTestingDist2D);
  }

  /// Test if hits at the track endpoinds do not stick out of TPC which they belong to.
  /// Here one can implement some configurable margin if needed for real data imeprfections.
  bool
  isContained(const pma::Track3D& trk, float margin = 0.0F) const
  {
    return (trk.FirstElement()->SameTPC(trk.front()->Point3D(), margin) &&
            trk.LastElement()->SameTPC(trk.back()->Point3D(), margin));
  }

  /// Build a track from two sets of hits from single TPC, hits should origin
  /// from at least two
  /// wire planes; number of segments used to create the track depends on the
  /// number of hits.
  pma::Track3D* buildTrack(const detinfo::DetectorPropertiesData& detProp,
                           const std::vector<art::Ptr<recob::Hit>>& hits_1,
                           const std::vector<art::Ptr<recob::Hit>>& hits_2 = {}) const;

  /// Build a track from sets of hits, multiple TPCs are OK (like taken from
  /// PFParticles),
  /// as far as hits origin from at least two wire planes.
  pma::Track3D* buildMultiTPCTrack(const detinfo::DetectorPropertiesData& clockData,
                                   const std::vector<art::Ptr<recob::Hit>>& hits) const;

  /// Build a shower segment from sets of hits and attached to the provided
  /// vertex.
  pma::Track3D* buildShowerSeg(const detinfo::DetectorPropertiesData& detProp,
                               const std::vector<art::Ptr<recob::Hit>>& hits,
                               const pma::Vector3D& vtx) const;

  /// Build a straight segment from two sets of hits (they should origin from
  /// two wire planes); method is intendet for short tracks or shower initial
  /// parts, where only a few hits per plane are available and there is no
  /// chance to see a curvature or any other features.
  pma::Track3D* buildSegment(const detinfo::DetectorPropertiesData& clockData,
                             const std::vector<art::Ptr<recob::Hit>>& hits_1,
                             const std::vector<art::Ptr<recob::Hit>>& hits_2 = {}) const;

  /// Build a straight segment from two sets of hits (they should origin from
  /// two wire planes), starting from a given point (like vertex known from
  /// another algorithm); method is intendet for short tracks or shower initial
  /// parts, where only a few hits per plane are available and there is no
  /// chance to see a curvature or any other features.
  pma::Track3D* buildSegment(const detinfo::DetectorPropertiesData& clockData,
                             const std::vector<art::Ptr<recob::Hit>>& hits_1,
                             const std::vector<art::Ptr<recob::Hit>>& hits_2,
                             const TVector3& point) const;

  /// Build a straight segment from set of hits (they should origin from two
  /// wire planes at least), starting from a given point.
  pma::Track3D* buildSegment(const detinfo::DetectorPropertiesData& detProp,
                             const std::vector<art::Ptr<recob::Hit>>& hits,
                             const TVector3& point) const;

  /// Get rid of small groups of hits around cascades; used to calculate cascade
  /// starting direction using the compact core cluster.
  void FilterOutSmallParts(const detinfo::DetectorPropertiesData& detProp,
                           double r2d,
                           const std::vector<art::Ptr<recob::Hit>>& hits_in,
                           std::vector<art::Ptr<recob::Hit>>& hits_out,
                           const TVector2& vtx2d) const;

  void RemoveNotEnabledHits(pma::Track3D& trk) const;

  /// Add more hits to an existing track, reoptimize, optionally add more nodes.
  pma::Track3D* extendTrack(const detinfo::DetectorPropertiesData& clockData,
                            const pma::Track3D& trk,
                            const std::vector<art::Ptr<recob::Hit>>& hits,
                            bool add_nodes) const;

  /// Add 3D reference points to clean endpoints of a track (both need to be in
  /// the same TPC).
  void guideEndpoints(const detinfo::DetectorPropertiesData& clockData,
                      pma::Track3D& trk,
                      const std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>>& hits) const;

  /// Add 3D reference points to clean endpoint of a track.
  void guideEndpoints(const detinfo::DetectorPropertiesData& clockData,
                      pma::Track3D& trk,
                      pma::Track3D::ETrackEnd endpoint,
                      const std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>>& hits) const;

  std::vector<pma::Hit3D*> trimTrackToVolume(pma::Track3D& trk, TVector3 p0, TVector3 p1) const;

  /// Flip tracks to get second as a continuation of first; returns false if not
  /// possible (tracks in reversed order).
  bool alignTracks(pma::Track3D& first, pma::Track3D& second) const;

  /// Add src to dst as it was its continuation; nodes of src are added to dst
  /// after its own nodes, hits of src are added to hits of dst, then dst is
  /// reoptimized.
  void mergeTracks(const detinfo::DetectorPropertiesData& detProp,
                   pma::Track3D& dst,
                   pma::Track3D& src,
                   bool reopt) const;

  /// Try to correct track direction of the stopping particle:
  ///   dir: kForward  - particle stop is at the end of the track;
  ///        kBackward - particle stop is at the beginning of the track;
  /// dQ/dx difference has to be above thr to actually flip the track;
  /// compares dQ/dx of n hits at each end of the track (default is based on the track length).
  void
  autoFlip(pma::Track3D& trk,
           pma::Track3D::EDirection dir = Track3D::kForward,
           double thr = 0.0,
           unsigned int n = 0) const
  {
    trk.AutoFlip(dir, thr, n);
  };

  /// Intendet to calculate dQ/dx in the initial part of EM cascade; collection
  /// view is used by default, but it works also with other projections.
  double selectInitialHits(pma::Track3D& trk,
                           unsigned int view = geo::kZ,
                           unsigned int* nused = 0) const;

private:
  // Helpers for guideEndpoints
  bool chkEndpointHits_(const detinfo::DetectorPropertiesData& detProp,
                        int wire,
                        int wdir,
                        double drift_x,
                        int view,
                        unsigned int tpc,
                        unsigned int cryo,
                        const pma::Track3D& trk,
                        const std::vector<art::Ptr<recob::Hit>>& hits) const;
  bool addEndpointRef_(const detinfo::DetectorPropertiesData& detProp,
                       pma::Track3D& trk,
                       const std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>>& hits,
                       std::pair<int, int> const* wires,
                       double const* xPos,
                       unsigned int tpc,
                       unsigned int cryo) const;

  // Helpers for FilterOutSmallParts
  bool GetCloseHits_(const detinfo::DetectorPropertiesData& detProp,
                     double r2d,
                     const std::vector<art::Ptr<recob::Hit>>& hits_in,
                     std::vector<size_t>& used,
                     std::vector<art::Ptr<recob::Hit>>& hits_out) const;

  bool Has_(const std::vector<size_t>& v, size_t idx) const;

  // Make segment shorter depending on mse
  void ShortenSeg_(const detinfo::DetectorPropertiesData& detProp,
                   pma::Track3D& trk,
                   const geo::TPCGeo& tpcgeom) const;

  // Control length of the track and number of hits which are still enabled
  bool TestTrk_(pma::Track3D& trk, const geo::TPCGeo& tpcgeom) const;

  // Calculate good number of segments depending on the number of hits.
  static size_t getSegCount_(size_t trk_size);

  // Parameters used in the algorithm

  double const fOptimizationEps; // relative change in the obj.function that
                                 // ends optimization, then next nodes are added
                                 // or track building is finished

  double const fFineTuningEps; // relative change in the obj.function that ends
                               // final tuning

  double const fTrkValidationDist2D; // max. distance [cm] used in the track
                                     // validation in the "third" plane
  double const fHitTestingDist2D;    // max. distance [cm] used in testing comp. of
                                     // hits with the track

  double const fMinTwoViewFraction; // min. length fraction covered with multiple 2D
                                    // view hits intertwinted with each other

  // Geometry and detector properties
  geo::GeometryCore const* fGeom;
};

#endif

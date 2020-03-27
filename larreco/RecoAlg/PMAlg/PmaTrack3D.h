/**
 *  @file   PmaTrack3D.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 *          Based on the "Precise 3D track reco..." AHEP (2013) 260820, with all the tricks that we
 *          developed later and with the work for the track-vertex topology optimization done at FNAL.
 *
 *          Progress:
 *             May-June 2015:  basic functionality for single (not-branching) 3D track optimization and dQ/dx.
 *             August 2015:    3D vertex finding and branching tracks optimization.
 */

#ifndef PmaTrack3D_h
#define PmaTrack3D_h

#include "TVector2.h"
#include "TVector3.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/PMAlg/PmaHit3D.h"
#include "larreco/RecoAlg/PMAlg/PmaNode3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace pma {
  class Segment3D;
  class Track3D;
}

class pma::Track3D {
public:
  enum ETrackEnd { kBegin = -1, kEnd = 1 };
  enum EDirection { kForward = -1, kBackward = 1 };
  enum ETag {
    kNotTagged = 0,
    kTrackLike = 0,
    kEmLike = 1,
    kStopping = 2,
    kCosmic = 4,

    kGeometry_YY = 0x000100,
    kGeometry_YZ = 0x000200,
    kGeometry_ZZ = 0x000300,
    kGeometry_XX = 0x000400,
    kGeometry_XY = 0x000500,
    kGeometry_XZ = 0x000600,

    kGeometry_Y = 0x001000,
    kGeometry_Z = 0x002000,
    kGeometry_X = 0x003000,

    kOutsideDrift_Partial = 0x010000,
    kOutsideDrift_Complete = 0x020000,
    kBeamIncompatible = 0x030000
  };
  ETag
  GetTag() const noexcept
  {
    return fTag;
  }
  bool
  HasTagFlag(ETag value) const noexcept
  {
    return (fTag & value);
  }
  void
  SetTagFlag(ETag value)
  {
    fTag = (ETag)(fTag | value);
  }

  Track3D();
  //  Track3D(const Track3D& src);
  ~Track3D();

  bool Initialize(detinfo::DetectorPropertiesData const& detProp, float initEndSegW = 0.05F);

  pma::Hit3D* release_at(size_t index);
  void
  push_back(pma::Hit3D* hit)
  {
    hit->fParent = this;
    fHits.push_back(hit);
  }
  bool push_back(const detinfo::DetectorPropertiesData& detProp, const art::Ptr<recob::Hit>& hit);
  bool erase(const art::Ptr<recob::Hit>& hit);

  pma::Hit3D* operator[](size_t index) { return fHits[index]; }
  pma::Hit3D const* operator[](size_t index) const { return fHits[index]; }
  pma::Hit3D const*
  front() const
  {
    return fHits.front();
  }
  pma::Hit3D const*
  back() const
  {
    return fHits.back();
  }
  size_t
  size() const
  {
    return fHits.size();
  }

  int index_of(const pma::Hit3D* hit) const;
  int index_of(const pma::Node3D* n) const;

  double
  Length(size_t step = 1) const
  {
    return Length(0, size() - 1, step);
  }
  double Length(size_t start, size_t stop, size_t step = 1) const;

  double Dist2(const TVector2& p2d, unsigned int view, unsigned int tpc, unsigned int cryo) const;
  double Dist2(const TVector3& p3d) const;

  /// Get trajectory direction at given hit index.
  pma::Vector3D GetDirection3D(size_t index) const;

  /// Add hits; does not update hit->node/seg assignments nor hit projection to
  /// track, so MakeProjection() and SortHits() should be called as needed.
  void AddHits(detinfo::DetectorPropertiesData const& detProp,
               const std::vector<art::Ptr<recob::Hit>>& hits);

  /// Remove hits; removes also hit->node/seg assignments.
  void RemoveHits(const std::vector<art::Ptr<recob::Hit>>& hits);

  unsigned int NHits(unsigned int view) const;
  unsigned int NEnabledHits(unsigned int view = geo::kUnknown) const;
  bool HasTwoViews(size_t nmin = 1) const;

  std::vector<unsigned int> TPCs() const;
  std::vector<unsigned int> Cryos() const;

  unsigned int
  FrontTPC() const
  {
    return fNodes.front()->TPC();
  }
  unsigned int
  FrontCryo() const
  {
    return fNodes.front()->Cryo();
  }

  unsigned int
  BackTPC() const
  {
    return fNodes.back()->TPC();
  }
  unsigned int
  BackCryo() const
  {
    return fNodes.back()->Cryo();
  }

  bool
  HasTPC(int tpc) const
  {
    for (auto n : fNodes)
      if (n->TPC() == tpc) return true;
    return false;
  }

  /// Rectangular region of the track 2D projection in view/tpc/cryo; first in
  /// the returned pair is (min_wire; min_drift), second is (max_wire;
  /// max_drift). Used for preselection of neighbouring hits in the track
  /// validation functions.
  std::pair<TVector2, TVector2> WireDriftRange(detinfo::DetectorPropertiesData const& detProp,
                                               unsigned int view,
                                               unsigned int tpc,
                                               unsigned int cryo) const;

  /// Invert the order of hits and vertices in the track, break other tracks if
  /// needed (new tracks are added to the allTracks vector). Returns true if
  /// successful or false if any of required track flips was not possible (e.g.
  /// resulting track would be composed of hits from a single 2D projection).
  bool Flip(const detinfo::DetectorPropertiesData& detProp, std::vector<pma::Track3D*>& allTracks);

  /// Invert the order of hits and vertices in the track, will fail on
  /// configuration that causes breaking another track.
  void Flip();

  /// Check if the track can be flipped without breaking any other track.
  bool CanFlip() const;

  void AutoFlip(pma::Track3D::EDirection dir, double thr = 0.0, unsigned int n = 0);
  bool AutoFlip(detinfo::DetectorPropertiesData const& detProp,
                std::vector<pma::Track3D*>& allTracks,
                pma::Track3D::EDirection dir,
                double thr = 0.0,
                unsigned int n = 0);

  /// MSE of 2D hits.
  double TestHitsMse(detinfo::DetectorPropertiesData const& detProp,
                     const std::vector<art::Ptr<recob::Hit>>& hits,
                     bool normalized = true) const; // normalize to the number of hits

  /// Count close 2D hits.
  unsigned int TestHits(detinfo::DetectorPropertiesData const& detProp,
                        const std::vector<art::Ptr<recob::Hit>>& hits,
                        double dist = 0.4) const; // max acceptable distance [cm]

  int NextHit(int index, unsigned int view = geo::kZ, bool inclDisabled = false) const;
  int PrevHit(int index, unsigned int view = geo::kZ, bool inclDisabled = false) const;

  /// Length of the track part associated with index'th hit. Calculated as a half distance to
  /// the preceding hit plus half distance to the subsequent hit. In case of the first (last)
  /// hit - missing part is estimated as 1/4 of the distance to the next (previous) hit.
  /// NOTE: only hits from a given view are considered; other hits are accounted for
  /// segment lengths but overall dx is calculated between hits in given view.
  double HitDxByView(size_t index, unsigned int view) const;

  /// Sequence of <hit_index, (wire, drift, X, Y, Z, dE, dx, range)> values for the track,
  /// hits tagged as outliers are skipped by default.
  /** Results are pushed into the dedx vector given in the function arguments:

            hit (segment middle if many hits) 2D projection in view:
              dedx[n][0] = wire;
              dedx[n][1] = drift;

            hit (segment middle if many hits) 3D position [cm]:
              dedx[n][2] = X;
              dedx[n][3] = Y;
              dedx[n][4] = Z;

              dedx[n][5] = dE [now ADC], energy assigned to the segment;

              dedx[n][6] = dx [cm], length of the segment.

              dedx[n][7] = range, total length to the track endpoint;

            Parameters:
              dedx  - vector to store results (empty at the begining);
              view  - view (U, V or Z) from which dedx is created;
              skip  - number of hits to skip at the begining (first hit has poorly estimated segment
                      length so it can be convenient to set skip=1 and handle first hit charge manually);
              inclDisabled - if true then artificial hits added with CompleteMissingWires() are used,
                             otherwise only true hits found in ADC are used.

            Return value: sum of ADC's of hits skipped at the begining. */
  double GetRawdEdxSequence(std::map<size_t, std::vector<double>>& dedx,
                            unsigned int view = geo::kZ,
                            unsigned int skip = 0,
                            bool inclDisabled = false) const;

  std::vector<float> DriftsOfWireIntersection(detinfo::DetectorPropertiesData const& detProp,
                                              unsigned int wire,
                                              unsigned int view) const;
  size_t CompleteMissingWires(detinfo::DetectorPropertiesData const& detProp, unsigned int view);

  void
  AddRefPoint(const TVector3& p)
  {
    fAssignedPoints.push_back(new TVector3(p));
  }
  void
  AddRefPoint(double x, double y, double z)
  {
    fAssignedPoints.push_back(new TVector3(x, y, z));
  }
  bool HasRefPoint(TVector3* p) const;

  /// MSE of hits weighted with hit amplidudes and wire plane coefficients.
  double GetMse(unsigned int view = geo::kUnknown) const;

  /// Objective function optimized in track reconstruction.
  double GetObjFunction(float penaltyFactor = 1.0F) const;

  /// Main optimization method.
  double Optimize(const detinfo::DetectorPropertiesData& detProp,
                  int nNodes = -1,
                  double eps = 0.01,
                  bool selAllHits = true,
                  bool setAllNodes = true,
                  size_t selSegHits = 0,
                  size_t selVtxHits = 0);

  void SortHitsInTree(bool skipFirst = false);
  void MakeProjectionInTree(bool skipFirst = false);
  bool UpdateParamsInTree(bool skipFirst, size_t& depth);
  double GetObjFnInTree(bool skipFirst = false);
  double TuneSinglePass(bool skipFirst = false);
  double TuneFullTree(double eps = 0.001, double gmax = 50.0);

  /// Adjust track tree position in the drift direction (when T0 is being
  /// corrected).
  void ApplyDriftShiftInTree(const detinfo::DetectorClocksData& clockData,
                             detinfo::DetectorPropertiesData const& detProp,
                             double dx,
                             bool skipFirst = false);
  /// Function to convert dx into dT0
  void SetT0FromDx(const detinfo::DetectorClocksData& clockData,
                   detinfo::DetectorPropertiesData const& detProp,
                   double dx);
  double
  GetT0() const
  {
    return fT0;
  }
  /// Check if the T0 has been set - enables us to distinguish between T0 set
  /// very close to zero or not set.
  bool
  HasT0() const noexcept
  {
    return fT0Flag;
  }

  /// Cut out tails with no hits assigned.
  void CleanupTails();

  /// Move the first/last Node3D to the first/last hit in the track;
  /// returns true if all OK, false if empty segments found.
  bool ShiftEndsToHits();

  std::vector<pma::Segment3D*> const&
  Segments() const noexcept
  {
    return fSegments;
  }

  pma::Segment3D* NextSegment(pma::Node3D* vtx) const;
  pma::Segment3D* PrevSegment(pma::Node3D* vtx) const;

  std::vector<pma::Node3D*> const&
  Nodes() const noexcept
  {
    return fNodes;
  }
  pma::Node3D*
  FirstElement() const
  {
    return fNodes.front();
  }
  pma::Node3D*
  LastElement() const
  {
    return fNodes.back();
  }

  void AddNode(pma::Node3D* node);
  void
  AddNode(detinfo::DetectorPropertiesData const& detProp,
          TVector3 const& p3d,
          unsigned int tpc,
          unsigned int cryo)
  {
    double ds = fNodes.empty() ? 0 : fNodes.back()->GetDriftShift();
    AddNode(new pma::Node3D(detProp, p3d, tpc, cryo, false, ds));
  }
  bool AddNode(detinfo::DetectorPropertiesData const& detProp);

  void InsertNode(detinfo::DetectorPropertiesData const& detProp,
                  TVector3 const& p3d,
                  size_t at_idx,
                  unsigned int tpc,
                  unsigned int cryo);
  bool RemoveNode(size_t idx);

  pma::Track3D* Split(detinfo::DetectorPropertiesData const& detProp,
                      size_t idx,
                      bool try_start_at_idx = true);

  bool AttachTo(pma::Node3D* vStart, bool noFlip = false);
  bool AttachBackTo(pma::Node3D* vStart);
  bool IsAttachedTo(pma::Track3D const* trk) const;

  /// Extend the track with everything from src, delete the src;
  void ExtendWith(pma::Track3D* src);

  pma::Track3D* GetRoot();
  bool GetBranches(std::vector<pma::Track3D const*>& branches, bool skipFirst = false) const;

  void MakeProjection();
  void UpdateProjection();
  void SortHits();

  unsigned int DisableSingleViewEnds();
  bool SelectHits(float fraction = 1.0F);
  bool SelectRndHits(size_t segmax, size_t vtxmax);
  bool SelectAllHits();

  float
  GetEndSegWeight() const noexcept
  {
    return fEndSegWeight;
  }
  void
  SetEndSegWeight(float value) noexcept
  {
    fEndSegWeight = value;
  }

  float
  GetPenalty() const noexcept
  {
    return fPenaltyFactor;
  }
  void
  SetPenalty(float value) noexcept
  {
    fPenaltyFactor = value;
  }

  unsigned int
  GetMaxHitsPerSeg() const noexcept
  {
    return fMaxHitsPerSeg;
  }
  void
  SetMaxHitsPerSeg(unsigned int value) noexcept
  {
    fMaxHitsPerSeg = value;
  }

private:
  void ClearNodes();
  void MakeFastProjection();

  bool AttachToSameTPC(pma::Node3D* vStart);
  bool AttachToOtherTPC(pma::Node3D* vStart);

  bool AttachBackToSameTPC(pma::Node3D* vStart);
  bool AttachBackToOtherTPC(pma::Node3D* vStart);

  void InternalFlip(std::vector<pma::Track3D*>& toSort);

  void UpdateHitsRadius();
  double AverageDist2() const;

  bool InitFromHits(detinfo::DetectorPropertiesData const& detProp,
                    int tpc,
                    int cryo,
                    float initEndSegW = 0.05F);
  bool InitFromRefPoints(detinfo::DetectorPropertiesData const& detProp, int tpc, int cryo);
  void InitFromMiddle(detinfo::DetectorPropertiesData const& detProp, int tpc, int cryo);

  pma::Track3D* GetNearestTrkInTree(const TVector3& p3d_cm, double& dist, bool skipFirst = false);
  pma::Track3D* GetNearestTrkInTree(const TVector2& p2d_cm,
                                    unsigned int view,
                                    unsigned int tpc,
                                    unsigned int cryo,
                                    double& dist,
                                    bool skipFirst = false);
  void ReassignHitsInTree(pma::Track3D* plRoot = nullptr);

  /// Distance to the nearest subsequent (dir = Track3D::kForward) or preceeding
  /// (dir = Track3D::kBackward) hit in given view. In case of last (first) hit
  /// in this view the half-distance in opposite direction is returned.
  /// Parameter secondDir is only for internal protection - please leave the
  /// default value.
  double HitDxByView(size_t index,
                     unsigned int view,
                     Track3D::EDirection dir,
                     bool secondDir = false) const;

  /// Calculate 3D position corresponding to 2D hit, return true if the 3D point
  /// is in the same TPC as the hit, false otherwise. Calculates also distance^2
  /// between the hit and 2D projection of the track. NOTE: results are
  /// meaningful only if the function returns true.
  bool GetUnconstrainedProj3D(detinfo::DetectorPropertiesData const& detProp,
                              art::Ptr<recob::Hit> hit,
                              TVector3& p3d,
                              double& dist2) const;

  void DeleteSegments();
  void RebuildSegments();
  bool SwapVertices(size_t v0, size_t v1);
  bool UpdateParams();

  bool CheckEndSegment(pma::Track3D::ETrackEnd endCode);

  std::vector<pma::Hit3D*> fHits;

  std::vector<TVector3*> fAssignedPoints;

  pma::Element3D* GetNearestElement(const TVector2& p2d,
                                    unsigned int view,
                                    int tpc = -1,
                                    bool skipFrontVtx = false,
                                    bool skipBackVtx = false) const;
  pma::Element3D* GetNearestElement(const TVector3& p3d) const;

  std::vector<pma::Node3D*> fNodes;
  std::vector<pma::Segment3D*> fSegments;

  unsigned int fMaxHitsPerSeg{70};
  float fPenaltyFactor{1.0F};
  float fMaxSegStopFactor{8.0F};

  unsigned int fSegStopValue{2};
  unsigned int fMinSegStop{2};
  unsigned int fMaxSegStop{2};

  float fSegStopFactor{0.2F};
  float fPenaltyValue{0.1F};
  float fEndSegWeight{0.05F};
  float fHitsRadius{1.0F};

  double fT0{};
  bool fT0Flag{false};

  ETag fTag{kNotTagged};
};

#endif

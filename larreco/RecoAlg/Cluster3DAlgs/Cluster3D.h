/**
 *
 *  @brief Definition of utility objects for use in the 3D clustering for LArSoft
 *
 *         The objects defined in this module are intended for internal use by the
 *         3D clustering (see Cluster3D_module.cc in larreco/ClusterFinder).
 *         These objects mostly contain volatile information and are not suitable
 *         for storage in the art event store
 *
 *  @author usher@slac.stanford.edu
 *
 */

#ifndef RECO_CLUSTER3D_H
#define RECO_CLUSTER3D_H

#include <iosfwd>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/DCEL.h"
namespace recob {
  class Hit;
}

// Eigen
#include <Eigen/Core>

namespace reco {

  // First define a container object to augment the sparse 2D hit information
  class ClusterHit2D {
  public:
    ClusterHit2D(); // Default constructor

  private:
    mutable unsigned m_statusBits; ///< Volatile status information of this 3D hit
    mutable float m_docaToAxis;    ///< DOCA of hit at POCA to associated cluster axis
    mutable float m_arcLenToPoca;  ///< arc length to POCA along cluster axis
    float m_xPosition;             ///< The x coordinate for this hit
    float m_timeTicks;             ///< The time (in ticks) for this hit
    geo::WireID m_wireID;          ///< Keep track this particular hit's wireID
    const recob::Hit* m_hit;       ///< Hit we are augmenting

  public:
    enum StatusBits {
      SHAREDINPAIR = 0x00080000,
      SHAREDINTRIPLET = 0x00040000,
      USEDINPAIR = 0x00008000,
      USEDINTRIPLET = 0x00004000,
      SHAREDINCLUSTER = 0x00000200,
      USEDINCLUSTER = 0x00000100,
      USED = 0x00000001
    };

    ClusterHit2D(unsigned statusBits,
                 float doca,
                 float poca,
                 float xPosition,
                 float timeTicks,
                 const geo::WireID& wireID,
                 const recob::Hit* recobHit);

    ClusterHit2D(const ClusterHit2D&);

    unsigned getStatusBits() const { return m_statusBits; }
    float getDocaToAxis() const { return m_docaToAxis; }
    float getArcLenToPoca() const { return m_arcLenToPoca; }
    float getXPosition() const { return m_xPosition; }
    float getTimeTicks() const { return m_timeTicks; }
    const geo::WireID& WireID() const { return m_wireID; }
    const recob::Hit* getHit() const { return m_hit; }

    void setStatusBit(unsigned bits) const { m_statusBits |= bits; }
    void clearStatusBits(unsigned bits) const { m_statusBits &= ~bits; }
    void setDocaToAxis(float doca) const { m_docaToAxis = doca; }
    void setArcLenToPoca(float poca) const { m_arcLenToPoca = poca; }

    void setHit(const recob::Hit* hit) { m_hit = hit; }

    friend std::ostream& operator<<(std::ostream& o, const ClusterHit2D& c);
    friend bool operator<(const ClusterHit2D& a, const ClusterHit2D& b);
  };

  using ClusterHit2DVec = std::vector<const reco::ClusterHit2D*>;

  // Now define an object with the recob::Hit information that will comprise the 3D cluster
  class ClusterHit3D {
  public:
    enum StatusBits {
      REJECTEDHIT = 0x80000000,     ///< Hit has been rejected for any reason
      SKELETONHIT = 0x10000000,     ///< Hit is a "skeleton" hit
      EDGEHIT = 0x20000000,         ///< Hit is an "edge" hit
      SEEDHIT = 0x40000000,         ///< Hit is part of Seed for track fits
      MADESPACEPOINT = 0x08000000,  ///< Hit has been made into Space Point
      CONVEXHULLVTX = 0x04000000,   ///< Point is on primary cluster convex hull
      EXTREMEPOINT = 0x02000000,    ///< Is a convex hull extreme point
      SKELETONPOSAVE = 0x00100000,  ///< Skeleton hit position averaged
      CLUSTERVISITED = 0x00008000,  ///< "visited" by a clustering algorithm
      CLUSTERNOISE = 0x00004000,    ///< Labelled "noise" by a clustering algorithm
      CLUSTERATTACHED = 0x00002000, ///< attached to a cluster
      CLUSTERSHARED = 0x00001000,   ///< 3D hit has 2D hits shared between clusters
      PATHCHECKED = 0x00000800,     ///< Path checking algorithm has seen this hit
      SELECTEDBYMST = 0x00000100,   ///< Hit has been used in Cluster Splitting MST
      PCAOUTLIER = 0x00000080,      ///< Hit labelled outlier in PCA
      HITINVIEW0 = 0x00000001,      ///< Hit contains 2D hit from view 0 (u plane)
      HITINVIEW1 = 0x00000002,      ///< Hit contains 2D hit from view 1 (v plane)
      HITINVIEW2 = 0x00000004       ///< Hit contains 2D hit from view 2 (w plane)
    };

    ClusterHit3D(); // Default constructor

    ClusterHit3D(size_t id,
                 unsigned int statusBits,
                 const Eigen::Vector3f& position,
                 float totalCharge,
                 float avePeakTime,
                 float deltaPeakTime,
                 float sigmaPeakTime,
                 float hitChiSquare,
                 float overlapFraction,
                 float chargeAsymmetry,
                 float docaToAxis,
                 float arclenToPoca,
                 const ClusterHit2DVec& hitVec,
                 const std::vector<float>& hitDelTSigVec,
                 const std::vector<geo::WireID>& wireIDVec);

    ClusterHit3D(const ClusterHit3D&);

    void initialize(size_t id,
                    unsigned int statusBits,
                    const Eigen::Vector3f& position,
                    float totalCharge,
                    float avePeakTime,
                    float deltaPeakTime,
                    float sigmaPeakTime,
                    float hitChiSquare,
                    float overlapFraction,
                    float chargeAsymmetry,
                    float docaToAxis,
                    float arclenToPoca,
                    const ClusterHit2DVec& hitVec,
                    const std::vector<float>& hitDelTSigVec,
                    const std::vector<geo::WireID>& wireIDVec);

    size_t getID() const { return fID; }
    unsigned int getStatusBits() const { return fStatusBits; }
    const Eigen::Vector3f getPosition() const { return fPosition; }
    float getX() const { return fPosition[0]; }
    float getY() const { return fPosition[1]; }
    float getZ() const { return fPosition[2]; }
    float getTotalCharge() const { return fTotalCharge; }
    float getAvePeakTime() const { return fAvePeakTime; }
    float getDeltaPeakTime() const { return fDeltaPeakTime; }
    float getSigmaPeakTime() const { return fSigmaPeakTime; }
    float getHitChiSquare() const { return fHitChiSquare; }
    float getOverlapFraction() const { return fOverlapFraction; }
    float getChargeAsymmetry() const { return fChargeAsymmetry; }
    float getDocaToAxis() const { return fDocaToAxis; }
    float getArclenToPoca() const { return fArclenToPoca; }
    const ClusterHit2DVec& getHits() const { return fHitVector; }
    const std::vector<float> getHitDelTSigVec() const { return fHitDelTSigVec; }
    const std::vector<geo::WireID>& getWireIDs() const { return fWireIDVector; }

    ClusterHit2DVec& getHits() { return fHitVector; }

    bool bitsAreSet(const unsigned int& bitsToCheck) const { return fStatusBits & bitsToCheck; }

    void setID(const size_t& id) const { fID = id; }
    void setStatusBit(unsigned bits) const { fStatusBits |= bits; }
    void clearStatusBits(unsigned bits) const { fStatusBits &= ~bits; }
    void setDocaToAxis(double doca) const { fDocaToAxis = doca; }
    void setArclenToPoca(double poca) const { fArclenToPoca = poca; }
    void setWireID(const geo::WireID& wid) const;

    void setPosition(const Eigen::Vector3f& pos) const { fPosition = pos; }

    const bool operator<(const reco::ClusterHit3D& other) const
    {
      if (fPosition[2] != other.fPosition[2])
        return fPosition[2] < other.fPosition[2];
      else
        return fPosition[0] < other.fPosition[0];
    }

    const bool operator==(const reco::ClusterHit3D& other) const { return fID == other.fID; }

    friend std::ostream& operator<<(std::ostream& o, const ClusterHit3D& c);
    //friend bool          operator <  (const ClusterHit3D & a, const ClusterHit3D & b);

  private:
    mutable size_t fID;                ///< "id" of this hit (useful for indexing)
    mutable unsigned int fStatusBits;  ///< Volatile status information of this 3D hit
    mutable Eigen::Vector3f fPosition; ///< position of this hit combination in world coordinates
    float fTotalCharge;                ///< Sum of charges of all associated recob::Hits
    float fAvePeakTime;                ///< Average peak time of all associated recob::Hits
    float fDeltaPeakTime;              ///< Largest delta peak time of associated recob::Hits
    float fSigmaPeakTime;              ///< Quad sum of peak time sigmas
    float fHitChiSquare;               ///< Hit ChiSquare relative to the average time
    float fOverlapFraction;            ///< Hit overlap fraction start/stop of triplet
    float fChargeAsymmetry;            ///< Assymetry of average of two closest to third charge
    mutable float fDocaToAxis;         ///< DOCA to the associated cluster axis
    mutable float fArclenToPoca;       ///< arc length along axis to DOCA point
    ClusterHit2DVec fHitVector;        ///< Hits comprising this 3D hit
    mutable std::vector<float> fHitDelTSigVec;      ///< Delta t of hit to matching pair / sig
    mutable std::vector<geo::WireID> fWireIDVector; ///< Wire ID's for the planes making up hit
  };

  // We also need to define a container for the output of the PCA Analysis
  class PrincipalComponents {
  public:
    using EigenValues = Eigen::Vector3f;
    using EigenVectors = Eigen::Matrix3f;

    PrincipalComponents();

  private:
    bool m_svdOK;                  ///< SVD Decomposition was successful
    int m_numHitsUsed;             ///< Number of hits in the decomposition
    EigenValues m_eigenValues;     ///< Eigen values from SVD decomposition
    EigenVectors m_eigenVectors;   ///< The three principle axes
    Eigen::Vector3f m_avePosition; ///< Average position of hits fed to PCA
    mutable double m_aveHitDoca;   ///< Average doca of hits used in PCA

  public:
    PrincipalComponents(bool ok,
                        int nHits,
                        const EigenValues& eigenValues,
                        const EigenVectors& eigenVecs,
                        const Eigen::Vector3f& avePos,
                        const float aveHitDoca = 9999.);

    bool getSvdOK() const { return m_svdOK; }
    int getNumHitsUsed() const { return m_numHitsUsed; }
    const EigenValues& getEigenValues() const { return m_eigenValues; }
    const EigenVectors& getEigenVectors() const { return m_eigenVectors; }
    const Eigen::Vector3f& getAvePosition() const { return m_avePosition; }
    const float getAveHitDoca() const { return m_aveHitDoca; }

    void flipAxis(size_t axis);
    void setAveHitDoca(double doca) const { m_aveHitDoca = doca; }

    friend std::ostream& operator<<(std::ostream& o, const PrincipalComponents& a);
    friend bool operator<(const PrincipalComponents& a, const PrincipalComponents& b);
  };

  class Cluster3D {
  public:
    Cluster3D(); ///Default constructor

  private:
    mutable unsigned m_statusBits;    ///< Status bits for the cluster
    PrincipalComponents m_pcaResults; ///< Output of the prinicipal componenets analysis
    float m_totalCharge;              ///< Total charge in the cluster
    float m_startPosition[3];         ///< "start" position for cluster (world coordinates)
    float m_endPosition[3];           ///< "end" position for cluster
    int m_clusterIdx;                 ///< ID for this cluster

  public:
    Cluster3D(unsigned statusBits,
              const PrincipalComponents& pcaResults,
              float totalCharge,
              const float* startPosition,
              const float* endPosition,
              int idx);

    unsigned getStatusBits() const { return m_statusBits; }
    const PrincipalComponents& getPcaResults() const { return m_pcaResults; }
    float getTotalCharge() const { return m_totalCharge; }
    const float* getStartPosition() const { return m_startPosition; }
    const float* getEndPosition() const { return m_endPosition; }
    int getClusterIdx() const { return m_clusterIdx; }

    void setStatusBit(unsigned bits) const { m_statusBits |= bits; }
    void clearStatusBits(unsigned bits) const { m_statusBits &= ~bits; }

    Cluster3D operator+(Cluster3D);
    friend std::ostream& operator<<(std::ostream& o, const Cluster3D& c);
    friend bool operator<(const Cluster3D& a, const Cluster3D& b);
  };

  /**
 *  @brief A utility class used in construction of 3D clusters
 */
  class RecobClusterParameters {
  public:
    RecobClusterParameters()
      : m_startTime(999999.)
      , m_sigmaStartTime(1.)
      , m_endTime(0.)
      , m_sigmaEndTime(1.)
      , m_totalCharge(0.)
      , m_startWire(9999999)
      , m_endWire(0)
      , m_plane(geo::PlaneID())
      , m_view(geo::kUnknown)
    {
      m_hitVector.clear();
    }

    void UpdateParameters(const reco::ClusterHit2D* hit);

    float m_startTime;
    float m_sigmaStartTime;
    float m_endTime;
    float m_sigmaEndTime;
    float m_totalCharge;
    unsigned int m_startWire;
    unsigned int m_endWire;
    geo::PlaneID m_plane;
    geo::View_t m_view;
    ClusterHit2DVec m_hitVector;
  };

  /**
 *  @brief export some data structure definitions
 */
  using Hit2DListPtr = std::list<const reco::ClusterHit2D*>;
  using HitPairListPtr = std::list<const reco::ClusterHit3D*>;
  using HitPairSetPtr = std::set<const reco::ClusterHit3D*>;
  using HitPairListPtrList = std::list<HitPairListPtr>;
  using HitPairClusterMap = std::map<int, HitPairListPtr>;
  using HitPairList = std::list<reco::ClusterHit3D>;
  //using HitPairList              = std::list<std::unique_ptr<reco::ClusterHit3D>>;

  using PCAHitPairClusterMapPair =
    std::pair<reco::PrincipalComponents, reco::HitPairClusterMap::iterator>;
  using PlaneToClusterParamsMap = std::map<size_t, RecobClusterParameters>;
  using EdgeTuple = std::tuple<const reco::ClusterHit3D*, const reco::ClusterHit3D*, double>;
  using EdgeList = std::list<EdgeTuple>;
  using Hit3DToEdgePair = std::pair<const reco::ClusterHit3D*, reco::EdgeList>;
  using Hit3DToEdgeMap = std::unordered_map<const reco::ClusterHit3D*, reco::EdgeList>;
  using Hit2DToHit3DListMap = std::unordered_map<const reco::ClusterHit2D*, reco::HitPairListPtr>;
  //using VertexPoint              = Eigen::Vector3f;
  //using VertexPointList          = std::list<Eigen::Vector3f>;

  using ProjectedPoint = std::
    tuple<float, float, const reco::ClusterHit3D*>; ///< Projected coordinates and pointer to hit
  using ProjectedPointList = std::list<ProjectedPoint>;
  using ConvexHullKinkTuple = std::
    tuple<ProjectedPoint, Eigen::Vector2f, Eigen::Vector2f>; ///< Point plus edges that point to it
  using ConvexHullKinkTupleList = std::list<ConvexHullKinkTuple>;

  /**
 *  @brief Define a container for working with the convex hull
 */
  class ConvexHull {
  public:
    ConvexHull()
    {
      fProjectedPointList.clear(), fConvexHullPointList.clear(), fConvexHullEdgeMap.clear(),
        fConvexHullEdgeList.clear(), fConvexHullExtremePoints.clear(),
        fConvexHullKinkPoints.clear();
    }

    void clear()
    {
      fProjectedPointList.clear(), fConvexHullPointList.clear(), fConvexHullEdgeMap.clear(),
        fConvexHullEdgeList.clear(), fConvexHullExtremePoints.clear(),
        fConvexHullKinkPoints.clear();
    }

    reco::ProjectedPointList& getProjectedPointList() { return fProjectedPointList; }
    reco::ProjectedPointList& getConvexHullPointList() { return fConvexHullPointList; }
    reco::Hit3DToEdgeMap& getConvexHullEdgeMap() { return fConvexHullEdgeMap; }
    reco::EdgeList& getConvexHullEdgeList() { return fConvexHullEdgeList; }
    reco::ProjectedPointList& getConvexHullExtremePoints() { return fConvexHullExtremePoints; }
    reco::ConvexHullKinkTupleList& getConvexHullKinkPoints() { return fConvexHullKinkPoints; }

  private:
    reco::ProjectedPointList
      fProjectedPointList; ///< The input set of points projected onto plane encompassed by the hull
    reco::ProjectedPointList fConvexHullPointList; ///< The points on the convex hull
    reco::Hit3DToEdgeMap fConvexHullEdgeMap;       ///< Map from 3D hit to associated edge
    reco::EdgeList fConvexHullEdgeList;            ///< An edge list translated back to 3D hits
    reco::ProjectedPointList
      fConvexHullExtremePoints; ///< The points furthest from each other on hull
    reco::ConvexHullKinkTupleList
      fConvexHullKinkPoints; ///< The points with large kinks along the convex hull
  };

  /**
 *  @brief Class wrapping the above and containing volatile information to characterize the cluster
 */
  class ClusterParameters;
  using ClusterParametersList = std::list<ClusterParameters>;

  class ClusterParameters {
  public:
    ClusterParameters()
    {
      fClusterParams.clear();
      fHitPairListPtr.clear();
      fHit2DToHit3DListMap.clear();
      fHit3DToEdgeMap.clear();
      fBestHitPairListPtr.clear();
      fBestEdgeList.clear();
      fConvexHull.clear();
      fFaceList.clear();
      fVertexList.clear();
      fHalfEdgeList.clear();
      fClusterParameters.clear();
    }

    ClusterParameters(reco::HitPairClusterMap::iterator& mapItr) : fHitPairListPtr(mapItr->second)
    {
      fClusterParams.clear();
      fHit2DToHit3DListMap.clear();
      fHit3DToEdgeMap.clear();
      fBestHitPairListPtr.clear();
      fBestEdgeList.clear();
      fConvexHull.clear();
      fFaceList.clear();
      fVertexList.clear();
      fHalfEdgeList.clear();
    }

    ClusterParameters(reco::HitPairListPtr& hitList) : fHitPairListPtr(hitList)
    {
      fClusterParams.clear();
      fHit2DToHit3DListMap.clear();
      fHit3DToEdgeMap.clear();
      fBestHitPairListPtr.clear();
      fBestEdgeList.clear();
      fConvexHull.clear();
      fFaceList.clear();
      fVertexList.clear();
      fHalfEdgeList.clear();
    }

    ClusterParametersList& daughterList() { return fClusterParameters; }

    void UpdateParameters(const reco::ClusterHit2D* hit)
    {
      fClusterParams[hit->WireID().Plane].UpdateParameters(hit);
    }

    void addHit3D(const reco::ClusterHit3D* hit3D)
    {
      fHitPairListPtr.emplace_back(hit3D);

      for (const auto& hit2D : hit3D->getHits())
        if (hit2D) fHit2DToHit3DListMap[hit2D].emplace_back(hit3D);
    }

    void fillHit2DToHit3DListMap()
    {
      for (const auto& hit3D : fHitPairListPtr) {
        for (const auto& hit2D : hit3D->getHits())
          if (hit2D) fHit2DToHit3DListMap[hit2D].emplace_back(hit3D);
      }
    }

    reco::PlaneToClusterParamsMap& getClusterParams() { return fClusterParams; }
    reco::Hit2DToHit3DListMap& getHit2DToHit3DListMap() { return fHit2DToHit3DListMap; }
    reco::HitPairListPtr& getHitPairListPtr() { return fHitPairListPtr; }
    reco::PrincipalComponents& getFullPCA() { return fFullPCA; }
    reco::PrincipalComponents& getSkeletonPCA() { return fSkeletonPCA; }
    reco::Hit3DToEdgeMap& getHit3DToEdgeMap() { return fHit3DToEdgeMap; }
    reco::HitPairListPtr& getBestHitPairListPtr() { return fBestHitPairListPtr; }
    reco::EdgeList& getBestEdgeList() { return fBestEdgeList; }
    reco::ConvexHull& getConvexHull() { return fConvexHull; }
    dcel2d::FaceList& getFaceList() { return fFaceList; }
    dcel2d::VertexList& getVertexList() { return fVertexList; }
    dcel2d::HalfEdgeList& getHalfEdgeList() { return fHalfEdgeList; }

    friend bool operator<(const ClusterParameters& a, const ClusterParameters& b)
    {
      return a.fHitPairListPtr.size() > b.fHitPairListPtr.size();
    }

  private:
    PlaneToClusterParamsMap fClusterParams;
    reco::HitPairListPtr fHitPairListPtr; // This contains the list of 3D hits in the cluster
    reco::Hit2DToHit3DListMap
      fHit2DToHit3DListMap;             // Provides a mapping between 2D hits and 3D hits they make
    reco::PrincipalComponents fFullPCA; // PCA run over full set of 3D hits
    reco::PrincipalComponents fSkeletonPCA; // PCA run over just the "skeleton" 3D hits
    reco::Hit3DToEdgeMap fHit3DToEdgeMap;
    reco::HitPairListPtr fBestHitPairListPtr;
    reco::EdgeList fBestEdgeList;       // This has become multiuse... really need to split it up
    reco::ConvexHull fConvexHull;       // Convex hull object
    dcel2d::FaceList fFaceList;         // Keeps track of "faces" from Voronoi Diagram
    dcel2d::VertexList fVertexList;     // Keeps track of "vertices" from Voronoi Diagram
    dcel2d::HalfEdgeList fHalfEdgeList; // Keeps track of "halfedges" from Voronoi Diagram
    ClusterParametersList fClusterParameters; // For possible daughter clusters
  };

  using ClusterToHitPairSetPair = std::pair<reco::ClusterParameters*, HitPairSetPtr>;
  using ClusterToHitPairSetMap = std::unordered_map<reco::ClusterParameters*, HitPairSetPtr>;
  using Hit2DToHit3DSetMap = std::unordered_map<const reco::ClusterHit2D*, HitPairSetPtr>;
  using Hit2DToClusterMap = std::unordered_map<const reco::ClusterHit2D*, ClusterToHitPairSetMap>;

}

#endif //RECO_CLUSTER3D_H

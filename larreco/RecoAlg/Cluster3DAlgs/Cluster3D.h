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
#include <vector>
#include <list>
#include <set>
#include <map>
#include <unordered_map>
#include <memory>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/DCEL.h"

// Eigen
#include <Eigen/Dense>

namespace reco {
    
// First define a container object to augment the sparse 2D hit information
class ClusterHit2D
{
public:
    
    ClusterHit2D();   // Default constructor
    
private:
    
    mutable unsigned  m_statusBits;   ///< Volatile status information of this 3D hit
    mutable float     m_docaToAxis;   ///< DOCA of hit at POCA to associated cluster axis
    mutable float     m_arcLenToPoca; ///< arc length to POCA along cluster axis
    float             m_xPosition;    ///< The x coordinate for this hit
    float             m_timeTicks;    ///< The time (in ticks) for this hit
    const recob::Hit* m_hit;          ///< Hit we are augmenting
    
public:
    
    enum StatusBits { SHAREDINPAIR    = 0x00080000,
                      SHAREDINTRIPLET = 0x00040000,
                      USEDINPAIR      = 0x00008000,
                      USEDINTRIPLET   = 0x00004000,
                      SHAREDINCLUSTER = 0x00000200,
                      USEDINCLUSTER   = 0x00000100,
                      USED            = 0x00000001
                    };
    
    ClusterHit2D(unsigned          statusBits,
                 float             doca,
                 float             poca,
                 float             xPosition,
                 float             timeTicks,
                 const recob::Hit& recobHit);
    
    unsigned          getStatusBits()   const {return m_statusBits;}
    float             getDocaToAxis()   const {return m_docaToAxis;}
    float             getArcLenToPoca() const {return m_arcLenToPoca;}
    float             getXPosition()    const {return m_xPosition;}
    float             getTimeTicks()    const {return m_timeTicks;}
    const recob::Hit& getHit()          const {return *m_hit;}
    
    void setStatusBit(unsigned bits)    const {m_statusBits   |= bits;}
    void clearStatusBits(unsigned bits) const {m_statusBits   &= ~bits;}
    void setDocaToAxis(float doca)      const {m_docaToAxis    = doca;}
    void setArcLenToPoca(float poca)    const {m_arcLenToPoca  = poca;}
    
    friend std::ostream& operator << (std::ostream& o, const ClusterHit2D& c);
    friend bool          operator <  (const ClusterHit2D & a, const ClusterHit2D & b);
    
};

using ClusterHit2DVec = std::vector<const reco::ClusterHit2D*>;

// Now define an object with the recob::Hit information that will comprise the 3D cluster
class ClusterHit3D
{
public:
    
    enum StatusBits { REJECTEDHIT     = 0x80000000,           ///< Hit has been rejected for any reason
                      SKELETONHIT     = 0x10000000,           ///< Hit is a "skeleton" hit
                      EDGEHIT         = 0x20000000,           ///< Hit is an "edge" hit
                      SEEDHIT         = 0x40000000,           ///< Hit is part of Seed for track fits
                      MADESPACEPOINT  = 0x08000000,           ///< Hit has been made into Space Point
                      SKELETONPOSAVE  = 0x00100000,           ///< Skeleton hit position averaged
                      CLUSTERVISITED  = 0x00008000,           ///< "visited" by a clustering algorithm
                      CLUSTERNOISE    = 0x00004000,           ///< Labelled "noise" by a clustering algorithm
                      CLUSTERATTACHED = 0x00002000,           ///< attached to a cluster
                      CLUSTERSHARED   = 0x00001000,           ///< 3D hit has 2D hits shared between clusters
                      PATHCHECKED     = 0x00000800,           ///< Path checking algorithm has seen this hit
                      SELECTEDBYMST   = 0x00000100,           ///< Hit has been used in Cluster Splitting MST
                      PCAOUTLIER      = 0x00000080,           ///< Hit labelled outlier in PCA
                      HITINVIEW0      = 0x00000001,           ///< Hit contains 2D hit from view 0 (u plane)
                      HITINVIEW1      = 0x00000002,           ///< Hit contains 2D hit from view 1 (v plane)
                      HITINVIEW2      = 0x00000004            ///< Hit contains 2D hit from view 2 (w plane)
                    };
        
    
    ClusterHit3D();   // Default constructor

    ClusterHit3D(size_t                          id,
                 unsigned int                    statusBits,
                 const float*                    position,
                 float                           totalCharge,
                 float                           avePeakTime,
                 float                           deltaPeakTime,
                 float                           sigmaPeakTime,
                 float                           hitChiSquare,
                 float                           docaToAxis,
                 float                           arclenToPoca,
                 const ClusterHit2DVec&          hitVec,
                 const std::vector<float>&       hitDelTSigVec,
                 const std::vector<geo::WireID>& wireIDVec);
    
    ClusterHit3D(const ClusterHit3D&);

    void initialize(size_t                          id,
                    unsigned int                    statusBits,
                    const float*                    position,
                    float                           totalCharge,
                    float                           avePeakTime,
                    float                           deltaPeakTime,
                    float                           sigmaPeakTime,
                    float                           hitChiSquare,
                    float                           docaToAxis,
                    float                           arclenToPoca,
                    const ClusterHit2DVec&          hitVec,
                    const std::vector<float>&       hitDelTSigVec,
                    const std::vector<geo::WireID>& wireIDVec);

    size_t                              getID()            const {return m_id;}
    unsigned int                        getStatusBits()    const {return m_statusBits;}
    const float*                        getPosition()      const {return m_position;}
    float                               getX()             const {return m_position[0];}
    float                               getY()             const {return m_position[1];}
    float                               getZ()             const {return m_position[2];}
    float                               getTotalCharge()   const {return m_totalCharge;}
    float                               getAvePeakTime()   const {return m_avePeakTime;}
    float                               getDeltaPeakTime() const {return m_deltaPeakTime;}
    float                               getSigmaPeakTime() const {return m_sigmaPeakTime;}
    float                               getHitChiSquare()  const {return m_hitChiSquare;}
    float                               getDocaToAxis()    const {return m_docaToAxis;}
    float                               getArclenToPoca()  const {return m_arclenToPoca;}
    const ClusterHit2DVec&              getHits()          const {return m_hitVector;}
    const std::vector<float>            getHitDelTSigVec() const {return m_hitDelTSigVec;}
    const std::vector<geo::WireID>&     getWireIDs()       const {return m_wireIDVector;}
    
    bool bitsAreSet(const unsigned int& bitsToCheck)       const {return m_statusBits & bitsToCheck;}

    void setID(const size_t& id)           const {m_id            = id;}
    void setStatusBit(unsigned bits)       const {m_statusBits   |= bits;}
    void clearStatusBits(unsigned bits)    const {m_statusBits   &= ~bits;}
    void setDocaToAxis(double doca)        const {m_docaToAxis    = doca;}
    void setArclenToPoca(double poca)      const {m_arclenToPoca  = poca;}
    void setWireID(const geo::WireID& wid) const;
    
    void setPosition(const float* pos) const {m_position[0] = pos[0]; m_position[1] = pos[1]; m_position[2] = pos[2];}

    const bool operator<(const reco::ClusterHit3D& other) const
    {
        if (m_position[2] != other.m_position[2]) return m_position[2] < other.m_position[2];
        else return m_position[0] < other.m_position[0];
    }

    const bool operator==(const reco::ClusterHit3D& other) const
    {
        return m_id == other.m_id;
    }
    
    friend std::ostream& operator << (std::ostream& o, const ClusterHit3D& c);
    //friend bool          operator <  (const ClusterHit3D & a, const ClusterHit3D & b);
    
private:
    
    mutable size_t                   m_id;                 ///< "id" of this hit (useful for indexing)
    mutable unsigned int             m_statusBits;         ///< Volatile status information of this 3D hit
    mutable float                    m_position[3];        ///< position of this hit combination in world coordinates
    float                            m_totalCharge;        ///< Sum of charges of all associated recob::Hits
    float                            m_avePeakTime;        ///< Average peak time of all associated recob::Hits
    float                            m_deltaPeakTime;      ///< Largest delta peak time of associated recob::Hits
    float                            m_sigmaPeakTime;      ///< Quad sum of peak time sigmas
    float                            m_hitChiSquare;       ///< Hit ChiSquare relative to the average time
    mutable float                    m_docaToAxis;         ///< DOCA to the associated cluster axis
    mutable float                    m_arclenToPoca;       ///< arc length along axis to DOCA point
    ClusterHit2DVec                  m_hitVector;          ///< Hits comprising this 3D hit
    mutable std::vector<float>       m_hitDelTSigVec;      ///< Delta t of hit to matching pair / sig
    mutable std::vector<geo::WireID> m_wireIDVector;       ///< Wire ID's for the planes making up hit
};
    
// We also need to define a container for the output of the PCA Analysis
class PrincipalComponents
{
public:

    typedef std::vector<std::vector<float>> EigenVectors;
    
    PrincipalComponents();
    
private:
    
    bool           m_svdOK;              ///< SVD Decomposition was successful
    int            m_numHitsUsed;        ///< Number of hits in the decomposition
    float          m_eigenValues[3];     ///< Eigen values from SVD decomposition
    EigenVectors   m_eigenVectors;       ///< The three principle axes
    float          m_avePosition[3];     ///< Average position of hits fed to PCA
    mutable double m_aveHitDoca;         ///< Average doca of hits used in PCA
    
public:
    
    PrincipalComponents(bool ok, int nHits, const float* eigenValues, const EigenVectors& eigenVecs, const float* avePos, const float aveHitDoca = 9999.);
    
    bool                getSvdOK()            const {return m_svdOK;}
    int                 getNumHitsUsed()      const {return m_numHitsUsed;}
    const float*        getEigenValues()      const {return m_eigenValues;}
    const EigenVectors& getEigenVectors()     const {return m_eigenVectors;}
    const float*        getAvePosition()      const {return m_avePosition;}
    const float         getAveHitDoca()       const {return m_aveHitDoca;}
    
    void                flipAxis(size_t axis);
    void                setAveHitDoca(double doca) const {m_aveHitDoca = doca;}
    
    friend std::ostream&  operator << (std::ostream & o, const PrincipalComponents& a);
    friend bool operator < (const PrincipalComponents& a, const PrincipalComponents& b);
    
};

class Cluster3D
{
public:

    Cluster3D();  ///Default constructor

private:

    mutable unsigned    m_statusBits;       ///< Status bits for the cluster
    PrincipalComponents m_pcaResults;       ///< Output of the prinicipal componenets analysis
    float               m_totalCharge;      ///< Total charge in the cluster
    float               m_startPosition[3]; ///< "start" position for cluster (world coordinates)
    float               m_endPosition[3];   ///< "end" position for cluster
    int                 m_clusterIdx;       ///< ID for this cluster
    
public:
    Cluster3D(unsigned                   statusBits,
              const PrincipalComponents& pcaResults,
              float                      totalCharge,
              const float*               startPosition,
              const float*               endPosition,
              int                        idx);
    
    unsigned                   getStatusBits()    const {return m_statusBits;}
    const PrincipalComponents& getPcaResults()    const {return m_pcaResults;}
    float                      getTotalCharge()   const {return m_totalCharge;}
    const float*               getStartPosition() const {return m_startPosition;}
    const float*               getEndPosition()   const {return m_endPosition;}
    int                        getClusterIdx()    const {return m_clusterIdx;}
    
    void setStatusBit(unsigned bits)    const {m_statusBits |= bits;}
    void clearStatusBits(unsigned bits) const {m_statusBits &= ~bits;}
    
    Cluster3D            operator +  (Cluster3D);
    friend std::ostream& operator << (std::ostream& o, const Cluster3D& c);
    friend bool          operator <  (const Cluster3D & a, const Cluster3D & b);
};

/**
 *  @brief A utility class used in construction of 3D clusters
 */
class RecobClusterParameters
{
public:
    RecobClusterParameters() : m_startTime(999999.),
    m_sigmaStartTime(1.),
    m_endTime(0.),
    m_sigmaEndTime(1.),
    m_totalCharge(0.),
    m_startWire(9999999),
    m_endWire(0),
    m_plane(100),
    m_view(geo::kUnknown)
    {
        m_hitVector.clear();
    }
    
    void UpdateParameters(const reco::ClusterHit2D* hit);
    
    float           m_startTime;
    float           m_sigmaStartTime;
    float           m_endTime;
    float           m_sigmaEndTime;
    float           m_totalCharge;
    unsigned int    m_startWire;
    unsigned int    m_endWire;
    unsigned int    m_plane;
    geo::View_t     m_view;
    ClusterHit2DVec m_hitVector;
};
    
    
/**
 *  @brief export some data structure definitions
 */
using Hit2DListPtr             = std::list<const reco::ClusterHit2D*>;
using HitPairListPtr           = std::list<const reco::ClusterHit3D*>;
using HitPairSetPtr            = std::set<const reco::ClusterHit3D*>;
using HitPairListPtrList       = std::list<HitPairListPtr>;
using HitPairClusterMap        = std::map<int, HitPairListPtr>;
using HitPairList              = std::list<std::unique_ptr<reco::ClusterHit3D>>;
    
using PCAHitPairClusterMapPair = std::pair<reco::PrincipalComponents, reco::HitPairClusterMap::iterator>;
using PlaneToClusterParamsMap  = std::map<size_t, RecobClusterParameters>;
using EdgeTuple                = std::tuple<const reco::ClusterHit3D*,const reco::ClusterHit3D*,double>;
using EdgeList                 = std::list<EdgeTuple>;
using Hit3DToEdgePair          = std::pair<const reco::ClusterHit3D*, reco::EdgeList>;
using Hit3DToEdgeMap           = std::unordered_map<const reco::ClusterHit3D*, reco::EdgeList>;
using Hit2DToHit3DListMap      = std::unordered_map<const reco::ClusterHit2D*, reco::HitPairListPtr>;
//using VertexPoint              = Eigen::Vector3f;
//using VertexPointList          = std::list<Eigen::Vector3f>;

/**
 *  @brief Class wrapping the above and containing volatile information to characterize the cluster
 */
class ClusterParameters;
using ClusterParametersList = std::list<ClusterParameters>;
    
class ClusterParameters
{
public:
    ClusterParameters()
    {
        m_clusterParams.clear();
        m_hitPairListPtr.clear();
        m_hit2DToHit3DListMap.clear();
        m_hit3DToEdgeMap.clear();
        m_bestHitPairListPtr.clear();
        m_bestEdgeList.clear();
        m_faceList.clear();
        m_vertexList.clear();
        m_halfEdgeList.clear();
        m_clusterParameters.clear();
    }
    
    ClusterParameters(reco::HitPairClusterMap::iterator& mapItr) : m_hitPairListPtr(mapItr->second)
    {
        m_clusterParams.clear();
        m_hit2DToHit3DListMap.clear();
        m_hit3DToEdgeMap.clear();
        m_bestHitPairListPtr.clear();
        m_bestEdgeList.clear();
        m_faceList.clear();
        m_vertexList.clear();
        m_halfEdgeList.clear();
    }
    
    ClusterParameters(reco::HitPairListPtr& hitList) : m_hitPairListPtr(hitList)
    {
        m_clusterParams.clear();
        m_hit2DToHit3DListMap.clear();
        m_hit3DToEdgeMap.clear();
        m_bestHitPairListPtr.clear();
        m_bestEdgeList.clear();
        m_faceList.clear();
        m_vertexList.clear();
        m_halfEdgeList.clear();
    }
    
    ClusterParametersList& daughterList() {return m_clusterParameters;}
    
    void UpdateParameters(const reco::ClusterHit2D* hit)
    {
        m_clusterParams[hit->getHit().WireID().Plane].UpdateParameters(hit);
    }
    
    void addHit3D(const reco::ClusterHit3D* hit3D)
    {
        m_hitPairListPtr.emplace_back(hit3D);
        
        for(const auto& hit2D : hit3D->getHits())
            if (hit2D) m_hit2DToHit3DListMap[hit2D].emplace_back(hit3D);
    }
    
    void fillHit2DToHit3DListMap()
    {
        for(const auto& hit3D : m_hitPairListPtr)
        {
            for(const auto& hit2D : hit3D->getHits())
                if (hit2D) m_hit2DToHit3DListMap[hit2D].emplace_back(hit3D);
        }
    }
    
    reco::PlaneToClusterParamsMap& getClusterParams()       {return m_clusterParams;}
    reco::Hit2DToHit3DListMap&     getHit2DToHit3DListMap() {return m_hit2DToHit3DListMap;}
    reco::HitPairListPtr&          getHitPairListPtr()      {return m_hitPairListPtr;}
    reco::PrincipalComponents&     getFullPCA()             {return m_fullPCA;}
    reco::PrincipalComponents&     getSkeletonPCA()         {return m_skeletonPCA;}
    reco::Hit3DToEdgeMap&          getHit3DToEdgeMap()      {return m_hit3DToEdgeMap;}
    reco::HitPairListPtr&          getBestHitPairListPtr()  {return m_bestHitPairListPtr;}
    reco::EdgeList&                getBestEdgeList()        {return m_bestEdgeList;}
    dcel2d::FaceList&              getFaceList()            {return m_faceList;}
    dcel2d::VertexList&            getVertexList()          {return m_vertexList;}
    dcel2d::HalfEdgeList&          getHalfEdgeList()        {return m_halfEdgeList;}

    friend bool operator < (const ClusterParameters &a, const ClusterParameters& b)
    {
        return a.m_hitPairListPtr.size() > b.m_hitPairListPtr.size();
    }

private:
    PlaneToClusterParamsMap   m_clusterParams;
    reco::HitPairListPtr      m_hitPairListPtr;      // This contains the list of 3D hits in the cluster
    reco::Hit2DToHit3DListMap m_hit2DToHit3DListMap; // Provides a mapping between 2D hits and 3D hits they make
    reco::PrincipalComponents m_fullPCA;             // PCA run over full set of 3D hits
    reco::PrincipalComponents m_skeletonPCA;         // PCA run over just the "skeleton" 3D hits
    reco::Hit3DToEdgeMap      m_hit3DToEdgeMap;
    reco::HitPairListPtr      m_bestHitPairListPtr;
    reco::EdgeList            m_bestEdgeList;
    dcel2d::FaceList          m_faceList;            // Keeps track of "faces" from Voronoi Diagram
    dcel2d::VertexList        m_vertexList;          // Keeps track of "vertices" from Voronoi Diagram
    dcel2d::HalfEdgeList      m_halfEdgeList;        // Keeps track of "halfedges" from Voronoi Diagram
    ClusterParametersList     m_clusterParameters;   // For possible daughter clusters
};

using ClusterToHitPairSetPair = std::pair<reco::ClusterParameters*,HitPairSetPtr>;
using ClusterToHitPairSetMap  = std::unordered_map<reco::ClusterParameters*,HitPairSetPtr>;
using Hit2DToHit3DSetMap      = std::unordered_map<const reco::ClusterHit2D*,HitPairSetPtr>;
using Hit2DToClusterMap       = std::unordered_map<const reco::ClusterHit2D*,ClusterToHitPairSetMap>;
    
}

#endif //RECO_CLUSTER3D_H

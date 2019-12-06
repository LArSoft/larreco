/**
 *  @file   DCEL2D.h
 *
 *  @brief  Definitions for a doubly connected edge list
 *          This will define a vertex, half edge and face
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef DCEL2D_h
#define DCEL2D_h

// std includes
#include <vector>
#include <list>
#include <algorithm>

// Eigen
#ifdef __clang__
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Dense>
#ifdef __clang__
#else
#pragma GCC diagnostic pop
#endif

//------------------------------------------------------------------------------------------------------------------------------------------

namespace reco
{
    class ClusterHit3D;
}

namespace dcel2d
{
class HalfEdge;  // Forward declaration

/**
 *  @brief Definitions used by the VoronoiDiagram algorithm
 */
using Point     = std::tuple<double,double,const reco::ClusterHit3D*>;
using PointList = std::list<Point>;
using Coords    = Eigen::Vector3f; //std::pair<double,double>;

class Vertex
{
    /**
     *  @brief  Vertex class definition for use in a doubly connected edge list
     *          a Vertex will contain the coordinates of the point it represents
     *          and a pointer to one of the half edges that emanates from it
     */
public:
    /**
     *  @brief  Constructor
     */
    Vertex() : fCoords(0.,0.,0.), fHalfEdge(NULL) {}

    Vertex(const double* coords, HalfEdge* half) : fHalfEdge(half)
    {
        setCoords(coords);
    }

    Vertex(const Coords& coords, HalfEdge* half) : fCoords(coords), fHalfEdge(half)
    {}

    const Coords&   getCoords()   const {return fCoords;}
    const HalfEdge* getHalfEdge() const {return fHalfEdge;}

    void setCoords(const double* coords)
    {
        fCoords[0] = coords[0];
        fCoords[1] = coords[1];
        fCoords[2] = coords[2];
    }

    void setCoords(const Coords& coords) {fCoords = coords;}

    void setHalfEdge(HalfEdge* half) {fHalfEdge = half;}

private:
    Coords    fCoords;    // x,y coordinates of this vertex
    HalfEdge* fHalfEdge;  // pointer to one of the half edges
};

class Face
{
    /**
     *  @brief  Face class definition for use in a doubly connected edge list
     *          A Face represents the area enclosed by a set of half edges and
     *          vertices. It simply needs to store a pointer to one of the half
     *          edges to be able to recover all
     */
public:
    /**
     *  @brief  Constructor
     */
//    Face() : fHalfEdge(NULL), fPoint(0.,0.,NULL) {}

    Face(HalfEdge* half, const Coords& coords, const reco::ClusterHit3D* clusterHit3D) :
        fHalfEdge(half),
        fConvexHull(false),
        fCoords(coords),
        fFaceArea(0.),
        fClusterHit3D(clusterHit3D)
    {}

    const HalfEdge*           getHalfEdge()     const {return fHalfEdge;}
    const bool                onConvexHull()    const {return fConvexHull;}
    const Coords&             getCoords()       const {return fCoords;}
    const double              getFaceArea()     const {return fFaceArea;}
    const reco::ClusterHit3D* getClusterHit3D() const {return fClusterHit3D;}

    void setHalfEdge(HalfEdge* half)  {fHalfEdge   = half;}
    void setOnConvexHull()            {fConvexHull = true;}
    void setFaceArea(double area)     {fFaceArea   = area;}

private:
    HalfEdge*                 fHalfEdge;     // pointer to one of the half edges
    mutable bool              fConvexHull;   // This face on convex hull
    Coords                    fCoords;       // projected coordinates of associated point
    double                    fFaceArea;     // The area of the face once constructed
    const reco::ClusterHit3D* fClusterHit3D; // The physical 3D hit this corresponds to
};

class HalfEdge
{
    /**
     *  @brief  HalfEdge class definition for use in a doubly connected edge list
     *          The half edge class represents one connection from the vertex it
     *          emanates from and pointing to the next vertex. It is the workhorse
     *          of the doubly connected edge list and contains all necessary
     *          pointers to enable completely traverse of the list.
     */
public:
    /**
     *  @brief  Constructor
     */
    HalfEdge() :
        m_targetVertex(NULL),
        m_face(NULL),
        m_twinHalfEdge(NULL),
        m_nextHalfEdge(NULL),
        m_lastHalfEdge(NULL)
    {}

    HalfEdge(Vertex*   vertex,
             Face*     face,
             HalfEdge* twin,
             HalfEdge* next,
             HalfEdge* last) :
        m_targetVertex(vertex),
        m_face(face),
        m_twinHalfEdge(twin),
        m_nextHalfEdge(next),
        m_lastHalfEdge(last)
    {}

    Vertex*   getTargetVertex() const {return m_targetVertex;}
    Face*     getFace()         const {return m_face;}
    HalfEdge* getTwinHalfEdge() const {return m_twinHalfEdge;}
    HalfEdge* getNextHalfEdge() const {return m_nextHalfEdge;}
    HalfEdge* getLastHalfEdge() const {return m_lastHalfEdge;}

    void setTargetVertex(Vertex* vertex) {m_targetVertex = vertex;}
    void setFace(Face* face)             {m_face         = face;}
    void setTwinHalfEdge(HalfEdge* twin) {m_twinHalfEdge = twin;}
    void setNextHalfEdge(HalfEdge* next) {m_nextHalfEdge = next;}
    void setLastHalfEdge(HalfEdge* last) {m_lastHalfEdge = last;}

private:
    Vertex*     m_targetVertex;  // Pointer to the vertex we point to
    Face*       m_face;          // Pointer to the face we are associated with
    HalfEdge*   m_twinHalfEdge;  // Pointer to the twin half edge (pointing opposite)
    HalfEdge*   m_nextHalfEdge;  // Pointer ot the next half edge
    HalfEdge*   m_lastHalfEdge;  // Pointer to the previous half edge
};

// Define containers to hold the above objects
using VertexList   = std::list<Vertex>;
using FaceList     = std::list<Face>;
using HalfEdgeList = std::list<HalfEdge>;

} // namespace lar_cluster3d
#endif

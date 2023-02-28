/**
 * @file   VoronoiDiagram_test.cxx
 * @brief  Unit test for the Voronoi Diagram code in cluster3d
 * @date   February 15, 2018
 * @author Tracy Usher (usher@slac.stanford.edu)
 *
 * Usage:
 *
 *     TBD....
 *     geometry_icarus_test  [ConfigurationFile [GeometryTestParameterSet]]
 *
 * By default, `GeometryTestParameterSet` is set to `"physics.analyzers.geotest"`.
 *
 */

// ICARUS libraries
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/Voronoi.h"

// LArSoft libraries

// utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

//------------------------------------------------------------------------------
//---  The test environment
//---

//------------------------------------------------------------------------------
//---  The tests
//---

/** ****************************************************************************
 * @brief Runs the test
 * @param argc number of arguments in argv
 * @param argv arguments to the function
 * @return number of detected errors (0 on success)
 * @throw cet::exception most of error situations throw
 *
 * The arguments in argv are:
 * 0. name of the executable ("Geometry_test")
 * 1. path to the FHiCL configuration file
 * 2. FHiCL path to the configuration of the geometry test
 *    (default: physics.analysers.geotest)
 * 3. FHiCL path to the configuration of the geometry
 *    (default: services.Geometry)
 *
 */
//------------------------------------------------------------------------------
int main(int argc, char const** argv)
{
  int nErrors(0);
  /*
    // This test program needs work before being released into the wild...
    // In the meantime, comment out the guts here.

    // Build a test point list
    dcel2d::PointList pointList;

    // Loop through hits and do projection to plane
//    for(const auto& hit3D : clusterParameters.getHitPairListPtr())
//    {
//        Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
//                                    hit3D->getPosition()[1] - pcaCenter(1),
//                                    hit3D->getPosition()[2] - pcaCenter(2));
//        Eigen::Vector3f pcaToHit = rotationMatrix * pcaToHitVec;
//
//        pointList.emplace_back(dcel2d::Point(pcaToHit(0),pcaToHit(1),hit3D));
//    }

    // Make a dummy 3D hit
    reco::ClusterHit3D clusterHit3D;

//    pointList.emplace_back(dcel2d::Point(-10.,  0., &clusterHit3D));
//    pointList.emplace_back(dcel2d::Point(  0., -5., &clusterHit3D));
//    pointList.emplace_back(dcel2d::Point(  0.,  5., &clusterHit3D));
//    pointList.emplace_back(dcel2d::Point( 10.,  0., &clusterHit3D));

    pointList.emplace_back(dcel2d::Point(-10.,  0., &clusterHit3D));
    pointList.emplace_back(dcel2d::Point( -5.,  0., &clusterHit3D));
    pointList.emplace_back(dcel2d::Point(  0.,  0., &clusterHit3D));
    pointList.emplace_back(dcel2d::Point(  5.,  0., &clusterHit3D));
    pointList.emplace_back(dcel2d::Point( 10.,  0., &clusterHit3D));

    // Sort the point vec by decreasing x, then increase y
    pointList.sort([](const auto& left, const auto& right){return (std::abs(std::get<0>(left) - std::get<0>(right)) > std::numeric_limits<float>::epsilon()) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);});

    // Get some useful containers
    dcel2d::FaceList          faceList;            // Keeps track of "faces" from Voronoi Diagram
    dcel2d::VertexList        vertexList;          // Keeps track of "vertices" from Voronoi Diagram
    dcel2d::HalfEdgeList      halfEdgeList;        // Keeps track of "halfedges" from Voronoi Diagram

    // Set up the voronoi diagram builder
    voronoi2d::VoronoiDiagram voronoiDiagram(halfEdgeList,vertexList,faceList);

    // And make the diagram
    voronoiDiagram.buildVoronoiDiagram(pointList);

    // Recover the voronoi diagram vertex list and the container to store them in
//    dcel2d::VertexList& vertexList = clusterParameters.getVertexList();

    // 4. And finally we cross fingers.
    if (nErrors > 0)
    {
        mf::LogError("VoronoiDiagram_test") << nErrors << " errors detected!";
    }
*/
  return nErrors;
} // main()

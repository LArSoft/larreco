/**
 *  @file   ClusterPathFinder_tool.cc
 * 
 *  @brief  art Tool for comparing clusters and merging those that are consistent
 * 
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

#include "larreco/RecoAlg/Cluster3DAlgs/IClusterModAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/ConvexHull/ConvexHull.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Voronoi/Voronoi.h"

// LArSoft includes
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// boost includes
#include <boost/range/adaptor/reversed.hpp>

// Eigen
#include <Eigen/Dense>

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {
    
class ClusterPathFinder : virtual public IClusterModAlg
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit ClusterPathFinder(const fhicl::ParameterSet&);
    
    /**
     *  @brief  Destructor
     */
    ~ClusterPathFinder();
    
    void configure(fhicl::ParameterSet const &pset) override;
    
    /**
     *  @brief Interface for initializing histograms if they are desired
     *         Note that the idea is to put hisgtograms in a subfolder
     *
     *  @param TFileDirectory - the folder to store the hists in
     */
    void initializeHistograms(art::TFileDirectory&) override;

    /**
     *  @brief Scan an input collection of clusters and modify those according
     *         to the specific implementing algorithm
     *
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    void ModifyClusters(reco::ClusterParametersList&) const override;
    
    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    float getTimeToExecute() const override {return m_timeToProcess;}
    
private:
    
    /**
     *  @brief Use PCA to try to find path in cluster
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    reco::ClusterParametersList::iterator breakIntoTinyBits(reco::ClusterParameters&              cluster,
                                                            reco::ClusterParametersList::iterator positionItr,
                                                            reco::ClusterParametersList&          outputClusterList,
                                                            int                                   level = 0) const;

    float closestApproach(const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          const Eigen::Vector3f&,
                          Eigen::Vector3f&,
                          Eigen::Vector3f&) const;
    
    using Point           = std::tuple<float,float,const reco::ClusterHit3D*>;
    using PointList       = std::list<Point>;
    using MinMaxPoints    = std::pair<Point,Point>;
    using MinMaxPointPair = std::pair<MinMaxPoints,MinMaxPoints>;

    void buildConvexHull(reco::ClusterParameters& clusterParameters, int level = 0)     const;
    void buildVoronoiDiagram(reco::ClusterParameters& clusterParameters, int level = 0) const;
    /**
     *  @brief Data members to follow
     */
    bool                                        m_enableMonitoring;      ///<
    size_t                                      m_minTinyClusterSize;    ///< Minimum size for a "tiny" cluster
    mutable float                               m_timeToProcess;         ///<
    
    geo::Geometry*                              m_geometry;              //< pointer to the Geometry service
    
    std::unique_ptr<lar_cluster3d::IClusterAlg> m_clusterAlg;            ///<  Algorithm to do 3D space point clustering
    PrincipalComponentsAlg                      m_pcaAlg;                // For running Principal Components Analysis
};

ClusterPathFinder::ClusterPathFinder(fhicl::ParameterSet const &pset) :
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterPathFinder::~ClusterPathFinder()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterPathFinder::configure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring   = pset.get<bool>  ("EnableMonitoring",  true  );
    m_minTinyClusterSize = pset.get<size_t>("MinTinyClusterSize",40);
    m_clusterAlg         = art::make_tool<lar_cluster3d::IClusterAlg>(pset.get<fhicl::ParameterSet>("ClusterAlg"));

    art::ServiceHandle<geo::Geometry> geometry;
    
    m_geometry = &*geometry;
    
    m_timeToProcess = 0.;
    
    return;
}
    
void ClusterPathFinder::initializeHistograms(art::TFileDirectory&)
{
    return;
}

void ClusterPathFinder::ModifyClusters(reco::ClusterParametersList& clusterParametersList) const
{
    /**
     *  @brief Top level interface for algorithm to consider pairs of clusters from the input
     *         list and determine if they are consistent with each other and, therefore, should
     *         be merged. This is done by looking at the PCA for each cluster and looking at the
     *         projection of the primary axis along the vector connecting their centers.
     */
    
    // Initial clustering is done, now trim the list and get output parameters
    cet::cpu_timer theClockBuildClusters;
    
    // Start clocks if requested
    if (m_enableMonitoring) theClockBuildClusters.start();
    
    int countClusters(0);

    // This is the loop over candidate 3D clusters
    // Note that it might be that the list of candidate clusters is modified by splitting
    // So we use the following construct to make sure we get all of them
    reco::ClusterParametersList::iterator clusterParametersListItr = clusterParametersList.begin();
    
    while(clusterParametersListItr != clusterParametersList.end())
    {
        // Dereference to get the cluster paramters
        reco::ClusterParameters& clusterParameters = *clusterParametersListItr;
        
        std::cout << "**> Looking at Cluster " << countClusters++ << std::endl;
        
        // It turns out that computing the convex hull surrounding the points in the 2D projection onto the
        // plane of largest spread in the PCA is a good way to break up the cluster... and we do it here since
        // we (currently) want this to be part of the standard output
        buildVoronoiDiagram(clusterParameters);

        // Make sure our cluster has enough hits...
        if (clusterParameters.getHitPairListPtr().size() > m_minTinyClusterSize)
        {
            // Get an interim cluster list
            reco::ClusterParametersList reclusteredParameters;
            
            // Call the main workhorse algorithm for building the local version of candidate 3D clusters
            m_clusterAlg->Cluster3DHits(clusterParameters.getHitPairListPtr(), reclusteredParameters);
            
            std::cout << ">>>>>>>>>>> Reclustered to " << reclusteredParameters.size() << " Clusters <<<<<<<<<<<<<<<" << std::endl;
            
            // Only process non-empty results
            if (!reclusteredParameters.empty())
            {
                // Loop over the reclustered set
                for (auto& cluster : reclusteredParameters)
                {
                    std::cout << "****> Calling breakIntoTinyBits" << std::endl;
                    
                    // Break our cluster into smaller elements...
                    breakIntoTinyBits(cluster, cluster.daughterList().end(), cluster.daughterList(), 4);
        
                    std::cout << "****> Broke Cluster with " << cluster.getHitPairListPtr().size() << " into " << cluster.daughterList().size() << " sub clusters";
                    for(auto& clus : cluster.daughterList()) std::cout << ", " << clus.getHitPairListPtr().size();
                    std::cout << std::endl;
            
                    // Add the daughters to the cluster
                    clusterParameters.daughterList().insert(clusterParameters.daughterList().end(),cluster); 
                }
            }
        }
    
        // Go to next cluster parameters object
        clusterParametersListItr++;
    }

    if (m_enableMonitoring)
    {
        theClockBuildClusters.stop();
        
        m_timeToProcess = theClockBuildClusters.accumulated_real_time();
    }
    
    mf::LogDebug("Cluster3D") << ">>>>> Cluster Path finding done" << std::endl;
    
    return;
}
    
reco::ClusterParametersList::iterator ClusterPathFinder::breakIntoTinyBits(reco::ClusterParameters&              clusterToBreak,
                                                                           reco::ClusterParametersList::iterator positionItr,
                                                                           reco::ClusterParametersList&          outputClusterList,
                                                                           int                                   level) const
{
    // This needs to be a recursive routine...
    // Idea is to take the input cluster and order 3D hits by arclength along PCA primary axis
    // If the cluster is above the minimum number of hits then we divide into two and call ourself
    // with the two halves. This means we form a new cluster with hits and PCA and then call ourself
    // If the cluster is below the minimum then we can't break any more, simply add this cluster to
    // the new output list.
    
    // set an indention string
    std::string pluses(level/2, '+');
    std::string indent(level/2, ' ');
    
    indent += pluses;
    
    reco::ClusterParametersList::iterator inputPositionItr = positionItr;
    
    // To make best use of this we'll also want the PCA for this cluster... so...
    // Recover the prime ingredients
    reco::PrincipalComponents& fullPCA     = clusterToBreak.getFullPCA();
    
    std::cout << indent << ">>> breakIntoTinyBits with " << clusterToBreak.getHitPairListPtr().size() << " input hits " << std::endl;
    
    // It turns out that computing the convex hull surrounding the points in the 2D projection onto the
    // plane of largest spread in the PCA is a good way to break up the cluster... and we do it here since
    // we (currently) want this to be part of the standard output
    buildConvexHull(clusterToBreak, level+2);

    bool storeCurrentCluster(true);
    int  minimumClusterSize(m_minTinyClusterSize);
    
    // Create a rough cut intended to tell us when we've reached the land of diminishing returns
    if (clusterToBreak.getBestEdgeList().size()   > 6 &&
        clusterToBreak.getHitPairListPtr().size() > size_t(2 * minimumClusterSize))
    {
        // Recover the list of 3D hits associated to this cluster
        reco::HitPairListPtr& clusHitPairVector = clusterToBreak.getHitPairListPtr();
        
        // Calculate the doca to the PCA primary axis for each 3D hit
        // Importantly, this also gives us the arclength along that axis to the hit
        m_pcaAlg.PCAAnalysis_calc3DDocas(clusHitPairVector, fullPCA);
        
        // Sort the hits along the PCA
        clusHitPairVector.sort([](const auto& left, const auto& right){return left->getArclenToPoca() < right->getArclenToPoca();});
        
        // Now we use the convex hull vertex points to form split points for breaking up the incoming cluster
        reco::EdgeList&                        bestEdgeList = clusterToBreak.getBestEdgeList();
        std::vector<const reco::ClusterHit3D*> vertexHitVec;
        
        std::cout << indent << "+> Breaking cluster, convex hull has " << bestEdgeList.size() << " edges to work with" << std::endl;
        
        for(const auto& edge : bestEdgeList)
        {
            vertexHitVec.push_back(std::get<0>(edge));
            vertexHitVec.push_back(std::get<1>(edge));
        }
        
        // Sort this vector, we aren't worried about duplicates right now...
        std::sort(vertexHitVec.begin(),vertexHitVec.end(),[](const auto& left, const auto& right){return left->getArclenToPoca() < right->getArclenToPoca();});
        
        // Now we create a list of pairs of iterators to the start and end of each subcluster
        using Hit3DItrPair   = std::pair<reco::HitPairListPtr::iterator,reco::HitPairListPtr::iterator>;
        using VertexPairList = std::list<Hit3DItrPair>;
        
        VertexPairList vertexPairList;
        reco::HitPairListPtr::iterator firstHitItr = clusHitPairVector.begin();
        
        for(const auto& hit3D : vertexHitVec)
        {
            reco::HitPairListPtr::iterator vertexItr = std::find(firstHitItr,clusHitPairVector.end(),hit3D);
            
            if (vertexItr == clusHitPairVector.end())
            {
                std::cout << indent << ">>>>>>>>>>>>>>>>> Hit not found in input list, cannot happen? <<<<<<<<<<<<<<<<<<<"  << std::endl;
                break;
            }
            
            std::cout << indent << "+> -- Distance from first to current vertex point: " << std::distance(firstHitItr,vertexItr) << " first: " << *firstHitItr << ", vertex: " << *vertexItr;
            
            // Require a minimum number of points...
            if (std::distance(firstHitItr,vertexItr) > minimumClusterSize)
            {
                vertexPairList.emplace_back(Hit3DItrPair(firstHitItr,vertexItr));
                firstHitItr = vertexItr;
                
                std::cout << " ++ made pair ";
            }
            
            std::cout << std::endl;
        }
        
        // Not done if there is distance from first to end of list
        if (std::distance(firstHitItr,clusHitPairVector.end()) > 0)
        {
            std::cout << indent << "+> loop over vertices done, remant distance: " << std::distance(firstHitItr,clusHitPairVector.end()) << std::endl;
            
            // In the event we don't have the minimum number of hits we simply extend the last pair
            if (!vertexPairList.empty() && std::distance(firstHitItr,clusHitPairVector.end()) < minimumClusterSize)
                vertexPairList.back().second = clusHitPairVector.end();
            else
                vertexPairList.emplace_back(Hit3DItrPair(firstHitItr,clusHitPairVector.end()));
        }
        
        std::cout << indent << "+> ---> breaking cluster into " << vertexPairList.size() << " subclusters" << std::endl;
        
        if (vertexPairList.size() > 1)
        {
            storeCurrentCluster = false;
            
            // Ok, now loop through our pairs
            for(auto& hit3DItrPair : vertexPairList)
            {
                reco::ClusterParameters clusterParams;
                reco::HitPairListPtr&   hitPairListPtr = clusterParams.getHitPairListPtr();
            
                std::cout << indent << "+>    -- building new cluster, size: " << std::distance(hit3DItrPair.first,hit3DItrPair.second) << std::endl;

                // size the container...
                hitPairListPtr.resize(std::distance(hit3DItrPair.first,hit3DItrPair.second));

                // and copy the hits into the container
                std::copy(hit3DItrPair.first,hit3DItrPair.second,hitPairListPtr.begin());
            
                // First stage of feature extraction runs here
                m_pcaAlg.PCAAnalysis_3D(hitPairListPtr, clusterParams.getFullPCA());
                
                // Recover the new fullPCA
                reco::PrincipalComponents& newFullPCA = clusterParams.getFullPCA();

                // Must have a valid pca
                if (newFullPCA.getSvdOK())
                {
                    std::cout << indent << "+>    -- >> cluster has a valid Full PCA" << std::endl;

                    // Need to check if the PCA direction has been reversed
                    Eigen::Vector3f fullPrimaryVec(fullPCA.getEigenVectors().row(0));
                    Eigen::Vector3f newPrimaryVec(newFullPCA.getEigenVectors().row(0));
                    
                    // If the PCA's are opposite the flip the axes
                    if (fullPrimaryVec.dot(newPrimaryVec) < 0.)
                    {
                        for(size_t vecIdx = 0; vecIdx < 3; vecIdx++) newFullPCA.flipAxis(vecIdx);
                    }

                    // Set the skeleton PCA to make sure it has some value
                    clusterParams.getSkeletonPCA() = clusterParams.getFullPCA();
                
                    positionItr = breakIntoTinyBits(clusterParams, positionItr, outputClusterList, level+4);
                }
            }
        }
    }

    // First question, are we done breaking?
    if (storeCurrentCluster)
    {
        // I think this is where we fill out the rest of the parameters?
        // Start by adding the 2D hits...
        // See if we can avoid duplicates by temporarily transferring to a set
        std::set<const reco::ClusterHit2D*> hitSet;
        
        // Loop through 3D hits to get a set of unique 2D hits
        for(const auto& hit3D : clusterToBreak.getHitPairListPtr())
        {
            for(const auto& hit2D : hit3D->getHits())
            {
                if (hit2D) hitSet.insert(hit2D);
            }
        }
        
        // Now add these to the new cluster
        for(const auto& hit2D : hitSet)
        {
            hit2D->setStatusBit(reco::ClusterHit2D::USED);
            clusterToBreak.UpdateParameters(hit2D);
        }
        
        std::cout << indent << "*********>>> storing new subcluster of size " << clusterToBreak.getHitPairListPtr().size() << std::endl;
        
        positionItr = outputClusterList.insert(positionItr,clusterToBreak);
        
        // The above points to the element, want the next element
        positionItr++;
    }
    else if (inputPositionItr == positionItr) std::cout << indent << "***** DID NOT STORE A CLUSTER *****" << std::endl;
    
    return positionItr;
}
    
void ClusterPathFinder::buildConvexHull(reco::ClusterParameters& clusterParameters, int level) const
{
    // set an indention string
    std::string minuses(level/2, '-');
    std::string indent(level/2, ' ');
    
    indent += minuses;

    // The plan is to build the enclosing 2D polygon around the points in the PCA plane of most spread for this cluster
    // To do so we need to start by building a list of 2D projections onto the plane of most spread...
    reco::PrincipalComponents& pca = clusterParameters.getFullPCA();

    // Recover the parameters from the Principal Components Analysis that we need to project and accumulate
    const Eigen::Vector3f& pcaCenter = pca.getAvePosition();

    //dcel2d::PointList pointList;
    using Point     = std::tuple<float,float,const reco::ClusterHit3D*>;
    using PointList = std::list<Point>;
    PointList pointList;

    // Loop through hits and do projection to plane
    for(const auto& hit3D : clusterParameters.getHitPairListPtr())
    {
        Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
                                    hit3D->getPosition()[1] - pcaCenter(1),
                                    hit3D->getPosition()[2] - pcaCenter(2));
        Eigen::Vector3f pcaToHit = pca.getEigenVectors() * pcaToHitVec;

        pointList.emplace_back(dcel2d::Point(pcaToHit(1),pcaToHit(2),hit3D));
    }
    
    // Sort the point vec by increasing x, then increase y
    pointList.sort([](const auto& left, const auto& right){return (std::abs(std::get<0>(left) - std::get<0>(right)) > std::numeric_limits<float>::epsilon()) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);});
    
    // containers for finding the "best" hull...
    std::vector<ConvexHull> convexHullVec;
    std::vector<PointList>  rejectedListVec;
    bool                    increaseDepth(pointList.size() > 5);
    float                   lastArea(std::numeric_limits<float>::max());
    
    while(increaseDepth)
    {
        // Get another convexHull container
        convexHullVec.push_back(ConvexHull(pointList));
        rejectedListVec.push_back(PointList());

        const ConvexHull& convexHull       = convexHullVec.back();
        PointList&        rejectedList     = rejectedListVec.back();
        const PointList&  convexHullPoints = convexHull.getConvexHull();
        
        increaseDepth = false;
        
        if (convexHull.getConvexHullArea() > 0.)
        {
            std::cout << indent << "-> built convex hull, 3D hits: " << pointList.size() << " with " << convexHullPoints.size() << " vertices" << ", area: " << convexHull.getConvexHullArea() << std::endl;
            std::cout << indent << "-> -Points:";
            for(const auto& point : convexHullPoints)
                std::cout << " (" << std::get<0>(point) << "," << std::get<1>(point) << ")";
            std::cout << std::endl;
            
            if (convexHullVec.size() < 2 || convexHull.getConvexHullArea() < 0.8 * lastArea)
            {
                for(auto& point : convexHullPoints)
                {
                    pointList.remove(point);
                    rejectedList.emplace_back(point);
                }
                lastArea      = convexHull.getConvexHullArea();
//                increaseDepth = true;
            }
        }
    }
    
    // do we have a valid convexHull?
    while(!convexHullVec.empty() && convexHullVec.back().getConvexHullArea() < 0.5)
    {
        convexHullVec.pop_back();
        rejectedListVec.pop_back();
    }

    // If we found the convex hull then build edges around the region
    if (!convexHullVec.empty())
    {
        size_t               nRejectedTotal(0);
        reco::HitPairListPtr hitPairListPtr = clusterParameters.getHitPairListPtr();
        
        for(const auto& rejectedList : rejectedListVec)
        {
            nRejectedTotal += rejectedList.size();
            
            for(const auto& rejectedPoint : rejectedList)
            {
                std::cout << indent << "-> -- Point is " << convexHullVec.back().findNearestDistance(rejectedPoint) << " from nearest edge" << std::endl;
                
                if (convexHullVec.back().findNearestDistance(rejectedPoint) > 0.5)
                    hitPairListPtr.remove(std::get<2>(rejectedPoint));
            }
        }
        
        std::cout << indent << "-> Removed " << nRejectedTotal << " leaving " << pointList.size() << "/" << hitPairListPtr.size() << " points" << std::endl;
        
        // Now add "edges" to the cluster to describe the convex hull (for the display)
        reco::Hit3DToEdgeMap& edgeMap      = clusterParameters.getHit3DToEdgeMap();
        reco::EdgeList&       bestEdgeList = clusterParameters.getBestEdgeList();

        Point lastPoint = convexHullVec.back().getConvexHull().front();
    
        for(auto& curPoint : convexHullVec.back().getConvexHull())
        {
            if (curPoint == lastPoint) continue;
        
            const reco::ClusterHit3D* lastPoint3D = std::get<2>(lastPoint);
            const reco::ClusterHit3D* curPoint3D  = std::get<2>(curPoint);
        
            float distBetweenPoints = (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) * (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0])
                                    + (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) * (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1])
                                    + (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]) * (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]);
        
            distBetweenPoints = std::sqrt(distBetweenPoints);
        
            reco::EdgeTuple edge(lastPoint3D,curPoint3D,distBetweenPoints);
        
            edgeMap[lastPoint3D].push_back(edge);
            edgeMap[curPoint3D].push_back(edge);
            bestEdgeList.emplace_back(edge);
        
            lastPoint = curPoint;
        }
    }

    return;
}

    
void ClusterPathFinder::buildVoronoiDiagram(reco::ClusterParameters& clusterParameters, int level) const
{
    // The plan is to build the enclosing 2D polygon around the points in the PCA plane of most spread for this cluster
    // To do so we need to start by building a list of 2D projections onto the plane of most spread...
    reco::PrincipalComponents& pca = clusterParameters.getFullPCA();
    
    // Recover the parameters from the Principal Components Analysis that we need to project and accumulate
    const Eigen::Vector3f& pcaCenter = pca.getAvePosition();
    dcel2d::PointList      pointList;
    
    // Loop through hits and do projection to plane
    for(const auto& hit3D : clusterParameters.getHitPairListPtr())
    {
        Eigen::Vector3f pcaToHitVec(hit3D->getPosition()[0] - pcaCenter(0),
                                    hit3D->getPosition()[1] - pcaCenter(1),
                                    hit3D->getPosition()[2] - pcaCenter(2));
        Eigen::Vector3f pcaToHit = pca.getEigenVectors() * pcaToHitVec;
        
        pointList.emplace_back(dcel2d::Point(pcaToHit(1),pcaToHit(2),hit3D));
    }
    
    // Sort the point vec by increasing x, then increase y
    pointList.sort([](const auto& left, const auto& right){return (std::abs(std::get<0>(left) - std::get<0>(right)) > std::numeric_limits<float>::epsilon()) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);});
    
    // Set up the voronoi diagram builder
    voronoi2d::VoronoiDiagram voronoiDiagram(clusterParameters.getHalfEdgeList(),clusterParameters.getVertexList(),clusterParameters.getFaceList());
    
    // And make the diagram
    voronoiDiagram.buildVoronoiDiagram(pointList);

    // Recover the voronoi diagram vertex list and the container to store them in
    dcel2d::VertexList& vertexList = clusterParameters.getVertexList();
    
    // Now get the inverse of the rotation matrix so we can get the vertex positions,
    // which lie in the plane of the two largest PCA axes, in the standard coordinate system
    Eigen::Matrix3f rotationMatrixInv = pca.getEigenVectors().inverse();
    
    // Translate and fill
    for(auto& vertex : vertexList)
    {
        Eigen::Vector3f coords = rotationMatrixInv * vertex.getCoords();
        
        coords += pcaCenter;
        
        vertex.setCoords(coords);
    }
    
    // Now do the Convex Hull
    // Now add "edges" to the cluster to describe the convex hull (for the display)
    reco::Hit3DToEdgeMap& edgeMap      = clusterParameters.getHit3DToEdgeMap();
    reco::EdgeList&       bestEdgeList = clusterParameters.getBestEdgeList();
    
//    const dcel2d::PointList& edgePoints = voronoiDiagram.getConvexHull();
//    PointList                localList;
//
//    for(const auto& edgePoint : edgePoints) localList.emplace_back(std::get<0>(edgePoint),std::get<1>(edgePoint),std::get<2>(edgePoint));
//
//    // Sort the point vec by increasing x, then increase y
//    localList.sort([](const auto& left, const auto& right){return (std::abs(std::get<0>(left) - std::get<0>(right)) > std::numeric_limits<float>::epsilon()) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);});
//
//    // Why do I need to do this?
//    ConvexHull convexHull(localList);
//
//    Point lastPoint = convexHull.getConvexHull().front();
    dcel2d::Point lastPoint = voronoiDiagram.getConvexHull().front();

//    std::cout << "@@@@>> Build convex hull, voronoi handed " << edgePoints.size() << " points, convexHull cut to " << convexHull.getConvexHull().size() << std::endl;
    
//    for(auto& curPoint : convexHull.getConvexHull())
    for(auto& curPoint : voronoiDiagram.getConvexHull())
    {
        if (curPoint == lastPoint) continue;
        
        const reco::ClusterHit3D* lastPoint3D = std::get<2>(lastPoint);
        const reco::ClusterHit3D* curPoint3D  = std::get<2>(curPoint);
        
        float distBetweenPoints = (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0]) * (curPoint3D->getPosition()[0] - lastPoint3D->getPosition()[0])
                                + (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1]) * (curPoint3D->getPosition()[1] - lastPoint3D->getPosition()[1])
                                + (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]) * (curPoint3D->getPosition()[2] - lastPoint3D->getPosition()[2]);
        
        distBetweenPoints = std::sqrt(distBetweenPoints);
        
        reco::EdgeTuple edge(lastPoint3D,curPoint3D,distBetweenPoints);
        
        edgeMap[lastPoint3D].push_back(edge);
        edgeMap[curPoint3D].push_back(edge);
        bestEdgeList.emplace_back(edge);
        
        lastPoint = curPoint;
    }

    std::cout << "****> vertexList containted " << vertexList.size() << " vertices" << std::endl;
    
    return;
}

    
float ClusterPathFinder::closestApproach(const Eigen::Vector3f& P0,
                                         const Eigen::Vector3f& u0,
                                         const Eigen::Vector3f& P1,
                                         const Eigen::Vector3f& u1,
                                         Eigen::Vector3f&       poca0,
                                         Eigen::Vector3f&       poca1) const
{
    // Technique is to compute the arclength to each point of closest approach
    Eigen::Vector3f w0 = P0 - P1;
    float a(1.);
    float b(u0.dot(u1));
    float c(1.);
    float d(u0.dot(w0));
    float e(u1.dot(w0));
    float den(a * c - b * b);
    
    float arcLen0 = (b * e - c * d) / den;
    float arcLen1 = (a * e - b * d) / den;
    
    poca0 = P0 + arcLen0 * u0;
    poca1 = P1 + arcLen1 * u1;
    
    return (poca0 - poca1).norm();
}
    
DEFINE_ART_CLASS_TOOL(ClusterPathFinder)
} // namespace lar_cluster3d

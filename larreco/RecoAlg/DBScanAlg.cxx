////////////////////////////////////////////////////////////////////////
//
// \file DBScanAlg.cxx
//
// kinga.partyka@yale.edu
//
//  This algorithm finds clusters of hits, they can be of arbitrary shape.You need to specify 2(3) parameters: 
// epsilon, epsilon2 and MinPoints as explained in the corresponding xml file.In my comments a 'point' reference 
// appears quite often. A 'point' is basically a simple hit which only contains wire and time information. This 
// algorithm is based on DBSCAN(Density Based Spatial Clustering of Applications with Noise): M. Ester, H.-P. Kriegel, 
// J. Sander, and X. Xu, A density-based algorithm for discovering clusters in large spatial databases with noise, 
// Second International Conference on Knowledge Discovery and Data Mining, pp. 226-231, AAAI Press. 1996. 
// ( Some of this code is from "Antonio Gulli's coding playground")  
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/LArProperties.h"
#include "lardata/Utilities/DetectorProperties.h"
#include "larreco/RecoAlg/DBScanAlg.h"
#include "lardata/RecoBase/Hit.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TH1.h"

//----------------------------------------------------------
// RStarTree stuff
//----------------------------------------------------------
class dbsPoint{
public:
  double x, y;
  double dx, dy;
  dbsPoint(double X=0.0, double Y=0.0, double dX=0.0, double dY=0.0)
    :x(X), y(Y), dx(dX), dy(dY){};
  BoundingBox bounds()const;
  void Expand(double DX, double DY){dx+=DX; dy+=DY;};
};

BoundingBox dbsPoint::bounds()const{
  BoundingBox bb;
  bb.edges[0].first  = x - std::abs(dx);
  bb.edges[0].second = x + std::abs(dx);
  
  bb.edges[1].first  = y - std::abs(dy);
  bb.edges[1].second = y + std::abs(dy);
  return bb;
}

//----------------------------------------------------------
// Set Visitor
//
// collects accepted leafs in the std::set sResult and the std::vector
// vResult
struct Visitor {
  unsigned int count;
  std::vector< unsigned int > vResult;
  std::set< unsigned int > sResult;
  const bool ContinueVisiting;
  Visitor() : count(0), vResult(), sResult(), ContinueVisiting(true) {};
  void operator()(const RTree::Leaf * const leaf){
    vResult.push_back(leaf->leaf);
    sResult.insert(leaf->leaf);
    count++;
  }
};

//----------------------------------------------------------
// Ellipse acceptor
//
// Roughly quivalent to what FindNeighbors was doing before, except
// that it doesn't handle dead wires
struct AcceptEllipse {
  const BoundingBox &m_bound;
  double r[2];
  double c[2]; ///< center of the bounding box
  double d[2]; ///< half 
  explicit AcceptEllipse(const BoundingBox&b, 
			 double r1, double r2):m_bound(b),r(),c(){
    r[0] = r1;
    r[1] = r2;
    c[0] = (m_bound.edges[0].second+m_bound.edges[0].first)/2.0;
    c[1] = (m_bound.edges[1].second+m_bound.edges[1].first)/2.0;
    d[0] = (m_bound.edges[0].second-m_bound.edges[0].first)/2.0;
    d[1] = (m_bound.edges[1].second-m_bound.edges[1].first)/2.0;
  }
  bool operator()(const RTree::Node * const node) const {
    // At the node level we use a rectangualr overlap condition
    return m_bound.overlaps(node->bound);
  }
  bool operator()(const RTree::Leaf * const leaf) const {
    // At the leaf level we have to more careful
    double C[2], D[2];
    C[0] = (leaf->bound.edges[0].second + leaf->bound.edges[0].first)/2.0;
    C[1] = (leaf->bound.edges[1].second + leaf->bound.edges[1].first)/2.0;
    D[0] = (leaf->bound.edges[0].second - leaf->bound.edges[0].first)/2.0;
    D[1] = (leaf->bound.edges[1].second - leaf->bound.edges[1].first)/2.0;
    double t=0;
    for (int i=0; i<2; ++i){
      // This is only approximate, it will accept a few classes of
      // non-overlapping ellipses
      t += ((c[i]-C[i])*(c[i]-C[i]))/((d[i]+D[i])*(d[i]+D[i]));
    }
    return (t < 1);
  }
private:
  static const BoundingBox EmptyBoundingBox; // for uninitialized bounds
  
  AcceptEllipse()
    :m_bound(EmptyBoundingBox),r(),c(){
    r[0] = r[1] = 1.0;
    c[0] = (m_bound.edges[0].second-m_bound.edges[0].first)/2.0;
    c[1] = (m_bound.edges[1].second-m_bound.edges[1].first)/2.0;
  }
};

const BoundingBox AcceptEllipse::EmptyBoundingBox;

//----------------------------------------------------------
// FindNeighbors acceptor
//
// Exactly equivalent to what FindNeighbors was doing before (assuming
// that points are preresented with no width in wire-space and with
// width in time-space
//
// We're going to make this work with nodal bounding boxes by
// accepting on center-inside the box or comparing to the nearest
// point on the edge using the point-to-point comparison function and
// a the maximum time-width.
struct AcceptFindNeighbors {
  const BoundingBox &fBound;
  double fEps[2];
  double fMaxWidth;
  double fWireDist;
  std::vector< unsigned int > &fBadWireSum;
  AcceptFindNeighbors(const BoundingBox&b, double eps, double eps2, 
		      double maxWidth, double wireDist,
		      std::vector< unsigned int > &badWireSum)
    :fBound(b), fEps(), fMaxWidth(maxWidth), fWireDist(wireDist)
    ,fBadWireSum(badWireSum)
  {
    fEps[0] = eps;
    fEps[1] = eps2;
  }
  // return a point-like box at the center of b
  BoundingBox center(const BoundingBox&b) const {
    BoundingBox c;
    c.edges[0].first = c.edges[0].second = ( b.edges[0].first +
					     b.edges[0].second )/2.0;
    c.edges[1].first = c.edges[1].second = ( b.edges[1].first +
					     b.edges[1].second )/2.0;
    return c;
  }
  BoundingBox center() const { return center(fBound); }
  // Here we implement the findNeighbor algorithm from point to point
  bool isNear(const BoundingBox &b) const { 
    // Precomupation of a few things...box centers, wire bridging
    // quantities, etc...
    double bCenter0 = center(b).edges[0].first;
    double bCenter1 = center(b).edges[1].first; 
    double tCenter0 = center().edges[0].first; // "t" is for test-point
    double tCenter1 = center().edges[1].first;
    // widths in the time direction
    double bWidth = std::abs(     b.edges[1].second -      b.edges[1].first);
    double tWidth = std::abs(fBound.edges[1].second - fBound.edges[1].first);
    // bad channel counting
    unsigned int wire1 = (unsigned int) (tCenter0/fWireDist + 0.5);
    unsigned int wire2 = (unsigned int) (bCenter0/fWireDist + 0.5);
    // Clamp the wize number to something resonably.
    ///\todo activating these should throw a warning or something
    if (wire1 < fBadWireSum.size()) wire1 = fBadWireSum.size();
    if (wire2 < fBadWireSum.size()) wire2 = fBadWireSum.size();
    // The getSimilarity[2] wirestobridge calculation is asymmetric,
    // but is plugged into the cache symmetrically.I am assuming that
    // this is OK because the wires that are hit cannot be bad.
    unsigned int wirestobridge = abs(fBadWireSum[wire1] - fBadWireSum[wire2]);
    double cmtobridge = wirestobridge*fWireDist;

    // getSimilarity()
    double sim = std::abs(tCenter0 - bCenter0) - cmtobridge;
    sim*=sim; // square it

    // getSimilarity2()
    if ( std::abs(tCenter0 - bCenter0) > 1e-10 ) { 
      cmtobridge *= std::abs((tCenter1-bCenter1)/(tCenter0-bCenter0));
    }
    double sim2 = std::abs(tCenter1 - bCenter1) - cmtobridge;
    sim2 *= sim2; // square it

    // getWidthFactor() 
    //     double k=0.13;  // see the comments on getWidthFactor
    //     double k=0.78;  // as I don't know where these came from
    //     double k=4.5;
    //     double k=1.96;
    double k=0.1;
    double WFactor = (exp(4.6*((tWidth*tWidth)+(bWidth*bWidth))))*k; // ??
    // We clamp WFactor on [ 1.0, 6.25 ]
    if (WFactor<1.0)  WFactor = 1.0;
    if (WFactor>6.25) WFactor = 6.25;

    // Now we implement the test...see FindNeighbors
    return (((sim )/(fEps[0]*fEps[0])          ) + 
	    ((sim2)/(fEps[1]*fEps[1]*(WFactor))) <= 1 );
  }
  BoundingBox nearestPoint(const BoundingBox&b) const {
    BoundingBox n;
    BoundingBox c = center();
    for (int i=0; i<2; ++i) {
      // The work for finding the nearest point is the same in both
      // dimensions
      if        ( b.edges[i].first  > c.edges[i].second  ) {
	// Our point is lower than the low edge of the box
	n.edges[i].first = n.edges[i].second = b.edges[i].first;
      } else if ( b.edges[0].second < c.edges[0].first ) {
	// Our point is higher than the high edge of the box
	n.edges[i].first = n.edges[i].second = b.edges[i].second;
      } else {
	// In this dimension our point lies within the boxes bounds
	n.edges[i].first = n.edges[i].second = c.edges[i].first;
      }
    }
    // Now give the time dimension a width
    n.edges[1].first  -= fMaxWidth/2.0;
    n.edges[1].second += fMaxWidth/2.0;
    return n;
  }
  bool operator()(const RTree::Node * const node) const {
    // if the our point overlaps the bounding box, accept immediately
    if (fBound.overlaps(node->bound)) return true;
    // No overlap, so compare to the nearest point on the bounding box
    // under the assumption that the maximum width applies for that
    // point
    return isNear(nearestPoint(node->bound));
  }
  bool operator()(const RTree::Leaf * const leaf) const {
    return isNear(leaf->bound);
  }
};

namespace cluster{
  const unsigned int kNO_CLUSTER    = UINT_MAX;
  const unsigned int kNOISE_CLUSTER = UINT_MAX-1;
}

//----------------------------------------------------------
// DBScanAlg stuff
//----------------------------------------------------------
cluster::DBScanAlg::DBScanAlg(fhicl::ParameterSet const& pset)
{
 this->reconfigure(pset); 
}

//----------------------------------------------------------
cluster::DBScanAlg::~DBScanAlg()
{
}

//----------------------------------------------------------
void cluster::DBScanAlg::reconfigure(fhicl::ParameterSet const& p)
{
  fEps            = p.get< double >("eps"   );
  fEps2           = p.get< double >("epstwo");
  fMinPts         = p.get< int    >("minPts");
  fClusterMethod  = p.get< int    >("Method");
  fDistanceMetric = p.get< int    >("Metric");
}

//----------------------------------------------------------
void cluster::DBScanAlg::InitScan(const std::vector< art::Ptr<recob::Hit> >& allhits, 
				  std::set<uint32_t>                   badChannels,
				  const std::vector<geo::WireID> & wireids)
{
  if (wireids.size()&&wireids.size()!=allhits.size()){
    throw cet::exception("DBScanAlg") << "allhits size = "<<allhits.size()<<" wireids size = "<<wireids.size()<<" do not match\n";
  }
  // clear all the data member vectors for the new set of hits
  fps.clear();
  fpointId_to_clusterId.clear();
  fnoise.clear();
  fvisited.clear();
  fsim.clear();
  fsim2.clear();
  fsim3.clear();
  fclusters.clear();
  fWirePitch.clear();

  fBadChannels = badChannels;
  fBadWireSum.clear();

  // Clear the RTree
  fRTree.Remove(RTree::AcceptAny(),RTree::RemoveLeaf());
  // and the bounds list
  fRect.clear();

  //------------------------------------------------------------------
  // Determine spacing between wires (different for each detector)
  ///get 2 first wires and find their spacing (wire_dist)

  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detp;
  art::ServiceHandle<geo::Geometry> geom;

  for(size_t p = 0; p < geom->Nplanes(); ++p)
    fWirePitch.push_back(geom->WirePitch(0,1,p));

  
  // Collect the bad wire list into a useful form
  if (fClusterMethod) { // Using the R*-tree
    fBadWireSum.resize(geom->Nchannels());
    unsigned int count=0;
    for (unsigned int i=0; i<fBadWireSum.size(); ++i) {
      count += fBadChannels.count(i);
      fBadWireSum[i] = count;
    }
  }

  // Collect the hits in a useful form,
  // and take note of the maximum time width
  fMaxWidth=0.0;
  for (unsigned int j = 0; j < allhits.size(); ++j){
    int dims = 3;//our point is defined by 3 elements:wire#,center of the hit, and the hit width
    std::vector<double> p(dims);
        
    double tickToDist = larp->DriftVelocity(larp->Efield(),larp->Temperature());
    tickToDist *= 1.e-3 * detp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    if (!wireids.size()) p[0] = (allhits[j]->WireID().Wire)*fWirePitch[allhits[j]->WireID().Plane];
    else p[0] = (wireids[j].Wire)*fWirePitch[allhits[j]->WireID().Plane];
    p[1] = allhits[j]->PeakTime()*tickToDist;
    p[2] = 2.*allhits[j]->RMS()*tickToDist;   //width of a hit in cm

    // check on the maximum width condition
    if ( p[2] > fMaxWidth ) fMaxWidth = p[2];
    
    fps.push_back(p);

    if (fClusterMethod) { // Using the R*-tree
      // Convert these same values into dbsPoints to feed into the R*-tree
      dbsPoint pp(p[0], p[1], 0.0, p[2]/2.0); // note dividing by two
      fRTree.Insert(j, pp.bounds());
      // Keep a parallel list already made up. We could use fps instead, but...
      fRect.push_back(pp);
    }
  }

  fpointId_to_clusterId.resize(fps.size(), kNO_CLUSTER); // Not zero as before!
  fnoise.resize(fps.size(), false);
  fvisited.resize(fps.size(), false);

  if (fClusterMethod) { // Using the R*-tree
    Visitor visitor = 
      fRTree.Query(RTree::AcceptAny(),Visitor());
    mf::LogInfo("DBscan") << "InitScan: hits RTree loaded with " 
			     << visitor.count << " items.";
  }
  mf::LogInfo("DBscan") << "InitScan: hits vector size is " << fps.size();

  return;
}

//----------------------------------------------------------
double cluster::DBScanAlg::getSimilarity(const std::vector<double> v1, const std::vector<double> v2){
  
   
  //for Euclidean distance comment everything out except this-->>>
  // return std::sqrt((v2[1]-v1[1])*(v2[1]-v1[1])+(v2[0]-v1[0])*(v2[0]-v1[0]));
  //------------------------------------------------------------------------
  // return std::abs( v2[0]-v1[0]); //for rectangle
  //---------------------------------------------------------------------- 
  //Manhattan distance:
  //return std::abs(v1[0]-v2[0])+std::abs(v1[1]-v2[1]);
  
  /// \todo this code assumes that all planes have the same wire pitch
  double wire_dist = fWirePitch[0];

  unsigned int wire1=(unsigned int)(v1[0]/wire_dist+0.5); //to make sure to get desired integer
  unsigned int wire2=(unsigned int)(v2[0]/wire_dist+0.5);
  int wirestobridge=0;

  if (wire1>wire2) {
    unsigned int wire = wire1;
    wire1 = wire2;
    wire2 = wire;
  }

  for(unsigned int i=wire1;i<wire2;i++){
    if(fBadChannels.find(i) != fBadChannels.end())
      wirestobridge++;
  }    
  
  double cmtobridge=wirestobridge*wire_dist;  
  //---------------------------------------------------------------------
  return (( std::abs(v2[0]-v1[0])-cmtobridge)*( std::abs(v2[0]-v1[0])-cmtobridge)); //for ellipse
}

//----------------------------------------------------------------
double cluster::DBScanAlg::getSimilarity2(const std::vector<double> v1, const std::vector<double> v2){

  //-------------------------------------------
  //return std::abs( v2[1]-v1[1]);//for rectangle
  //------------------------------------------

  /// \todo this code assumes all planes have the same wire pitch
  double wire_dist = fWirePitch[0];

  unsigned int wire1=(unsigned int)(v1[0]/wire_dist+0.5); //to make sure to get desired integer
  unsigned int wire2=(unsigned int)(v2[0]/wire_dist+0.5);
  int wirestobridge=0;

  if (wire1>wire2) {
    unsigned int wire = wire1;
    wire1 = wire2;
    wire2 = wire;
  }

  for(unsigned int i=wire1;i<wire2;i++){
    if(fBadChannels.find(i) != fBadChannels.end())
      wirestobridge++;
  }    
  
  double cmtobridge=wirestobridge*wire_dist;  
  
  if (std::abs(v2[0]-v1[0])>1e-10){
    cmtobridge *= std::abs((v2[1]-v1[1])/(v2[0]-v1[0]));
  }
  else cmtobridge = 0;

  return (( std::abs(v2[1]-v1[1])-cmtobridge)*( std::abs(v2[1]-v1[1])-cmtobridge));//for ellipse
  
  
}

//----------------------------------------------------------------
double cluster::DBScanAlg::getWidthFactor(const std::vector<double> v1, const std::vector<double> v2){
 
  //double k=0.13; //this number was determined by looking at flat muon hits' widths. 
                   //The average width of these hits in cm is 0.505, so 4*2*(w1^2)=2.04 
                   //where w1=w2=0.505, e^2.044= 7.69. In order not to change the distance 
                   //in time direction of the ellipse we want to make it equal to 1 for 
                   //these hits. Thus the k factor is k=1/7.69=0.13//for coeff=4

  //double k=0.78;
  //..................................................
  double k = 0.1;//for 4.5 coeff
  double WFactor = (exp(4.6*(( v1[2]*v1[2])+( v2[2]*v2[2]))))*k;
  //........................................................
  //Let's try something different:
  // double k=1.96;
  // double WFactor=(( v1[2]*v1[2])+( v2[2]*v2[2]))*k;
  if(WFactor > 1){
    if(WFactor < 6.25) return WFactor;  //remember that we are increasing the distance in 
                                        //eps2 as std::sqrt of this number (i.e std::sqrt(6.25))
    else return 6.25;
   
  }
  else return 1.0;  
}

//----------------------------------------------------------------
//\todo this is O(n) in the number of hits, while the high performance
//      claimed for DBSCAN relies on it being O(log n)!
std::vector<unsigned int> cluster::DBScanAlg::findNeighbors( unsigned int pid, 
							  double threshold,
							  double threshold2) {
  std::vector<unsigned int> ne;
  
  for ( int unsigned j=0; j < fsim.size(); j++){
    if((pid != j ) 
       && (((fsim[pid][j])/ (threshold*threshold))
	   + ((fsim2[pid][j])/ (threshold2*threshold2*(fsim3[pid][j]))))<1){ //ellipse
      ne.push_back(j);
    }
  }// end loop over fsim
  
  return ne;
}

//-----------------------------------------------------------------
void cluster::DBScanAlg::computeSimilarity()
{
  int size = fps.size();
  fsim.resize(size, std::vector<double>(size));
  for ( int i=0; i < size; i++){
    for ( int j=i+1; j < size; j++){
      fsim[j] [i] = fsim[i][ j] = getSimilarity(fps[i], fps[j]);
    }
  }
}

//------------------------------------------------------------------
void cluster::DBScanAlg::computeSimilarity2()
{
  int size = fps.size();
  fsim2.resize(size, std::vector<double>(size));
  for ( int i=0; i < size; i++){
    for ( int j=i+1; j < size; j++){
      fsim2[j] [i] = fsim2[i][ j] = getSimilarity2(fps[i], fps[j]);
    }
  }
}

//------------------------------------------------------------------
void cluster::DBScanAlg::computeWidthFactor()
{
  int size = fps.size();
  fsim3.resize(size, std::vector<double>(size));
       
  for ( int i=0; i < size; i++){
    for ( int j=i+1; j < size; j++){
      fsim3[j] [i] = fsim3[i][ j] = getWidthFactor(fps[i], fps[j]);
    }
  }
}


//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
// Run the selected clustering algorithm
void cluster::DBScanAlg::run_cluster() {
  switch(fClusterMethod) {
  case 2: 
    return run_dbscan_cluster();
  case 1:
    return run_FN_cluster();
  default:
    computeSimilarity();  // watch out for this, they are *slow*
    computeSimilarity2(); // "
    computeWidthFactor(); // "
    return run_FN_naive_cluster();
  }
}

//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
//
//  DWM's implementation of DBScanAlg as much like the paper as possible
void cluster::DBScanAlg::run_dbscan_cluster() {
  unsigned int cid = 0;
  // foreach pid
  for (size_t pid = 0; pid < fps.size(); pid++){
    // not already visited
    if (fpointId_to_clusterId[pid] == kNO_CLUSTER) {
      if ( ExpandCluster(pid,cid) ) {
	cid++;
      }
    } // if (!visited
  } // for
  //  END DBSCAN  

  // Construct clusters, count noise, etc..
  int noise = 0;
  fclusters.resize(cid);
  for(size_t y = 0; y < fpointId_to_clusterId.size(); ++y){
    if (fpointId_to_clusterId[y] == kNO_CLUSTER) {
      // This shouldn't happen...all points should be clasified by now!
      mf::LogWarning("DBscan") << "Unclassified point!";
    } 
    else if (fpointId_to_clusterId[y]==kNOISE_CLUSTER) {
      ++noise;
    } 
    else {
      unsigned int c = fpointId_to_clusterId[y];
      if (c >= cid) {
	mf::LogWarning("DBscan") << "Point in cluster " << c 
			      << " when only " << cid 
			      << " clusters wer found [0-" << cid-1
				 << "]";
      }
      fclusters[c].push_back(y);
    }
  }  
  mf::LogInfo("DBscan") << "DWM (R*-tree): Found " 
			   << cid << " clusters...";
  for (unsigned int c = 0; c < cid; ++c){
    mf::LogVerbatim("DBscan") << "\t" << "Cluster " << c << ":\t" 
			     << fclusters[c].size();
  }
  mf::LogVerbatim("DBscan") << "\t" << "...and " << noise << " noise points.";
}

//----------------------------------------------------------------
// Find the neighbos of the given point
std::set<unsigned int> cluster::DBScanAlg::RegionQuery(unsigned int point){
  dbsPoint region(fRect[point]);
  Visitor visitor = 
//     fRTree.Query(RTree::AcceptOverlapping(region.bounds()),Visitor());
//    fRTree.Query(AcceptEllipse(region.bounds(),fEps,fEps2),Visitor());
  fRTree.Query(AcceptFindNeighbors(region.bounds(),
				   fEps,fEps2,
				   fMaxWidth,fWirePitch[0],//\todo
				   fBadWireSum),           //assumes
	       Visitor());				   //equal
							   //pitch
  return visitor.sResult;
}
//----------------------------------------------------------------
// Find the neighbos of the given point
std::vector<unsigned int> cluster::DBScanAlg::RegionQuery_vector(unsigned int point){
  dbsPoint region(fRect[point]);
  Visitor visitor = 
//     fRTree.Query(RTree::AcceptOverlapping(region.bounds()),Visitor());
//    fRTree.Query(AcceptEllipse(region.bounds(),fEps,fEps2),Visitor());
  fRTree.Query(AcceptFindNeighbors(region.bounds(),
				   fEps,fEps2,
				   fMaxWidth,fWirePitch[0],//\todo
				   fBadWireSum),           //assumes
	       Visitor());				   //equal
							   //pitch
  std::vector<unsigned int> &v = visitor.vResult;
  // find neighbors insures that the called point is not in the
  // returned and this is intended as a drop-in replacement, so insure
  // this condition
  v.erase(std::remove(v.begin(), v.end(), point),v.end());
  return v;
}


//----------------------------------------------------------------
// Try to make a new cluster on the basis of point
bool cluster::DBScanAlg::ExpandCluster(unsigned int point,
				       unsigned int clusterID)
{
  /* GetSetOfPoints for point*/
  //std::vector<unsigned int> ne = findNeighbors(point, fEps,fEps2);
  std::set< unsigned int > seeds = RegionQuery(point);

  // not enough support -> mark as noise
  if (seeds.size() < fMinPts){
    fpointId_to_clusterId[point] = kNOISE_CLUSTER;
    return false;
  } else {
    // Add to the currecnt cluster
    fpointId_to_clusterId[point]=clusterID;
    for (std::set<unsigned int>::iterator itr = seeds.begin(); itr != seeds.end(); itr++){
      fpointId_to_clusterId[*itr]=clusterID;
    }
    seeds.erase(point);
    while (!seeds.empty()) {
      unsigned int currentP = *(seeds.begin());
      std::set< unsigned int > result = RegionQuery(currentP);
      
      if (result.size() >= fMinPts){
	for (std::set<unsigned int>::iterator itr = result.begin();
	     itr != result.end();
	     itr++){
	  unsigned int resultP = *itr;
	  // not already assigned to a cluster
	  if (fpointId_to_clusterId[resultP] == kNO_CLUSTER ||
	      fpointId_to_clusterId[resultP] == kNOISE_CLUSTER ){
	    if (fpointId_to_clusterId[resultP] == kNO_CLUSTER ) {
	      seeds.insert(resultP);
	    }
	    fpointId_to_clusterId[resultP]=clusterID;
	  } // unclassified or noise
	} // for
      } // enough support
      seeds.erase(currentP);
    } // while
    return true;
  }
}

//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
//
// The original findNeignbor based code converted to use a R*-tree,
// but not rearranged
void cluster::DBScanAlg::run_FN_cluster() 
{

  unsigned int cid = 0;
  // foreach pid
  for (size_t pid = 0; pid < fps.size(); pid++){
    // not already visited
    if (!fvisited[pid]){  
      
      fvisited[pid] = true;
      // get the neighbors
      //std::vector<unsigned int> ne = findNeighbors(pid, fEps,fEps2);
      std::vector<unsigned int> ne = RegionQuery_vector(pid);

      // not enough support -> mark as noise
      if (ne.size() < fMinPts){
	fnoise[pid] = true;
      }     
      else{
	// Add p to current cluster
	
	std::vector<unsigned int> c;              // a new cluster
	
	c.push_back(pid);   	// assign pid to cluster
	fpointId_to_clusterId[pid]=cid; 
	// go to neighbors
	for (size_t i = 0; i < ne.size(); ++i){
	  unsigned int nPid = ne[i];
	  
	  // not already visited
	  if (!fvisited[nPid]){
	    fvisited[nPid] = true;
	    // go to neighbors
	    //std::vector<unsigned int> ne1 = findNeighbors(nPid, fEps, fEps2);
	    std::vector<unsigned int> ne1 = RegionQuery_vector(nPid);
	    // enough support
	    if (ne1.size() >= fMinPts){
		       	
	      // join
	      
	      for(size_t i = 0; i < ne1.size(); ++i){
		// join neighbord
		ne.push_back(ne1[i]); 
	      }
	    }
	  }
		
	  // not already assigned to a cluster
	  //if (!fpointId_to_clusterId[nPid]){
	  if (fpointId_to_clusterId[nPid] == kNO_CLUSTER ){
	    c.push_back(nPid);
	    fpointId_to_clusterId[nPid]=cid;
	  }
	}
	    
	fclusters.push_back(c);
	
	
	cid++;
      }
    } // if (!visited
  } // for
  

  int noise=0;
  //no_hits=fnoise.size();

  for(size_t y = 0;y < fpointId_to_clusterId.size(); ++y){
    //if  (fpointId_to_clusterId[y]==0) noise++;
    if (fpointId_to_clusterId[y]==kNO_CLUSTER) noise++;
  }  
  mf::LogInfo("DBscan") << "FindNeighbors (R*-tree): Found " 
			   << cid << " clusters...";
  for (unsigned int c = 0; c < cid; ++c){
    mf::LogVerbatim("DBscan") << "\t" << "Cluster " << c << ":\t" 
			     << fclusters[c].size();
  }
  mf::LogVerbatim("DBscan") << "\t" << "...and " << noise << " noise points.";

}

//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
//
// The original findNeighrbor-based code. 
void cluster::DBScanAlg::run_FN_naive_cluster() 
{

  unsigned int cid = 0;
  // foreach pid
  for (size_t pid = 0; pid < fps.size(); ++pid){
    // not already visited
    if (!fvisited[pid]){  
      
      fvisited[pid] = true;
      // get the neighbors
      std::vector<unsigned int> ne = findNeighbors(pid, fEps,fEps2);
      
      // not enough support -> mark as noise
      if (ne.size() < fMinPts){
	fnoise[pid] = true;
      }     
      else{
	// Add p to current cluster
	
	std::vector<unsigned int> c;              // a new cluster
	
	c.push_back(pid);   	// assign pid to cluster
	fpointId_to_clusterId[pid] = cid; 
	// go to neighbors
	for (size_t i = 0; i < ne.size(); ++i){
	  unsigned int nPid = ne[i];
	  
	  // not already visited
	  if (!fvisited[nPid]){
	    fvisited[nPid] = true;
	    // go to neighbors
	    std::vector<unsigned int> ne1 = findNeighbors(nPid, fEps, fEps2);
	    // enough support
	    if (ne1.size() >= fMinPts){
		       	
	      // join
	      
	      for(unsigned int i=0;i<ne1.size();i++){
		// join neighbord
		ne.push_back(ne1[i]); 
	      }
	    }
	  }
		
	  // not already assigned to a cluster
	  //if (!fpointId_to_clusterId[nPid]){
	  if (fpointId_to_clusterId[nPid] == kNO_CLUSTER){
	    c.push_back(nPid);
	    fpointId_to_clusterId[nPid]=cid;
	  }
	}
	
	fclusters.push_back(c);
	
	
	cid++;
      }
    } // if (!visited
  } // for
  
  
  int noise = 0;
  //no_hits=fnoise.size();

  for(size_t y = 0; y < fpointId_to_clusterId.size(); ++y){
    //if  (fpointId_to_clusterId[y]==0) noise++;
    if  (fpointId_to_clusterId[y] == kNO_CLUSTER) ++noise;
  }
  mf::LogInfo("DBscan") << "FindNeighbors (naive): Found " << cid 
			   << " clusters...";
  for (unsigned int c = 0; c < cid; ++c){
    mf::LogVerbatim("DBscan") << "\t" << "Cluster " << c << ":\t" 
			     << fclusters[c].size() << " points";
  }
  mf::LogVerbatim("DBscan") << "\t" << "...and " << noise << " noise points.";
  
}

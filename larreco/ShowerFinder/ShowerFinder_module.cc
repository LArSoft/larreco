////////////////////////////////////////////////////////////////////////
// ShowerFinder class
//
// \author roxanne.guenette@yale.edu
//
////////////////////////////////////////////////////////////////////////

#include <string>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <math.h>

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"


// ROOT
#include "TMath.h"


///shower finding
namespace shwf {

  class ShowerFinder : public art::EDProducer {

  public:

    explicit ShowerFinder(fhicl::ParameterSet const&);
    virtual ~ShowerFinder();

    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt);

  private:

    std::string fVertexModuleLabel;         ///< label of module finding 2D endpoint
    std::string fClusterModuleLabel;        ///< label of module finding clusters
    std::string fHoughLineModuleLabel;      ///< label of module finding hough line
    std::string fVertexStrengthModuleLabel; ///< label of module finding 2D endpoint
    double      fRcone;                     ///< radious of cone for method
    double      fLcone;                     ///< length of the cone

  protected:


  }; // class ShowerFinder


}

namespace shwf{

  //-------------------------------------------------
  ShowerFinder::ShowerFinder(fhicl::ParameterSet const& pset)  :
    EDProducer{pset}
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Shower>                >();
    produces< art::Assns<recob::Shower, recob::Cluster> >();
    produces< art::Assns<recob::Shower, recob::Hit>     >();
  }


  //-------------------------------------------------
  ShowerFinder::~ShowerFinder()
  {
  }

  //-------------------------------------------------
  void ShowerFinder::reconfigure(fhicl::ParameterSet const& pset)
  {
    fVertexModuleLabel        = pset.get<std::string > ("VertexModuleLabel");
    fClusterModuleLabel       = pset.get<std::string > ("ClusterModuleLabel");
    fHoughLineModuleLabel     = pset.get<std::string > ("HoughLineModuleLabel");
    fVertexStrengthModuleLabel= pset.get<std::string > ("VertexStrengthModuleLabel");
    fRcone                    = pset.get<double      > ("Rcone");
    fLcone                    = pset.get<double      > ("Lcone");

    return;
  }

  //
  //-------------------------------------------------
  /// \todo This method appears to produce a recob::Cluster really as it is
  /// \todo a collection of 2D clusters from single planes
  void ShowerFinder::produce(art::Event& evt)
  {

    std::unique_ptr<std::vector<recob::Shower> > showercol(new std::vector<recob::Shower>);
    std::unique_ptr< art::Assns<recob::Shower, recob::Cluster> > cassn(new art::Assns<recob::Shower, recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Shower, recob::Hit>     > hassn(new art::Assns<recob::Shower, recob::Hit>);

    // Read in the vertex List object(s).
    art::Handle< std::vector<recob::EndPoint2D> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);

    // Read in the hough line List object(s).
    art::Handle< std::vector<recob::Cluster> > houghListHandle;
    evt.getByLabel(fHoughLineModuleLabel,houghListHandle);

    // Read in the cluster List object(s).
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);

    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);

    // Read in the vertex Strength List object(s).
    art::Handle< std::vector<recob::EndPoint2D> > vertexStrengthListHandle;
    evt.getByLabel(fVertexStrengthModuleLabel,vertexStrengthListHandle);

    std::vector<size_t> protoShowers; //vector of indices of clusters associated to a cone

    std::vector< art::Ptr<recob::Hit> > clusterhits; //hits in the cluster

    art::ServiceHandle<geo::Geometry const> geom;

    //This vector will contain all strong and strongest vertices
    art::PtrVector<recob::EndPoint2D> vertSel;

    //This loop is going over all the vertices in the event
    //and is interested in ONLY strong and strongest vertices.
    mf::LogInfo("ShowerFinder") << "Vertex STRENGTH list size = " << vertexStrengthListHandle->size()
				<< " AND vertices:" << vertexListHandle->size()
				<< "\nCLUSTER list size = " << clusterListHandle->size()
				<< " AND Hough: :" << houghListHandle->size();

    for(size_t iv = 0; iv < vertexListHandle->size(); ++iv){
      art::Ptr<recob::EndPoint2D> vertex(vertexListHandle, iv);
      //std::cout << "Vertex " << iv << " :  str = " << vertex->ID() << std::endl;
      //if(vertex->Strength() == 4 || vertex->Strength() == 3){
      if(vertex->ID() == 1 || vertex->Strength() == 3){ //only use Strongest and strong
	vertSel.push_back(vertex);
      }
      else continue;
    }

    //Definition of the geometry of the cone (which is basically a triangle)
    double scan_angle = 0; //angle of the scan steps
    double xa_cone = 0; // x coordinate of the cone's apex (wire number)
    double ya_cone = 0; // y coordinate of the cone's apex (drift time)
    double x1_cone = 0; // x coordinate of the cone's top right point (wire number)
    double y1_cone = 0; // y coordinate of the cone's top right point (drift time)
    double x2_cone = 0; // x coordinate of the cone's top left point (wire number)
    double y2_cone = 0; // y coordinate of the cone's top left point  (drift time)

    //The length of the side of the cone
    double fScone = std::sqrt( (fRcone*fRcone) + (fLcone*fLcone));

    // Opening angle of the cone (defined from input parameters)
    double cone_angle = (TMath::ATan(fRcone / fLcone)) / 2.0;
    mf::LogInfo("ShowerFinder") << "Cone Opening Angle: " << (180.0*cone_angle)/TMath::Pi();
    double compl_angle = 0;

    unsigned int n_scan =1 + (int)(TMath::Pi() / (2.0*cone_angle));
    mf::LogInfo("ShowerFinder") << "N scan: " << n_scan;

    double x_hit = 0; //x coordinate of hit
    double y_hit = 0; //y coordinate of hit

    int hits_cluster_counter = 0; //count the number of hits in a cluster that is inside a cone
    //int hits_cluster_Total = 0; //The total number of hits in a cluster

    // For EVERY vertex, the algorithm is going to scan the plane to find clusters
    // contained in the scanning cones

    for(size_t ivert = 0; ivert < vertSel.size(); ++ivert){

      mf::LogInfo("ShowerFinder") << "Number of STRONG vertices = " << vertSel.size();

      //get the coordinates of the vertex for the summit of the cone
      xa_cone = vertSel[ivert]->WireID().Wire;   //for update to EndPoint2D ... WK 4/22/13
      ya_cone = vertSel[ivert]->DriftTime();

      mf::LogInfo("ShowerFinder") << "Vertex at: (" << xa_cone << ", " << ya_cone << ")";

      //Beginning of the scan!
      for(unsigned int iscan = 0; iscan < n_scan; ++iscan){

	mf::LogInfo("ShowerFinder") << ">>>> Start SCAN: " << iscan;

	//define the scan anlge
	scan_angle = (TMath::Pi()/2.0) - (iscan*(2.0*cone_angle));

	mf::LogInfo("ShowerFinder") << "Scan Angle: " << (180.*scan_angle)/TMath::Pi();

	//get the complementary angle for geometry puurposes
	compl_angle = scan_angle - cone_angle;

	//Calculate the coordinates of the top right corner of the cone
	x1_cone = xa_cone + fScone*(std::cos(compl_angle));
	y1_cone = ya_cone + fScone*(std::sin(compl_angle));

	//Calculate the coordinates of the top left corner of the cone
	x2_cone = xa_cone + fScone*(std::cos(scan_angle + cone_angle));
	y2_cone = ya_cone + fScone*(std::sin(scan_angle + cone_angle));

	//Looking if a cluster is in this cone (loop over all hits of all clusters)
	protoShowers.clear();
	for(size_t iclust = 0; iclust < clusterListHandle->size(); ++iclust){

	//  art::Ptr<recob::Cluster> clust(clusterListHandle, iclust);
	  recob::Cluster const& clust = clusterListHandle->at(iclust);

	  //Get the hits vector from the cluster
	  clusterhits = fmh.at(iclust);
	  if(clusterhits.size() == 0) continue;

	  //Loop over ALL hits in the cluster. Looking if the cluster's
	  // hit is comprised in the cone
	  for(size_t ihits = 0; ihits < clusterhits.size(); ++ihits){


	    x_hit = clusterhits[ihits]->WireID().Wire;   //for update to EndPoint2D ... WK 4/22/13
	    y_hit = clusterhits[ihits]->PeakTime();

	    // Check in hits is INSIDE cone

	    //define the 2 line equations:
	    if(y_hit <= ((y2_cone - ya_cone)/(x2_cone - xa_cone))*x_hit + ya_cone &&
	       y_hit >= ((y1_cone - ya_cone)/(x1_cone - xa_cone))*x_hit + ya_cone){
	      hits_cluster_counter++;
	    }

	  }//end hits loop

	  //If there is more than 50% if the cluster INSIDE the cone, this is a protoshower
	  if(clusterhits.size() == 0) continue;
	  if(((double)hits_cluster_counter / (double)clusterhits.size()) >= 0.5){
	    mf::LogInfo("ShowerFinder") << "GOT A SHOWER!!!  in scan " << iscan
					<< "  cluster: " << iclust << " : " << clust.ID();

	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    /// \todo NEED TO TAKE OUT THE HOUGH LINES FROM THE PROTOSHOWERS!!!!!
	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    protoShowers.push_back(iclust);
	  }
	  clusterhits.clear();
	  hits_cluster_counter = 0;

	} //end cluster loop

	if(protoShowers.empty()) continue;

	// loop over hits in the protoShowers to determine the total charge of the shower
	double totalCharge = 0.;

	for(size_t p = 0; p < protoShowers.size(); ++p){
	  const size_t psIndex = protoShowers[p];
	  for (art::Ptr<recob::Hit> const& hit: fmh.at(psIndex))
	    if(hit->SignalType() == geo::kCollection) totalCharge += hit->Integral();
	}

	/// \todo really need to determine the values of the arguments of the recob::Shower ctor
	// fill with bogus values for now
	//double dcosVtx[3]    = { util::kBogusD };
	//double dcosVtxErr[3] = { util::kBogusD };
	//double maxTransWidth[2] = { util::kBogusD };
	//double distMaxWidth = util::kBogusD;

	//showercol->push_back( recob::Shower(dcosVtx, dcosVtxErr, maxTransWidth, totalCharge, distMaxWidth) );
	showercol->push_back(recob::Shower());

	// associate the shower with its clusters
	util::CreateAssn(*this, evt, *cassn,
	  showercol->size() - 1, protoShowers.begin(), protoShowers.end());

	// get the hits associated with each cluster and associate those with the shower
	for(size_t p = 0; p < protoShowers.size(); ++p){
	  const size_t psIndex = protoShowers[p];
	  std::vector< art::Ptr<recob::Hit> > hits = fmh.at(psIndex);
	  util::CreateAssn(*this, evt, *showercol, hits, *hassn);
	}

      } //end scan loop
    } //end vertices loop

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //NEED TO SEPARATE THE SHOWERS FROM THE DIFFERENT VERTEX!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    mf::LogInfo("ShowerFinder") << "--->  Recorded shower =  "<< showercol->size();
      /// \todo check if protoshower from further vertex is also contained in vertex nearer...
      //if the shower is stand alone ok, else, erase the next one
      //shower.SetID(is);
      //shower.SetVertexCoord(xa_cone, ya_cone);

    vertSel.clear();

    evt.put(std::move(showercol));
    evt.put(std::move(cassn));
    evt.put(std::move(hassn));

    return;
  } // end of produce

} // end of namespace


namespace shwf{

  DEFINE_ART_MODULE(ShowerFinder)

} //end of shower namespace

////////////////////////////////////////////////////////////////////////
//
// VertexFinder2D class
//
// tjyang@fnal.gov
//
// This algorithm is designed to reconstruct the vertices using the
// 2D cluster information
// 
// This is Preliminary Work and needs modifications
// ////////////////////////////////////////////////////////////////////////
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Vertex.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "TMath.h"
#include "TH1D.h"
#include "TVectorD.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "TVector3.h"

struct CluLen{
  int index;
  float length;
};

bool myfunction (CluLen c1, CluLen c2) { return (c1.length>c2.length);}

struct SortByWire {
  bool operator() (art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2) const { 
    return h1->Channel() < h2->Channel();
  }
};

///vertex reconstruction
namespace vertex {
   
 class VertexFinder2D :  public art::EDProducer {
    
  public:
    
    explicit VertexFinder2D(fhicl::ParameterSet const& pset); 
    virtual ~VertexFinder2D();        
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

    
    void produce(art::Event& evt);

  private:
    
    TH1D *dtIC;
  
    std::string fClusterModuleLabel;

  };
    
}

namespace vertex{

//-----------------------------------------------------------------------------
  VertexFinder2D::VertexFinder2D(fhicl::ParameterSet const& pset)
  {  
    this->reconfigure(pset);    
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::EndPoint2D, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
  }
//-----------------------------------------------------------------------------
  VertexFinder2D::~VertexFinder2D()
  {
  }

  //---------------------------------------------------------------------------
  void VertexFinder2D::reconfigure(fhicl::ParameterSet const& p) 
  {
    fClusterModuleLabel  = p.get< std::string >("ClusterModuleLabel");
    return;
  }
  //-------------------------------------------------------------------------
  void VertexFinder2D::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    dtIC = tfs->make<TH1D>("dtIC","It0-Ct0",100,-5,5);
    dtIC->Sumw2();
  }

// //-----------------------------------------------------------------------------
  void VertexFinder2D::produce(art::Event& evt)
  {

    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;
    // define TPC parameters
    TString tpcName = geom->GetLArTPCVolumeName();
    
    double YC =  (geom->DetHalfHeight())*2.;

    // wire angle with respect to the vertical direction
    double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; 
    
    // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
    double timetick = detprop->SamplingRate()*1.e-3; //time sample in us
    double presamplings = detprop->TriggerOffset(); //trigger offset

    double wire_pitch   = geom->WirePitch(0,1,0); //wire pitch in cm
    double Efield_drift = larprop->Efield();      // Electric Field in the drift region in kV/cm
    double Temperature  = larprop->Temperature(); // LAr Temperature in K
    
    //drift velocity in the drift region (cm/us)
    double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature); 

    //time sample (cm) 
    double timepitch = driftvelocity*timetick; 
    
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);

    art::PtrVector<recob::Cluster> clusters;
    for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }

    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);

    //Point to a collection of vertices to output.
    std::unique_ptr<std::vector<recob::Vertex> >                 vcol(new std::vector<recob::Vertex>);          //3D vertex
    std::unique_ptr<std::vector<recob::EndPoint2D> >             epcol(new std::vector<recob::EndPoint2D>);  //2D vertex
    std::unique_ptr< art::Assns<recob::EndPoint2D, recob::Hit> > assnep(new art::Assns<recob::EndPoint2D, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> >  assnsh(new art::Assns<recob::Vertex, recob::Shower>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track> >   assntr(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit> >     assnh(new art::Assns<recob::Vertex, recob::Hit>);

    // nplanes here is really being used as a proxy for the 
    // number of views in the detector
    int nplanes = geom->Views().size();
	
    std::vector< std::vector<int> > Cls(nplanes); //index to clusters in each view
    std::vector< std::vector<CluLen> > clulens(nplanes);

    std::vector<double> dtdwstart;
      
    //loop over clusters
    for(size_t iclu = 0; iclu < clusters.size(); ++iclu){
      
      float w0 = clusters[iclu]->StartWire();
      float w1 = clusters[iclu]->EndWire();
      float t0 = clusters[iclu]->StartTick();
      float t1 = clusters[iclu]->EndTick();
//      t0 -= detprop->GetXTicksOffset(clusters[iclu]->View(),0,0);
//      t1 -= detprop->GetXTicksOffset(clusters[iclu]->View(),0,0);

      CluLen clulen;
      clulen.index = iclu;
      clulen.length = sqrt(pow((w0-w1)*wire_pitch,2)+pow(detprop->ConvertTicksToX(t0,clusters[iclu]->View(),0,0)-detprop->ConvertTicksToX(t1,clusters[iclu]->View(),0,0),2));

      switch(clusters[iclu]->View()){
	
      case geo::kU :
	clulens[0].push_back(clulen);
	break;
      case geo::kV :
	clulens[1].push_back(clulen);
	break;
      case geo::kZ :
	clulens[2].push_back(clulen);
	break;
      default :
	break;
      }

      std::vector<double> wires;
      std::vector<double> times;
      
      std::vector< art::Ptr<recob::Hit> > hit = fmh.at(iclu);
      std::sort(hit.begin(), hit.end(), SortByWire());
      int n = 0;
      for(size_t i = 0; i < hit.size(); ++i){
	wires.push_back(hit[i]->WireID().Wire);
	times.push_back(hit[i]->PeakTime());
	++n;
      }
      if(n>=2){
	TGraph *the2Dtrack = new TGraph(std::min(10,n),&wires[0],&times[0]);           
	try{
	  the2Dtrack->Fit("pol1","Q");
	  TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
	  double par[2];
	  pol1->GetParameters(par);
	  //std::cout<<iclu<<" "<<par[1]<<" "<<clusters[iclu]->dTdW()<<std::endl;
	  dtdwstart.push_back(par[1]);
	}
	catch(...){
	  mf::LogWarning("VertexFinder2D") << "Fitter failed";
	  delete the2Dtrack;
	  dtdwstart.push_back(std::tan(clusters[iclu]->StartAngle()));
	  continue;
	}
	delete the2Dtrack;
      }
      else dtdwstart.push_back(std::tan(clusters[iclu]->StartAngle()));
    
    }
    
    //sort clusters based on 2D length
    for (size_t i = 0; i<clulens.size(); ++i){
      std::sort (clulens[i].begin(),clulens[i].end(), myfunction);
      for (size_t j = 0; j<clulens[i].size(); ++j){
	Cls[i].push_back(clulens[i][j].index);
      }
    }

    std::vector< std::vector<int> > cluvtx(nplanes);
    std::vector<double> vtx_w;
    std::vector<double> vtx_t;
    
    for (int i = 0; i < nplanes; ++i){
      if (Cls[i].size() >= 1){
	//at least one cluster
	//find the longest two clusters
	int c1 = -1;
	int c2 = -1;
	double ww0 = -999;
	double wb1 = -999;
	double we1 = -999;
	double wb2 = -999;
	double we2 = -999;
	double tt1 = -999;
	double tt2 = -999;
	double dtdw1 = -999;
	double dtdw2 = -999;
	double lclu1 = -999;
	double lclu2 = -999;
	for (unsigned j = 0; j<Cls[i].size(); ++j){
	  double lclu = std::sqrt(pow((clusters[Cls[i][j]]->StartWire()-clusters[Cls[i][j]]->EndWire())*13.5,2)
				  +pow(clusters[Cls[i][j]]->StartTick()-clusters[Cls[i][j]]->EndTick(),2));
	  bool rev = false;
	  bool deltaraylike = false;
	  bool enoughhits = false;
	  if (c1 != -1){
	    double wb = clusters[Cls[i][j]]->StartWire();
	    double we = clusters[Cls[i][j]]->EndWire();
	    double tt = clusters[Cls[i][j]]->StartTick();
	    double dtdw = dtdwstart[Cls[i][j]];
	    int nhits = fmh.at(Cls[i][j]).size();
	    ww0 = (tt-tt1+dtdw1*wb1-dtdw*wb)/(dtdw1-dtdw);
	    if (std::abs(wb1-ww0) > std::abs(we1-ww0)) rev = true;//reverse cluster dir
	    if ((!rev && ww0 > wb1+15)||(rev && ww0 < we1-15)) deltaraylike = true;
	    if (((!rev && ww0 > wb1+10)||(rev && ww0 < we1-10)) && nhits < 5) deltaraylike = true;
	    if (wb > wb1+20 && nhits < 20) deltaraylike = true;
	    if (wb > wb1+50 && nhits < 20) deltaraylike = true;
	    if (wb > wb1+8 && TMath::Abs(dtdw1-dtdw) < 0.15) deltaraylike = true;
	    if (std::abs(wb-wb1) > 30 && std::abs(we-we1) > 30) deltaraylike = true;
	    if (std::abs(tt-tt1) > 100) deltaraylike = true; //not really deltaray, but isolated cluster
	    //make sure there are enough hits in the cluster
	    //at leaset 2 hits if goes horizentally, at leaset 4 hits if goes vertically
	    double alpha = std::atan(dtdw);
	    if (nhits >= int(2+3*(1-std::abs(std::cos(alpha))))) enoughhits = true;
	    if (nhits < 5 && (ww0 < wb1-20 || ww0 > we1+20)) enoughhits = false;
	    
	  }
	  //do not replace the second cluster if the 3rd cluster is not consistent with the existing 2
	  bool replace = true;
	  if (c1 != -1 && c2 != -1){
	    double wb = clusters[Cls[i][j]]->StartWire();
	    double we = clusters[Cls[i][j]]->EndWire();
	    ww0 = (tt2-tt1+dtdw1*wb1-dtdw2*wb2)/(dtdw1-dtdw2);
	    if ((std::abs(ww0-wb1) < 10 || std::abs(ww0-we1) < 10) &&
		(std::abs(ww0-wb2) < 10 || std::abs(ww0-we2) < 10)){
	      if (std::abs(ww0-wb) > 15 && std::abs(ww0-we) > 15) replace = false;
	    }
	    //std::cout<<c1<<" "<<c2<<" "<<ww0<<" "<<wb1<<" "<<wb2<<" "<<wb<<" "<<we<<std::endl;
	  }
	  if (lclu1 < lclu){
	    if (c1 != -1 && !deltaraylike && enoughhits){
	      lclu2 = lclu1;
	      c2 = c1;
	      wb2 = wb1;
	      we2 = we1;
	      tt2 = tt1;
	      dtdw2 = dtdw1;
	    }
	    lclu1 = lclu;
	    c1 = Cls[i][j];
	    wb1 = clusters[Cls[i][j]]->StartWire();
	    we1 = clusters[Cls[i][j]]->EndWire();
	    tt1 = clusters[Cls[i][j]]->StartTick();
	    if (wb1>we1){
	      wb1 = clusters[Cls[i][j]]->EndWire();
	      we1 = clusters[Cls[i][j]]->StartWire();
	      tt1 = clusters[Cls[i][j]]->EndTick();
	    }
	    dtdw1 = dtdwstart[Cls[i][j]];
	  }
	  else if (lclu2 < lclu){
	    if (!deltaraylike && enoughhits && replace){
	      lclu2 = lclu;
	      c2 = Cls[i][j];
	      wb2 = clusters[Cls[i][j]]->StartWire();
	      we2 = clusters[Cls[i][j]]->EndWire();
	      tt2 = clusters[Cls[i][j]]->StartTick();
	      dtdw2 = dtdwstart[Cls[i][j]];
	    }
	  }
	}
	if (c1 != -1 && c2 != -1){
	  cluvtx[i].push_back(c1);
	  cluvtx[i].push_back(c2);
	  
	  double w1 = clusters[c1]->StartWire();
	  double t1 = clusters[c1]->StartTick();
	  if (clusters[c1]->StartWire()>clusters[c1]->EndWire()){
	    w1 = clusters[c1]->EndWire();
	    t1 = clusters[c1]->EndTick();
	  }
	  double k1 = dtdwstart[c1];
	  double w2 = clusters[c2]->StartWire();
	  double t2 = clusters[c2]->StartTick();
	  if (clusters[c2]->StartWire()>clusters[c2]->EndWire()){
	    w1 = clusters[c2]->EndWire();
	    t1 = clusters[c2]->EndTick();
	  }
	  double k2 = dtdwstart[c2];
//	  std::cout<<c1<<" "<<w1<<" "<<t1<<" "<<k1<<" "<<std::endl;
//	  std::cout<<c2<<" "<<w2<<" "<<t2<<" "<<k2<<" "<<std::endl;
	  //calculate the vertex
	  if (std::abs(k1-k2) < 0.5){
	    vtx_w.push_back(w1);
	    vtx_t.push_back(t1);
	  }
	  else{
	    double t0 = (k1*k2*(w1-w2)+k1*t2-k2*t1)/(k1-k2);
	    double w0 = (t2-t1+k1*w1-k2*w2)/(k1-k2);
	    vtx_w.push_back(w0);
	    vtx_t.push_back(t0);
	  }
	}
	else if (Cls[i].size() >= 1){
	  if (c1 != -1){
	    cluvtx[i].push_back(c1);
	    vtx_w.push_back(wb1);
	    vtx_t.push_back(tt1);
	  }
	  else{
	    cluvtx[i].push_back(Cls[i][0]);
	    vtx_w.push_back(clusters[Cls[i][0]]->StartWire());
	    vtx_t.push_back(clusters[Cls[i][0]]->StartTick());
	  }
	}
	//save 2D vertex
	// make an empty art::PtrVector of hits
	/// \todo should really get the actual vector of hits corresponding to end point
	/// \todo for now will get all hits from the current cluster
	std::vector< art::Ptr<recob::Hit> > hits = fmh.at(Cls[i][0]);
	double totalQ = 0.;
	for(size_t h = 0; h < hits.size(); ++h) totalQ += hits[h]->Integral();
	
	geo::WireID wireID(hits[0]->WireID().Cryostat,
			   hits[0]->WireID().TPC,
			   hits[0]->WireID().Plane,
			   (unsigned int)vtx_w.back());  //for update to EndPoint2D ... WK 4/22/13
	
	recob::EndPoint2D vertex(vtx_t.back(),
				 wireID, //for update to EndPoint2D ... WK 4/22/13
				 1,
				 epcol->size(),
				 clusters[Cls[i][0]]->View(),
				 totalQ);
	epcol->push_back(vertex);
	
	util::CreateAssn(*this, evt, *epcol, hits, *assnep);
	
      }
      else{
	//no cluster found
	vtx_w.push_back(-1);
	vtx_t.push_back(-1);
      }
    }
    //std::cout<<vtx_w[0]<<" "<<vtx_t[0]<<" "<<vtx_w[1]<<" "<<vtx_t[1]<<std::endl;
    
    Double_t vtxcoord[3];
    if (Cls[0].size()>0&&Cls[1].size()>0){//ignore w view
      double Iw0 = (vtx_w[0]+3.95)*wire_pitch;
      double Cw0 = (vtx_w[1]+1.84)*wire_pitch;
      
      double It0 = vtx_t[0] - presamplings;
      It0 *= timepitch;
      double Ct0 = vtx_t[1] - presamplings ;
      Ct0 *= timepitch;
      vtxcoord[0] = detprop->ConvertTicksToX(vtx_t[1],1,0,0);
      vtxcoord[1] = (Cw0-Iw0)/(2.*std::sin(Angle));
      vtxcoord[2] = (Cw0+Iw0)/(2.*std::cos(Angle))-YC/2.*std::tan(Angle);
      
      double yy,zz;       
      if(vtx_w[0]>=0&&vtx_w[0]<=239&&vtx_w[1]>=0&&vtx_w[1]<=239){
	if(geom->ChannelsIntersect(geom->PlaneWireToChannel(0,(int)((Iw0/wire_pitch)-3.95)),      
				   geom->PlaneWireToChannel(1,(int)((Cw0/wire_pitch)-1.84)),
				   yy,zz)){
	  //channelsintersect provides a slightly more accurate set of y and z coordinates. 
	  // use channelsintersect in case the wires in question do cross.
	  vtxcoord[1] = yy;
	  vtxcoord[2] = zz;
	}
	else{
	  vtxcoord[0] = -99999;
	  vtxcoord[1] = -99999;
	  vtxcoord[2] = -99999;
	}
      }	
      dtIC->Fill(It0-Ct0);
    }
    else{
      vtxcoord[0] = -99999;
      vtxcoord[1] = -99999;
      vtxcoord[2] = -99999;
    }
    
    /// \todo need to actually make tracks and showers to go into 3D vertex
    /// \todo currently just passing empty collections to the ctor
    art::PtrVector<recob::Track> vTracks_vec;
    art::PtrVector<recob::Shower> vShowers_vec;
    
    recob::Vertex the3Dvertex(vtxcoord, vcol->size());
    vcol->push_back(the3Dvertex);
    
    if(vShowers_vec.size() > 0){
      util::CreateAssn(*this, evt, *vcol, vShowers_vec, *assnsh);
      // get the hits associated with each track and associate those with the vertex
      ///\todo uncomment following lines when the shower vector actually contains showers from the art::Event
      // 	  for(size_t p = 0; p < vShowers_vec.size(); ++p){
      // 	    std::vector< art::Ptr<recob::Hit> > hits = fms.at(p);
      // 	    util::CreateAssn(*this, evt, *vcol, hits, *assnh);
      // 	  }
    }
    
    if(vTracks_vec.size() > 0){
      util::CreateAssn(*this, evt, *vcol, vTracks_vec, *assntr);
      // get the hits associated with each track and associate those with the vertex
      ///\todo uncomment following lines when the track vector actually contains tracks from the art::Event
      // 	  for(size_t p = 0; p < vTracks_vec.size(); ++p){
      // 	    std::vector< art::Ptr<recob::Hit> > hits = fmt.at(p);
      // 	    util::CreateAssn(*this, evt, *vcol, hits, *assnh);
      // 	  }
    }
    
    LOG_VERBATIM("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    LOG_VERBATIM("Summary") << "VertexFinder2D Summary:";
    for(size_t i = 0; i<epcol->size(); ++i) LOG_VERBATIM("Summary") << epcol->at(i) ;
    for(size_t i = 0; i<vcol->size(); ++i) LOG_VERBATIM("Summary") << vcol->at(i) ;
    
    evt.put(std::move(epcol));
    evt.put(std::move(vcol));
    evt.put(std::move(assnep));
    evt.put(std::move(assntr));
    evt.put(std::move(assnsh));
    evt.put(std::move(assnh));

  } // end of produce
} // end of vertex namespace

// //-----------------------------------------------------------------------------



namespace vertex{

  DEFINE_ART_MODULE(VertexFinder2D)

} // end of vertex namespace

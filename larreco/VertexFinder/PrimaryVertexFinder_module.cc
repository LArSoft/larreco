////////////////////////////////////////////////////////////////////////
//
// PrimaryVertexFinder class
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <iomanip>
#include <math.h>
#include <algorithm>
#include <vector>
#include <string>

#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"

///vertex reconstruction
namespace vertex {

 class PrimaryVertexFinder :  public art::EDProducer {

  public:

    explicit PrimaryVertexFinder(fhicl::ParameterSet const& pset);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);


    void produce(art::Event& evt);

  private:


    std::string fTrackModuleLabel;
    double      fVertexWindow;
    double      StartPointSeperation(recob::SpacePoint sp1, recob::SpacePoint sp2);
    bool        IsInVertexCollection(int a, std::vector<std::vector<int> > vertex_collection);
    int         IndexInVertexCollection(int a, int b, std::vector<std::vector<int> > vertex_collection);
    bool        IsInNewVertex(int a, std::vector<int> newvertex);
    double      gammavalue(TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2);
    double      alphavalue(double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2);
    double      MinDist(double alpha, double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2);
    TVector3    PointOnExtendedTrack(double alphagamma, TVector3 startpoint,  TVector3 dircos);
    TH2F*       fNoTracks;
    TH1F*       fLength_1stTrack;
    TH1F*       fLength_2ndTrack;
    TH1F*       fLength_3rdTrack;
    TH1F*       fLength_4thTrack;
    TH1F*       fLength_5thTrack;

  };

}

//------------------------------------------------------------------------------
bool sort_pred2(const std::pair<art::Ptr<recob::Track>,double>& left, const std::pair<art::Ptr<recob::Track>,double>& right)
{
  return left.second < right.second;
}

namespace vertex{

  //-----------------------------------------------------------------------------
  PrimaryVertexFinder::PrimaryVertexFinder(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
  {
    this->reconfigure(pset);
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Vertex, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
  }

  //---------------------------------------------------------------------------
  void PrimaryVertexFinder::reconfigure(fhicl::ParameterSet const& p)
  {
    fTrackModuleLabel  = p.get< std::string >("TrackModuleLabel");
    fVertexWindow      = p.get<double     >  ("VertexWindow");
  }

  //-------------------------------------------------------------------------
  void PrimaryVertexFinder::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    //    fNoVertices= tfs->make<TH2F>("fNoVertices", ";Event No; No of vertices", 100,0, 100, 30, 0, 30);
     fNoTracks= tfs->make<TH2F>("fNoTracks", ";Event No; No of Tracks", 10,0, 10, 10, 0, 10);
     fLength_1stTrack = tfs->make<TH1F>("fLength_Track1", "Muon Track Length", 100,0,100);
     fLength_2ndTrack = tfs->make<TH1F>("fLength_Track2", "2nd Track Length", 100,0,100);
     fLength_3rdTrack = tfs->make<TH1F>("fLength_Track3", "3rd Track Length", 100,0,100);
     fLength_4thTrack = tfs->make<TH1F>("fLength_Track4", "4th Track Length", 100,0,100);
     fLength_5thTrack = tfs->make<TH1F>("fLength_Track5", "5th Track Length", 100,0,100);
  }

// //-----------------------------------------------------------------------------
  void PrimaryVertexFinder::produce(art::Event& evt)
  {

    mf::LogInfo("PrimaryVertexFinder") << "------------------------------------------------------------------------------";

    //   std::cout << "run    : " << evt.Header().Run() << std::endl;
    //   std::cout << "subrun : " << evt.Header().Subrun() << std::endl;
    //std::cout << "event  : " << evt.Header().Event() << std::endl;

    mf::LogInfo("PrimaryVertexFinder") << "event  : " << evt.id().event();


    art::ServiceHandle<geo::Geometry const> geom;

    //mf::LogInfo("PrimaryVertexFinder") << "I am in Primary vertex finder " << std::endl;

    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrackModuleLabel,trackListHandle);

    //Point to a collection of vertices to output.
    std::unique_ptr< std::vector<recob::Vertex> > vcol(new std::vector<recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit> > vhassn(new art::Assns<recob::Vertex, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track> > vtassn(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> > vsassn(new art::Assns<recob::Vertex, recob::Shower>);


    std::vector<recob::Track> const& trkIn = *trackListHandle;

    mf::LogInfo("PrimaryVertexFinder") << "number of tracks in this event = " << trkIn.size();
    fNoTracks->Fill(evt.id().event(),trkIn.size());

    std::vector<recob::SpacePoint> startpoints_vec; // first space point of each track

    std::vector <TVector3> startvec;
    TVector3 startXYZ;

    std::vector <TVector3> endvec;
    TVector3 endXYZ;

    std::vector <TVector3> dircosvec;
    TVector3 dircosXYZ;

    std::vector< std::pair<art::Ptr<recob::Track>, double> > trackpair;

    art::FindMany<recob::SpacePoint> TrackSpacePoints
      (trackListHandle, evt, fTrackModuleLabel);

    for(unsigned int i = 0; i<trkIn.size(); ++i){
      recob::Track::Point_t start, end;
      std::tie(start, end) = trkIn[i].Extent();
      startXYZ.SetXYZ(start.X(),start.Y(),start.Z());
      endXYZ.SetXYZ(end.X(),end.Y(),end.Z());


      double length = (endXYZ-startXYZ).Mag();// (endvec[i]-startvec[i]).Mag();
      //mf::LogInfo("PrimaryVertexFinder") << "Track length calculated = " << length << std::endl;
      trackpair.push_back(std::pair<art::Ptr<recob::Track>,double>({ trackListHandle, i },length));
    }

    for(size_t i = 0; i<trackpair.size(); ++i){
      mf::LogInfo("PrimaryVertexFinder") << "track id is  = " << (trackpair[i].first)->ID()
					 << " track length = " << (trackpair[i].second);
    }

    std::sort(trackpair.rbegin(), trackpair.rend(), sort_pred2);

    mf::LogInfo("PrimaryVertexFinder") << "AFTER SORTING ";
    for(size_t i = 0; i < trackpair.size(); ++i){
      mf::LogInfo("PrimaryVertexFinder") << "track id is  = " << (trackpair[i].first)->ID()
					 << " track length = " << (trackpair[i].second);
    }

    if(trackpair.size()>0)
    fLength_1stTrack->Fill(trackpair[0].second);

    if(trackpair.size()>1)
    fLength_2ndTrack->Fill(trackpair[1].second);

    if(trackpair.size()>2)
    fLength_3rdTrack->Fill(trackpair[2].second);

    if(trackpair.size()>3)
    fLength_4thTrack->Fill(trackpair[3].second);

    if(trackpair.size()>4)
    fLength_5thTrack->Fill(trackpair[4].second);

    for(size_t j = 0; j < trackpair.size(); ++j) { //loop over tracks
      art::Ptr<recob::Track> const& track = trackpair[j].first;

      // the index of this track in the query is the same as its position in
      // the data product:
      std::vector<recob::SpacePoint const*> const& spacepoints
        = TrackSpacePoints.at(track.key());

      startXYZ  = trackpair[j].first->Vertex<TVector3>();
      endXYZ    = trackpair[j].first->End<TVector3>();
      dircosXYZ = trackpair[j].first->VertexDirection<TVector3>();

      startvec.push_back(startXYZ);
      endvec.push_back(endXYZ);
      dircosvec.push_back(dircosXYZ);

      mf::LogInfo("PrimaryVertexFinder") << "PrimaryVertexFinder got "<< spacepoints.size()
					 <<" 3D spacepoint(s) from Track3Dreco.cxx";

      // save the first SpacePoint of each Track... from now the SpacePoint ID represents the Track ID!!
      startpoints_vec.emplace_back(
        spacepoints[0]->XYZ(), spacepoints[0]->ErrXYZ(),
        spacepoints[0]->Chisq(), startpoints_vec.size()
        );

    }// loop over tracks

    for(size_t i = 0; i < startvec.size();  ++i){ //trackpair.size()
      mf::LogInfo("PrimaryVertexFinder") << "Tvector3 start point SORTED = ";
      startvec[i].Print();
    }
    for(size_t i = 0; i < dircosvec.size(); ++i){ //trackpair.size()
      mf::LogInfo("PrimaryVertexFinder") << "Tvector3 dir cos SORTED = ";
      dircosvec[i].Print();
    }

    std::vector<std::vector<int> > vertex_collection_int;
    std::vector <std::vector <TVector3> > vertexcand_vec;

    for (unsigned int i=0; i<trackpair.size(); ++i){
      for (unsigned int j=i+1; j<trackpair.size(); ++j){
	mf::LogInfo("PrimaryVertexFinder") << "distance between " << i << " and " << j
					   << " = "
					   << StartPointSeperation(startpoints_vec[i], startpoints_vec[j]);
	double GAMMA = gammavalue(startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	double ALPHA = alphavalue(GAMMA, startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	double MINDIST = MinDist(ALPHA, GAMMA, startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	mf::LogInfo("PrimaryVertexFinder") << "alpha = " << ALPHA << " gamma = "
					   << GAMMA << " MINIMUM DISTANCE = " << MINDIST;

	TVector3 TRACK1POINT = PointOnExtendedTrack(ALPHA, startvec[i], dircosvec[i]);
	TVector3 TRACK2POINT = PointOnExtendedTrack(GAMMA, startvec[j], dircosvec[j]);

	mf::LogInfo("PrimaryVertexFinder") << "POINTS ON THE TRACKS ARE:: ";
	TRACK1POINT.Print();
	TRACK2POINT.Print();

	//if(StartPointSeperation(startpoints_vec[i], startpoints_vec[j])<fVertexWindow){ ///// correct this
	//if(MINDIST<2 && trackpair[i].second >30 && trackpair[j].second >30){
	if(MINDIST < fVertexWindow && ((TRACK1POINT-startvec[i]).Mag()) < fVertexWindow){

	  if((!IsInVertexCollection(i, vertex_collection_int)) && (!IsInVertexCollection(j, vertex_collection_int))){
	    std::vector<int> newvertex_int;
	    std::vector <TVector3> vertexcand;
	    newvertex_int.push_back(i);
	    newvertex_int.push_back(j);
	    vertex_collection_int.push_back(newvertex_int);
	    //newvertex.clear();
	    vertexcand.push_back(TRACK1POINT);
	    vertexcand.push_back(TRACK2POINT);
	    vertexcand_vec.push_back(vertexcand);
	  }
	  else{
	    int index = IndexInVertexCollection(i, j, vertex_collection_int);
	    //mf::LogInfo("PrimaryVertexFinder") << "index where a new vertex will be added = " << index << std::endl;
	    if(!IsInNewVertex(i, vertex_collection_int[index])){
	      vertex_collection_int[index].push_back(i);
	      vertexcand_vec[index].push_back(TRACK1POINT); //need to fix for delta rays
	    }
	    if(!IsInNewVertex(j, vertex_collection_int[index])){
	      vertex_collection_int[index].push_back(j);
	      vertexcand_vec[index].push_back(TRACK2POINT); //need to fix for delta rays
	    }
	  }
	}// end else
      }
    }


    //now add the unmatched track IDs to the collection
    for(size_t i = 0; i < trackpair.size(); ++i){
      if(!IsInVertexCollection(i, vertex_collection_int)){
	//if(trackpair[i].second>30){
	std::vector<int> temp;
	std::vector <TVector3> temp1;
	temp.push_back(i);
	temp1.push_back(startvec[i]);
	vertex_collection_int.push_back(temp);
	vertexcand_vec.push_back(temp1);
	//}
      }
    }

    // indices (in their data products) of tracks and showers connected to the vertex
    std::vector<size_t> vTrackIndices, vShowerIndices;

    // find the hits of all the tracks
    art::FindManyP<recob::Hit> TrackHits(trackListHandle, evt, fTrackModuleLabel);

    // find the hits of all the showers
  //  art::FindManyP<recob::Hit> ShowerHits(showerListHandle, evt, fShowerModuleLabel);
    ///\todo replace with the real query when this module is updated to look for showers too
    art::FindManyP<recob::Hit> ShowerHits(std::vector<art::Ptr<recob::Shower>>(), evt, fTrackModuleLabel);

    for(size_t i = 0; i < vertex_collection_int.size(); ++i){
      double x = 0.;
      double y = 0.;
      double z = 0.;
      int elemsize = 0.;
      for(std::vector<int>::iterator itr = vertex_collection_int[i].begin(); itr < vertex_collection_int[i].end(); ++itr){
        mf::LogInfo("PrimaryVertexFinder") << "vector elements at index " << i << " are " << *itr
                                           << "\ntrack original ID = " << (trackpair[*itr].first)->ID();
        // save the index in the data product of this track
        vTrackIndices.push_back(trackpair[*itr].first.key());
      }
      mf::LogInfo("PrimaryVertexFinder") << "------------";


      for(std::vector<TVector3>::iterator itr = vertexcand_vec[i].begin(); itr < vertexcand_vec[i].end(); ++itr){
	//calculate sum of x, y and z of a vertex
	x += (*itr).X();
	y += (*itr).Y();
	z += (*itr).Z();
	elemsize = vertexcand_vec[i].size();
      }

      double avgx = x/elemsize;
      double avgy = y/elemsize;
      double avgz = z/elemsize;

      Double_t vtxcoord[3];
      vtxcoord[0] = avgx;
      vtxcoord[1] = avgy;
      vtxcoord[2] = avgz;

      recob::Vertex the3Dvertex(vtxcoord, vcol->size());
      vcol->push_back(the3Dvertex);

      if(!vTrackIndices.empty()){
        // associate the tracks and their hits with the vertex
        util::CreateAssn(*this, evt, *vtassn,
          vcol->size() - 1, vTrackIndices.begin(), vTrackIndices.end());
        for(size_t tIndex: vTrackIndices) {
          std::vector<art::Ptr<recob::Hit>> const& hits = TrackHits.at(tIndex);
          util::CreateAssn(*this, evt, *vcol, hits, *vhassn);
        }
        vTrackIndices.clear();
      } // if tracks

      if(!vShowerIndices.empty()){
        // associate the showers and their hits with the vertex
        util::CreateAssn(*this, evt, *vsassn,
          vcol->size() - 1, vShowerIndices.begin(), vShowerIndices.end());
        for(size_t sIndex: vShowerIndices){
          std::vector<art::Ptr<recob::Hit>> const& hits = ShowerHits.at(sIndex);
          util::CreateAssn(*this, evt, *vcol, hits, *vhassn);
        }
        vShowerIndices.clear();
      } // if showers


    }// end loop over vertex_collection_ind

    MF_LOG_VERBATIM("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    MF_LOG_VERBATIM("Summary") << "PrimaryVertexFinder Summary:";
    for(size_t i = 0; i < vcol->size(); ++i) MF_LOG_VERBATIM("Summary") << vcol->at(i) ;

    evt.put(std::move(vcol));
    evt.put(std::move(vtassn));
    evt.put(std::move(vhassn));
    evt.put(std::move(vsassn));

  } // end of produce
} // end of vertex namespace

// //-----------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::StartPointSeperation(recob::SpacePoint sp1, recob::SpacePoint sp2)
{
  double x= (sp2.XYZ()[0])-(sp1.XYZ()[0]);
  double y= (sp2.XYZ()[1])-(sp1.XYZ()[1]);
  double z= (sp2.XYZ()[2])-(sp1.XYZ()[2]);
  double distance = std::sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  return distance;
}
// //---------------------------------------------------------------------------------
bool vertex::PrimaryVertexFinder::IsInVertexCollection(int a, std::vector<std::vector<int> > vertex_collection)
{
  int flag = 0;

  for(unsigned int i = 0; i < vertex_collection.size() ; i++){
    for(std::vector<int>::iterator itr = vertex_collection[i].begin(); itr < vertex_collection[i].end(); ++itr){
      if (a == *itr){
	flag = 1;
	break;
      }
    }
  }
  if(flag==1)
    return true;
  return false;
}
// //------------------------------------------------------------------------------
int vertex::PrimaryVertexFinder::IndexInVertexCollection(int a, int b, std::vector<std::vector<int> > vertex_collection)
{
  int index = -1;
  for(unsigned int i = 0; i < vertex_collection.size() ; i++){
    for(std::vector<int>::iterator itr = vertex_collection[i].begin(); itr < vertex_collection[i].end(); ++itr){
      if (a == *itr || b == *itr)
	index = i;
    }
  }
  return index;
}
// //------------------------------------------------------------------------------
bool vertex::PrimaryVertexFinder::IsInNewVertex(int a, std::vector<int> newvertex)
{
  int flag = 0;
  for(unsigned int i = 0; i < newvertex.size() ; i++){
    if (a == newvertex[i]){
      flag = 1;
      break;
    }
  }

  if(flag==1)
    return true;
  return false;
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::gammavalue(TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  double gamma = ((startpoint1*dircos2)-(startpoint2*dircos2)+((dircos1*dircos2)*(startpoint2*dircos1))-((dircos1*dircos2)*(startpoint1*dircos1)))/(1-((dircos1*dircos2)*(dircos1*dircos2)));

  return gamma;
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::alphavalue(double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  double alpha = (gamma*(dircos1*dircos2)) + (startpoint2*dircos1) - (startpoint1*dircos1);

  return alpha;
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::MinDist(double alpha, double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  TVector3 mindis_vector = startpoint1 - startpoint2 + alpha*dircos1 - gamma*dircos2;
  double mindis = mindis_vector.Mag();
  return mindis;
}
// //------------------------------------------------------------------------------
TVector3 vertex::PrimaryVertexFinder::PointOnExtendedTrack(double alphagamma, TVector3 startpoint,  TVector3 dircos)
{
  TVector3 PointOnExtendedTrack = startpoint + (alphagamma * dircos);
  return PointOnExtendedTrack;
}
// //------------------------------------------------------------------------------


namespace vertex{

  DEFINE_ART_MODULE(PrimaryVertexFinder)

} // end of vertex namespace

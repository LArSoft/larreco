////////////////////////////////////////////////////////////////////
// Implementation of the EMShower algorithm
//
// Forms EM showers from clusters and associated tracks.
// Also provides methods for finding the vertex and further
// properties of the shower.
//
// Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
////////////////////////////////////////////////////////////////////

#ifndef EMShowerAlg_hxx
#define EMShowerAlg_hxx

// framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Core/FindManyP.h"

// larsoft
#include "Utilities/DetectorProperties.h"
#include "AnalysisAlg/CalorimetryAlg.h"
#include "RecoAlg/ShowerEnergyAlg.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Shower.h"
#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

// C++
#include <iostream>
#include <map>

// ROOT
#include "TVector2.h"
#include "TMath.h"

namespace shower {
  class EMShowerAlg;
  class HitPosition;
}

class shower::EMShowerAlg {
public:

  EMShowerAlg(fhicl::ParameterSet const& pset);
  void AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
				  art::FindManyP<recob::Hit> const& fmh,
				  art::FindManyP<recob::Track> const& fmt,
				  std::vector<int> const& clustersToIgnore,
				  std::map<int,std::vector<int> >& clusterToTracks,
				  std::map<int,std::vector<int> >& trackToClusters);
  void AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
				  art::FindManyP<recob::Hit> const& fmh,
				  art::FindManyP<recob::Track> const& fmt,
				  std::map<int,std::vector<int> >& clusterToTracks,
				  std::map<int,std::vector<int> >& trackToClusters);
  void CheckShowerPlanes(std::vector<std::vector<int> > const& initialShowers,
			 std::vector<int>& clustersToIgnore,
			 std::vector<art::Ptr<recob::Cluster> > const& clusters,
			 art::FindManyP<recob::Hit> const& fmh);
  void FindShowers(std::map<int,std::vector<int> > const& trackToClusters, std::vector<std::vector<int> >& showers);
  void FindInitialTrack(art::PtrVector<recob::Hit> const& hits, art::Ptr<recob::Track>& initialTrack, std::vector<art::Ptr<recob::Hit> >& initialTrackHits);
  recob::Shower MakeShower(art::PtrVector<recob::Hit> const& hits,
			   art::Ptr<recob::Track> const& initialTrack,
			   std::vector<art::Ptr<recob::Hit> > const& initialTrackHits);

private:

  std::vector<art::Ptr<recob::Hit> > FindInitialTrack(std::vector<art::Ptr<recob::Hit> > const& shower);
  void ProjectVertexIn2D(TVector3 const& vertex,
			 std::map<int,std::vector<art::Ptr<recob::Hit> > >& trackHitsMap,
			 std::map<int,std::pair<std::vector<art::Ptr<recob::Hit> >,std::vector<art::Ptr<recob::Hit> > > > const& trackHitsBothEndsMap);
  void FindShowerStartDirection(art::Ptr<recob::Track> const& vertexTrack, std::map<int,TVector2> const& showerCentreMap, TVector3& showerVertex, TVector3& showerDirection);
  art::Ptr<recob::Track> FindVertexTrack(std::map<int,art::Ptr<recob::Hit> > const& vertexMap,
					 std::map<int,art::Ptr<recob::Track> > const& trackMap,
					 std::map<int,std::vector<art::Ptr<recob::Hit> > > const& trackHitsMap);
  std::vector<art::Ptr<recob::Hit> > FindTrack(std::vector<art::Ptr<recob::Hit> > const& shower,
					       TVector2 const& start,
					       TVector2 const& end);
  void FindShowerEnds(std::vector<art::Ptr<recob::Hit> > const& shower,
		      art::Ptr<recob::Hit>& end1,
		      art::Ptr<recob::Hit>& end2);
  double FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, art::Ptr<recob::Track> const& track);
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);
  TVector2 HitPosition(art::Ptr<recob::Hit> const& hit);
  TVector2 HitPosition(TVector2 const& pos, geo::PlaneID planeID);
  double GlobalWire(geo::WireID wireID);
  TVector2 Project3DPointOntoPlane(TVector3 const& point, unsigned int plane);

  // Parameters
  double fMinTrackLength;
  double fdEdxTrackLength;

  // Services used by this class
  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;
  art::ServiceHandle<art::TFileService> tfs;

  // Algs used by this class
  shower::ShowerEnergyAlg fShowerEnergyAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  pma::ProjectionMatchingAlg fProjectionMatchingAlg;

};

class shower::HitPosition {
 public:

  TVector2 WireTick;
  TVector2 Cm;

  // default constructor
  HitPosition();
  // contructors
  HitPosition(TVector2 wiretick, TVector2 cm) {
    WireTick = wiretick;
    Cm = cm;
  }
  HitPosition(TVector2 wiretick, geo::PlaneID planeID) {
    WireTick = wiretick;
    Cm = ConvertWireTickToCm(wiretick, planeID);
  }

  TVector2 ConvertWireTickToCm(TVector2 wiretick, geo::PlaneID planeID) {
    return TVector2(wiretick.X() * fGeom->WirePitch(planeID),
		    fDetProp->ConvertTicksToX(wiretick.Y(), planeID));
  }

 private:

  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

};

#endif

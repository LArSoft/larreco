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
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/AnalysisAlg/CalorimetryAlg.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/Shower.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

// C++
#include <iostream>
#include <map>
#include <iterator>

// ROOT
#include "TVector2.h"
#include "TMath.h"

//temp
#include "TGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TString.h"

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
  void FindInitialTrack(art::PtrVector<recob::Hit> const& hits,
			std::unique_ptr<recob::Track>& initialTrack,
			std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits);
  void FindInitialTrack(art::PtrVector<recob::Hit> const& hits,
			std::unique_ptr<recob::Track>& initialTrack,
			std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits,
			art::FindManyP<recob::Cluster> const& fmc, int plane);
  recob::Shower MakeShower(art::PtrVector<recob::Hit> const& hits,
			   std::unique_ptr<recob::Track> const& initialTrack,
			   std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialTrackHits);
  recob::Shower MakeShower(art::PtrVector<recob::Hit> const& hits,
			   art::Ptr<recob::Vertex> const& vertex,
			   int & iok);
  void FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >const& showerHits,
			    art::Ptr<recob::Vertex> const& vertex,
			    std::vector<art::Ptr<recob::Hit> >& trackHits);
  std::unique_ptr<recob::Track> ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
					       std::vector<art::Ptr<recob::Hit> > const& track2);
  std::unique_ptr<recob::Track> ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
					       std::vector<art::Ptr<recob::Hit> > const& track2,
					       std::map<geo::PlaneID,TVector2> const& showerCentreMap);
  std::map<int,std::vector<art::Ptr<recob::Hit> > > FindShowerStart(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap);
  double OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
			 std::vector<art::Ptr<recob::Hit> >& orderedShower,
			 art::FindManyP<recob::Cluster> const& fmc);
  std::vector<art::Ptr<recob::Hit> > OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower);
  void OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
			 std::vector<art::Ptr<recob::Hit> >& orderedShower,
			 art::Ptr<recob::Vertex> const& vertex);
  double OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
					      std::vector<art::Ptr<recob::Hit> >& showerHits);
  TVector3 Construct3DPoint(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2);

  Int_t WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm);

  bool isCleanShower(std::vector<art::Ptr<recob::Hit> > const& hits);

private:

  std::vector<int> IdentifyBadPlanes(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap);
  std::vector<int> IdentifyBadPlanes(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap,
				     std::map<int,double> const& goodnessOfOrderMap);
  std::map<int,std::vector<art::Ptr<recob::Hit> > > CheckShowerHits(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap,
								    std::map<int,double> const& goodnessOfOrderMap);
  std::unique_ptr<recob::Track> MakeInitialTrack(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap,
						 std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerCentreMap);
  std::vector<double> GetShowerDirectionProperties(std::vector<art::Ptr<recob::Hit> > const& showerHits, TVector2 const& direction, std::string end);
  std::vector<art::Ptr<recob::Hit> > FindOrderOfHits(std::vector<art::Ptr<recob::Hit> > const& hits);
  double FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, std::unique_ptr<recob::Track> const& track);
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);
  TVector2 HitPosition(art::Ptr<recob::Hit> const& hit);
  TVector2 HitPosition(TVector2 const& pos, geo::PlaneID planeID);
  double GlobalWire(geo::WireID wireID);
  TVector2 Project3DPointOntoPlane(TVector3 const& point, geo::PlaneID planeID);

  // Parameters
  double fMinTrackLength;
  double fdEdxTrackLength;
  // Parameters to fit wire vs time
  unsigned int         fNfitpass;
  std::vector<unsigned int>     fNfithits;
  std::vector<double>  fToler;

  // Services used by this class
  art::ServiceHandle<geo::Geometry> fGeom;
  detinfo::DetectorProperties const* fDetProp;
  art::ServiceHandle<art::TFileService> tfs;

  // Algs used by this class
  shower::ShowerEnergyAlg fShowerEnergyAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  pma::ProjectionMatchingAlg fProjectionMatchingAlg;

  bool debug = false;

};

class shower::HitPosition {
 public:

  TVector2 WireTick;
  TVector2 Cm;

  // default constructor
  HitPosition();
  // contructors
  HitPosition(TVector2 wiretick, TVector2 cm): HitPosition()
  {
    WireTick = wiretick;
    Cm = cm;
  }
  HitPosition(TVector2 wiretick, geo::PlaneID planeID): HitPosition() {
    WireTick = wiretick;
    Cm = ConvertWireTickToCm(wiretick, planeID);
  }

  TVector2 ConvertWireTickToCm(TVector2 wiretick, geo::PlaneID planeID) {
    return TVector2(wiretick.X() * fGeom->WirePitch(planeID),
		    fDetProp->ConvertTicksToX(wiretick.Y(), planeID));
  }

 private:

  geo::GeometryCore const* fGeom = nullptr;
  detinfo::DetectorProperties const* fDetProp = nullptr;

};

#endif

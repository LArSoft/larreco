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

// LArSoft includes
#include "art/Framework/Core/FindManyP.h"
#include "AnalysisAlg/CalorimetryAlg.h"
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
#include "lardata/DataProviders/IDetectorProperties.h"

// C++
#include <iostream>
#include <map>

// ROOT
#include "TVector2.h"
#include "TMath.h"

namespace shower {
  class EMShowerAlg;
}

class shower::EMShowerAlg {
public:
  EMShowerAlg();

  void MakeShowers(std::map<int,std::vector<int> > const& trackToClusters, std::vector<std::vector<int> >& showers);
  void FindVertexTrack(art::Ptr<recob::Track>& vertexTrack, std::map<int,art::Ptr<recob::Hit> > const& vertexMap, std::map<int,art::Ptr<recob::Track> > const& trackMap, std::map<int,std::vector<int> > const& trackHitsMap);
  void FindShowerProperties(art::PtrVector<recob::Hit> const& hits, art::FindManyP<recob::Track> const& fmt, calo::CalorimetryAlg const& calo,
			    TVector3& direction, TVector3& directionError, TVector3& vertex, TVector3& vertexError,
			    std::vector<double>& totalEnergy, std::vector<double>& totalEnergyError, std::vector<double>& dEdx, std::vector<double>& dEdxError,
			    int& bestPlane);
  void FindShowerStartDirection(art::Ptr<recob::Track> const& vertexTrack, std::map<int,TVector2> const& showerCentreMap, TVector3& showerVertex, TVector3& showerDirection);
  double FinddEdx(art::PtrVector<recob::Hit> const& shower, art::Ptr<recob::Track> const& track, calo::CalorimetryAlg const& calo, geo::View_t const& view, std::vector<int> const& trackHits);
  double FindTotalEnergy(art::PtrVector<recob::Hit> const& hits, int plane);
  void FindShowerEnds(art::PtrVector<recob::Hit> const& shower, TVector2 const& centre, art::Ptr<recob::Hit>& end1, art::Ptr<recob::Hit>& end2);
  art::Ptr<recob::Hit> FindVertex(art::PtrVector<recob::Hit> const& shower, art::Ptr<recob::Hit> const& end1, art::Ptr<recob::Hit> const& end2);
  void FindVertex(art::PtrVector<recob::Hit> const& shower, TVector2 const& end1, TVector2 const& end2, std::vector<int>& trackHits);
  std::vector<int> FindTrack(art::PtrVector<recob::Hit> const& shower, TVector2 const& start, TVector2 const& end);
  void FindTrack(art::PtrVector<recob::Hit> const& shower, std::map<double,int> const& hitToEnd, std::vector<int>& trackHits);
  void FindTrack(TVector2 const& start, TVector2 const& end, std::map<int,std::vector<int> > const& hitWires, std::vector<int>& trackHits);
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);
  double GlobalWire(geo::WireID wireID);
  int FindTrackID(art::Ptr<recob::Hit> const& hit);
  int FindTrueTrack(std::vector<art::Ptr<recob::Hit> > const& showerHits);

private:

  art::ServiceHandle<geo::Geometry> fGeom;
  dataprov::IDetectorProperties const* fDetProp;

};

#endif

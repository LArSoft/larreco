////////////////////////////////////////////////////////////////////////
// Class:       TCShower
// Plugin Type: producer (art v2_11_02)
// File:        TCShower_module.cc
//
// Generated at Fri Jun  8 14:55:04 2018 by Rory Fitzpatrick using cetskelgen
// from cetlib version v3_03_01.
// 
// Contact: roryfitz@umich.edu
// 
// module produces showers by selecting tracks surround by many 
// showerLike trajectories as defined by trajcluster with negative
// cluster IDs 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"

#include <memory>

bool compare(const art::Ptr<recob::Track>& l, const art::Ptr<recob::Track>& r) //(2)
{
  return l->Vertex().Z() > r->Vertex().Z();
}

namespace shower {
  class TCShower;
}

class shower::TCShower : public art::EDProducer {
public:
  explicit TCShower(fhicl::ParameterSet const & p);

  TCShower(TCShower const &) = delete;
  TCShower(TCShower &&) = delete;
  TCShower & operator = (TCShower const &) = delete;
  TCShower & operator = (TCShower &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  int goodHit(art::Ptr<recob::Hit>, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2);
  
  int goodHit(art::Ptr<recob::Hit>, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2, int& pull);

  std::vector<TVector3> getSecondPoint(art::Ptr<recob::Track> thistrack);

  std::string fClusterModuleLabel;
  std::string fTrackModuleLabel;
  std::string fHitModuleLabel;

};

// -----------------------------------------------------

shower::TCShower::TCShower(fhicl::ParameterSet const & pset) : 
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel", "trajcluster" ) ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel", "pmtrack" ) ),
  fHitModuleLabel           (pset.get< std::string >("HitModuleLabel", "trajcluster" ) ) {

  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
}

// -----------------------------------------------------

void shower::TCShower::produce(art::Event & evt) {
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitShowerAssociations(new art::Assns<recob::Shower, recob::Hit>);

  // hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  // clusters
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  // tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

  // detector properties
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;

  // get associations
  art::FindManyP<recob::Hit> cls_fm(clusterListHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Hit> trk_fm(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Track> hit_fm(hitListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster> hitcls_fm(hitListHandle, evt, fClusterModuleLabel);

  // sort tracks
  std::vector<art::Ptr<recob::Track> > tracklistL20;
  std::vector<art::Ptr<recob::Track> > tracklistG20;

  std::vector<art::Ptr<recob::Track> > tracklistSorted;
  std::vector<size_t> trackIndicesSorted;

  for (size_t i = 0; i < tracklist.size(); ++i) {
    double length = tracklist[i]->Length();

    if (length < 20) tracklistL20.push_back(tracklist[i]);
    else tracklistG20.push_back(tracklist[i]);
  }

  std::sort(tracklistL20.begin(), tracklistL20.end(), compare);
  std::sort(tracklistG20.begin(), tracklistG20.end(), compare);
  
  // add all points to tracklistSorted
  for (size_t i = 0; i < tracklistL20.size(); ++i) {
    tracklistSorted.push_back(tracklistL20[i]);
  }
  for (size_t i = 0; i < tracklistG20.size(); ++i) {
    tracklistSorted.push_back(tracklistG20[i]);
  }

  std::reverse(tracklistSorted.begin(), tracklistSorted.end());

  for (size_t i = 0; i < tracklistSorted.size(); ++i) {
    for (size_t j = 0; j < tracklist.size(); ++j) {
      if (tracklistSorted[i]->ID() == tracklist[j]->ID()) trackIndicesSorted.push_back(j);
    } // loop through tracks
  } // loop through sorted tracks

  for (size_t i = 0; i < tracklistSorted.size(); ++i) {
    std::cout << tracklistSorted[i]->ID() << " " << tracklistSorted[i]->Length() << " " << tracklistSorted[i]->Vertex().Z() << std::endl; 
  }

  tracklist = tracklistSorted;

  for (size_t i = 0; i < tracklist.size(); ++i) {

    int tolerance = 100; // how many shower like cluster you need to define a shower
    double pullTolerance = 0.6; // hits should be evenly distributed around the track
    double maxDist = 10; // how far a shower like cluster can be from the track
    double minDistVert = 15; // exclude tracks near the vertex
    
    if (tracklist[i]->Length() < 20) continue;

    // adjust tolerances for short tracks
    if (tracklist[i]->Length() < 50) {
      tolerance = 50;
      pullTolerance = 0.9;
    }

    std::vector< art::Ptr<recob::Hit> > trk_hitlist = trk_fm.at(trackIndicesSorted[i]);
    std::vector< art::Ptr<recob::Hit> > showerHits;

    // add track hits to shower
    for (size_t ii = 0; ii < trk_hitlist.size(); ++ii) {
      showerHits.push_back(trk_hitlist[ii]);
    } // loop over track hits

    int nShowerHits = 0;
    double showerHitPull = 0;
    bool showerCandidate = false;

    TVector3 trkStart = tracklist[i]->Vertex();
    TVector3 trkPt2; // a second point along the track
    recob::Track::Point_t trkPt2temp  = tracklist[i]->TrajectoryPoint(20).position;
    trkPt2[0] = trkPt2temp.X();
    trkPt2[1] = trkPt2temp.Y();
    trkPt2[2] = trkPt2temp.Z();

    recob::Track::Point_t trkStarttemp  = tracklist[i]->TrajectoryPoint(5).position;
    trkStart[0] = trkStarttemp.X();
    trkStart[1] = trkStarttemp.Y();
    trkStart[2] = trkStarttemp.Z();

    // track vertex
    std::map<geo::PlaneID, double> trk_tick1;
    std::map<geo::PlaneID, double> trk_wire1;

    // second track point
    std::map<geo::PlaneID, double> trk_tick2;
    std::map<geo::PlaneID, double> trk_wire2;

    for (auto iPlane = geom->begin_plane_id(); iPlane != geom->end_plane_id(); ++iPlane){
      trk_tick1[*iPlane] = detprop->ConvertXToTicks(trkStart[0], *iPlane);
      trk_wire1[*iPlane] = geom->WireCoordinate(trkStart[1], trkStart[2], *iPlane);
      trk_tick2[*iPlane] = detprop->ConvertXToTicks(trkPt2[0], *iPlane);
      trk_wire2[*iPlane] = geom->WireCoordinate(trkPt2[1], trkPt2[2], *iPlane);
    }
    
    for (size_t j = 0; j < clusterlist.size(); ++j) {

      std::vector< art::Ptr<recob::Hit> > cls_hitlist = cls_fm.at(j); 

      if (clusterlist[j]->ID() > 0 && cls_hitlist.size() > 10) continue;

      bool isGoodCluster = false; // true if the hit belongs to a cluster that should be added to the shower
      
      for (size_t jj = 0; jj < cls_hitlist.size(); ++jj) {
	// don't count hits associated with the track
	std::vector< art::Ptr<recob::Track> > hit_trklist = hit_fm.at(cls_hitlist[jj].key());
	if (hit_trklist.size()) {
	  if (hit_trklist[0]->ID() == tracklist[i]->ID()) continue;
	}

	int isGoodHit = goodHit(cls_hitlist[jj], maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2);

	if (isGoodHit == -1) {
	  isGoodCluster = false;
	  break;
	}
	else if (isGoodHit == 1) {
	  isGoodCluster = true;
	}

      } // loop over hits in cluster

      // add hits to shower
      if (isGoodCluster) {
	for (size_t jj = 0; jj < cls_hitlist.size(); ++jj) {
	  nShowerHits++;

	  int showerHitPullAdd = 0;
	  goodHit(cls_hitlist[jj], maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2, showerHitPullAdd);
	  showerHitPull += showerHitPullAdd;
	  
	  showerHits.push_back(cls_hitlist[jj]);
	} // loop over hits in cluster
      } // cluster contains hit close to track
      
    } // loop over clusters

    showerHitPull /= nShowerHits;
    if (nShowerHits > tolerance && std::abs(showerHitPull) < pullTolerance) {
      showerCandidate = true;
      
      // loop over hits to find those that aren't associated with any clusters
      for (size_t k = 0; k < hitlist.size(); ++k) {
	std::vector< art::Ptr<recob::Cluster> > hit_clslist = hitcls_fm.at(k);
	if (hit_clslist.size()) continue;
	int isGoodHit = goodHit(hitlist[k], maxDist*2, minDistVert*2, trk_wire1, trk_tick1, trk_wire2, trk_tick2);
	if (isGoodHit == 1) showerHits.push_back(hitlist[k]);
      } // loop over hits


      // TODO: Add missing tracks to shower
      /*
      TVector3 parentDir = trkPt2 - trkStart;

      // loop over tracks to see if any fall within the shower
      for (size_t k = 0; k < tracklist.size(); ++k) {
	if (k == i) continue; // don't check current track

	TVector3 trkStartOther = tracklist[k]->Vertex();
	TVector3 trkPt2Other; // a second point along the track                                                            
	recob::Track::Point_t trkPt2temp  = tracklist[k]->TrajectoryPoint(5).position;
	trkPt2Other[0] = trkPt2temp.X();
	trkPt2Other[1] = trkPt2temp.Y();
	trkPt2Other[2] = trkPt2temp.Z();
	
	TVector3 otherDir = trkPt2Other - trkStartOther;

	double ang = parentDir.Angle(otherDir);

	std::vector< art::Ptr<recob::Hit> > trk_hitlistOther = trk_fm.at(k);

	int nhits = 0;
	int ngoodhits = 0;

	maxDist = 10;
	if (tracklist[k]->Length() > 30) maxDist /= 2;

	// are the track hits close?
	for (size_t kk = 0; kk < trk_hitlistOther.size(); ++kk) {
	  nhits++;
	  int isGoodHit = goodHit(trk_hitlist[kk], maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2);

	  if (isGoodHit == -1){
	    ngoodhits = 0;
	    break;
	  }
	  else if (isGoodHit == 1) {
	    ngoodhits++;
	  }
	} // loop over track hits   
	
	double fracGood = (double)ngoodhits/nhits;

	double x1 = tracklist[i]->End().X();
	double y1 = tracklist[i]->End().Y();
	double z1 = tracklist[i]->End().Z();

	double x2 = tracklist[k]->End().X();
	double y2 = tracklist[k]->End().Y();
	double z2 = tracklist[k]->End().Z();

	double dist = sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2) );

	if (dist > 20) 
	  if (z2 - z1 < 0) continue;

	bool isGoodTrack = (ang < 0.5 && fracGood > 0.4 && nhits > 50) || (ang < 1 && fracGood > 0.6 && nhits < 20);

	if (tracklist[k]->Length() > 30) isGoodTrack = (ang < 0.5 && fracGood > 0.7);

	//	  std::cout << isGoodTrack << " " << tracklist[k]->ID() << " " << tracklist[k]->Length() << " " << ang  << " " << fracGood << " " << nhits << " " << dist << std::endl;

	if (isGoodTrack) {

	  for (size_t kk = 0; kk < trk_hitlistOther.size(); ++kk) {
	    //	    showerHits.push_back(trk_hitlistOther[kk]); 
	  } // loop over hits to add them to shower
	} // good track
	
      } // loop over tracks
      */
  
    } // decide if shower

    TVector3 dcosVtxErr;
    TVector3 xyzErr;

    std::vector<double> totalEnergy(2);
    std::vector<double> totalEnergyErr(2);

    std::vector<double> dEdx(2);
    std::vector<double> dEdxErr(2);

    if (showerCandidate) {

      std::cout << "track ID " << tracklist[i]->ID() << " " << tracklist[i]->Length() << " " << showerHitPull<< " " << nShowerHits << std::endl;
      showers->push_back(recob::Shower(trkPt2-trkStart, dcosVtxErr, tracklist[i]->Vertex(), xyzErr, totalEnergy, totalEnergyErr, dEdx, dEdx, 0, 0));
      showers->back().set_id(showers->size()-1);
      
      util::CreateAssn(*this, evt, *(showers.get()), showerHits, *(hitShowerAssociations.get()) );

      break;

    }

  } // loop over tracks

  evt.put(std::move(showers));
  evt.put(std::move(hitShowerAssociations));

} // produce

// -----------------------------------------------------
// return -1 if hit is too close to track vertex or has a wide opening angle
// return 1 if hit is close to the shower axis
// return 0 otherwise

int shower::TCShower::goodHit(art::Ptr<recob::Hit> hit, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2){

  int pull = 0;
  return goodHit(hit, maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2, pull);

} // goodHit

// -----------------------------------------------------
// return -1 if hit is too close to track vertex or has a wide opening angle
// return 1 if hit is close to the shower axis
// return 0 otherwise

int shower::TCShower::goodHit(art::Ptr<recob::Hit> hit, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2, int& pull){

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;

  double wirePitch = geom->WirePitch(hit->WireID());
  double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
  tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns                                
  double UnitsPerTick = tickToDist / wirePitch;

  double x0 = hit->WireID().Wire;
  double y0 = hit->PeakTime() * UnitsPerTick;

  double x1 = trk_wire1[hit->WireID()];
  double y1 = trk_tick1[hit->WireID()] * UnitsPerTick;

  double x2 = trk_wire2[hit->WireID()];
  double y2 = trk_tick2[hit->WireID()] * UnitsPerTick;

  double distToVert = std::sqrt( pow(x0 - x1, 2) + pow(y0 - y1, 2) );
  if (distToVert < minDistVert) return -1;

  // exclude cluster if it's "behind" the vertex                                                                      
  double a = std::sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) );
  double b = std::sqrt( pow(x0 - x1, 2) + pow(y0 - y1, 2) );
  double c = std::sqrt( pow(x0 - x2, 2) + pow(y0 - y2, 2) );
  double costheta = -( pow(c,2) - pow(a,2) - pow(b,2) ) / (2 * a * b);
  if (costheta < 0) return -1;

  double dist = std::abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/std::sqrt( pow((y2-y1), 2) + pow((x2-x1), 2) );

  if (dist < maxDist) {
    if ( ( (y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1) > 0) pull = 1;
    else pull = -1; 
    return 1;
  }

  return 0;

} // goodHit

// -----------------------------------------------------
// Performs a linear regression on the first n trajectory points
// of a track and a point along the line.

std::vector<TVector3> shower::TCShower::getSecondPoint(art::Ptr<recob::Track> thistrack) {

  size_t starttrajpoint = 0;
  size_t ntrajpoints = 10 + starttrajpoint;

  TMatrix Y = TMatrix((int) ntrajpoints, 1);
  TMatrix X = TMatrix((int) ntrajpoints, 3);

  for (size_t i = starttrajpoint; i < ntrajpoints; ++i) {

    Y(i, 0) = thistrack->TrajectoryPoint(i).position.Z();
    X(i, 0) = 1;
    X(i, 1) = thistrack->TrajectoryPoint(i).position.X();
    X(i, 2) = thistrack->TrajectoryPoint(i).position.Y();

  } // loop over trajectory points

  TMatrix XT = X;
  XT.T();

  TMatrix XTY = XT * Y;
  TMatrix XTX = XT * X;
  TMatrix Z = XTX.Invert() * XTY;

  TVector3 trkPt1;
  TVector3 trkPt2;

  double x = thistrack->TrajectoryPoint(10).position.X();
  double y = thistrack->TrajectoryPoint(10).position.Y();
  double z = Z(0,0) + x * Z(1,0) + y * Z(2,0); 

  trkPt2[0] = x;
  trkPt2[1] = y;
  trkPt2[2] = z;

  std::vector<TVector3> trkPts;

  trkPts.push_back(trkPt1);
  trkPts.push_back(trkPt2);

  return trkPts;

} // getSecondPoint

// -----------------------------------------------------

void shower::TCShower::beginJob() {

  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(shower::TCShower)

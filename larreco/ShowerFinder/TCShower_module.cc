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
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TCShowerAlg.h"

#include "TH1F.h"

#include <memory>

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
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  int getShowersWithSlices(art::Event & evt, art::Ptr<recob::Slice> thisslice);
  int getShowersWithoutSlices(art::Event & evt);

  shower::TCShowerAlg fTCAlg;

  std::string fClusterModuleLabel;
  std::string fTrackModuleLabel;
  std::string fHitModuleLabel;
  std::string fSliceModuleLabel;
  std::string fVertexModuleLabel;
  std::string fCalorimetryModuleLabel;

};

// -----------------------------------------------------

shower::TCShower::TCShower(fhicl::ParameterSet const & pset) : 
  fTCAlg(pset.get< fhicl::ParameterSet >("TCAlg") ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel", "trajcluster" ) ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel", "trajclusterKalmanTrack" ) ),
  fHitModuleLabel           (pset.get< std::string >("HitModuleLabel", "trajcluster" ) ),
  fSliceModuleLabel         (pset.get< std::string >("SliceModuleLabel", "dbcluster3d" ) ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel", "trajcluster" ) ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel")  ) {

  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Slice, recob::Shower> >();
}

// -----------------------------------------------------

void shower::TCShower::produce(art::Event & evt) {
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitShowerAssociations(new art::Assns<recob::Shower, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Slice, recob::Shower> > sliceShowerAssociations(new art::Assns<recob::Slice, recob::Shower>);

  // slices
  art::Handle< std::vector<recob::Slice> > sliceListHandle;
  std::vector<art::Ptr<recob::Slice> > slicelist;
  if (evt.getByLabel(fSliceModuleLabel,sliceListHandle))
    art::fill_ptr_vector(slicelist, sliceListHandle);

  int foundShower = -1;

  if (slicelist.size()) { // use slices
    for (size_t i = 0; i < slicelist.size(); ++i) {
      std::cout << "---------- slice " << i << " ----------" << std::endl;

      foundShower = getShowersWithSlices(evt, slicelist[i]);

      if (foundShower) {
	std::cout << "FOUND SHOWER " << foundShower << std::endl;
	std::cout << "shower hits " << fTCAlg.showerHits.size() << std::endl;

	showers->push_back(recob::Shower(fTCAlg.shwDir, fTCAlg.dcosVtxErr, fTCAlg.shwvtx, fTCAlg.xyzErr, fTCAlg.totalEnergy, fTCAlg.totalEnergyErr, fTCAlg.dEdx, fTCAlg.dEdxErr, fTCAlg.bestplane, 0));
	showers->back().set_id(showers->size()-1);

	util::CreateAssn(*this, evt, *(showers.get()), fTCAlg.showerHits, *(hitShowerAssociations.get()) );
	util::CreateAssn(*this, evt, *showers, slicelist[i], *sliceShowerAssociations );
      }

    } // loop through slices
  } // with slices
  else {
    foundShower = getShowersWithoutSlices(evt);

    if (foundShower) {
      showers->push_back(recob::Shower(fTCAlg.shwDir, fTCAlg.dcosVtxErr, fTCAlg.shwvtx, fTCAlg.xyzErr, fTCAlg.totalEnergy, fTCAlg.totalEnergyErr, fTCAlg.dEdx, fTCAlg.dEdxErr, fTCAlg.bestplane, 0));
      showers->back().set_id(showers->size()-1);
    
      util::CreateAssn(*this, evt, *(showers.get()), fTCAlg.showerHits, *(hitShowerAssociations.get()) );

    }

  } // no slices

  evt.put(std::move(showers));
  evt.put(std::move(hitShowerAssociations));
  evt.put(std::move(sliceShowerAssociations));

} // produce

// -----------------------------------------------------
int shower::TCShower::getShowersWithSlices(art::Event & evt, art::Ptr<recob::Slice> thisslice) {
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitModuleLabel,hitListHandle);

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  /*
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel, trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

  for (size_t i = 0; i < tracklist.size(); ++i) {
    std::cout << "track id " << tracklist[i]->ID() << std::endl;
  }
  */
  art::Handle< std::vector<recob::Slice> > sliceListHandle;
  evt.getByLabel(fSliceModuleLabel,sliceListHandle);

  art::Handle< std::vector<recob::Vertex> > vtxListHandle;
  evt.getByLabel(fVertexModuleLabel,vtxListHandle);

  //  art::Handle< std::vector<recob::EndPoint2D> > vx2ListHandle;
  //  evt.getByLabel(fVertexModuleLabel, vx2ListHandle);

  art::Handle< std::vector<recob::PFParticle> > pfpListHandle;
  evt.getByLabel(fHitModuleLabel,pfpListHandle);

  art::FindManyP<recob::Hit> hitslice_fm(sliceListHandle, evt, fHitModuleLabel);
  art::FindManyP<recob::PFParticle> pfpslice_fm(sliceListHandle, evt, fHitModuleLabel);
  art::FindManyP<recob::Cluster> clsslice_fm(sliceListHandle, evt, fHitModuleLabel);
  art::FindManyP<recob::Cluster> clspfp_fm(pfpListHandle, evt, fHitModuleLabel);
  art::FindManyP<recob::Vertex> vtxpfp_fm(pfpListHandle, evt, fVertexModuleLabel);
  //  art::FindManyP<recob::EndPoint2D> vx2cls_fm(clusterListHandle, evt, fClusterModuleLabel);

  std::vector<art::Ptr<recob::Hit> > hitlist;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  std::vector<art::Ptr<recob::Vertex> > vertexlist;
  std::vector<art::Ptr<recob::EndPoint2D> > vx2list;
  
  // get all hits with hit-slice association
  hitlist = hitslice_fm.at(thisslice.key());

  // get all clusters with cluster-slice association
  clusterlist = clsslice_fm.at(thisslice.key());
  /*
  for (size_t i = 0; i < clusterlist.size(); ++i) {
    std::vector<art::Ptr<recob::EndPoint2D> > eplist = vx2cls_fm.at(clusterlist[i].key());
    
    for (size_t j = 0; j < eplist.size(); ++j) {
      bool addToList = true;
      for (size_t k = 0; k < vx2list.size(); ++k) {
	if (eplist[j]->ID() == vx2list[k]->ID() ) {
	  addToList = false;
	  break;
	}
      }  

      if (addToList) vx2list.push_back(eplist[j]);

    }
  } // loop through clusterlist
  */
  std::vector<art::Ptr<recob::PFParticle> > pfplist = pfpslice_fm.at(thisslice.key());

  for (size_t i = 0; i < pfplist.size(); ++i) {
    //    std::vector<art::Ptr<recob::Track> > thistracklist = trkpfp_fm.at(pfplist[i].key());    
    std::vector<art::Ptr<recob::Vertex> > thisvtxlist = vtxpfp_fm.at(pfplist[i].key());    
    // get all verticies with slice-pfparticle, pfparticle-vertex
     for (size_t j = 0; j < thisvtxlist.size(); ++j) {
      vertexlist.push_back(thisvtxlist[j]);
    } // loop through tracks
   } // loop through pfparticles
 

  // get associations
  art::FindManyP<recob::Hit> cls_fm(clusterListHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::PFParticle> hit_fm(hitListHandle, evt, fHitModuleLabel);
  art::FindManyP<recob::Cluster> hitcls_fm(hitListHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Track> trkpfp_fm(pfpListHandle, evt, fTrackModuleLabel);

  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

  return fTCAlg.makeShowers(pfplist, vertexlist, clusterlist, hitlist, cls_fm, clspfp_fm, vtxpfp_fm, hit_fm, hitcls_fm, trkpfp_fm, fmcal);

}

// -----------------------------------------------------
int shower::TCShower::getShowersWithoutSlices(art::Event & evt) {

  // pfparticles
  art::Handle< std::vector<recob::PFParticle> > pfpListHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfplist;
  if (evt.getByLabel(fHitModuleLabel,pfpListHandle))
    art::fill_ptr_vector(pfplist, pfpListHandle);

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::Handle< std::vector<recob::Vertex> > vtxListHandle;
  std::vector<art::Ptr<recob::Vertex> > vertexlist;
  if (evt.getByLabel(fVertexModuleLabel,vtxListHandle))
    art::fill_ptr_vector(vertexlist, vtxListHandle);

  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);

  // get associations                                                                                                      
  art::FindManyP<recob::Cluster> clspfp_fm(pfpListHandle, evt, fHitModuleLabel);
  art::FindManyP<recob::Vertex> vtxpfp_fm(pfpListHandle, evt, fVertexModuleLabel);
  art::FindManyP<recob::Hit> cls_fm(clusterListHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::PFParticle> hit_fm(hitListHandle, evt, fHitModuleLabel);
  art::FindManyP<recob::Cluster> hitcls_fm(hitListHandle, evt, fClusterModuleLabel);

  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

  //  return fTCAlg.makeShowers(pfplist, vertexlist, clusterlist, hitlist, cls_fm, clspfp_fm, vtxpfp_fm, hit_fm, hitcls_fm, fmcal);

  return 0;

}

// -----------------------------------------------------
void shower::TCShower::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
}

// -----------------------------------------------------

DEFINE_ART_MODULE(shower::TCShower)

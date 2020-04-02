////////////////////////////////////////////////////////////////////////
// Class:       SCECorrection
// Plugin Type: producer (art v3_04_00)
// File:        SCECorrection_module.cc
//
// Generated at Sun Mar 22 09:23:33 2020 by Edward Tyley using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

namespace sce {
  class SCECorrection;
}


class sce::SCECorrection : public art::EDProducer {
  public:
    explicit SCECorrection(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    SCECorrection(SCECorrection const&) = delete;
    SCECorrection(SCECorrection&&) = delete;
    SCECorrection& operator=(SCECorrection const&) = delete;
    SCECorrection& operator=(SCECorrection&&) = delete;

    // Required functions.
    void produce(art::Event& evt) override;

  private:

    // Declare member data here.
    detinfo::DetectorProperties const* fDetProp;
    geo::GeometryCore const* fGeom;
    spacecharge::SpaceCharge const* fSCE;

    bool fCorrectNoT0Tag;
    bool fCorrectSCE;

    std::string fPFPLabel;
    std::string fTrackLabel;
    std::vector<std::string> fT0Labels;
    std::vector<bool> fT0LabelsCorrectT0;

    geo::Vector_t applyT0Shift(const double& t0, const geo::TPCID tpcId) const;

    std::map<art::Ptr<anab::T0>, bool> getSliceT0s(
        const art::Event& evt,
        const std::vector<art::Ptr<recob::PFParticle> > slicePFPs,
        const art::Handle<std::vector<recob::PFParticle> >& pfpHandle,
        const art::Handle<std::vector<recob::Track> >& trackHandle,
        const art::FindManyP<recob::Track> fmPFPTrack) const;

    std::pair<art::Ptr<anab::T0>, bool> getSliceBestT0(
        std::map<art::Ptr<anab::T0>, bool> sliceT0CorrectMap);
};

sce::SCECorrection::SCECorrection(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
  fGeom(lar::providerFrom<geo::Geometry>()),
  fSCE(lar::providerFrom<spacecharge::SpaceChargeService>())
{
  fCorrectNoT0Tag     = p.get<bool>("CorrectNoT0Tag");
  fCorrectSCE         = p.get<bool>("CorrectSCE");
  fPFPLabel           = p.get<std::string>("PFPLabel");
  fTrackLabel         = p.get<std::string>("TrackLabel");
  fT0Labels           = p.get<std::vector<std::string> >("T0Labels");
  fT0LabelsCorrectT0  = p.get<std::vector<bool> >("T0LabelsCorrectT0");

  produces<std::vector<anab::T0> >();
  produces<std::vector<recob::Slice> >();
  produces<std::vector<recob::PFParticle> >();
  produces<std::vector<recob::SpacePoint> >();
  produces<std::vector<recob::Cluster> >();
  produces<std::vector<recob::Vertex> >();
  produces<std::vector<larpandoraobj::PFParticleMetadata> >();

  // produces<art::Assns<anab::T0, recob::Slice> >();
  produces<art::Assns<anab::T0, recob::PFParticle> >();
  produces<art::Assns<recob::Slice, recob::Hit>>();

  produces<art::Assns<recob::PFParticle, recob::Slice> >();
  produces<art::Assns<recob::PFParticle, recob::SpacePoint> >();
  produces<art::Assns<recob::PFParticle, recob::Vertex> >();
  produces<art::Assns<recob::PFParticle, recob::Cluster> >();
  produces<art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> >();
  produces<art::Assns<recob::SpacePoint, recob::Hit> >();
  produces<art::Assns<recob::Cluster, recob::Hit> >();
}

void sce::SCECorrection::produce(art::Event& evt)
{
  // Implementation of required member function here.
  // auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto t0Collection      = std::make_unique<std::vector<anab::T0> >();
  auto pfpCollection     = std::make_unique<std::vector<recob::PFParticle> >();
  auto clusterCollection = std::make_unique<std::vector<recob::Cluster> >();
  auto spCollection      = std::make_unique<std::vector<recob::SpacePoint> >();
  auto vtxCollection     = std::make_unique<std::vector<recob::Vertex> >();
  auto sliceCollection   = std::make_unique<std::vector<recob::Slice> >();
  auto pfpMetaCollection = std::make_unique<std::vector<larpandoraobj::PFParticleMetadata> > ();

  // auto t0SliceAssn       = std::make_unique<art::Assns<anab::T0, recob::Slice> >();
  auto t0PFPAssn         = std::make_unique<art::Assns<anab::T0, recob::PFParticle> >();
  auto slcHitAssn        = std::make_unique<art::Assns<recob::Slice, recob::Hit> >();
  auto pfpSliceAssn      = std::make_unique<art::Assns<recob::PFParticle, recob::Slice> >();
  auto pfpVtxAssn        = std::make_unique<art::Assns<recob::PFParticle, recob::Vertex> >();
  auto pfpSPAssn         = std::make_unique<art::Assns<recob::PFParticle, recob::SpacePoint> >();
  auto pfpClusterAssn    = std::make_unique<art::Assns<recob::PFParticle, recob::Cluster> >();
  auto pfpMetaAssn       = std::make_unique<art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> >();
  auto spHitAssn         = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit> >();
  auto clusterHitAssn    = std::make_unique<art::Assns<recob::Cluster, recob::Hit> >();

  art::PtrMaker<anab::T0>                          t0PtrMaker{evt};
  art::PtrMaker<recob::PFParticle>                 pfpPtrMaker{evt};
  art::PtrMaker<recob::Cluster>                    clusterPtrMaker{evt};
  art::PtrMaker<recob::Vertex>                     vtxPtrMaker{evt};
  art::PtrMaker<recob::Slice>                      slicePtrMaker{evt};
  art::PtrMaker<recob::SpacePoint>                 spPtrMaker{evt};
  art::PtrMaker<larpandoraobj::PFParticleMetadata> pfpMetaPtrMaker{evt};

  // Get all the slices in the event
  art::Handle<std::vector<recob::Slice> > sliceHandle;
  std::vector<art::Ptr<recob::Slice> > allSlices;
  if (evt.getByLabel(fPFPLabel, sliceHandle))
    art::fill_ptr_vector(allSlices, sliceHandle);

  // Get all the Clusters in the event
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  std::vector<art::Ptr<recob::Cluster> > allClusters;
  if (evt.getByLabel(fPFPLabel, clusterHandle))
    art::fill_ptr_vector(allClusters, clusterHandle);

  // Get all the SpacePoints in the event
  art::Handle<std::vector<recob::SpacePoint> > spHandle;
  std::vector<art::Ptr<recob::SpacePoint> > allSpacePoints;
  if (evt.getByLabel(fPFPLabel, spHandle))
    art::fill_ptr_vector(allSpacePoints, spHandle);

  // Get all the PFParticles in the event
  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> > allPFParticles;
  if (evt.getByLabel(fPFPLabel, pfpHandle))
    art::fill_ptr_vector(allPFParticles, pfpHandle);

  // Get all the Tracks in the event
  art::Handle<std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > allTracks;
  if (evt.getByLabel(fTrackLabel, trackHandle))
    art::fill_ptr_vector(allTracks, trackHandle);

  art::FindManyP<recob::PFParticle>                 fmSlicePFP(sliceHandle, evt, fPFPLabel);
  art::FindManyP<recob::Track>                      fmPFPTrack(pfpHandle, evt, fTrackLabel);
  art::FindManyP<recob::SpacePoint>                 fmPFPSP(pfpHandle, evt, fPFPLabel);
  art::FindManyP<recob::Cluster>                    fmPFPCluster(pfpHandle, evt, fPFPLabel);
  art::FindManyP<recob::Vertex>                     fmPFPVertex(pfpHandle, evt, fPFPLabel);
  art::FindManyP<recob::Hit>                        fmClusterHit(clusterHandle, evt, fPFPLabel);
  art::FindManyP<recob::Hit>                        fmSPHit(spHandle, evt, fPFPLabel);
  art::FindManyP<recob::Hit>                        fmSlcHit(sliceHandle, evt, fPFPLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta(pfpHandle, evt, fPFPLabel);

  if (!fmPFPTrack.isValid()){
    std::cout<<"PFPTrack assn not valid"<<std::endl;
    return;
  }

  // For each slice, get all the PFPs and tracks and check for T0 tags
  // std::cout<<"Test: Slices: "<<allSlices.size()<<std::endl;
  for (auto const& slice: allSlices){

    //Cretae a new slice
    recob::Slice newSlice(*slice);
    sliceCollection->push_back(newSlice);
    art::Ptr<recob::Slice> newSlicePtr = slicePtrMaker(sliceCollection->size()-1);

    // Get the pfps and hits associated to the slice
    std::vector<art::Ptr<recob::PFParticle> > slicePFPs = fmSlicePFP.at(slice.key());

    std::map<art::Ptr<anab::T0>, bool> sliceT0CorrectMap = getSliceT0s(
        evt, slicePFPs, pfpHandle, trackHandle, fmPFPTrack);

    std::pair<art::Ptr<anab::T0>, bool> sliceT0CorrectPair = getSliceBestT0(sliceT0CorrectMap);

    if (sliceT0CorrectPair.first.isNull() && !fCorrectNoT0Tag){
      continue;
    }

    art::Ptr<anab::T0> newT0Ptr;
    double t0Offset(0);
    if (!sliceT0CorrectPair.first.isNull()){
      // Calculate the shift we need to apply for the t0
      t0Offset = fDetProp->DriftVelocity() * sliceT0CorrectPair.first->Time() / 1e3;
      // Create a new T0
      t0Collection->push_back(*sliceT0CorrectPair.first);
      newT0Ptr = t0PtrMaker(t0Collection->size()-1);
      // t0SliceAssn->addSingle(newT0Ptr, newSlicePtr);
    }

    // Make an association with the new slice and the old hits
    const std::vector<art::Ptr<recob::Hit> > sliceHits = fmSlcHit.at(slice.key());
    for (const art::Ptr<recob::Hit> &hitPtr: sliceHits) {
      slcHitAssn->addSingle(newSlicePtr, hitPtr);
    }

    // Correct all PFPs in the slice
    for (auto const& pfp: slicePFPs){

      // Create new PFPs and associate them to the slice
      recob::PFParticle newPFP(*pfp);
      pfpCollection->push_back(newPFP);
      art::Ptr<recob::PFParticle> newPFPPtr = pfpPtrMaker(pfpCollection->size()-1);
      pfpSliceAssn->addSingle(newPFPPtr, newSlicePtr);

      if (!newT0Ptr.isNull()){
        t0PFPAssn->addSingle(newT0Ptr, newPFPPtr);
      }

      // Get the vertex associated to the PFP
      std::vector<art::Ptr<recob::Vertex> > pfpVertexs = fmPFPVertex.at(pfp.key());
      std::vector<art::Ptr<recob::SpacePoint> > pfpSPs = fmPFPSP.at(pfp.key());
      for (auto const& pfpVertex: pfpVertexs){

        geo::Point_t vtxPos(pfpVertex->position());
        //Find the closest SP to the vertex
        // If the PFP has no space points, look in the whole event
        std::vector<art::Ptr<recob::SpacePoint> vtxSPs = pfpSPs.size() ? pfpSPs: allSpacePoints;

        double minVtxSPDist = std::numeric_limits<double>::max();
        art::Ptr<recob::SpacePoint> spPtr;
        for (auto const& sp: vtxSPs){
          geo::Point_t spPos{sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2]};
          geo::Vector_t vtxSPDiff = vtxPos - spPos;
          if (vtxSPDiff.Mag2() < minVtxSPDist){
            spPtr = sp;
            minVtxSPDist = vtxSPDiff.Mag2();
          }
        }

        if (spPtr.isNull())
          continue;

        // Get the hit and TPC Id associated to closest SP
        art::Ptr<recob::Hit> spHitPtr = fmSPHit.at(spPtr.key()).front();
        geo::TPCID tpcId = spHitPtr->WireID().asTPCID();

        if (!sliceT0CorrectPair.first.isNull() && sliceT0CorrectPair.second){
          geo::Vector_t posOffset = applyT0Shift(t0Offset, tpcId);
          vtxPos += posOffset;
        }

        if(fCorrectSCE && fSCE->EnableCalSpatialSCE()){
          geo::Vector_t posOffset = fSCE->GetCalPosOffsets(vtxPos, tpcId.TPC);
          vtxPos += posOffset;
        }

        // Create a new vertex and associate it to the PFP
        recob::Vertex newVtx(vtxPos, pfpVertex->covariance(), pfpVertex->chi2(),
            pfpVertex->ndof(), pfpVertex->ID());
        vtxCollection->push_back(newVtx);
        art::Ptr<recob::Vertex> newVtxPtr = vtxPtrMaker(vtxCollection->size()-1);
        pfpVtxAssn->addSingle(newPFPPtr, newVtxPtr);
      }

      for (auto const& sp: pfpSPs){

        //Get the spacepoint position in a nicer form
        geo::Point_t spPos{sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2]};

        // Get the hit so we know what TPC the sp was in
        // N.B. We can't use SP position to infer the TPC as it could be
        // shifted into another TPC
        art::Ptr<recob::Hit> spHitPtr = fmSPHit.at(sp.key()).front();
        geo::TPCID tpcId = spHitPtr->WireID().asTPCID();

        if (!sliceT0CorrectPair.first.isNull() && sliceT0CorrectPair.second){
          geo::Vector_t posOffset = applyT0Shift(t0Offset, tpcId);
          spPos += posOffset;
        }

        if(fCorrectSCE && fSCE->EnableCalSpatialSCE()){
          geo::Vector_t posOffset = fSCE->GetCalPosOffsets(spPos, tpcId.TPC);
          spPos += posOffset;
        }

        // Create new spacepoint and associate it to the pfp and hit
        Double32_t spXYZ[3] = {spPos.X(), spPos.Y(), spPos.Z()};
        recob::SpacePoint correctedSP(spXYZ, sp->ErrXYZ(), sp->Chisq(), sp->ID());

        spCollection->push_back(correctedSP);
        art::Ptr<recob::SpacePoint> spPtr = spPtrMaker(spCollection->size()-1);
        pfpSPAssn->addSingle(newPFPPtr, spPtr);
        spHitAssn->addSingle(spPtr, spHitPtr);
      } // pspSPs

      // Create new clusters and associations
      std::vector<art::Ptr<recob::Cluster> > pfpClusters = fmPFPCluster.at(pfp.key());
      for (auto const& pfpCluster: pfpClusters){
        recob::Cluster newCluster(*pfpCluster);
        clusterCollection->push_back(newCluster);
        art::Ptr<recob::Cluster> newClusterPtr = clusterPtrMaker(clusterCollection->size()-1);

        std::vector<art::Ptr<recob::Hit> > clusterHits = fmClusterHit.at(pfpCluster.key());
        pfpClusterAssn->addSingle(newPFPPtr, newClusterPtr);
        for (auto const& clusterHit: clusterHits){
          clusterHitAssn->addSingle(newClusterPtr, clusterHit);
        }
      }

      // Create new PFParticle Metadata objects and associations
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetas = fmPFPMeta.at(pfp.key());
      for (const art::Ptr<larpandoraobj::PFParticleMetadata> &pfpMeta: pfpMetas) {
        larpandoraobj::PFParticleMetadata newPFPMeta(*pfpMeta);
        pfpMetaCollection->push_back(newPFPMeta);
        art::Ptr<larpandoraobj::PFParticleMetadata> newPFPMetaPtr =
          pfpMetaPtrMaker(pfpMetaCollection->size()-1);
        pfpMetaAssn->addSingle(newPFPPtr, newPFPMetaPtr);
      }
    } // slicePFPs
  } // slice

  // std::cout<<"Test: SPs: "<<spCollection->size()<<std::endl;

  // Order the PFPs so the id matched the position in the vector
  std::sort(pfpCollection->begin(), pfpCollection->end(),
      [](recob::PFParticle lhs, recob::PFParticle rhs) {return lhs.Self() < rhs.Self();});

  // Put all the things we just produced into the event
  evt.put(std::move(t0Collection));
  evt.put(std::move(sliceCollection));
  evt.put(std::move(clusterCollection));
  evt.put(std::move(pfpCollection));
  evt.put(std::move(spCollection));
  evt.put(std::move(vtxCollection));
  evt.put(std::move(pfpMetaCollection));

  evt.put(std::move(t0PFPAssn));
  // evt.put(std::move(t0SliceAssn));
  evt.put(std::move(slcHitAssn));
  evt.put(std::move(pfpSPAssn));
  evt.put(std::move(spHitAssn));
  evt.put(std::move(pfpVtxAssn));
  evt.put(std::move(pfpSliceAssn));
  evt.put(std::move(pfpClusterAssn));
  evt.put(std::move(clusterHitAssn));
  evt.put(std::move(pfpMetaAssn));
}

geo::Vector_t sce::SCECorrection::applyT0Shift(const double& t0Offset, const geo::TPCID tpcId) const{

  const geo::TPCGeo& tpcGeo = fGeom->GetElement(tpcId);
  int driftDirection = tpcGeo.DetectDriftDirection();

  switch (std::abs(driftDirection)) {
    case 1: return geo::Vector_t{t0Offset*driftDirection, 0, 0};
    case 2: return geo::Vector_t{0, t0Offset*driftDirection, 0};
    case 3: return geo::Vector_t{0, 0, t0Offset*driftDirection};
    default: return geo::Vector_t{0,0,0}; //TODO: make an exception
  }
}

std::map<art::Ptr<anab::T0>, bool> sce::SCECorrection::getSliceT0s (
    const art::Event& evt,
    const std::vector<art::Ptr<recob::PFParticle> > slicePFPs,
    const art::Handle<std::vector<recob::PFParticle> >& pfpHandle,
    const art::Handle<std::vector<recob::Track> >& trackHandle,
    const art::FindManyP<recob::Track> fmPFPTrack) const{

  std::map<art::Ptr<anab::T0>, bool> pfpT0CorrectMap;
  // Loop over all of the PFPs in the slice
  for (auto const& pfp: slicePFPs){

    // Loop over all of the T0 labels
    // We will take the first label to have a T0, so the order matters
    for (unsigned int i=0; i<fT0Labels.size(); i++){
      std::string t0Label = fT0Labels.at(i);

      // Get the T0
      art::FindManyP<anab::T0> fmPFPT0(pfpHandle, evt, t0Label);
      if (fmPFPT0.isValid()){
        std::vector<art::Ptr<anab::T0> > pfpT0s = fmPFPT0.at(pfp.key());
        if (pfpT0s.size()==1){
          pfpT0CorrectMap[pfpT0s.front()] = fT0LabelsCorrectT0.at(i);
          break;
        }
      }
      // If not, Check the track associated to the PFP
      std::vector<art::Ptr<recob::Track> > pfpTracks = fmPFPTrack.at(pfp.key());
      if (pfpTracks.size()!=1){
        continue;
      }
      art::Ptr<recob::Track> pfpTrack = pfpTracks.front();

      // Check if the track has a T0
      art::FindManyP<anab::T0> fmTrackT0(trackHandle, evt, t0Label);
      if (fmTrackT0.isValid()){
        std::vector<art::Ptr<anab::T0> > trackT0s = fmTrackT0.at(pfpTrack.key());
        if (trackT0s.size()==1){
          pfpT0CorrectMap[trackT0s.front()] = fT0LabelsCorrectT0.at(i);
          break;
        }
      }
    } // fT0Labels
  } // slicePFPs
  return pfpT0CorrectMap;
}

std::pair<art::Ptr<anab::T0>, bool> sce::SCECorrection::getSliceBestT0(
    std::map<art::Ptr<anab::T0>, bool> sliceT0CorrectMap){

  if (!sliceT0CorrectMap.size()){
    return std::pair<art::Ptr<anab::T0>, bool>();
  }

  double minT0 = std::numeric_limits<double>::max();
  std::pair<art::Ptr<anab::T0>, bool> sliceT0CorrectPair;
  for (auto const& sliceT0CorrectIter: sliceT0CorrectMap) {
    double t0Time = abs(sliceT0CorrectIter.first->Time());
    if (t0Time < minT0){
      minT0 = t0Time;
      sliceT0CorrectPair = sliceT0CorrectIter;
    }
  }
  return sliceT0CorrectPair;
}

DEFINE_ART_MODULE(sce::SCECorrection)

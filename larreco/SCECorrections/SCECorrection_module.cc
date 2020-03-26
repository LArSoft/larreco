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
    void produce(art::Event& e) override;

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

    geo::Vector_t correctT0(const geo::Point_t& point, const double& t0,
        const geo::TPCID tpcId) const;
};

sce::SCECorrection::SCECorrection(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
  fGeom(lar::providerFrom<geo::Geometry>()),
  fSCE(lar::providerFrom<spacecharge::SpaceChargeService>())
{
  fCorrectNoT0Tag    = p.get<bool>("CorrectNoT0Tag");
  fCorrectSCE        = p.get<bool>("CorrectSCE");
  fPFPLabel          = p.get<std::string>("PFPLabel");
  fTrackLabel        = p.get<std::string>("TrackLabel");
  fT0Labels          = p.get<std::vector<std::string> >("T0Labels");
  fT0LabelsCorrectT0 = p.get<std::vector<bool> >("T0LabelsCorrectT0");

  produces<std::vector<anab::T0> >();
  produces<std::vector<recob::Slice> >();
  produces<std::vector<recob::PFParticle> >();
  produces<std::vector<recob::SpacePoint> >();
  produces<std::vector<recob::Cluster> >();
  produces<std::vector<recob::Vertex> >();

  produces<art::Assns<anab::T0, recob::Slice> >();
  produces<art::Assns<anab::T0, recob::PFParticle> >();
  produces<art::Assns<recob::PFParticle, recob::Slice> >();
  produces<art::Assns<recob::PFParticle, recob::SpacePoint> >();
  produces<art::Assns<recob::PFParticle, recob::Vertex> >();
  produces<art::Assns<recob::PFParticle, recob::Cluster> >();
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

  auto t0SliceAssn       = std::make_unique<art::Assns<anab::T0, recob::Slice> >();
  auto t0PFPAssn         = std::make_unique<art::Assns<anab::T0, recob::PFParticle> >();
  auto pfpSliceAssn      = std::make_unique<art::Assns<recob::PFParticle, recob::Slice> >();
  auto pfpVtxAssn        = std::make_unique<art::Assns<recob::PFParticle, recob::Vertex> >();
  auto pfpSPAssn         = std::make_unique<art::Assns<recob::PFParticle, recob::SpacePoint> >();
  auto pfpClusterAssn    = std::make_unique<art::Assns<recob::PFParticle, recob::Cluster> >();
  auto spHitAssn         = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit> >();
  auto clusterHitAssn    = std::make_unique<art::Assns<recob::Cluster, recob::Hit> >();

  art::PtrMaker<anab::T0>          t0PtrMaker{evt};
  art::PtrMaker<recob::PFParticle> pfpPtrMaker{evt};
  art::PtrMaker<recob::Cluster>    clusterPtrMaker{evt};
  art::PtrMaker<recob::Vertex>     vtxPtrMaker{evt};
  art::PtrMaker<recob::Slice>      slicePtrMaker{evt};
  art::PtrMaker<recob::SpacePoint> spPtrMaker{evt};

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

  art::FindManyP<recob::PFParticle> fmSlicePFP(sliceHandle, evt, fPFPLabel);
  art::FindManyP<recob::Track>      fmPFPTrack(pfpHandle, evt, fTrackLabel);
  art::FindManyP<recob::SpacePoint> fmPFPSP(pfpHandle, evt, fPFPLabel);
  art::FindManyP<recob::Cluster>    fmPFPCluster(pfpHandle, evt, fPFPLabel);
  art::FindManyP<recob::Vertex>     fmPFPVertex(pfpHandle, evt, fPFPLabel);
  art::FindManyP<recob::Hit>        fmClusterHit(clusterHandle, evt, fPFPLabel);
  art::FindManyP<recob::Hit>        fmSPHit(spHandle, evt, fPFPLabel);

  if (!fmPFPTrack.isValid()){
    std::cout<<"PFPTrack assn not valid"<<std::endl;
    return;
  }

  // For each slice, get all the PFPs and tracks and check for T0 tags
  // std::cout<<"Test: Slices: "<<allSlices.size()<<std::endl;
  for (auto const& slice: allSlices){

    recob::Slice newSlice(*slice);
    sliceCollection->push_back(newSlice);
    art::Ptr<recob::Slice> newSlicePtr = slicePtrMaker(sliceCollection->size()-1);

    // Set the default T0 to 0 (Trigger time)
    std::vector<art::Ptr< anab::T0> > sliceT0s;
    std::vector<unsigned int> sliceT0Labels;
    std::vector<art::Ptr<recob::PFParticle> > slicePFPs = fmSlicePFP.at(slice.key());

    // Loop over all of the PFPs in the slice
    // std::cout<<"Test: Slice PPFs: "<<slicePFPs.size()<<std::endl;
    for (auto const& pfp: slicePFPs){

      // Loop over all of the T0 labels
      // We will take the first label to have a T0, so the order matters
      for (unsigned int i=0; i<fT0Labels.size(); i++){
        std::string t0Label = fT0Labels.at(i);

        art::FindManyP<anab::T0> fmPFPT0(pfpHandle, evt, t0Label);
        // Check if the PFP has a T0
        if (fmPFPT0.isValid()){
          std::vector<art::Ptr<anab::T0> > pfpT0s = fmPFPT0.at(pfp.key());
          if (pfpT0s.size()==1){
            sliceT0s.push_back(pfpT0s.front());
            sliceT0Labels.push_back(i);
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
            sliceT0s.push_back(trackT0s.front());
            sliceT0Labels.push_back(i);
            break;
          }
        }
      } // fT0Labels
    } // slicePFPs

    // If there are no T0 tags in the slice and fCorrectNoT0Tag is false, skip the slice
    // If there is 1 t0 tagged PFP, use that t0
    // If there is >1 t0 tagges PFP, take the t0 closest to the trigger (0)
    // std::cout<<"Test: Slice T0s: "<<pfpT0Times.size()<<std::endl;
    art::Ptr<anab::T0> sliceT0;
    unsigned int sliceT0Label(0);
    if (!sliceT0s.size() && !fCorrectNoT0Tag){
      continue;
    } else if (sliceT0s.size()==1){
      sliceT0 = sliceT0s.front();
      sliceT0Label = sliceT0Labels.front();
    } else {
      double sliceT0Time = std::numeric_limits<double>::max();
      for (unsigned int i=0; i<sliceT0s.size(); i++){
        art::Ptr<anab::T0> sliceT0Ptr = sliceT0s.at(i);
        if (abs(sliceT0Ptr->Time())<abs(sliceT0Time)){
          sliceT0 = sliceT0Ptr;
          sliceT0Label = sliceT0Labels.at(i);
          sliceT0Time =  sliceT0Ptr->Time();
        }
      } // sliceT0Times
    }

    bool t0Set = !sliceT0.isNull();
    art::Ptr<anab::T0> sliceT0Ptr;
    if (t0Set){
      t0Collection->push_back(*sliceT0);
      sliceT0Ptr = t0PtrMaker(t0Collection->size()-1);
      t0SliceAssn->addSingle(sliceT0Ptr, newSlicePtr);
    }

    // Create a Map of PFParticles to allo of their spacepoints
    for (auto const& pfp: slicePFPs){

      recob::PFParticle newPFP(*pfp);
      pfpCollection->push_back(newPFP);
      art::Ptr<recob::PFParticle> newPFPPtr = pfpPtrMaker(pfpCollection->size()-1);
      pfpSliceAssn->addSingle(newPFPPtr, newSlicePtr);

      if (t0Set && !sliceT0Ptr.isNull()){
        t0PFPAssn->addSingle(sliceT0Ptr, newPFPPtr);
      }

      std::vector<art::Ptr<recob::Vertex> > pfpVertexs = fmPFPVertex.at(pfp.key());
      std::vector<art::Ptr<recob::SpacePoint> > pfpSPs = fmPFPSP.at(pfp.key());
      for (auto const& pfpVertex: pfpVertexs){

        geo::Point_t vtxPos(pfpVertex->position());
        //Find the closest SP to the vertex
        double minVtxSPDist = std::numeric_limits<double>::max();
        art::Ptr<recob::SpacePoint> spPtr;
        for (auto const& sp: allSpacePoints){ //TODO: Make slice SPs
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

        if (t0Set && fT0LabelsCorrectT0.at(sliceT0Label)){
          geo::Vector_t posOffset = correctT0(vtxPos, sliceT0->Time(), tpcId);
          vtxPos += posOffset;
        }

        if(fCorrectSCE && fSCE->EnableCalSpatialSCE()){
          geo::Vector_t posOffset = fSCE->GetCalPosOffsets(vtxPos, tpcId.TPC);
          vtxPos += posOffset;
        }

        recob::Vertex newVtx(vtxPos, pfpVertex->covariance(), pfpVertex->chi2(),
            pfpVertex->ndof(), pfpVertex->ID());
        vtxCollection->push_back(newVtx);
        art::Ptr<recob::Vertex> newVtxPtr = vtxPtrMaker(vtxCollection->size()-1);

        pfpVtxAssn->addSingle(newPFPPtr, newVtxPtr);
      }

      for (auto const& sp: pfpSPs){

        //Get the spacepoint position
        geo::Point_t spPos{sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2]};

        art::Ptr<recob::Hit> spHitPtr = fmSPHit.at(sp.key()).front();
        geo::TPCID tpcId = spHitPtr->WireID().asTPCID();
        //TODO: make a switch for drift direction
        // Apply T0 Correction if required
        if (t0Set && fT0LabelsCorrectT0.at(sliceT0Label)){
          geo::Vector_t posOffset = correctT0(spPos, sliceT0->Time(), tpcId);
          spPos += posOffset;
        }

        if(fCorrectSCE && fSCE->EnableCalSpatialSCE()){
          geo::Vector_t posOffset = fSCE->GetCalPosOffsets(spPos, tpcId.TPC);
          spPos += posOffset;
        }

        Double32_t spXYZ[3] = {spPos.X(), spPos.Y(), spPos.Z()};
        recob::SpacePoint correctedSP(spXYZ, sp->ErrXYZ(), sp->Chisq(), sp->ID());

        spCollection->push_back(correctedSP);
        art::Ptr<recob::SpacePoint> spPtr = spPtrMaker(spCollection->size()-1);
        pfpSPAssn->addSingle(newPFPPtr, spPtr);

        spHitAssn->addSingle(spPtr, spHitPtr);
      } // pspSPs
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
    } // slicePFPs
  } // slice

  // std::cout<<"Test: SPs: "<<spCollection->size()<<std::endl;

  evt.put(std::move(t0Collection));
  evt.put(std::move(sliceCollection));
  evt.put(std::move(clusterCollection));
  evt.put(std::move(pfpCollection));
  evt.put(std::move(spCollection));
  evt.put(std::move(vtxCollection));

  evt.put(std::move(t0PFPAssn));
  evt.put(std::move(t0SliceAssn));
  evt.put(std::move(pfpSPAssn));
  evt.put(std::move(spHitAssn));
  evt.put(std::move(pfpVtxAssn));
  evt.put(std::move(pfpSliceAssn));
  evt.put(std::move(pfpClusterAssn));
  evt.put(std::move(clusterHitAssn));
}

geo::Vector_t sce::SCECorrection::correctT0(const geo::Point_t& point, const double& t0,
    const geo::TPCID tpcId) const{

  double t0Offset = fDetProp->DriftVelocity() * t0 / 1e3;

  const geo::TPCGeo& tpcGeo = fGeom->GetElement(tpcId);
  int driftDirection = tpcGeo.DetectDriftDirection();

  switch (std::abs(driftDirection)) {
    case 1: return geo::Vector_t{t0Offset*driftDirection, 0, 0};
    case 2: return geo::Vector_t{0, t0Offset*driftDirection, 0};
    case 3: return geo::Vector_t{0, 0, t0Offset*driftDirection};
    default: return geo::Vector_t{0,0,0}; //TODO: make an exception
  }
}

DEFINE_ART_MODULE(sce::SCECorrection)

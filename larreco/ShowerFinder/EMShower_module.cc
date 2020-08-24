////////////////////////////////////////////////////////////////////////////
// Class:       EMShower
// Module Type: producer
// File:        EMShower_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
//
// Module to make EM showers.
// Takes the output from cluster finding and track finding and combines
// information to make complete 3D shower.
//
// See DUNE-DocDB 1369 (public) for a detailed description.
////////////////////////////////////////////////////////////////////////////

// Framework includes:
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib/pow.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/ArtDataHelper/ToElement.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/EMShowerAlg.h"

// ROOT includes
#include "RtypesCore.h"
#include "TVector3.h"

// C++ STL includes
#include <algorithm>
#include <array>
#include <float.h>
#include <iostream>
#include <map>
#include <math.h>
#include <memory>
#include <stddef.h>
#include <string>
#include <vector>

#include "range/v3/view.hpp"

namespace {
  template <typename T, typename U>
  struct AddMany {
    AddMany(art::Ptr<T> const& ptr, art::Assns<T, U>& assns) : ptr_{ptr}, assns_{assns} {}

    void
    operator()(art::Ptr<U> const& b) const
    {
      assns_.addSingle(ptr_, b);
    }

    art::Ptr<T> const ptr_;
    art::Assns<T, U>& assns_;
  };
}

using lar::to_element;

namespace shower {
  class EMShower;
}

class shower::EMShower : public art::EDProducer {
public:
  EMShower(fhicl::ParameterSet const& pset);

private:
  void produce(art::Event& evt);

  art::InputTag const fHitsModuleLabel;
  art::InputTag const fClusterModuleLabel;
  art::InputTag const fTrackModuleLabel;
  art::InputTag const fPFParticleModuleLabel;
  art::InputTag const fVertexModuleLabel;
  art::InputTag const fCNNEMModuleLabel;
  int const fShower;
  int const fPlane;
  int const fDebug;
  EMShowerAlg const fEMShowerAlg;
  bool const fSaveNonCompleteShowers;
  bool const fFindBadPlanes;
  bool const fMakeSpacePoints;
  bool const fUseCNNtoIDEMPFP;
  bool const fUseCNNtoIDEMHit;
  double const fMinTrackLikeScore;

  art::ServiceHandle<geo::Geometry const> fGeom;
};

shower::EMShower::EMShower(fhicl::ParameterSet const& pset)
  : EDProducer{pset}
  , fHitsModuleLabel{pset.get<art::InputTag>("HitsModuleLabel")}
  , fClusterModuleLabel{pset.get<art::InputTag>("ClusterModuleLabel")}
  , fTrackModuleLabel{pset.get<art::InputTag>("TrackModuleLabel")}
  , fPFParticleModuleLabel{pset.get<art::InputTag>("PFParticleModuleLabel", "")}
  , fVertexModuleLabel{pset.get<art::InputTag>("VertexModuleLabel", "")}
  , fCNNEMModuleLabel{pset.get<art::InputTag>("CNNEMModuleLabel", "")}
  , fShower{pset.get<int>("Shower", -1)}
  , fPlane{pset.get<int>("Plane", -1)}
  , fDebug{pset.get<int>("Debug", 0)}
  , fEMShowerAlg{pset.get<fhicl::ParameterSet>("EMShowerAlg"), fDebug}
  , fSaveNonCompleteShowers{pset.get<bool>("SaveNonCompleteShowers")}
  , fFindBadPlanes{pset.get<bool>("FindBadPlanes")}
  , fMakeSpacePoints{pset.get<bool>("MakeSpacePoints")}
  , fUseCNNtoIDEMPFP{pset.get<bool>("UseCNNtoIDEMPFP")}
  , fUseCNNtoIDEMHit{pset.get<bool>("UseCNNtoIDEMHit")}
  , fMinTrackLikeScore{pset.get<double>("MinTrackLikeScore")}
{
  produces<std::vector<recob::Shower>>();
  produces<std::vector<recob::SpacePoint>>();
  produces<art::Assns<recob::Shower, recob::Hit>>();
  produces<art::Assns<recob::Shower, recob::Cluster>>();
  produces<art::Assns<recob::Shower, recob::Track>>();
  produces<art::Assns<recob::Shower, recob::SpacePoint>>();
  produces<art::Assns<recob::SpacePoint, recob::Hit>>();
}

void
shower::EMShower::produce(art::Event& evt)
{
  // Output -- showers and associations with hits and clusters
  auto showers = std::make_unique<std::vector<recob::Shower>>();
  auto spacePoints = std::make_unique<std::vector<recob::SpacePoint>>();
  auto clusterAssociations = std::make_unique<art::Assns<recob::Shower, recob::Cluster>>();
  auto hitShowerAssociations = std::make_unique<art::Assns<recob::Shower, recob::Hit>>();
  auto trackAssociations = std::make_unique<art::Assns<recob::Shower, recob::Track>>();
  auto spShowerAssociations = std::make_unique<art::Assns<recob::Shower, recob::SpacePoint>>();
  auto hitSpAssociations = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();

  // Event has hits, tracks and clusters found already

  // Hits
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hits;
  if (evt.getByLabel(fHitsModuleLabel, hitHandle)) art::fill_ptr_vector(hits, hitHandle);

  // Tracks
  art::Handle<std::vector<recob::Track>> trackHandle;
  std::vector<art::Ptr<recob::Track>> tracks;
  if (evt.getByLabel(fTrackModuleLabel, trackHandle)) art::fill_ptr_vector(tracks, trackHandle);

  // Clusters
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  std::vector<art::Ptr<recob::Cluster>> clusters;
  if (evt.getByLabel(fClusterModuleLabel, clusterHandle))
    art::fill_ptr_vector(clusters, clusterHandle);

  // PFParticles
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfps;
  if (evt.getByLabel(fPFParticleModuleLabel, pfpHandle)) art::fill_ptr_vector(pfps, pfpHandle);

  // PFParticles
  art::Handle<std::vector<recob::Vertex>> vtxHandle;
  std::vector<art::Ptr<recob::Vertex>> vertices;
  if (evt.getByLabel(fVertexModuleLabel, vtxHandle)) art::fill_ptr_vector(vertices, vtxHandle);

  // Associations
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Track> fmt(hitHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmsp(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster> fmc(hitHandle, evt, fHitsModuleLabel);

  // Make showers
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
  std::vector<std::vector<int>> newShowers;
  std::vector<unsigned int> pfParticles;

  std::map<int, std::vector<int>> clusterToTracks;
  std::map<int, std::vector<int>> trackToClusters;

  if (!pfpHandle.isValid()) {

    // Map between tracks and clusters
    fEMShowerAlg.AssociateClustersAndTracks(clusters, fmh, fmt, clusterToTracks, trackToClusters);

    // Make initial showers
    std::vector<std::vector<int>> initialShowers = fEMShowerAlg.FindShowers(trackToClusters);

    // Deal with views in which 2D reconstruction failed
    std::vector<int> clustersToIgnore;
    if (fFindBadPlanes)
      clustersToIgnore = fEMShowerAlg.CheckShowerPlanes(initialShowers, clusters, fmh);

    if (clustersToIgnore.size() > 0) {
      clusterToTracks.clear();
      trackToClusters.clear();
      fEMShowerAlg.AssociateClustersAndTracks(
        clusters, fmh, fmt, clustersToIgnore, clusterToTracks, trackToClusters);
      newShowers = fEMShowerAlg.FindShowers(trackToClusters);
    }
    else
      newShowers = initialShowers;
  }

  else {

    // Use pfparticle information
    art::FindManyP<recob::Cluster> fmcp(pfpHandle, evt, fPFParticleModuleLabel);
    for (size_t ipfp = 0; ipfp < pfps.size(); ++ipfp) {
      art::Ptr<recob::PFParticle> pfp = pfps[ipfp];
      if (fCNNEMModuleLabel != "" && fUseCNNtoIDEMPFP) { //use CNN to identify EM pfparticle
        auto hitResults = anab::MVAReader<recob::Hit, 4>::create(evt, fCNNEMModuleLabel);
        if (!hitResults) {
          throw cet::exception("EMShower") << "Cannot get MVA results from " << fCNNEMModuleLabel;
        }
        int trkLikeIdx = hitResults->getIndex("track");
        int emLikeIdx = hitResults->getIndex("em");
        if ((trkLikeIdx < 0) || (emLikeIdx < 0)) {
          throw cet::exception("EMShower") << "No em/track labeled columns in MVA data products.";
        }
        if (fmcp.isValid()) { //find clusters
          std::vector<art::Ptr<recob::Hit>> pfphits;
          std::vector<art::Ptr<recob::Cluster>> clus = fmcp.at(ipfp);
          for (size_t iclu = 0; iclu < clus.size(); ++iclu) {
            std::vector<art::Ptr<recob::Hit>> ClusterHits = fmh.at(clus[iclu].key());
            pfphits.insert(pfphits.end(), ClusterHits.begin(), ClusterHits.end());
          }
          if (pfphits.size()) { //find hits
            auto vout = hitResults->getOutput(pfphits);
            double trk_like = -1, trk_or_em = vout[trkLikeIdx] + vout[emLikeIdx];
            if (trk_or_em > 0) {
              trk_like = vout[trkLikeIdx] / trk_or_em;
              if (trk_like < fMinTrackLikeScore) { //EM like
                std::vector<int> clusters;
                for (size_t iclu = 0; iclu < clus.size(); ++iclu) {
                  clusters.push_back(clus[iclu].key());
                }
                if (clusters.size()) {
                  newShowers.push_back(clusters);
                  pfParticles.push_back(ipfp);
                }
              }
            }
          }
        }
        else {
          throw cet::exception("EMShower") << "Cannot get associated cluster for PFParticle "
                                           << fPFParticleModuleLabel.encode() << "[" << ipfp << "]";
        }
      }
      else if (pfp->PdgCode() == 11) { //shower particle
        if (fmcp.isValid()) {
          std::vector<int> clusters;
          std::vector<art::Ptr<recob::Cluster>> clus = fmcp.at(ipfp);
          for (size_t iclu = 0; iclu < clus.size(); ++iclu) {
            clusters.push_back(clus[iclu].key());
          }
          if (clusters.size()) {
            newShowers.push_back(clusters);
            pfParticles.push_back(ipfp);
          }
        }
      }
    }
  }

  // Make output larsoft products
  int showerNum = 0;
  for (auto const& newShower : newShowers) {

    if (showerNum != fShower and fShower != -1) continue;

    // New shower
    if (fDebug > 0) std::cout << "\n\nStart shower " << showerNum << '\n';

    // New associations
    art::PtrVector<recob::Hit> showerHits;
    art::PtrVector<recob::Cluster> showerClusters;
    art::PtrVector<recob::Track> showerTracks;
    art::PtrVector<recob::SpacePoint> showerSpacePoints_p;

    std::vector<int> associatedTracks;

    // Make showers and associations
    for (int const showerCluster : newShower) {

      // Clusters
      art::Ptr<recob::Cluster> const cluster = clusters.at(showerCluster);
      showerClusters.push_back(cluster);

      // Hits
      std::vector<art::Ptr<recob::Hit>> const& showerClusterHits = fmh.at(cluster.key());
      if (fCNNEMModuleLabel != "" && fUseCNNtoIDEMHit) { // use CNN to identify EM hits
        auto hitResults = anab::MVAReader<recob::Hit, 4>::create(evt, fCNNEMModuleLabel);
        if (!hitResults) {
          throw cet::exception("EMShower")
            << "Cannot get MVA results from " << fCNNEMModuleLabel.encode();
        }
        int trkLikeIdx = hitResults->getIndex("track");
        int emLikeIdx = hitResults->getIndex("em");
        if (trkLikeIdx < 0 || emLikeIdx < 0) {
          throw cet::exception("EMShower") << "No em/track labeled columns in MVA data products.";
        }
        for (auto const& showerHit : showerClusterHits) {
          auto vout = hitResults->getOutput(showerHit);
          double trk_like = -1, trk_or_em = vout[trkLikeIdx] + vout[emLikeIdx];
          if (trk_or_em > 0) {
            trk_like = vout[trkLikeIdx] / trk_or_em;
            if (trk_like < fMinTrackLikeScore) { // EM like
              showerHits.push_back(showerHit);
            }
          }
        }
      }
      else {
        for (auto const& showerClusterHit : showerClusterHits)
          showerHits.push_back(showerClusterHit);
      }
      // Tracks
      if (!pfpHandle.isValid()) { // Only do this for non-pfparticle mode
        for (int const clusterTrack : clusterToTracks.at(showerCluster))
          if (not cet::search_all(associatedTracks, clusterTrack))
            associatedTracks.push_back(clusterTrack);
      }
    }

    if (!pfpHandle.isValid()) { // For non-pfparticles, get space points from tracks
      // Tracks and space points
      for (int const trackIndex : associatedTracks) {
        showerTracks.push_back(tracks.at(trackIndex));
      }
    }
    else { // For pfparticles, get space points from hits
      art::FindManyP<recob::SpacePoint> fmspp(showerHits, evt, fPFParticleModuleLabel);
      if (fmspp.isValid()) {
        for (size_t ihit = 0; ihit < showerHits.size(); ++ihit) {
          for (auto const& spPtr : fmspp.at(ihit))
            showerSpacePoints_p.push_back(spPtr);
        }
      }
    }

    if (!pfpHandle.isValid()) {

      // First, order the hits into the correct shower order in each plane
      if (fDebug > 1)
        std::cout << " ------------------ Ordering shower hits --------------------\n";
      std::map<int, std::vector<art::Ptr<recob::Hit>>> showerHitsMap =
        fEMShowerAlg.OrderShowerHits(detProp, showerHits, fPlane);
      if (fDebug > 1)
        std::cout << " ------------------ End ordering shower hits "
                     "--------------------\n";

      // Find the track at the start of the shower
      std::unique_ptr<recob::Track> initialTrack;
      std::map<int, std::vector<art::Ptr<recob::Hit>>> initialTrackHits;
      fEMShowerAlg.FindInitialTrack(detProp, showerHitsMap, initialTrack, initialTrackHits, fPlane);

      // Make space points
      std::vector<std::vector<art::Ptr<recob::Hit>>> hitAssns;
      std::vector<recob::SpacePoint> showerSpacePoints;
      if (fMakeSpacePoints)
        showerSpacePoints = fEMShowerAlg.MakeSpacePoints(detProp, showerHitsMap, hitAssns);
      else {
        for (auto const& trkPtr : showerTracks) {
          for (auto const& trackSpacePoint :
               fmsp.at(trkPtr.key()) | ranges::view::transform(to_element)) {
            showerSpacePoints.push_back(trackSpacePoint);
            hitAssns.push_back(std::vector<art::Ptr<recob::Hit>>());
          }
        }
      }

      // Save space points
      art::PtrMaker<recob::SpacePoint> const make_space_point_ptr{evt};
      size_t firstSpacePoint = spacePoints->size(), nSpacePoint = 0;
      for (auto const& ssp : showerSpacePoints) {
        spacePoints->emplace_back(ssp.XYZ(), ssp.ErrXYZ(), ssp.Chisq(), spacePoints->size());
        auto const index = spacePoints->size() - 1;
        auto const space_point_ptr = make_space_point_ptr(index);
        cet::for_all(hitAssns.at(nSpacePoint), AddMany{space_point_ptr, *hitSpAssociations});
      }
      auto const lastSpacePoint = spacePoints->size();

      // Make shower object and associations
      recob::Shower shower =
        fEMShowerAlg.MakeShower(clockData, detProp, showerHits, initialTrack, initialTrackHits);
      shower.set_id(showerNum);
      if (fSaveNonCompleteShowers or
          (!fSaveNonCompleteShowers and shower.ShowerStart() != TVector3{})) {
        showers->push_back(shower);

        auto const shower_ptr = art::PtrMaker<recob::Shower>{evt}(showers->size() - 1);
        cet::for_all(showerHits, AddMany{shower_ptr, *hitShowerAssociations});
        cet::for_all(showerClusters, AddMany{shower_ptr, *clusterAssociations});
        cet::for_all(showerTracks, AddMany{shower_ptr, *trackAssociations});
        for (size_t i = firstSpacePoint; i < lastSpacePoint; ++i) {
          spShowerAssociations->addSingle(shower_ptr, make_space_point_ptr(i));
        }
      }
      else
        mf::LogInfo("EMShower") << "Discarding shower " << showerNum
                                << " due to incompleteness (SaveNonCompleteShowers == false)";
    }

    else { // pfParticle

      if (vertices.size()) {
        // found the most upstream vertex
        using recob::tracking::Point_t;
        Point_t nuvtx{0, 0, DBL_MAX};
        for (auto const& vtx : vertices) {
          auto const pos = vtx->position();
          if (pos.Z() < nuvtx.Z()) { nuvtx = pos; }
        }

        Point_t shwvtx{0, 0, 0};
        double mindist2 = DBL_MAX;
        for (auto const& sp : showerSpacePoints_p | ranges::view::transform(to_element)) {
          double const dist2 = cet::sum_of_squares(
            nuvtx.X() - sp.XYZ()[0], nuvtx.Y() - sp.XYZ()[1], nuvtx.Z() - sp.XYZ()[2]);
          if (dist2 < mindist2) {
            mindist2 = dist2;
            shwvtx.SetXYZ(sp.XYZ()[0], sp.XYZ()[1], sp.XYZ()[2]);
          }
        }

        art::Ptr<recob::Vertex> bestvtx;
        mindist2 = DBL_MAX;
        for (auto const& vtx : vertices) {
          auto const pos = vtx->position();
          double const dist2 =
            cet::sum_of_squares(pos.X() - shwvtx.X(), pos.Y() - shwvtx.Y(), pos.Z() - shwvtx.Z());
          if (dist2 < mindist2) {
            mindist2 = dist2;
            bestvtx = vtx;
          }
        }

        int iok = 0;
        recob::Shower shower =
          fEMShowerAlg.MakeShower(clockData, detProp, showerHits, bestvtx, iok);
        if (iok == 0) {
          showers->push_back(shower);
          auto const index = showers->size() - 1;
          showers->back().set_id(index);

          auto const shower_ptr = art::PtrMaker<recob::Shower>{evt}(index);
          cet::for_all(showerHits, AddMany{shower_ptr, *hitShowerAssociations});
          cet::for_all(showerClusters, AddMany{shower_ptr, *clusterAssociations});
          cet::for_all(showerTracks, AddMany{shower_ptr, *trackAssociations});
          cet::for_all(showerSpacePoints_p, AddMany{shower_ptr, *spShowerAssociations});
        }
      }
    }
  }

  // Put in event
  evt.put(std::move(showers));
  evt.put(std::move(spacePoints));
  evt.put(std::move(hitShowerAssociations));
  evt.put(std::move(clusterAssociations));
  evt.put(std::move(trackAssociations));
  evt.put(std::move(spShowerAssociations));
  evt.put(std::move(hitSpAssociations));
}

DEFINE_ART_MODULE(shower::EMShower)

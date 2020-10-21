/**
 * @file   TrajCluster_module.cc
 * @brief  Cluster finder using trajectories
 * @author Bruce Baller (baller@fnal.gov)
 *
*
 */

// C/C++ standard libraries
#include <string>

// Framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/TrajClusterAlg.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/PFPUtils.h"

//root includes
#include "TTree.h"

// ... more includes in the implementation section

namespace cluster {
  /**
   * @brief Produces clusters by the TrajCluster algorithm
   *
   * Configuration parameters
   * -------------------------
   *
   * - *HitFinderModuleLabel* (InputTag, mandatory): label of the hits to be
   *   used as input (usually the label of the producing module is enough)
   * - *TrajClusterAlg* (parameter set, mandatory): full configuration for
   *   TrajClusterAlg algorithm
   *
   */
  class TrajCluster : public art::EDProducer {
  public:
    explicit TrajCluster(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;
    void beginJob() override;
    void endJob() override;

    tca::TrajClusterAlg fTCAlg; // define TrajClusterAlg object
    void GetHits(const std::vector<recob::Hit>& inputHits,
                 const geo::TPCID& tpcid,
                 std::vector<std::vector<unsigned int>>& tpcHits);
    void GetHits(const std::vector<recob::Hit>& inputHits,
                 const geo::TPCID& tpcid,
                 const std::vector<recob::Slice>& inputSlices,
                 art::FindManyP<recob::Hit>& hitFromSlc,
                 std::vector<std::vector<unsigned int>>& tpcHits,
                 std::vector<int>& slcIDs);
    art::InputTag fHitModuleLabel;
    art::InputTag fSliceModuleLabel;
    art::InputTag fSpacePointModuleLabel;
    art::InputTag fSpacePointHitAssnLabel;

    unsigned int fMaxSliceHits;
    bool fDoWireAssns;
    bool fDoRawDigitAssns;
    bool fSaveAll2DVertices;
    bool fMakeTracks;
    bool fMakeTrackSpacePoints;
    unsigned int fEventsProcessed;
    unsigned int fNumTracks;
    unsigned int fNumClusters;
  }; // class TrajCluster

} // namespace cluster

//******************************************************************************
//*** implementation
//***

// C/C++ standard libraries
#include <memory> // std::move()

// Framework libraries
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

//LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/ArtDataHelper/HitCreator.h" // recob::HitCollectionAssociator
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"

namespace cluster {

  struct HitLoc {
    unsigned int index; // index of this entry in a sort vector
    unsigned int ctp;   // encoded Cryostat, TPC and Plane
    unsigned int wire;
    int tick;         // hit StartTick using typedef int TDCtick_t in RawTypes.h
    short localIndex; // defined in Hit.h
  };

  //----------------------------------------------------------------------------
  bool
  SortHits(HitLoc const& h1, HitLoc const& h2)
  {
    // sort by hit location (Cryostat, TPC, Plane, Wire, StartTick, hit LocalIndex)
    if (h1.ctp != h2.ctp) return h1.ctp < h2.ctp;
    if (h1.wire != h2.wire) return h1.wire < h2.wire;
    if (h1.tick != h2.tick) return h1.tick < h2.tick;
    return h1.localIndex < h2.localIndex;
  } // SortHits

  //----------------------------------------------------------------------------
  TrajCluster::TrajCluster(fhicl::ParameterSet const& pset)
    : EDProducer{pset}, fTCAlg{pset.get<fhicl::ParameterSet>("TrajClusterAlg")}
  {
    fHitModuleLabel = "NA";
    if (pset.has_key("HitModuleLabel")) fHitModuleLabel = pset.get<art::InputTag>("HitModuleLabel");
    fSliceModuleLabel = "NA";
    if (pset.has_key("SliceModuleLabel"))
      fSliceModuleLabel = pset.get<art::InputTag>("SliceModuleLabel");
    fMaxSliceHits = UINT_MAX;
    if (pset.has_key("MaxSliceHits")) fMaxSliceHits = pset.get<unsigned int>("MaxSliceHits");
    fSpacePointModuleLabel = "NA";
    if (pset.has_key("SpacePointModuleLabel"))
      fSpacePointModuleLabel = pset.get<art::InputTag>("SpacePointModuleLabel");
    fSpacePointHitAssnLabel = "NA";
    if (pset.has_key("SpacePointHitAssnLabel"))
      fSpacePointHitAssnLabel = pset.get<art::InputTag>("SpacePointHitAssnLabel");
    fDoWireAssns = pset.get<bool>("DoWireAssns", true);
    fDoRawDigitAssns = pset.get<bool>("DoRawDigitAssns", true);
    fSaveAll2DVertices = false;
    if(pset.has_key("SaveAll2DVertices")) fSaveAll2DVertices = pset.get<bool>("SaveAll2DVertices");
    fMakeTrackSpacePoints = true;
    if(pset.has_key("MakeTrackSpacePoints")) fMakeTrackSpacePoints = pset.get<bool>("MakeTrackSpacePoints");
    fMakeTracks = true;
    if(pset.has_key("MakeTracks")) fMakeTracks = pset.get<bool>("MakeTracks");

    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionAssociator::declare_products(producesCollector(),"",fDoWireAssns,fDoRawDigitAssns);

    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< std::vector<recob::Seed> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();

    produces< std::vector<recob::PFParticle> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::Vertex> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();
    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Track> >();

    produces< std::vector<recob::Track> >();
    produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

    produces< std::vector<recob::SpacePoint> >();

    produces< art::Assns<recob::Slice, recob::Cluster> >();
    produces< art::Assns<recob::Slice, recob::PFParticle> >();
    produces< art::Assns<recob::Slice, recob::Hit> >();

    // www: declear/create SpacePoint and association between SpacePoint and Hits from TrajCluster (Hit->SpacePoint)
    produces<art::Assns<recob::SpacePoint, recob::Hit>>();
  } // TrajCluster::TrajCluster()

  //----------------------------------------------------------------------------
  void
  TrajCluster::beginJob()
  {
    fEventsProcessed = 0;
    fNumClusters = 0;
    fNumTracks = 0;
  }

  //----------------------------------------------------------------------------
  void
  TrajCluster::endJob()
  {
    // Print out end of job statistics
    if(!tca::tcc.modes[tca::kModeDebug]) return;
    std::vector<unsigned int> const& fAlgModCount = fTCAlg.GetAlgModCount();
    std::vector<std::string> const& fAlgBitNames = fTCAlg.GetAlgBitNames();
    if (fAlgBitNames.size() != fAlgModCount.size()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<"Produced "<<fNumTracks<<" tracks and "<<fNumClusters<<" clusters in "<<fEventsProcessed<<" events\n";
    unsigned short icol = 0;
    for (unsigned short ib = 0; ib < fAlgModCount.size(); ++ib) {
      if (ib == tca::kKilled) continue;
      myprt << std::left << std::setw(18) << fAlgBitNames[ib] << std::right << std::setw(10)
            << fAlgModCount[ib] << " ";
      ++icol;
      if (icol == 4) {
        myprt << "\n";
        icol = 0;
      }
    } // ib
  }   // endJob

  //----------------------------------------------------------------------------
  void
  TrajCluster::produce(art::Event& evt)
  {
    // Get a single hit collection from a HitsModuleLabel or multiple sets of "sliced" hits
    // (aka clusters of hits that are close to each other in 3D) from a SliceModuleLabel.
    // A pointer to the full hit collection is passed to TrajClusterAlg. The hits that are
    // in each slice are reconstructed to find 2D trajectories (that become clusters),
    // 2D vertices (EndPoint2D), 3D vertices and PFParticles. These data products
    // are then collected and written to the event. Each slice is considered as an independent
    // collection of hits with the additional requirement that all hits in a slice reside in
    // one TPC

    // pointers to the slices in the event
    std::vector<art::Ptr<recob::Slice>> slices;
    std::vector<int> slcIDs;
    unsigned int nInputHits = 0;
 
    // get a reference to the Hit collection
    auto inputHits = art::Handle<std::vector<recob::Hit>>();
    if (!evt.getByLabel(fHitModuleLabel, inputHits))
      throw cet::exception("TrajClusterModule")
        << "Failed to get a handle to hit collection '" << fHitModuleLabel.label() << "'\n";
    nInputHits = (*inputHits).size();
    if (!fTCAlg.SetInputHits(*inputHits, evt.run(), evt.event()))
      throw cet::exception("TrajClusterModule")
        << "Failed to process hits from '" << fHitModuleLabel.label() << "'\n";
    // Try to determine the source of the hit collection using the assumption that it was
    // derived from gaushit. If this is successful, pass the handle to TrajClusterAlg to
    // recover hits that were incorrectly removed by disambiguation (DUNE)
    if (fHitModuleLabel != "gaushit") {
      auto sourceHits = art::Handle<std::vector<recob::Hit>>();
      art::InputTag sourceModuleLabel("gaushit");
      if (evt.getByLabel(sourceModuleLabel, sourceHits)) fTCAlg.SetSourceHits(*sourceHits);
    } // look for gaushit collection
 
    // get an optional reference to the Slice collection
    auto inputSlices = art::Handle<std::vector<recob::Slice>>();
    if (fSliceModuleLabel != "NA") {
      fTCAlg.ExpectSlicedHits();
      if (!evt.getByLabel(fSliceModuleLabel, inputSlices))
        throw cet::exception("TrajClusterModule") << "Failed to get a inputSlices";
    } // fSliceModuleLabel specified

    // get an optional reference to the SpacePoint collection
    auto InputSpts = art::Handle<std::vector<recob::SpacePoint>>();
    if (fSpacePointModuleLabel != "NA") {
      if (!evt.getByLabel(fSpacePointModuleLabel, InputSpts))
        throw cet::exception("TrajClusterModule") << "Failed to get a handle to SpacePoints\n";
      tca::evt.sptHits.resize((*InputSpts).size(), {{UINT_MAX, UINT_MAX, UINT_MAX}});
      art::FindManyP<recob::Hit> hitsFromSpt(InputSpts, evt, fSpacePointHitAssnLabel);
      // TrajClusterAlg doesn't use the SpacePoint positions (only the assns to hits) but pass it
      // anyway in case it is useful
      fTCAlg.SetInputSpts(*InputSpts);
      if (!hitsFromSpt.isValid())
        throw cet::exception("TrajClusterModule")
          << "Failed to get a handle to SpacePoint -> Hit assns\n";
      // ensure that the assn is to the inputHit collection
      auto& firstHit = hitsFromSpt.at(0)[0];
      if (firstHit.id() != inputHits.id())
        throw cet::exception("TrajClusterModule")
          << "The SpacePoint -> Hit assn doesn't reference the input hit collection\n";
      tca::evt.sptHits.resize((*InputSpts).size(), {{UINT_MAX, UINT_MAX, UINT_MAX}});
      for (unsigned int isp = 0; isp < (*InputSpts).size(); ++isp) {
        auto& hits = hitsFromSpt.at(isp);
        for (unsigned short iht = 0; iht < hits.size(); ++iht) {
          unsigned short plane = hits[iht]->WireID().Plane;
          tca::evt.sptHits[isp][plane] = hits[iht].key();
        } // iht
      }   // isp
    }     // fSpacePointModuleLabel specified

    if (nInputHits > 0) {

      // look for a debug hit
      if (tca::tcc.dbgStp) {
        tca::debug.Hit = UINT_MAX;
        for (unsigned int iht = 0; iht < (*inputHits).size(); ++iht) {
          auto& hit = (*inputHits)[iht];
          if ((int)hit.WireID().TPC == tca::debug.TPC &&
              (int)hit.WireID().Plane == tca::debug.Plane &&
              (int)hit.WireID().Wire == tca::debug.Wire &&
              hit.PeakTime() > tca::debug.Tick - 10 && hit.PeakTime() < tca::debug.Tick + 10) {
            std::cout << "TCM found debug hit " << iht;
            std::cout << " RMS " << hit.RMS();
            std::cout << " Multiplicity " << hit.Multiplicity();
            std::cout << " GoodnessOfFit " << hit.GoodnessOfFit();
            std::cout << "\n";
            tca::debug.Hit = iht;
            break;
          } // Look for debug hit
        }   // iht
      }     // tca::tcc.dbgStp

      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
      auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
      auto const* geom = lar::providerFrom<geo::Geometry>();
      for (const auto& tpcid : geom->IterateTPCIDs()) {
        // ignore protoDUNE dummy TPCs
        if (geom->TPC(tpcid).DriftDistance() < 25.0) continue;
        // a vector for the subset of hits in each slice in a TPC
        //    slice      hits in this tpc
        std::vector<std::vector<unsigned int>> sltpcHits;
        if (inputSlices.isValid()) {
          // get hits in this TPC and slice
          art::FindManyP<recob::Hit> hitFromSlc(inputSlices, evt, fSliceModuleLabel);
          GetHits(*inputHits, tpcid, *inputSlices, hitFromSlc, sltpcHits, slcIDs);
        }
        else {
          // get hits in this TPC
          // All hits are in one "fake" slice
          GetHits(*inputHits, tpcid, sltpcHits);
          slcIDs.resize(1);
          slcIDs[0] = 1;
        }
        if (sltpcHits.empty()) continue;
        for (unsigned int isl = 0; isl < sltpcHits.size(); ++isl) {
          auto& tpcHits = sltpcHits[isl];
          if (tpcHits.empty()) continue;
          if(tca::tcc.modes[tca::kModeDebug] && tpcHits.size() > 10) 
          std::cout<<"Found "<<tpcHits.size()<<" hits in TPC "<<tpcid.TPC<<"\n";
          // only reconstruct slices with MC-matched hits?
          // sort the slice hits by Cryostat, TPC, Wire, Plane, Start Tick and LocalIndex.
          // This assumes that hits with larger LocalIndex are at larger Tick.
          std::vector<HitLoc> sortVec(tpcHits.size());
          for (unsigned int indx = 0; indx < tpcHits.size(); ++indx) {
            auto& hit = (*inputHits)[tpcHits[indx]];
            sortVec[indx].index = indx;
            sortVec[indx].ctp = tca::EncodeCTP(hit.WireID());
            sortVec[indx].wire = hit.WireID().Wire;
            sortVec[indx].tick = hit.StartTick();
            sortVec[indx].localIndex = hit.LocalIndex();
          } // iht
          std::sort(sortVec.begin(), sortVec.end(), SortHits);
          std::vector tmp = tpcHits;
          for (unsigned int ii = 0; ii < tpcHits.size(); ++ii)
            tpcHits[ii] = tmp[sortVec[ii].index];
          // clear the temp vector
          tmp.resize(0);
          sortVec.resize(0);
          fTCAlg.RunTrajClusterAlg(clockData, detProp, tpcHits, slcIDs[isl]);
        } // isl
      }   // TPC
      // stitch PFParticles between TPCs, create PFP start vertices, etc
      fTCAlg.FinishEvent(detProp);
    } // nInputHits > 0

    // Vectors to hold all data products that will go into the event
    std::vector<recob::Hit> hitCol; // output hit collection
    std::vector<recob::Cluster> clsCol;
    std::vector<recob::PFParticle> pfpCol;
    std::vector<recob::Track> trkCol;
    std::vector<recob::SpacePoint> sptCol;
    std::vector<recob::Vertex> vx3Col;
    std::vector<recob::EndPoint2D> vx2Col;
    std::vector<recob::Seed> sedCol;
    // a vector to correlate inputHits with output hits
    std::vector<unsigned int> newIndex(nInputHits, UINT_MAX);

    // assns for those data products
    // Cluster -> ...
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>> cls_hit_assn(
      new art::Assns<recob::Cluster, recob::Hit>);
    // unsigned short is the end to which a vertex is attached
    std::unique_ptr<art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>> cls_vx2_assn(
      new art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>> cls_vx3_assn(
      new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);
    // PFParticle -> ...
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Cluster>>
      pfp_cls_assn(new art::Assns<recob::PFParticle, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Track>>
      pfp_trk_assn(new art::Assns<recob::PFParticle, recob::Track>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Vertex>>
      pfp_vx3_assn(new art::Assns<recob::PFParticle, recob::Vertex>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Seed>>
      pfp_sed_assn(new art::Assns<recob::PFParticle, recob::Seed>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::SpacePoint>>
      pfp_spt_assn(new art::Assns<recob::PFParticle, recob::SpacePoint>);
    // Track -> ...
    std::unique_ptr<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>
      trk_hit_meta_assn(new art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>);
    // Slice -> ...
    std::unique_ptr<art::Assns<recob::Slice, recob::Cluster>> slc_cls_assn(
      new art::Assns<recob::Slice, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Slice, recob::PFParticle>> slc_pfp_assn(
      new art::Assns<recob::Slice, recob::PFParticle>);
    std::unique_ptr<art::Assns<recob::Slice, recob::Hit>> slc_hit_assn(
      new art::Assns<recob::Slice, recob::Hit>);
    // www: Hit -> SpacePoint
    std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit>>
      spt_hit_assn(new art::Assns<recob::SpacePoint, recob::Hit>);

    // temp struct to get the index of a 2D (or 3D vertex) into vx2Col (or vx3Col)
    // given a slice index and a vertex ID (not UID)
    struct slcVxStruct {
      unsigned short slIndx;
      int ID;
      unsigned short vxColIndx;
    };
    std::vector<slcVxStruct> vx2StrList;
    // vector to map 3V UID -> ID in each sub-slice
    std::vector<slcVxStruct> vx3StrList;

    if (nInputHits > 0) {
      unsigned short nSlices = fTCAlg.GetSlicesSize();
      // define a hit collection begin index to pass to CreateAssn for each cluster
      unsigned int hitColBeginIndex = 0;
      for (unsigned short isl = 0; isl < nSlices; ++isl) {
        unsigned short slcIndex = 0;
        if (!slices.empty()) {
          for (slcIndex = 0; slcIndex < slices.size(); ++slcIndex)
            if (slices[slcIndex]->ID() == slcIDs[isl]) break;
          if (slcIndex == slices.size()) continue;
        }
        auto& slc = fTCAlg.GetSlice(isl);
        // See if there was a serious reconstruction failure that made the sub-slice invalid
        if (!slc.isValid) continue;
        // make EndPoint2Ds
        for (auto& vx2 : slc.vtxs) {
          if (vx2.ID <= 0) continue;
          // skip complete 2D vertices?
          if (!fSaveAll2DVertices && vx2.Vx3ID != 0) continue;
          unsigned int vtxID = vx2.UID;
          unsigned int wire = std::nearbyint(vx2.Pos[0]);
          geo::PlaneID plID = tca::DecodeCTP(vx2.CTP);
          geo::WireID wID = geo::WireID(plID.Cryostat, plID.TPC, plID.Plane, wire);
          geo::View_t view = tca::tcc.geom->View(wID);
          vx2Col.emplace_back((double)vx2.Pos[1] / tca::tcc.unitsPerTick, // Time
                              wID,                                        // WireID
                              vx2.Score,                                  // strength = score
                              vtxID,                                      // ID
                              view,                                       // View
                              0); // total charge - not relevant

          // fill the mapping struct
          slcVxStruct tmp;
          tmp.slIndx = isl;
          tmp.ID = vx2.ID;
          tmp.vxColIndx = vx2Col.size() - 1;
          vx2StrList.push_back(tmp);

        } // vx2
        // make Vertices
        for (auto& vx3 : slc.vtx3s) {
          if (vx3.ID <= 0) continue;
          unsigned int vtxID = vx3.UID;
          double xyz[3];
          xyz[0] = vx3.X;
          xyz[1] = vx3.Y;
          xyz[2] = vx3.Z;
          vx3Col.emplace_back(xyz, vtxID);

          // fill the mapping struct
          slcVxStruct tmp;
          tmp.slIndx = isl;
          tmp.ID = vx3.ID;
          tmp.vxColIndx = vx3Col.size() - 1;
          vx3StrList.push_back(tmp);
        } // vx3
        // Convert the tjs to clusters
        for (auto& tj : slc.tjs) {
          if (tj.AlgMod[tca::kKilled]) continue;
          hitColBeginIndex = hitCol.size();
          for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
            auto& tp = tj.Pts[ipt];
            if (tp.Chg <= 0) continue;
            // index of inputHits indices  of hits used in one TP
            std::vector<unsigned int> tpHits;
            for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
              if (!tp.UseHit[ii]) continue;
              if (tp.Hits[ii] > slc.slHits.size() - 1) { break; } // bad slHits index
              unsigned int allHitsIndex = slc.slHits[tp.Hits[ii]].allHitsIndex;
              if (allHitsIndex > nInputHits - 1) { break; } // bad allHitsIndex
              tpHits.push_back(allHitsIndex);
              if (newIndex[allHitsIndex] != UINT_MAX) {
                std::cout << "Bad Slice " << isl << " tp.Hits " << tp.Hits[ii] << " allHitsIndex "
                          << allHitsIndex;
                std::cout << " old newIndex " << newIndex[allHitsIndex];
                auto& oldhit = (*inputHits)[allHitsIndex];
                std::cout << " old " << oldhit.WireID().Plane << ":" << oldhit.WireID().Wire << ":"
                          << (int)oldhit.PeakTime();
                auto& newhit = hitCol[newIndex[allHitsIndex]];
                std::cout << " new " << newhit.WireID().Plane << ":" << newhit.WireID().Wire << ":"
                          << (int)newhit.PeakTime();
                std::cout << " hitCol size " << hitCol.size();
                std::cout << "\n";
                break;
              }
            } // ii
            // Let the alg define the hit either by merging multiple hits or by a simple copy
            // of a single hit from inputHits
            // Merge hits in the TP that are on the same wire or create hits on multiple wires
            // and update the old hits -> new hits assn (newIndex)
            if (tj.AlgMod[tca::kHaloTj]) {
              // dressed muon - don't merge hits
              for (auto iht : tpHits) {
                hitCol.push_back((*inputHits)[iht]);
                newIndex[iht] = hitCol.size() - 1;
              } // iht
            }
            else {
              fTCAlg.MergeTPHits(tpHits, hitCol, newIndex);
            }
          } // tp
          if (hitCol.empty()) continue;
          // Sum the charge and make the associations
          float sumChg = 0;
          float sumADC = 0;
          for (unsigned int indx = hitColBeginIndex; indx < hitCol.size(); ++indx) {
            auto& hit = hitCol[indx];
            sumChg += hit.Integral();
            sumADC += hit.SummedADC();
            if (!slices.empty() &&
                !util::CreateAssn(*this, evt, hitCol, slices[slcIndex], *slc_hit_assn, indx)) {
              throw art::Exception(art::errors::ProductRegistrationFailure)
                << "Failed to associate hits with Slice";
            }
          } // indx
          geo::View_t view = hitCol[hitColBeginIndex].View();
          auto& firstTP = tj.Pts[tj.EndPt[0]];
          auto& lastTP = tj.Pts[tj.EndPt[1]];
          int clsID = tj.UID;
          if (tj.PDGCode == 11) clsID = -clsID;
          // dressed muon - give the halo cluster the same ID as the parent
          if (tj.AlgMod[tca::kHaloTj]) clsID = -tj.ParentID;
          unsigned int nclhits = hitCol.size() - hitColBeginIndex + 1;
          clsCol.emplace_back(firstTP.Pos[0],                         // Start wire
                              0,                                      // sigma start wire
                              firstTP.Pos[1] / tca::tcc.unitsPerTick, // start tick
                              0,                                      // sigma start tick
                              firstTP.AveChg,                         // start charge
                              firstTP.Ang,                            // start angle
                              0,             // start opening angle (0 for line-like clusters)
                              lastTP.Pos[0], // end wire
                              0,             // sigma end wire
                              lastTP.Pos[1] / tca::tcc.unitsPerTick, // end tick
                              0,                                     // sigma end tick
                              lastTP.AveChg,                         // end charge
                              lastTP.Ang,                            // end angle
                              0,       // end opening angle (0 for line-like clusters)
                              sumChg,  // integral
                              0,       // sigma integral
                              sumADC,  // summed ADC
                              0,       // sigma summed ADC
                              nclhits, // n hits
                              0,       // wires over hits
                              0,       // width (0 for line-like clusters)
                              clsID,   // ID from TrajClusterAlg
                              view,    // view
                              tca::DecodeCTP(tj.CTP), // planeID
                              recob::Cluster::Sentry  // sentry
          );
          if (!util::CreateAssn(
                *this, evt, clsCol, hitCol, *cls_hit_assn, hitColBeginIndex, hitCol.size())) {
            throw art::Exception(art::errors::ProductRegistrationFailure)
              << "Failed to associate hits with cluster ID " << tj.UID;
          } // exception
          // make Slice -> cluster assn
          if (!slices.empty()) {
            if (!util::CreateAssn(*this, evt, clsCol, slices[slcIndex], *slc_cls_assn)) {
              throw art::Exception(art::errors::ProductRegistrationFailure)
                << "Failed to associate slice with PFParticle";
            } // exception
          }   // slices exist
          // Make cluster -> 2V and cluster -> 3V assns
          for (unsigned short end = 0; end < 2; ++end) {
            if (tj.VtxID[end] <= 0) continue;
            for (auto& vx2str : vx2StrList) {
              if (vx2str.slIndx != isl) continue;
              if (vx2str.ID != tj.VtxID[end]) continue;
              if (!util::CreateAssnD(
                    *this, evt, *cls_vx2_assn, clsCol.size() - 1, vx2str.vxColIndx, end)) {
                throw art::Exception(art::errors::ProductRegistrationFailure)
                  << "Failed to associate cluster " << tj.UID << " with EndPoint2D";
              } // exception
              auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
              if (vx2.Vx3ID > 0) {
                std::cout<<" isl "<<isl;
                std::cout<<" UT"<<tj.UID<<" end "<<end<<" 2V"<<vx2.ID<<" -> 3V"<<vx2.Vx3ID<<"\n";
                for (auto vx3str : vx3StrList) {
                  std::cout<<"  vx3str "<<vx3str.slIndx<<" "<<vx3str.ID<<"\n";
                  if (vx3str.slIndx != isl) continue;
                  if (vx3str.ID != vx2.Vx3ID) continue;
                  if (!util::CreateAssnD(
                        *this, evt, *cls_vx3_assn, clsCol.size() - 1, vx3str.vxColIndx, end)) {
                    throw art::Exception(art::errors::ProductRegistrationFailure)
                      << "Failed to associate cluster " << tj.UID << " with Vertex";
                  } // exception
                  std::cout<<"Associated UT"<<tj.UID<<" -> 2V"<<vx2.ID<<" -> 3V"<<vx2.Vx3ID<<"\n";
                  break;
                } // vx3str
              }   // vx2.Vx3ID > 0
              break;
            } // vx2str
          }   // end
        }     // tj (aka cluster)
      }     // slice isl

      // Add PFParticles now that clsCol is filled
      for (unsigned short isl = 0; isl < nSlices; ++isl) {
        unsigned short slcIndex = 0;
        if (!slices.empty()) {
          for (slcIndex = 0; slcIndex < slices.size(); ++slcIndex)
            if (slices[slcIndex]->ID() == slcIDs[isl]) break;
          if (slcIndex == slices.size()) continue;
        }
        auto& slc = fTCAlg.GetSlice(isl);
        // See if there was a serious reconstruction failure that made the slice invalid
        if (!slc.isValid) continue;
        // make PFParticles
        for (size_t ipfp = 0; ipfp < slc.pfps.size(); ++ipfp) {
          auto& pfp = slc.pfps[ipfp];
          if (pfp.ID <= 0) continue;
          // parents and daughters are indexed within a slice so find the index offset in pfpCol
          size_t self = pfpCol.size();
          size_t offset = self - ipfp;
          size_t parentIndex = UINT_MAX;
          if (pfp.ParentUID > 0) parentIndex = pfp.ParentUID + offset - 1;
          std::vector<size_t> dtrIndices(pfp.DtrUIDs.size());
          for (unsigned short idtr = 0; idtr < pfp.DtrUIDs.size(); ++idtr)
            dtrIndices[idtr] = pfp.DtrUIDs[idtr] + offset - 1;
          pfpCol.emplace_back(pfp.PDGCode, self, parentIndex, dtrIndices);
          auto tp3d0 = EndTP3D(pfp, 0);
          auto pos = tp3d0.Pos;
          auto dir = tp3d0.Dir;
          double sp[] = {pos[0], pos[1], pos[2]};
          double sd[] = {dir[0], dir[1], dir[2]};
          double spe[] = {0., 0., 0.};
          double sde[] = {0., 0., 0.};
          sedCol.emplace_back(sp, sd, spe, sde);
          // PFParticle -> clusters
          std::vector<unsigned int> clsIndices;
          for (auto tuid : pfp.TjUIDs) {
            unsigned int clsIndex = 0;
            for (clsIndex = 0; clsIndex < clsCol.size(); ++clsIndex)
              if (abs(clsCol[clsIndex].ID()) == tuid) break;
            if (clsIndex == clsCol.size()) continue;
            clsIndices.push_back(clsIndex);
          } // tjid
          if (!util::CreateAssn(*this,
                                evt,
                                *pfp_cls_assn,
                                pfpCol.size() - 1,
                                clsIndices.begin(),
                                clsIndices.end())) {
            throw art::Exception(art::errors::ProductRegistrationFailure)
              << "Failed to associate clusters with PFParticle";
          } // exception
          // PFParticle -> Vertex
          if (pfp.Vx3ID[0] > 0) {
            for (auto vx3str : vx3StrList) {
              if (vx3str.slIndx != isl) continue;
              if (vx3str.ID != pfp.Vx3ID[0]) continue;
              std::vector<unsigned short> indx(1, vx3str.vxColIndx);
              if (!util::CreateAssn(
                    *this, evt, *pfp_vx3_assn, pfpCol.size() - 1, indx.begin(), indx.end())) {
                throw art::Exception(art::errors::ProductRegistrationFailure)
                  << "Failed to associate PFParticle " << pfp.UID << " with Vertex";
              } // exception
              break;
            } // vx3Index
          }   // start vertex exists
          // PFParticle -> Seed
          if (!sedCol.empty()) {
            if (!util::CreateAssn(*this,
                                  evt,
                                  pfpCol,
                                  sedCol,
                                  *pfp_sed_assn,
                                  sedCol.size() - 1,
                                  sedCol.size(),
                                  pfpCol.size() - 1)) {
              throw art::Exception(art::errors::ProductRegistrationFailure)
                << "Failed to associate seed with PFParticle";
            } // exception
          }   // seeds exist
          // PFParticle -> Slice
          if (!slices.empty()) {
            if (!util::CreateAssn(*this, evt, pfpCol, slices[slcIndex], *slc_pfp_assn)) {
              throw art::Exception(art::errors::ProductRegistrationFailure)
                << "Failed to associate slice with PFParticle";
            } // exception
          }   // slices exist
          // PFParticle -> Track (after making the track)
          if(fMakeTracks) {
            recob::Track trk;
            std::vector<unsigned int> trkHits;
            fTCAlg.MakeTrackFromPFP(pfp, newIndex, trk, trkHits);
            if(!trkHits.empty()) {
              trkCol.push_back(trk);
              // PFParticle -> Track
              if(!util::CreateAssn(*this, evt, pfpCol, trkCol, *pfp_trk_assn, trkCol.size()-1, trkCol.size())) {
                  throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate Track with PFParticle";
              } // fail
              // Track -> Hit + meta data assn
              for(unsigned int iht = 0; iht < trkHits.size(); ++iht) {
                recob::TrackHitMeta metadata(iht,-1);
                if(!util::CreateAssnD(*this, evt, *trk_hit_meta_assn, trkCol.size()-1, trkHits[iht], metadata)) {
                  throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate Hits with Track";
                }
              } // iht
            } // !trkHits.empty()
            if(fMakeTrackSpacePoints) {
              std::vector<recob::SpacePoint> spts;
              // each SpacePoint is associated with one hit
              std::vector<unsigned int> sptsHit;
              fTCAlg.MakeSpacePointsFromPFP(pfp, newIndex, spts, sptsHit);
              if(!spts.empty()) {
                for(unsigned int isp = 0; isp < spts.size(); ++isp) {
                  sptCol.push_back(spts[isp]);
                  // PFParticle -> SpacePoint
                  if(!util::CreateAssn(*this, evt, pfpCol, sptCol, *pfp_spt_assn, sptCol.size()-1, sptCol.size())) {
                    throw art::Exception(art::errors::ProductRegistrationFailure)
                      << "Failed to associate SpacePoint with Track";
                  } // exception
                  // SpacePoint -> Hit
                  if(!util::CreateAssn(*this, evt, sptCol, hitCol, *spt_hit_assn, sptCol.size()-1, sptCol.size(), sptsHit[isp])) {
                    throw art::Exception(art::errors::ProductRegistrationFailure)
                      << "Failed to associate Hit with SpacePoint";
                  } // exception
                } // isp
              } // !spts.empty()
            } // fMakeTrackSpacePoints
          } // fMakeTracks
        }   // ipfp
      }     // isl

      fNumClusters += clsCol.size();
      fNumTracks += trkCol.size();

      // add the hits that weren't used in any slice to hitCol unless this is a
      // special debugging mode and would be a waste of time
      if (!slices.empty() && tca::tcc.recoSlice == 0) {
        auto inputSlices = evt.getValidHandle<std::vector<recob::Slice>>(fSliceModuleLabel);
        art::FindManyP<recob::Hit> hitFromSlc(inputSlices, evt, fSliceModuleLabel);
        for (unsigned int allHitsIndex = 0; allHitsIndex < nInputHits; ++allHitsIndex) {
          if (newIndex[allHitsIndex] != UINT_MAX) continue;
          std::vector<unsigned int> oneHit(1, allHitsIndex);
          fTCAlg.MergeTPHits(oneHit, hitCol, newIndex);
          // find out which slice it is in
          bool gotit = false;
          for (size_t isl = 0; isl < slices.size(); ++isl) {
            auto& hit_in_slc = hitFromSlc.at(isl);
            for (auto& hit : hit_in_slc) {
              if (hit.key() != allHitsIndex) continue;
              gotit = true;
              // Slice -> Hit assn
              if (!util::CreateAssn(*this, evt, hitCol, slices[isl], *slc_hit_assn)) {
                throw art::Exception(art::errors::ProductRegistrationFailure)
                  << "Failed to associate old Hit with Slice";
              } // exception
              break;
            } // hit
            if (gotit) break;
          } // isl
        }   // allHitsIndex
      }     // slices exist
      else {
        // no recob::Slices. Just copy the unused hits
        for (unsigned int allHitsIndex = 0; allHitsIndex < nInputHits; ++allHitsIndex) {
          if (newIndex[allHitsIndex] != UINT_MAX) continue;
          std::vector<unsigned int> oneHit(1, allHitsIndex);
          fTCAlg.MergeTPHits(oneHit, hitCol, newIndex);
        } // allHitsIndex
      }   // recob::Slices
    }     // input hits exist

    // www: find spacepoint from hits (inputHits) through SpacePoint->Hit assns, then create association between spacepoint and trajcluster hits (here, hits in hitCol)
    if (nInputHits > 0) {
      // www: expecting to find spacepoint from hits (inputHits): SpacePoint->Hit assns
      if (fSpacePointModuleLabel != "NA") {
        art::FindManyP<recob::SpacePoint> spFromHit(inputHits, evt, fSpacePointModuleLabel);
        // www: using sp from hit
        for (unsigned int allHitsIndex = 0; allHitsIndex < nInputHits; ++allHitsIndex) {
          if (newIndex[allHitsIndex] == UINT_MAX)
            continue; // skip hits not used in slice (not TrajCluster hits)
          auto& sp_from_hit = spFromHit.at(allHitsIndex);
          for (auto& sp : sp_from_hit) {
            // SpacePoint -> Hit assn
            if(!util::CreateAssn(*this, evt, hitCol, sp, *spt_hit_assn, newIndex[allHitsIndex])) {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate new Hit with SpacePoint";
            } // exception
          }   // sp
        }     // allHitsIndex
      }       // fSpacePointModuleLabel != "NA"
    }         // nInputHits > 0

    // clear the alg data structures
    fTCAlg.ClearResults();

    // convert vectors to unique_ptrs
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>(std::move(hitCol)));
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(clsCol)));
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(vx2Col)));
    std::unique_ptr<std::vector<recob::Vertex> > v3col(new std::vector<recob::Vertex>(std::move(vx3Col)));
    std::unique_ptr<std::vector<recob::PFParticle> > pcol(new std::vector<recob::PFParticle>(std::move(pfpCol)));
    std::unique_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>(std::move(trkCol)));
    std::unique_ptr<std::vector<recob::SpacePoint> > spcol(new std::vector<recob::SpacePoint>(std::move(sptCol)));
    std::unique_ptr<std::vector<recob::Seed> > sdcol(new std::vector<recob::Seed>(std::move(sedCol)));

    // move the cluster collection and the associations into the event:
    if (fHitModuleLabel != "NA") {
      recob::HitRefinerAssociator shcol(evt, fHitModuleLabel, fDoWireAssns, fDoRawDigitAssns);
      shcol.use_hits(std::move(hcol));
      shcol.put_into(evt);
    }
    else {
      recob::HitRefinerAssociator shcol(evt, fSliceModuleLabel, fDoWireAssns, fDoRawDigitAssns);
      shcol.use_hits(std::move(hcol));
      shcol.put_into(evt);
    }
    evt.put(std::move(ccol));
    evt.put(std::move(cls_hit_assn));
    evt.put(std::move(v2col));
    evt.put(std::move(v3col));
    evt.put(std::move(sdcol));
    evt.put(std::move(cls_vx2_assn));
    evt.put(std::move(cls_vx3_assn));
    evt.put(std::move(pcol));
    evt.put(std::move(pfp_cls_assn));
    evt.put(std::move(pfp_vx3_assn));
    evt.put(std::move(pfp_sed_assn));
    evt.put(std::move(pfp_spt_assn));
    evt.put(std::move(pfp_trk_assn));
    evt.put(std::move(tcol));
    evt.put(std::move(spcol));
    evt.put(std::move(trk_hit_meta_assn));
    evt.put(std::move(slc_cls_assn));
    evt.put(std::move(slc_pfp_assn));
    evt.put(std::move(slc_hit_assn));
    evt.put(std::move(spt_hit_assn));
  } // TrajCluster::produce()

  ////////////////////////////////////////////////
  void
  TrajCluster::GetHits(const std::vector<recob::Hit>& inputHits,
                       const geo::TPCID& tpcid,
                       std::vector<std::vector<unsigned int>>& tpcHits)
  {
    // Put hits in this TPC into a single "slice", unless a special debugging mode is specified to
    // only reconstruct hits that are MC-matched
    unsigned int tpc = tpcid.TPC;
    tpcHits.resize(1);
    for (size_t iht = 0; iht < inputHits.size(); ++iht) {
      auto& hit = inputHits[iht];
      if (hit.WireID().TPC == tpc) tpcHits[0].push_back(iht);
    }
  } // GetHits

  ////////////////////////////////////////////////
  void
  TrajCluster::GetHits(const std::vector<recob::Hit>& inputHits,
                       const geo::TPCID& tpcid,
                       const std::vector<recob::Slice>& inputSlices,
                       art::FindManyP<recob::Hit>& hitFromSlc,
                       std::vector<std::vector<unsigned int>>& tpcHits,
                       std::vector<int>& slcIDs)
  {
    // Put the hits in all slices into tpcHits in this TPC
    tpcHits.clear();
    slcIDs.clear();
    if (!hitFromSlc.isValid()) return;

    unsigned int tpc = tpcid.TPC;

    for (size_t isl = 0; isl < inputSlices.size(); ++isl) {
      auto& hit_in_slc = hitFromSlc.at(isl);
      if (hit_in_slc.size() < 3) continue;
      int slcID = inputSlices[isl].ID();
      for (auto& hit : hit_in_slc) {
        if (hit->WireID().TPC != tpc) continue;
        unsigned short indx = 0;
        for (indx = 0; indx < slcIDs.size(); ++indx)
          if (slcID == slcIDs[indx]) break;
        if (indx == slcIDs.size()) {
          slcIDs.push_back(slcID);
          tpcHits.resize(tpcHits.size() + 1);
        }
        tpcHits[indx].push_back(hit.key());
      } // hit
    }   // isl

   } // GetHits

  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(TrajCluster)

} // namespace cluster

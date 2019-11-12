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
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"

//LArSoft includes
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/TrajClusterAlg.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/PFPUtils.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

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
  class TrajCluster: public art::EDProducer {
  public:
    explicit TrajCluster(fhicl::ParameterSet const & pset);

  private:
    void produce(art::Event & evt) override;
    void beginJob() override;
    void endJob() override;


    tca::TrajClusterAlg fTCAlg; // define TrajClusterAlg object
    TTree* showertree;
    void GetHits(const std::vector<recob::Hit>& inputHits,
                 const geo::TPCID& tpcid, std::vector<std::vector<unsigned int>>& tpcHits);
    void GetHits(const std::vector<recob::Hit>& inputHits, 
                 const geo::TPCID& tpcid,
                 const std::vector<recob::Slice>& inputSlices,
                 art::FindManyP<recob::Hit>& hitFromSlc,
                 std::vector<std::vector<unsigned int>>& tpcHits,
                 std::vector<int>& slcIDs);
    void FillMCPList(art::Event & evt, 
                     art::InputTag& fHitTruthModuleLabel, 
                     art::Handle<std::vector<recob::Hit>> & inputHits);


    art::InputTag fHitModuleLabel;
    art::InputTag fSliceModuleLabel;
    art::InputTag fHitTruthModuleLabel;
    art::InputTag fSpacePointModuleLabel;
    art::InputTag fSpacePointHitAssnLabel;

    unsigned int fMaxSliceHits;    
    bool fDoWireAssns;
    bool fDoRawDigitAssns;
    bool fSaveAll2DVertices;
  }; // class TrajCluster

} // namespace cluster

//******************************************************************************
//*** implementation
//***

// C/C++ standard libraries
#include <memory> // std::move()

// Framework libraries
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/ArtDataHelper/HitCreator.h" // recob::HitCollectionAssociator
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Seed.h"


namespace cluster {

  struct HitLoc {
    unsigned int index; // index of this entry in a sort vector
    unsigned int ctp;   // encoded Cryostat, TPC and Plane
    unsigned int wire;
    int tick;           // hit StartTick using typedef int TDCtick_t in RawTypes.h
    short localIndex;   // defined in Hit.h
  };

  //----------------------------------------------------------------------------
  bool SortHits(HitLoc const& h1, HitLoc const& h2)
  {
    // sort by hit location (Cryostat, TPC, Plane, Wire, StartTick, hit LocalIndex)
    if(h1.ctp != h2.ctp) return h1.ctp < h2.ctp;
    if(h1.wire != h2.wire) return h1.wire < h2.wire;
    if(h1.tick != h2.tick) return h1.tick < h2.tick;
    return h1.localIndex < h2.localIndex;
  } // SortHits

  //----------------------------------------------------------------------------
  TrajCluster::TrajCluster(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fTCAlg{pset.get< fhicl::ParameterSet >("TrajClusterAlg")}
  {
    fHitModuleLabel = "NA";
    if(pset.has_key("HitModuleLabel")) fHitModuleLabel = pset.get<art::InputTag>("HitModuleLabel");
    fSliceModuleLabel = "NA";
    if(pset.has_key("SliceModuleLabel")) fSliceModuleLabel = pset.get<art::InputTag>("SliceModuleLabel");
    fHitTruthModuleLabel = "NA";
    if(pset.has_key("HitTruthModuleLabel")) fHitTruthModuleLabel = pset.get<art::InputTag>("HitTruthModuleLabel");
    fMaxSliceHits = UINT_MAX;
    if(pset.has_key("MaxSliceHits")) fMaxSliceHits = pset.get<unsigned int>("MaxSliceHits");
    fSpacePointModuleLabel = "NA";
    if(pset.has_key("SpacePointModuleLabel")) fSpacePointModuleLabel = pset.get<art::InputTag>("SpacePointModuleLabel");
    fSpacePointHitAssnLabel = "NA";
    if(pset.has_key("SpacePointHitAssnLabel")) fSpacePointHitAssnLabel = pset.get<art::InputTag>("SpacePointHitAssnLabel");
    fDoWireAssns = pset.get<bool>("DoWireAssns",true);
    fDoRawDigitAssns = pset.get<bool>("DoRawDigitAssns",true);
    fSaveAll2DVertices = false;
    if(pset.has_key("SaveAll2DVertices")) fSaveAll2DVertices = pset.get<bool>("SaveAll2DVertices");

    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionAssociator::declare_products(producesCollector(),"",fDoWireAssns,fDoRawDigitAssns);

    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< std::vector<recob::Seed> >();
    produces< std::vector<recob::Shower> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();
    produces< art::Assns<recob::Shower, recob::Hit> >();

    produces< std::vector<recob::PFParticle> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::Shower> >();
    produces< art::Assns<recob::PFParticle, recob::Vertex> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();

    produces< art::Assns<recob::Slice, recob::Cluster> >();
    produces< art::Assns<recob::Slice, recob::PFParticle> >();
    produces< art::Assns<recob::Slice, recob::Hit> >();

    produces< std::vector<anab::CosmicTag>>();
    produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();

    // www: declear/create SpacePoint and association between SpacePoint and Hits from TrajCluster (Hit->SpacePoint)
    produces< art::Assns<recob::SpacePoint, recob::Hit> >();
  } // TrajCluster::TrajCluster()

  //----------------------------------------------------------------------------
  void TrajCluster::beginJob()
  {
    art::ServiceHandle<art::TFileService const> tfs;

    showertree = tfs->make<TTree>("showervarstree", "showerVarsTree");
    fTCAlg.DefineShTree(showertree);
//    crtree = tfs->make<TTree>("crtree", "Cosmic removal variables");
//    fTCAlg.DefineCRTree(crtree);
  }

  //----------------------------------------------------------------------------
  void TrajCluster::endJob()
  {
    std::vector<unsigned int> const& fAlgModCount = fTCAlg.GetAlgModCount();
    std::vector<std::string> const& fAlgBitNames = fTCAlg.GetAlgBitNames();
    if(fAlgBitNames.size() != fAlgModCount.size()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<"TrajCluster algorithm counts\n";
    unsigned short icol = 0;
    for(unsigned short ib = 0; ib < fAlgModCount.size(); ++ib) {
      if(ib == tca::kKilled) continue;
      myprt<<std::left<<std::setw(18)<<fAlgBitNames[ib]<<std::right<<std::setw(10)<<fAlgModCount[ib]<<" ";
      ++icol;
      if(icol == 4) { myprt<<"\n"; icol = 0; }
    } // ib
  } // endJob

  //----------------------------------------------------------------------------
  void TrajCluster::produce(art::Event & evt)
  {
    // Get a single hit collection from a HitsModuleLabel or multiple sets of "sliced" hits
    // (aka clusters of hits that are close to each other in 3D) from a SliceModuleLabel.
    // A pointer to the full hit collection is passed to TrajClusterAlg. The hits that are
    // in each slice are reconstructed to find 2D trajectories (that become clusters),
    // 2D vertices (EndPoint2D), 3D vertices, PFParticles and Showers. These data products
    // are then collected and written to the event. Each slice is considered as an independent
    // collection of hits with the additional requirement that all hits in a slice reside in
    // one TPC

    // pointers to the slices in the event
    std::vector<art::Ptr<recob::Slice>> slices;
    std::vector<int> slcIDs;
    unsigned int nInputHits = 0;
    
    // get a reference to the Hit collection
    auto inputHits = art::Handle<std::vector<recob::Hit>>();
    if(!evt.getByLabel(fHitModuleLabel, inputHits)) throw cet::exception("TrajClusterModule")<<"Failed to get a handle to hit collection '"<<fHitModuleLabel.label()<<"'\n";
    nInputHits = (*inputHits).size();
    if(!fTCAlg.SetInputHits(*inputHits, evt.run(), evt.event())) throw cet::exception("TrajClusterModule")<<"Failed to process hits from '"<<fHitModuleLabel.label()<<"'\n";
    
    // get an optional reference to the Slice collection
    auto inputSlices = art::Handle<std::vector<recob::Slice>>();
    if(fSliceModuleLabel != "NA") {
      fTCAlg.ExpectSlicedHits();
      if(!evt.getByLabel(fSliceModuleLabel, inputSlices)) throw cet::exception("TrajClusterModule")<<"Failed to get a inputSlices";
    } // fSliceModuleLabel specified

    // get an optional reference to the SpacePoint collection
    auto InputSpts = art::Handle<std::vector<recob::SpacePoint>>();
    if(fSpacePointModuleLabel != "NA") {
      if(!evt.getByLabel(fSpacePointModuleLabel, InputSpts)) throw cet::exception("TrajClusterModule")<<"Failed to get a handle to SpacePoints\n";
      fTCAlg.SetInputSpts(*InputSpts);
      // Size the Hit -> SpacePoint assn vector
      tca::evt.allHitsSptIndex.resize(nInputHits, UINT_MAX);
      art::FindManyP<recob::Hit> hitsFromSpt (InputSpts, evt, fSpacePointHitAssnLabel);
      if(!hitsFromSpt.isValid()) throw cet::exception("TrajClusterModule")<<"Failed to get a handle to SpacePoint -> Hit assns\n";
      for(unsigned int isp = 0; isp < (*InputSpts).size(); ++isp) {
        auto &hits = hitsFromSpt.at(isp);
        for(auto& hit : hits) tca::evt.allHitsSptIndex[hit.key()] = isp;
      } // isp
    } // fSpacePointModuleLabel specified
    
    // load MCParticles?
    if(!evt.isRealData() && tca::tcc.matchTruth[0] >= 0) FillMCPList(evt, fHitTruthModuleLabel, inputHits);
    
    if(nInputHits > 0) {
      auto const* geom = lar::providerFrom<geo::Geometry>();
      // a list of TPCs that will be considered when comparing with MC
      std::vector<unsigned int> tpcList;
      for(const auto& tpcid : geom->IterateTPCIDs()) {
        // only reconstruct hits in a user-selected TPC in debug mode
        if(tca::tcc.modes[tca::kDebug] && tca::tcc.recoTPC >= 0 && (short)tpcid.TPC != tca::tcc.recoTPC) continue;
        // a vector for the subset of hits in each slice in a TPC
        //    slice      hits in this tpc
        std::vector<std::vector<unsigned int>> sltpcHits;
        if(inputSlices.isValid()) {
          // get hits in this TPC and slice
          art::FindManyP<recob::Hit> hitFromSlc(inputSlices, evt, fSliceModuleLabel);
          GetHits(*inputHits, tpcid, *inputSlices, hitFromSlc, sltpcHits, slcIDs);
        } else {
          // get hits in this TPC
          // All hits are in one "fake" slice
          GetHits(*inputHits, tpcid, sltpcHits);
          slcIDs.resize(1);
          slcIDs[0] = 1;
        }
        if(sltpcHits.empty()) continue;
        for(unsigned short isl = 0; isl < sltpcHits.size(); ++isl) {
          auto& tpcHits = sltpcHits[isl];
          // sort the slice hits by Cryostat, TPC, Wire, Plane, Start Tick and LocalIndex.
          // This assumes that hits with larger LocalIndex are at larger Tick.
          std::vector<HitLoc> sortVec(tpcHits.size());
          for(unsigned int indx = 0; indx < tpcHits.size(); ++indx) {
            auto& hit = (*inputHits)[tpcHits[indx]];
            sortVec[indx].index = indx;
            sortVec[indx].ctp = tca::EncodeCTP(hit.WireID());
            sortVec[indx].wire = hit.WireID().Wire;
            sortVec[indx].tick = hit.StartTick();
            sortVec[indx].localIndex = hit.LocalIndex();
          } // iht
          std::sort(sortVec.begin(), sortVec.end(), SortHits);
          std::vector tmp = tpcHits;
          for(unsigned int ii = 0; ii < tpcHits.size(); ++ii) tpcHits[ii] = tmp[sortVec[ii].index];
          // clear the temp vector
          tmp.resize(0);
          sortVec.resize(0);
          // look for a debug hit
          if(tca::tcc.dbgStp) {
            tca::debug.Hit = UINT_MAX;
            for(unsigned short indx = 0; indx < tpcHits.size(); ++indx) {
              auto& hit = (*inputHits)[tpcHits[indx]];
              if((int)hit.WireID().TPC == tca::debug.TPC &&
                 (int)hit.WireID().Plane == tca::debug.Plane &&
                 (int)hit.WireID().Wire == tca::debug.Wire &&
                 hit.PeakTime() > tca::debug.Tick - 10  && hit.PeakTime() < tca::debug.Tick + 10) {
                std::cout<<"Debug hit "<<tpcHits[indx]<<" found in slice ID "<<slcIDs[isl];
                std::cout<<" RMS "<<hit.RMS();
                std::cout<<" Multiplicity "<<hit.Multiplicity();
                std::cout<<" GoodnessOfFit "<<hit.GoodnessOfFit();
                std::cout<<"\n";
                tca::debug.Hit = tpcHits[indx];
                break;
              } // Look for debug hit
            } // iht
          } // tca::tcc.dbgStp
          fTCAlg.RunTrajClusterAlg(tpcHits, slcIDs[isl]);
          // this is only used for MC truth matching
          tpcList.push_back(tpcid.TPC);
        } // isl
      } // TPC
      // stitch PFParticles between TPCs, create PFP start vertices, etc
      fTCAlg.FinishEvent();
      if(!evt.isRealData()) fTCAlg.fTM.MatchTruth(tpcList);
      if(tca::tcc.matchTruth[0] >= 0) fTCAlg.fTM.PrintResults(evt.event());
      if(tca::tcc.dbgSummary) tca::PrintAll("TCM");
    } // nInputHits > 0

    // Vectors to hold all data products that will go into the event
    std::vector<recob::Hit> hitCol;       // output hit collection
    std::vector<recob::Cluster> clsCol;
    std::vector<recob::PFParticle> pfpCol;
    std::vector<recob::Vertex> vx3Col;
    std::vector<recob::EndPoint2D> vx2Col;
    std::vector<recob::Seed> sedCol;
    std::vector<recob::Shower> shwCol;
    std::vector<anab::CosmicTag> ctCol;
    // a vector to correlate inputHits with output hits
    std::vector<unsigned int> newIndex(nInputHits, UINT_MAX);

    // assns for those data products
    // Cluster -> ...
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>>
      cls_hit_assn(new art::Assns<recob::Cluster, recob::Hit>);
    // unsigned short is the end to which a vertex is attached
    std::unique_ptr<art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>>
      cls_vx2_assn(new art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>>
      cls_vx3_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);
    // Shower -> ...
    std::unique_ptr<art::Assns<recob::Shower, recob::Hit>>
      shwr_hit_assn(new art::Assns<recob::Shower, recob::Hit>);
    // PFParticle -> ...
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Cluster>>
      pfp_cls_assn(new art::Assns<recob::PFParticle, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Shower>>
      pfp_shwr_assn(new art::Assns<recob::PFParticle, recob::Shower>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Vertex>>
      pfp_vx3_assn(new art::Assns<recob::PFParticle, recob::Vertex>);
    std::unique_ptr<art::Assns<recob::PFParticle, anab::CosmicTag>>
      pfp_cos_assn(new art::Assns<recob::PFParticle, anab::CosmicTag>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Seed>>
      pfp_sed_assn(new art::Assns<recob::PFParticle, recob::Seed>);
    // Slice -> ...
    std::unique_ptr<art::Assns<recob::Slice, recob::Cluster>>
      slc_cls_assn(new art::Assns<recob::Slice, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Slice, recob::PFParticle>>
      slc_pfp_assn(new art::Assns<recob::Slice, recob::PFParticle>);
    std::unique_ptr<art::Assns<recob::Slice, recob::Hit>>
      slc_hit_assn(new art::Assns<recob::Slice, recob::Hit>);
    // www: Hit -> SpacePoint
    std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit>>
      sp_hit_assn(new art::Assns<recob::SpacePoint, recob::Hit>);

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

    if(nInputHits > 0) {
      unsigned short nSlices = fTCAlg.GetSlicesSize();
      // define a hit collection begin index to pass to CreateAssn for each cluster
      unsigned int hitColBeginIndex = 0;
      for(unsigned short isl = 0; isl < nSlices; ++isl) {
        unsigned short slcIndex = 0;
        if(!slices.empty()) {
          for(slcIndex = 0; slcIndex < slices.size(); ++slcIndex) if(slices[slcIndex]->ID() == slcIDs[isl]) break;
          if(slcIndex == slices.size()) continue;
        }
        auto& slc = fTCAlg.GetSlice(isl);
        // See if there was a serious reconstruction failure that made the sub-slice invalid
        if(!slc.isValid) continue;
        // make EndPoint2Ds
        for(auto& vx2 : slc.vtxs) {
          if(vx2.ID <= 0) continue;
          // skip complete 2D vertices?
          if(!fSaveAll2DVertices && vx2.Vx3ID != 0) continue;
          unsigned int vtxID = vx2.UID;
          unsigned int wire = std::nearbyint(vx2.Pos[0]);
          geo::PlaneID plID = tca::DecodeCTP(vx2.CTP);
          geo::WireID wID = geo::WireID(plID.Cryostat, plID.TPC, plID.Plane, wire);
          geo::View_t view = tca::tcc.geom->View(wID);
          vx2Col.emplace_back((double)vx2.Pos[1]/tca::tcc.unitsPerTick,  // Time
                              wID,                  // WireID
                              vx2.Score,            // strength = score
                              vtxID,                // ID
                              view,                 // View
                              0);                   // total charge - not relevant

	  // fill the mapping struct
          slcVxStruct tmp;
          tmp.slIndx = isl;
          tmp.ID = vx2.ID;
          tmp.vxColIndx = vx2Col.size() - 1;
          vx2StrList.push_back(tmp);

        } // vx2
        // make Vertices
        for(auto& vx3 : slc.vtx3s) {
          if(vx3.ID <= 0) continue;
          // ignore incomplete vertices
          if(vx3.Wire >= 0) continue;
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
        for(auto& tj : slc.tjs) {
          if(tj.AlgMod[tca::kKilled]) continue;
          hitColBeginIndex = hitCol.size();
          for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
            auto& tp = tj.Pts[ipt];
            if(tp.Chg <= 0) continue;
            // index of inputHits indices  of hits used in one TP
            std::vector<unsigned int> tpHits;
            for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
              if(!tp.UseHit[ii]) continue;
              if(tp.Hits[ii] > slc.slHits.size() - 1) {
                break;
              } // bad slHits index
              unsigned int allHitsIndex = slc.slHits[tp.Hits[ii]].allHitsIndex;
              if(allHitsIndex > nInputHits - 1) {
                break;
              } // bad allHitsIndex
              tpHits.push_back(allHitsIndex);
              if(newIndex[allHitsIndex] != UINT_MAX) {
                std::cout<<"Bad Slice "<<isl<<" tp.Hits "<<tp.Hits[ii]<<" allHitsIndex "<<allHitsIndex;
                std::cout<<" old newIndex "<<newIndex[allHitsIndex];
                auto& oldhit = (*inputHits)[allHitsIndex];
                std::cout<<" old "<<oldhit.WireID().Plane<<":"<<oldhit.WireID().Wire<<":"<<(int)oldhit.PeakTime();
                auto& newhit = hitCol[newIndex[allHitsIndex]];
                std::cout<<" new "<<newhit.WireID().Plane<<":"<<newhit.WireID().Wire<<":"<<(int)newhit.PeakTime();
                std::cout<<" hitCol size "<<hitCol.size();
                std::cout<<"\n";
                break;
              }
            } // ii
            // Let the alg define the hit either by merging multiple hits or by a simple copy
            // of a single hit from inputHits
            // Merge hits in the TP that are on the same wire or create hits on multiple wires
            // and update the old hits -> new hits assn (newIndex)
            if(tj.AlgMod[tca::kHaloTj]) {
              // dressed muon - don't merge hits
              for(auto iht : tpHits) {
                hitCol.push_back((*inputHits)[iht]);
                newIndex[iht] = hitCol.size() - 1;
              } // iht
            } else {
              fTCAlg.MergeTPHits(tpHits, hitCol, newIndex);
            }
          } // tp
          if(hitCol.empty()) continue;
          // Sum the charge and make the associations
          float sumChg = 0;
          float sumADC = 0;
          for(unsigned int indx =  hitColBeginIndex; indx < hitCol.size(); ++indx) {
            auto& hit = hitCol[indx];
            sumChg += hit.Integral();
            sumADC += hit.SummedADC();
            if(!slices.empty() && !util::CreateAssn(*this, evt, hitCol, slices[slcIndex], *slc_hit_assn, indx)) {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate hits with Slice";
            }
          } // indx
          geo::View_t view = hitCol[hitColBeginIndex].View();
          auto& firstTP = tj.Pts[tj.EndPt[0]];
          auto& lastTP = tj.Pts[tj.EndPt[1]];
          int clsID = tj.UID;
          if(tj.AlgMod[tca::kShowerLike]) clsID = -clsID;
          // dressed muon - give the halo cluster the same ID as the parent
          if(tj.AlgMod[tca::kHaloTj]) clsID = -tj.ParentID;
          unsigned int nclhits = hitCol.size() - hitColBeginIndex + 1;
          clsCol.emplace_back(
                              firstTP.Pos[0],         // Start wire
                              0,                      // sigma start wire
                              firstTP.Pos[1]/tca::tcc.unitsPerTick,         // start tick
                              0,                      // sigma start tick
                              firstTP.AveChg,         // start charge
                              firstTP.Ang,            // start angle
                              0,                      // start opening angle (0 for line-like clusters)
                              lastTP.Pos[0],          // end wire
                              0,                      // sigma end wire
                              lastTP.Pos[1]/tca::tcc.unitsPerTick,           // end tick
                              0,                      // sigma end tick
                              lastTP.AveChg,          // end charge
                              lastTP.Ang,             // end angle
                              0,                      // end opening angle (0 for line-like clusters)
                              sumChg,                 // integral
                              0,                      // sigma integral
                              sumADC,                 // summed ADC
                              0,                      // sigma summed ADC
                              nclhits,                // n hits
                              0,                      // wires over hits
                              0,                      // width (0 for line-like clusters)
                              clsID,                  // ID from TrajClusterAlg
                              view,                   // view
                              tca::DecodeCTP(tj.CTP), // planeID
                              recob::Cluster::Sentry  // sentry
                              );
          if(!util::CreateAssn(*this, evt, clsCol, hitCol, *cls_hit_assn, hitColBeginIndex, hitCol.size()))
          {
            throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate hits with cluster ID "<<tj.UID;
          } // exception
          // make Slice -> cluster assn
          if(!slices.empty()) {
            if(!util::CreateAssn(*this, evt, clsCol, slices[slcIndex], *slc_cls_assn))
            {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate slice with PFParticle";
            } // exception
          } // slices exist
          // Make cluster -> 2V and cluster -> 3V assns
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] <= 0) continue;
            for(auto& vx2str : vx2StrList) {
              if(vx2str.slIndx != isl) continue;
              if(vx2str.ID != tj.VtxID[end]) continue;
              if(!util::CreateAssnD(*this, evt, *cls_vx2_assn, clsCol.size() - 1, vx2str.vxColIndx, end)) {
                throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate cluster "<<tj.UID<<" with EndPoint2D";
              } // exception
              auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
              if(vx2.Vx3ID > 0) {
                for(auto vx3str : vx3StrList) {
                  if(vx3str.slIndx != isl) continue;
                  if(vx3str.ID != vx2.Vx3ID) continue;
                  if(!util::CreateAssnD(*this, evt, *cls_vx3_assn, clsCol.size() - 1, vx3str.vxColIndx, end))
                  {
                    throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate cluster "<<tj.UID<<" with Vertex";
                  } // exception
                  break;
                } // vx3str
              } // vx2.Vx3ID > 0
              break;
            } // vx2str
          } // end
        } // tj (aka cluster)
        
        // make Showers
        for(auto& ss3 : slc.showers) {
          if(ss3.ID <= 0) continue;
          recob::Shower shower;
          shower.set_id(ss3.UID);
          shower.set_total_energy(ss3.Energy);
          shower.set_total_energy_err(ss3.EnergyErr);
          shower.set_total_MIPenergy(ss3.MIPEnergy);
          shower.set_total_MIPenergy_err(ss3.MIPEnergyErr);
          shower.set_total_best_plane(ss3.BestPlane);
          TVector3 dir = {ss3.Dir[0], ss3.Dir[1], ss3.Dir[2]};
          shower.set_direction(dir);
          TVector3 dirErr = {ss3.DirErr[0], ss3.DirErr[1], ss3.DirErr[2]};
          shower.set_direction_err(dirErr);
          TVector3 pos = {ss3.Start[0], ss3.Start[1], ss3.Start[2]};
          shower.set_start_point(pos);
          TVector3 posErr = {ss3.StartErr[0], ss3.StartErr[1], ss3.StartErr[2]};
          shower.set_start_point_err(posErr);
          shower.set_dedx(ss3.dEdx);
          shower.set_dedx_err(ss3.dEdxErr);
          shower.set_length(ss3.Len);
          shower.set_open_angle(ss3.OpenAngle);
          shwCol.push_back(shower);
          // make the shower - hit association
          std::vector<unsigned int> shwHits(ss3.Hits.size());
          for(unsigned int iht = 0; iht < ss3.Hits.size(); ++iht) shwHits[iht] = newIndex[ss3.Hits[iht]];
          if(!util::CreateAssn(*this, evt, *shwr_hit_assn, shwCol.size()-1, shwHits.begin(), shwHits.end()))
          {
            throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate hits with Shower";
          } // exception
        } // ss3
      } // slice isl


      // Add PFParticles now that clsCol is filled
      for(unsigned short isl = 0; isl < nSlices; ++isl) {
        unsigned short slcIndex = 0;
        if(!slices.empty()) {
          for(slcIndex = 0; slcIndex < slices.size(); ++slcIndex) if(slices[slcIndex]->ID() == slcIDs[isl]) break;
          if(slcIndex == slices.size()) continue;
        }
        auto& slc = fTCAlg.GetSlice(isl);
        // See if there was a serious reconstruction failure that made the slice invalid
        if(!slc.isValid) continue;
        // make PFParticles
        for(size_t ipfp = 0; ipfp < slc.pfps.size(); ++ipfp) {
          auto& pfp = slc.pfps[ipfp];
          if(pfp.ID <= 0) continue;
          // parents and daughters are indexed within a slice so find the index offset in pfpCol
          size_t self = pfpCol.size();
          size_t offset = self - ipfp;
          size_t parentIndex = UINT_MAX;
          if(pfp.ParentUID > 0) parentIndex = pfp.ParentUID + offset - 1;
          std::vector<size_t> dtrIndices(pfp.DtrUIDs.size());
          for(unsigned short idtr = 0; idtr < pfp.DtrUIDs.size(); ++idtr) dtrIndices[idtr] = pfp.DtrUIDs[idtr] + offset - 1;
          pfpCol.emplace_back(pfp.PDGCode, self, parentIndex, dtrIndices);
          auto pos = PosAtEnd(pfp, 0);
          auto dir = DirAtEnd(pfp, 0);
          double sp[] = {pos[0],pos[1],pos[2]};
          double sd[] = {dir[0],dir[1],dir[2]};
          double spe[] = {0.,0.,0.};
          double sde[] = {0.,0.,0.};
          sedCol.emplace_back(sp,sd,spe,sde);
          // PFParticle -> clusters
          std::vector<unsigned int> clsIndices;
          for(auto tuid : pfp.TjUIDs) {
            unsigned int clsIndex = 0;
            for(clsIndex = 0; clsIndex < clsCol.size(); ++clsIndex) if(abs(clsCol[clsIndex].ID()) == tuid) break;
            if(clsIndex == clsCol.size()) {
              std::cout<<"TrajCluster module invalid P"<<pfp.UID<<" -> T"<<tuid<<" -> cluster index \n";
              continue;
            }
            clsIndices.push_back(clsIndex);
          } // tjid
          if(!util::CreateAssn(*this, evt, *pfp_cls_assn, pfpCol.size()-1, clsIndices.begin(), clsIndices.end()))
          {
            throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate clusters with PFParticle";
          } // exception
          // PFParticle -> Vertex
          if(pfp.Vx3ID[0] > 0) {
            for(auto vx3str : vx3StrList) {
              if(vx3str.slIndx != isl) continue;
              if(vx3str.ID != pfp.Vx3ID[0]) continue;
              std::vector<unsigned short> indx(1, vx3str.vxColIndx);
              if(!util::CreateAssn(*this, evt, *pfp_vx3_assn, pfpCol.size() - 1, indx.begin(), indx.end()))
              {
                throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate PFParticle "<<pfp.UID<<" with Vertex";
              } // exception
              break;
            } // vx3Index
          } // start vertex exists
          // PFParticle -> Seed
          if(!sedCol.empty()) {
            if(!util::CreateAssn(*this, evt, pfpCol, sedCol, *pfp_sed_assn, sedCol.size()-1, sedCol.size(), pfpCol.size()-1))
            {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate seed with PFParticle";
            } // exception
          } // seeds exist
          // PFParticle -> Slice
          if(!slices.empty()) {
            if(!util::CreateAssn(*this, evt, pfpCol, slices[slcIndex], *slc_pfp_assn))
            {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate slice with PFParticle";
            } // exception
          } // slices exist
          // PFParticle -> Shower
          if(pfp.PDGCode == 1111) {
            std::vector<unsigned short> shwIndex(1, 0);
            for(auto& ss3 : slc.showers) {
              if(ss3.ID <= 0) continue;
              if(ss3.PFPIndex == ipfp) break;
              ++shwIndex[0];
            } // ss3
            if(shwIndex[0] < shwCol.size()) {
              if(!util::CreateAssn(*this, evt, *pfp_shwr_assn, pfpCol.size()-1, shwIndex.begin(), shwIndex.end()))
              {
                throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate shower with PFParticle";
              } // exception
            } // valid shwIndex
          } // pfp -> Shower
          // PFParticle cosmic tag
          if(tca::tcc.modes[tca::kTagCosmics]) {
            std::vector<float> tempPt1, tempPt2;
            tempPt1.push_back(-999);
            tempPt1.push_back(-999);
            tempPt1.push_back(-999);
            tempPt2.push_back(-999);
            tempPt2.push_back(-999);
            tempPt2.push_back(-999);
            ctCol.emplace_back(tempPt1, tempPt2, pfp.CosmicScore, anab::CosmicTagID_t::kNotTagged);
            if (!util::CreateAssn(*this, evt, pfpCol, ctCol, *pfp_cos_assn, ctCol.size()-1, ctCol.size())){
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate CosmicTag with PFParticle";
            }
          } // cosmic tag
        } // ipfp
      } // isl

      // add the hits that weren't used in any slice to hitCol unless this is a
      // special debugging mode and would be a waste of time
      if(!slices.empty() && tca::tcc.recoSlice == 0) {
        auto inputSlices = evt.getValidHandle<std::vector<recob::Slice>>(fSliceModuleLabel);
        art::FindManyP<recob::Hit> hitFromSlc(inputSlices, evt, fSliceModuleLabel);
        for(unsigned int allHitsIndex = 0; allHitsIndex < nInputHits; ++allHitsIndex) {
          if(newIndex[allHitsIndex] != UINT_MAX) continue;
          std::vector<unsigned int> oneHit(1, allHitsIndex);
          fTCAlg.MergeTPHits(oneHit, hitCol, newIndex);
          // find out which slice it is in
          bool gotit = false;
          for(size_t isl = 0; isl < slices.size(); ++isl) {
            auto& hit_in_slc = hitFromSlc.at(isl);
            for(auto& hit : hit_in_slc) {
              if(hit.key() != allHitsIndex) continue;
              gotit = true;
              // Slice -> Hit assn
              if(!util::CreateAssn(*this, evt, hitCol, slices[isl], *slc_hit_assn))
              {
                throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate old Hit with Slice";
              } // exception
              break;
            } // hit
            if(gotit) break;
          } // isl
        } // allHitsIndex
      } // slices exist
      else {
        // no recob::Slices. Just copy the unused hits
        for(unsigned int allHitsIndex = 0; allHitsIndex < nInputHits; ++allHitsIndex) {
          if(newIndex[allHitsIndex] != UINT_MAX) continue;
          std::vector<unsigned int> oneHit(1, allHitsIndex);
          fTCAlg.MergeTPHits(oneHit, hitCol, newIndex);
        } // allHitsIndex
      } // recob::Slices
    } // input hits exist

    // www: find spacepoint from hits (inputHits) through SpacePoint->Hit assns, then create association between spacepoint and trajcluster hits (here, hits in hitCol)
    if (nInputHits > 0) {
      // www: expecting to find spacepoint from hits (inputHits): SpacePoint->Hit assns
      if (fSpacePointModuleLabel != "NA") {
        art::FindManyP<recob::SpacePoint> spFromHit (inputHits, evt, fSpacePointModuleLabel);
        // www: using sp from hit
        for (unsigned int allHitsIndex = 0; allHitsIndex < nInputHits; ++allHitsIndex) {
          if (newIndex[allHitsIndex] == UINT_MAX) continue; // skip hits not used in slice (not TrajCluster hits)
          auto & sp_from_hit = spFromHit.at(allHitsIndex);
          for (auto& sp : sp_from_hit) {
            // SpacePoint -> Hit assn
            if(!util::CreateAssn(*this, evt, hitCol, sp, *sp_hit_assn, newIndex[allHitsIndex])) {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate new Hit with SpacePoint";
            } // exception
          } // sp
        } // allHitsIndex
      } // fSpacePointModuleLabel != "NA"
    } // nInputHits > 0

    // clear the alg data structures
    fTCAlg.ClearResults();

    // convert vectors to unique_ptrs
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>(std::move(hitCol)));
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(clsCol)));
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(vx2Col)));
    std::unique_ptr<std::vector<recob::Vertex> > v3col(new std::vector<recob::Vertex>(std::move(vx3Col)));
    std::unique_ptr<std::vector<recob::PFParticle> > pcol(new std::vector<recob::PFParticle>(std::move(pfpCol)));
    std::unique_ptr<std::vector<recob::Seed> > sdcol(new std::vector<recob::Seed>(std::move(sedCol)));
    std::unique_ptr<std::vector<recob::Shower> > scol(new std::vector<recob::Shower>(std::move(shwCol)));
    std::unique_ptr<std::vector<anab::CosmicTag>> ctgcol(new std::vector<anab::CosmicTag>(std::move(ctCol)));

    // move the cluster collection and the associations into the event:
    if(fHitModuleLabel != "NA") {
      recob::HitRefinerAssociator shcol(evt, fHitModuleLabel, fDoWireAssns, fDoRawDigitAssns);
      shcol.use_hits(std::move(hcol));
      shcol.put_into(evt);
    } else {
      recob::HitRefinerAssociator shcol(evt, fSliceModuleLabel, fDoWireAssns, fDoRawDigitAssns);
      shcol.use_hits(std::move(hcol));
      shcol.put_into(evt);
    }
    evt.put(std::move(ccol));
    evt.put(std::move(cls_hit_assn));
    evt.put(std::move(v2col));
    evt.put(std::move(v3col));
    evt.put(std::move(scol));
    evt.put(std::move(sdcol));
    evt.put(std::move(shwr_hit_assn));
    evt.put(std::move(cls_vx2_assn));
    evt.put(std::move(cls_vx3_assn));
    evt.put(std::move(pcol));
    evt.put(std::move(pfp_cls_assn));
    evt.put(std::move(pfp_shwr_assn));
    evt.put(std::move(pfp_vx3_assn));
    evt.put(std::move(pfp_sed_assn));
    evt.put(std::move(slc_cls_assn));
    evt.put(std::move(slc_pfp_assn));
    evt.put(std::move(slc_hit_assn));
    evt.put(std::move(ctgcol));
    evt.put(std::move(pfp_cos_assn));
    evt.put(std::move(sp_hit_assn)); // www: association between sp and hit (trjaclust)
  } // TrajCluster::produce()

  ////////////////////////////////////////////////
  void TrajCluster::GetHits(const std::vector<recob::Hit>& inputHits, 
                            const geo::TPCID& tpcid, 
                            std::vector<std::vector<unsigned int>>& tpcHits)
  {
    // Put hits in this TPC into a single "slice"
    unsigned int tpc = tpcid.TPC;
    tpcHits.resize(1);
    for(unsigned int iht = 0; iht < inputHits.size(); ++iht) {
      auto& hit = inputHits[iht];
      if(hit.WireID().TPC == tpc) tpcHits[0].push_back(iht);
    }
  } // GetHits

  
  ////////////////////////////////////////////////
  void TrajCluster::GetHits(const std::vector<recob::Hit>& inputHits, 
                            const geo::TPCID& tpcid,
                            const std::vector<recob::Slice>& inputSlices,
                            art::FindManyP<recob::Hit>& hitFromSlc,
                            std::vector<std::vector<unsigned int>>& tpcHits,
                            std::vector<int>& slcIDs)
  {
    // Put the hits in all slices into tpcHits in this TPC
    tpcHits.clear();
    slcIDs.clear();
    if(!hitFromSlc.isValid()) return;

    unsigned int tpc = tpcid.TPC;
    
    for(size_t isl = 0; isl < inputSlices.size(); ++isl) {
      auto& hit_in_slc = hitFromSlc.at(isl);
      if(hit_in_slc.size() < 3) continue;
      int slcID = inputSlices[isl].ID();
      for(auto& hit : hit_in_slc) {
        if(hit->WireID().TPC != tpc) continue;
        unsigned short indx = 0;
        for(indx = 0; indx < slcIDs.size(); ++indx) if(slcID == slcIDs[indx]) break;
        if(indx == slcIDs.size()) {
          slcIDs.push_back(slcID);
          tpcHits.resize(tpcHits.size() + 1);
        }
        tpcHits[indx].push_back(hit.key());
      } // hit
    } // isl

  } // GetHits

  ////////////////////////////////////////////////
  void TrajCluster::FillMCPList(art::Event & evt, 
                                art::InputTag& fHitTruthModuleLabel, 
                                art::Handle<std::vector<recob::Hit>> & inputHits)
  {
    // pass a reference to the MCParticle collection to TrajClusterAlg
    auto mcpHandle = art::Handle<std::vector<simb::MCParticle>>();
    if(!evt.getByLabel("largeant", mcpHandle)) throw cet::exception("TrajClusterModule")<<"Failed to get a handle to MCParticles\n";
    fTCAlg.SetMCPHandle(*mcpHandle);
    // size the hit -> MCParticle match vector
    tca::evt.allHitsMCPIndex.resize((*inputHits).size(), UINT_MAX);
    // TODO: Add a check here to ensure that a neutrino vertex exists inside any TPC
    // when checking neutrino reconstruction performance.
    // create a list of MCParticles of interest
    // save MCParticles that have the desired MCTruth origin using
    // the Origin_t typedef enum: kUnknown, kBeamNeutrino, kCosmicRay, kSuperNovaNeutrino, kSingleParticle
    simb::Origin_t origin = (simb::Origin_t)tca::tcc.matchTruth[0];
    // or save them all
    bool anySource = (origin == simb::kUnknown);
    // only reconstruct slices that have hits matched to the desired MC origin?
//    if(tca::tcc.matchTruth.size() > 4 && tca::tcc.matchTruth[4] > 0) requireSliceMCTruthMatch = true;
    // get the assns
    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(inputHits, evt, fHitTruthModuleLabel);
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    sim::ParticleList const& plist = pi_serv->ParticleList();
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      auto& p = (*ipart).second;
      int trackID = p->TrackId();
      const art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
      int KE = 1000 * (p->E() - p->Mass());
      if(!anySource && theTruth->Origin() != origin) continue;
      if(tca::tcc.matchTruth[1] > 1 && KE > 10 && p->Process() == "primary") {
        std::cout<<"TCM: mcp Origin "<<theTruth->Origin()
        <<" pdg "<<p->PdgCode()
        <<std::setw(7)<<KE
        <<" "<<p->Process()
        <<"\n";
      }
    } // ipart
    std::vector<art::Ptr<simb::MCParticle>> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
    for(unsigned int iht = 0; iht < (*inputHits).size(); ++iht) {
      particle_vec.clear(); match_vec.clear();
      try{ particles_per_hit.get(iht, particle_vec, match_vec); }
      catch(...) {
        std::cout<<"BackTrackerHitMatchingData not found\n";
        break;
      }
      if(particle_vec.empty()) continue;
      int trackID = 0;
      for(unsigned short im = 0; im < match_vec.size(); ++im) {
        if(match_vec[im]->ideFraction < 0.5) continue;
        trackID = particle_vec[im]->TrackId();
        break;
      } // im
      if(trackID == 0) continue;
      // see if this is a MCParticle that should be tracked
      const art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
      if(!anySource && theTruth->Origin() != origin) continue;
      // get the index
      for(unsigned int indx = 0; indx < (*mcpHandle).size(); ++indx) {
        auto& mcp = (*mcpHandle)[indx];
        if(mcp.TrackId() != trackID) continue;
        tca::evt.allHitsMCPIndex[iht] = indx;
        break;
      } // indx
    } // iht
  } // FillMCPList
  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(TrajCluster)

} // namespace cluster

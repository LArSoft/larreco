/**
 * @file   TrajCluster_module.cc
 * @brief  Cluster finder using trajectories
 * @author Bruce Baller (baller@fnal.gov)
 * 
*
 */


// C/C++ standard libraries
#include <string>
#include <utility> // std::unique_ptr<>

// Framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Optional/TFileService.h"

//LArSoft includes
#include "larreco/RecoAlg/TrajClusterAlg.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

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
    
    void reconfigure(fhicl::ParameterSet const & pset) ;
    void produce(art::Event & evt) override;
    void beginJob() override;
    void endJob() override;
    
  private:
//    bool SortHits(HitLoc const& h1, HitLoc const& h2);

    std::unique_ptr<tca::TrajClusterAlg> fTCAlg; // define TrajClusterAlg object
    TTree* showertree;
//    TTree* crtree;

    art::InputTag fHitModuleLabel;
    art::InputTag fSlicerModuleLabel;
    
    bool fDoWireAssns;
    bool fDoRawDigitAssns;
    
  }; // class TrajCluster
  
} // namespace cluster

//******************************************************************************
//*** implementation
//***

// C/C++ standard libraries
#include <vector>
#include <memory> // std::move()

// Framework libraries
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"

//LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/ArtDataHelper/HitCreator.h" // recob::HitCollectionAssociator
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"


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
  void TrajCluster::reconfigure(fhicl::ParameterSet const & pset)
  {
    // this trick avoids double configuration on construction
    if (fTCAlg)
      fTCAlg->reconfigure(pset.get< fhicl::ParameterSet >("TrajClusterAlg"));
    else {
      fTCAlg.reset(new tca::TrajClusterAlg(pset.get< fhicl::ParameterSet >("TrajClusterAlg")));
    }

    fHitModuleLabel = pset.get<art::InputTag>("HitModuleLabel");
    fSlicerModuleLabel = "NA";
    // When SlicerModuleLabel is defined, slices of the hit collection will be passed to TrajClusterAlg
    // separately and the results combined at the end
    if(pset.has_key("SlicerModuleLabel")) fSlicerModuleLabel = pset.get<art::InputTag>("SlicerModuleLabel");
    fDoWireAssns = pset.get<bool>("DoWireAssns",true);
    fDoRawDigitAssns = pset.get<bool>("DoRawDigitAssns",true);

  } // TrajCluster::reconfigure()
  
  //----------------------------------------------------------------------------
  TrajCluster::TrajCluster(fhicl::ParameterSet const& pset) {
    
    reconfigure(pset);
    
    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionAssociator::declare_products(*this,"",fDoWireAssns,fDoRawDigitAssns);
    
    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< std::vector<recob::Shower> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();
    produces< art::Assns<recob::Shower, recob::Hit> >();
    
    produces< std::vector<recob::PFParticle> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::Shower> >();
    produces< art::Assns<recob::PFParticle, recob::Vertex> >();

    produces< std::vector<anab::CosmicTag>>();
    produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
  } // TrajCluster::TrajCluster()

  //----------------------------------------------------------------------------
  void TrajCluster::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    showertree = tfs->make<TTree>("showervarstree", "showerVarsTree");
    fTCAlg->DefineShTree(showertree);
//    crtree = tfs->make<TTree>("crtree", "Cosmic removal variables");
//    fTCAlg->DefineCRTree(crtree);
  }
  
  //----------------------------------------------------------------------------
  void TrajCluster::endJob()
  {
    std::vector<unsigned int> const& fAlgModCount = fTCAlg->GetAlgModCount();
    std::vector<std::string> const& fAlgBitNames = fTCAlg->GetAlgBitNames();
    if(fAlgBitNames.size() != fAlgModCount.size()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<"TrajCluster algorithm counts\n";
    unsigned short icol = 0;
    for(unsigned short ib = 0; ib < fAlgModCount.size(); ++ib) {
      if(ib == tca::kKilled) continue;
      myprt<<std::left<<std::setw(16)<<fAlgBitNames[ib]<<std::right<<std::setw(10)<<fAlgModCount[ib]<<" ";
      ++icol;
      if(icol == 4) { myprt<<"\n"; icol = 0; }
    } // ib
  } // endJob
  
  //----------------------------------------------------------------------------
  void TrajCluster::produce(art::Event & evt)
  {
    
    // pass a pointer to the full hit collection to TrajClusterAlg
    auto inputHits = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
    // This is a pointer to a vector of recob::Hits that exist in the event. The hits
    // are not copied.
    fTCAlg->SetInputHits(*inputHits);
    unsigned int nInputHits = (*inputHits).size();
    
    // Next define a vector of indices into inputHits (= evt.fAllHits in TrajClusterAlg) 
    // for each slice for hits associated with 3D-matched PFParticles that were found 
    // with simple 3D clustering (else just the full collection)
    std::vector<std::vector<unsigned int>> slHitsVec;
    if(fSlicerModuleLabel != "NA") {
      // Expecting to find PFParticles (aka slices) -> Clusters -> Hits
//      auto clusHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fSlicerModuleLabel);
//      auto pfpsHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(fSlicerModuleLabel);
    } else {
      // There was no pre-processing of the hits to define logical slices
      // so just consider all hits as one slice
      slHitsVec.resize(1);
      slHitsVec[0].resize(nInputHits);
      for(unsigned int iht = 0; iht < nInputHits; ++iht) slHitsVec[0][iht] = iht;
    } // no input PFParticles

    // First sort the hits in each slice and then reconstruct
    for(auto& slhits : slHitsVec) {
      // sort the slice hits by Cryostat, TPC, Wire, Plane, Start Tick and LocalIndex.
      // This assumes that hits with larger LocalIndex are at larger Tick.
      std::vector<HitLoc> sortVec(slhits.size());
      bool badHit = false;
      for(unsigned short indx = 0; indx < slhits.size(); ++indx) {
        if(slhits[indx] > nInputHits - 1) {
          badHit = true;
          break;
        }
        auto& hit = (*inputHits)[slhits[indx]];
        sortVec[indx].index = indx;
        sortVec[indx].ctp = tca::EncodeCTP(hit.WireID());
        sortVec[indx].wire = hit.WireID().Wire;
        sortVec[indx].tick = hit.StartTick();
        sortVec[indx].localIndex = hit.LocalIndex();
      } // iht
      if(badHit) {
        std::cout<<"TrajCluster found an invalid slice reference to the input hit collection. Ignoring this slice.\n";
        continue;
      }
      std::sort(sortVec.begin(), sortVec.end(), SortHits);
      // make a temporary copy of slhits
      auto tmp = slhits;
      // put slhits into proper order
      for(unsigned short ii = 0; ii < slhits.size(); ++ii) slhits[ii] = tmp[sortVec[ii].index];
      // release the memory for the copy
      tmp.resize(0);
      // reconstruct using the hits in this slice. The data products are stored internally in
      // TrajCluster data formats.
      fTCAlg->RunTrajClusterAlg(slhits);
    } // slhit
    
    // Vectors to hold all data products that will go into the event
    std::vector<recob::Hit> hitCol;       // Hit collection after merging
    std::vector<recob::Cluster> clsCol;
    std::vector<recob::PFParticle> pfpCol;
    std::vector<recob::Vertex> vx3Col;
    std::vector<recob::EndPoint2D> vx2Col;
    std::vector<recob::Shower> shwCol;
    std::vector<anab::CosmicTag> ctCol;
    // a vector for the ID of a cluster in which a hit in hitCol is used. 
    std::vector<int> inClus;
    // a vector of bools to flag hits in the inputHit collection either used in a
    // cluster or declared invalid (used) while merging. Hits in the inputHits collection
    // that are !hitUsed are simply copied from the inputHit collection to hitCol.
    std::vector<bool> hitUsed(nInputHits, false);
    
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
      pfp_vtx_assn(new art::Assns<recob::PFParticle, recob::Vertex>);
    std::unique_ptr<art::Assns<recob::PFParticle, anab::CosmicTag>>
      pfp_cos_assn(new art::Assns<recob::PFParticle, anab::CosmicTag>);

    // 
    int clsID = 0;
//    int vx2ID = 0;
//    int vx3ID = 0;
//    int pfpID = 0;
    unsigned short nSlices = fTCAlg->GetSlicesSize();
    for(unsigned short isl = 0; isl < nSlices; ++isl) {
      auto& slc = fTCAlg->GetSlice(isl);
      // See if there was a serious reconstruction failure that made the slice invalid
      if(!slc.isValid) continue;
      // Convert the tjs to clusters
      bool badSlice = false;
      for(auto& tj : slc.tjs) {
        if(tj.AlgMod[tca::kKilled]) continue;
        ++clsID;
        std::cout<<" got clsID "<<clsID<<"\n";
        for(auto& tp : tj.Pts) {
          if(tp.Chg <= 0) continue;
          std::vector<unsigned int> tpHits;
          for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
            if(!tp.UseHit[ii]) continue;
            unsigned int allHitsIndex = slc.slHits[tp.Hits[ii]].allHitsIndex;
            // consistency check
            if(hitUsed[allHitsIndex]) {
              std::cout<<"TrajCluster module trying to use an already-used hit\n";
              badSlice = true;
              break;
            } // hitUsed
            tpHits.push_back(allHitsIndex);
          } // ii
          if(badSlice) break;
          recob::Hit newHit = fTCAlg->MergeTPHits(tpHits);
          // Let the alg define the hit either by merging or by a simple copy
          if(newHit.Channel() == raw::InvalidChannelID) {
            std::cout<<"TrajCluster module failed merging hits\n";
            badSlice = true;
            break;
          } // MergeTPHits failed
          // add it to the new hits collection
          hitCol.push_back(newHit);
          // make the hit -> cluster association
          inClus.push_back(clsID);
          // flag the inputHits used
          for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
            if(!tp.UseHit[ii]) continue;
            unsigned int allHitsIndex = slc.slHits[tp.Hits[ii]].allHitsIndex;
            hitUsed[allHitsIndex] = true;
          } // ii
        } // tp
        if(badSlice) break;
      } // tj
      if(badSlice) continue;
    } // isl

    // clear the slices vector
    fTCAlg->ClearResults();
/*
    // The YieldHits function merges hits that are used in each trajectory point, creates (and returns)
    // a new hit collection and creates a cluster hit collection with
    // the correct associations. The clusters are retrieved with the GetClusters
//    std::unique_ptr<std::vector<recob::Hit>> newHits = std::make_unique<std::vector<recob::Hit>>(fTCAlg->YieldHits());

//    std::vector<tca::ClusterStore> const& Clusters = fTCAlg->GetClusters(0);
    
    // make EndPoints (aka 2D vertices)
    std::vector<tca::VtxStore> const& EndPts = fTCAlg->GetEndPoints(0);
    art::ServiceHandle<geo::Geometry> geom;
    for(tca::VtxStore const& vtx2: EndPts) {
      if(vtx2.ID == 0) continue;
      unsigned int vtxID = vtx2.ID;
      unsigned int wire = std::nearbyint(vtx2.Pos[0]);
      geo::PlaneID plID = tca::DecodeCTP(vtx2.CTP);
      geo::WireID wID = geo::WireID(plID.Cryostat, plID.TPC, plID.Plane, wire);
      geo::View_t view = geom->View(wID);
      sv2col.emplace_back((double)vtx2.Pos[1],  // Time
                          wID,                  // WireID
                          0,                    // strength - not relevant
                          vtxID,                // ID
                          view,                 // View
                          0);                   // total charge - not relevant
    } // Endpoints

    // make 3D vertices
    std::vector<tca::Vtx3Store> const& Vertices = fTCAlg->GetVertices(0);
    double xyz[3] = {0, 0, 0};
    for(tca::Vtx3Store const& vtx3: Vertices) {
      // ignore incomplete vertices or obsolete
      if(vtx3.Wire >= 0) continue;
      if(vtx3.ID == 0) continue;
      unsigned int vtxID = vtx3.ID;
      xyz[0] = vtx3.X;
      xyz[1] = vtx3.Y;
      xyz[2] = vtx3.Z;
      sv3col.emplace_back(xyz, vtxID);
    } // 3D vertices
    
    // make the clusters and associations
    std::vector<art::Ptr<recob::Hit>> clusterHits;
    float sumChg, sumADC;
    std::vector<size_t> dtrIndices;
    unsigned short clsID = 0;
    for(size_t icl = 0; icl < Clusters.size(); ++icl) {
      tca::ClusterStore const& clstr = Clusters[icl];
      ++clsID;
      geo::PlaneID planeID = tca::DecodeCTP(clstr.CTP);
      unsigned short plane = planeID.Plane;
      unsigned short nclhits = clstr.tclhits.size();

      sumChg = 0;
      sumADC = 0;
      for(unsigned short itt = 0; itt < nclhits; ++itt) {
        unsigned int iht = clstr.tclhits[itt];
        recob::Hit const& hit = (*newHits)[iht];
        sumChg += hit.Integral();
        sumADC += hit.SummedADC();
      } // itt
      geo::View_t view = (*newHits)[clstr.tclhits[0]].View();

      sccol.emplace_back(
          clstr.BeginWir,  // Start wire
          0,                      // sigma start wire
          clstr.BeginTim,         // start tick
          0,                      // sigma start tick
          clstr.BeginChg,         // start charge
          clstr.BeginAng,         // start angle
          0,                      // start opening angle (0 for line-like clusters)
          clstr.EndWir,    // end wire
          0,                      // sigma end wire
          clstr.EndTim,           // end tick
          0,                      // sigma end tick
          clstr.EndChg,           // end charge
          clstr.EndAng,           // end angle
          0,                      // end opening angle (0 for line-like clusters)
          sumChg,                 // integral
          0,                      // sigma integral
          sumADC,                 // summed ADC
          0,                      // sigma summed ADC
          nclhits,                // n hits
          0,                      // wires over hits
          0,                      // width (0 for line-like clusters)
          clstr.ID,               // ID from TrajClusterAlg
          view,                   // view
          planeID,                // plane
          recob::Cluster::Sentry  // sentry
          );
      // make the cluster - hit association
      if(!util::CreateAssn(*this, evt, *cls_hit_assn, sccol.size()-1, clstr.tclhits.begin(), clstr.tclhits.end()))
      {
        throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate hits with cluster ID "<<clstr.ID;
      } // exception

      // make the cluster - endpoint associations
      unsigned short end;
      if(clstr.BeginVtx >= 0) {
        end = 0;
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(tca::Vtx3Store const& vtx3 : Vertices) {
          // ignore killed vertices
          if(vtx3.ID == 0) continue;
          // ignore incomplete vertices
          if(vtx3.Wire > 0) continue;
          if(vtx3.Vx2ID[plane] == 0) continue;
          if(vtx3.Vx2ID[plane] == clstr.BeginVtx) {
            if(!util::CreateAssnD(*this, evt, *cls_vtx_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate cluster "<<icl<<" with vertex";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
      if(clstr.EndVtx >= 0) {
        end = 1;
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(tca::Vtx3Store const& vtx3 : Vertices) {
          // ignore killed vertices
          if(vtx3.ID == 0) continue;
          // ignore incomplete vertices
          if(vtx3.Wire >= 0) continue;
          if(vtx3.Vx2ID[plane] == 0) continue;
          if(vtx3.Vx2ID[plane] == clstr.EndVtx) {
            if(!util::CreateAssnD(*this, evt, *cls_vtx_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate cluster ID "<<clsID<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
    } // icl
    
    // Get the list of PFParticles. These are a subset of the set of 3D matches of trajectory hits
    std::vector<tca::PFPStruct> pfpList = fTCAlg->GetPFParticles(0);
    // Make showers
    std::vector<unsigned int> shwrIndices(pfpList.size(),UINT_MAX);
//    unsigned short nshower = fTCAlg->GetShowerStructSize();
    unsigned short nshower = 0;
    for(unsigned short ish = 0; ish < nshower; ++ish) {
//      tca::ShowerStruct3D const& ss3 = fTCAlg->GetShowerStruct(ish);
      if(ss3.ID == 0) continue;
      recob::Shower shower;
      shower.set_id(ish + 1);
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
      sscol.push_back(shower);
      if(ss3.PFPIndex < shwrIndices.size()) {
        shwrIndices[ss3.PFPIndex] = ish;
      } else {
//        std::cout<<"Invalid PFPIndex "<<ss3.PFPIndex<<" "<<shwrIndices.size()<<"\n";
      }
      // make the shower - hit association
      if(!util::CreateAssn(*this, evt, *shwr_hit_assn, sscol.size()-1, ss3.Hits.begin(), ss3.Hits.end()))
      {
        throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate hits with Shower";
      } // exception
    } // ish
    // get each of the match vector elements and construct the PFParticle
    for(size_t ipfp = 0; ipfp < pfpList.size(); ++ipfp) {
      auto& pfp = pfpList[ipfp];
      if(pfp.ID == 0) continue;
      // ignore special PFParticles (e.g. truth photons)
      if(pfp.PDGCode == 22) continue;
      size_t parentIndex = pfp.ID - 1;
      std::vector<size_t> dtrIndices(pfp.DtrIDs.size());
      for(unsigned short idtr = 0; idtr < pfp.DtrIDs.size(); ++idtr) dtrIndices[idtr] = pfp.DtrIDs[idtr] - 1;
      spcol.emplace_back(pfp.PDGCode, ipfp, parentIndex, dtrIndices);
      // make a list of clusters that are associated with this PFParticle. Trace the association
      // through the trajectories that 
      std::vector<unsigned int> clsIndices;
      for(auto tjid : pfp.TjIDs) {
        if(tjid == 0) {
          std::cout<<"TC: Found an invalid tj ID "<<tjid<<" in P"<<pfp.ID;
          continue;
        }
//        unsigned int clsIndex = fTCAlg->GetTjClusterIndex(tjid);
        unsigned int clsIndex = 0;
        if(clsIndex > Clusters.size() - 1) {
          std::cout<<"Retrieved an invalid cluster index for PFParticle "<<pfp.ID<<" TjID "<<tjid<<". Ignoring it...\n";
          clsIndices.clear();
          break;
        }
        clsIndices.push_back(clsIndex);
      } // tjid
      if(pfp.Vx3ID[0] > (int)Vertices.size()) std::cout<<"TC module: Bad Vtx3DIndex = "<<pfp.Vx3ID[0]<<" size "<<Vertices.size()<<"\n";
      
      // PFParticle - Cluster associations
      if(!util::CreateAssn(*this, evt, *pfp_cls_assn, spcol.size()-1, clsIndices.begin(), clsIndices.end()))
      {
        throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate clusters with PFParticle";
      } // exception
      // PFParticle - Vertex association
      std::vector<unsigned int> vtmp(1);
      // Translate the 3D vertex index into the index of complete 3D vertices that have been put into sv3col
      unsigned short vtxIndex = 0;
      for(unsigned short iv = 0; iv < Vertices.size(); ++iv) {
        if(Vertices[iv].ID == 0) continue;
        if(Vertices[iv].Wire >= 0) continue;
        if(pfp.Vx3ID[0] == Vertices[iv].ID) {
          vtmp[0] = vtxIndex;
          if(!util::CreateAssn(*this, evt, *pfp_vtx_assn, spcol.size()-1, vtmp.begin(), vtmp.end())) 
          {
            throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate vertex with PFParticle";
          }
          break;
        }
        ++vtxIndex;
      } // iv
      // PFParticle - Shower associations
      if (shwrIndices[ipfp]<UINT_MAX) {
        if(!util::CreateAssn(*this, evt, *pfp_shwr_assn, spcol.size()-1, shwrIndices.begin()+ipfp, shwrIndices.begin()+ipfp+1))
        {
          throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate showers with PFParticle";
        } // exception
      }
      // PFParticle - CosmicTag association
      if (fTCAlg->GetTJS().TagCosmics){
        std::vector<float> tempPt1, tempPt2;
        tempPt1.push_back(-999);
        tempPt1.push_back(-999);
        tempPt1.push_back(-999);
        tempPt2.push_back(-999);
        tempPt2.push_back(-999);
        tempPt2.push_back(-999);
        ctcol.emplace_back(tempPt1, tempPt2, pfp.CosmicScore, anab::CosmicTagID_t::kNotTagged);
        if (!util::CreateAssn(*this, evt, spcol, ctcol, *pfp_cos_assn, ctcol.size()-1, ctcol.size())){
          throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate CosmicTag with PFParticle";
        }
      } // TagCosmics
    } // ipfp
*/
    // convert vectors to unique_ptrs
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>(std::move(hitCol)));
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(clsCol)));
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(vx2Col)));
    std::unique_ptr<std::vector<recob::Vertex> > v3col(new std::vector<recob::Vertex>(std::move(vx3Col)));
    std::unique_ptr<std::vector<recob::PFParticle> > pcol(new std::vector<recob::PFParticle>(std::move(pfpCol)));
    std::unique_ptr<std::vector<recob::Shower> > scol(new std::vector<recob::Shower>(std::move(shwCol)));
    std::unique_ptr<std::vector<anab::CosmicTag>> ctgcol(new std::vector<anab::CosmicTag>(std::move(ctCol)));


    // move the cluster collection and the associations into the event:
    recob::HitRefinerAssociator shcol(*this, evt, fHitModuleLabel, fDoWireAssns, fDoRawDigitAssns);
    shcol.use_hits(std::move(hcol));
    shcol.put_into(evt);
    evt.put(std::move(ccol));
    evt.put(std::move(cls_hit_assn));
    evt.put(std::move(v2col));
    evt.put(std::move(v3col));
    evt.put(std::move(scol));
    evt.put(std::move(shwr_hit_assn));
    evt.put(std::move(cls_vx2_assn));
    evt.put(std::move(cls_vx3_assn));
    evt.put(std::move(pcol));
    evt.put(std::move(pfp_cls_assn));
    evt.put(std::move(pfp_shwr_assn));
    evt.put(std::move(pfp_vtx_assn));
    evt.put(std::move(ctgcol));
    evt.put(std::move(pfp_cos_assn));
  } // TrajCluster::produce()
  
  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(TrajCluster)
  
} // namespace cluster


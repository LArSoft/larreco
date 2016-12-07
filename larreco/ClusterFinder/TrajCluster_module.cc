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

//LArSoft includes
#include "larreco/RecoAlg/TrajClusterAlg.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "lardataobj/RecoBase/PFParticle.h"

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
    virtual ~TrajCluster() = default;
    
    void reconfigure(fhicl::ParameterSet const & pset) override;
    void produce(art::Event & evt) override;
    void endJob();
    
  private:
    std::unique_ptr<tca::TrajClusterAlg> fTCAlg; // define TrajClusterAlg object
    
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
#include "lardata/RecoBaseArt/HitCreator.h" // recob::HitCollectionAssociator
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"

namespace cluster {
  
  //----------------------------------------------------------------------------
  void TrajCluster::reconfigure(fhicl::ParameterSet const & pset)
  {
    // this trick avoids double configuration on construction
    if (fTCAlg)
      fTCAlg->reconfigure(pset.get< fhicl::ParameterSet >("TrajClusterAlg"));
    else {
      fTCAlg.reset(new tca::TrajClusterAlg(pset.get< fhicl::ParameterSet >("TrajClusterAlg")));
    }
  } // TrajCluster::reconfigure()
  
  //----------------------------------------------------------------------------
  TrajCluster::TrajCluster(fhicl::ParameterSet const& pset) {
    
    reconfigure(pset);
    
    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionAssociator::declare_products(*this);
    
    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();
    produces< std::vector<recob::PFParticle> >();
    produces< art::Assns<recob::Cluster, recob::PFParticle> >();
  } // TrajCluster::TrajCluster()
  
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

    // look for clusters in all planes
    fTCAlg->RunTrajClusterAlg(evt);
    
    std::unique_ptr<std::vector<recob::Hit>> newHits = std::make_unique<std::vector<recob::Hit>>(std::move(fTCAlg->YieldHits()));

    std::vector<recob::Cluster> sccol;
    std::vector<recob::PFParticle> spcol;
    std::vector<recob::Vertex> sv3col;
    std::vector<recob::EndPoint2D> sv2col;

    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>>
        hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>> 
        cv_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::PFParticle>>
        cp_assn(new art::Assns<recob::Cluster, recob::PFParticle>);

    std::vector<tca::ClusterStore> const& Clusters = fTCAlg->GetClusters();
    
    // make EndPoints (aka 2D vertices)
    std::vector<tca::VtxStore> const& EndPts = fTCAlg->GetEndPoints();
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int vtxID = 0;
    for(tca::VtxStore const& vtx2: EndPts) {
      if(vtx2.NTraj == 0) continue;
      ++vtxID;
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
    // convert 2D Vertex vector to unique_ptrs
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(sv2col)));

    // make 3D vertices
    std::vector<tca::Vtx3Store> const& Vertices = fTCAlg->GetVertices();
    double xyz[3] = {0, 0, 0};
    vtxID = 0;
    for(tca::Vtx3Store const& vtx3: Vertices) {
      // ignore incomplete vertices
      if(vtx3.Ptr2D[0] < 0) continue;
      if(vtx3.Ptr2D[1] < 0) continue;
      if(vtx3.Ptr2D[2] < 0) continue;
      ++vtxID;
      xyz[0] = vtx3.X;
      xyz[1] = vtx3.Y;
      xyz[2] = vtx3.Z;
      sv3col.emplace_back(xyz, vtxID);
    } // 3D vertices
    // convert Vertex vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Vertex> > v3col(new std::vector<recob::Vertex>(std::move(sv3col)));
    
    // make the clusters and associations
    std::vector<art::Ptr<recob::Hit>> clusterHits;
    float sumChg, sumADC;
    std::vector<size_t> dtrIndices;
    unsigned short clsID = 0;
    for(size_t icl = 0; icl < Clusters.size(); ++icl) {
      tca::ClusterStore const& clstr = Clusters[icl];
//      std::cout<<"cls "<<clstr.ID<<" "<<(int)clstr.BeginWir<<":"<<(int)clstr.BeginTim<<" "<<(int)clstr.EndWir<<":"<<(int)clstr.EndTim<<"\n";
      if(clstr.ID < 0) continue;
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
      if(!util::CreateAssn(*this, evt, *hc_assn, sccol.size()-1, clstr.tclhits.begin(), clstr.tclhits.end()))
      {
        throw art::Exception(art::errors::InsertFailure)<<"Failed to associate hits with cluster ID "<<clstr.ID;
      } // exception

      // make the cluster - endpoint associations
      unsigned short end;
      if(clstr.BeginVtx >= 0) {
        end = 0;
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(tca::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Wire < 0) continue;
          if(vtx3.Ptr2D[plane] < 0) continue;
          if(vtx3.Ptr2D[plane] == clstr.BeginVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with vertex";
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
        for(tca::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Ptr2D[0] < 0) continue;
          if(vtx3.Ptr2D[1] < 0) continue;
          if(vtx3.Ptr2D[2] < 0) continue;
          if(vtx3.Ptr2D[plane] == clstr.EndVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster ID "<<clsID<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
      // make PFParticles with ID = cluster ID, parent index = 0, no daughters
      // There is one PFParticle for each cluster
      size_t parent = clstr.ParentCluster;
      if(parent == USHRT_MAX) parent = recob::PFParticle::kPFParticlePrimary;
      dtrIndices.clear();
      for(unsigned short jcl = 0; jcl < Clusters.size(); ++jcl)
        if(Clusters[jcl].ParentCluster == icl) dtrIndices.push_back(jcl);
      spcol.emplace_back((int)clstr.PDGCode, icl, parent, dtrIndices);
      // cluster - PFParticle association
      size_t cEnd = sccol.size();
      size_t cStart = cEnd - 1;
      if(!util::CreateAssn(*this, evt, sccol, spcol, *cp_assn, cStart, cEnd, sccol.size()-1))
      {
        throw art::Exception(art::errors::InsertFailure)
        <<"Failed to associate cluster ID "<<clsID<<" with PFParticle";
      } // exception
    } // icl
    
//    std::cout<<"module clusters "<<Clusters.size()<<" sccol size "<<sccol.size()<<"\n";

    // convert cluster vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));
    std::unique_ptr<std::vector<recob::PFParticle> > pcol(new std::vector<recob::PFParticle>(std::move(spcol)));


    // clean up
    fTCAlg->ClearResults();

    // move the cluster collection and the associations into the event:
    art::InputTag hitModuleLabel = fTCAlg->GetHitFinderModuleLabel();
    recob::HitRefinerAssociator shcol(*this, evt, hitModuleLabel);
    shcol.use_hits(std::move(newHits));
    shcol.put_into(evt);
    evt.put(std::move(ccol));
    evt.put(std::move(hc_assn));
    evt.put(std::move(v2col));
    evt.put(std::move(v3col));
    evt.put(std::move(cv_assn));
    evt.put(std::move(pcol));
    evt.put(std::move(cp_assn));

  } // TrajCluster::produce()
  
  
  
  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(TrajCluster)
  
} // namespace cluster


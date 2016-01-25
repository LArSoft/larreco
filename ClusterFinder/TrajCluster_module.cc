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
#include "art/Utilities/InputTag.h"

//LArSoft includes
#include "RecoAlg/TrajClusterAlg.h"

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
      
    private:
      std::unique_ptr<TrajClusterAlg> fTCAlg; // define TrajClusterAlg object
      
//      art::InputTag fHitFinderModuleLabel; ///< label of module producing input hits
    
  }; // class TrajCluster
  
} // namespace cluster

//******************************************************************************
//*** implementation
//***

// C/C++ standard libraries
#include <vector>
#include <memory> // std::move()

// Framework libraries
#include "art/Utilities/Exception.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Assns.h"

//LArSoft includes
#include "SimpleTypesAndConstants/geo_types.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Vertex.h"
#include "RecoBaseArt/HitCreator.h" // recob::HitCollectionAssociator
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"


namespace cluster {
  
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
    //produces< art::Assns<recob::Cluster, recob::PFParticle, unsigned short> >();
  } // TrajCluster::TrajCluster()
  
  
  //----------------------------------------------------------------------------
  void TrajCluster::reconfigure(fhicl::ParameterSet const & pset)
  {
//    fHitFinderModuleLabel = pset.get<art::InputTag>("HitFinderModuleLabel");
    
    // this trick avoids double configuration on construction
    if (fTCAlg)
      fTCAlg->reconfigure(pset.get< fhicl::ParameterSet >("TrajClusterAlg"));
    else {
      fTCAlg.reset(new TrajClusterAlg
        (pset.get< fhicl::ParameterSet >("TrajClusterAlg")));
    }
  } // TrajCluster::reconfigure()
  
  //----------------------------------------------------------------------------
  void TrajCluster::produce(art::Event & evt)
  {
/*
    //Get the hits for this event:
    art::Handle< std::vector<recob::Hit> > hitVecHandle;
    evt.getByLabel(fHitFinderModuleLabel, hitVecHandle);
    std::vector<art::Ptr<recob::Hit> > allhits;
    allhits.resize(hitVecHandle->size());
    //wrap the hits in art::Ptrs to pass to the Alg
    for (unsigned int iHit = 0; iHit < allhits.size(); iHit++)
      allhits[iHit] = art::Ptr< recob::Hit>(hitVecHandle, iHit);
*/
    // look for clusters in all planes
    fTCAlg->RunTrajClusterAlg(evt);
    
//    std::unique_ptr<std::vector<recob::Hit>> FinalHits
//      (new std::vector<recob::Hit>(std::move(fTCAlg->YieldHits())));
    
    // shcol contains the hit collection
    // and its associations to wires and raw digits;
    // we get the association to raw digits through wire associations
//    recob::HitRefinerAssociator shcol(*this, evt, fHitFinderModuleLabel);
    std::vector<recob::Cluster> sccol;
    std::vector<recob::Vertex> sv3col;
    std::vector<recob::EndPoint2D> sv2col;

    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > 
        hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>> 
        cv_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);

//    std::vector<TrajClusterAlg::ClusterStore> const& Clusters = fTCAlg->GetClusters();
    
//    std::vector<short> const& inClus = fTCAlg->GetinClus();

// Consistency check
/*
    for(unsigned int icl = 0; icl < Clusters.size(); ++icl) {
        TrajClusterAlg::ClusterStore const& clstr = Clusters[icl];
        if(clstr.ID < 0) continue;
        geo::PlaneID planeID = TrajClusterAlg::DecodeCTP(clstr.CTP);
        unsigned short plane = planeID.Plane;
        for(unsigned short ii = 0; ii < clstr.tclhits.size(); ++ii) {
          unsigned int iht = clstr.tclhits[ii];
          recob::Hit const& theHit = FinalHits->at(iht);
          if(theHit.WireID().Plane != plane) {
            mf::LogError("TrajCluster")<<"Cluster-hit plane mis-match "<<theHit.WireID().Plane<<" "<<plane
            <<" in cluster "<<clstr.ID<<" WT "<<clstr.BeginWir<<":"<<(int)clstr.BeginTim<<" cluster CTP "<<clstr.CTP;
            return;
          }
          if(inClus[iht] != clstr.ID) {
            mf::LogError("TrajCluster") << "InClus mis-match " << inClus[iht]
            << " ID " << clstr.ID << " in cluster ID " << clstr.ID<<" cluster ProcCode "<<clstr.ProcCode;;
            return;
          }
        } // ii
      } // icl
*/
    
/*
    // make EndPoints (aka 2D vertices)
    std::vector<TrajClusterAlg::VtxStore> const& EndPts = fTCAlg->GetEndPoints();
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int vtxID = 0, end, wire;
    for(TrajClusterAlg::VtxStore const& vtx2: EndPts) {
      if(vtx2.NClusters == 0) continue;
      ++vtxID;
      wire = (0.5 + vtx2.Wire);
      geo::PlaneID plID = TrajClusterAlg::DecodeCTP(vtx2.CTP);
      geo::WireID wID = geo::WireID(plID.Cryostat, plID.TPC, plID.Plane, wire);
      geo::View_t view = geom->View(wID);
      sv2col.emplace_back((double)vtx2.Time,    // Time
                          wID,                  // WireID
                          0,                    // strength - not relevant
                          vtxID,                // ID
                          view,                 // View
                          0);                   // total charge - not relevant
    } // Endpoints
    // convert 2D Vertex vector to unique_ptrs
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(sv2col)));

    // make 3D vertices
    std::vector<TrajClusterAlg::Vtx3Store> const& Vertices = fTCAlg->GetVertices();
    double xyz[3] = {0, 0, 0};
    vtxID = 0;
    for(TrajClusterAlg::Vtx3Store const& vtx3: Vertices) {
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
    float sumChg, sumADC;
    unsigned int clsID = 0, nclhits;
    for(unsigned int icl = 0; icl < Clusters.size(); ++icl) {
      TrajClusterAlg::ClusterStore const& clstr = Clusters[icl];
      if(clstr.ID < 0) continue;
      ++clsID;
      sumChg = 0;
      sumADC = 0;
      geo::PlaneID planeID = TrajClusterAlg::DecodeCTP(clstr.CTP);
      unsigned short plane = planeID.Plane;
      nclhits = clstr.tclhits.size();
      std::vector<unsigned int> clsHitIndices;

      // correct the hit indices to refer to the valid hits that were just added
      for(unsigned int itt = 0; itt < nclhits; ++itt) {
        unsigned int iht = clstr.tclhits[itt];
        recob::Hit const& hit = FinalHits->at(iht);
        sumChg += hit.Integral();
        sumADC += hit.SummedADC();
      } // itt

      // get the wire, plane from a hit
      unsigned int iht = clstr.tclhits[0];
      
//      geo::View_t view = FinalHits->at(iht).View();
      geo::View_t view = kU;
      sccol.emplace_back(
          (float)clstr.BeginWir,  // Start wire
          0,                      // sigma start wire
          clstr.BeginTim,         // start tick
          0,                      // sigma start tick
          clstr.BeginChg,         // start charge
          clstr.BeginAng,         // start angle
          0,                      // start opening angle (0 for line-like clusters)
          (float)clstr.EndWir,    // end wire
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
          clsID,                  // ID
          view,                   // view
          planeID,                // plane
          recob::Cluster::Sentry  // sentry
          );
      // make the cluster - hit association
      if(!util::CreateAssn(
        *this, evt, *hc_assn, sccol.size()-1, clstr.tclhits.begin(), clstr.tclhits.end())
        )
      {
        throw art::Exception(art::errors::InsertFailure)
          <<"Failed to associate hit "<<iht<<" with cluster "<<icl;
      } // exception
      // make the cluster - endpoint associations
      if(clstr.BeginVtx >= 0) {
        end = 0;
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(TrajClusterAlg::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Ptr2D[0] < 0) continue;
          if(vtx3.Ptr2D[1] < 0) continue;
          if(vtx3.Ptr2D[2] < 0) continue;
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
        for(TrajClusterAlg::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Ptr2D[0] < 0) continue;
          if(vtx3.Ptr2D[1] < 0) continue;
          if(vtx3.Ptr2D[2] < 0) continue;
          if(vtx3.Ptr2D[plane] == clstr.EndVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
    } // icl
*/
    // convert cluster vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));

//    shcol.use_hits(std::move(FinalHits));
    
    // clean up
    fTCAlg->ClearResults();

    // move the hit collection and the associations into the event:
//    shcol.put_into(evt);
    evt.put(std::move(ccol));
    evt.put(std::move(hc_assn));
//    evt.put(std::move(v2col));
//    evt.put(std::move(v3col));
    evt.put(std::move(cv_assn));

  } // TrajCluster::produce()
  
  
  
  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(TrajCluster)
  
} // namespace cluster


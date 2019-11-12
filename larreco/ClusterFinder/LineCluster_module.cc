/**
 * @file   LineCluster_module.cc
 * @brief  Cluster finder using cluster crawler algorithm
 * @author Bruce Baller (bballer@fnal.gov)
 *
 * Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod
 * from cetpkgsupport v1_02_00.
 */

// Framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

//LArSoft includes
#include "larreco/RecoAlg/ClusterCrawlerAlg.h"

// ... more includes in the implementation section

namespace cluster {
  /**
   * @brief Produces clusters by ClusterCrawler algorithm
   *
   * Configuration parameters
   * -------------------------
   *
   * - *HitFinderModuleLabel* (InputTag, mandatory): label of the hits to be
   *   used as input (usually the label of the producing module is enough)
   * - *ClusterCrawlerAlg* (parameter set, mandatory): full configuration for
   *   ClusterCrawlerAlg algorithm
   *
   */
  class LineCluster: public art::EDProducer {

    public:
      explicit LineCluster(fhicl::ParameterSet const & pset);

    private:
      void produce(art::Event & evt) override;

      ClusterCrawlerAlg fCCAlg; // define ClusterCrawlerAlg object

      art::InputTag fHitFinderLabel; ///< label of module producing input hits

      bool fDoWireAssns;
      bool fDoRawDigitAssns;

  }; // class LineCluster

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

//LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/ArtDataHelper/HitCreator.h" // recob::HitCollectionAssociator

namespace cluster {

  //----------------------------------------------------------------------------
  LineCluster::LineCluster(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fCCAlg{pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg")}
  {
    fHitFinderLabel = pset.get<art::InputTag>("HitFinderModuleLabel");
    fDoWireAssns = pset.get<bool>("DoWireAssns",true);
    fDoRawDigitAssns = pset.get<bool>("DoRawDigitAssns",false);

    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionAssociator::declare_products(producesCollector(),"",fDoWireAssns,fDoRawDigitAssns);

    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();
    produces< art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short> >();
  } // LineCluster::LineCluster()


  //----------------------------------------------------------------------------
  void LineCluster::produce(art::Event & evt)
  {
    // fetch the wires needed by CCHitFinder

    // make this accessible to ClusterCrawler_module
    art::ValidHandle< std::vector<recob::Hit>> hitVecHandle
     = evt.getValidHandle<std::vector<recob::Hit>>(fHitFinderLabel);

    // look for clusters in all planes
    fCCAlg.RunCrawler(*hitVecHandle);

    std::unique_ptr<std::vector<recob::Hit>> FinalHits
      (new std::vector<recob::Hit>(std::move(fCCAlg.YieldHits())));

    // shcol contains the hit collection
    // and its associations to wires and raw digits;
    // we get the association to raw digits through wire associations
    recob::HitRefinerAssociator shcol(evt, fHitFinderLabel, fDoWireAssns, fDoRawDigitAssns);
    std::vector<recob::Cluster> sccol;
    std::vector<recob::Vertex> sv3col;
    std::vector<recob::EndPoint2D> sv2col;

    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> >
        hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>>
        cv_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>>
       cep_assn(new art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>);

    std::vector<ClusterCrawlerAlg::ClusterStore> const& Clusters = fCCAlg.GetClusters();


// Consistency check
/*
    std::vector<short> const& inClus = fCCAlg.GetinClus();
    for(unsigned int icl = 0; icl < Clusters.size(); ++icl) {
        ClusterCrawlerAlg::ClusterStore const& clstr = Clusters[icl];
        if(clstr.ID < 0) continue;
        geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
        unsigned short plane = planeID.Plane;
        for(unsigned short ii = 0; ii < clstr.tclhits.size(); ++ii) {
          unsigned int iht = clstr.tclhits[ii];
          recob::Hit const& theHit = FinalHits->at(iht);
          if(theHit.WireID().Plane != plane) {
            mf::LogError("LineCluster")<<"Cluster-hit plane mis-match "<<theHit.WireID().Plane<<" "<<plane
            <<" in cluster "<<clstr.ID<<" WT "<<clstr.BeginWir<<":"<<(int)clstr.BeginTim<<" cluster CTP "<<clstr.CTP;
            return;
          }
          if(inClus[iht] != clstr.ID) {
            mf::LogError("LineCluster") << "InClus mis-match " << inClus[iht]
            << " ID " << clstr.ID << " in cluster ID " << clstr.ID<<" cluster ProcCode "<<clstr.ProcCode;;
            return;
          }
        } // ii
      } // icl
*/
    // make EndPoints (aka 2D vertices)
    std::vector<ClusterCrawlerAlg::VtxStore> const& EndPts = fCCAlg.GetEndPoints();
    std::vector<unsigned int> indxToIndx(EndPts.size());
    art::ServiceHandle<geo::Geometry const> geom;
    unsigned short vtxID = 0, end, wire, ivx;
    for(ivx = 0; ivx < EndPts.size(); ++ivx) {
      if(EndPts[ivx].NClusters == 0) continue;
      indxToIndx[ivx] = vtxID;
       ++vtxID;
//      std::cout<<"EndPt "<<ivx<<" vtxID "<<vtxID<<"\n";
      wire = (0.5 + EndPts[ivx].Wire);
      geo::PlaneID plID = ClusterCrawlerAlg::DecodeCTP(EndPts[ivx].CTP);
      geo::WireID wID = geo::WireID(plID.Cryostat, plID.TPC, plID.Plane, wire);
      geo::View_t view = geom->View(wID);
      sv2col.emplace_back((double)EndPts[ivx].Time,    // Time
                          wID,                  // WireID
                          0,                    // strength - not relevant
                          vtxID,                // ID
                          view,                 // View
                          0);                   // total charge - not relevant
    } // iv
    // convert 2D Vertex vector to unique_ptrs
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(sv2col)));

    // make 3D vertices
    std::vector<ClusterCrawlerAlg::Vtx3Store> const& Vertices = fCCAlg.GetVertices();
    double xyz[3] = {0, 0, 0};
    vtxID = 0;
    for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertices) {
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
      ClusterCrawlerAlg::ClusterStore const& clstr = Clusters[icl];
      if(clstr.ID < 0) continue;
      ++clsID;
      sumChg = 0;
      sumADC = 0;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
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

      geo::View_t view = FinalHits->at(iht).View();
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
        throw art::Exception(art::errors::ProductRegistrationFailure)
          <<"Failed to associate hit "<<iht<<" with cluster "<<icl;
      } // exception
      // make the cluster - EndPoint2D and Vertex associations
      if(clstr.BeginVtx >= 0) {
        end = 0;
//        std::cout<<clstr.ID<<" clsID "<<clsID<<" Begin vtx "<<clstr.BeginVtx<<" vtxID "<<indxToIndx[clstr.BeginVtx]<<"\n";
        if(!util::CreateAssnD(*this, evt, *cep_assn, clsID - 1, indxToIndx[clstr.BeginVtx], end))
        {
          throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate cluster "<<clsID<<" with EndPoint2D "<<clstr.BeginVtx;
        } // exception
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Ptr2D[0] < 0) continue;
          if(vtx3.Ptr2D[1] < 0) continue;
          if(vtx3.Ptr2D[2] < 0) continue;
          if(vtx3.Ptr2D[plane] == clstr.BeginVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::ProductRegistrationFailure)
                <<"Failed to associate cluster "<<icl<<" with vertex";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
      if(clstr.EndVtx >= 0) {
        end = 1;
//        std::cout<<clstr.ID<<" clsID "<<clsID<<" End   vtx "<<clstr.EndVtx<<" vtxID "<<indxToIndx[clstr.EndVtx]<<"\n";
        if(!util::CreateAssnD(*this, evt, *cep_assn, clsID - 1, indxToIndx[clstr.EndVtx], end))
        {
          throw art::Exception(art::errors::ProductRegistrationFailure)<<"Failed to associate cluster "<<clsID<<" with EndPoint2D "<<clstr.BeginVtx;
        } // exception
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Ptr2D[0] < 0) continue;
          if(vtx3.Ptr2D[1] < 0) continue;
          if(vtx3.Ptr2D[2] < 0) continue;
          if(vtx3.Ptr2D[plane] == clstr.EndVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::ProductRegistrationFailure)
                <<"Failed to associate cluster "<<icl<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
    } // icl

    // convert cluster vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));

    shcol.use_hits(std::move(FinalHits));

    // clean up
    fCCAlg.ClearResults();

    // move the hit collection and the associations into the event:
    shcol.put_into(evt);
    evt.put(std::move(ccol));
    evt.put(std::move(hc_assn));
    evt.put(std::move(v2col));
    evt.put(std::move(v3col));
    evt.put(std::move(cv_assn));
    evt.put(std::move(cep_assn));

  } // LineCluster::produce()



  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(LineCluster)

} // namespace cluster

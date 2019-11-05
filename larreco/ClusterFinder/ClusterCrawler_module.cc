////////////////////////////////////////////////////////////////////////
// Class:       ClusterCrawler
// Module Type: producer
// File:        ClusterCrawler_module.cc
//
// Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod
// from cetpkgsupport v1_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>

//LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h" // recob::HitCollectionAssociator
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/CCHitFinderAlg.h"
#include "larreco/RecoAlg/ClusterCrawlerAlg.h"

namespace cluster {
  class ClusterCrawler;
}

class cluster::ClusterCrawler : public art::EDProducer {

public:
  explicit ClusterCrawler(fhicl::ParameterSet const & pset);

private:
  void produce(art::Event & evt) override;

  hit::CCHitFinderAlg fCCHFAlg; // define CCHitFinderAlg object
  ClusterCrawlerAlg fCCAlg; // define ClusterCrawlerAlg object
  std::string fCalDataModuleLabel; ///< label of module producing input wires
};


namespace cluster {

  ClusterCrawler::ClusterCrawler(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    fCCHFAlg           (pset.get<fhicl::ParameterSet>("CCHitFinderAlg")),
    fCCAlg             (pset.get<fhicl::ParameterSet>("ClusterCrawlerAlg")),
    fCalDataModuleLabel(pset.get<std::string>("CalDataModuleLabel"))
  {
    mf::LogWarning("ClusterCrawler") <<
      "\nClusterCrawler module has been deprecated and will be removed."
      "\nIt is now replaced by HitFinder and LineCluster modules."
      ;

    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionAssociator::declare_products(producesCollector());

    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();
  }

  void ClusterCrawler::produce(art::Event & evt)
  {
    // fetch the wires needed by CCHitFinder

    // make this accessible to ClusterCrawler_module
    art::ValidHandle< std::vector<recob::Wire>> wireVecHandle
      = evt.getValidHandle<std::vector<recob::Wire>>(fCalDataModuleLabel);

    // find hits in all planes
    fCCHFAlg.RunCCHitFinder(*wireVecHandle);

    // extract the result of the algorithm (it's moved)
    std::vector<recob::Hit> FirstHits = fCCHFAlg.YieldHits();

    // look for clusters in all planes
    fCCAlg.RunCrawler(FirstHits);

    std::unique_ptr<std::vector<recob::Hit>> FinalHits
      (new std::vector<recob::Hit>(std::move(fCCAlg.YieldHits())));

    art::ServiceHandle<geo::Geometry const> geo;

    // shcol contains the hit collection
    // and its associations to wires and raw digits;
    // we get the association to raw digits through wire associations
    recob::HitCollectionAssociator shcol(evt, fCalDataModuleLabel, true);
    std::vector<recob::Cluster> sccol;
    std::vector<recob::Vertex> sv3col;

    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> >
      hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>>
      cv_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);


    std::vector<ClusterCrawlerAlg::ClusterStore> const& tcl = fCCAlg.GetClusters();

    std::vector<short> const& inClus = fCCAlg.GetinClus();

    // Consistency check
    for(unsigned int icl = 0; icl < tcl.size(); ++icl) {
      ClusterCrawlerAlg::ClusterStore const& clstr = tcl[icl];
      if(clstr.ID < 0) continue;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
      unsigned short plane = planeID.Plane;
      for(unsigned short ii = 0; ii < clstr.tclhits.size(); ++ii) {
        unsigned int iht = clstr.tclhits[ii];
        recob::Hit const& theHit = FinalHits->at(iht);
        if(theHit.WireID().Plane != plane) {
          std::cout<<"CC: cluster-hit plane mis-match "<<theHit.WireID().Plane<<" "<<plane
                   <<" in cluster "<<clstr.ID<<" WT "<<clstr.BeginWir<<":"<<(int)clstr.BeginTim<<"\n";
          return;
        }
        if(inClus[iht] != clstr.ID) {
          std::cout << "CC: InClus mis-match " << inClus[iht]
                    << " ID " << clstr.ID << " in cluster " << icl << "\n";
          return;
        }
      } // ii
    } // icl

    // make 3D vertices
    std::vector<ClusterCrawlerAlg::Vtx3Store> const& Vertices
      = fCCAlg.GetVertices();

    double xyz[3] = {0, 0, 0};
    unsigned int vtxID = 0, end;
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
    for(unsigned int icl = 0; icl < tcl.size(); ++icl) {
      ClusterCrawlerAlg::ClusterStore const& clstr = tcl[icl];
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
      // make the cluster - endpoint associations
      if(clstr.BeginVtx >= 0) {
        end = 0;
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
    evt.put(std::move(v3col));
    evt.put(std::move(cv_assn));

  } // produce
} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawler)

}

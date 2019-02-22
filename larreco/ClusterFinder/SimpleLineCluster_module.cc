////////////////////////////////////////////////////////////////////////
// Class:       SimpleLineCluster
// Plugin Type: producer (art v2_11_02)
// File:        SimpleLineCluster_module.cc
//
// Generated at Fri Jun 22 19:50:13 2018 by Tingjun Yang using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


//LArSoft includes
#include "larreco/RecoAlg/CCHitFinderAlg.h"
#include "larreco/RecoAlg/ClusterCrawlerAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include <memory>

namespace cluster {
  class SimpleLineCluster;
}


class cluster::SimpleLineCluster : public art::EDProducer {
public:
  explicit SimpleLineCluster(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleLineCluster(SimpleLineCluster const &) = delete;
  SimpleLineCluster(SimpleLineCluster &&) = delete;
  SimpleLineCluster & operator = (SimpleLineCluster const &) = delete;
  SimpleLineCluster & operator = (SimpleLineCluster &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::unique_ptr<ClusterCrawlerAlg> fCCAlg; // define ClusterCrawlerAlg object
  
  art::InputTag fHitFinderLabel; ///< label of module producing input hits
  art::InputTag fPFParticleLabel; ///< label of module producing input pfparticles
  
};


cluster::SimpleLineCluster::SimpleLineCluster(fhicl::ParameterSet const & p)
  : EDProducer{p}
  , fHitFinderLabel(p.get<art::InputTag>("HitFinderModuleLabel"))
  , fPFParticleLabel(p.get<art::InputTag>("PFParticleModuleLabel"))
{
  // this trick avoids double configuration on construction
  if (fCCAlg)
    fCCAlg->reconfigure(p.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
  else {
    fCCAlg.reset(new ClusterCrawlerAlg
                 (p.get< fhicl::ParameterSet >("ClusterCrawlerAlg")));
  }

  recob::HitCollectionAssociator::declare_products(*this);
  produces< std::vector<recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
}

void cluster::SimpleLineCluster::produce(art::Event & evt)
{
  //Retrieve data products
  auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitFinderLabel);
  auto clusHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fPFParticleLabel);
  auto pfpsHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(fPFParticleLabel);
  
  // all hits in the collection
  std::vector< art::Ptr<recob::Hit> > hits;
  art::fill_ptr_vector(hits, hitsHandle);
  
  // all pfps in the collection
  std::vector< art::Ptr<recob::PFParticle> > pfps;
  art::fill_ptr_vector(pfps, pfpsHandle);

  art::FindManyP < recob::Cluster > cluFromPfp(pfpsHandle, evt, fPFParticleLabel);
  art::FindManyP < recob::Hit > hitFromClu(clusHandle, evt, fPFParticleLabel);

  // also get the associated wires and raw digits;
  // we assume they have been created by the same module as the hits
  art::FindOneP<raw::RawDigit> channelHitRawDigits(hitsHandle, evt, fHitFinderLabel);
  art::FindOneP<recob::Wire>   channelHitWires    (hitsHandle, evt, fHitFinderLabel);

  std::vector<recob::Hit> allHits;
  std::vector<ClusterCrawlerAlg::ClusterStore> allClusters;
  std::vector<unsigned int> clusterhitindex;

  for (size_t i = 0; i < pfps.size(); ++i){
    std::vector<recob::Hit> pfphits;
    auto& clus_in_pfp = cluFromPfp.at(i);
    for (auto& clu : clus_in_pfp){
      auto& hits_in_clu = hitFromClu.at(clu.key());
      for (auto& hit : hits_in_clu) pfphits.push_back(*hit);
    }
    fCCAlg->RunCrawler(pfphits);

    //copy hits
    size_t nhits = allHits.size();
    for (auto & hit : fCCAlg->GetHits()){
      allHits.push_back(hit);
    }

    auto & Clusters = fCCAlg->GetClusters();
    for(auto &clstr : Clusters){
      allClusters.push_back(clstr);
      clusterhitindex.push_back(nhits);
    }

    fCCAlg->ClearResults();
  }
  

  std::vector<recob::Cluster> sccol;

  std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > 
    hc_assn(new art::Assns<recob::Cluster, recob::Hit>);

  // make the clusters and associations
  float sumChg, sumADC;
  unsigned int clsID = 0, nclhits;
  for(unsigned int icl = 0; icl < allClusters.size(); ++icl) {
    ClusterCrawlerAlg::ClusterStore const& clstr = allClusters[icl];
    if(clstr.ID < 0) continue;
    ++clsID;
    sumChg = 0;
    sumADC = 0;
    geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
    //unsigned short plane = planeID.Plane;
    nclhits = clstr.tclhits.size();
    std::vector<unsigned int> clsHitIndices;
    // correct the hit indices to refer to the valid hits that were just added
    for(unsigned int itt = 0; itt < nclhits; ++itt) {
      unsigned int iht = clstr.tclhits[itt];
      recob::Hit const& hit = allHits[iht];
      sumChg += hit.Integral();
      sumADC += hit.SummedADC();
    } // itt
      // get the wire, plane from a hit
    unsigned int iht = clstr.tclhits[0];
    
    geo::View_t view = allHits[iht].View();
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
    std::vector<size_t> indices;
    for (auto & ih : clstr.tclhits){
      indices.push_back(ih+clusterhitindex[icl]);
    }
    if(!util::CreateAssn(
                         *this, evt, sccol, allHits, *hc_assn, indices)
       )
      {
        throw art::Exception(art::errors::ProductRegistrationFailure)
          <<"Failed to associate hit "<<iht<<" with cluster "<<icl;
      } // exception
  }

  // convert cluster vector to unique_ptrs
  std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));

  std::unique_ptr<std::vector<recob::Hit>> FinalHits
    (new std::vector<recob::Hit>(std::move(allHits)));

  recob::HitRefinerAssociator shcol(*this, evt, fHitFinderLabel, 
                                    channelHitWires.isValid(), 
                                    channelHitRawDigits.isValid());

  shcol.use_hits(std::move(FinalHits));

  // move the hit collection and the associations into the event:
  shcol.put_into(evt);
  evt.put(std::move(ccol));
  evt.put(std::move(hc_assn));
  

}

DEFINE_ART_MODULE(cluster::SimpleLineCluster)

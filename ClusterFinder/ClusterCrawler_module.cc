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
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>

//LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "Utilities/AssociationUtil.h"
// #include "Filters/ChannelFilter.h"
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
#include "RecoAlg/CCHitRefinerAlg.h"


namespace cluster {
  class ClusterCrawler;
}

class cluster::ClusterCrawler : public art::EDProducer {

  public:
    explicit ClusterCrawler(fhicl::ParameterSet const & pset);
    virtual ~ClusterCrawler();

    void reconfigure(fhicl::ParameterSet const & pset) override;
    void produce(art::Event & evt) override;
    void beginJob();

  private:
    CCHitFinderAlg fCCHFAlg; // define CCHitFinderAlg object
    ClusterCrawlerAlg fCCAlg; // define ClusterCrawlerAlg object
    CCHitRefinerAlg fCCHRAlg; // define CCHitRefinerAlg object    
};


namespace cluster {

  ClusterCrawler::ClusterCrawler(fhicl::ParameterSet const& pset) :
    fCCHFAlg(pset.get< fhicl::ParameterSet >("CCHitFinderAlg" )),
    fCCAlg(  pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg")),
    fCCHRAlg(pset.get< fhicl::ParameterSet >("CCHitRefinerAlg" ))
  {  
    this->reconfigure(pset);
    produces< std::vector<recob::Hit> >();
    produces< std::vector<recob::Cluster> >();  
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< std::vector<recob::EndPoint2D> >();
  }

  ClusterCrawler::~ClusterCrawler()
  {
  }

  void ClusterCrawler::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCCAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
    fCCHFAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitFinderAlg"));
    fCCHRAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitRefinerAlg"));
  }
  
  void ClusterCrawler::beginJob(){
  }
  
  void ClusterCrawler::produce(art::Event & evt)
  {

    // find hits in all planes
    fCCHFAlg.RunCCHitFinder(evt);

    // look for clusters in all planes
    fCCAlg.RunCrawler(fCCHFAlg.allhits);
  
    // refine the hits, set cluster Begin/End
    fCCHRAlg.RunCCHitRefiner(fCCHFAlg.allhits, fCCHFAlg.hitcuts, 
                               fCCAlg.tcl, fCCAlg.vtx, fCCAlg);

    art::ServiceHandle<geo::Geometry> geo;
    
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::vector<recob::Hit> shcol;
    std::vector<recob::Cluster> sccol;
    std::vector<recob::EndPoint2D> svcol;

    // put clusters and hits into std::vectors
    unsigned short nclus = 0;
    unsigned short hitcnt = 0;
    for(size_t icl = 0; icl < fCCAlg.tcl.size(); ++icl) {
      ClusterCrawlerAlg::ClusterStore clstr = fCCAlg.tcl[icl];
      if(clstr.ID < 0) continue;
      // start cluster numbering at 1
      ++nclus;
      // make the hits on this cluster
      double totalQ = 0.;
      unsigned short firsthit = hitcnt;
      for(unsigned short itt = 0; itt < clstr.tclhits.size(); ++itt) {
        unsigned short iht = clstr.tclhits[itt];
        if(iht > fCCHFAlg.allhits.size() - 1) {
          mf::LogError("ClusterCrawler")<<"Bad hit index "<<iht;
          return;
        }
        CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
        // consistency check
        if(theHit.InClus != clstr.ID) {
          mf::LogError("ClusterCrawler")<<"Using bad hit in cluster "<<clstr.ID
            <<" hit ID "<<iht<<" InClus "<<theHit.InClus;
          return;
        }
        art::Ptr<recob::Wire> theWire = theHit.Wire;
        uint32_t channel = theWire->Channel();
        // get the Wire ID from the channel
        std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
        if(!wids[0].isValid) {
          mf::LogError("ClusterCrawler")<<"Invalid Wire ID "<<theWire<<" "<<channel;
          return;
        }
        recob::Hit hit(theHit.Wire,  wids[0],
              (double) theHit.Time - theHit.RMS, 0.,
              (double) theHit.Time + theHit.RMS, 0.,
              (double) theHit.Time, theHit.TimeErr,
              (double) theHit.Charge, theHit.ChargeErr,
              (double) theHit.Amplitude, theHit.AmplitudeErr,
              (int)    theHit.numHits, 
              (double) theHit.ChiDOF);
        shcol.push_back(hit);
        ++hitcnt;
      } // itt
      // get the view from a hit on the cluster
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[clstr.tclhits[0]];
      art::Ptr<recob::Wire> theWire = theHit.Wire;
      uint32_t channel = theWire->Channel();
      
      // create the recob::Cluster directly in the vector
      sccol.emplace_back((double)clstr.BeginWir, 0.,
                         (double)clstr.BeginTim, 0.,
                         (double)clstr.EndWir, 0.,
                         (double)clstr.EndTim, 0.,
                         (double)clstr.EndSlp, (double)clstr.EndSlpErr,
                         -999.,0.,
                         totalQ,
                         geo->View(channel),
                         nclus,
                         geo->ChannelToWire(channel).front().planeID() // get plane from channel geometry
                         );
      
      // associate the hits to this cluster
      util::CreateAssn(*this, evt, sccol, shcol, *hc_assn, firsthit, hitcnt);
    } // cluster iterator
//  std::cout<<"# clusters "<<nclus<<" # Hits in clusters "<<hitcnt;
    
    // make hits that are not associated with any cluster
    hitcnt = 0;
    for(unsigned short iht = 0; iht < fCCHFAlg.allhits.size(); ++iht) {
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
      // obsolete or used hit?
      if(theHit.InClus != 0) continue;
      ++hitcnt;
      art::Ptr<recob::Wire> theWire = theHit.Wire;
      uint32_t channel = theWire->Channel();
      // get the Wire ID from the channel
      std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
      if(!wids[0].isValid) {
        mf::LogError("ClusterCrawler")<<"Invalid Wire ID "<<theWire<<" "<<channel;
      }
      recob::Hit hit(theHit.Wire,  wids[0],
            (double) theHit.Time - theHit.RMS, 0.,
            (double) theHit.Time + theHit.RMS, 0.,
            (double) theHit.Time , theHit.TimeErr,
            (double) theHit.Charge , theHit.ChargeErr,
            (double) theHit.Amplitude , theHit.AmplitudeErr,
            (int)    theHit.numHits, 
            (double) theHit.ChiDOF);
      shcol.push_back(hit);
    }
//  std::cout<<" # Hits NOT in clusters "<<hitcnt;
    
    // convert to unique_ptrs
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);
    for(unsigned short iht = 0; iht < shcol.size(); ++iht) {
      hcol->push_back(shcol[iht]);
    }
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));

    // deal with cluster-EndPoint2D assns later (if necessary/desired)
    std::unique_ptr<std::vector<recob::EndPoint2D> > vcol(new std::vector<recob::EndPoint2D>);

    // make the vertex collection
//  std::cout<<" # vertices "<<fCCAlg.vtx.size()<<std::endl;
    for(unsigned short iv = 0; iv < fCCAlg.vtx.size(); iv++) {
      ClusterCrawlerAlg::VtxStore vtx = fCCAlg.vtx[iv];
      if(vtx.Wght <= 0) continue;
      if(vtx.CTP > 2) {
        mf::LogError("ClusterCrawler")<<"Bad vtx CTP "<<vtx.CTP;
        continue;
      }
      unsigned int cstat = vtx.CTP / 100;
      unsigned int tpc = (vtx.CTP - 100 * cstat) / 10;
      unsigned int plane = vtx.CTP - 100 * cstat - 10 * tpc;
      unsigned int wire = vtx.Wire;
      if(wire > geo->Nwires(plane) - 1) {
        mf::LogError("ClusterCrawler")<<"Bad vtx wire "<<wire<<" plane "
          <<plane<<" vtx # "<<iv;
        continue;
      }
      uint32_t channel = geo->PlaneWireToChannel(plane, wire, tpc, cstat);
      // get the Wire ID from the channel
      std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
      if(!wids[0].isValid) {
        mf::LogError("ClusterCrawler")<<"Invalid Wire ID "<<plane<<" "<<wire<<" "<<tpc<<" "<<cstat;
        continue;
      }
      recob::EndPoint2D myvtx((double)vtx.Time, wids[0], (double)vtx.Wght,
        (int)iv, geo->View(channel), 0.);
      vcol->push_back(myvtx);
    } // iv


    // clean up
    fCCHFAlg.allhits.clear();
    fCCAlg.tcl.clear();
    fCCAlg.vtx.clear();

    evt.put(std::move(hcol));
    evt.put(std::move(ccol));
    evt.put(std::move(hc_assn));
    evt.put(std::move(vcol));

  } // produce
} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawler)
  
} 


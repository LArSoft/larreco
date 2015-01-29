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
#include "art/Framework/Core/FindOneP.h"

#include <vector>
#include <memory> // std::move
#include <utility> // std::pair<>, std::unique_ptr<>

//LArSoft includes
#include "SimpleTypesAndConstants/geo_types.h"
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBaseArt/HitCreator.h"
#include "RecoBase/Vertex.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
//#include "RecoAlg/CCHitRefinerAlg.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"


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
//    CCHitRefinerAlg fCCHRAlg; // define CCHitRefinerAlg object    
};


namespace cluster {

  ClusterCrawler::ClusterCrawler(fhicl::ParameterSet const& pset) :
    fCCHFAlg(pset.get< fhicl::ParameterSet >("CCHitFinderAlg" )),
    fCCAlg(  pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"))
//    fCCHRAlg(pset.get< fhicl::ParameterSet >("CCHitRefinerAlg" ))
  {  
    this->reconfigure(pset);
    
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
    
    produces< std::vector<recob::Cluster> >();  
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< std::vector<recob::Vertex> >();
  }

  ClusterCrawler::~ClusterCrawler()
  {
  }

  void ClusterCrawler::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCCAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
    fCCHFAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitFinderAlg"));
//    fCCHRAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitRefinerAlg"));
  }
  
  void ClusterCrawler::beginJob(){
  }
  
  // used for sorting clusters by length
  typedef std::pair<unsigned int, unsigned int> mypair;
  bool SortByLen(const mypair& L, const mypair& R) {return (L.first > R.first);}
  
  void ClusterCrawler::produce(art::Event & evt)
  {

    // find hits in all planes
    fCCHFAlg.RunCCHitFinder(evt);

    // look for clusters in all planes
    fCCAlg.RunCrawler(fCCHFAlg.allhits);
  
/*
    // refine the hits, set cluster Begin/End
    fCCHRAlg.RunCCHitRefiner(fCCHFAlg.allhits, fCCHFAlg.hitcuts, 
                               fCCAlg.tcl, fCCAlg.vtx, fCCAlg);
*/

    art::ServiceHandle<geo::Geometry> geo;
    
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    // shcol contains the hit collection
    // and its associations to wires and raw digits
    recob::HitCollectionCreator shcol(*this, evt);
    std::vector<recob::Cluster> sccol;

    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here; configuration would happen here,
    // but we are using the default configuration for that algorithm
    ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;
    
    // learn from the algorithm which wires it used,
    // and get both wires and their association with raw digits
    std::string WireCreatorModuleLabel = fCCHFAlg.CalDataModuleLabel();
    art::ValidHandle< std::vector<recob::Wire>> wireVecHandle
     = evt.getValidHandle<std::vector<recob::Wire>>(WireCreatorModuleLabel);
    art::FindOneP<raw::RawDigit> WireToRawDigit
      (wireVecHandle, evt, WireCreatorModuleLabel);
    
    // put clusters and hits into std::vectors
    unsigned short nclus = 0;
    unsigned short hitcnt = 0;
    unsigned short icl, firsthit, itt, iht;
    double qtot;
    for(icl = 0; icl < fCCAlg.tcl.size(); ++icl) {
      ClusterCrawlerAlg::ClusterStore clstr = fCCAlg.tcl[icl];
      if(clstr.ID < 0) continue;
      // start cluster numbering at 1
      ++nclus;
      qtot = 0;
      // make the hits on this cluster
      firsthit = hitcnt;
      for(itt = 0; itt < clstr.tclhits.size(); ++itt) {
        iht = clstr.tclhits[itt];
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
        art::Ptr<recob::Wire> const& theWire = theHit.Wire;
        art::Ptr<raw::RawDigit> const& theRawDigit
          = WireToRawDigit.at(theWire.key());
        raw::ChannelID_t channel = theWire->Channel();
        // get the Wire ID from the channel
        std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
        if(!wids[0].isValid) {
          mf::LogError("ClusterCrawler")<<"Invalid Wire ID "<<theWire<<" "<<channel;
          return;
        }
        recob::HitCreator hit(
          *theWire,                  // wire reference
          wids[0],                   // wire ID
          theHit.LoTime,             // start_tick 
          theHit.HiTime,             // end_tick 
          theHit.RMS,                // rms
          theHit.Time,               // peak_time
          theHit.TimeErr,            // sigma_peak_time
          theHit.Amplitude,          // peak_amplitude
          theHit.AmplitudeErr,       // sigma_peak_amplitude
          theHit.Charge,             // hit_integral
          theHit.ChargeErr,          // hit_sigma_integral
          theHit.ADCSum,             // summed ADC
          theHit.numHits,            // multiplicity
          iht - theHit.LoHitID,      // local_index
          theHit.ChiDOF,             // goodness_of_fit
          theHit.DOF,                // dof
          std::vector<float>()       // signal FIXME
          );
        shcol.emplace_back(hit.move(), theWire, theRawDigit);
        ++hitcnt;
        qtot += theHit.Charge;
      } // itt
      // get the view from a hit on the cluster
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[clstr.tclhits[0]];
      art::Ptr<recob::Wire> theWire = theHit.Wire;
      
      // Stuff 2D vertex info into unused Cluster variables to help
      // associate the Begin cluster - vertex 2D and End cluster - vertex 2D
      double clBeginEP2Index = -1;
      double clBeginVtxIndex = -1;
      double clEndEP2Index = -1;
      double clEndVtxIndex = -1;
      short iv2 = 0;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
      unsigned short cPlane = planeID.Plane;
      if(!planeID.isValid) {
        mf::LogError("ClusterCrawler")<<"Invalid PlaneID from CTP "<<clstr.CTP;
        return;
      } // !planeID.isValid
      if(clstr.BeginVtx >= 0) {
        iv2 = clstr.BeginVtx;
        clBeginEP2Index = iv2;
        // See if this 2D vertex is associated with a 3D vertex
        for(unsigned short iv3 = 0; iv3 < fCCAlg.vtx3.size(); ++iv3) {
          // ignore incomplete vertices
          if(fCCAlg.vtx3[iv3].Ptr2D[0] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[1] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[2] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[cPlane] == iv2) {
            clBeginVtxIndex = iv3;
            break;
          }
        } // iv3
      } // clstr.BeginVtx >= 0
      if(clstr.EndVtx >= 0) {
        iv2 = clstr.EndVtx;
        clEndEP2Index = iv2;
        for(unsigned short iv3 = 0; iv3 < fCCAlg.vtx3.size(); ++iv3) {
          // ignore incomplete vertices
          if(fCCAlg.vtx3[iv3].Ptr2D[0] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[1] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[2] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[cPlane] == iv2) {
            clEndVtxIndex = iv3;
            break;
          }
        } // iv3
      } // clstr.EndVtx >= 0
      
      // feed the algorithm with all the cluster hits
      // (they are still in the HitCollectionCreator, but we can peek at them)
      ClusterParamAlgo.ImportHits(shcol.peek());
      
      // create the recob::Cluster directly in the vector;
      // the value of some of the quantities is defined by the track-like nature
      // imposed on the cluster by the algorithm.
      sccol.emplace_back(
      /* FIXME this is the old code
                         (double)clstr.BeginWir, clBeginEP2Index,
                         (double)clstr.BeginTim, clBeginVtxIndex,
                         (double)clstr.EndWir, clEndEP2Index,
                         (double)clstr.EndTim, clEndVtxIndex,
                         (double)clstr.EndSlp, (double)clstr.BeginSlp,
                         (double)clstr.BeginChg, (double)clstr.EndChg,
                         qtot,
                         geo->View(channel),
                         nclus,
                         planeID
//                         ClusterCrawlerAlg::DecodeCTP(clstr.CTP)
        */
        (float) clstr.BeginWir,                    // start_wire
        0.5,                                       // sigma_start_wire FIXME
        clstr.BeginTim,                            // start_time
        0.0,                                       // sigma_start_tick FIXME
        ClusterParamAlgo.StartCharge().value(),    // start_charge     FIXME
        ClusterParamAlgo.StartAngle().value(),     // start_angle      FIXME
        0.,                                        // start_opening (by definition)
        (float) clstr.EndWir,                      // end_wire
        0.5,                                       // sigma_end_wire   FIXME
        clstr.EndTim,                              // end_time
        0.0,                                       // sigma_end_tick   FIXME
        ClusterParamAlgo.EndCharge().value(),      // end_charge       FIXME
        ClusterParamAlgo.EndAngle().value(),       // end_angle        FIXME
        0.,                                        // end_opening (by definition)
        ClusterParamAlgo.Integral().value(),       // integral
        ClusterParamAlgo.IntegralStdDev().value(), // integral_stddev
        ClusterParamAlgo.SummedADC().value(),      // summedADC
        ClusterParamAlgo.SummedADCStdDev().value(),// summedADC_stddev
        ClusterParamAlgo.NHits(),                  // n_hits
        ClusterParamAlgo.NWiresOverNHits(),        // wires_over_hits
        0.,                                        // width (by definition)
        nclus,                                     // ID
        theWire->View(),                           // view
        planeID,                                   // planeID
        recob::Cluster::Sentry                     // sentry
        );
      
      // associate the hits to this cluster
      util::CreateAssn(*this, evt, sccol, shcol.peek(), *hc_assn, firsthit, hitcnt);
    } // cluster iterator
    
    // make hits that are not associated with any cluster
    hitcnt = 0;
    for(iht = 0; iht < fCCHFAlg.allhits.size(); ++iht) {
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
      // obsolete or used hit?
      if(theHit.InClus != 0) continue;
      ++hitcnt;
      art::Ptr<recob::Wire> theWire = theHit.Wire;
      art::Ptr<raw::RawDigit> const& theRawDigit
        = WireToRawDigit.at(theWire.key());
      raw::ChannelID_t channel = theWire->Channel();
      // get the Wire ID from the channel
      std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
      if(!wids[0].isValid) {
        mf::LogError("ClusterCrawler")<<"Invalid Wire ID "<<theWire<<" "<<channel;
      }
      recob::HitCreator hit(
        *theWire,                  // wire reference
        wids[0],                   // wire ID
        theHit.LoTime,             // start_tick
        theHit.HiTime,             // end_tick
        theHit.RMS,                // rms
        theHit.Time,               // peak_time
        theHit.TimeErr,            // sigma_peak_time
        theHit.Amplitude,          // peak_amplitude
        theHit.AmplitudeErr,       // sigma_peak_amplitude
        theHit.Charge,             // hit_integral
        theHit.ChargeErr,          // hit_sigma_integral
        theHit.ADCSum,             // summedADC FIXME
        theHit.numHits,            // multiplicity
        iht - theHit.LoHitID,      // local_index FIXME
        theHit.ChiDOF,             // goodness_of_fit
        theHit.DOF,                // dof
        std::vector<float>()       // signal FIXME
        );
      shcol.emplace_back(hit.move(), theWire, theRawDigit);
    } // for unassociated hits
    
    // convert to unique_ptrs
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));

    // 2D and 3D vertex collections
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>);
    std::unique_ptr<std::vector<recob::Vertex> > v3col(new std::vector<recob::Vertex>);

    // make the 2D vertex collection
    for(unsigned short iv = 0; iv < fCCAlg.vtx.size(); iv++) {
      ClusterCrawlerAlg::VtxStore aVtx = fCCAlg.vtx[iv];
      if(aVtx.Wght <= 0) continue;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(aVtx.CTP);
      unsigned int wire = (0.5 + aVtx.Wire);
      if(wire > geo->Nwires(planeID.Plane) - 1) {
        mf::LogError("ClusterCrawler")<<"Bad vtx wire "<<wire<<" plane "
          <<planeID.Plane<<" vtx # "<<iv;
        continue;
      }
      raw::ChannelID_t channel = geo->PlaneWireToChannel(planeID.Plane, wire, planeID.TPC, planeID.Cryostat);
      // get the Wire ID from the channel
      std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
      if(!wids[0].isValid) {
        mf::LogError("ClusterCrawler")<<"Invalid Wire ID "<<planeID.Plane<<" "<<wire<<" "<<planeID.TPC<<" "<<planeID.Cryostat;
        continue;
      }
      // stuff the (float) wire coordinate into the (not-useful) endpoint charge
      recob::EndPoint2D myvtx((double)aVtx.Time, wids[0], (double)aVtx.Wght,
        (int)iv, geo->View(channel), aVtx.Wire);
      v2col->push_back(myvtx);
    } // iv
    
    // make the 3D vertex collection
    double xyz[3] = {0, 0, 0};
    for(unsigned short iv = 0; iv < fCCAlg.vtx3.size(); iv++) {
      ClusterCrawlerAlg::Vtx3Store vtx3 = fCCAlg.vtx3[iv];
      // ignore incomplete vertices
      if(vtx3.Ptr2D[0] < 0) continue;
      if(vtx3.Ptr2D[1] < 0) continue;
      if(vtx3.Ptr2D[2] < 0) continue;
      xyz[0] = vtx3.X;
      xyz[1] = vtx3.Y;
      xyz[2] = vtx3.Z;
      recob::Vertex myvtx(xyz, iv);
      v3col->push_back(myvtx);
    }

    // clean up
    fCCHFAlg.allhits.clear();
    fCCAlg.tcl.clear();
    fCCAlg.vtx.clear();
    fCCAlg.vtx3.clear();

    // move the hit collection and the associations into the event:
    shcol.put_into(evt);
    evt.put(std::move(ccol));
    evt.put(std::move(hc_assn));
    evt.put(std::move(v2col));
    evt.put(std::move(v3col));

  } // produce
} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawler)
  
} 


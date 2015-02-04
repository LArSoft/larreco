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
};


namespace cluster {

  ClusterCrawler::ClusterCrawler(fhicl::ParameterSet const& pset) :
    fCCHFAlg(pset.get< fhicl::ParameterSet >("CCHitFinderAlg" )),
    fCCAlg(  pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"))
  {  
    this->reconfigure(pset);
    
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
    
    produces< std::vector<recob::Cluster> >();  
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::EndPoint2D, recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Vertex, recob::Cluster> >();
  }

  ClusterCrawler::~ClusterCrawler()
  {
  }

  void ClusterCrawler::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCCAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
    fCCHFAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitFinderAlg"));
  }
  
  void ClusterCrawler::beginJob(){
  }
  
  void ClusterCrawler::produce(art::Event & evt)
  {

    // find hits in all planes
    fCCHFAlg.RunCCHitFinder(evt);

    // look for clusters in all planes
    fCCAlg.RunCrawler(fCCHFAlg.allhits);
  
    art::ServiceHandle<geo::Geometry> geo;
    
    // shcol contains the hit collection
    // and its associations to wires and raw digits
    recob::HitCollectionCreator shcol(*this, evt);
    std::vector<recob::Cluster> sccol;
    std::vector<recob::EndPoint2D> sv2col;
    std::vector<recob::Vertex> sv3col;

    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > 
        hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>> 
       cep_assn(new art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>> 
        cv_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);
    
    // learn from the algorithm which wires it used,
    // and get both wires and their association with raw digits
    std::string WireCreatorModuleLabel = fCCHFAlg.CalDataModuleLabel();
    art::ValidHandle< std::vector<recob::Wire>> wireVecHandle
     = evt.getValidHandle<std::vector<recob::Wire>>(WireCreatorModuleLabel);
    art::FindOneP<raw::RawDigit> WireToRawDigit
      (wireVecHandle, evt, WireCreatorModuleLabel);

    // Temp vector for indexing valid (i.e. not obsolete) hits
    std::vector<int> vhit(fCCHFAlg.allhits.size(), -1);
    unsigned int icl, iht, itt, nvhit = 0;
    for(iht = 0; iht < fCCHFAlg.allhits.size(); ++iht) {
      if(fCCHFAlg.allhits[iht].InClus < 0) continue;
      vhit[iht] = nvhit;
      ++nvhit;
      // store the hit
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
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
      recob::Hit hit(
        channel,                   // channel
        theHit.LoTime,             // start_tick 
        theHit.HiTime,             // end_tick 
        theHit.Time,               // peak_time
        theHit.TimeErr,            // sigma_peak_time
        theHit.RMS,                // rms
        theHit.Amplitude,          // peak_amplitude
        theHit.AmplitudeErr,       // sigma_peak_amplitude
        theHit.Charge,             // hit_integral
        theHit.ADCSum,             // summed ADC
        0,                         // sigma summed ADC
        theHit.numHits,            // multiplicity
        iht - theHit.LoHitID,      // local_index FIXME LOHITID...
        theHit.ChiDOF,             // goodness_of_fit
        theHit.DOF,                // dof
        theWire->View(),           // view
        geo->SignalType(channel),  // ...
        wids[0]                    // wire ID
        );
      shcol.emplace_back(std::move(hit), theWire, theRawDigit);
    } // iht
    
    // make the endpoints
    // Temp vector for indexing valid (i.e. not obsolete) endpoints
    std::vector<int> vep2(fCCAlg.vtx.size(), -1);
    unsigned int ep2ID = 0;
    unsigned short end;
    for(unsigned int iep = 0; iep < fCCAlg.vtx.size(); iep++) {
      ClusterCrawlerAlg::VtxStore& aVtx = fCCAlg.vtx[iep];
      if(aVtx.Wght <= 0) continue;
      vep2[iep] = ep2ID;
      ++ep2ID;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(aVtx.CTP);
      unsigned int wire = (0.5 + aVtx.Wire);
      if(wire > geo->Nwires(planeID.Plane) - 1) {
        mf::LogError("ClusterCrawler")<<"Bad vtx wire "<<wire<<" plane "
          <<planeID.Plane<<" vtx # "<<iep;
        continue;
      }
      raw::ChannelID_t channel = geo->PlaneWireToChannel(planeID.Plane, wire, planeID.TPC, planeID.Cryostat);
      // get the Wire ID from the channel
      std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
      if(!wids[0].isValid) {
        mf::LogError("ClusterCrawler")<<"Invalid Wire ID "<<planeID.Plane
          <<" "<<wire<<" "<<planeID.TPC<<" "<<planeID.Cryostat;
        continue;
      }
      sv2col.emplace_back(
        (double)aVtx.Time,     // time
        wids[0],              // wire ID
        (double)aVtx.Wght,    // strength - whatever that is...
        (int)ep2ID,           // endpoint ID
        geo->View(channel),   // view
        (double)aVtx.Wire     // total charge --> wire number in double precision
      );
    } // iep
    // convert EndPoint2D vector to unique_ptrs
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(sv2col)));

    // make 3D vertices
    double xyz[3] = {0, 0, 0};
    unsigned int vtxID = 0;
    for(unsigned int iv = 0; iv < fCCAlg.vtx3.size(); iv++) {
      ClusterCrawlerAlg::Vtx3Store& vtx3 = fCCAlg.vtx3[iv];
      // ignore incomplete vertices
      if(vtx3.Ptr2D[0] < 0) continue;
      if(vtx3.Ptr2D[1] < 0) continue;
      if(vtx3.Ptr2D[2] < 0) continue;
      ++vtxID;
      xyz[0] = vtx3.X;
      xyz[1] = vtx3.Y;
      xyz[2] = vtx3.Z;
      sv3col.emplace_back(xyz, vtxID);
    } // iv
    // convert Vertex vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Vertex> > v3col(new std::vector<recob::Vertex>(std::move(sv3col)));
    
    // make the clusters and associations
    float sumChg, sumADC, wovrh;
    unsigned int clsID = 0, nclhits, iep;
    for(icl = 0; icl < fCCAlg.tcl.size(); ++icl) {
      ClusterCrawlerAlg::ClusterStore& clstr = fCCAlg.tcl[icl];
      if(clstr.ID < 0) continue;
      ++clsID;
      sumChg = 0;
      sumADC = 0;
      nclhits = clstr.tclhits.size();
      // correct the hit indices to refer to the valid hits that were just added
      for(itt = 0; itt < nclhits; ++itt) {
        iht = clstr.tclhits[itt];
        if(iht > fCCHFAlg.allhits.size() - 1 || vhit[iht] < 0) {
          mf::LogError("ClusterCrawler")<<"Bad hit index "<<iht;
          return;
        } // bad hit index
        // set iht to the valid hit index
        iht = vhit[iht];
        clstr.tclhits[itt] = iht;
        sumChg += fCCHFAlg.allhits[iht].Charge;
        sumADC += fCCHFAlg.allhits[iht].ADCSum;
      } // itt
      // get the wire, plane from a hit
      iht = clstr.tclhits[0];
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
      art::Ptr<recob::Wire> const& theWire = theHit.Wire;
      raw::ChannelID_t channel = theWire->Channel();
      std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
      wovrh = (float)(clstr.BeginWir - clstr.EndWir) / (float)nclhits;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
      unsigned short plane = planeID.Plane;
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
          wovrh,                  // wires over hits
          0,                      // width (0 for line-like clusters)
          clsID,                   // ID
          geo->View(channel),     // view
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
        iep = vep2[clstr.BeginVtx];
        end = 0;
        if(iep < 0) {
          mf::LogError("ClusterCrawler")<<"Invalid cluster -> endpoint association ";
          return;
        } // iep < 0
        if(!util::CreateAssnD(*this, evt, *cep_assn, clsID - 1, iep, end))
        {
          throw art::Exception(art::errors::InsertFailure)
            <<"Failed to associate cluster "<<icl<<" with endpoint";
        } // exception
        // See if this endpoint is associated with a 3D vertex
        unsigned short iv = 0;
        for(unsigned short iv3 = 0; iv3 < fCCAlg.vtx3.size(); ++iv3) {
          // ignore incomplete vertices
          if(fCCAlg.vtx3[iv3].Ptr2D[0] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[1] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[2] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[plane] == clstr.BeginVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, iv3, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++iv;
        } // iv3
      } // clstr.BeginVtx >= 0
      if(clstr.EndVtx >= 0) {
        iep = vep2[clstr.EndVtx];
        end = 1;
        if(iep < 0) {
          mf::LogError("ClusterCrawler")<<"Invalid cluster -> endpoint association ";
          return;
        } // iep < 0
        if(!util::CreateAssnD(*this, evt, *cep_assn, clsID - 1, iep, end))
        {
          throw art::Exception(art::errors::InsertFailure)
            <<"Failed to associate cluster "<<icl<<" with endpoint";
        } // exception
        // See if this endpoint is associated with a 3D vertex
        unsigned short iv = 0;
        for(unsigned short iv3 = 0; iv3 < fCCAlg.vtx3.size(); ++iv3) {
          // ignore incomplete vertices
          if(fCCAlg.vtx3[iv3].Ptr2D[0] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[1] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[2] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[plane] == clstr.EndVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, iv3, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++iv;
        } // iv3
      } // clstr.BeginVtx >= 0
    } // icl
    
    // convert cluster vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));

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
    evt.put(std::move(cep_assn));
    evt.put(std::move(v3col));
    evt.put(std::move(cv_assn));

  } // produce
} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawler)
  
} 


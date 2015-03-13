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
#include <algorithm> // std::max()
#include <functional> // std::mem_fn()
#include <memory> // std::move
#include <utility> // std::pair<>, std::unique_ptr<>
#include <limits> // std::numeric_limits<>

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
#include "RecoBaseArt/HitCreator.h" // recob::HitCollectionCreator
#include "RecoBase/Vertex.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/MakeIndex.h"
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
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
    std::string fCalDataModuleLabel; ///< label of module producing input wires
};


namespace cluster {

  ClusterCrawler::ClusterCrawler(fhicl::ParameterSet const& pset) :
    fCCHFAlg           (pset.get<fhicl::ParameterSet>("CCHitFinderAlg")),
    fCCAlg             (pset.get<fhicl::ParameterSet>("ClusterCrawlerAlg")),
    fCalDataModuleLabel(pset.get<std::string>("CalDataModuleLabel"))
  {
    this->reconfigure(pset);
    
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
    
    produces< std::vector<recob::Cluster> >();  
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();
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
    // fetch the wires needed by CCHitFinder

    // make this accessible to ClusterCrawler_module
    art::ValidHandle< std::vector<recob::Wire>> wireVecHandle
     = evt.getValidHandle<std::vector<recob::Wire>>(fCalDataModuleLabel);

    // find hits in all planes
    fCCHFAlg.RunCCHitFinder(*wireVecHandle);
    
    // extract the result of the algorithm (it's moved)
    std::vector<recob::Hit> FirstHits = fCCHFAlg.YieldHits();

    // look for clusters in all planes
    fCCAlg.RunCrawler(fCCHFAlg.allhits, FirstHits);
    
    // access to the algorithm results
    ClusterCrawlerAlg::HitInCluster_t const& HitInCluster
      = fCCAlg.GetHitInCluster();
    
    std::vector<recob::Hit> FinalHits = fCCAlg.YieldHits();
    
    art::ServiceHandle<geo::Geometry> geo;
    
    // shcol contains the hit collection
    // and its associations to wires and raw digits
    recob::HitCollectionCreator shcol(*this, evt);
    std::vector<recob::Cluster> sccol;
    std::vector<recob::Vertex> sv3col;

    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > 
        hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>> 
        cv_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);
    
    // learn from the algorithm which wires it used,
    // and get both wires and their association with raw digits
    art::FindOneP<raw::RawDigit> WireToRawDigit
      (wireVecHandle, evt, fCalDataModuleLabel);

    unsigned int icl, iht, itt;

// Consistency check
  for(icl = 0; icl < fCCAlg.tcl.size(); ++icl) {
    ClusterCrawlerAlg::ClusterStore& clstr = fCCAlg.tcl[icl];
    if(clstr.ID < 0) continue;
    geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
    unsigned short plane = planeID.Plane;
    for(unsigned short ii = 0; ii < clstr.tclhits.size(); ++ii) {
      iht = clstr.tclhits[ii];
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
      if(theHit.WirID.Plane != plane) {
        std::cout<<"CC: cluster-hit plane mis-match "<<theHit.WirID.Plane<<" "<<plane
        <<" in cluster "<<clstr.ID<<" WT "<<clstr.BeginWir<<":"<<(int)clstr.BeginTim<<"\n";
        return;
      }
      if(HitInCluster[iht] != clstr.ID) {
        std::cout << "CC: InClus mis-match " << HitInCluster[iht]
          << " ID " << clstr.ID << " in cluster " << icl << "\n";
        return;
      }
    } // ii
  } // icl

    // fill a map of wire index vs. channel number
    std::vector<size_t> WireMap
      = util::MakeIndex(*wireVecHandle, std::mem_fn(&recob::Wire::Channel));
    
    // Temp vector for indexing valid (i.e. not obsolete) hits
    std::vector<int> vhit(fCCHFAlg.allhits.size(), -1);
    unsigned int locIndex, nvhit = 0;
    for(iht = 0; iht < fCCHFAlg.allhits.size(); ++iht) {
      if(!HitInCluster.isPresent(iht)) continue;
      vhit[iht] = nvhit;
      ++nvhit;
      // store the hit
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
      raw::ChannelID_t channel = theHit.Wire->Channel();
      size_t WireKey = WireMap.at(channel);
      if (WireKey == std::numeric_limits<size_t>::max()) {
        throw art::Exception(art::errors::LogicError)
          << "There is no wire for channel " << channel << ", yet hit #"
          << iht << " was created from it!";
      }
      const art::Ptr<recob::Wire> theWire(wireVecHandle, WireKey);
      art::Ptr<raw::RawDigit> const& theRawDigit = WireToRawDigit.at(theWire.key());
      if (theRawDigit->Channel() != channel) {
        throw art::Exception(art::errors::LogicError)
          << "I am trying to associate hit #" << iht << " on channel "
          << channel << " with a raw digit on channel "
          << theRawDigit->Channel() << "!!!";
      }
      if (theWire->Channel() != channel) {
        throw art::Exception(art::errors::LogicError)
          << "I am trying to associate hit #" << iht << " on channel "
          << channel << " with a wire on channel "
          << theWire->Channel() << "!!!";
      }
      // find the local index in a hit multiplet
      locIndex = iht - theHit.LoHitID;
      recob::Hit hit(
        channel,                   // channel
        theHit.LoTime,             // start_tick 
        theHit.HiTime,             // end_tick 
        theHit.Time,               // peak_time
        theHit.TimeErr,            // sigma_peak_time
        theHit.RMS,                // rms
        theHit.Amplitude,          // peak_amplitude
        theHit.AmplitudeErr,       // sigma_peak_amplitude
        theHit.ADCSum,             // summed ADC
        theHit.Charge,             // hit_integral
        theHit.ChargeErr,          // sigma hit_integral
        theHit.numHits,            // multiplicity
        locIndex,                  // local_index
        theHit.ChiDOF,             // goodness_of_fit
        theHit.DOF,                // dof
        theWire->View(),           // view
        geo->SignalType(channel),  // ...
        theHit.WirID               // wire ID
        );
      
      shcol.emplace_back(std::move(hit), theWire, theRawDigit);
    } // iht

    // make 3D vertices
    double xyz[3] = {0, 0, 0};
    unsigned int vtxID = 0, end;
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
    float sumChg, sumADC;
    unsigned int clsID = 0, nclhits;
    for(icl = 0; icl < fCCAlg.tcl.size(); ++icl) {
      ClusterCrawlerAlg::ClusterStore& clstr = fCCAlg.tcl[icl];
      if(clstr.ID < 0) continue;
      ++clsID;
      sumChg = 0;
      sumADC = 0;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
      unsigned short plane = planeID.Plane;
      nclhits = clstr.tclhits.size();
      std::vector<unsigned int> clsHitIndices;
      // correct the hit indices to refer to the valid hits that were just added
      for(itt = 0; itt < nclhits; ++itt) {
        iht = clstr.tclhits[itt];
        if(iht > fCCHFAlg.allhits.size() - 1 || vhit[iht] < 0 || 
           vhit[iht] > (int)fCCHFAlg.allhits.size() - 1) 
          throw cet::exception("ClusterCrawler")<<"Bad hit index "<<iht;
        sumChg += fCCHFAlg.allhits[iht].Charge;
        sumADC += fCCHFAlg.allhits[iht].ADCSum;
        clsHitIndices.push_back(vhit[iht]);
      } // itt
      // get the wire, plane from a hit
      iht = clstr.tclhits[0];
      CCHitFinderAlg::CCHit& theHit = fCCHFAlg.allhits[iht];
      size_t WireKey = WireMap.at(theHit.Wire->Channel());
      const art::Ptr<recob::Wire> theWire(wireVecHandle, WireKey);
      raw::ChannelID_t channel = theWire->Channel();
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
          clsID,                   // ID
          geo->View(channel),     // view
          planeID,                // plane
          recob::Cluster::Sentry  // sentry
          );
      // make the cluster - hit association
      if(!util::CreateAssn(
        *this, evt, *hc_assn, sccol.size()-1, clsHitIndices.begin(), clsHitIndices.end())
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
        for(unsigned short iv3 = 0; iv3 < fCCAlg.vtx3.size(); ++iv3) {
          // ignore incomplete vertices
          if(fCCAlg.vtx3[iv3].Ptr2D[0] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[1] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[2] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[plane] == clstr.BeginVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with vertex";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // iv3
      } // clstr.BeginVtx >= 0
      if(clstr.EndVtx >= 0) {
        end = 1;
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(unsigned short iv3 = 0; iv3 < fCCAlg.vtx3.size(); ++iv3) {
          // ignore incomplete vertices
          if(fCCAlg.vtx3[iv3].Ptr2D[0] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[1] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[2] < 0) continue;
          if(fCCAlg.vtx3[iv3].Ptr2D[plane] == clstr.EndVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
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
    evt.put(std::move(v3col));
    evt.put(std::move(cv_assn));

  } // produce
} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawler)
  
} 


////////////////////////////////////////////////////////////////////////
//
//  CCTrackMaker
//
//  Make tracks using ClusterCrawler clusters and vertex info
//
//  baller@fnal.gov, October 2014
//  May 2015: Major re-write
//
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/TrackTrajectoryAlg.h"
#include "larreco/RecoAlg/VertexFitAlg.h"

struct CluLen {
  int index;
  int length;
};

bool greaterThan(CluLen c1, CluLen c2)
{
  return (c1.length > c2.length);
}
bool lessThan(CluLen c1, CluLen c2)
{
  return (c1.length < c2.length);
}

namespace trkf {

  class CCTrackMaker : public art::EDProducer {
  public:
    explicit CCTrackMaker(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    std::string fHitModuleLabel;
    std::string fClusterModuleLabel;
    std::string fVertexModuleLabel;

    // services
    art::ServiceHandle<geo::Geometry const> geom;
    geo::WireReadoutGeom const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

    TrackTrajectoryAlg fTrackTrajectoryAlg;
    VertexFitAlg fVertexFitAlg;

    // Track matching parameters
    unsigned short algIndex;
    std::vector<short> fMatchAlgs;
    std::vector<float> fXMatchErr;
    std::vector<float> fAngleMatchErr;
    std::vector<float> fChgAsymFactor;
    std::vector<float> fMatchMinLen;
    std::vector<bool> fMakeAlgTracks;

    // Cluster merging parameters
    float fMaxDAng;
    float fChainMaxdX;
    float fChainVtxAng;
    float fMergeChgAsym;
    float fMaxMergeError;
    float fMergeErrorCut;

    float fChgWindow;
    float fWirePitch;
    // cosmic ray tagging
    float fFiducialCut;
    float fDeltaRayCut;

    bool fMakePFPs;

    // vertex fitting
    unsigned short fNVtxTrkHitsFit;
    float fHitFitErrFac;

    // temp
    bool fuBCode;

    // debugging inputs
    short fDebugAlg;
    short fDebugPlane;
    short fDebugCluster;
    bool fPrintAllClusters;
    bool prt;

    // hit indexing info
    std::array<unsigned int, 3> firstWire;
    std::array<unsigned int, 3> lastWire;
    std::array<unsigned int, 3> firstHit;
    std::array<unsigned int, 3> lastHit;
    std::array<std::vector<std::pair<int, int>>, 3> WireHitRange;

    std::vector<art::Ptr<recob::Hit>> allhits;

    // Cluster parameters
    struct clPar {
      std::array<float, 2> Wire;       // Begin/End Wire
      std::array<float, 2> X;          // Begin/End X
      std::array<short, 2> Time;       // Begin/End Time
      std::array<float, 2> Slope;      // Begin/End slope
      std::array<float, 2> Charge;     // Begin/End charge
      std::array<float, 2> ChgNear;    // Charge near the cluster at each end
      std::array<float, 2> Angle;      // Begin/End angle (radians)
      std::array<short, 2> Dir;        // Product of end * slope
      std::array<short, 2> VtxIndex;   // Vertex index
      std::array<short, 2> mVtxIndex;  // "Maybe" Vertex index
      std::array<short, 2> BrkIndex;   // Broken cluster index
      std::array<float, 2> MergeError; // Broken cluster merge error (See MakeClusterChains)
      unsigned short EvtIndex;         // index of the cluster in clusterlist
      short InTrack;                   // cluster -> chain index
      unsigned short Length;           // cluster length (wires)
      float TotChg;                    // total charge of cluster (or series of clusters)
    };
    // vector of cluster parameters in each plane
    std::array<std::vector<clPar>, 3> cls;

    // cluster chain parameters
    struct ClsChainPar {
      std::array<float, 2> Wire;     // Begin/End Wire
      std::array<float, 2> X;        // Begin/End X
      std::array<float, 2> Time;     // Begin/End Time
      std::array<float, 2> Slope;    // Begin/End slope
      std::array<float, 2> Angle;    // Begin/End angle (radians)
      std::array<short, 2> VtxIndex; // Vertex index
      std::array<short, 2> Dir;
      std::array<float, 2> ChgNear;   // Charge near the cluster at each end
      std::array<short, 2> mBrkIndex; // a "Maybe" cluster index
      unsigned short Length;
      float TotChg;
      std::vector<unsigned short> ClsIndex;
      std::vector<unsigned short> Order;
      short InTrack; // cluster -> track ID (-1 if none, 0 if under construction)
    };
    // vector of cluster parameters in each plane
    std::array<std::vector<ClsChainPar>, 3> clsChain;

    // 3D Vertex info
    struct vtxPar {
      short ID;
      unsigned short EvtIndex;
      float X;
      float Y;
      float Z;
      std::array<unsigned short, 3> nClusInPln;
      bool Neutrino;
    };

    std::vector<vtxPar> vtx;

    // cluster indices assigned to one vertex. Filled in VtxMatch
    std::array<std::vector<unsigned short>, 3> vxCls;

    struct TrkPar {
      short ID;
      unsigned short Proc; // 1 = VtxMatch, 2 = ...
      std::array<std::vector<art::Ptr<recob::Hit>>, 3> TrkHits;
      std::vector<TVector3> TrjPos;
      std::vector<TVector3> TrjDir;
      std::array<short, 2> VtxIndex;
      std::vector<unsigned short> ClsEvtIndices;
      float Length;
      short ChgOrder;
      short MomID;
      std::array<bool, 2> EndInTPC;
      std::array<bool, 2> GoodEnd; // set true if there are hits in all planes at the end point
      std::vector<short> DtrID;
      short PDGCode;
    };
    std::vector<TrkPar> trk;

    // Array of pointers to hits in each plane for one track
    std::array<std::vector<art::Ptr<recob::Hit>>, 3> trkHits;
    // and for one seed
    std::array<std::vector<art::Ptr<recob::Hit>>, 3> seedHits;
    // relative charge normalization between planes
    std::array<float, 3> ChgNorm;

    // Vector of PFParticle -> track IDs. The first element will be
    // track ID = 0 indicating that this is a neutrino PFParticle and
    // there is no associated track
    std::vector<unsigned short> pfpToTrkID;

    // characterize the match between clusters in 2 or 3 planes
    struct MatchPars {
      std::array<short, 3> Cls;
      std::array<unsigned short, 3> End;
      std::array<float, 3> Chg;
      short Vtx;
      float dWir; // wire difference at the matching end
      float dAng; // angle difference at the matching end
      float dX;   // X difference
      float Err;  // Wire,Angle,Time match error
      short oVtx;
      float odWir; // wire difference at the other end
      float odAng; // angle difference at the other end
      float odX;   // time difference at the other end
      float dSP;   // space point difference
      float oErr;  // dAngle dX match error
    };
    // vector of many match combinations
    std::vector<MatchPars> matcomb;

    void PrintClusters(detinfo::DetectorPropertiesData const& detProp,
                       geo::TPCID const& tpcid) const;

    void PrintTracks() const;

    void MakeClusterChains(detinfo::DetectorPropertiesData const& detProp,
                           art::FindManyP<recob::Hit> const& fmCluHits,
                           geo::TPCID const& tpcid);
    float dXClTraj(art::FindManyP<recob::Hit> const& fmCluHits,
                   unsigned short ipl,
                   unsigned short icl1,
                   unsigned short end1);
    void FillChgNear(detinfo::DetectorPropertiesData const& detProp, geo::TPCID const& tpcid);
    void FillWireHitRange(geo::TPCID const& tpcid);

    // Find clusters that point to vertices but do not have a
    // cluster-vertex association made by ClusterCrawler
    void FindMaybeVertices(geo::TPCID const& tpcid);
    // match clusters associated with vertices
    void VtxMatch(detinfo::DetectorPropertiesData const& detProp,
                  art::FindManyP<recob::Hit> const& fmCluHits,
                  geo::TPCID const& tpcid);
    // match clusters in all planes
    void PlnMatch(detinfo::DetectorPropertiesData const& detProp, geo::TPCID const& tpcid);
    // match clusters in all planes
    void AngMatch(art::FindManyP<recob::Hit> const& fmCluHits);

    // Make the track/vertex and mother/daughter relationships
    void MakeFamily(geo::TPCID const& tpcid);
    void TagCosmics(geo::TPCID const& tpcid);

    void FitVertices(detinfo::DetectorPropertiesData const& detProp, geo::TPCID const& tpcid);

    // fill the end matching parameters in the MatchPars struct
    void FillEndMatch(detinfo::DetectorPropertiesData const& detProp,
                      geo::TPCID const& tpcid,
                      MatchPars& match);
    // 2D version
    void FillEndMatch2(MatchPars& match);

    float ChargeAsym(std::array<float, 3>& mChg);

    bool FindMissingCluster(unsigned short kpl,
                            short& kcl,
                            unsigned short& kend,
                            float kWir,
                            float kX,
                            float okWir,
                            float okX);

    bool DupMatch(MatchPars& match, unsigned short nplanes);

    void SortMatches(detinfo::DetectorPropertiesData const& detProp,
                     art::FindManyP<recob::Hit> const& fmCluHits,
                     unsigned short procCode,
                     geo::TPCID const& tpcid);

    // fill the trkHits array using information.
    void FillTrkHits(art::FindManyP<recob::Hit> const& fmCluHits,
                     unsigned short imat,
                     unsigned short nplanes);

    // store the track in the trk vector
    void StoreTrack(detinfo::DetectorPropertiesData const& detProp,
                    art::FindManyP<recob::Hit> const& fmCluHits,
                    unsigned short imat,
                    unsigned short procCode,
                    geo::TPCID const& tpcid);

    // returns the charge along the line between (wire1, time1) and (wire2, time2)
    float ChargeNear(geo::PlaneID const& planeid,
                     unsigned short wire1,
                     float time1,
                     unsigned short wire2,
                     float time2);

    // inflate cuts at large angle
    float AngleFactor(float slope);

  }; // class CCTrackMaker

  //-------------------------------------------------
  CCTrackMaker::CCTrackMaker(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fVertexModuleLabel = pset.get<std::string>("VertexModuleLabel");
    // track matching
    fMatchAlgs = pset.get<std::vector<short>>("MatchAlgs");
    fXMatchErr = pset.get<std::vector<float>>("XMatchErr");
    fAngleMatchErr = pset.get<std::vector<float>>("AngleMatchErr");
    fChgAsymFactor = pset.get<std::vector<float>>("ChgAsymFactor");
    fMatchMinLen = pset.get<std::vector<float>>("MatchMinLen");
    fMakeAlgTracks = pset.get<std::vector<bool>>("MakeAlgTracks");
    // Cluster merging
    fMaxDAng = pset.get<float>("MaxDAng");
    fChainMaxdX = pset.get<float>("ChainMaxdX");
    fChainVtxAng = pset.get<float>("ChainVtxAng");
    fMergeChgAsym = pset.get<float>("MergeChgAsym");
    // Cosmic ray tagging
    fFiducialCut = pset.get<float>("FiducialCut");
    fDeltaRayCut = pset.get<float>("DeltaRayCut");
    // make PFParticles
    fMakePFPs = pset.get<bool>("MakePFPs");
    // vertex fitting
    fNVtxTrkHitsFit = pset.get<unsigned short>("NVtxTrkHitsFit");
    fHitFitErrFac = pset.get<float>("HitFitErrFac");
    // uB code
    fuBCode = pset.get<bool>("uBCode");
    // debugging inputs
    fDebugAlg = pset.get<short>("DebugAlg");
    fDebugPlane = pset.get<short>("DebugPlane");
    fDebugCluster = pset.get<short>("DebugCluster");
    fPrintAllClusters = pset.get<bool>("PrintAllClusters");

    // Consistency check
    if (fMatchAlgs.size() > fXMatchErr.size() || fMatchAlgs.size() > fAngleMatchErr.size() ||
        fMatchAlgs.size() > fChgAsymFactor.size() || fMatchAlgs.size() > fMatchMinLen.size() ||
        fMatchAlgs.size() > fMakeAlgTracks.size()) {
      mf::LogError("CCTM") << "Incompatible fcl input vector sizes";
      return;
    }
    // Reality check
    for (unsigned short ii = 0; ii < fMatchAlgs.size(); ++ii) {
      if (fAngleMatchErr[ii] <= 0 || fXMatchErr[ii] <= 0) {
        mf::LogError("CCTM") << "Invalid matching parameters " << fAngleMatchErr[ii] << " "
                             << fXMatchErr[ii];
        return;
      }
    } // ii

    produces<std::vector<recob::PFParticle>>();
    produces<art::Assns<recob::PFParticle, recob::Track>>();
    produces<art::Assns<recob::PFParticle, recob::Cluster>>();
    produces<art::Assns<recob::PFParticle, recob::Seed>>();
    produces<art::Assns<recob::PFParticle, recob::Vertex>>();
    produces<std::vector<recob::Vertex>>();
    produces<std::vector<recob::Track>>();
    produces<art::Assns<recob::Track, recob::Hit>>();
    produces<std::vector<recob::Seed>>();
  }

  //------------------------------------------------------------------------------------//
  void CCTrackMaker::produce(art::Event& evt)
  {
    fWirePitch = wireReadoutGeom.Plane({0, 0, 0}).WirePitch();

    fChgWindow = 40; // window (ticks) for identifying shower-like clusters

    auto tcol = std::make_unique<std::vector<recob::Track>>();
    auto thassn = std::make_unique<art::Assns<recob::Track, recob::Hit>>();
    auto vcol = std::make_unique<std::vector<recob::Vertex>>();
    auto pcol = std::make_unique<std::vector<recob::PFParticle>>();
    auto ptassn = std::make_unique<art::Assns<recob::PFParticle, recob::Track>>();
    auto pcassn = std::make_unique<art::Assns<recob::PFParticle, recob::Cluster>>();
    auto psassn = std::make_unique<art::Assns<recob::PFParticle, recob::Seed>>();
    auto pvassn = std::make_unique<art::Assns<recob::PFParticle, recob::Vertex>>();

    // seed collection
    auto scol = std::make_unique<std::vector<recob::Seed>>();
    auto shassn = std::make_unique<art::Assns<recob::Seed, recob::Hit>>();

    // all hits
    art::Handle<std::vector<recob::Hit>> allhitsListHandle;
    // cluster list
    art::Handle<std::vector<recob::Cluster>> clusterListHandle;
    std::vector<art::Ptr<recob::Cluster>> clusterlist;
    // ClusterCrawler Vertices
    art::Handle<std::vector<recob::Vertex>> VtxListHandle;
    std::vector<art::Ptr<recob::Vertex>> vtxlist;

    // get Hits
    allhits.clear();
    if (evt.getByLabel(fHitModuleLabel, allhitsListHandle))
      art::fill_ptr_vector(allhits, allhitsListHandle);

    // get Clusters
    if (evt.getByLabel(fClusterModuleLabel, clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);
    if (clusterlist.size() == 0) return;
    // get cluster - hit associations
    art::FindManyP<recob::Hit> fmCluHits(clusterListHandle, evt, fClusterModuleLabel);

    // get Vertices
    if (evt.getByLabel(fVertexModuleLabel, VtxListHandle))
      art::fill_ptr_vector(vtxlist, VtxListHandle);
    art::FindManyP<recob::Cluster, unsigned short> fmVtxCls(VtxListHandle, evt, fVertexModuleLabel);

    std::vector<CluLen> clulens;

    unsigned short itr, tID, tIndex;

    // maximum error (see MakeClusterChains) for considering clusters broken
    fMaxMergeError = 30;
    // Cut on the error for a definitely broken cluster. Clusters with fMergeErrorCut < MergeError < fMaxMergeError
    // are possibly broken clusters but we will consider these when matching between planes
    fMergeErrorCut = 10;

    std::vector<art::Ptr<recob::Hit>> tmpHits;
    std::vector<art::Ptr<recob::Cluster>> tmpCls;
    std::vector<art::Ptr<recob::Vertex>> tmpVtx;

    // vector for PFParticle constructor
    std::vector<size_t> dtrIndices;

    // some vectors for recob::Seed
    double sPos[3], sDir[3];
    double sErr[3] = {0, 0, 0};

    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    // check consistency between clusters and associated hits
    std::vector<art::Ptr<recob::Hit>> clusterhits;
    for (std::size_t icl = 0; icl < clusterlist.size(); ++icl) {
      auto const ipl = clusterlist[icl]->Plane().Plane;
      clusterhits = fmCluHits.at(icl);
      if (clusterhits[0]->WireID().Wire != std::nearbyint(clusterlist[icl]->EndWire())) {
        std::cout << "CCTM Cluster-Hit End wire mis-match " << clusterhits[0]->WireID().Wire
                  << " vs " << std::nearbyint(clusterlist[icl]->EndWire()) << " Bail out! \n";
        return;
      }
      for (unsigned short iht = 0; iht < clusterhits.size(); ++iht) {
        if (clusterhits[iht]->WireID().Plane != ipl) {
          std::cout << "CCTM Cluster-Hit plane mis-match " << ipl << " vs "
                    << clusterhits[iht]->WireID().Plane << " on hit " << iht << " Bail out! \n";
          return;
        } // hit-cluster plane mis-match
      }   // iht
    }     // icl
    // end check consistency

    vtx.clear();
    trk.clear();
    for (auto const& tpcid : geom->Iterate<geo::TPCID>()) {
      unsigned int const nplanes = wireReadoutGeom.Nplanes(tpcid);
      if (nplanes > 3) continue;
      for (std::size_t ipl = 0; ipl < 3; ++ipl) {
        cls[ipl].clear();
        clsChain[ipl].clear();
        trkHits[ipl].clear();
      } // ipl
      // FillWireHitRange also calculates the charge in each plane
      FillWireHitRange(tpcid);
      for (auto const& planeid : wireReadoutGeom.Iterate<geo::PlaneID>(tpcid)) {
        clulens.clear();
        // sort clusters by increasing End wire number
        for (std::size_t icl = 0; icl < clusterlist.size(); ++icl) {
          if (clusterlist[icl]->Plane() != planeid) continue;
          CluLen clulen;
          clulen.index = icl;
          clulen.length = clusterlist[icl]->EndWire();
          clulens.push_back(clulen);
        }
        if (clulens.empty()) continue;
        // sort clusters
        std::sort(clulens.begin(), clulens.end(), lessThan);
        if (clulens.empty()) continue;
        for (unsigned short ii = 0; ii < clulens.size(); ++ii) {
          const unsigned short icl = clulens[ii].index;
          clPar clstr;
          clstr.EvtIndex = icl;
          recob::Cluster const& cluster = *(clusterlist[icl]);
          // Begin info -> end index 1 (DS)
          clstr.Wire[1] = cluster.StartWire();
          clstr.Time[1] = cluster.StartTick();
          clstr.X[1] = (float)detProp.ConvertTicksToX(cluster.StartTick(), planeid);
          clstr.Angle[1] = cluster.StartAngle();
          clstr.Slope[1] = std::tan(cluster.StartAngle());
          clstr.Dir[1] = 0;
          clstr.Charge[1] = ChgNorm[planeid.Plane] * cluster.StartCharge();
          // this will be filled later
          clstr.ChgNear[1] = 0;
          clstr.VtxIndex[1] = -1;
          clstr.mVtxIndex[1] = -1;
          clstr.BrkIndex[1] = -1;
          clstr.MergeError[1] = fMaxMergeError;
          // End info -> end index 0 (US)
          clstr.Wire[0] = cluster.EndWire();
          clstr.Time[0] = cluster.EndTick();
          clstr.X[0] = (float)detProp.ConvertTicksToX(cluster.EndTick(), planeid);
          clstr.Angle[0] = cluster.EndAngle();
          clstr.Slope[0] = std::tan(cluster.EndAngle());
          clstr.Dir[0] = 0;
          if (clstr.Time[1] > clstr.Time[0]) {
            clstr.Dir[0] = 1;
            clstr.Dir[1] = -1;
          }
          else {
            clstr.Dir[0] = -1;
            clstr.Dir[1] = 1;
          }
          clstr.Charge[0] = ChgNorm[planeid.Plane] * cluster.EndCharge();
          // this will be filled later
          clstr.ChgNear[1] = 0;
          clstr.VtxIndex[0] = -1;
          clstr.mVtxIndex[0] = -1;
          clstr.BrkIndex[0] = -1;
          clstr.MergeError[0] = fMaxMergeError;
          // other info
          clstr.InTrack = -1;
          clstr.Length = (unsigned short)(0.5 + clstr.Wire[1] - clstr.Wire[0]);
          clstr.TotChg = ChgNorm[planeid.Plane] * cluster.Integral();
          if (clstr.TotChg <= 0) clstr.TotChg = 1;
          clusterhits = fmCluHits.at(icl);
          if (clusterhits.size() == 0) {
            mf::LogError("CCTM") << "No associated cluster hits " << icl;
            continue;
          }
          // correct charge for missing cluster hits
          clstr.TotChg *= clstr.Length / (float)clusterhits.size();
          cls[planeid.Plane].push_back(clstr);
        } // ii (icl)
      }   // planeid

      FillChgNear(detProp, tpcid);

      // and finally the vertices
      double xyz[3];
      for (unsigned short ivx = 0; ivx < vtxlist.size(); ++ivx) {
        vtxPar aVtx;
        aVtx.ID = ivx + 1;
        aVtx.EvtIndex = ivx;
        vtxlist[ivx]->XYZ(xyz);
        aVtx.X = xyz[0];
        aVtx.Y = xyz[1];
        aVtx.Z = xyz[2];
        aVtx.nClusInPln[0] = 0;
        aVtx.nClusInPln[1] = 0;
        aVtx.nClusInPln[2] = 0;
        std::vector<art::Ptr<recob::Cluster>> const& vtxCls = fmVtxCls.at(ivx);
        std::vector<const unsigned short*> const& vtxClsEnd = fmVtxCls.data(ivx);
        for (unsigned short vcass = 0; vcass < vtxCls.size(); ++vcass) {
          auto const icl = vtxCls[vcass].key();
          // the cluster plane
          auto const ipl = vtxCls[vcass]->Plane().Plane;
          auto end = *vtxClsEnd[vcass];
          if (end > 1)
            throw cet::exception("CCTM")
              << "Invalid end data from vertex - cluster association" << end;
          bool gotit = false;
          // CCTM end is opposite of CC end
          end = 1 - end;
          for (unsigned short jcl = 0; jcl < cls[ipl].size(); ++jcl) {
            if (cls[ipl][jcl].EvtIndex == icl) {
              cls[ipl][jcl].VtxIndex[end] = ivx;
              ++aVtx.nClusInPln[ipl];
              gotit = true;
              break;
            } // index check
          }   // jcl
          if (!gotit)
            throw cet::exception("CCTM") << "Bad index from vertex - cluster association" << icl;
        } // icl
        vtx.push_back(aVtx);
      } // ivx
      // Find broken clusters
      MakeClusterChains(detProp, fmCluHits, tpcid);
      FindMaybeVertices(tpcid);

      // call algorithms in the specified order
      matcomb.clear();
      for (algIndex = 0; algIndex < fMatchAlgs.size(); ++algIndex) {
        if (fMatchAlgs[algIndex] == 1) {
          prt = (fDebugAlg == 1);
          VtxMatch(detProp, fmCluHits, tpcid);
          if (fMakeAlgTracks[algIndex]) SortMatches(detProp, fmCluHits, 1, tpcid);
        }
        else if (fMatchAlgs[algIndex] == 2) {
          prt = (fDebugAlg == 2);
          PlnMatch(detProp, tpcid);
          if (fMakeAlgTracks[algIndex]) SortMatches(detProp, fmCluHits, 2, tpcid);
        }
        if (prt) PrintClusters(detProp, tpcid);
      } // algIndex
      prt = false;
      pfpToTrkID.clear();
      // Determine the vertex/track hierarchy
      if (fMakePFPs) {
        TagCosmics(tpcid);
        MakeFamily(tpcid);
      }
      FitVertices(detProp, tpcid);

      // Make PFParticles -> tracks
      for (unsigned short ipf = 0; ipf < pfpToTrkID.size(); ++ipf) {
        // trackID of this PFP
        tID = pfpToTrkID[ipf];
        if (tID > trk.size() + 1) {
          mf::LogWarning("CCTM") << "Bad track ID";
          continue;
        }
        dtrIndices.clear();
        // load the daughter PFP indices
        mf::LogVerbatim("CCTM") << "PFParticle " << ipf << " tID " << tID;
        for (unsigned short jpf = 0; jpf < pfpToTrkID.size(); ++jpf) {
          itr = pfpToTrkID[jpf] - 1; // convert from track ID to track index
          if (trk[itr].MomID == tID) dtrIndices.push_back(jpf);
          if (trk[itr].MomID == tID)
            mf::LogVerbatim("CCTM") << " dtr jpf " << jpf << " itr " << itr;
        } // jpf
        unsigned short parentIndex = USHRT_MAX;
        if (tID == 0) {
          // make neutrino PFP USHRT_MAX == primary PFP
          recob::PFParticle pfp(14, ipf, parentIndex, dtrIndices);
          pcol->emplace_back(std::move(pfp));
          for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
            if (!vtx[ivx].Neutrino) continue;
            // make the vertex
            xyz[0] = vtx[ivx].X;
            xyz[1] = vtx[ivx].Y;
            xyz[2] = vtx[ivx].Z;
            size_t vStart = vcol->size();
            recob::Vertex vertex(xyz, vtx[ivx].ID);
            vcol->emplace_back(std::move(vertex));
            size_t vEnd = vcol->size();
            // PFParticle - vertex association
            util::CreateAssn(evt, *pcol, *vcol, *pvassn, vStart, vEnd);
            vtx[ivx].ID = -vtx[ivx].ID;
            break;
          } // ivx
        }
        else if (tID > 0) {
          // make non-neutrino PFP. Need to find the parent PFP index
          // Track index of this PFP
          tIndex = tID - 1;
          // trackID of this PFP parent
          for (unsigned short ii = 0; ii < pfpToTrkID.size(); ++ii) {
            if (pfpToTrkID[ii] == trk[tIndex].MomID) {
              parentIndex = ii;
              break;
            }
          } // ii
          // PFParticle
          mf::LogVerbatim("CCTM") << "daughters tID " << tID << " pdg " << trk[tIndex].PDGCode
                                  << " ipf " << ipf << " parentIndex " << parentIndex
                                  << " dtr size " << dtrIndices.size();
          recob::PFParticle pfp(trk[tIndex].PDGCode, ipf, parentIndex, dtrIndices);
          pcol->emplace_back(std::move(pfp));
          // track
          // make the track
          size_t tStart = tcol->size();
          recob::Track track(
            recob::TrackTrajectory(recob::tracking::convertCollToPoint(trk[tIndex].TrjPos),
                                   recob::tracking::convertCollToVector(trk[tIndex].TrjDir),
                                   recob::Track::Flags_t(trk[itr].TrjPos.size()),
                                   false),
            0,
            -1.,
            0,
            recob::tracking::SMatrixSym55(),
            recob::tracking::SMatrixSym55(),
            tID);
          tcol->emplace_back(std::move(track));
          size_t tEnd = tcol->size();
          // PFParticle - track association
          util::CreateAssn(evt, *pcol, *tcol, *ptassn, tStart, tEnd);
          // flag this track as already put in the event
          trk[tIndex].ID = -trk[tIndex].ID;
          // track -> hits association
          tmpHits.clear();
          for (unsigned int ipl = 0; ipl < nplanes; ++ipl)
            tmpHits.insert(
              tmpHits.end(), trk[tIndex].TrkHits[ipl].begin(), trk[tIndex].TrkHits[ipl].end());
          util::CreateAssn(evt, *tcol, tmpHits, *thassn);
          // Find seed hits and the end of the track that is best
          unsigned int end = 0;
          unsigned short itj = 0;
          if (end > 0) itj = trk[tIndex].TrjPos.size() - 1;
          for (unsigned short ii = 0; ii < 3; ++ii) {
            sPos[ii] = trk[tIndex].TrjPos[itj](ii);
            sDir[ii] = trk[tIndex].TrjDir[itj](ii);
          } // ii
          size_t sStart = scol->size();
          recob::Seed seed(sPos, sDir, sErr, sErr);
          scol->emplace_back(std::move(seed));
          size_t sEnd = scol->size();
          // PFP-seed association
          util::CreateAssn(evt, *pcol, *scol, *psassn, sStart, sEnd);
          // Seed-hit association
          tmpHits.clear();
          for (unsigned int ipl = 0; ipl < nplanes; ++ipl)
            tmpHits.insert(tmpHits.end(), seedHits[ipl].begin(), seedHits[ipl].end());
          util::CreateAssn(evt, *scol, tmpHits, *shassn);
          // cluster association
          // PFP-cluster association
          tmpCls.clear();
          for (unsigned short ii = 0; ii < trk[tIndex].ClsEvtIndices.size(); ++ii) {
            auto const icl = trk[tIndex].ClsEvtIndices[ii];
            tmpCls.push_back(clusterlist[icl]);
          } // ii
          util::CreateAssn(evt, *pcol, tmpCls, *pcassn);
        } // non-neutrino PFP
      }   // ipf
      // make non-PFP tracks
      for (unsigned short itr = 0; itr < trk.size(); ++itr) {
        // ignore already saved tracks
        if (trk[itr].ID < 0) continue;
        recob::Track track(
          recob::TrackTrajectory(recob::tracking::convertCollToPoint(trk[itr].TrjPos),
                                 recob::tracking::convertCollToVector(trk[itr].TrjDir),
                                 recob::Track::Flags_t(trk[itr].TrjPos.size()),
                                 false),
          0,
          -1.,
          0,
          recob::tracking::SMatrixSym55(),
          recob::tracking::SMatrixSym55(),
          trk[itr].ID);
        tcol->emplace_back(std::move(track));
        tmpHits.clear();
        for (unsigned int ipl = 0; ipl < nplanes; ++ipl)
          tmpHits.insert(tmpHits.end(), trk[itr].TrkHits[ipl].begin(), trk[itr].TrkHits[ipl].end());
        util::CreateAssn(evt, *tcol, tmpHits, *thassn);
      } // itr
      // make remnant vertices
      for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
        if (vtx[ivx].ID < 0) continue;
        xyz[0] = vtx[ivx].X;
        xyz[1] = vtx[ivx].Y;
        xyz[2] = vtx[ivx].Z;
        recob::Vertex vertex(xyz, vtx[ivx].ID);
        vcol->emplace_back(std::move(vertex));
      }
      if (fDebugAlg > 0) PrintTracks();

      double orphanLen = 0;
      for (unsigned int ipl = 0; ipl < nplanes; ++ipl) {
        for (std::size_t icl = 0; icl < cls[ipl].size(); ++icl) {
          if (cls[ipl][icl].Length > 40 && cls[ipl][icl].InTrack < 0) {
            orphanLen += cls[ipl][icl].Length;
            // unused cluster
            mf::LogVerbatim("CCTM")
              << "Orphan long cluster " << ipl << ":" << icl << ":" << cls[ipl][icl].Wire[0] << ":"
              << (int)cls[ipl][icl].Time[0] << " length " << cls[ipl][icl].Length;
          }
        } // icl

        cls[ipl].clear();
        clsChain[ipl].clear();
        std::cout << "Total orphan length " << orphanLen << "\n";
        trkHits[ipl].clear();
        seedHits[ipl].clear();
        vxCls[ipl].clear();
      } // ipl
    }   // tpcs

    evt.put(std::move(pcol));
    evt.put(std::move(ptassn));
    evt.put(std::move(pcassn));
    evt.put(std::move(pvassn));
    evt.put(std::move(psassn));
    evt.put(std::move(tcol));
    evt.put(std::move(thassn));
    evt.put(std::move(scol));
    evt.put(std::move(vcol));

    // final cleanup
    vtx.clear();
    trk.clear();
    allhits.clear();
    matcomb.clear();
    pfpToTrkID.clear();

  } // produce

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FitVertices(detinfo::DetectorPropertiesData const& detProp,
                                 geo::TPCID const& tpcid)
  {
    std::vector<std::vector<geo::WireID>> hitWID;
    std::vector<std::vector<double>> hitX;
    std::vector<std::vector<double>> hitXErr;
    TVector3 pos, posErr;
    std::vector<TVector3> trkDir;
    std::vector<TVector3> trkDirErr;
    float ChiDOF;

    if (fNVtxTrkHitsFit == 0) return;

    unsigned short indx, indx2, iht, nHitsFit;

    for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      if (!vtx[ivx].Neutrino) continue;
      hitWID.clear();
      hitX.clear();
      hitXErr.clear();
      trkDir.clear();
      // find the tracks associated with this vertex
      unsigned int thePln, theTPC, theCst;
      for (unsigned short itk = 0; itk < trk.size(); ++itk) {
        for (unsigned short end = 0; end < 2; ++end) {
          if (trk[itk].VtxIndex[end] != ivx) continue;
          unsigned short itj = 0;
          if (end == 1) itj = trk[itk].TrjPos.size() - 1;
          // increase the size of the hit vectors by 1 for this track
          indx = hitX.size();
          hitWID.resize(indx + 1);
          hitX.resize(indx + 1);
          hitXErr.resize(indx + 1);
          trkDir.resize(indx + 1);
          trkDir[indx] = trk[itk].TrjDir[itj];
          auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
          for (unsigned short ipl = 0; ipl < nplanes; ++ipl) {
            if (trk[itk].TrkHits[ipl].size() < 2) continue;
            // make slots for the hits on this track in this plane
            nHitsFit = trk[itk].TrkHits[ipl].size();
            if (nHitsFit > fNVtxTrkHitsFit) nHitsFit = fNVtxTrkHitsFit;
            indx2 = hitWID[indx].size();
            hitWID[indx].resize(indx2 + nHitsFit);
            hitX[indx].resize(indx2 + nHitsFit);
            hitXErr[indx].resize(indx2 + nHitsFit);
            for (unsigned short ii = 0; ii < nHitsFit; ++ii) {
              if (end == 0) { iht = ii; }
              else {
                iht = trk[itk].TrkHits[ipl].size() - ii - 1;
              }
              hitWID[indx][indx2 + ii] = trk[itk].TrkHits[ipl][iht]->WireID();
              thePln = trk[itk].TrkHits[ipl][iht]->WireID().Plane;
              theTPC = trk[itk].TrkHits[ipl][iht]->WireID().TPC;
              theCst = trk[itk].TrkHits[ipl][iht]->WireID().Cryostat;
              hitX[indx][indx2 + ii] = detProp.ConvertTicksToX(
                trk[itk].TrkHits[ipl][iht]->PeakTime(), thePln, theTPC, theCst);
              hitXErr[indx][indx2 + ii] = fHitFitErrFac * trk[itk].TrkHits[ipl][iht]->RMS();
            } // ii
          }   // ipl
        }     // end
      }       // itk
      if (hitX.size() < 2) {
        mf::LogVerbatim("CCTM") << "Not enough hits to fit vtx " << ivx;
        continue;
      } // hitX.size() < 2
      pos(0) = vtx[ivx].X;
      pos(1) = vtx[ivx].Y;
      pos(2) = vtx[ivx].Z;
      fVertexFitAlg.VertexFit(hitWID, hitX, hitXErr, pos, posErr, trkDir, trkDirErr, ChiDOF);
      if (ChiDOF > 3000) continue;
      // update the vertex position
      vtx[ivx].X = pos(0);
      vtx[ivx].Y = pos(1);
      vtx[ivx].Z = pos(2);
      // and the track trajectory
      unsigned short fitTrk = 0;
      for (unsigned short itk = 0; itk < trk.size(); ++itk) {
        for (unsigned short end = 0; end < 2; ++end) {
          if (trk[itk].VtxIndex[end] != ivx) continue;
          unsigned short itj = 0;
          if (end == 1) itj = trk[itk].TrjPos.size() - 1;
          trk[itk].TrjDir[itj] = trkDir[fitTrk];
          ++fitTrk;
        } // end
      }   // itk
    }     // ivx
  }       // FitVertices

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillChgNear(detinfo::DetectorPropertiesData const& detProp,
                                 geo::TPCID const& tpcid)
  {
    // fills the CHgNear array

    unsigned short wire, nwires, indx;
    float dir, ctime, cx, chgWinLo, chgWinHi;
    float cnear;

    for (auto const& planeid : wireReadoutGeom.Iterate<geo::PlaneID>(tpcid)) {
      auto const ipl = planeid.Plane;
      for (unsigned short icl = 0; icl < cls[ipl].size(); ++icl) {
        // find the nearby charge at the US and DS ends
        nwires = cls[ipl][icl].Length / 2;
        if (nwires < 2) continue;
        // maximum of 30 wires for long clusters
        if (nwires > 30) nwires = 30;
        for (unsigned short end = 0; end < 2; ++end) {
          cnear = 0;
          // direction for adding/subtracting wire numbers
          dir = 1 - 2 * end;
          for (unsigned short w = 0; w < nwires; ++w) {
            wire = cls[ipl][icl].Wire[end] + dir * w;
            cx = cls[ipl][icl].X[end] + dir * w * cls[ipl][icl].Slope[end] * fWirePitch;
            ctime = detProp.ConvertXToTicks(cx, planeid);
            chgWinLo = ctime - fChgWindow;
            chgWinHi = ctime + fChgWindow;
            indx = wire - firstWire[ipl];
            if (WireHitRange[ipl][indx].first < 0) continue;
            unsigned int firhit = WireHitRange[ipl][indx].first;
            unsigned int lashit = WireHitRange[ipl][indx].second;
            for (unsigned int hit = firhit; hit < lashit; ++hit) {
              if (hit > allhits.size() - 1) {
                mf::LogError("CCTM")
                  << "FillChgNear bad lashit " << lashit << " size " << allhits.size() << "\n";
                continue;
              }
              if (allhits[hit]->PeakTime() < chgWinLo) continue;
              if (allhits[hit]->PeakTime() > chgWinHi) continue;
              cnear += allhits[hit]->Integral();
            } // hit
          }   // w
          cnear /= (float)(nwires - 1);
          if (cnear > cls[ipl][icl].Charge[end]) {
            cls[ipl][icl].ChgNear[end] = ChgNorm[ipl] * cnear / cls[ipl][icl].Charge[end];
          }
          else {
            cls[ipl][icl].ChgNear[end] = 0;
          }
        } // end
      }   // icl
    }     // ipl

  } // FillChgNear

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::MakeFamily(geo::TPCID const& tpcid)
  {
    // define the track/vertex and mother/daughter relationships

    unsigned short ivx, itr, ipl, ii, jtr;
    unsigned short nus, nds, nuhs, ndhs;
    float longUSTrk, longDSTrk, qual;

    // min distance^2 between a neutrino vertex candidate and a through
    // going muon
    float tgMuonCut2 = 9;

    // struct for neutrino vertex candidates
    struct NuVtx {
      unsigned short VtxIndex;
      unsigned short nUSTk;
      unsigned short nDSTk;
      unsigned short nUSHit;
      unsigned short nDSHit;
      float longUSTrk;
      float longDSTrk;
      float Qual;
    };
    std::vector<NuVtx> nuVtxCand;

    NuVtx aNuVtx;

    // analyze all of the vertices
    float best = 999, dx, dy, dz, dr;
    short imbest = -1;
    bool skipVtx;
    unsigned short itj;
    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    for (ivx = 0; ivx < vtx.size(); ++ivx) {
      vtx[ivx].Neutrino = false;
      nus = 0;
      nds = 0;
      nuhs = 0;
      ndhs = 0;
      longUSTrk = 0;
      longDSTrk = 0;
      skipVtx = false;
      // skip vertices that are close to through-going muons
      for (itr = 0; itr < trk.size(); ++itr) {
        if (trk[itr].ID < 0) continue;
        if (trk[itr].PDGCode != 13) continue;
        for (itj = 0; itj < trk[itr].TrjPos.size(); ++itj) {
          dx = trk[itr].TrjPos[itj](0) - vtx[ivx].X;
          dy = trk[itr].TrjPos[itj](1) - vtx[ivx].Y;
          dz = trk[itr].TrjPos[itj](2) - vtx[ivx].Z;
          dr = dx * dx + dy * dy + dz * dz;
          if (dr < tgMuonCut2) {
            skipVtx = true;
            break;
          }
          if (skipVtx) break;
        } // itj
        if (skipVtx) break;
      } // itr
      if (skipVtx) continue;
      for (itr = 0; itr < trk.size(); ++itr) {
        if (trk[itr].ID < 0) continue;
        if (trk[itr].VtxIndex[0] == ivx) {
          // DS-going track
          ++nds;
          if (trk[itr].Length > longDSTrk) longDSTrk = trk[itr].Length;
          for (ipl = 0; ipl < nplanes; ++ipl)
            ndhs += trk[itr].TrkHits[ipl].size();
        } // DS-going track
        // Reject this vertex as a neutrino candidate if the track being
        // considered has a starting vertex
        if (trk[itr].VtxIndex[1] == ivx && trk[itr].VtxIndex[0] >= 0) {
          skipVtx = true;
          break;
        } // trk[itr].VtxIndex[0] > 0
        if (trk[itr].VtxIndex[1] == ivx && trk[itr].VtxIndex[0] < 0) {
          // US-going track w no US vertex
          ++nus;
          if (trk[itr].Length > longUSTrk) longUSTrk = trk[itr].Length;
          for (ipl = 0; ipl < nplanes; ++ipl)
            nuhs += trk[itr].TrkHits[ipl].size();
        } // US-going track
      }   // itr
      if (skipVtx) continue;
      if (nds == 0) continue;
      qual = 1 / (float)nds;
      qual /= (float)ndhs;
      if (nus > 0) qual *= (float)nuhs / (float)ndhs;
      if (qual < best) {
        best = qual;
        imbest = ivx;
      }
      if (nds > 0 && longDSTrk > 5) {
        // at least one longish track going DS
        aNuVtx.VtxIndex = ivx;
        aNuVtx.nUSTk = nus;
        aNuVtx.nDSTk = nds;
        aNuVtx.nUSHit = nuhs;
        aNuVtx.nDSHit = ndhs;
        aNuVtx.longUSTrk = longUSTrk;
        aNuVtx.longDSTrk = longDSTrk;
        aNuVtx.Qual = qual;
        nuVtxCand.push_back(aNuVtx);
      }
    } // ivx
    if (imbest < 0) return;

    // Found the neutrino interaction vertex
    ivx = imbest;
    vtx[ivx].Neutrino = true;
    // Put a 0 into the PFParticle -> track vector to identify a neutrino
    // track with no trajectory or hits
    pfpToTrkID.push_back(0);

    // list of DS-going current generation daughters so we can fix next
    // generation daughter assignments. This code section assumes that there
    // are no decays or 2ndry interactions with US-going primary tracks.
    std::vector<unsigned short> dtrGen;
    std::vector<unsigned short> dtrNextGen;
    for (itr = 0; itr < trk.size(); ++itr) {
      if (trk[itr].ID < 0) continue;
      if (trk[itr].VtxIndex[0] == ivx) {
        // DS-going primary track
        trk[itr].MomID = 0;
        // call every track coming from a neutrino vertex a proton
        trk[itr].PDGCode = 2212;
        pfpToTrkID.push_back(trk[itr].ID);
        dtrGen.push_back(itr);
      } // DS-going primary track
      if (trk[itr].VtxIndex[1] == ivx) {
        // US-going primary track
        trk[itr].MomID = 0;
        // call every track coming from a neutrino vertex a proton
        trk[itr].PDGCode = 2212;
        pfpToTrkID.push_back(trk[itr].ID);
        // reverse the track trajectory
        std::reverse(trk[itr].TrjPos.begin(), trk[itr].TrjPos.end());
        for (ii = 0; ii < trk[itr].TrjDir.size(); ++ii)
          trk[itr].TrjDir[ii] = -trk[itr].TrjDir[ii];
      } // DS-going primary track
    }   // itr

    if (dtrGen.size() == 0) return;

    unsigned short tmp, indx;
    unsigned short nit = 0;

    // follow daughters through all generations (< 10). Daughter tracks
    // may go US or DS
    while (nit < 10) {
      ++nit;
      dtrNextGen.clear();
      // Look for next generation daughters
      for (ii = 0; ii < dtrGen.size(); ++ii) {
        itr = dtrGen[ii];
        if (trk[itr].VtxIndex[1] >= 0) {
          // found a DS vertex
          ivx = trk[itr].VtxIndex[1];
          // look for a track associated with this vertex
          for (jtr = 0; jtr < trk.size(); ++jtr) {
            if (jtr == itr) continue;
            if (trk[jtr].VtxIndex[0] == ivx) {
              // DS-going track
              indx = trk[itr].DtrID.size();
              trk[itr].DtrID.resize(indx + 1);
              trk[itr].DtrID[indx] = jtr;
              trk[jtr].MomID = trk[itr].ID;
              // call all daughters pions
              trk[jtr].PDGCode = 211;
              dtrNextGen.push_back(jtr);
              pfpToTrkID.push_back(trk[jtr].ID);
            } // DS-going track
            if (trk[jtr].VtxIndex[1] == ivx) {
              // US-going track
              indx = trk[itr].DtrID.size();
              trk[itr].DtrID.resize(indx + 1);
              trk[itr].DtrID[indx] = jtr;
              trk[jtr].MomID = trk[itr].ID;
              // call all daughters pions
              trk[jtr].PDGCode = 211;
              dtrNextGen.push_back(jtr);
              pfpToTrkID.push_back(trk[jtr].ID);
              // reverse the track trajectory
              std::reverse(trk[jtr].TrjPos.begin(), trk[jtr].TrjPos.end());
              for (unsigned short jj = 0; jj < trk[jtr].TrjDir.size(); ++jj)
                trk[jtr].TrjDir[jj] = -trk[jtr].TrjDir[jj];
              // interchange the trk - vtx assignments
              tmp = trk[jtr].VtxIndex[0];
              trk[jtr].VtxIndex[0] = trk[jtr].VtxIndex[1];
              trk[jtr].VtxIndex[1] = tmp;
            } // DS-going track
          }   // jtr
        }     // trk[itr].VtxIndex[0] >= 0
      }       // ii (itr)
      // break out if no next gen daughters found
      if (dtrNextGen.size() == 0) break;
      dtrGen = dtrNextGen;
    } // nit

  } // MakeFamily

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::TagCosmics(geo::TPCID const& tpcid)
  {
    // Make cosmic ray PFParticles
    unsigned short ipf, itj;
    bool skipit = true;

    // Y,Z limits of the detector
    const geo::TPCGeo& thetpc = geom->TPC(tpcid);
    auto const world = thetpc.GetCenter();

    float XLo = world.X() - thetpc.HalfWidth() + fFiducialCut;
    float XHi = world.X() + thetpc.HalfWidth() - fFiducialCut;
    float YLo = world.Y() - thetpc.HalfHeight() + fFiducialCut;
    float YHi = world.Y() + thetpc.HalfHeight() - fFiducialCut;
    float ZLo = world.Z() - thetpc.Length() / 2 + fFiducialCut;
    float ZHi = world.Z() + thetpc.Length() / 2 - fFiducialCut;

    bool startsIn, endsIn;

    for (unsigned short itk = 0; itk < trk.size(); ++itk) {
      // ignore already used tracks
      if (trk[itk].ID < 0) continue;
      // ignore short tracks (< 10 cm)
      if (trk[itk].Length < 10) continue;
      // check for already identified PFPs
      skipit = false;
      for (ipf = 0; ipf < pfpToTrkID.size(); ++ipf) {
        if (pfpToTrkID[ipf] == trk[itk].ID) {
          skipit = true;
          break;
        }
      } // ipf
      if (skipit) continue;
      startsIn = true;
      if (trk[itk].TrjPos[0](0) < XLo || trk[itk].TrjPos[0](0) > XHi) startsIn = false;
      if (trk[itk].TrjPos[0](1) < YLo || trk[itk].TrjPos[0](1) > YHi) startsIn = false;
      if (trk[itk].TrjPos[0](2) < ZLo || trk[itk].TrjPos[0](2) > ZHi) startsIn = false;
      if (startsIn) continue;
      endsIn = true;
      itj = trk[itk].TrjPos.size() - 1;
      if (trk[itk].TrjPos[itj](0) < XLo || trk[itk].TrjPos[itj](0) > XHi) endsIn = false;
      if (trk[itk].TrjPos[itj](1) < YLo || trk[itk].TrjPos[itj](1) > YHi) endsIn = false;
      if (trk[itk].TrjPos[itj](2) < ZLo || trk[itk].TrjPos[itj](2) > ZHi) endsIn = false;
      if (endsIn) continue;
      // call it a cosmic muon
      trk[itk].PDGCode = 13;
      pfpToTrkID.push_back(trk[itk].ID);
    } // itk

    if (fDeltaRayCut <= 0) return;

    for (unsigned short itk = 0; itk < trk.size(); ++itk) {
      // find a tagged cosmic ray
      if (trk[itk].PDGCode != 13) continue;

    } // itk

  } // TagCosmics

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::VtxMatch(detinfo::DetectorPropertiesData const& detProp,
                              art::FindManyP<recob::Hit> const& fmCluHits,
                              geo::TPCID const& tpcid)
  {
    // Use vertex assignments to match clusters
    unsigned short ivx, ii, ipl, icl, jj, jpl, jcl, kk, kpl, kcl;
    short idir, iend, jdir, jend, kdir, kend, ioend;

    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    for (ivx = 0; ivx < vtx.size(); ++ivx) {

      // list of cluster chains associated with this vertex in each plane
      for (ipl = 0; ipl < nplanes; ++ipl) {
        vxCls[ipl].clear();
        for (icl = 0; icl < clsChain[ipl].size(); ++icl) {
          if (clsChain[ipl][icl].InTrack >= 0) continue;
          for (iend = 0; iend < 2; ++iend) {
            if (clsChain[ipl][icl].VtxIndex[iend] == vtx[ivx].EvtIndex) vxCls[ipl].push_back(icl);
          } // end
        }   // icl
      }     // ipl

      if (prt) {
        mf::LogVerbatim myprt("CCTM");
        myprt << "VtxMatch: Vertex ID " << vtx[ivx].EvtIndex << "\n";
        for (ipl = 0; ipl < nplanes; ++ipl) {
          myprt << "ipl " << ipl << " cls";
          for (unsigned short ii = 0; ii < vxCls[ipl].size(); ++ii)
            myprt << " " << vxCls[ipl][ii];
          myprt << "\n";
        } // ipl
      }   // prt
      // match between planes, requiring clusters matched to this vertex
      iend = 0;
      jend = 0;
      bool gotkcl;
      float totErr;
      for (ipl = 0; ipl < nplanes; ++ipl) {
        if (nplanes == 2 && ipl > 0) continue;
        for (ii = 0; ii < vxCls[ipl].size(); ++ii) {
          icl = vxCls[ipl][ii];
          // ignore used clusters
          if (clsChain[ipl][icl].InTrack >= 0) continue;
          jpl = (ipl + 1) % nplanes;
          kpl = (jpl + 1) % nplanes;
          for (jj = 0; jj < vxCls[jpl].size(); ++jj) {
            jcl = vxCls[jpl][jj];
            if (clsChain[jpl][jcl].InTrack >= 0) continue;
            for (iend = 0; iend < 2; ++iend) {
              if (clsChain[ipl][icl].VtxIndex[iend] != vtx[ivx].EvtIndex) continue;
              ioend = 1 - iend;
              idir = clsChain[ipl][icl].Dir[iend];
              for (jend = 0; jend < 2; ++jend) {
                if (clsChain[jpl][jcl].VtxIndex[jend] != vtx[ivx].EvtIndex) continue;
                jdir = clsChain[jpl][jcl].Dir[jend];
                if (idir != 0 && jdir != 0 && idir != jdir) continue;
                // ignore outrageously bad other end X matches
                if (fabs(clsChain[jpl][jcl].X[1 - jend] - clsChain[ipl][icl].X[ioend]) > 50)
                  continue;
                MatchPars match;
                match.Cls[ipl] = icl;
                match.End[ipl] = iend;
                match.Cls[jpl] = jcl;
                match.End[jpl] = jend;
                match.Vtx = ivx;
                match.oVtx = -1;
                // set large so that DupMatch doesn't get confused when called before FillEndMatch
                match.Err = 1E6;
                match.oErr = 1E6;
                if (nplanes == 2) {
                  FillEndMatch2(match);
                  if (prt)
                    mf::LogVerbatim("CCTM")
                      << "FillEndMatch2: Err " << match.Err << " oErr " << match.oErr;
                  if (match.Err + match.oErr > 100) continue;
                  if (DupMatch(match, nplanes)) continue;
                  matcomb.push_back(match);
                  continue;
                }
                match.Cls[kpl] = -1;
                match.End[kpl] = 0;
                if (prt)
                  mf::LogVerbatim("CCTM")
                    << "VtxMatch: check " << ipl << ":" << icl << ":" << iend << " and " << jpl
                    << ":" << jcl << ":" << jend << " for cluster in kpl " << kpl;
                gotkcl = false;
                for (kk = 0; kk < vxCls[kpl].size(); ++kk) {
                  kcl = vxCls[kpl][kk];
                  if (clsChain[kpl][kcl].InTrack >= 0) continue;
                  for (kend = 0; kend < 2; ++kend) {
                    kdir = clsChain[kpl][kcl].Dir[kend];
                    if (idir != 0 && kdir != 0 && idir != kdir) continue;
                    if (clsChain[kpl][kcl].VtxIndex[kend] != vtx[ivx].EvtIndex) continue;
                    // rough check of other end match
                    // TODO: SHOWER-LIKE CLUSTER CHECK
                    match.Cls[kpl] = kcl;
                    match.End[kpl] = kend;
                    // first call to ignore redundant matches
                    if (DupMatch(match, nplanes)) continue;
                    FillEndMatch(detProp, tpcid, match);
                    // ignore if no signal at the other end
                    if (match.Chg[kpl] <= 0) continue;
                    if (match.Err + match.oErr > 100) continue;
                    // second call to keep matches with better error
                    if (DupMatch(match, nplanes)) continue;
                    matcomb.push_back(match);
                    gotkcl = true;
                    //                    break;
                  } // kend
                }   // kk -> kcl
                if (gotkcl) continue;
                // look for a cluster that missed the vertex assignment
                float best = 10;
                short kbst = -1;
                unsigned short kbend = 0;
                if (prt)
                  mf::LogVerbatim("CCTM") << "VtxMatch: look for missed cluster chain in kpl";
                for (kcl = 0; kcl < clsChain[kpl].size(); ++kcl) {
                  if (clsChain[kpl][kcl].InTrack >= 0) continue;
                  for (kend = 0; kend < 2; ++kend) {
                    kdir = clsChain[kpl][kcl].Dir[kend];
                    if (idir != 0 && kdir != 0 && idir != kdir) continue;
                    if (clsChain[kpl][kcl].VtxIndex[kend] >= 0) continue;
                    // make a rough dX cut at the match end
                    if (fabs(clsChain[kpl][kcl].X[kend] - vtx[ivx].X) > 5) continue;
                    // and at the other end
                    if (fabs(clsChain[kpl][kcl].X[1 - kend] - clsChain[ipl][icl].X[ioend]) > 50)
                      continue;
                    // check the error
                    match.Cls[kpl] = kcl;
                    match.End[kpl] = kend;
                    if (DupMatch(match, nplanes)) continue;
                    FillEndMatch(detProp, tpcid, match);
                    totErr = match.Err + match.oErr;
                    if (prt) {
                      mf::LogVerbatim myprt("CCTM");
                      myprt << "VtxMatch: Chk missing cluster match ";
                      for (unsigned short ii = 0; ii < nplanes; ++ii)
                        myprt << " " << ii << ":" << match.Cls[ii] << ":" << match.End[ii];
                      myprt << " Err " << match.Err << "\n";
                    }
                    if (totErr > 100) continue;
                    if (totErr < best) {
                      best = totErr;
                      kbst = kcl;
                      kbend = kend;
                    }
                  } // kend
                }   // kcl
                if (kbst >= 0) {
                  // found a decent match
                  match.Cls[kpl] = kbst;
                  match.End[kpl] = kbend;
                  FillEndMatch(detProp, tpcid, match);
                  matcomb.push_back(match);
                  // assign the vertex to this cluster
                  clsChain[kpl][kbst].VtxIndex[kbend] = ivx;
                  // and update vxCls
                  vxCls[kpl].push_back(kbst);
                }
                else {
                  // Try a 2 plane match if a 3 plane match didn't work
                  match.Cls[kpl] = -1;
                  match.End[kpl] = 0;
                  if (DupMatch(match, nplanes)) continue;
                  FillEndMatch(detProp, tpcid, match);
                  if (match.Err + match.oErr < 100) matcomb.push_back(match);
                }
              } // jend
            }   // iend
          }     // jj
        }       // ii -> icl
      }         // ipl

      if (matcomb.size() == 0) continue;
      SortMatches(detProp, fmCluHits, 1, tpcid);

    } // ivx

    for (ipl = 0; ipl < 3; ++ipl)
      vxCls[ipl].clear();

  } // VtxMatch

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FindMaybeVertices(geo::TPCID const& tpcid)
  {
    // Project clusters to vertices and fill mVtxIndex. No requirement is
    // made that charge exists on the line between the Begin (End) of the
    // cluster and the vertex

    if (vtx.empty()) return;

    for (auto const& plane : wireReadoutGeom.Iterate<geo::PlaneGeo>(tpcid)) {
      auto const& planeID = plane.ID();
      unsigned int const ipl = planeID.Cryostat;
      for (unsigned int icl = 0; icl < cls[ipl].size(); ++icl) {
        for (unsigned int end = 0; end < 2u; ++end) {
          // ignore already attached clusters
          if (cls[ipl][icl].VtxIndex[end] >= 0) continue;
          short ibstvx = -1;
          float best = 1.;
          // index of the other end
          unsigned int oend = 1 - end;
          for (std::size_t ivx = 0; ivx < vtx.size(); ++ivx) {
            // ignore if the other end is attached to this vertex (which can happen with short clusters)
            if (cls[ipl][icl].VtxIndex[oend] == static_cast<short>(ivx)) continue;
            float const dWire = plane.WireCoordinate(geo::Point_t{0, vtx[ivx].Y, vtx[ivx].Z}) -
                                cls[ipl][icl].Wire[end];
            if (end == 0) {
              if (dWire > 30 || dWire < -2) continue;
            }
            else {
              if (dWire < -30 || dWire > 2) continue;
            }
            // project the cluster to the vertex wire
            float const dX = fabs(cls[ipl][icl].X[end] +
                                  cls[ipl][icl].Slope[end] * fWirePitch * dWire - vtx[ivx].X);
            if (dX < best) {
              best = dX;
              ibstvx = ivx;
            }
          } // ivx
          if (ibstvx >= 0) {
            // attach
            cls[ipl][icl].VtxIndex[end] = ibstvx;
            cls[ipl][icl].mVtxIndex[end] = ibstvx;
          }
        } // end
      }   // icl
    }     // ipl

  } // FindMaybeVertices

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::MakeClusterChains(detinfo::DetectorPropertiesData const& detProp,
                                       art::FindManyP<recob::Hit> const& fmCluHits,
                                       geo::TPCID const& tpcid)
  {
    unsigned short icl, icl1, icl2;
    float dw, dx, dWCut, dw1Max, dw2Max;
    float dA, dA2, dACut = fMaxDAng, chgAsymCut;
    float dXCut, chgasym, mrgErr;
    // long straight clusters
    bool ls1, ls2;
    bool gotprt = false;
    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    for (unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      if (cls[ipl].size() > 1) {
        for (icl1 = 0; icl1 < cls[ipl].size() - 1; ++icl1) {
          prt = (fDebugAlg == 666 && ipl == fDebugPlane && icl1 == fDebugCluster);
          if (prt) gotprt = true;
          // maximum delta Wire overlap is 10% of the total length
          dw1Max = 0.6 * cls[ipl][icl1].Length;
          ls1 = (cls[ipl][icl1].Length > 100 &&
                 fabs(cls[ipl][icl1].Angle[0] - cls[ipl][icl1].Angle[1]) < 0.04);
          for (icl2 = icl1 + 1; icl2 < cls[ipl].size(); ++icl2) {
            ls2 = (cls[ipl][icl2].Length > 100 &&
                   fabs(cls[ipl][icl2].Angle[0] - cls[ipl][icl2].Angle[1]) < 0.04);
            dw2Max = 0.6 * cls[ipl][icl2].Length;
            // set overlap cut to be the shorter of the two
            dWCut = dw1Max;
            if (dw2Max < dWCut) dWCut = dw2Max;
            // but not exceeding 20 for very long clusters
            if (dWCut > 100) dWCut = 100;
            if (dWCut < 2) dWCut = 2;
            chgAsymCut = fMergeChgAsym;
            // Compare end 1 of icl1 with end 0 of icl2

            if (prt)
              mf::LogVerbatim("CCTM")
                << "MCC P:C:W icl1 " << ipl << ":" << icl1 << ":" << cls[ipl][icl1].Wire[1]
                << " vtx " << cls[ipl][icl1].VtxIndex[1] << " ls1 " << ls1 << " icl2 " << ipl << ":"
                << icl2 << ":" << cls[ipl][icl2].Wire[0] << " vtx " << cls[ipl][icl2].VtxIndex[0]
                << " ls2 " << ls2 << " dWCut " << dWCut;
            if (std::abs(cls[ipl][icl1].Wire[1] - cls[ipl][icl2].Wire[0]) > dWCut) continue;
            // ignore if the clusters begin/end on the same wire
            //            if(cls[ipl][icl1].Wire[1]  == cls[ipl][icl2].Wire[0]) continue;
            // or if the angle exceeds the cut
            float af = AngleFactor(cls[ipl][icl1].Slope[1]);
            dACut = fMaxDAng * af;
            dXCut = fChainMaxdX * 5 * af;
            dA = fabs(cls[ipl][icl1].Angle[1] - cls[ipl][icl2].Angle[0]);
            // compare the match angle at the opposite ends of the clusters.
            // May have a bad end/begin angle if there is a delta-ray in the middle
            dA2 = fabs(cls[ipl][icl1].Angle[0] - cls[ipl][icl2].Angle[1]);

            if (prt)
              mf::LogVerbatim("CCTM")
                << " dA " << dA << " dA2 " << dA2 << " DACut " << dACut << " dXCut " << dXCut;

            if (dA2 < dA) dA = dA2;
            // ignore vertices that have only two associated clusters and
            // the angle is small
            if (dA < fChainVtxAng && cls[ipl][icl1].VtxIndex[1] >= 0) {
              dw = fWirePitch * (cls[ipl][icl2].Wire[0] - cls[ipl][icl1].Wire[1]);
              dx = cls[ipl][icl1].X[1] + cls[ipl][icl1].Slope[1] * dw * fWirePitch -
                   cls[ipl][icl2].X[0];
              unsigned short ivx = cls[ipl][icl1].VtxIndex[1];
              if (vtx[ivx].nClusInPln[ipl] == 2 && fabs(dx) < 1) {
                cls[ipl][icl1].VtxIndex[1] = -2;
                cls[ipl][icl2].VtxIndex[0] = -2;
                vtx[ivx].nClusInPln[ipl] = 0;
                if (prt) mf::LogVerbatim("CCTM") << " clobbered vertex " << ivx;
              } // vertices match
            }   // dA < 0.1 && ...

            // don't merge if a vertex exists at these ends
            if (cls[ipl][icl1].VtxIndex[1] >= 0) continue;
            if (cls[ipl][icl2].VtxIndex[0] >= 0) continue;

            // expand the angle match cut for clusters that appear to be stopping
            if (cls[ipl][icl2].Wire[0] - cls[ipl][icl1].Wire[1] < 3 &&
                (cls[ipl][icl1].Length < 3 || cls[ipl][icl2].Length < 3)) {
              if (prt) mf::LogVerbatim("CCTM") << "Stopping cluster";
              dACut *= 1.5;
              chgAsymCut *= 1.5;
              dXCut *= 3;
            } // stopping US cluster

            // find the angle made by the endpoints of the clusters
            dw = fWirePitch * (cls[ipl][icl2].Wire[0] - cls[ipl][icl1].Wire[1]);
            if (dw != 0) {
              dx = cls[ipl][icl2].X[0] - cls[ipl][icl1].X[1];
              float dA3 = std::abs(atan(dx / dw) - cls[ipl][icl1].Angle[1]);
              if (prt) mf::LogVerbatim("CCTM") << " dA3 " << dA3;
              if (dA3 > dA) dA = dA3;
            }

            // angle matching
            if (dA > dACut) continue;

            if (prt)
              mf::LogVerbatim("CCTM")
                << " rough dX " << fabs(cls[ipl][icl1].X[1] - cls[ipl][icl2].X[0]) << " cut = 20";

            // make a rough dX cut
            if (fabs(cls[ipl][icl1].X[1] - cls[ipl][icl2].X[0]) > 20) continue;

            // handle cosmic ray clusters that are broken at delta rays
            if (ls1 || ls2) {
              // tighter angle cuts but no charge cut
              if (dA > fChainVtxAng) continue;
            }
            else {
              chgasym = fabs(cls[ipl][icl1].Charge[1] - cls[ipl][icl2].Charge[0]);
              chgasym /= cls[ipl][icl1].Charge[1] + cls[ipl][icl2].Charge[0];
              if (prt) mf::LogVerbatim("CCTM") << " chgasym " << chgasym << " cut " << chgAsymCut;
              if (chgasym > chgAsymCut) continue;
            } // ls1 || ls2
            // project the longer cluster to the end of the shorter one
            if (cls[ipl][icl1].Length > cls[ipl][icl2].Length) {
              dw = fWirePitch * (cls[ipl][icl2].Wire[0] - cls[ipl][icl1].Wire[1]);
              dx = cls[ipl][icl1].X[1] + cls[ipl][icl1].Slope[1] * dw * fWirePitch -
                   cls[ipl][icl2].X[0];
            }
            else {
              dw = fWirePitch * (cls[ipl][icl1].Wire[1] - cls[ipl][icl2].Wire[0]);
              dx = cls[ipl][icl2].X[0] + cls[ipl][icl2].Slope[0] * dw * fWirePitch -
                   cls[ipl][icl1].X[1];
            }

            // handle overlapping clusters
            if (dA2 < 0.01 && abs(dx) > dXCut && dx < -1) {
              dx = dXClTraj(fmCluHits, ipl, icl1, 1);
              if (prt) mf::LogVerbatim("CCTM") << " new dx from dXClTraj " << dx;
            }

            if (prt)
              mf::LogVerbatim("CCTM")
                << " X0 " << cls[ipl][icl1].X[1] << " slp " << cls[ipl][icl1].Slope[1] << " dw "
                << dw << " oX " << cls[ipl][icl2].X[0] << " dx " << dx << " cut " << dXCut;

            if (fabs(dx) > dXCut) continue;

            // calculate a merge error that will be used to adjudicate between multiple merge attempts AngleFactor
            float xerr = dx / dXCut;
            float aerr = dA / dACut;
            mrgErr = xerr * xerr + aerr * aerr;

            if (prt)
              mf::LogVerbatim("CCTM")
                << "icl1 mrgErr " << mrgErr << " MergeError " << cls[ipl][icl1].MergeError[1]
                << " icl2 MergeError " << cls[ipl][icl2].MergeError[0];

            // this merge better than a previous one?
            if (mrgErr > cls[ipl][icl1].MergeError[1]) continue;
            if (mrgErr > cls[ipl][icl2].MergeError[0]) continue;

            // un-merge icl1 - this should always be true but check anyway
            if (cls[ipl][icl1].BrkIndex[1] >= 0) {
              unsigned short ocl = cls[ipl][icl1].BrkIndex[1];
              if (prt) mf::LogVerbatim("CCTM") << "clobber old icl1 BrkIndex " << ocl;
              if (cls[ipl][ocl].BrkIndex[0] == icl1) {
                cls[ipl][ocl].BrkIndex[0] = -1;
                cls[ipl][ocl].MergeError[0] = fMaxMergeError;
              }
              if (cls[ipl][ocl].BrkIndex[1] == icl1) {
                cls[ipl][ocl].BrkIndex[1] = -1;
                cls[ipl][ocl].MergeError[1] = fMaxMergeError;
              }
            } // cls[ipl][icl1].BrkIndex[1] >= 0
            cls[ipl][icl1].BrkIndex[1] = icl2;
            cls[ipl][icl1].MergeError[1] = mrgErr;

            // un-merge icl2
            if (cls[ipl][icl2].BrkIndex[0] >= 0) {
              unsigned short ocl = cls[ipl][icl2].BrkIndex[0];
              if (prt) mf::LogVerbatim("CCTM") << "clobber old icl2 BrkIndex " << ocl;
              if (cls[ipl][ocl].BrkIndex[0] == icl1) {
                cls[ipl][ocl].BrkIndex[0] = -1;
                cls[ipl][ocl].MergeError[0] = fMaxMergeError;
              }
              if (cls[ipl][ocl].BrkIndex[1] == icl1) {
                cls[ipl][ocl].BrkIndex[1] = -1;
                cls[ipl][ocl].MergeError[1] = fMaxMergeError;
              }
            } // cls[ipl][icl2].BrkIndex[0] >= 0
            cls[ipl][icl2].BrkIndex[0] = icl1;
            cls[ipl][icl2].MergeError[0] = mrgErr;
            if (prt) mf::LogVerbatim("CCTM") << " merge " << icl1 << " and " << icl2;

          } // icl2
        }   // icl1

        // look for broken clusters in which have a C shape similar to a Begin-Begin vertex. The clusters
        // will have large and opposite sign angles
        bool gotone;
        for (icl1 = 0; icl1 < cls[ipl].size() - 1; ++icl1) {
          gotone = false;
          for (icl2 = icl1 + 1; icl2 < cls[ipl].size(); ++icl2) {
            // check both ends
            for (unsigned short end = 0; end < 2; ++end) {
              // Ignore already identified broken clusters
              if (cls[ipl][icl1].BrkIndex[end] >= 0) continue;
              if (cls[ipl][icl2].BrkIndex[end] >= 0) continue;
              // require a large angle cluster
              if (fabs(cls[ipl][icl1].Angle[end]) < 1) continue;
              // and a second large angle cluster
              if (fabs(cls[ipl][icl2].Angle[end]) < 1) continue;
              if (prt)
                mf::LogVerbatim("CCTM")
                  << "BrokenC: clusters " << cls[ipl][icl1].Wire[end] << ":"
                  << (int)cls[ipl][icl1].Time[end] << " " << cls[ipl][icl2].Wire[end] << ":"
                  << (int)cls[ipl][icl2].Time[end] << " angles " << cls[ipl][icl1].Angle[end] << " "
                  << cls[ipl][icl2].Angle[end] << " dWire "
                  << fabs(cls[ipl][icl1].Wire[end] - cls[ipl][icl2].Wire[end]);
              if (fabs(cls[ipl][icl1].Wire[end] - cls[ipl][icl2].Wire[end]) > 5) continue;
              // This is really crude but maybe OK
              // project 1 -> 2
              float dsl = cls[ipl][icl2].Slope[end] - cls[ipl][icl1].Slope[end];
              float fvw =
                (cls[ipl][icl1].X[end] - cls[ipl][icl1].Wire[end] * cls[ipl][icl1].Slope[end] -
                 cls[ipl][icl2].X[end] + cls[ipl][icl2].Wire[end] * cls[ipl][icl2].Slope[end]) /
                dsl;
              if (prt) mf::LogVerbatim("CCTM") << " fvw " << fvw;
              if (fabs(cls[ipl][icl1].Wire[end] - fvw) > 4) continue;
              if (fabs(cls[ipl][icl2].Wire[end] - fvw) > 4) continue;
              cls[ipl][icl1].BrkIndex[end] = icl2;
              // TODO This could use some improvement if necessary
              cls[ipl][icl1].MergeError[end] = 1;
              cls[ipl][icl2].BrkIndex[end] = icl1;
              cls[ipl][icl2].MergeError[end] = 1;
              gotone = true;
              dx = fabs(cls[ipl][icl1].X[end] - cls[ipl][icl2].X[end]);
              if (prt)
                mf::LogVerbatim("CCTM")
                  << "BrokenC: icl1:W " << icl1 << ":" << cls[ipl][icl1].Wire[end] << " icl2:W "
                  << icl2 << ":" << cls[ipl][icl2].Wire[end] << " end " << end << " dx " << dx;
            } // end
            if (gotone) break;
          } // icl2
        }   // icl1

      } // cls[ipl].size() > 1

      // follow mother-daughter broken clusters and put them in the cluster chain array
      unsigned short mom, momBrkEnd, dtrBrkEnd, nit;
      short dtr;

      std::vector<bool> gotcl(cls[ipl].size());
      for (icl = 0; icl < cls[ipl].size(); ++icl)
        gotcl[icl] = false;
      if (prt)
        mf::LogVerbatim("CCTM") << "ipl " << ipl << " cls.size() " << cls[ipl].size() << "\n";

      std::vector<unsigned short> sCluster;
      std::vector<unsigned short> sOrder;
      for (icl = 0; icl < cls[ipl].size(); ++icl) {
        sCluster.clear();
        sOrder.clear();
        if (gotcl[icl]) continue;
        // don't start with a cluster broken at both ends
        if (cls[ipl][icl].BrkIndex[0] >= 0 && cls[ipl][icl].BrkIndex[1] >= 0) continue;
        for (unsigned int end = 0; end < 2; ++end) {
          if (cls[ipl][icl].BrkIndex[end] < 0) continue;
          if (cls[ipl][icl].MergeError[end] > fMergeErrorCut) continue;
          gotcl[icl] = true;
          mom = icl;
          // end where the mom is broken
          momBrkEnd = end;
          sCluster.push_back(mom);
          if (momBrkEnd == 1) {
            // typical case - broken at the DS end
            sOrder.push_back(0);
          }
          else {
            // broken at the US end
            sOrder.push_back(1);
          }
          dtr = cls[ipl][icl].BrkIndex[end];
          nit = 0;
          while (dtr >= 0 && dtr < (short)cls[ipl].size() && nit < cls[ipl].size()) {
            // determine which end of the dtr should be attached to mom
            for (dtrBrkEnd = 0; dtrBrkEnd < 2; ++dtrBrkEnd)
              if (cls[ipl][dtr].BrkIndex[dtrBrkEnd] == mom) break;
            if (dtrBrkEnd == 2) {
              gotcl[icl] = false;
              break;
            }
            // check for reasonable merge error
            if (cls[ipl][dtr].MergeError[dtrBrkEnd] < fMergeErrorCut) {
              sCluster.push_back(dtr);
              sOrder.push_back(dtrBrkEnd);
              gotcl[dtr] = true;
            }
            ++nit;
            // set up to check the new mom
            mom = dtr;
            // at the other end
            momBrkEnd = 1 - dtrBrkEnd;
            // with the new dtr
            dtr = cls[ipl][mom].BrkIndex[momBrkEnd];
          } // dtr >= 0 ...
          if (dtrBrkEnd == 2) continue;
        } // end

        if (!gotcl[icl]) {
          // a single unbroken cluster
          sCluster.push_back(icl);
          sOrder.push_back(0);
        }

        if (sCluster.size() == 0) {
          mf::LogError("CCTM") << "MakeClusterChains error in plane " << ipl << " cluster " << icl;
          return;
        }

        ClsChainPar ccp;
        // fill the struct parameters assuming this is a cluster chain containing only one cluster
        unsigned short jcl = sCluster[0];
        if (jcl > cls[ipl].size()) std::cout << "oops MCC\n";
        unsigned short oend;
        for (unsigned int end = 0; end < 2; ++end) {
          oend = end;
          if (sOrder[0] > 0) oend = 1 - end;
          ccp.Wire[end] = cls[ipl][jcl].Wire[oend];
          ccp.Time[end] = cls[ipl][jcl].Time[oend];
          ccp.X[end] = cls[ipl][jcl].X[oend];
          ccp.Slope[end] = cls[ipl][jcl].Slope[oend];
          ccp.Angle[end] = cls[ipl][jcl].Angle[oend];
          ccp.Dir[end] = cls[ipl][icl].Dir[oend];
          ccp.VtxIndex[end] = cls[ipl][jcl].VtxIndex[oend];
          ccp.ChgNear[end] = cls[ipl][jcl].ChgNear[oend];
          ccp.mBrkIndex[end] = cls[ipl][jcl].BrkIndex[oend];
        } // end
        ccp.Length = cls[ipl][icl].Length;
        ccp.TotChg = cls[ipl][icl].TotChg;
        ccp.InTrack = -1;

        for (unsigned short ii = 1; ii < sCluster.size(); ++ii) {
          jcl = sCluster[ii];
          if (jcl > cls[ipl].size()) std::cout << "oops MCC\n";
          // end is the end where the break is being mended
          unsigned int end = sOrder[ii];
          if (end > 1) std::cout << "oops2 MCC\n";
          oend = 1 - end;
          // update the parameters at the other end of the chain
          ccp.Wire[1] = cls[ipl][jcl].Wire[oend];
          ccp.Time[1] = cls[ipl][jcl].Time[oend];
          ccp.X[1] = cls[ipl][jcl].X[oend];
          ccp.Slope[1] = cls[ipl][jcl].Slope[oend];
          ccp.Angle[1] = cls[ipl][jcl].Angle[oend];
          ccp.Dir[1] = cls[ipl][jcl].Dir[oend];
          ccp.VtxIndex[1] = cls[ipl][jcl].VtxIndex[oend];
          ccp.ChgNear[1] = cls[ipl][jcl].ChgNear[oend];
          ccp.mBrkIndex[1] = cls[ipl][jcl].BrkIndex[oend];
          ccp.Length += cls[ipl][jcl].Length;
          ccp.TotChg += cls[ipl][jcl].TotChg;
        } // ii
        ccp.ClsIndex = sCluster;
        ccp.Order = sOrder;
        // redo the direction
        if (ccp.Time[1] > ccp.Time[0]) {
          ccp.Dir[0] = 1;
          ccp.Dir[1] = -1;
        }
        else {
          ccp.Dir[0] = -1;
          ccp.Dir[1] = 1;
        }
        clsChain[ipl].push_back(ccp);

      } // icl

      // re-index mBrkIndex to point to a cluster chain, not a cluster
      unsigned short brkCls;
      bool gotit;
      for (unsigned short ccl = 0; ccl < clsChain[ipl].size(); ++ccl) {
        for (unsigned short end = 0; end < 2; ++end) {
          if (clsChain[ipl][ccl].mBrkIndex[end] < 0) continue;
          brkCls = clsChain[ipl][ccl].mBrkIndex[end];
          gotit = false;
          // find this cluster index in a cluster chain
          for (unsigned short ccl2 = 0; ccl2 < clsChain[ipl].size(); ++ccl2) {
            if (ccl2 == ccl) continue;
            if (std::find(clsChain[ipl][ccl2].ClsIndex.begin(),
                          clsChain[ipl][ccl2].ClsIndex.end(),
                          brkCls) == clsChain[ipl][ccl2].ClsIndex.end())
              continue;
            // found it
            clsChain[ipl][ccl].mBrkIndex[end] = ccl2;
            gotit = true;
            break;
          } // ccl2
          if (!gotit)
            mf::LogError("CCTM") << "MCC: Cluster chain " << ccl << " end " << end
                                 << " Failed to find brkCls " << brkCls << " in plane " << ipl;
        } // end
      }   // ccl

    } // ipl

    if (gotprt) PrintClusters(detProp, tpcid);
    prt = false;

  } // MakeClusterChains

  ///////////////////////////////////////////////////////////////////////
  float CCTrackMaker::dXClTraj(art::FindManyP<recob::Hit> const& fmCluHits,
                               unsigned short ipl,
                               unsigned short icl1,
                               unsigned short end1)
  {
    // project cluster icl1 at end1 to find the best intersection with icl2
    float dw, dx, best = 999;
    std::vector<art::Ptr<recob::Hit>> clusterhits = fmCluHits.at(cls[ipl][icl1].EvtIndex);
    for (unsigned short hit = 0; hit < clusterhits.size(); ++hit) {
      dw = clusterhits[hit]->WireID().Wire - cls[ipl][icl1].Wire[end1];
      dx = fabs(cls[ipl][icl1].Time[end1] + dw * fWirePitch * cls[ipl][icl1].Slope[end1] -
                clusterhits[hit]->PeakTime());
      if (dx < best) best = dx;
      if (dx < 0.01) break;
    } // hit
    return best;
  } // dXClTraj

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::StoreTrack(detinfo::DetectorPropertiesData const& detProp,
                                art::FindManyP<recob::Hit> const& fmCluHits,
                                unsigned short imat,
                                unsigned short procCode,
                                geo::TPCID const& tpcid)
  {
    // store the current "under construction" track in the trk vector

    TrkPar newtrk;

    if (imat > matcomb.size() - 1) {
      mf::LogError("CCTM") << "Bad imat in StoreTrack";
      return;
    }

    // ensure there are at least 2 hits in at least 2 planes
    unsigned short nhitinpl = 0;
    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    for (unsigned short ipl = 0; ipl < nplanes; ++ipl)
      if (trkHits[ipl].size() > 1) ++nhitinpl;
    if (nhitinpl < 2) {
      mf::LogError("CCTM") << "StoreTrack: Not enough hits in each plane\n";
      return;
    }
    if (prt)
      mf::LogVerbatim("CCTM") << "In StoreTrack: matcomb " << imat << " cluster chains "
                              << matcomb[imat].Cls[0] << " " << matcomb[imat].Cls[1] << " "
                              << matcomb[imat].Cls[2];

    // Track hit vectors for fitting the trajectory
    std::array<std::vector<geo::WireID>, 3> trkWID;
    std::array<std::vector<double>, 3> trkX;
    std::array<std::vector<double>, 3> trkXErr;

    // track trajectory for a track
    std::vector<TVector3> trkPos;
    std::vector<TVector3> trkDir;

    newtrk.ID = trk.size() + 1;
    newtrk.Proc = procCode;
    newtrk.TrkHits = trkHits;
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    newtrk.VtxIndex = {{-1, -1}};
    newtrk.ChgOrder = 0;
    newtrk.MomID = -1;
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    newtrk.EndInTPC = {{false, false}};
    newtrk.GoodEnd = {{false, false}};
    newtrk.DtrID = {0};
    newtrk.PDGCode = -1;

    unsigned short icl, iht;

    if (prt)
      mf::LogVerbatim("CCTM") << "CCTM: Make traj for track " << newtrk.ID << " procCode "
                              << procCode << " nhits in planes " << trkHits[0].size() << " "
                              << trkHits[1].size() << " " << trkHits[2].size();
    // make the track trajectory
    if (nplanes == 2) {
      trkWID[2].resize(0);
      trkX[2].resize(0);
      trkXErr[2].resize(0);
    }
    for (auto const& planeid : wireReadoutGeom.Iterate<geo::PlaneID>(tpcid)) {
      auto const ipl = planeid.Plane;
      trkWID[ipl].resize(trkHits[ipl].size());
      trkX[ipl].resize(trkHits[ipl].size());
      trkXErr[ipl].resize(trkHits[ipl].size());
      for (iht = 0; iht < trkHits[ipl].size(); ++iht) {
        trkWID[ipl][iht] = trkHits[ipl][iht]->WireID();
        trkX[ipl][iht] = detProp.ConvertTicksToX(trkHits[ipl][iht]->PeakTime(), planeid);
        trkXErr[ipl][iht] =
          fHitFitErrFac * trkHits[ipl][iht]->RMS() * trkHits[ipl][iht]->Multiplicity();
      } // iht
    }   // ipl
    fTrackTrajectoryAlg.TrackTrajectory(trkWID, trkX, trkXErr, trkPos, trkDir);
    if (trkPos.size() < 2) {
      mf::LogError("CCTM") << "StoreTrack: No trajectory points on failed track " << newtrk.ID
                           << " in StoreTrack: matcomb " << imat << " cluster chains "
                           << matcomb[imat].Cls[0] << " " << matcomb[imat].Cls[1] << " "
                           << matcomb[imat].Cls[2];
      // make a garbage trajectory
      trkPos.resize(2);
      trkPos[1](2) = 1;
      trkDir.resize(2);
      trkDir[1](2) = 1;
    }
    newtrk.TrjPos = trkPos;
    newtrk.TrjDir = trkDir;

    if (prt) mf::LogVerbatim("CCTM") << " number of traj points " << trkPos.size();

    // determine if each end is good in the sense that there are hits in each plane
    // that are consistent in time and are presumed to form a good 3D space point
    unsigned short nClose, indx, jndx;
    float xErr;
    for (unsigned int end = 0; end < 2; ++end) {
      nClose = 0;
      for (unsigned short ipl = 0; ipl < nplanes - 1; ++ipl) {
        if (trkX[ipl].size() == 0) continue;
        for (unsigned short jpl = ipl + 1; jpl < nplanes; ++jpl) {
          if (trkX[jpl].size() == 0) continue;
          if (end == 0) {
            indx = 0;
            jndx = 0;
          }
          else {
            indx = trkXErr[ipl].size() - 1;
            jndx = trkXErr[jpl].size() - 1;
          }
          xErr = 3 * (trkXErr[ipl][indx] + trkXErr[jpl][jndx]);
          if (std::abs(trkX[ipl][indx] - trkX[jpl][jndx]) < xErr) ++nClose;
        } // jpl
      }   // ipl
      if (nClose == nplanes) newtrk.GoodEnd[end] = true;
    } // end

    // set trajectory end points to a vertex if one exists
    unsigned short ivx, itj, ccl;
    float dx, dy, dz, dr0, dr1;
    unsigned short attachEnd;
    for (unsigned int end = 0; end < 2; ++end) {
      ivx = USHRT_MAX;
      if (end == 0 && matcomb[imat].Vtx >= 0) ivx = matcomb[imat].Vtx;
      if (end == 1 && matcomb[imat].oVtx >= 0) ivx = matcomb[imat].oVtx;
      if (ivx == USHRT_MAX) continue;
      // determine the proper end using the TrjPos order and brute force
      itj = 0;
      dx = vtx[ivx].X - newtrk.TrjPos[itj](0);
      dy = vtx[ivx].Y - newtrk.TrjPos[itj](1);
      dz = vtx[ivx].Z - newtrk.TrjPos[itj](2);
      dr0 = dx * dx + dy * dy + dz * dz;
      itj = newtrk.TrjPos.size() - 1;
      dx = vtx[ivx].X - newtrk.TrjPos[itj](0);
      dy = vtx[ivx].Y - newtrk.TrjPos[itj](1);
      dz = vtx[ivx].Z - newtrk.TrjPos[itj](2);
      dr1 = dx * dx + dy * dy + dz * dz;
      attachEnd = 1;
      if (dr0 < dr1) {
        itj = 0;
        attachEnd = 0;
        // a really bad match to the vertex
        if (dr0 > 5) return;
      }
      else {
        // a really bad match to the vertex
        if (dr1 > 5) return;
      }
      newtrk.TrjPos[itj](0) = vtx[ivx].X;
      newtrk.TrjPos[itj](1) = vtx[ivx].Y;
      newtrk.TrjPos[itj](2) = vtx[ivx].Z;
      newtrk.VtxIndex[attachEnd] = ivx;
      // correct the trajectory direction
      TVector3 dir;
      if (itj == 0) {
        dir = newtrk.TrjPos[1] - newtrk.TrjPos[0];
        newtrk.TrjDir[0] = dir.Unit();
      }
      else {
        dir = newtrk.TrjPos[itj - 1] - newtrk.TrjPos[itj];
        newtrk.TrjDir[itj] = dir.Unit();
      }
    } // end

    if (newtrk.VtxIndex[0] >= 0 && newtrk.VtxIndex[0] == newtrk.VtxIndex[1]) {
      mf::LogError("CCTM")
        << "StoreTrack: Trying to attach a vertex to both ends of a track. imat = " << imat;
      return;
    }

    // calculate the length
    newtrk.Length = 0;
    float norm;
    double X, Y, Z;
    for (unsigned short itj = 1; itj < newtrk.TrjPos.size(); ++itj) {
      X = newtrk.TrjPos[itj](0) - newtrk.TrjPos[itj - 1](0);
      Y = newtrk.TrjPos[itj](1) - newtrk.TrjPos[itj - 1](1);
      Z = newtrk.TrjPos[itj](2) - newtrk.TrjPos[itj - 1](2);
      norm = sqrt(X * X + Y * Y + Z * Z);
      newtrk.Length += norm;
    }

    // store the cluster -> track assignment
    newtrk.ClsEvtIndices.clear();
    for (unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      if (matcomb[imat].Cls[ipl] < 0) continue;
      ccl = matcomb[imat].Cls[ipl];
      if (ccl > clsChain[ipl].size()) std::cout << "oops StoreTrack\n";
      clsChain[ipl][ccl].InTrack = newtrk.ID;
      for (unsigned short icc = 0; icc < clsChain[ipl][ccl].ClsIndex.size(); ++icc) {
        icl = clsChain[ipl][ccl].ClsIndex[icc];
        if (icl > cls[ipl].size()) std::cout << "oops StoreTrack\n";
        cls[ipl][icl].InTrack = newtrk.ID;
        if (cls[ipl][icl].EvtIndex > fmCluHits.size() - 1) {
          std::cout << "ooops2 store track EvtIndex " << cls[ipl][icl].EvtIndex << " size "
                    << fmCluHits.size() << " icl " << icl << "\n";
          continue;
        }
        newtrk.ClsEvtIndices.push_back(cls[ipl][icl].EvtIndex);
      } // icc
    }   // ipl

    if (prt) mf::LogVerbatim("CCTM") << " track ID " << newtrk.ID << " stored in StoreTrack";

    trk.push_back(newtrk);
  } // StoreTrack

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::PlnMatch(detinfo::DetectorPropertiesData const& detProp,
                              geo::TPCID const& tpcid)
  {
    // Match clusters in all planes
    bool ignoreSign;
    float kSlp, kAng, kX, kWir, okWir;
    short idir, ioend, jdir, joend, kdir;

    double yp, zp;
    auto const& first_tpc = geom->TPC({0, 0});
    float tpcSizeY = first_tpc.HalfWidth();
    float tpcSizeZ = first_tpc.Length();

    float dxcut = 2;
    float dxkcut;
    float dwcut = 6;
    if (fuBCode) {
      dxcut = 20;
      dwcut = 60;
    }

    // temp array for making a rough charge asymmetry cut
    std::array<float, 3> mchg;
    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    for (unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      geo::PlaneID const iplane_id{tpcid, ipl};
      for (unsigned short icl = 0; icl < clsChain[ipl].size(); ++icl) {
        if (clsChain[ipl][icl].InTrack >= 0) continue;
        // skip short clusters
        if (clsChain[ipl][icl].Length < fMatchMinLen[algIndex]) continue;
        unsigned short jpl = (ipl + 1) % nplanes;
        unsigned short kpl = (jpl + 1) % nplanes;
        geo::PlaneID const jplane_id{tpcid, jpl};
        for (unsigned short jcl = 0; jcl < clsChain[jpl].size(); ++jcl) {
          if (clsChain[jpl][jcl].InTrack >= 0) continue;
          // skip short clusters
          if (clsChain[jpl][jcl].Length < fMatchMinLen[algIndex]) continue;
          // make first charge asymmetry cut
          mchg[0] = clsChain[ipl][icl].TotChg;
          mchg[1] = clsChain[jpl][jcl].TotChg;
          mchg[2] = mchg[1];
          if (fChgAsymFactor[algIndex] > 0 && ChargeAsym(mchg) > 0.5) continue;
          for (unsigned short iend = 0; iend < 2; ++iend) {
            idir = clsChain[ipl][icl].Dir[iend];
            for (unsigned short jend = 0; jend < 2; ++jend) {
              jdir = clsChain[jpl][jcl].Dir[jend];
              if (idir != 0 && jdir != 0 && idir != jdir) continue;
              // make an X cut
              if (fabs(clsChain[ipl][icl].X[iend] - clsChain[jpl][jcl].X[jend]) > dxcut) continue;
              ioend = 1 - iend;
              joend = 1 - jend;
              // Find the expected third (k) plane parameters
              kSlp = wireReadoutGeom.ThirdPlaneSlope(
                ipl, clsChain[ipl][icl].Slope[iend], jpl, clsChain[jpl][jcl].Slope[jend], tpcid);
              kAng = atan(kSlp);

              // Ensure the match end is within the TPC
              auto const intersection =
                wireReadoutGeom
                  .WireIDsIntersect(
                    geo::WireID{iplane_id, (unsigned int)(0.5 + clsChain[ipl][icl].Wire[iend])},
                    geo::WireID{jplane_id, (unsigned int)(0.5 + clsChain[jpl][jcl].Wire[jend])})
                  .value_or(geo::WireIDIntersection::invalid());
              yp = intersection.y;
              zp = intersection.z;
              if (yp > tpcSizeY || yp < -tpcSizeY) continue;
              if (zp < 0 || zp > tpcSizeZ) continue;

              kX = 0.5 * (clsChain[ipl][icl].X[iend] + clsChain[jpl][jcl].X[jend]);
              auto const& kplane = wireReadoutGeom.Plane({tpcid, kpl});
              kWir = kplane.WireCoordinate(geo::Point_t{0, yp, zp});

              // now look at the other end
              auto const ointersection =
                wireReadoutGeom
                  .WireIDsIntersect(
                    geo::WireID{iplane_id, (unsigned int)(0.5 + clsChain[ipl][icl].Wire[ioend])},
                    geo::WireID{jplane_id, (unsigned int)(0.5 + clsChain[jpl][jcl].Wire[joend])})
                  .value_or(geo::WireIDIntersection::invalid());
              yp = ointersection.y;
              zp = ointersection.z;
              if (yp > tpcSizeY || yp < -tpcSizeY) continue;
              if (zp < 0 || zp > tpcSizeZ) continue;
              okWir = kplane.WireCoordinate(geo::Point_t{0, yp, zp});

              if (prt)
                mf::LogVerbatim("CCTM")
                  << "PlnMatch: chk i " << ipl << ":" << icl << ":" << iend << " idir " << idir
                  << " X " << clsChain[ipl][icl].X[iend] << " j " << jpl << ":" << jcl << ":"
                  << jend << " jdir " << jdir << " X " << clsChain[jpl][jcl].X[jend];

              if (prt)
                mf::LogVerbatim("CCTM")
                  << "PlnMatch: chk j " << ipl << ":" << icl << ":" << iend << " " << jpl << ":"
                  << jcl << ":" << jend << " iSlp " << std::setprecision(2)
                  << clsChain[ipl][icl].Slope[iend] << " jSlp " << std::setprecision(2)
                  << clsChain[jpl][jcl].Slope[jend] << " kWir " << (int)kWir << " okWir "
                  << (int)okWir << " kSlp " << std::setprecision(2) << kSlp << " kAng "
                  << std::setprecision(2) << kAng << " kX " << std::setprecision(1) << kX;

              // handle the case near pi/2, where the errors on large slopes
              // could result in a wrong-sign kAng
              ignoreSign = (fabs(kSlp) > 1.5);
              if (ignoreSign) kAng = fabs(kAng);
              dxkcut = dxcut * AngleFactor(kSlp);
              bool gotkcl = false;
              for (unsigned short kcl = 0; kcl < clsChain[kpl].size(); ++kcl) {
                if (clsChain[kpl][kcl].InTrack >= 0) continue;
                // make second charge asymmetry cut
                mchg[0] = clsChain[ipl][icl].TotChg;
                mchg[1] = clsChain[jpl][jcl].TotChg;
                mchg[2] = clsChain[kpl][kcl].TotChg;
                if (fChgAsymFactor[algIndex] > 0 && ChargeAsym(mchg) > 0.5) continue;
                for (unsigned short kend = 0; kend < 2; ++kend) {
                  kdir = clsChain[kpl][kcl].Dir[kend];
                  if (idir != 0 && kdir != 0 && idir != kdir) continue;
                  if (prt)
                    mf::LogVerbatim("CCTM")
                      << " kcl " << kcl << " kend " << kend << " dx "
                      << std::abs(clsChain[kpl][kcl].X[kend] - kX) << " dxkcut " << dxkcut;
                  if (std::abs(clsChain[kpl][kcl].X[kend] - kX) > dxkcut) continue;
                  // rough dWire cut
                  if (prt)
                    mf::LogVerbatim("CCTM")
                      << " kcl " << kcl << " kend " << kend << " dw "
                      << (clsChain[kpl][kcl].Wire[kend] - kWir) << " ignoreSign " << ignoreSign;
                  if (fabs(clsChain[kpl][kcl].Wire[kend] - kWir) > dwcut) continue;
                  if (prt) mf::LogVerbatim("CCTM") << " chk k " << kpl << ":" << kcl << ":" << kend;
                  MatchPars match;
                  match.Cls[ipl] = icl;
                  match.End[ipl] = iend;
                  match.Cls[jpl] = jcl;
                  match.End[jpl] = jend;
                  match.Cls[kpl] = kcl;
                  match.End[kpl] = kend;
                  match.Err = 100;
                  if (DupMatch(match, nplanes)) continue;
                  match.Chg[ipl] = 0;
                  match.Chg[jpl] = 0;
                  match.Chg[kpl] = 0;
                  match.Vtx = clsChain[ipl][icl].VtxIndex[iend];
                  match.oVtx = -1;
                  FillEndMatch(detProp, tpcid, match);
                  if (prt)
                    mf::LogVerbatim("CCTM") << " PlnMatch: match k " << kpl << ":" << match.Cls[kpl]
                                            << ":" << match.End[kpl] << " oChg " << match.Chg[kpl]
                                            << " mErr " << match.Err << " oErr " << match.oErr;
                  if (match.Chg[kpl] == 0) continue;
                  if (match.Err > 10 || match.oErr > 10) continue;
                  if (prt) mf::LogVerbatim("CCTM") << " dup? ";
                  if (DupMatch(match, nplanes)) continue;
                  matcomb.push_back(match);
                  gotkcl = true;
                } // kend
              }   // kcl
              if (prt) mf::LogVerbatim("CCTM") << " PlnMatch: gotkcl " << gotkcl;
              if (!gotkcl) {
                // make a 2-plane match and try again
                MatchPars match;
                match.Cls[ipl] = icl;
                match.End[ipl] = iend;
                match.Cls[jpl] = jcl;
                match.End[jpl] = jend;
                match.Cls[kpl] = -1;
                match.End[kpl] = 0;
                match.Err = 100;
                if (DupMatch(match, nplanes)) continue;
                match.Chg[ipl] = 0;
                match.Chg[jpl] = 0;
                match.Chg[kpl] = 0;
                match.Vtx = clsChain[ipl][icl].VtxIndex[iend];
                match.oVtx = -1;
                FillEndMatch(detProp, tpcid, match);
                if (prt)
                  mf::LogVerbatim("CCTM")
                    << " Tried 2-plane match"
                    << " k " << kpl << ":" << match.Cls[kpl] << ":" << match.End[kpl] << " Chg "
                    << match.Chg[kpl] << " Err " << match.Err << " match.oErr " << match.oErr;
                if (match.Chg[kpl] <= 0) continue;
                if (match.Err > 10 || match.oErr > 10) continue;
                matcomb.push_back(match);
              } // !gotkcl
            }   // jend
          }     // iend
        }       // jcl
      }         // icl
    }           // ipl
  }             // PlnMatch

  ///////////////////////////////////////////////////////////////////////
  bool CCTrackMaker::DupMatch(MatchPars& match, unsigned short const nplanes)
  {
    unsigned short nMatCl, nMiss;
    float toterr = match.Err + match.oErr;
    for (unsigned int imat = 0; imat < matcomb.size(); ++imat) {
      // check for exact matches
      if (match.Cls[0] == matcomb[imat].Cls[0] && match.Cls[1] == matcomb[imat].Cls[1] &&
          match.Cls[2] == matcomb[imat].Cls[2]) {

        // compare the error
        if (toterr < matcomb[imat].Err + matcomb[imat].oErr) {
          // keep the better one
          matcomb[imat].End[0] = match.End[0];
          matcomb[imat].End[1] = match.End[1];
          matcomb[imat].End[2] = match.End[2];
          matcomb[imat].Vtx = match.Vtx;
          matcomb[imat].dWir = match.dWir;
          matcomb[imat].dAng = match.dAng;
          matcomb[imat].dX = match.dX;
          matcomb[imat].Err = match.Err;
          matcomb[imat].oVtx = match.oVtx;
          matcomb[imat].odWir = match.odWir;
          matcomb[imat].odAng = match.odAng;
          matcomb[imat].odX = match.odX;
          matcomb[imat].oErr = match.oErr;
        }
        return true;
      } // test
      // check for a 3-plane match vs 2-plane match
      nMatCl = 0;
      nMiss = 0;
      for (unsigned short ipl = 0; ipl < nplanes; ++ipl) {
        if (match.Cls[ipl] >= 0) {
          if (match.Cls[ipl] == matcomb[imat].Cls[ipl] &&
              (match.End[0] == matcomb[imat].End[0] || match.End[1] == matcomb[imat].End[1]))
            ++nMatCl;
        }
        else {
          ++nMiss;
        }
      } // ipl
      if (nMatCl == 2 && nMiss == 1) return true;
    } // imat
    return false;
  } // DupMatch

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::SortMatches(detinfo::DetectorPropertiesData const& detProp,
                                 art::FindManyP<recob::Hit> const& fmCluHits,
                                 unsigned short const procCode,
                                 geo::TPCID const& tpcid)
  {
    // sort cluster matches by increasing total match error. Find the minimum total error of all
    // cluster match combinations and make tracks from them
    CluLen merr;
    std::vector<CluLen> materr;
    unsigned int ii, im;

    if (matcomb.size() == 0) return;

    // sort by decreasing error
    for (ii = 0; ii < matcomb.size(); ++ii) {
      merr.index = ii;
      merr.length = matcomb[ii].Err + matcomb[ii].oErr;
      materr.push_back(merr);
    } // ii
    std::sort(materr.begin(), materr.end(), lessThan);

    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    if (prt) {
      mf::LogVerbatim myprt("CCTM");
      myprt << "SortMatches\n";
      myprt << "   ii    im  Vx   Err     dW     dA     dX  oVx   oErr    odW   odA    odX   Asym  "
               " icl   jcl   kcl \n";
      for (ii = 0; ii < materr.size(); ++ii) {
        im = materr[ii].index;
        float asym =
          fabs(matcomb[im].Chg[0] - matcomb[im].Chg[1]) / (matcomb[im].Chg[0] + matcomb[im].Chg[1]);
        asym *=
          fabs(matcomb[im].Chg[1] - matcomb[im].Chg[2]) / (matcomb[im].Chg[1] + matcomb[im].Chg[2]);
        myprt << std::fixed << std::right << std::setw(5) << ii << std::setw(5) << im
              << std::setw(4) << matcomb[im].Vtx << std::setw(7) << std::setprecision(2)
              << matcomb[im].Err << std::setw(7) << std::setprecision(1) << matcomb[im].dWir
              << std::setw(7) << std::setprecision(2) << matcomb[im].dAng << std::setw(7)
              << std::setprecision(2) << matcomb[im].dX << std::setw(4) << matcomb[im].oVtx
              << std::setw(7) << std::setprecision(2) << matcomb[im].oErr << std::setw(7)
              << std::setprecision(1) << matcomb[im].odWir << std::setw(7) << std::setprecision(2)
              << matcomb[im].odAng << std::setw(7) << std::setprecision(2) << matcomb[im].odX
              << std::setw(7) << std::setprecision(3) << asym << " 0:" << matcomb[im].Cls[0] << ":"
              << matcomb[im].End[0] << " 1:" << matcomb[im].Cls[1] << ":" << matcomb[im].End[1];
        if (nplanes > 2) myprt << " 2:" << matcomb[im].Cls[2] << ":" << matcomb[im].End[2];
        myprt << "\n";
      } // ii
    }   // prt

    // define an array to ensure clusters are only used once
    std::array<std::vector<bool>, 3> pclUsed;
    unsigned short ipl;
    for (ipl = 0; ipl < nplanes; ++ipl) {
      pclUsed[ipl].resize(clsChain[ipl].size());
    }

    // count the total number of clusters and length used in matcomb
    unsigned short matcombTotCl = 0;
    float matcombTotLen = 0;
    unsigned short icl;
    for (ii = 0; ii < matcomb.size(); ++ii) {
      for (ipl = 0; ipl < nplanes; ++ipl) {
        if (matcomb[ii].Cls[ipl] < 0) continue;
        icl = matcomb[ii].Cls[ipl];
        ++matcombTotCl;
        matcombTotLen += clsChain[ipl][icl].Length;
      }
    }

    if (prt)
      mf::LogVerbatim("CCTM") << "Number of clusters to match " << matcombTotCl << " total length "
                              << matcombTotLen;

    if (matcombTotLen <= 0) {
      mf::LogError("CCTM") << "SortMatches: bad matcomb total length " << matcombTotLen;
      return;
    }

    // vector of matcomb indices of unique cluster matches
    std::vector<unsigned short> matIndex;
    // vector of matcomb indices of unique cluster matches that have the best total error
    std::vector<unsigned short> bestMatIndex;
    float totLen, totErr, bestTotErr = 9999;
    // start with the best match
    unsigned short jj, jm, nused, jcl;
    // fraction of the length of all clustters in matcomb that are used in a match
    float fracLen;

    for (ii = 0; ii < materr.size(); ++ii) {
      im = materr[ii].index;
      matIndex.clear();
      // skip really bad matches
      if (matcomb[im].Err > bestTotErr) continue;
      totLen = 0;
      // initialize pclUsed and flag the clusters in this match
      for (ipl = 0; ipl < nplanes; ++ipl) {
        // initialize to no clusters used
        std::fill(pclUsed[ipl].begin(), pclUsed[ipl].end(), false);
        // check for 2 plane match
        if (matcomb[im].Cls[ipl] < 0) continue;
        icl = matcomb[im].Cls[ipl];
        pclUsed[ipl][icl] = true;
        totLen += clsChain[ipl][icl].Length;
      } // ipl
      // Initialize the error sum
      totErr = matcomb[im].Err;
      // Save the index
      matIndex.push_back(im);
      // look for matches in the rest of the list that are not already matched.
      for (jj = 0; jj < materr.size(); ++jj) {
        if (jj == ii) continue;
        jm = materr[jj].index;
        // skip really bad matches
        if (matcomb[jm].Err > bestTotErr) continue;
        // check for non-unique cluster indices
        nused = 0;
        for (ipl = 0; ipl < nplanes; ++ipl) {
          if (matcomb[jm].Cls[ipl] < 0) continue;
          jcl = matcomb[jm].Cls[ipl];
          if (pclUsed[ipl][jcl]) ++nused;
          // This cluster chain was used in a previous match
          if (nused > 0) break;
          totLen += clsChain[ipl][jcl].Length;
        } // ipl
        // at least one of the clusters in this match have been used
        if (nused != 0) continue;
        // found a match with an unmatched set of clusters. Update the total error and flag them
        totErr += matcomb[jm].Err;
        matIndex.push_back(jm);
        // Flag the clusters used and see if all of them are used
        for (ipl = 0; ipl < nplanes; ++ipl) {
          if (matcomb[jm].Cls[ipl] < 0) continue;
          jcl = matcomb[jm].Cls[ipl];
          pclUsed[ipl][jcl] = true;
        } // ipl
      }   // jm
      if (totLen == 0) continue;
      nused = 0;
      for (ipl = 0; ipl < nplanes; ++ipl) {
        for (unsigned short indx = 0; indx < pclUsed[ipl].size(); ++indx)
          if (pclUsed[ipl][indx]) ++nused;
      } // ipl
      if (totLen > matcombTotLen) std::cout << "Oops " << totLen << " " << matcombTotLen << "\n";
      // weight the total error by the total length of all clusters
      fracLen = totLen / matcombTotLen;
      totErr /= fracLen;
      if (prt) {
        mf::LogVerbatim myprt("CCTM");
        myprt << "match " << im << " totErr " << totErr << " nused " << nused << " fracLen "
              << fracLen << " totLen " << totLen << " mat: ";
        for (unsigned short indx = 0; indx < matIndex.size(); ++indx)
          myprt << " " << matIndex[indx];
      } // prt
      // check for more used clusters and a better total error
      if (totErr < bestTotErr) {
        bestTotErr = totErr;
        bestMatIndex = matIndex;
        if (nused == matcombTotCl) break;
        if (prt) {
          mf::LogVerbatim myprt("CCTM");
          myprt << "bestTotErr " << bestTotErr << " nused " << nused << " matcombTotCl "
                << matcombTotCl << " mat: ";
          for (unsigned short indx = 0; indx < bestMatIndex.size(); ++indx)
            myprt << " " << bestMatIndex[indx];
        } // prt
        // stop looking if we have found everything
        if (fracLen > 0.999) break;
      } // totErr < bestTotErr
    }   // im

    if (bestTotErr > 9000) return;

    for (ii = 0; ii < bestMatIndex.size(); ++ii) {
      im = bestMatIndex[ii];
      FillTrkHits(fmCluHits, im, nplanes);
      // look for missing clusters?
      // store this track with processor code 1
      StoreTrack(detProp, fmCluHits, im, procCode, tpcid);
    } // ii

  } // SortMatches

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillEndMatch2(MatchPars& match)
  {
    // 2D version of FillEndMatch

    match.Err = 100;
    match.oErr = 100;
    match.Chg[2] = 0;
    match.dWir = 0;
    match.dAng = 0;
    match.odWir = 0;
    match.odAng = 0;

    unsigned short ipl = 0;
    unsigned short jpl = 1;

    if (match.Cls[0] < 0 || match.Cls[1] < 0) return;

    unsigned short icl = match.Cls[ipl];
    unsigned short iend = match.End[ipl];
    match.Chg[ipl] = clsChain[ipl][icl].TotChg;
    // cluster i match end
    float miX = clsChain[ipl][icl].X[iend];
    // cluster i other end
    unsigned short oiend = 1 - iend;
    float oiX = clsChain[ipl][icl].X[oiend];

    unsigned short jcl = match.Cls[jpl];
    unsigned short jend = match.End[jpl];
    match.Chg[jpl] = clsChain[jpl][jcl].TotChg;
    // cluster j match end
    float mjX = clsChain[jpl][jcl].X[jend];
    // cluster j other end
    unsigned short ojend = 1 - jend;
    float ojX = clsChain[jpl][jcl].X[ojend];

    // look for a match end vertex match
    match.Vtx = -1;
    if (clsChain[ipl][icl].VtxIndex[iend] >= 0 &&
        clsChain[ipl][icl].VtxIndex[iend] == clsChain[jpl][jcl].VtxIndex[jend]) {
      match.Vtx = clsChain[ipl][icl].VtxIndex[iend];
      miX = vtx[match.Vtx].X;
      mjX = vtx[match.Vtx].X;
    }

    // look for an other end vertex match
    match.oVtx = -1;
    if (clsChain[ipl][icl].VtxIndex[oiend] >= 0 &&
        clsChain[ipl][icl].VtxIndex[oiend] == clsChain[jpl][jcl].VtxIndex[ojend]) {
      match.oVtx = clsChain[ipl][icl].VtxIndex[oiend];
      oiX = vtx[match.oVtx].X;
      ojX = vtx[match.oVtx].X;
    }

    // find the charge asymmetry
    float chgAsym = 1;
    if (fChgAsymFactor[algIndex] > 0) {
      chgAsym = fabs(match.Chg[ipl] - match.Chg[jpl]) / (match.Chg[ipl] + match.Chg[jpl]);
      if (chgAsym > 0.5) return;
      chgAsym = 1 + fChgAsymFactor[algIndex] * chgAsym;
    }

    // find the error at the match end
    float maxSlp = fabs(clsChain[ipl][icl].Slope[iend]);
    if (fabs(clsChain[jpl][jcl].Slope[jend]) > maxSlp)
      maxSlp = fabs(clsChain[jpl][jcl].Slope[jend]);
    float sigmaX = fXMatchErr[algIndex] + std::max(maxSlp, (float)20);
    match.dX = fabs(miX - mjX);
    match.Err = chgAsym * match.dX / sigmaX;

    // find the error at the other end
    maxSlp = fabs(clsChain[ipl][icl].Slope[oiend]);
    if (fabs(clsChain[jpl][jcl].Slope[ojend]) > maxSlp)
      maxSlp = fabs(clsChain[jpl][jcl].Slope[ojend]);
    sigmaX = fXMatchErr[algIndex] + std::max(maxSlp, (float)20);
    match.odX = fabs(oiX - ojX);
    match.oErr = chgAsym * match.odX / sigmaX;

    if (prt)
      mf::LogVerbatim("CCTM") << "FEM2: m " << ipl << ":" << icl << ":" << iend << " miX " << miX
                              << " - " << jpl << ":" << jcl << ":" << jend << " mjX " << mjX
                              << " match.dX " << match.dX << " match.Err " << match.Err
                              << " chgAsym " << chgAsym << " o "
                              << " oiX " << oiX << " ojX " << ojX << " match.odX " << match.odX
                              << " match.oErr " << match.oErr << "\n";

  } // FillEndMatch2

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillEndMatch(detinfo::DetectorPropertiesData const& detProp,
                                  geo::TPCID const& tpcid,
                                  MatchPars& match)
  {
    // fill the matching parameters for this cluster match. The calling routine
    // should set the match end vertex ID (if applicable) as well as the
    // cluster IDs and matching ends in each plane. Note that the matching variables
    // Note that dWir, dAng and dTim are not filled if there is a vertex (match.Vtx >= 0).
    // Likewise, odWir, odAng and odX are not filled if there is a vertex match
    // at the other end
    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);

    if (nplanes == 2) {
      FillEndMatch2(match);
      return;
    }

    std::array<short, 3> mVtx;
    std::array<short, 3> oVtx;
    std::array<float, 3> oWir;
    std::array<float, 3> oSlp;
    std::array<float, 3> oAng;
    std::array<float, 3> oX;

    std::array<float, 3> mChg;

    unsigned short ii, ipl, iend, jpl, jend, kpl, kend, oend;
    short icl, jcl, kcl;

    for (ipl = 0; ipl < 3; ++ipl) {
      mVtx[ipl] = -1;
      oVtx[ipl] = -1;
      oWir[ipl] = -66;
      oSlp[ipl] = -66;
      oAng[ipl] = -66;
      oX[ipl] = -66;
      mChg[ipl] = -1;
    } // ipl

    // initialize parameters that shouldn't have been set by the calling routine
    match.dWir = 0;
    match.dAng = 0;
    match.dX = 0;
    match.Err = 100;
    match.odWir = 0;
    match.odAng = 0;
    match.odX = 0;
    match.oErr = 100;
    match.oVtx = -1;

    if (prt) {
      mf::LogVerbatim myprt("CCTM");
      myprt << "FEM ";
      for (ipl = 0; ipl < nplanes; ++ipl) {
        myprt << " " << ipl << ":" << match.Cls[ipl] << ":" << match.End[ipl];
      }
    }

    short missingPlane = -1;
    unsigned short nClInPln = 0;
    // number of vertex matches at each end
    short aVtx = -1;
    unsigned short novxmat = 0;
    short aoVtx = -1;
    unsigned short nvxmat = 0;
    unsigned short nShortCl = 0;
    // fill the other end parameters in each plane
    for (ipl = 0; ipl < nplanes; ++ipl) {
      if (match.Cls[ipl] < 0) {
        missingPlane = ipl;
        continue;
      }
      ++nClInPln;
      icl = match.Cls[ipl];
      match.Chg[ipl] = clsChain[ipl][icl].TotChg;
      mChg[ipl] = clsChain[ipl][icl].TotChg;
      iend = match.End[ipl];
      mVtx[ipl] = clsChain[ipl][icl].VtxIndex[iend];
      if (clsChain[ipl][icl].Length < 6) ++nShortCl;
      if (mVtx[ipl] >= 0) {
        if (aVtx < 0) aVtx = mVtx[ipl];
        if (mVtx[ipl] == aVtx) ++nvxmat;
      }
      if (prt)
        mf::LogVerbatim("CCTM") << "FEM: m " << ipl << ":" << icl << ":" << iend << " Vtx "
                                << mVtx[ipl] << "  Wir " << clsChain[ipl][icl].Wire[iend]
                                << std::fixed << std::setprecision(3) << " Slp "
                                << clsChain[ipl][icl].Slope[iend] << std::fixed
                                << std::setprecision(1) << " X " << clsChain[ipl][icl].X[iend];

      oend = 1 - iend;
      oWir[ipl] = clsChain[ipl][icl].Wire[oend];
      oAng[ipl] = clsChain[ipl][icl].Angle[oend];
      oSlp[ipl] = clsChain[ipl][icl].Slope[oend];
      oX[ipl] = clsChain[ipl][icl].X[oend];
      oVtx[ipl] = clsChain[ipl][icl].VtxIndex[oend];
      if (oVtx[ipl] >= 0) {
        if (aoVtx < 0) aoVtx = oVtx[ipl];
        if (oVtx[ipl] == aoVtx) ++novxmat;
      }

      if (prt)
        mf::LogVerbatim("CCTM") << "     o " << ipl << ":" << icl << ":" << oend << " oVtx "
                                << oVtx[ipl] << " oWir " << oWir[ipl] << std::fixed
                                << std::setprecision(3) << " oSlp " << oSlp[ipl] << std::fixed
                                << std::setprecision(1) << " oX " << oX[ipl] << " Chg "
                                << (int)mChg[ipl];

    } // ipl

    bool isShort = (nShortCl > 1);

    if (nClInPln < 2) {
      mf::LogWarning("CCTM") << "Not enough matched planes supplied";
      return;
    }

    if (prt)
      mf::LogVerbatim("CCTM") << "FEM: Vtx m " << aVtx << " count " << nvxmat << " o " << aoVtx
                              << " count " << novxmat << " missingPlane " << missingPlane
                              << " nClInPln " << nClInPln;

    // perfect match
    if (nvxmat == 3 && novxmat == 3) {
      match.Vtx = aVtx;
      match.Err = 0;
      match.oVtx = aoVtx;
      match.oErr = 0;
      return;
    }

    // 2-plane vertex match?
    // factors applied to error = 1 (no vtx), 0.5 (2 pln vtx), 0.33 (3 pln vtx)
    float vxFactor = 1;
    float ovxFactor = 1;
    if (nClInPln == 3) {
      // a cluster in all 3 planes
      if (nvxmat == 3) {
        // and all vertex assignments agree at the match end
        match.Vtx = aVtx;
        vxFactor = 0.33;
      }
      if (novxmat == 3) {
        // and all vertex assignments agree at the other end
        match.oVtx = aoVtx;
        ovxFactor = 0.33;
      }
    }
    else {
      // a cluster in 2 planes
      if (nvxmat == 2) {
        match.Vtx = aVtx;
        vxFactor = 0.5;
      }
      if (novxmat == 2) {
        match.oVtx = aoVtx;
        ovxFactor = 0.5;
      }
    } // nClInPln

    // find wire, X and Time at both ends

    // Find the "other end" matching parameters with
    // two cases: a 3-plane match or a 2-plane match
    // and with/without an other end vertex

    float kWir, okWir;
    float kSlp, okSlp, kAng, okAng, okX, kX, kTim, okTim;

    if (nClInPln == 3) {
      ipl = 0;
      jpl = 1;
      kpl = 2;
    }
    else {
      // 2-plane match
      kpl = missingPlane;
      if (kpl == 0) {
        ipl = 1;
        jpl = 2;
      }
      else if (kpl == 1) {
        ipl = 2;
        jpl = 0;
      }
      else {
        ipl = 0;
        jpl = 1;
      } // kpl test
    }   // missing plane
    iend = match.End[ipl];
    jend = match.End[jpl];
    icl = match.Cls[ipl];
    jcl = match.Cls[jpl];
    if (nplanes > 2) {
      kcl = match.Cls[kpl];
      kend = match.End[kpl];
    }

    geo::PlaneID const iplaneid{tpcid, ipl};
    geo::PlaneID const jplaneid{tpcid, jpl};
    geo::PlaneID const kplaneid{tpcid, kpl};
    auto const& kplane = wireReadoutGeom.Plane(kplaneid);
    /////////// Wire, Angle, X and Time at the Other end
    okSlp = wireReadoutGeom.ThirdPlaneSlope(ipl, oSlp[ipl], jpl, oSlp[jpl], tpcid);
    okAng = atan(okSlp);
    // handle the case near pi/2, where the errors on large slopes could result in
    // a wrong-sign kAng
    bool ignoreSign = (fabs(okSlp) > 10);
    if (ignoreSign) okAng = fabs(okAng);
    if (match.oVtx >= 0) {
      // a vertex exists at the other end
      okWir = kplane.WireCoordinate(geo::Point_t{0, vtx[match.oVtx].Y, vtx[match.oVtx].Z});
      okX = vtx[match.oVtx].X;
    }
    else {
      // no vertex at the other end
      auto const [ypos, zpos, _] =
        wireReadoutGeom
          .WireIDsIntersect(geo::WireID(iplaneid, oWir[ipl]), geo::WireID(jplaneid, oWir[jpl]))
          .value_or(geo::WireIDIntersection::invalid());
      okWir = 0.5 + kplane.WireCoordinate(geo::Point_t{0, ypos, zpos});
      okX = 0.5 * (oX[ipl] + oX[jpl]);
    }
    okTim = detProp.ConvertXToTicks(okX, kplaneid);
    if (prt)
      mf::LogVerbatim("CCTM") << "FEM: oEnd"
                              << " kpl " << kpl << " okSlp " << okSlp << " okAng " << okAng
                              << " okWir " << (int)okWir << " okX " << okX;

    /////////// Wire, Angle, X and Time at the Match end

    kSlp = wireReadoutGeom.ThirdPlaneSlope(
      ipl, clsChain[ipl][icl].Slope[iend], jpl, clsChain[jpl][jcl].Slope[jend], kplaneid);
    kAng = atan(kSlp);
    if (ignoreSign) kAng = fabs(kAng);
    if (match.Vtx >= 0) {
      if (vtx.size() == 0 || (unsigned int)match.Vtx > vtx.size() - 1) {
        mf::LogError("CCTM") << "FEM: Bad match.Vtx " << match.Vtx << " vtx size " << vtx.size();
        return;
      }
      // a vertex exists at the match end
      kWir = kplane.WireCoordinate(geo::Point_t{0, vtx[match.Vtx].Y, vtx[match.Vtx].Z});
      kX = vtx[match.Vtx].X;
    }
    else {
      // no vertex at the match end
      auto const [ypos, zpos, _] =
        wireReadoutGeom
          .WireIDsIntersect(geo::WireID(iplaneid, clsChain[ipl][icl].Wire[iend]),
                            geo::WireID(jplaneid, clsChain[jpl][jcl].Wire[jend]))
          .value_or(geo::WireIDIntersection::invalid());
      kWir = 0.5 + kplane.WireCoordinate(geo::Point_t{0, ypos, zpos});
      kX = 0.5 * (clsChain[ipl][icl].X[iend] + clsChain[jpl][jcl].X[jend]);
    }
    kTim = detProp.ConvertXToTicks(kX, kplaneid);
    if (prt)
      mf::LogVerbatim("CCTM") << "FEM: mEnd"
                              << " kpl " << kpl << " kSlp " << kSlp << " kAng " << kAng << " kX "
                              << kX;

    // try to find a 3-plane match using this information
    if (nClInPln < 3 && FindMissingCluster(kpl, kcl, kend, kWir, kX, okWir, okX)) {
      nClInPln = 3;
      // update local variables
      match.Cls[kpl] = kcl;
      match.End[kpl] = kend;
      match.Chg[kpl] = clsChain[kpl][kcl].TotChg;
      mChg[kpl] = clsChain[kpl][kcl].TotChg;
      oend = 1 - kend;
      oWir[kpl] = clsChain[kpl][kcl].Wire[oend];
      oX[kpl] = clsChain[kpl][kcl].X[oend];
      oAng[kpl] = clsChain[kpl][kcl].Angle[oend];
      oSlp[kpl] = clsChain[kpl][kcl].Slope[oend];
    } // FindMissingCluster

    // decide whether to continue with a 2-plane match. The distance between match and other end should
    // be large enough to create a cluster
    if (nClInPln == 2 && fabs(okWir - kWir) > 3) return;

    // Calculate the cluster charge asymmetry. This factor will be applied
    // to the error of the end matches
    float chgAsym = 1;
    // Get the charge in the plane without a matching cluster
    if (nClInPln < 3 && mChg[missingPlane] <= 0) {
      if (missingPlane != kpl)
        mf::LogError("CCTM") << "FEM bad missingPlane " << missingPlane << " " << kpl << "\n";
      mChg[kpl] = ChargeNear(kplaneid, (unsigned short)kWir, kTim, (unsigned short)okWir, okTim);
      match.Chg[kpl] = mChg[kpl];
      if (prt)
        mf::LogVerbatim("CCTM") << "FEM: Missing cluster in " << kpl << " ChargeNear " << (int)kWir
                                << ":" << (int)kTim << " " << (int)okWir << ":" << (int)okTim
                                << " chg " << mChg[kpl];
      if (mChg[kpl] <= 0) return;
    }

    if (fChgAsymFactor[algIndex] > 0) {
      chgAsym = ChargeAsym(mChg);
      if (chgAsym > 0.5) return;
      chgAsym = 1 + fChgAsymFactor[algIndex] * chgAsym;
    }

    if (prt) mf::LogVerbatim("CCTM") << "FEM: charge asymmetry factor " << chgAsym;
    float sigmaX, sigmaA;
    float da, dx, dw;

    /////////// Matching error at the Match end
    // check for vertex consistency at the match end
    aVtx = -1;
    bool allPlnVtxMatch = false;
    if (nClInPln == 3) {
      unsigned short nmvtx = 0;
      for (ii = 0; ii < nplanes; ++ii) {
        if (mVtx[ii] >= 0) {
          if (aVtx < 0) aVtx = mVtx[ii];
          ++nmvtx;
        }
      } // ii
      // same vertex in all planes
      if (nmvtx) allPlnVtxMatch = true;
    } // nClInPln

    // inflate the X match error to allow for missing one wire on the end of a cluster
    sigmaX = fXMatchErr[algIndex] + std::max(kSlp, (float)20);
    sigmaA = fAngleMatchErr[algIndex] * AngleFactor(kSlp);
    if (prt)
      mf::LogVerbatim("CCTM") << "bb " << algIndex << "  " << fXMatchErr[algIndex] << " "
                              << fAngleMatchErr[algIndex] << " kslp " << kSlp;

    if (nClInPln == 3) {
      kcl = match.Cls[kpl];
      kend = match.End[kpl];
      dw = kWir - clsChain[kpl][kcl].Wire[kend];
      match.dWir = dw;
      if (fabs(match.dWir) > 100) return;
      if (match.Vtx >= 0) { match.dX = kX - clsChain[kpl][kcl].X[kend]; }
      else {
        match.dX = std::abs(clsChain[ipl][icl].X[iend] - clsChain[jpl][jcl].X[jend]) +
                   std::abs(clsChain[ipl][icl].X[iend] - clsChain[kpl][kcl].X[kend]);
      }
      if (prt) mf::LogVerbatim("CCTM") << " dw " << dw << " dx " << match.dX;
      // TODO: Angle matching has problems with short clusters
      if (!isShort) {
        if (ignoreSign) { match.dAng = kAng - fabs(clsChain[kpl][kcl].Angle[kend]); }
        else {
          match.dAng = kAng - clsChain[kpl][kcl].Angle[kend];
        }
      } // !isShort
      da = fabs(match.dAng) / sigmaA;
      dx = fabs(match.dX) / sigmaX;
      if (allPlnVtxMatch) {
        // matched vertex. Use angle for match error
        match.Err = vxFactor * chgAsym * da / 3;
        if (prt) mf::LogVerbatim("CCTM") << " 3-pln w Vtx  match.Err " << match.Err;
      }
      else {
        dw /= 2;
        // divide by 9
        match.Err = vxFactor * chgAsym * sqrt(dx * dx + da * da + dw * dw) / 9;
        if (prt) mf::LogVerbatim("CCTM") << " 3-pln match.Err " << match.Err;
      }
    }
    else {
      // 2-plane match
      match.dWir = -1;
      match.dAng = -1;
      match.dX = clsChain[ipl][icl].X[iend] - clsChain[jpl][jcl].X[jend];
      // degrade error by 3 for 2-plane matches
      match.Err = 3 + vxFactor * chgAsym * fabs(match.dX) / sigmaX;
      if (prt) mf::LogVerbatim("CCTM") << " 2-pln Err " << match.Err;
    } // !(nClInPln == 3)

    /////////// Matching error at the Other end
    if (nClInPln == 3) {
      // A cluster in all 3 planes
      dw = okWir - oWir[kpl];
      match.odWir = dw;
      if (match.oVtx >= 0) { match.odX = okX - oX[kpl]; }
      else {
        match.odX = std::abs(oX[ipl] - oX[jpl]) + std::abs(oX[ipl] - oX[kpl]);
      }
      if (prt)
        mf::LogVerbatim("CCTM") << " odw " << match.odWir << " odx " << match.odX << " sigmaX "
                                << sigmaX;
      // TODO: CHECK FOR SHOWER-LIKE CLUSTER OTHER END MATCH
      if (!isShort) {
        if (ignoreSign) { match.odAng = okAng - fabs(oAng[kpl]); }
        else {
          match.odAng = okAng - oAng[kpl];
        }
      } // !isShort
      da = match.odAng / sigmaA;
      dx = fabs(match.odX) / sigmaX;
      // error for wire number match
      dw /= 2;
      // divide by number of planes with clusters * 3 for dx, da and dw
      match.oErr = ovxFactor * chgAsym * sqrt(dx * dx + da * da + dw * dw) / 9;
      if (prt) mf::LogVerbatim("CCTM") << " 3-pln match.oErr " << match.oErr;
    }
    else {
      // Only 2 clusters in 3 planes
      match.odX = (oX[ipl] - oX[jpl]) / sigmaX;
      match.oErr = 3 + ovxFactor * chgAsym * fabs(match.odX);
      if (prt) mf::LogVerbatim("CCTM") << " 2-pln match.oErr " << match.oErr;
    }
  } // FillEndMatch

  ///////////////////////////////////////////////////////////////////////
  bool CCTrackMaker::FindMissingCluster(unsigned short kpl,
                                        short& kcl,
                                        unsigned short& kend,
                                        float kWir,
                                        float kX,
                                        float okWir,
                                        float okX)
  {
    // try to attach a missing cluster to the cluster chain kcl. kend is the "match end"

    unsigned short okend;
    float dxcut;

    if (kcl >= 0) return false;

    // Look for a missing cluster with loose cuts
    float kslp = fabs((okX - kX) / (okWir - kWir));
    if (kslp > 20) kslp = 20;
    // expand dX cut assuming there is a missing hit on the end of a cluster => 1 wire
    dxcut = 3 * fXMatchErr[algIndex] + kslp;
    unsigned short nfound = 0;
    unsigned short foundCl = 0, foundEnd = 0;
    for (unsigned short ccl = 0; ccl < clsChain[kpl].size(); ++ccl) {
      if (clsChain[kpl][ccl].InTrack >= 0) continue;
      // require a match at both ends
      for (unsigned short end = 0; end < 2; ++end) {
        okend = 1 - end;
        if (fabs(clsChain[kpl][ccl].Wire[end] - kWir) > 4) continue;
        if (fabs(clsChain[kpl][ccl].Wire[okend] - okWir) > 4) continue;
        // require at least one end to match
        if (fabs(clsChain[kpl][ccl].X[end] - kX) > dxcut &&
            fabs(clsChain[kpl][ccl].X[okend] - okX) > dxcut)
          continue;
        ++nfound;
        foundCl = ccl;
        foundEnd = end;
      } // end
    }   // ccl
    if (nfound == 0) return false;
    if (nfound > 1) {
      mf::LogVerbatim("CCTM") << "FindMissingCluster: Found too many matches. Write some code "
                              << nfound;
      return false;
    }
    kcl = foundCl;
    kend = foundEnd;
    return true;

  } // FindMissingCluster

  ///////////////////////////////////////////////////////////////////////
  float CCTrackMaker::ChargeAsym(std::array<float, 3>& mChg)
  {
    // find charge asymmetry between the cluster with the highest (lowest)
    // charge
    float big = 0, small = 1.E9;
    for (unsigned short ii = 0; ii < 3; ++ii) {
      if (mChg[ii] < small) small = mChg[ii];
      if (mChg[ii] > big) big = mChg[ii];
    }
    // chgAsym varies between 0 (perfect charge match) and 1 (dreadfull)
    return (big - small) / (big + small);
  } // CalculateChargeAsym

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillTrkHits(art::FindManyP<recob::Hit> const& fmCluHits,
                                 unsigned short imat,
                                 unsigned short nplanes)
  {
    // Fills the trkHits vector using cluster hits associated with the match combo imat

    unsigned short ipl;

    for (ipl = 0; ipl < 3; ++ipl)
      trkHits[ipl].clear();

    if (imat > matcomb.size()) return;

    unsigned short indx;
    std::vector<art::Ptr<recob::Hit>> clusterhits;
    unsigned short icc, ccl, icl, ecl, iht, ii;
    short endOrder, fillOrder;

    if (prt)
      mf::LogVerbatim("CCTM") << "In FillTrkHits: matcomb " << imat << " cluster chains "
                              << matcomb[imat].Cls[0] << " " << matcomb[imat].Cls[1] << " "
                              << matcomb[imat].Cls[2];

    for (ipl = 0; ipl < nplanes; ++ipl) {
      if (matcomb[imat].Cls[ipl] < 0) continue;
      // ccl is the cluster chain index
      ccl = matcomb[imat].Cls[ipl];
      // endOrder = 1 for normal order (hits added from US to DS), and -1 for reverse order
      endOrder = 1 - 2 * matcomb[imat].End[ipl];
      // re-order the sequence of cluster indices for reverse order
      if (endOrder < 0) {
        std::reverse(clsChain[ipl][ccl].ClsIndex.begin(), clsChain[ipl][ccl].ClsIndex.end());
        std::reverse(clsChain[ipl][ccl].Order.begin(), clsChain[ipl][ccl].Order.end());
        for (ii = 0; ii < clsChain[ipl][ccl].Order.size(); ++ii)
          clsChain[ipl][ccl].Order[ii] = 1 - clsChain[ipl][ccl].Order[ii];
      }
      if (ccl > clsChain[ipl].size() - 1) {
        mf::LogError("CCTM") << "Bad cluster chain index " << ccl << " in plane " << ipl;
        continue;
      }
      // loop over all the clusters in the chain
      for (icc = 0; icc < clsChain[ipl][ccl].ClsIndex.size(); ++icc) {
        icl = clsChain[ipl][ccl].ClsIndex[icc];
        if (icl > fmCluHits.size() - 1) {
          std::cout << "oops in FTH " << icl << " clsChain size " << clsChain[ipl].size() << "\n";
          exit(1);
        }
        ecl = cls[ipl][icl].EvtIndex;
        if (ecl > fmCluHits.size() - 1) {
          std::cout << "FTH bad EvtIndex " << ecl << " fmCluHits size " << fmCluHits.size() << "\n";
          continue;
        }
        clusterhits = fmCluHits.at(ecl);
        if (clusterhits.size() == 0) {
          std::cout << "FTH no cluster hits for EvtIndex " << cls[ipl][icl].EvtIndex << "\n";
          continue;
        }
        indx = trkHits[ipl].size();
        trkHits[ipl].resize(indx + clusterhits.size());
        // ensure the hit fill ordering is consistent
        fillOrder = 1 - 2 * clsChain[ipl][ccl].Order[icc];
        if (fillOrder == 1) {
          for (iht = 0; iht < clusterhits.size(); ++iht) {
            if (indx + iht > trkHits[ipl].size() - 1) std::cout << "FTH oops3\n";
            trkHits[ipl][indx + iht] = clusterhits[iht];
          }
        }
        else {
          for (ii = 0; ii < clusterhits.size(); ++ii) {
            iht = clusterhits.size() - 1 - ii;
            if (indx + ii > trkHits[ipl].size() - 1) std::cout << "FTH oops4\n";
            trkHits[ipl][indx + ii] = clusterhits[iht];
          } // ii
        }
      } // icc
      ii = trkHits[ipl].size() - 1;
      if (prt)
        mf::LogVerbatim("CCTM") << "plane " << ipl << " first p " << trkHits[ipl][0]->WireID().Plane
                                << " w " << trkHits[ipl][0]->WireID().Wire << ":"
                                << (int)trkHits[ipl][0]->PeakTime() << " last p "
                                << trkHits[ipl][ii]->WireID().Plane << " w "
                                << trkHits[ipl][ii]->WireID().Wire << ":"
                                << (int)trkHits[ipl][ii]->PeakTime();
    } // ipl

    // TODO Check the ends of trkHits to see if there are missing hits that should have been included
    // in a cluster

  } // FillTrkHits

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::PrintTracks() const
  {
    mf::LogVerbatim myprt("CCTM");
    myprt << "********* PrintTracks \n";
    myprt << "vtx  Index    X      Y      Z\n";
    for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      myprt << std::right << std::setw(4) << ivx << std::setw(4) << vtx[ivx].EvtIndex;
      myprt << std::fixed;
      myprt << std::right << std::setw(10) << std::setprecision(1) << vtx[ivx].X;
      myprt << std::right << std::setw(7) << std::setprecision(1) << vtx[ivx].Y;
      myprt << std::right << std::setw(7) << std::setprecision(1) << vtx[ivx].Z;
      if (vtx[ivx].Neutrino) myprt << " Neutrino vertex";
      myprt << "\n";
    } // ivx

    myprt << ">>>>>>>>>> Tracks \n";
    myprt << "trk  ID  Proc nht nTrj  sX     sY     sZ     eX     eY     eZ  sVx eVx sGd eGd "
             "ChgOrd  dirZ Mom PDG     ClsIndices\n";
    for (unsigned short itr = 0; itr < trk.size(); ++itr) {
      myprt << std::right << std::setw(3) << itr << std::setw(4) << trk[itr].ID;
      myprt << std::right << std::setw(5) << std::setw(4) << trk[itr].Proc;
      unsigned short nht = 0;
      for (unsigned short ii = 0; ii < 3; ++ii)
        nht += trk[itr].TrkHits[ii].size();
      myprt << std::right << std::setw(5) << nht;
      myprt << std::setw(4) << trk[itr].TrjPos.size();
      myprt << std::fixed;
      myprt << std::right << std::setw(7) << std::setprecision(1) << trk[itr].TrjPos[0](0);
      myprt << std::right << std::setw(7) << std::setprecision(1) << trk[itr].TrjPos[0](1);
      myprt << std::right << std::setw(7) << std::setprecision(1) << trk[itr].TrjPos[0](2);
      unsigned short itj = trk[itr].TrjPos.size() - 1;
      myprt << std::right << std::setw(7) << std::setprecision(1) << trk[itr].TrjPos[itj](0);
      myprt << std::right << std::setw(7) << std::setprecision(1) << trk[itr].TrjPos[itj](1);
      myprt << std::right << std::setw(7) << std::setprecision(1) << trk[itr].TrjPos[itj](2);
      myprt << std::setw(4) << trk[itr].VtxIndex[0] << std::setw(4) << trk[itr].VtxIndex[1];
      myprt << std::setw(4) << trk[itr].GoodEnd[0];
      myprt << std::setw(4) << trk[itr].GoodEnd[1];
      myprt << std::setw(4) << trk[itr].ChgOrder;
      myprt << std::right << std::setw(10) << std::setprecision(3) << trk[itr].TrjDir[itj](2);
      myprt << std::right << std::setw(4) << trk[itr].MomID;
      myprt << std::right << std::setw(5) << trk[itr].PDGCode << "   ";
      for (unsigned short ii = 0; ii < trk[itr].ClsEvtIndices.size(); ++ii)
        myprt << " " << trk[itr].ClsEvtIndices[ii];
      myprt << "\n";
    } // itr

  } // PrintTracks

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::PrintClusters(detinfo::DetectorPropertiesData const& detProp,
                                   geo::TPCID const& tpcid) const
  {

    unsigned short iTime;
    mf::LogVerbatim myprt("CCTM");
    myprt << "******* PrintClusters *********  Num_Clusters_in             Wire:Time\n";
    myprt << "vtx  Index    X       Y       Z  Pln0 Pln1 Pln2          Pln0    Pln1    Pln2\n";
    for (unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      myprt << std::right << std::setw(3) << ivx << std::setw(7) << ivx;
      myprt << std::fixed;
      myprt << std::right << std::setw(7) << std::setprecision(1) << vtx[ivx].X;
      myprt << std::right << std::setw(7) << std::setprecision(1) << vtx[ivx].Y;
      myprt << std::right << std::setw(7) << std::setprecision(1) << vtx[ivx].Z;
      myprt << std::right << std::setw(5) << vtx[ivx].nClusInPln[0];
      myprt << std::right << std::setw(5) << vtx[ivx].nClusInPln[1];
      myprt << std::right << std::setw(5) << vtx[ivx].nClusInPln[2];
      myprt << "    ";
      for (auto const& plane : wireReadoutGeom.Iterate<geo::PlaneGeo>(tpcid)) {
        int time = 0.5 + detProp.ConvertXToTicks(vtx[ivx].X, plane.ID());
        int wire = plane.WireCoordinate(geo::Point_t{0, vtx[ivx].Y, vtx[ivx].Z});
        myprt << std::right << std::setw(7) << wire << ":" << time;
      }

      myprt << "\n";
    } // ivx

    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    for (unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      myprt << ">>>>>>>>>> Cluster chains in Plane " << ipl << "\n";
      myprt << "ipl   ccl  Len    Chg    W0:T0     Ang0 Dir0  Vx0  mBk0   W1:T1     Ang1 Dir1  Vx1 "
               " mBk1  InTk    cls:Order \n";
      for (unsigned short ccl = 0; ccl < clsChain[ipl].size(); ++ccl) {
        myprt << std::right << std::setw(3) << ipl;
        myprt << std::right << std::setw(5) << ccl;
        myprt << std::right << std::setw(6) << clsChain[ipl][ccl].Length;
        myprt << std::right << std::setw(8) << (int)clsChain[ipl][ccl].TotChg;
        for (unsigned short end = 0; end < 2; ++end) {
          iTime = clsChain[ipl][ccl].Time[end];
          myprt << std::right << std::setw(5) << (int)clsChain[ipl][ccl].Wire[end] << ":"
                << std::setprecision(1) << iTime;
          if (iTime < 10) { myprt << "   "; }
          else if (iTime < 100) {
            myprt << "  ";
          }
          else if (iTime < 1000)
            myprt << " ";
          myprt << std::right << std::setw(7) << std::setprecision(2)
                << clsChain[ipl][ccl].Angle[end];
          myprt << std::right << std::setw(5) << clsChain[ipl][ccl].Dir[end];
          myprt << std::right << std::setw(5) << clsChain[ipl][ccl].VtxIndex[end];
          myprt << std::fixed << std::right << std::setw(6) << std::setprecision(1)
                << clsChain[ipl][ccl].mBrkIndex[end];
        }
        myprt << std::right << std::setw(7) << clsChain[ipl][ccl].InTrack;
        myprt << "   ";
        for (unsigned short ii = 0; ii < clsChain[ipl][ccl].ClsIndex.size(); ++ii)
          myprt << " " << clsChain[ipl][ccl].ClsIndex[ii] << ":" << clsChain[ipl][ccl].Order[ii];
        myprt << "\n";
      } // ccl
      if (fPrintAllClusters) {
        myprt << ">>>>>>>>>> Clusters in Plane " << ipl << "\n";
        myprt << "ipl  icl  Evt   Len     Chg   W0:T0     Ang0 Dir0  Vx0  CN0   W1:T1      Ang1  "
                 "Dir1  Vx1  CN1 InTk Brk0 MrgEr0 Brk1 MrgEr1\n";
        for (unsigned short icl = 0; icl < cls[ipl].size(); ++icl) {
          myprt << std::right << std::setw(3) << ipl;
          myprt << std::right << std::setw(5) << icl;
          myprt << std::right << std::setw(5) << cls[ipl][icl].EvtIndex;
          myprt << std::right << std::setw(6) << cls[ipl][icl].Length;
          myprt << std::right << std::setw(8) << (int)cls[ipl][icl].TotChg;
          for (unsigned short end = 0; end < 2; ++end) {
            iTime = cls[ipl][icl].Time[end];
            myprt << std::right << std::setw(5) << (int)cls[ipl][icl].Wire[end] << ":" << iTime;
            if (iTime < 10) { myprt << "   "; }
            else if (iTime < 100) {
              myprt << "  ";
            }
            else if (iTime < 1000)
              myprt << " ";
            myprt << std::right << std::setw(7) << std::setprecision(2) << cls[ipl][icl].Angle[end];
            myprt << std::right << std::setw(5) << cls[ipl][icl].Dir[end];
            myprt << std::right << std::setw(5) << cls[ipl][icl].VtxIndex[end];
            myprt << std::fixed << std::right << std::setw(5) << std::setprecision(1)
                  << cls[ipl][icl].ChgNear[end];
          }
          myprt << std::fixed;
          myprt << std::right << std::setw(5) << cls[ipl][icl].InTrack;
          myprt << std::right << std::setw(5) << (int)cls[ipl][icl].BrkIndex[0];
          myprt << std::right << std::setw(7) << std::setprecision(1)
                << cls[ipl][icl].MergeError[0];
          myprt << std::right << std::setw(5) << (int)cls[ipl][icl].BrkIndex[1];
          myprt << std::right << std::setw(7) << std::setprecision(1)
                << cls[ipl][icl].MergeError[1];
          myprt << "\n";
        } // icl
      }   // fPrintAllClusters
    }     // ipl
  }       // PrintClusters

  ///////////////////////////////////////////////////////////////////////
  float CCTrackMaker::AngleFactor(float slope)
  {
    float slp = fabs(slope);
    if (slp > 10.) slp = 30.;
    // return a value between 1 and 46
    return 1 + 0.05 * slp * slp;
  } // AngleFactor

  ///////////////////////////////////////////////////////////////////////
  float CCTrackMaker::ChargeNear(geo::PlaneID const& planeid,
                                 unsigned short wire1,
                                 float time1,
                                 unsigned short wire2,
                                 float time2)
  {
    // returns the hit charge along a line between (wire1, time1) and
    // (wire2, time2)

    // put in increasing wire order (wire2 > wire1)
    unsigned short w1 = wire1;
    unsigned short w2 = wire2;
    double t1 = time1;
    double t2 = time2;
    double slp, prtime;
    if (w1 == w2) { slp = 0; }
    else {
      if (w1 > w2) {
        w1 = wire2;
        w2 = wire1;
        t1 = time2;
        t2 = time1;
      }
      slp = (t2 - t1) / (w2 - w1);
    }

    unsigned short wire;

    float chg = 0;
    for (unsigned short hit = 0; hit < allhits.size(); ++hit) {
      if (allhits[hit]->WireID().asPlaneID() != planeid) continue;
      wire = allhits[hit]->WireID().Wire;
      if (wire < w1) continue;
      if (wire > w2) continue;
      prtime = t1 + (wire - w1) * slp;
      if (prtime > allhits[hit]->PeakTimePlusRMS(3)) continue;
      if (prtime < allhits[hit]->PeakTimeMinusRMS(3)) continue;
      chg += ChgNorm[planeid.Plane] * allhits[hit]->Integral();
    } // hit
    return chg;
  } // ChargeNear

  ///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillWireHitRange(geo::TPCID const& tpcid)
  {
    // fills the WireHitRange vector. Slightly modified version of the one in ClusterCrawlerAlg

    unsigned short ipl;

    // initialize everything
    for (ipl = 0; ipl < 3; ++ipl) {
      firstWire[ipl] = INT_MAX;
      lastWire[ipl] = 0;
      firstHit[ipl] = INT_MAX;
      lastHit[ipl] = 0;
      WireHitRange[ipl].clear();
      ChgNorm[ipl] = 0;
    }

    // find the first and last wire with a hit
    unsigned short oldipl = 0;
    for (unsigned int hit = 0; hit < allhits.size(); ++hit) {
      if (static_cast<geo::TPCID const&>(allhits[hit]->WireID()) != tpcid) continue;
      ipl = allhits[hit]->WireID().Plane;
      if (allhits[hit]->WireID().Wire >
          wireReadoutGeom.Nwires(allhits[hit]->WireID().asPlaneID())) {
        if (lastWire[ipl] < firstWire[ipl]) {
          mf::LogError("CCTM") << "Invalid WireID().Wire " << allhits[hit]->WireID().Wire;
          return;
        }
      }
      //      ChgNorm[ipl] += allhits[hit]->Integral();
      if (ipl < oldipl) {
        mf::LogError("CCTM") << "Hits are not in increasing-plane order\n";
        return;
      }
      oldipl = ipl;
      if (firstHit[ipl] == INT_MAX) {
        firstHit[ipl] = hit;
        firstWire[ipl] = allhits[hit]->WireID().Wire;
      }
      lastHit[ipl] = hit;
      lastWire[ipl] = allhits[hit]->WireID().Wire;
    } // hit

    // xxx
    auto const nplanes = wireReadoutGeom.Nplanes(tpcid);
    for (ipl = 0; ipl < nplanes; ++ipl) {
      if (lastWire[ipl] < firstWire[ipl]) {
        mf::LogError("CCTM") << "Invalid first/last wire in plane " << ipl;
        return;
      }
    } // ipl

    // normalize charge in induction planes to the collection plane
    //    for(ipl = 0; ipl < nplanes; ++ipl) ChgNorm[ipl] = ChgNorm[nplanes - 1] / ChgNorm[ipl];
    for (ipl = 0; ipl < nplanes; ++ipl)
      ChgNorm[ipl] = 1;

    // now we can define the WireHitRange vector.
    int sflag, nwires, wire;
    unsigned int indx, thisWire, thisHit, lastFirstHit;
    //uint32_t chan;
    for (ipl = 0; ipl < nplanes; ++ipl) {
      nwires = lastWire[ipl] - firstWire[ipl] + 1;
      WireHitRange[ipl].resize(nwires);
      // start by defining the "no hits on wire" condition
      sflag = -2;
      for (wire = 0; wire < nwires; ++wire)
        WireHitRange[ipl][wire] = std::make_pair(sflag, sflag);
      // overwrite with the "dead wires" condition
      sflag = -1;
      // next overwrite with the index of the first/last hit on each wire
      lastWire[ipl] = firstWire[ipl];
      thisHit = firstHit[ipl];
      lastFirstHit = firstHit[ipl];
      for (unsigned int hit = firstHit[ipl]; hit <= lastHit[ipl]; ++hit) {
        thisWire = allhits[hit]->WireID().Wire;
        if (thisWire > lastWire[ipl]) {
          indx = lastWire[ipl] - firstWire[ipl];
          int tmp1 = lastFirstHit;
          int tmp2 = thisHit;
          WireHitRange[ipl][indx] = std::make_pair(tmp1, tmp2);
          lastWire[ipl] = thisWire;
          lastFirstHit = thisHit;
        }
        else if (thisWire < lastWire[ipl]) {
          mf::LogError("CCTM") << "Hit not in proper order in plane " << ipl;
          exit(1);
        }
        ++thisHit;
      } // hit
      // define the last wire
      indx = lastWire[ipl] - firstWire[ipl];
      int tmp1 = lastFirstHit;
      ++lastHit[ipl];
      int tmp2 = lastHit[ipl];
      WireHitRange[ipl][indx] = std::make_pair(tmp1, tmp2);
      // add one to lastWire and lastHit for more standard indexing
      ++lastWire[ipl];
    } // ipl

    // error checking
    for (ipl = 0; ipl < nplanes; ++ipl) {
      if (firstWire[ipl] < INT_MAX) continue;
      if (lastWire[ipl] > 0) continue;
      if (firstHit[ipl] < INT_MAX) continue;
      if (lastHit[ipl] > 0) continue;
      std::cout << "FWHR problem\n";
      exit(1);
    } // ipl

    unsigned int nht = 0;
    std::vector<bool> hchk(allhits.size());
    for (unsigned int ii = 0; ii < hchk.size(); ++ii)
      hchk[ii] = false;
    for (ipl = 0; ipl < nplanes; ++ipl) {
      for (unsigned int w = firstWire[ipl]; w < lastWire[ipl]; ++w) {
        indx = w - firstWire[ipl];
        if (indx > lastWire[ipl]) {
          std::cout << "FWHR bad index " << indx << "\n";
          exit(1);
        }
        // no hit on this wire
        if (WireHitRange[ipl][indx].first < 0) continue;
        unsigned int firhit = WireHitRange[ipl][indx].first;
        unsigned int lashit = WireHitRange[ipl][indx].second;
        for (unsigned int hit = firhit; hit < lashit; ++hit) {
          ++nht;
          if (hit > allhits.size() - 1) {
            std::cout << "FWHR: Bad hchk " << hit << " size " << allhits.size() << "\n";
            continue;
          }
          hchk[hit] = true;
          if (allhits[hit]->WireID().Plane != ipl || allhits[hit]->WireID().Wire != w) {
            std::cout << "FWHR bad plane " << allhits[hit]->WireID().Plane << " " << ipl
                      << " or wire " << allhits[hit]->WireID().Wire << " " << w << "\n";
            exit(1);
          }
        } // hit
      }   // w
    }     // ipl

    if (nht != allhits.size()) {
      std::cout << "FWHR hit count problem " << nht << " " << allhits.size() << "\n";
      for (unsigned int ii = 0; ii < hchk.size(); ++ii)
        if (!hchk[ii])
          std::cout << " " << ii << " " << allhits[ii]->WireID().Plane << " "
                    << allhits[ii]->WireID().Wire << " " << (int)allhits[ii]->PeakTime() << "\n";
      exit(1);
    }

  } // FillWireHitRange
  DEFINE_ART_MODULE(CCTrackMaker)

} // namespace

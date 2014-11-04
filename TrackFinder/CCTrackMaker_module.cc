////////////////////////////////////////////////////////////////////////
//
//  CCTrackMaker
//
//  Make 3D tracks using ClusterCrawler clusters and vertex info
// 
//  baller@fnal.gov, October 2014
// 
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
//#include "art/Framework/Services/Optional/TFileService.h" 
//#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Vertex.h"
// #include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/TrackTrajectoryAlg.h"
#include "RecoAlg/LinFitAlg.h"


struct CluLen{
  int index;
  float length;
};

bool greaterThan (CluLen c1, CluLen c2) { return (c1.length > c2.length);}
bool lessThan (CluLen c1, CluLen c2) { return (c1.length < c2.length);}

namespace trkf {

  class CCTrackMaker : public art::EDProducer {
    
  public:
    
    explicit CCTrackMaker(fhicl::ParameterSet const& pset);
    virtual ~CCTrackMaker();
    
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:

    std::string     fHitModuleLabel; 
    std::string     fClusterModuleLabel; 
    std::string     fEndPoint2DModuleLabel; 
    std::string     fVertexModuleLabel;
    
    // all hits
    art::Handle< std::vector<recob::Hit> > allhitsListHandle;
    std::vector<art::Ptr<recob::Hit>> allhits;
    // cluster list
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    std::vector<art::Ptr<recob::Cluster>> clusterlist;
    // hits associated with clusters
    std::vector<art::Ptr<recob::Hit>> clusterhits;
    // ClusterCrawler EndPoints
    art::Handle< std::vector<recob::EndPoint2D> > EP2ListHandle;
    std::vector<art::Ptr<recob::EndPoint2D>> ep2list;
    // ClusterCrawler Vertices
    art::Handle< std::vector<recob::Vertex> > VtxListHandle;
    std::vector<art::Ptr<recob::Vertex>> vtxlist;

    // services
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    TrackTrajectoryAlg fTrackTrajectoryAlg;
    LinFitAlg fLinFitAlg;

    // characterize the match between clusters in 2 or 3 planes
    struct MatchPars {
      std::array<unsigned short, 3> Pln;
      std::array<short, 3> Cls;
      std::array<unsigned short, 3> End;
      std::array<float, 3> Chg;
      short Vtx;
      float dWir;   // wire difference at the matching end
      float dAng;   // angle difference at the matching end
      float dX;     // X difference
      float RMS;    // Wire,Angle,Time match rms
      short oVtx;
      float odWir;  // wire difference at the other end
      float odAng;  // angle difference at the other end
      float odX;  // time difference at the other end
      float oRMS;   // dAngle dX match rms
    };
    
    void PrintStructs();
    
    void PrintTracks();
    
    void FindBustedClusters();

    // match clusters associated with vertices
    void VtxMatch(art::Event& evt,
                art::FindManyP<recob::Hit> const& fmCluHits, 
                std::unique_ptr<std::vector<recob::Track>>& tcol,
                std::unique_ptr<art::Assns<recob::Track, recob::Hit>>& thassn);

    // match clusters in all planes
    void PlnMatch(art::Event& evt,
                art::FindManyP<recob::Hit> const& fmCluHits, 
                std::unique_ptr<std::vector<recob::Track>>& tcol,
                std::unique_ptr<art::Assns<recob::Track, recob::Hit>>& thassn);

    
    // Make the track/vertex and mother/daughter relationships
    void MakeFamily();

    // fill the end matching parameters in the MatchPars struct
    void FillEndMatch(MatchPars& match, bool prt, mf::LogVerbatim& myprt);

    bool DupMatch(MatchPars& match, std::vector<MatchPars>& matcomb);
    
    void SortMatches(std::vector<MatchPars>& matcomb,
                     art::FindManyP<recob::Hit> const& fmCluHits, 
                     unsigned short procCode,
                     mf::LogVerbatim& myprt);

    // fill the trkHits array using information
    void FillTrkHits(art::FindManyP<recob::Hit> const& fmCluHits, bool& success);

    void FillClmat(unsigned short pln, short cls, unsigned short end);

    // store the track in the trk vector
    void StoreTrack(art::FindManyP<recob::Hit> const& fmCluHits, 
      unsigned short procCode,
      bool& prt, mf::LogVerbatim& myprt);
    
    // determine the order of hits in the trkHits vector using charge
    short ChgHitOrder(bool& prt, mf::LogVerbatim& myprt);
    
    // Updates the spl array (index of plane numbers sorted by decreasing
    // number of hits in trkHits array)
    void SortPlanes(std::array< std::vector<art::Ptr<recob::Hit>>, 3>& trkHits);
    std::array<unsigned short, 3> spl;
    
    // returns the charge along the line between (wire1, time1) and (wire2, time2)
    float ChargeNear(unsigned short ipl, unsigned short wire1, double time1,
                     unsigned short wire2, double time2);
    
    // inflate cuts at large angle
    float AngleFactor(float slope);

    short fMatchAlgs;
    float fMaxDAng;
    // conversion from dT/dW to dX/dU
    float fTickToDist;
    float fWirePitch;
    float fScaleF;
    float fXMatchErr;
    float fAngleMatchErr;
    float fChgAsymFactor;
    
    short fDebugPlane;
    short fDebugCluster;
    bool prt;
    
    unsigned short nplanes;
    unsigned int cstat;
    unsigned int tpc;

    // Cluster parameters
    struct clPar{
      std::array<float, 2> Wire; // Begin/End Wire
//      std::array<float, 2> Time; // Begin/End Time
      std::array<float, 2> X; // Begin/End X
      std::array<short, 2> Time; // Begin/End Time
      std::array<float, 2> Slope; // Begin/End slope (dT/dW)
      std::array<float, 2> Charge; // Begin/End charge
      std::array<float, 2> Angle; // Begin/End angle (radians)
      std::array<float, 2> Dir; // Product of end * slope
      std::array<short, 2> EP2Index; // EndPoint2D index
      std::array<short, 2> VtxIndex; // Vertex index
      std::array<short, 2> BrkIndex; // Broken cluster index
      std::array<short, 2> mBrkIndex; // Multi broken cluster index
      short EvtIndex;   // index of the cluster in clusterlist
      short InTrack;          // cluster -> track ID (-1 if none, 0 if under construction)
      unsigned short Length;           // cluster length (wires)
      float TotChg;     // total charge of cluster (or series of clusters)
    };
    // vector of cluster parameters in each plane
    std::array<std::vector<clPar>, 3> cls;
    
    // array of cluster matches in each plane
    std::array<std::vector<short>, 3> clmat;
    // begin and end vertex indices
    std::array<short, 2> clmatVxID;
    // and the end match ordering
    std::array<unsigned short, 3> clend;
    
    // EndPoint2D parameters 
    struct EP2Par{
      float Wir;
      float X;
    };
    std::array<std::vector<EP2Par>, 3> ep2;
    
    // 3D Vertex info
    struct vtxPar{
      unsigned short ID;
      float X;
      float Y;
      float Z;
      std::array<std::vector<short>, 3> clIndex; // index to clusters in each plane (<0 = DS)
    };
    
    std::vector<vtxPar> vtx;
    
    struct TrkPar{
      unsigned short ID;
      unsigned short Proc; // 1 = VtxMatch, 2 = ...
      std::array< std::vector<art::Ptr<recob::Hit>>, 3> trkHits;
      std::vector<TVector3> trjPos;
      std::vector<TVector3> trjDir;
      std::array<short, 2> VtxIndex;
      short ChgOrder;
      short MomID;
      short nDtr;
      std::array<short, 20> DtrID;
    };
    std::vector<TrkPar> trk;
    
    // track trajectory for a track under construction
    std::vector<TVector3> trkPos;
    std::vector<TVector3> trkDir;
    
    // Array of pointers to hits in each plane for one track
    std::array< std::vector<art::Ptr<recob::Hit>>, 3> trkHits;
    

  }; // class CCTrackMaker

  //-------------------------------------------------
  CCTrackMaker::CCTrackMaker(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<recob::Track>                        >();
//    produces< art::Assns<recob::Track,      recob::Cluster>    >();
//    produces< art::Assns<recob::Track,      recob::SpacePoint> >();
//    produces< art::Assns<recob::SpacePoint, recob::Hit>        >();
    produces< art::Assns<recob::Track,      recob::Hit>        >();
  }

  //-------------------------------------------------
  void CCTrackMaker::reconfigure(fhicl::ParameterSet const& pset)
  {
    fHitModuleLabel         = pset.get< std::string >("HitModuleLabel");
    fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
    fEndPoint2DModuleLabel  = pset.get< std::string >("EndPoint2DModuleLabel");
    fVertexModuleLabel      = pset.get< std::string >("VertexModuleLabel");
    fMatchAlgs              = pset.get< short >("MatchAlgs");
    fMaxDAng                = pset.get< float >("MaxDAng");
    fXMatchErr              = pset.get< float >("XMatchErr");
    fAngleMatchErr          = pset.get< float >("AngleMatchErr");
    fChgAsymFactor          = pset.get< float >("ChgAsymFactor");
    fDebugPlane             = pset.get< short >("DebugPlane");
    fDebugCluster           = pset.get< short >("DebugCluster");
  }

  //-------------------------------------------------
  CCTrackMaker::~CCTrackMaker()
  {
  }

  //-------------------------------------------------
  void CCTrackMaker::beginJob()
  {
  }

  //-------------------------------------------------
  void CCTrackMaker::endJob()
  {
  }

  //------------------------------------------------------------------------------------//
  void CCTrackMaker::produce(art::Event& evt)
  {

    fWirePitch = geom->WirePitch();
    fTickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
    fTickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    // scale factor to convert dT/dW into angle (radians)
    fScaleF = fTickToDist / fWirePitch;


    std::unique_ptr<std::vector<recob::Track>> tcol (new std::vector<recob::Track>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit> > thassn (new art::Assns<recob::Track, recob::Hit>);

/*
    std::unique_ptr<std::vector<recob::SpacePoint> > 	        spcol(new std::vector<recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::Cluster> >    tcassn(new art::Assns<recob::Track, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> >   shassn(new art::Assns<recob::SpacePoint, recob::Hit>);
*/

    // get Hits
    if (evt.getByLabel(fHitModuleLabel, allhitsListHandle))
      art::fill_ptr_vector(allhits, allhitsListHandle);

    // get Clusters
    if (evt.getByLabel(fClusterModuleLabel, clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);
    art::FindManyP<recob::Hit> fmCluHits(clusterListHandle, evt, fClusterModuleLabel);
    if(clusterlist.size() == 0) return;

    // get ClusterCrawler EndPoints
    if (evt.getByLabel(fEndPoint2DModuleLabel, EP2ListHandle))
      art::fill_ptr_vector(ep2list, EP2ListHandle);
    // get ClusterCrawler Vertices
    if (evt.getByLabel(fVertexModuleLabel, VtxListHandle))
      art::fill_ptr_vector(vtxlist, VtxListHandle);
    
    unsigned short ipl, icl, end, iep, ivx;
    // some junk vectors to satisfy the recob::Track constructor
    std::vector< std::vector<double> > dQdx;
    std::vector<double> mom(2, util::kBogusD);
    std::vector<art::Ptr<recob::Hit> > tmpHits;    
    for(cstat = 0; cstat < geom->Ncryostats(); ++cstat) {
      for(tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc) {
        nplanes = geom->Cryostat(cstat).TPC(tpc).Nplanes();
        std::vector<double> tmp;
        double dW, dX;
        // Sort clusters by length using Tingjun's code
        std::vector<std::vector<CluLen>> clulens(nplanes);
        for(icl = 0; icl < clusterlist.size(); ++icl) {
          if(clusterlist[icl]->Plane().Cryostat != cstat) continue;
          if(clusterlist[icl]->Plane().TPC != tpc) continue;
          ipl = clusterlist[icl]->Plane().Plane;
          dW = fWirePitch * (clusterlist[icl]->EndPos()[0] - clusterlist[icl]->StartPos()[0]);
          dX = detprop->ConvertTicksToX(clusterlist[icl]->EndPos()[1], ipl, tpc, cstat) -
               detprop->ConvertTicksToX(clusterlist[icl]->StartPos()[1], ipl, tpc, cstat);
          CluLen clulen;
          clulen.index = icl;
          clulen.length = sqrt(dW * dW + dX * dX);
          clulens[ipl].push_back(clulen);
        } // icl
        
        // put the cluster information into the struct
        for(ipl = 0; ipl < nplanes; ++ipl) {
          // sort clusters in this plane
          std::sort (clulens[ipl].begin(),clulens[ipl].end(), greaterThan);
          if(clulens[ipl].size() == 0) continue;
          for(unsigned short ii = 0; ii < clulens[ipl].size(); ++ii) {
            icl = clulens[ipl][ii].index;
            clPar clstr;
            clstr.EvtIndex = icl;
            // Begin info
            tmp = clusterlist[icl]->StartPos();
            clstr.Wire[0] = tmp[0];
            clstr.Time[0] = tmp[1];
            clstr.X[0] = (float)detprop->ConvertTicksToX(tmp[1], ipl, tpc, cstat);
            clstr.Slope[0] = fScaleF * clusterlist[icl]->SigmadTdW();
            clstr.Charge[0] = clusterlist[icl]->dQdW();
            clstr.Angle[0] = atan(clstr.Slope[0]);
            clstr.Dir[0] = -1 * (2*(clstr.Slope[0]>0)-1);
            tmp = clusterlist[icl]->SigmaStartPos();
            clstr.EP2Index[0] = (short)tmp[0];
            clstr.VtxIndex[0] = (short)tmp[1];
            clstr.BrkIndex[0] = -1;
            clstr.mBrkIndex[0] = -1;
            // End info
            tmp = clusterlist[icl]->EndPos();
            clstr.Wire[1] = tmp[0];
            clstr.Time[1] = tmp[1];
            clstr.X[1] = (float)detprop->ConvertTicksToX(tmp[1], ipl, tpc, cstat);
            clstr.Slope[1] = fScaleF * clusterlist[icl]->dTdW();
            clstr.Charge[1] = clusterlist[icl]->SigmadQdW();
            clstr.Angle[1] = atan(clstr.Slope[1]);
            clstr.Dir[1] = 1 * (2*(clstr.Slope[1]>0)-1);
            tmp = clusterlist[icl]->SigmaEndPos();
            clstr.EP2Index[1] = (short)tmp[0];
            clstr.VtxIndex[1] = (short)tmp[1];
            clstr.BrkIndex[1] = -1;
            clstr.mBrkIndex[1] = -1;
            // other info
            clstr.InTrack = -1;
            clstr.Length = (unsigned short)(0.5 + clstr.Wire[0] - clstr.Wire[1]);
            clstr.TotChg = clusterlist[icl]->Charge();
            if(clstr.TotChg <= 0) clstr.TotChg = 1;
            cls[ipl].push_back(clstr);
          } // ii (icl)
        } // ipl
        // now get the endpoints
        for(iep = 0; iep < ep2list.size(); ++iep) {
          for(ipl = 0; ipl < nplanes; ++ipl) {
            if(ep2list[iep]->WireID().Plane == ipl) {
              EP2Par endpt;
              endpt.Wir = (float)ep2list[iep]->Charge();
              endpt.X = (float)detprop->ConvertTicksToX(ep2list[iep]->DriftTime(), ipl, tpc, cstat);
              ep2[ipl].push_back(endpt);
              break;
            } // ep2list[iep]->WireID().Plane == ipl
          } // ipl
        } // ivx
        // and finally the vertices
        double xyz[3];
        for(ivx = 0; ivx < vtxlist.size(); ++ivx) {
          vtxPar aVtx;
          aVtx.ID = vtxlist[ivx]->ID();
          vtxlist[ivx]->XYZ(xyz);
          aVtx.X = xyz[0];
          aVtx.Y = xyz[1];
          aVtx.Z = xyz[2];
          // make the list of pointers to clusters in each plane
          for(ipl = 0; ipl < nplanes; ++ipl) {
            if(cls[ipl].size() == 0) continue;
            for(icl = 0; icl < cls[ipl].size(); ++icl) {
              for(end = 0; end < 2; ++end) {
                if(cls[ipl][icl].VtxIndex[end] == aVtx.ID) {
                  aVtx.clIndex[ipl].push_back((short)icl);
                  // point to the index of the vertex - not the vertex ID
                  cls[ipl][icl].VtxIndex[end] = ivx;
                }
              } // end
            } // icl
          } // ipl
          vtx.push_back(aVtx);
        } // ivx
        // Find broken clusters
        FindBustedClusters();
        if(fMatchAlgs & 1) VtxMatch(evt, fmCluHits, tcol, thassn);
        if(fMatchAlgs & 2) PlnMatch(evt, fmCluHits, tcol, thassn);
        // Determine the vertex/track/PFParticle(?) hierarchy
//        MakeFamily();
        for(unsigned short itr = 0; itr < trk.size(); ++itr) {
          tcol->push_back(recob::Track(trk[itr].trjPos, trk[itr].trjDir, dQdx, mom, trk[itr].ID));
          tmpHits.clear();
          for(ipl = 0; ipl < nplanes; ++ipl) 
            tmpHits.insert(tmpHits.end(), trk[itr].trkHits[ipl].begin(), trk[itr].trkHits[ipl].end());
          util::CreateAssn(*this, evt, *tcol, tmpHits, *thassn);
        } // itr
        if(fDebugPlane >= 0) {
          PrintStructs();
          PrintTracks();
        }
        for(ipl = 0; ipl < nplanes; ++ipl) {
          cls[ipl].clear();
          ep2[ipl].clear();
        }
        vtx.clear();
        trk.clear();
      } // tpc
    } // cstat
    
    evt.put(std::move(tcol));
    evt.put(std::move(thassn));

  } // produce

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::VtxMatch(art::Event& evt,
                art::FindManyP<recob::Hit> const& fmCluHits, 
                std::unique_ptr<std::vector<recob::Track>>& tcol,
                std::unique_ptr<art::Assns<recob::Track, recob::Hit>>& thassn)
  { // Use vertex assignments to match clusters
    unsigned short ivx, ii, ipl, icl, jj, jpl, jcl, kk, kpl, kcl;
    short iend, jend, kend;

    mf::LogVerbatim myprt("CCTM");
    
    // vector of many match combinations
    std::vector<MatchPars> matcomb;
    
    bool prt = (fDebugPlane == 1);
    MatchPars match;

    for(ivx = 0; ivx < vtx.size(); ++ivx) {
    // list of clusters associated with this vertex in each plane
      std::array<std::vector<unsigned short>, 3> vxtk;
      for(ipl = 0; ipl < nplanes; ++ipl) {
        for(icl = 0; icl< cls[ipl].size(); ++icl) {
          for(iend = 0; iend < 2; ++iend) {
            if(cls[ipl][icl].VtxIndex[iend] == vtx[ivx].ID) 
                    vxtk[ipl].push_back(icl);
          } // end
        } // icl
      } // ipl
      
  if(prt) myprt<<"VtxMatch: Vertex ID "<<vtx[ivx].ID<<"\n";
      for(ipl = 0; ipl < nplanes; ++ipl) {
  if(prt) myprt<<"ipl "<<ipl<<" cls";
        for(unsigned short ii = 0; ii < vxtk[ipl].size(); ++ii) {
  if(prt) myprt<<" "<<vxtk[ipl][ii];
        } // ii
  if(prt) myprt<<"\n";
      } // ipl

      // match between planes, requiring clusters matched to this vertex
      matcomb.clear();
      iend = 0; jend = 0;
      bool gotkcl;
      for(ipl = 0; ipl < nplanes; ++ipl) {
        for(ii = 0; ii < vxtk[ipl].size(); ++ii) {
          icl = vxtk[ipl][ii];
          if(cls[ipl][icl].InTrack >= 0) continue;
          jpl = (ipl + 1) % nplanes;
          kpl = (jpl + 1) % nplanes;
          for(jj = 0; jj < vxtk[jpl].size(); ++jj) {
            jcl = vxtk[jpl][jj];
            if(cls[jpl][jcl].InTrack >= 0) continue;
            for(iend = 1; iend >= 0; --iend) {
              if(cls[ipl][icl].VtxIndex[iend] != vtx[ivx].ID) continue;
              for(jend = 1; jend >= 0; --jend) {
                if(cls[jpl][jcl].VtxIndex[jend] != vtx[ivx].ID) continue;
                if(cls[ipl][icl].Dir[iend] != cls[jpl][jcl].Dir[jend]) continue;
                match.Pln[ipl] = ipl; match.Cls[ipl] = icl; match.End[ipl] = iend;
                match.Pln[jpl] = jpl; match.Cls[jpl] = jcl; match.End[jpl] = jend;
                match.Vtx = ivx;
                gotkcl = false;
                for(kk = 0; kk < vxtk[kpl].size(); ++kk) {
                  kcl = vxtk[kpl][kk];
                  if(cls[kpl][kcl].InTrack >= 0) continue;
                  for(kend = 1; kend >= 0; --kend) {
                    if(cls[kpl][kcl].Dir[kend] != cls[ipl][icl].Dir[iend]) continue;
                    if(cls[kpl][kcl].VtxIndex[kend] != vtx[ivx].ID) continue;
                    match.Pln[kpl] = kpl; match.Cls[kpl] = kcl; match.End[kpl] = kend;
                    FillEndMatch(match, prt, myprt);
                    if(DupMatch(match, matcomb)) continue;
                    // ignore if no signal at the other end
                    if(match.Chg[kpl] <= 0) continue;
                    matcomb.push_back(match);
                    gotkcl = true;
                  } // kend
                } // kk -> kcl
                if(!gotkcl) {
                  // Try a 2 plane match if a 3 plane match didn't work
                  match.Cls[kpl] = -1; match.End[kpl] = 0;
                  FillEndMatch(match, prt, myprt);
                  matcomb.push_back(match);
                }
              } // jend
            } // iend
          } // jj
        } // ii -> icl
      } // ipl
      if(matcomb.size() == 0) continue;

      SortMatches(matcomb, fmCluHits, 1, myprt);

    } // ivx
    
  } // VtxMatch

/*
///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::VtxMatch(art::Event& evt,
                art::FindManyP<recob::Hit> const& fmCluHits, 
                std::unique_ptr<std::vector<recob::Track>>& tcol,
                std::unique_ptr<art::Assns<recob::Track, recob::Hit>>& thassn)
  { // Use vertex assignments to match clusters
    unsigned short ivx, ii, ipl, icl, jj, jpl, jcl, kk, kpl, kcl;
    bool ignoreSign;
    float best, da, kSlp, kAng;
    short kbest, kbestend, iend, jend, kend;

    mf::LogVerbatim myprt("CCTM");
    
    // vector of many match combinations
    std::vector<MatchPars> matcomb;
    
    bool prt = (fDebugPlane == 1);
    MatchPars match;

    for(ivx = 0; ivx < vtx.size(); ++ivx) {
    // list of clusters associated with this vertex in each plane
      std::array<std::vector<unsigned short>, 3> vxtk;
      for(ipl = 0; ipl < nplanes; ++ipl) {
        for(icl = 0; icl< cls[ipl].size(); ++icl) {
          for(iend = 0; iend < 2; ++iend) {
            if(cls[ipl][icl].VtxIndex[iend] == vtx[ivx].ID) 
                    vxtk[ipl].push_back(icl);
          } // end
        } // icl
      } // ipl
      
  if(prt) myprt<<"VtxMatch: Vertex ID "<<vtx[ivx].ID<<"\n";
      for(ipl = 0; ipl < nplanes; ++ipl) {
  if(prt) myprt<<"ipl "<<ipl<<" cls";
        for(unsigned short ii = 0; ii < vxtk[ipl].size(); ++ii) {
  if(prt) myprt<<" "<<vxtk[ipl][ii];
        } // ii
  if(prt) myprt<<"\n";
      } // ipl

      // match between planes, requiring clusters matched to this vertex
      // in all 3 planes
      matcomb.clear();
      iend = 0; jend = 0;
      for(ipl = 0; ipl < nplanes; ++ipl) {
        for(ii = 0; ii < vxtk[ipl].size(); ++ii) {
          icl = vxtk[ipl][ii];
          if(cls[ipl][icl].InTrack >= 0) continue;
          jpl = (ipl + 1) % nplanes;
          kpl = (jpl + 1) % nplanes;
          for(jj = 0; jj < vxtk[jpl].size(); ++jj) {
            jcl = vxtk[jpl][jj];
            if(cls[jpl][jcl].InTrack >= 0) continue;
            for(iend = 1; iend >= 0; --iend) {
              if(cls[ipl][icl].VtxIndex[iend] != vtx[ivx].ID) continue;
              for(jend = 1; jend >= 0; --jend) {
                if(cls[jpl][jcl].VtxIndex[jend] != vtx[ivx].ID) continue;
                if(cls[ipl][icl].Dir[iend] != cls[jpl][jcl].Dir[jend]) continue;
                kSlp = geom->ThirdPlaneSlope(ipl, cls[ipl][icl].Slope[iend],
                                 jpl, cls[jpl][jcl].Slope[jend], tpc, cstat);
                kAng = atan(fScaleF * kSlp);
    if(prt) myprt<<"chk i "<<ipl<<":"<<icl<<":"<<iend
      <<" "<<jpl<<":"<<jcl<<":"<<jend
      <<" iSlp "<<cls[ipl][icl].Slope[iend]<<" jSlp "<<cls[jpl][jcl].Slope[jend]
      <<" kSlp "<<kSlp<<" kAng "<<kAng<<"\n";
                // handle the case near pi/2, where the errors on large slopes 
                // could result in a wrong-sign kAng
                ignoreSign = (fabs(kSlp) > 1.5);
                if(ignoreSign) kAng = fabs(kAng);
                best = 0.4;
                kbest = -1;
                kbestend = 0;
                for(kk = 0; kk < vxtk[kpl].size(); ++kk) {
                  kcl = vxtk[kpl][kk];
                  if(cls[kpl][kcl].InTrack >= 0) continue;
                  for(kend = 1; kend >= 0; --kend) {
                    if(cls[kpl][kcl].Dir[kend] != cls[ipl][icl].Dir[iend]) continue;
                    if(cls[kpl][kcl].VtxIndex[kend] != vtx[ivx].ID) continue;
                    if(ignoreSign) {
                      da = fabs(fabs(cls[kpl][kcl].Angle[kend]) - kAng);
                    } else {
                      da = fabs(cls[kpl][kcl].Angle[kend] - kAng);
                    }
    if(prt) myprt<<" chk k "<<kpl<<":"<<kcl<<":"<<kend<<" da "<<da<<"\n";
                    if(da < best) {
                      best = da;
                      kbest = kcl;
                      kbestend = kend;
                    } // best
                  } // kend
                } // kk -> kcl
                match.Pln[ipl] = ipl; match.Cls[ipl] = icl; match.End[ipl] = iend;
                match.Pln[jpl] = jpl; match.Cls[jpl] = jcl; match.End[jpl] = jend;
                match.Vtx = ivx; 
                // the other variables are initialized in FillEndMatch
                if(kbest >= 0) {
                  match.Pln[kpl] = kpl; match.Cls[kpl] = kbest; match.End[kpl] = kbestend;
                } else {
                  match.Pln[kpl] = kpl; match.Cls[kpl] = -1; match.End[kpl] = 0;
                }
                FillEndMatch(match, prt, myprt);
                if(DupMatch(match, matcomb)) continue;
  if(prt) myprt<<"VtxMatch: chk i "<<ipl<<":"<<icl<<":"<<iend
    <<" "<<jpl<<":"<<jcl<<":"<<jend
    <<" iSlp "<<std::setprecision(2)<<cls[ipl][icl].Slope[iend]
    <<" jSlp "<<std::setprecision(2)<<cls[jpl][jcl].Slope[jend]
    <<" kWir "<<(int)kWir<<" okWir "<<(int)okWir
    <<" kSlp "<<std::setprecision(2)<<kSlp
    <<" kAng "<<std::setprecision(2)<<kAng<<"\n";
                // ignore if no signal at the other end
                if(match.Chg[kpl] <= 0) continue;
                // Drop the kpl cluster if the match is poor
                if(match.RMS + match.oRMS > 5) {
                  match.Cls[kpl] = -1; match.End[kpl] = 0;
                  FillEndMatch(match, prt, myprt);
                  if(match.RMS + match.oRMS > 10) continue;
                }
                matcomb.push_back(match);
              } // jend
            } // iend
          } // jj
        } // ii -> icl
      } // ipl
      if(matcomb.size() == 0) continue;

      SortMatches(matcomb, fmCluHits, 1, myprt);

    } // ivx
    
  } // VtxMatch
*/
///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FindBustedClusters()
  {
    // Find em and fill BrkIndex and mBrkIndex


    // re-purpose for sorting by wire number
    CluLen wnum;
    std::vector<CluLen> wireNum;
    
    unsigned short ipl, icl, icl1, icl2, end, oend;
    short dw, dx, dWCut, dw1Max, dw2Max;
    float dA, dACut = fMaxDAng;
    float dXCut, chgrat;
    // long straight clusters
    bool ls1, ls2;
    for(ipl = 0; ipl < nplanes; ++ipl) {
      if(cls[ipl].size() == 0) continue;
      for(icl1 = 0; icl1 < cls[ipl].size() - 1; ++icl1) {
        // maximum delta Wire overlap is 10% of the total length
        dw1Max = 0.1 * cls[ipl][icl1].Length;
        ls1 = (cls[ipl][icl1].Length > 100 && 
              fabs(cls[ipl][icl1].Angle[0] - cls[ipl][icl1].Angle[1]) < 0.04);
        for(icl2 = icl1 + 1; icl2 < cls[ipl].size(); ++icl2) {
          ls2 = (cls[ipl][icl2].Length > 100 && 
                fabs(cls[ipl][icl2].Angle[0] - cls[ipl][icl2].Angle[1]) < 0.04);
          dw2Max = 0.1 * cls[ipl][icl2].Length;
          // set overlap cut to be the shorter of the two
          dWCut = dw1Max;
          if(dw2Max < dWCut) dWCut = dw2Max;
          // but not exceeding 10 for very long clusters
          if(dWCut > 10) dWCut = 10;
          if(dWCut == 0) dWCut = 2;
          // check both ends
          for(end = 0; end < 2; ++end) {
            oend = 1 - end;
            // ignore if cluster ends at a vertex
            if(cls[ipl][icl1].VtxIndex[end] != -1) continue;
            if(cls[ipl][icl2].VtxIndex[oend] != -1) continue;
            // ignore if dWire exceeds the cut
            if(abs(cls[ipl][icl1].Wire[end]  
                 - cls[ipl][icl2].Wire[oend]) > dWCut) continue;
            // or if the angle exceeds the cut
            if(fabs(cls[ipl][icl1].Angle[end]) < 0.8) {
              dACut = fMaxDAng;
              dXCut = 4;
            } else {
              dACut = 2 * fMaxDAng;
              dXCut = 16;
            }
            dA = fabs(cls[ipl][icl1].Angle[end] - cls[ipl][icl2].Angle[oend]);
            // angle matching
            if(dA > dACut) continue;
            // make a rough dTime cut
            if(abs(cls[ipl][icl1].X[end] 
                 - cls[ipl][icl2].X[oend]) > 16) continue;
            // charge ratio cut
            // handle cosmic ray clusters that are broken at delta rays
            if(ls1 || ls2) {
              // tighter angle cuts but no charge cut
              if(dA > 0.04) continue;
            } else {
              chgrat = cls[ipl][icl1].Charge[end] / cls[ipl][icl2].Charge[oend];
              if(chgrat < 0.7 || chgrat > 1.5) continue;
            } // ls1 || ls2
            // project to DS end of icl2 and make a tighter time cut
            dw = fWirePitch * (cls[ipl][icl2].Wire[oend] - 
                               cls[ipl][icl1].Wire[end]);
            dx = cls[ipl][icl1].X[end] + cls[ipl][icl1].Slope[end] * dw 
               - cls[ipl][icl2].X[oend];
            if(abs(dx) > dXCut) continue;
            cls[ipl][icl1].BrkIndex[end] = icl2;
            cls[ipl][icl2].BrkIndex[oend] = icl1;
          } // end
       } // icl2
      } // icl1

      // now define the "mother" break index that points between the first and
      // last cluster. This allows easier referencing when there are more than
      // two broken clusters in a track. Sort by wire number - low to high. The
      // mother cluster is the one that is furthest DS
      wireNum.resize(cls[ipl].size());
      for(icl = 0; icl < cls[ipl].size(); ++icl) {
        wnum.index = icl;
        wnum.length = cls[ipl][icl].Wire[0];
        wireNum[icl] = wnum;
      } // icl
      std::sort (wireNum.begin(),wireNum.end(), greaterThan);
      
      short dtr;
      unsigned short nit, i1;
      std::vector<unsigned short> used;
      float chg;
      for(i1 = 0; i1 < cls[ipl].size(); ++i1) {
        icl = wireNum[i1].index;
        used.push_back(icl1);
        // index of broken cluster US of icl
        dtr = cls[ipl][icl].BrkIndex[1];
        nit = 0;
        chg = cls[ipl][icl].TotChg;
        while(dtr >= 0 && dtr < (short)cls[ipl].size() &&
              nit < cls[ipl].size()) {
          if(dtr > 0) {
            if(std::find(used.begin(), used.end(), dtr) != used.end()) break;
            used.push_back(dtr);
          } // dtr
          if(cls[ipl][dtr].BrkIndex[1] < 0) break;
          dtr = cls[ipl][dtr].BrkIndex[1];
          chg += cls[ipl][dtr].TotChg;
          ++nit;
        } // dtr >= 0 && 
        if(dtr >= 0 && dtr != cls[ipl][icl].BrkIndex[1]) {
          cls[ipl][icl].mBrkIndex[1] = dtr;
          cls[ipl][icl].TotChg = chg;
          cls[ipl][dtr].mBrkIndex[0] = icl;
          cls[ipl][dtr].TotChg = chg;
        }
      } // i1

    } // ipl
  } // FindBustedClusters

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::MakeFamily()
  {
    // define the track/vertex and mother/daughter relationships

/*
    // check proimity of tracks to vertices
    float drMax = 1, dx, dy, dr;
    std::array<float, 2> dz;
    
    unsigned short ivx, itr, end, itj;
    for(ivx = 0; ivx < vtx.size(); ++ivx) {
      for(itr = 0; itr < trk.size(); ++itr) {
        // determine which end is closer
        dz[0] = fabs(trk[itr].trjPos[0](2) - vtx[ivx].Z);
        itj = trk[itr].trjPos.size() - 1;
        dz[1] = fabs(trk[itr].trjPos[itj](2) - vtx[ivx].Z);
        // rough distance cut
        if(dz[0] > 3 && dz[1] > 3) continue;
        if(dz[0] < dz[1]) {
          end = 0;
          itj = 0;
        } else {
          end = 1;
          itj = trk[itr].trjPos.size() - 1;
        }
        dx = trk[itr].trjPos[itj](0) - vtx[ivx].X;
        dy = trk[itr].trjPos[itj](1) - vtx[ivx].Y;
        // give higher weight to X
        dr = sqrt(4*dx*dx + dy*dy + dz[end]*dz[end]) / 6;
        if(dr < drMax) {
          trk[itr].VtxIndex[end] = ivx;
          // set the trajectory point to the vertex position
          trk[itr].trjPos[itj](0) = vtx[ivx].X;
          trk[itr].trjPos[itj](1) = vtx[ivx].Y;
          trk[itr].trjPos[itj](2) = vtx[ivx].Z;
        }
      } // itr
    } // ivx

    // find the mothers
    unsigned short jtr, cnt;
    short mom;
    for(itr = 0; itr < trk.size(); ++itr) {
      if(trk[itr].ChgOrder != 0) {
        // track direction determined by charge
        if(trk[itr].VtxIndex[0] >= 0) {
          // There is a vertex at the START of this track.
          // Count the number of tracks that END at this vertex
          ivx = trk[itr].VtxIndex[0];
          cnt = 0;
          mom = -1;
          for(jtr = 0; jtr < trk.size(); ++jtr) {
            if(trk[jtr].VtxIndex[1] == ivx) {
              ++cnt;
              mom = jtr;
            }
          } // jtr
          if(cnt == 1) {
            // only one track. It must be the mother
            trk[itr].MomID = jtr;
          } // cnt == 1
        } // trk[itr].VtxIndex[0] >= 0
        // now consider the other end
      } // trk[itr].ChgOrder != 0
    } // itr
*/    
  } // MakeFamily


///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::SortPlanes(std::array< std::vector<art::Ptr<recob::Hit>>, 3>& trkHits)
  {
    // Updates the spl array (index of plane numbers sorted by decreasing
    // number of hits in trkHits array)

    // sort the planes by decreasing number of hits
    CluLen clulen;
    // index of planes sorted by decreasing number of hits
    std::vector<CluLen> tmp(3);
    unsigned short ipl;
    for(ipl = 0; ipl < nplanes; ++ipl) {
      clulen.index = ipl;
      clulen.length = trkHits[ipl].size();
      tmp[ipl] = clulen;
    }
    std::sort (tmp.begin(),tmp.end(), greaterThan);
    
    for(ipl = 0; ipl < nplanes; ++ipl) spl[ipl] = tmp[ipl].index;

  } // SortPlanes

///////////////////////////////////////////////////////////////////////
  short CCTrackMaker::ChgHitOrder(bool& prt, mf::LogVerbatim& myprt)
  {
    // returns a number (-1, 0 or 1) that represents the hit charge is
    // increasing (+1) or decreasing (-1) in the trkHits vector
    
    // charge ratios
    std::array<float, 3> chgrat;
    float chgbeg, chgend;

    SortPlanes(trkHits);
    unsigned short ip, ipl, icl, iend, bend;
    ipl = spl[1];
    if(trkHits[ipl].size() < 3) {
      mf::LogWarning("CCTM")<<"ChgHitOrder: trkHits too short "<<trkHits[ipl].size()<<"\n";
      return 0;
    }
    
    if(prt) myprt<<"ChgHitOrder \n";
    
    for(ip = 0; ip < nplanes; ++ip) {
      ipl = spl[ip];
      chgrat[ipl] = 1;
      if(clmat[ipl].size() == 0) continue;
      if(clmat[ipl][0] < 0) continue;
      icl = clmat[ipl][0];
      iend = clend[ipl];
      chgbeg = cls[ipl][icl].Charge[iend];
      // charge at the other end, allowing for broken clusters
      bend = clmat[ipl].size() - 1;
      icl = clmat[ipl][bend];
      chgend = cls[ipl][icl].Charge[1 - iend];
      chgrat[ipl] = chgbeg / chgend;
  if(prt) myprt<<"ChgHitOrder: "<<ipl<<" chgbeg "<<chgbeg<<" chgend "<<chgend
    <<" chgrat "<<chgrat[ipl]<<"\n";
    } // ip
    
    // Use charge ratio if 2 of the 3 planes have a large charge ratio
    unsigned short nbeg = 0, nend = 0;
    for(ip = 0; ip < nplanes; ++ip) {
      ipl = spl[ip];
      if(chgrat[ipl] > 1.1 ) ++nend;
      if(chgrat[ipl] < 0.9 ) ++nbeg;
    } // ipl
    
    if(nbeg > 1) return 1;
    if(nend > 1) return -1;
    
    return 0;

  } // ChgHitOrder

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::StoreTrack(art::FindManyP<recob::Hit> const& fmCluHits, 
    unsigned short procCode, bool& prt, mf::LogVerbatim& myprt)
  {
    // store the current "under construction" track in the trk vector
    
    bool success;
    FillTrkHits(fmCluHits, success);
  if(prt && !success) myprt<<"Failed in FillTrkHits\n";
    if(!success) return;
    
    TrkPar newtrk;
    
    // Track hit X and WireIDs in each plane
    std::array<std::vector<std::pair<double, geo::WireID>>,3> trkXW;
    // Track hit charge ...
    std::array<std::vector<double>,3> trkChg;
    
    newtrk.ID = trk.size() + 1;
    newtrk.Proc = procCode;
    newtrk.trkHits = trkHits;
    
    unsigned short ipl, iht;
    double xx;

    // make the track trajectory
    for(ipl = 0; ipl < 3; ++ipl) {
      trkXW[ipl].resize(trkHits[ipl].size());
      trkChg[ipl].resize(trkHits[ipl].size());
      for(iht = 0; iht < trkHits[ipl].size(); ++iht) {
        xx = detprop->ConvertTicksToX(trkHits[ipl][iht]->PeakTime(), 
                                ipl, tpc, cstat);
        trkXW[ipl][iht] = std::make_pair(xx, trkHits[ipl][iht]->WireID());
        trkChg[ipl][iht] = trkHits[ipl][iht]->Charge();
      } // iht
    } // ip
    fTrackTrajectoryAlg.TrackTrajectory(trkXW, trkPos, trkDir, trkChg);
    // failure occurred?
    if(trkPos.size() < 2 || trkPos.size() != trkDir.size()) {
      TVector3 dum;
      while(trkPos.size() < 3) {
        trkPos.push_back(dum);
        trkDir.push_back(dum);
      }
    } // bad trajectory size
    newtrk.trjPos = trkPos;
    newtrk.trjDir = trkDir;

    // Use the average end positions as the first and last trajectory point
    unsigned short icl, end, ihit, jhit, jpl, ivx;
    double X, aveX, naveX;
    double Y, Z, aveY, aveZ, naveYZ;
    unsigned short itj = 0;
    for(end = 0; end < 2; ++end) {
      if(end > 0) itj = newtrk.trjPos.size() - 1;
      aveX = 0;
      naveX = 0;
      for(ipl = 0; ipl < 3; ++ipl) {
        if(trkHits[ipl].size() == 0) continue;
        // index of the first or last hit
        if(end == 0) {
          ihit = 0;
        } else {
          ihit = trkHits[ipl].size() - 1;
        }
        X = detprop->ConvertTicksToX(trkHits[ipl][ihit]->PeakTime(), 
                                    ipl, tpc, cstat);
        aveX += X;
        ++naveX;
      } // ipl
      newtrk.trjPos[itj](0) = aveX / naveX;
      // Find the average Y,Z
      aveY = 0, aveZ = 0, naveYZ = 0;
      for(ipl = 0; ipl < 2; ++ipl) {
        if(trkHits[ipl].size() == 0) continue;
        for(jpl = ipl + 1; jpl < 3; ++jpl) {
          if(trkHits[jpl].size() == 0) continue;
          // index of the first or last hit
          if(end == 0) {
            ihit = 0;
            jhit = 0;
          } else {
            ihit = trkHits[ipl].size() - 1;
            jhit = trkHits[jpl].size() - 1;
          }
          geom->IntersectionPoint((0.5+trkHits[ipl][ihit]->WireID().Wire),
                                  (0.5+trkHits[jpl][jhit]->WireID().Wire),
                                  ipl, jpl, cstat, tpc, Y, Z);
          aveY += Y;
          aveZ += Z;
          ++naveYZ;
          if(naveYZ == 2) break;
        } // jpl
        if(naveYZ == 2) break;
      } // ipl
      if(naveYZ > 0) {
        if(end == 0) {
          newtrk.trjPos[0](1) = aveY / naveYZ;
          newtrk.trjPos[0](2) = aveZ / naveYZ;
          // correct the direction
          TVector3 tdir = newtrk.trjPos[1] - newtrk.trjPos[0];
          newtrk.trjDir[itj] = tdir.Unit();
        } else {
          newtrk.trjPos[itj](1) = aveY / naveYZ;
          newtrk.trjPos[itj](2) = aveZ / naveYZ;
          // correct the direction
          TVector3 tdir = newtrk.trjPos[itj] - newtrk.trjPos[itj-1];
          newtrk.trjDir[itj] = tdir.Unit();
        }
      } // naveYZ
    } // end

    // set vertex assignments and update the trajectory ends
    newtrk.VtxIndex[0] = clmatVxID[0];
    newtrk.VtxIndex[1] = clmatVxID[1];
    for(unsigned short end = 0; end < 2; ++end) {
      if(newtrk.VtxIndex[end] >= 0) {
        itj = 0;
        if(end == 1) itj = newtrk.trjPos.size() - 1;
        ivx = newtrk.VtxIndex[end];
        newtrk.trjPos[itj](0) = vtx[ivx].X;
        newtrk.trjPos[itj](1) = vtx[ivx].Y;
        newtrk.trjPos[itj](2) = vtx[ivx].Z;
      } //
    } // end
    newtrk.ChgOrder = ChgHitOrder(prt, myprt);
    // reverse trajectory according to charge?
    if(newtrk.ChgOrder < 0) {
      std::reverse(newtrk.trjPos.begin(), newtrk.trjPos.end());
      std::reverse(newtrk.trjDir.begin(), newtrk.trjDir.end());
      for(itj = 0; itj < newtrk.trjDir.size(); ++itj) {
        for(unsigned short ii = 0; ii < 3; ++ii) 
          newtrk.trjDir[itj](ii) = -newtrk.trjDir[itj](ii);
      } // itj
      short tmp = newtrk.VtxIndex[0];
      newtrk.VtxIndex[0] = newtrk.VtxIndex[1];
      newtrk.VtxIndex[1] = tmp;
    } // reverse
    newtrk.MomID = -1;
    newtrk.nDtr = 0;
    for(unsigned short idtr = 0; idtr < 20; ++idtr) newtrk.DtrID[idtr] = 0;
    
    // make the cluster-> track assignment
    for(ipl = 0; ipl < nplanes; ++ipl) {
      for(icl = 0; icl < cls[ipl].size(); ++icl) 
        if(cls[ipl][icl].InTrack == 0) cls[ipl][icl].InTrack = newtrk.ID;
    } // ipl
    
  if(prt) myprt<<" track ID "<<newtrk.ID<<" stored in StoreTrack \n";
    trk.push_back(newtrk);
    
  } // StoreTrack

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::PlnMatch(art::Event& evt,
                art::FindManyP<recob::Hit> const& fmCluHits, 
                std::unique_ptr<std::vector<recob::Track>>& tcol,
                std::unique_ptr<art::Assns<recob::Track, recob::Hit>>& thassn)
  {
    // Match clusters in all planes
    unsigned short ipl, icl, jpl, jcl, kpl, kcl, ivx;
    bool ignoreSign;
    float da, kSlp, kAng, kX, kWir, okWir;
    short iend, ioend, jend, joend, kend;
    
    double yp, zp;
    float tpcSizeY = geom->DetHalfWidth();
    float tpcSizeZ = geom->DetLength();
    
    float angCut = 5 * fAngleMatchErr;

    mf::LogVerbatim myprt("CCTM");
    
    MatchPars match;
    // vector of many match combinations
    std::vector<MatchPars> matcomb;
    
    bool prt = false;

    // match between planes, requiring clusters matched to this vertex
    // in all 3 planes
    matcomb.clear();
    iend = 0; jend = 0;
    for(ipl = 0; ipl < nplanes; ++ipl) {
      for(icl = 0; icl < cls[ipl].size(); ++icl) {
        if(cls[ipl][icl].InTrack >= 0) continue;
        jpl = (ipl + 1) % nplanes;
        kpl = (jpl + 1) % nplanes;
        prt = (ipl == fDebugPlane && icl == fDebugCluster);
        for(jcl = 0; jcl < cls[jpl].size(); ++jcl) {
          if(cls[jpl][jcl].InTrack >= 0) continue;
          for(iend = 1; iend >= 0; --iend) {
            // don't try to match if the cluster is broken on this end
            if(cls[ipl][icl].BrkIndex[iend] >= 0) continue;
            for(jend = 1; jend >= 0; --jend) {
              if(cls[jpl][jcl].BrkIndex[jend] >= 0) continue;
              if(cls[ipl][icl].Dir[iend] != cls[jpl][jcl].Dir[jend]) continue;
              // make a rough X cut
              if(fabs(cls[ipl][icl].X[iend] - 
                      cls[jpl][jcl].X[jend]) > 40) continue;
              ioend = 1 - iend; joend = 1 - jend;
              // Find the expected third (k) plane parameters
              kSlp = geom->ThirdPlaneSlope(ipl, cls[ipl][icl].Slope[iend],
                                           jpl, cls[jpl][jcl].Slope[jend],
                                           tpc, cstat);
              kAng = atan(kSlp);
              // Look at the match end
              if(cls[ipl][icl].VtxIndex[iend] >= 0 &&
                 cls[ipl][icl].VtxIndex[iend] == cls[jpl][jcl].VtxIndex[jend]) {
                // vertex match
                ivx = cls[ipl][icl].VtxIndex[iend];
                kX = vtx[ivx].X;
                kWir = geom->WireCoordinate(vtx[ivx].Y, vtx[ivx].Z, kpl, tpc, cstat);
              } else {
                // No Vertex. Ensure the match end is within the TPC
                geom->IntersectionPoint((0.5+cls[ipl][icl].Wire[iend]), 
                                        (0.5+cls[jpl][jcl].Wire[jend]),
                                        ipl, jpl, cstat, tpc, yp, zp);
                if(yp > tpcSizeY || yp < -tpcSizeY) continue;
                if(zp < 0 || zp > tpcSizeZ) continue;
                kX = 0.5 * (cls[ipl][icl].X[iend] + cls[jpl][jcl].X[jend]);
                kWir = geom->WireCoordinate(yp, zp, kpl, tpc, cstat);
              }
              // now look at the other end
              if(cls[ipl][icl].VtxIndex[ioend] >= 0 &&
                 cls[ipl][icl].VtxIndex[ioend] == cls[jpl][jcl].VtxIndex[joend]) {
                ivx = cls[ipl][icl].VtxIndex[ioend];
                okWir = geom->WireCoordinate(vtx[ivx].Y, vtx[ivx].Z, kpl, tpc, cstat);
              } else {
                // no vertex. Ensure the other end is within the TPC
                geom->IntersectionPoint((0.5+cls[ipl][icl].Wire[ioend]), 
                                        (0.5+cls[jpl][jcl].Wire[joend]),
                                        ipl, jpl, cstat, tpc, yp, zp);
                if(yp > tpcSizeY || yp < -tpcSizeY) continue;
                if(zp < 0 || zp > tpcSizeZ) continue;
                okWir = geom->WireCoordinate(yp, zp, kpl, tpc, cstat);
              }
  if(prt) myprt<<"PlnMatch: chk i "<<ipl<<":"<<icl<<":"<<iend
    <<" "<<jpl<<":"<<jcl<<":"<<jend
    <<" iSlp "<<std::setprecision(2)<<cls[ipl][icl].Slope[iend]
    <<" jSlp "<<std::setprecision(2)<<cls[jpl][jcl].Slope[jend]
    <<" kWir "<<(int)kWir<<" okWir "<<(int)okWir
    <<" kSlp "<<std::setprecision(2)<<kSlp
    <<" kAng "<<std::setprecision(2)<<kAng<<"\n";
              // expect a cluster in kpl?
              if(fabs(kWir - okWir) > 4) {
                // handle the case near pi/2, where the errors on large slopes 
                // could result in a wrong-sign kAng
                ignoreSign = (fabs(kSlp) > 1.5);
                if(ignoreSign) kAng = fabs(kAng);
                for(kcl = 0; kcl < cls[kpl].size(); ++kcl) {
                  if(cls[kpl][kcl].InTrack >= 0) continue;
                  for(kend = 1; kend >= 0; --kend) {
                    if(cls[kpl][kcl].BrkIndex[kend] >= 0) continue;
                    if(cls[kpl][kcl].Dir[kend] != cls[ipl][icl].Dir[iend]) continue;
                    // rough dX cut
                    if(fabs(cls[kpl][kcl].X[kend] - kX) > 10) continue;
                    // rough dWire cut
                    if(fabs(cls[kpl][kcl].Wire[kend] - kWir) > 10) continue;
                    if(ignoreSign) {
                      da = fabs(fabs(cls[kpl][kcl].Angle[kend]) - kAng);
                    } else {
                      da = fabs(cls[kpl][kcl].Angle[kend] - kAng);
                    }
  if(prt) myprt<<" chk k "<<kpl<<":"<<kcl<<":"<<kend<<" da "<<da<<"\n";
                    if(da < angCut) {
                      match.Pln[ipl] = ipl; match.Cls[ipl] = icl; match.End[ipl] = iend;
                      match.Pln[jpl] = jpl; match.Cls[jpl] = jcl; match.End[jpl] = jend;
                      match.Pln[kpl] = kpl; match.Cls[kpl] = kcl; match.End[kpl] = kend;
                      match.Chg[ipl] =   0; match.Chg[jpl] =   0; match.Chg[kpl] = 0;
                      match.Vtx = cls[ipl][icl].VtxIndex[iend]; 
                      FillEndMatch(match, prt, myprt);
  if(prt) myprt<<" dup? "<<"\n";
                      if(DupMatch(match, matcomb)) continue;
  if(prt) myprt<<" Match"
    <<" k "<<match.Pln[2]<<":"<<match.Cls[2]<<":"<<match.End[2]
    <<" oChg "<<match.Chg[kpl]<<" RMS "<<match.RMS + match.oRMS<<"\n";
                      if(match.Chg[kpl] == 0) continue;
                      // drop the kpl cluster if the match is poor
                      if(match.RMS + match.oRMS > 5) {
                        match.Cls[kpl] = -1; match.End[kpl] = 0;
                        FillEndMatch(match, prt, myprt);
                        if(match.RMS + match.oRMS > 10) continue;
                      }
                      matcomb.push_back(match);
                    } // da < angCut
                  } // kend
                } // kcl
              } else {
                // No cluster expected
                match.Pln[ipl] = ipl; match.Cls[ipl] = icl; match.End[ipl] = iend;
                match.Pln[jpl] = jpl; match.Cls[jpl] = jcl; match.End[jpl] = jend;
                match.Pln[kpl] = kpl; match.Cls[kpl] =  -1; match.End[kpl] = 0;
                if(DupMatch(match, matcomb)) continue;
                match.Chg[ipl] = 0; match.Chg[jpl] = 0; match.Chg[kpl] = 0;
                match.Vtx = cls[ipl][icl].VtxIndex[iend]; 
                FillEndMatch(match, prt, myprt);
  if(prt) myprt<<" Short"
//    <<" i "<<match.Pln[0]<<":"<<match.Cls[0]<<":"<<match.End[0]
//    <<" j "<<match.Pln[1]<<":"<<match.Cls[1]<<":"<<match.End[1]
    <<" k "<<match.Pln[2]<<":"<<match.Cls[2]<<":"<<match.End[2]
    <<" Chg "<<match.Chg[kpl]<<" RMS "<<match.RMS + match.oRMS<<"\n";
                if(match.Chg[kpl] <= 0) continue;
                if(match.RMS + match.oRMS > 10) continue;
                matcomb.push_back(match);
              }
            } // jend
          } // iend
        } // jcl
      } // icl
    } // ipl

    if(matcomb.size() == 0) return;
    
    SortMatches(matcomb, fmCluHits, 2, myprt);

  } // PlnMatch

///////////////////////////////////////////////////////////////////////
  bool CCTrackMaker::DupMatch(MatchPars& match, std::vector<MatchPars>& matcomb)
  {
    for(unsigned short imat = 0; imat < matcomb.size(); ++imat) {
      if(match.Cls[0] == matcomb[imat].Cls[0] && 
         match.Cls[1] == matcomb[imat].Cls[1] &&
         match.Cls[2] == matcomb[imat].Cls[2]) {
/*
  if((match.RMS + match.oRMS) != (matcomb[imat].RMS + matcomb[imat].oRMS )) 
    std::cout<<"different rms "<<imat<<" "<<match.RMS + match.oRMS
    <<" "<<matcomb[imat].RMS + matcomb[imat].oRMS
    <<"\n";
*/
        // compare the rms
        if((match.RMS + match.oRMS) < (matcomb[imat].RMS + matcomb[imat].oRMS )) {
          // keep the better one
          matcomb[imat].End[0] = match.End[0];
          matcomb[imat].End[1] = match.End[1];
          matcomb[imat].End[2] = match.End[2]; 
          matcomb[imat].dWir = match.dWir;
          matcomb[imat].dAng = match.dAng;
          matcomb[imat].dX = match.dX;
          matcomb[imat].RMS = match.RMS;
          matcomb[imat].oVtx = match.oVtx;
          matcomb[imat].odWir = match.odWir;
          matcomb[imat].odAng = match.odAng;
          matcomb[imat].odX = match.odX;
          matcomb[imat].oRMS = match.oRMS;
        }
        return true;
      } // test
    } // imat
    return false;
  } // DupMatch

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::SortMatches(std::vector<MatchPars>& matcomb,
           art::FindManyP<recob::Hit> const& fmCluHits, 
           unsigned short procCode, mf::LogVerbatim& myprt)
  {
    // sort by increasing total match RMS using both ends
    CluLen mrms;
    std::vector<CluLen> matrms;
    unsigned short ii, imat, ipl, icl;
    
    bool prt = false;

    for(ii = 0; ii < matcomb.size(); ++ii) {
      mrms.index = ii;
      mrms.length = matcomb[ii].RMS + matcomb[ii].oRMS;
      matrms.push_back(mrms);
    } // ii
    std::sort(matrms.begin(), matrms.end(), lessThan);
    myprt<<"SortMatches\n";
    myprt<<" ii  Vx   SRMS   RMS   dW     dA     dX oVx   oRMS  odW   odA    odX   Asym   icl   jcl   kcl \n";
    for(ii = 0; ii < matcomb.size(); ++ii) {
      imat = matrms[ii].index;
      float asym = fabs(matcomb[imat].Chg[0] - matcomb[imat].Chg[1]) / 
                       (matcomb[imat].Chg[0] + matcomb[imat].Chg[1]);
      asym *= fabs(matcomb[imat].Chg[1] - matcomb[imat].Chg[2]) / 
                  (matcomb[imat].Chg[1] + matcomb[imat].Chg[2]);
      myprt<<std::fixed<<std::right
        <<std::setw(3)<<ii
        <<std::setw(4)<<matcomb[imat].Vtx
        <<std::setw(7)<<std::setprecision(2)<<matcomb[imat].RMS + matcomb[imat].oRMS
        <<std::setw(7)<<std::setprecision(2)<<matcomb[imat].RMS
        <<std::setw(5)<<std::setprecision(1)<<matcomb[imat].dWir
        <<std::setw(7)<<std::setprecision(2)<<matcomb[imat].dAng
        <<std::setw(7)<<std::setprecision(2)<<matcomb[imat].dX
        <<std::setw(3)<<matcomb[imat].oVtx
        <<std::setw(7)<<std::setprecision(2)<<matcomb[imat].oRMS
        <<std::setw(5)<<std::setprecision(1)<<matcomb[imat].odWir
        <<std::setw(7)<<std::setprecision(2)<<matcomb[imat].odAng
        <<std::setw(7)<<std::setprecision(2)<<matcomb[imat].odX
        <<std::setw(7)<<std::setprecision(3)<<asym
        <<" "<<matcomb[imat].Pln[0]<<":"<<matcomb[imat].Cls[0]<<":"<<matcomb[imat].End[0]
        <<" "<<matcomb[imat].Pln[1]<<":"<<matcomb[imat].Cls[1]<<":"<<matcomb[imat].End[1]
        <<" "<<matcomb[imat].Pln[2]<<":"<<matcomb[imat].Cls[2]<<":"<<matcomb[imat].End[2]
        <<"\n";
    } // ii
    myprt<<"\n";
      // put the clusters for each match combo into clmat
      bool skipIt = false;
    for(ii = 0; ii < matcomb.size(); ++ii) {
      imat = matrms[ii].index;
      skipIt = false;
      for(ipl = 0; ipl < 3; ++ipl) {
        clmat[ipl].clear();
        if(matcomb[imat].Pln[ipl] < 0 || matcomb[imat].Cls[ipl] < 0) continue;
        // check for cluster used in previous matches
        icl = matcomb[imat].Cls[ipl];
        if(cls[ipl][icl].InTrack >= 0) {
//  if(prt) myprt<<" skipIt: "<<ii<<" cluster used "<<icl<<"\n";
          skipIt = true;
          break;
        } // 
        FillClmat(matcomb[imat].Pln[ipl],matcomb[imat].Cls[ipl],matcomb[imat].End[ipl]);
      } // ipl
      if(skipIt) continue;
  if(prt) {
    for(ipl = 0; ipl < 3; ++ipl) {
      myprt<<"ipl "<<ipl<<" clmat ";
      for(unsigned short iii = 0; iii < clmat[ipl].size(); ++iii) myprt<<" "<<clmat[ipl][iii];
      myprt<<"\n";
    }
  } // prt
      clmatVxID[0] = matcomb[imat].Vtx;
      clmatVxID[1] = matcomb[imat].oVtx;
      StoreTrack(fmCluHits, procCode, prt, myprt);
    } // ii
    
  } // SortMatches

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillEndMatch(MatchPars& match, bool prt, mf::LogVerbatim& myprt)
  {
    // fill the matching parameters for this cluster match. The calling routine
    // should set the match end vertex ID (if applicable) as well as the 
    // cluster IDs and matching ends in each plane. Note that the matching variables
    // Note that dWir, dAng and dTim are not filled if there is a vertex (match.Vtx >= 0).
    // Likewise, odWir, odAng and odX are not filled if there is a vertex match
    // at the other end
    
    std::array<short, 3> mVtx;
    std::array<short, 3> oVtx;
    std::array<short, 3> oWir;
    std::array<float, 3> oSlp;
    std::array<float, 3> oAng;
    std::array<float, 3> oX;
    
    std::array<float, 3> mChg;
    
    unsigned short ii, ipl, icl, iend, jpl, jcl, jend, kpl, kcl, kend, oend, ocl;
    
    for(ipl = 0; ipl < 3; ++ipl) {
      mVtx[ipl] = -1; oVtx[ipl] = -1;
      oWir[ipl] = -66; oSlp[ipl] = -66; oAng[ipl] = -66; oX[ipl] = -66;
      mChg[ipl] = -1;
    } // ipl
    
    // initialize parameters that shouldn't have been set by the calling routine
    match.dWir = 0;  match.dAng = 0;  match.dX = 0;  match.RMS = 100;
    match.odWir = 0; match.odAng = 0; match.odX = 0; match.oRMS = 100;
    match.oVtx = -1;

  if(prt) {
    myprt<<"FEM ";
    for(ipl = 0; ipl < nplanes; ++ipl) {
      myprt<<" "<<match.Pln[ipl]<<":"<<match.Cls[ipl]<<":"<<match.End[ipl];
    }
    myprt<<"\n";
  }
    
    short missingPlane = -1;
    unsigned short nClInPln = 0;
    // number of vertex matches at each end
    short aVtx = -1;
    unsigned short novxmat = 0;
    short aoVtx = -1;
    unsigned short nvxmat = 0;
    // fill the other end parameters in each plane
    for(ipl = 0; ipl < nplanes; ++ipl) {
      if(match.Cls[ipl] < 0) {
        missingPlane = ipl;
        continue;
      }
      ++nClInPln;
      icl = match.Cls[ipl];
      match.Chg[ipl] = cls[ipl][icl].TotChg;
      mChg[ipl] = cls[ipl][icl].TotChg;
      iend = match.End[ipl];
      mVtx[ipl] = cls[ipl][icl].VtxIndex[iend];
      if(mVtx[ipl] >= 0) {
        if(aVtx < 0) aVtx = mVtx[ipl];
        if(mVtx[ipl] == aVtx) ++nvxmat;
      }
  if(prt) myprt<<"FEM: m "<<ipl<<":"<<icl<<":"<<iend
    <<" Vtx "<<mVtx[ipl]<<" Wir "<<cls[ipl][icl].Wire[iend]
    <<std::fixed<<std::setprecision(1)
    <<" Slp "<<cls[ipl][icl].Slope[iend]
    <<" X "<<cls[ipl][icl].X[iend];
      oend = 1 - iend;
      ocl = icl;
      if(cls[ipl][icl].mBrkIndex[oend] >= 0) {
        // cluster at the other end of a chain of broken clusters
        ocl = cls[ipl][icl].mBrkIndex[oend];
      } else if(cls[ipl][icl].BrkIndex[oend] >= 0) {
        ocl = cls[ipl][icl].BrkIndex[oend];
      }
      oVtx[ipl] = cls[ipl][ocl].VtxIndex[oend];
      if(oVtx[ipl] >= 0) {
        if(aoVtx < 0) aoVtx = oVtx[ipl];
        if(oVtx[ipl] == aoVtx) ++novxmat;
      }
      oWir[ipl] = cls[ipl][ocl].Wire[oend];
      oAng[ipl] = cls[ipl][ocl].Angle[oend];
      oSlp[ipl] = cls[ipl][ocl].Slope[oend];
      oX[ipl] = cls[ipl][ocl].X[oend];
  if(prt) myprt<<" o "<<ipl<<":"<<ocl<<":"<<oend
    <<" oVtx "<<oVtx[ipl]<<" oWir "<<oWir[ipl]
    <<std::fixed<<std::setprecision(1)
    <<" oSlp "<<oSlp[ipl]
    <<" oX "<<oX[ipl]<<" Chg "<<(int)mChg[ipl]<<"\n";
    } // ipl
    
    if(nClInPln < 2) {
      mf::LogWarning("CCTM")<<"Not enough matched planes supplied";
      return;
    }

  if(prt) myprt<<"FEM: Vtx m "<<aVtx<<" count "<<nvxmat
    <<" o "<<aoVtx<<" count "<<novxmat
    <<" missingPlane "<<missingPlane
    <<" nClInPln "<<nClInPln<<"\n";

    // perfect match
    if(nvxmat == 3 && novxmat == 3) {
      match.Vtx  =  aVtx; match.RMS  = 0;
      match.oVtx = aoVtx; match.oRMS = 0;
      return;
    }
    
    // 2-plane vertex match?
    // factors applied to RMS = 1 (no vtx), 0.5 (2 pln vtx), 0.33 (3 pln vtx)
    float vxFactor = 1;
    float ovxFactor = 1;
    if(nClInPln == 3) {
      // a cluster in all 3 planes
      if(nvxmat == 3) {
        // and all vertex assignments agree at the match end
        match.Vtx  =  aVtx; vxFactor = 0.33;
      }
      if(novxmat == 3) {
        // and all vertex assignments agree at the other end
        match.oVtx  =  aoVtx; ovxFactor = 0.33;
      }
    } else {
      // a cluster in 2 planes
      if(nvxmat == 2) {
        match.Vtx  =  aVtx; vxFactor = 0.5;
      }
      if(novxmat == 2) {
        match.oVtx  =  aoVtx; ovxFactor = 0.5;
      }
    } // nClInPln

    // find wire, X and Time at both ends

    // Find the "other end" matching parameters with
    // two cases: a 3-plane match or a 2-plane match
    // and with/without an other end vertex
    
    double ypos, zpos;
    unsigned short kWir, okWir;
    float kSlp, okSlp, kAng, okAng, okX, kX, kTim, okTim;

    if(nClInPln == 3) {
      ipl = 0; jpl = 1; kpl = 2;
    } else {
      // 2-plane match 
      kpl = missingPlane;
      if(kpl == 0) {
        ipl = 1; jpl = 2;
      } else if(kpl == 1) {
        ipl = 2; jpl = 0;
      } else {
        ipl = 0; jpl = 1;
      } // kpl test
    } // missing plane
    iend = match.End[ipl]; jend = match.End[jpl]; kend = match.End[kpl];
    icl  = match.Cls[ipl];  jcl = match.Cls[jpl];  kcl = match.Cls[kpl];
    
/////////// Wire, Angle, X and Time at the Other end
    okSlp = geom->ThirdPlaneSlope(ipl, oSlp[ipl], jpl, oSlp[jpl], tpc, cstat);
    okAng = atan(okSlp);
    // handle the case near pi/2, where the errors on large slopes could result in 
    // a wrong-sign kAng
    bool ignoreSign = (fabs(okSlp) > 10);
    if(ignoreSign) okAng = fabs(okAng);
    if(match.oVtx >= 0) {
      // a vertex exists at the other end
      okWir = geom->WireCoordinate(vtx[match.oVtx].Y, vtx[match.oVtx].Z, kpl, tpc, cstat);
      okX = vtx[match.oVtx].X;
    } else {
      // no vertex at the other end
      geom->IntersectionPoint(oWir[ipl], oWir[jpl],ipl, jpl, cstat, tpc, ypos, zpos);
      okWir = (0.5 + geom->WireCoordinate(ypos, zpos, kpl, tpc, cstat));
      okX = 0.5 * (oX[ipl] + oX[jpl]);
    }
    okTim = detprop->ConvertXToTicks(okX, kpl, tpc, cstat);
  if(prt) myprt<<"FEM: oEnd"
    <<" kpl "<<kpl<<" okSlp "<<okSlp
    <<" okAng "<<okAng<<" okWir "<<okWir<<" okX "<<okX<<"\n";

/////////// Wire, Angle, X and Time at the Match end
    kSlp = geom->ThirdPlaneSlope(ipl, cls[ipl][icl].Slope[iend],
                                 jpl, cls[jpl][jcl].Slope[jend], tpc, cstat);
    kAng = atan(kSlp);
    if(ignoreSign) kAng = fabs(kAng);
    if(match.Vtx >= 0) {
      // a vertex exists at the match end
      kWir = geom->WireCoordinate(vtx[match.Vtx].Y, vtx[match.Vtx].Z, kpl, tpc, cstat);
      kX = vtx[match.Vtx].X;
    } else {
      // no vertex at the match end
      geom->IntersectionPoint(cls[ipl][icl].Wire[iend],
                              cls[jpl][jcl].Wire[jend],
                              ipl, jpl, cstat, tpc, ypos, zpos);
      kWir = (0.5 + geom->WireCoordinate(ypos, zpos, kpl, tpc, cstat));
      kX = 0.5 * (cls[ipl][icl].X[iend] + cls[jpl][jcl].X[jend]);
    }
    kTim = detprop->ConvertXToTicks(kX, kpl, tpc, cstat);
  if(prt) myprt<<"FEM: mEnd"
    <<" kpl "<<kpl<<" kSlp "<<kSlp
    <<" kAng "<<kAng<<" kX "<<kX<<"\n";

    // Calculate the cluster charge asymmetry. This factor will be applied
    // to the rms of the end matches
    float chgAsym = 1;
    // Get the charge in the plane without a matching cluster
    if(nClInPln < 3 && mChg[missingPlane] <= 0) {
      if(missingPlane != kpl) std::cout<<"Whoops "<<missingPlane<<" "<<kpl<<"\n";
      mChg[kpl] = ChargeNear(kpl, kWir, kTim, okWir, okTim);
      match.Chg[kpl] = mChg[kpl];
  if(prt) myprt<<"FEM: ChargeNear "<<kWir<<":"<<(int)kTim
  <<" "<<okWir<<":"<<(int)okTim<<" chg "<<mChg[kpl]<<"\n";
      if(mChg[kpl] <= 0) return;
    }

    // asymmetry factor varies between 1 and 2
    chgAsym = fabs(mChg[0] - mChg[1]) / (mChg[0] + mChg[1]);
    chgAsym *= fabs(mChg[1] - mChg[2]) / (mChg[1] + mChg[2]);
    chgAsym = 1 + fChgAsymFactor * chgAsym;

  if(prt) myprt<<"FEM: charge asymmetry factor "<<chgAsym<<"\n";

    float sigmaX = fXMatchErr * AngleFactor(okSlp);
    float sigmaA = fAngleMatchErr * AngleFactor(okSlp);
    float da, dx, dw;
    
/////////// Matching rms at the Other end
    if(nClInPln == 3) {
      // A cluster in all 3 planes
      dw = okWir - oWir[kpl];
      match.odWir = fabs(dw);
  if(prt && match.odWir > 100) myprt<<"Bad odWir: okWir = "<<okWir
    <<" oWir"<<oWir[kpl]<<"\n";
      if(match.odWir > 100) return;
      match.odX = fabs(okX - oX[kpl]);
  if(prt) myprt<<" odw "<<match.odWir<<" odx "<<match.odX<<"\n";
      if(ignoreSign) {
        match.odAng = fabs(okAng - fabs(oAng[kpl]));
      } else {
        match.odAng = fabs(okAng - oAng[kpl]);
      }
      da = match.odAng / sigmaA;
      dx = match.odX / sigmaX;
      // error for wire number match
      dw /= 2;
      // divide by number of planes with clusters * 3 for dx, da and dw
      match.oRMS = ovxFactor * chgAsym * sqrt(dx*dx + da*da + dw*dw) / 9;
  if(prt) myprt<<" 3-pln oRMS "<<match.oRMS<<"\n";
    } else {
      // Only 2 clusters in 3 planes
      match.oRMS = ovxFactor * chgAsym * 3;
  if(prt) myprt<<" 2-pln oRMS "<<match.oRMS<<"\n";
    }
    
/////////// Matching rms at the Match end
    // check for vertex consistency at the match end
    aVtx = -1;
    bool allPlnVtxMatch = false;
    if(nClInPln == 3) {
      unsigned short nmvtx = 0;
      for(ii = 0; ii < nplanes; ++ii) {
        if(mVtx[ii] >= 0) {
          if(aVtx < 0) aVtx = mVtx[ii];
          ++nmvtx;
        }
      } // ii
      // same vertex in all planes
      if(nmvtx ) allPlnVtxMatch = true;
    } // nClInPln
    sigmaX = fXMatchErr * AngleFactor(kSlp);
    sigmaA = fAngleMatchErr * AngleFactor(kSlp);
    
    if(nClInPln == 3) {
      kcl = match.Cls[kpl];
      kend = match.End[kpl];
      dw = kWir - cls[kpl][kcl].Wire[kend];
      match.dWir = fabs(dw);
  if(prt && match.dWir > 100) myprt<<"Bad dWir: kWir = "<<kWir
    <<" cls Wir"<<cls[kpl][kcl].Wire[kend]<<"\n";
      if(match.dWir > 100) return;
      match.dX = fabs(kX - cls[kpl][kcl].X[kend]);
  if(prt) myprt<<" dw "<<dw<<" dx "<<match.dX<<"\n";
      if(ignoreSign) {
        match.dAng = fabs(kAng - fabs(cls[kpl][kcl].Angle[kend]));
      } else {
        match.dAng = fabs(kAng - cls[kpl][kcl].Angle[kend]);
      }
      da = match.dAng / sigmaA;
      dx = match.dX / sigmaX;
      if(allPlnVtxMatch) {
        // matched vertex. Use angle for match rms
        match.RMS = vxFactor * chgAsym * da / 3;
  if(prt) myprt<<" 3-pln w Vtx  RMS "<<match.RMS<<"\n";
      } else {
        dw /= 2;
        // divide by 9
        match.RMS = vxFactor * chgAsym * sqrt(dx*dx + da*da + dw*dw) / 9;
  if(prt) myprt<<" 3-pln RMS "<<match.RMS<<"\n";
      }
    } else {
      // 2-plane match
      match.dWir = -1;
      match.dAng = -1;
      match.dX = fabs(cls[ipl][icl].X[iend] - cls[jpl][jcl].X[jend]);
      match.RMS = vxFactor * chgAsym * match.dX / (2 * sigmaX);
  if(prt) myprt<<" 2-pln RMS "<<match.RMS<<"\n";
    } // !(nClInPln == 3)

  } // FillEndMatch

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillClmat(unsigned short pln, short cind, unsigned short end)
  {
  
    short oend = 1 - end;
    short dtr;
    
    clmat[pln].clear();
    
    if(cind < 0) return;
    
    clend[pln] = end;
    clmat[pln].push_back(cind);
    
    unsigned short icl = cind;

    // simple unbroken cluster
    if(cls[pln][icl].BrkIndex[oend] < 0) return;

    if(cls[pln][icl].mBrkIndex[oend] >= 0) {
      // multi-broken cluster. Follow the chain of breaks
      dtr = cls[pln][icl].BrkIndex[oend];
      unsigned short nit = 0;
      while(dtr >= 0 && (unsigned short)dtr < cls[pln].size()) {
        clmat[pln].push_back(dtr);
        if(cls[pln][dtr].BrkIndex[oend] < 0) break;
        dtr = cls[pln][dtr].BrkIndex[oend];
        ++nit;
        if(nit > 20) return;
      } // dtr >= 0
    } else {
      // a simple broken cluster
      dtr = cls[pln][icl].BrkIndex[oend];
      clmat[pln].push_back(dtr);
    }

    // reverse the cluster index order?
    if(end == 1) std::reverse(clmat[pln].begin(),clmat[pln].end());
    
  } // FillClmat

///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::FillTrkHits(art::FindManyP<recob::Hit> const& fmCluHits,
    bool& success)
  {
    // put hits associated with clusters identified by the clmat array
    unsigned short ipl, icl, im, iht;
    
    success = false;
    
    for(ipl = 0; ipl < 3; ++ipl) trkHits[ipl].clear();

    for(ipl = 0; ipl < 3; ++ipl) {
      if(clmat[ipl].size() == 0) continue;
      // enter hits in the proper order for each cluster in clmat
      if(clend[ipl] == 0) {
        for(im = 0; im < clmat[ipl].size(); ++im) {
          if(clmat[ipl][im] < 0) continue;
          icl = clmat[ipl][im];
          if(cls[ipl][icl].InTrack != -1) {
            mf::LogWarning("CCTM")<<"FillTrkHits: trying to re-use cluster "
              <<icl<<" in plane "<<ipl;
            return;
          } // sanity check
          cls[ipl][icl].InTrack = 0;
          clusterhits = fmCluHits.at(cls[ipl][icl].EvtIndex);
          for(iht = 0; iht < clusterhits.size(); ++iht) {
            trkHits[ipl].push_back(clusterhits[iht]);
          } // iht
        } // imat
      } else {
        for(im = 0; im < clmat[ipl].size(); ++im) {
          if(clmat[ipl][im] < 0) continue;
          icl = clmat[ipl][clmat[ipl].size()-1-im];
          if(cls[ipl][icl].InTrack != -1) {
            mf::LogWarning("CCTM")<<"FillTrkHits: trying to re-use cluster "
              <<icl<<" in plane "<<ipl;
            return;
          } // sanity check
          cls[ipl][icl].InTrack = 0;
          clusterhits = fmCluHits.at(cls[ipl][icl].EvtIndex);
          for(iht = 0; iht < clusterhits.size(); ++iht) {
            trkHits[ipl].push_back(clusterhits[clusterhits.size()-1-iht]);
          } // iht
        } // imat
      }
    } // ii
    
    success = true;
    
  } // FillTrkHits
  
///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::PrintTracks()
  {
    mf::LogVerbatim myprt("CCTrackmaker");
    myprt<<"********* PrintTracks \n";
    myprt<<"vtx  ID    X      Y      Z\n";
    for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      myprt<<std::right<<std::setw(3)<<ivx<<std::setw(4)<<vtx[ivx].ID;
      myprt<<std::fixed;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<vtx[ivx].X;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<vtx[ivx].Y;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<vtx[ivx].Z;
      myprt<<"\n";
    } // ivx
    
    myprt<<">>>>>>>>>> Tracks \n";
    myprt<<"trk  ID  Proc nht nTrj  sX     sY     sZ     eX     eY     eZ  sVx eVx ChgOrder dirZ Mom\n";
    for(unsigned short itr = 0; itr < trk.size(); ++itr) {
      myprt<<std::right<<std::setw(3)<<itr<<std::setw(4)<<trk[itr].ID;
      myprt<<std::right<<std::setw(5)<<std::setw(4)<<trk[itr].Proc;
      unsigned short nht = 0;
      for(unsigned short ii = 0; ii < 3; ++ii) nht += trk[itr].trkHits[ii].size();
      myprt<<std::right<<std::setw(5)<<nht;
      myprt<<std::setw(4)<<trk[itr].trjPos.size();
      myprt<<std::fixed;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<trk[itr].trjPos[0](0);
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<trk[itr].trjPos[0](1);
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<trk[itr].trjPos[0](2);
      unsigned short itj = trk[itr].trjPos.size() - 1;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<trk[itr].trjPos[itj](0);
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<trk[itr].trjPos[itj](1);
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<trk[itr].trjPos[itj](2);
      myprt<<std::setw(4)<<trk[itr].VtxIndex[0]<<std::setw(4)<<trk[itr].VtxIndex[1];
      myprt<<std::setw(4)<<trk[itr].ChgOrder;
      myprt<<std::right<<std::setw(10)<<std::setprecision(3)<<trk[itr].trjDir[itj](2);
      myprt<<std::right<<std::setw(4)<<trk[itr].MomID;
      myprt<<"\n";
    } // itr
    
  } // PrintTracks


///////////////////////////////////////////////////////////////////////
  void CCTrackMaker::PrintStructs()
  {

    mf::LogVerbatim myprt("CCTrackmaker");
    myprt<<"********* PrintStructs nplanes = "<<nplanes<<"\n";
    myprt<<"vtx  ID    X       Y       Z   clusters \n";
    for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      myprt<<std::right<<std::setw(3)<<ivx<<std::setw(4)<<vtx[ivx].ID;
      myprt<<std::fixed;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<vtx[ivx].X;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<vtx[ivx].Y;
      myprt<<std::right<<std::setw(7)<<std::setprecision(1)<<vtx[ivx].Z;
      for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
        myprt<<" pln"<<ipl;
        for(unsigned short icl = 0; icl < vtx[ivx].clIndex[ipl].size(); ++icl) {
          myprt<<" "<<vtx[ivx].clIndex[ipl][icl];
        }
      } // ipl
      myprt<<"\n";
    } // ivx
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
/*
      myprt<<">>>>>>>>>> EndPoints in Plane "<<ipl<<"\n";
      myprt<<"iep    W   X \n";
      for(unsigned short iep = 0; iep < ep2[ipl].size(); ++iep) {
        myprt<<std::right<<std::setw(3)<<iep;
        myprt<<std::right<<std::setw(5)<<(int)ep2[ipl][iep].Wir;
        myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(1)<<ep2[ipl][iep].X<<"\n";
      } // iep
*/
      myprt<<">>>>>>>>>> Clusters in Plane "<<ipl<<"\n";
      myprt<<"ipl icl ID  wLen  BBk  EBk BmBk EmBK  BVx  EVx InTk     Chg   BW:T     BAng   EW:T     EAng\n";
      for(unsigned short icl = 0; icl < cls[ipl].size(); ++icl) {
        myprt<<std::right<<std::setw(3)<<ipl;
        myprt<<std::right<<std::setw(4)<<icl;
        myprt<<std::right<<std::setw(4)<<cls[ipl][icl].EvtIndex;
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].Length;
        myprt<<std::fixed;
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].BrkIndex[0];
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].BrkIndex[1];
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].mBrkIndex[0];
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].mBrkIndex[1];
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].VtxIndex[0];
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].VtxIndex[1];
        myprt<<std::right<<std::setw(5)<<cls[ipl][icl].InTrack;
        myprt<<std::right<<std::setw(8)<<(int)cls[ipl][icl].TotChg;
        for(unsigned short end = 0; end < 2; ++end) {
          myprt<<std::right<<std::setw(5)<<(int)cls[ipl][icl].Wire[end]
            <<":"<<std::setprecision(1)<<cls[ipl][icl].Time[end];
          myprt<<std::right<<std::setw(6)<<std::setprecision(2)<<cls[ipl][icl].Angle[end];
//          myprt<<std::right<<std::setw(5)<<std::setprecision(0)<<cls[ipl][icl].Charge[end];
        }
        myprt<<"\n";
      } // icl
    } // ipl
  } // PrintStructs

///////////////////////////////////////////////////////////////////////
  float CCTrackMaker::AngleFactor(float slope)
  {
    float slp = fabs(slope);
    if(slp > 10.) slp = 10.;
    // return a value between 1 and 6
    return 1 + 0.05 * slp * slp;
  } // AngleFactor
  
///////////////////////////////////////////////////////////////////////
  float CCTrackMaker::ChargeNear(unsigned short ipl, unsigned short wire1, double time1,
          unsigned short wire2, double time2)
  {
    // returns the hit charge along a line between (wire1, time1) and
    // (wire2, time2)
    
    // put in increasing wire order (wire2 > wire1)
    if(wire1 > wire2) {
      unsigned short itmp = wire1;
      wire1 = wire2;
      wire2 = itmp;
      double tmp = time1;
      time1 = time2;
      time2 = tmp;
    } // wire
    
    double slp = 0, prtime;
    if(wire2 > wire1) slp = (time2 - time1) / (wire2 - wire1);
    unsigned short wire;

    float chg = 0;
    for(unsigned short hit = 0; hit < allhits.size(); ++hit) {
      if(allhits[hit]->WireID().Cryostat != cstat) continue;
      if(allhits[hit]->WireID().TPC != tpc) continue;
      if(allhits[hit]->WireID().Plane != ipl) continue;
      wire = allhits[hit]->WireID().Wire;
      if(wire < wire1) continue;
      if(wire > wire2) continue;
      prtime = time1 + (wire - wire1) * slp;
      if(prtime > allhits[hit]->EndTime() + 20) continue;
      if(prtime < allhits[hit]->StartTime() - 20) continue;
      chg += allhits[hit]->Charge();
    } // hit
    return chg;
  } // ChargeNear
  
  
  DEFINE_ART_MODULE(CCTrackMaker)

} // namespace

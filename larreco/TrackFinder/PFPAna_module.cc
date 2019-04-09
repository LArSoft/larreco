////////////////////////////////////////////////////////////////////////
//
// PFPAna class
//
// Bruce Baller
//
////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <TH1F.h>
#include <TProfile.h>
#include <vector>
#include <string>
#include <array>

//Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/Utilities/AssociationUtil.h"


#include "art/Framework/Core/EDAnalyzer.h"



class TH1F;
class TH2F;

namespace pfpf {

   
  class PFPAna : public art::EDAnalyzer {

  public:
          
    explicit PFPAna(fhicl::ParameterSet const& pset); 
    virtual ~PFPAna();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

  private:
    TH1F* fNClusters;
    TH1F* fNHitInCluster;
    // Cosmic Rays
    TH1F* fCREP2;
    TH1F* fCRE;
    TH1F* fCRP;
    // Neutrino interactions
    TH1F* fNuKE_elec;
    TH1F* fNuKE_muon;
    TH1F* fNuKE_pion;
    TH1F* fNuKE_kaon;
    TH1F* fNuKE_prot;
  //  TH1F* fNuEP2;
    TH1F* fNuEP2_elec;
    TH1F* fNuEP2_muon;
    TH1F* fNuEP2_pion;
    TH1F* fNuEP2_kaon;
    TH1F* fNuEP2_prot;
    TH1F* fNuE_elec;
    TH1F* fNuE_muon;
    TH1F* fNuE_pion;
    TH1F* fNuE_kaon;
    TH1F* fNuE_prot;
    TH1F* fNuP_elec;
    TH1F* fNuP_muon;
    TH1F* fNuP_pion;
    TH1F* fNuP_kaon;
    TH1F* fNuP_prot;
    
    TH1F* fNuVtx_dx;
    TH1F* fNuVtx_dy;
    TH1F* fNuVtx_dz;
    
    TProfile* fNuEP2_KE_elec;
    TProfile* fNuEP2_KE_muon;
    TProfile* fNuEP2_KE_pion;
    TProfile* fNuEP2_KE_kaon;
    TProfile* fNuEP2_KE_prot;
    
    std::string fHitsModuleLabel;
    std::string fClusterModuleLabel;
    std::string fTrackModuleLabel;
    std::string fPFParticleModuleLabel;
    std::string fVertexModuleLabel;
    std::vector<float> fElecKERange;
    std::vector<float> fMuonKERange;
    std::vector<float> fPionKERange;
    std::vector<float> fKaonKERange;
    std::vector<float> fProtKERange;
    short fTrackWeightOption;
    bool  fMergeDaughters;
  //  float fMergeAngleCut;
    bool  fSkipCosmics;
    short fPrintLevel;
      	 
  }; // class PFPAna


  //--------------------------------------------------------------------
  PFPAna::PFPAna(fhicl::ParameterSet const& pset)  
    : EDAnalyzer(pset)
    , fHitsModuleLabel      (pset.get< std::string > ("HitsModuleLabel"))
    , fClusterModuleLabel   (pset.get< std::string > ("ClusterModuleLabel"))
    , fTrackModuleLabel     (pset.get< std::string > ("TrackModuleLabel"))
    , fPFParticleModuleLabel   (pset.get< std::string > ("PFParticleModuleLabel"))
    , fVertexModuleLabel    (pset.get< std::string > ("VertexModuleLabel"))
    , fElecKERange    (pset.get< std::vector<float>> ("ElecKERange"))
    , fMuonKERange    (pset.get< std::vector<float>> ("MuonKERange"))
    , fPionKERange    (pset.get< std::vector<float>> ("PionKERange"))
    , fKaonKERange    (pset.get< std::vector<float>> ("KaonKERange"))
    , fProtKERange    (pset.get< std::vector<float>> ("ProtKERange"))
    , fTrackWeightOption    (pset.get< short >       ("TrackWeightOption"))
    , fMergeDaughters       (pset.get< bool >        ("MergeDaughters"))
    , fSkipCosmics          (pset.get< bool >        ("SkipCosmics"))
    , fPrintLevel           (pset.get< short >       ("PrintLevel"))
  {
  
  
  }
  
  //------------------------------------------------------------------
  PFPAna::~PFPAna()
  {
  
  }
  
  void PFPAna::beginJob()
  {
  
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
  
    fNClusters=tfs->make<TH1F>("fNoClustersInEvent","Number of Clusters", 40,0 ,400);
    fNHitInCluster = tfs->make<TH1F>("fNHitInCluster","NHitInCluster",100,0,100);

    if(!fSkipCosmics) {
      // Cosmic ray histos
      fCREP2 = tfs->make<TH1F>("CREP2","CREP2",50,0,1);
      fCRE   = tfs->make<TH1F>("CRE","CR Efficiency",50,0,1);
      fCRP   = tfs->make<TH1F>("CRP","CR Purity",50,0,1);
    }
    // Neutrino Int histos
    fNuKE_elec = tfs->make<TH1F>("NuKE_elec","NuKE electron",100,0,4000);
    fNuKE_muon = tfs->make<TH1F>("NuKE_muon","NuKE muon",100,0,4000);
    fNuKE_pion = tfs->make<TH1F>("NuKE_pion","NuKE pion",100,0,4000);
    fNuKE_kaon = tfs->make<TH1F>("NuKE_kaon","NuKE kaon",100,0,4000);
    fNuKE_prot = tfs->make<TH1F>("NuKE_prot","NuKE proton",100,0,4000);
    
    fNuEP2_elec = tfs->make<TH1F>("NuEP2_elec","NuEP2 electron",50,0,1);
    fNuEP2_muon = tfs->make<TH1F>("NuEP2_muon","NuEP2 muon",50,0,1);
    fNuEP2_pion = tfs->make<TH1F>("NuEP2_pion","NuEP2 pion",50,0,1);
    fNuEP2_kaon = tfs->make<TH1F>("NuEP2_kaon","NuEP2 kaon",50,0,1);
    fNuEP2_prot = tfs->make<TH1F>("NuEP2_prot","NuEP2 proton",50,0,1);

    fNuE_elec = tfs->make<TH1F>("NuE_elec","Nu Efficiency electron",50,0,1);
    fNuE_muon = tfs->make<TH1F>("NuE_muon","Nu Efficiency muon",50,0,1);
    fNuE_pion = tfs->make<TH1F>("NuE_pion","Nu Efficiency pion",50,0,1);
    fNuE_kaon = tfs->make<TH1F>("NuE_kaon","Nu Efficiency kaon",50,0,1);
    fNuE_prot = tfs->make<TH1F>("NuE_prot","Nu Efficiency proton",50,0,1);

    fNuP_elec = tfs->make<TH1F>("NuP_elec","Nu Purity electron",50,0,1);
    fNuP_muon = tfs->make<TH1F>("NuP_muon","Nu Purity muon",50,0,1);
    fNuP_pion = tfs->make<TH1F>("NuP_pion","Nu Purity pion",50,0,1);
    fNuP_kaon = tfs->make<TH1F>("NuP_kaon","Nu Purity kaon",50,0,1);
    fNuP_prot = tfs->make<TH1F>("NuP_prot","Nu Purity proton",50,0,1);

    // True - Reco vertex difference
    fNuVtx_dx = tfs->make<TH1F>("Vtx dx","Vtx dx",80,-10,10);
    fNuVtx_dy = tfs->make<TH1F>("Vtx dy","Vtx dy",80,-10,10);
    fNuVtx_dz = tfs->make<TH1F>("Vtx dz","Vtx dz",80,-10,10);
    
    fNuEP2_KE_elec = tfs->make<TProfile>("NuEP2_KE_elec","NuEP2 electron vs KE",20,0,2000);
    fNuEP2_KE_muon = tfs->make<TProfile>("NuEP2_KE_muon","NuEP2 muon vs KE",20,0,2000);
    fNuEP2_KE_pion = tfs->make<TProfile>("NuEP2_KE_pion","NuEP2 pion vs KE",20,0,2000);
    fNuEP2_KE_kaon = tfs->make<TProfile>("NuEP2_KE_kaon","NuEP2 kaon vs KE",20,0,2000);
    fNuEP2_KE_prot = tfs->make<TProfile>("NuEP2_KE_prot","NuEP2 proton vs KE",20,0,2000);
  
  }
  
  void PFPAna::analyze(const art::Event& evt)
  {
    
    // code stolen from TrackAna_module.cc
    art::ServiceHandle<geo::Geometry>      geom;  
    if(geom->Nplanes() > 3) return;
    
    // get all hits in the event
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    evt.getByLabel(fHitsModuleLabel, hitListHandle);
    std::vector< art::Ptr<recob::Hit> > allhits;
    art::fill_ptr_vector(allhits, hitListHandle);
    if(allhits.size() == 0) return;

    // get clusters and cluster-hit associations
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);
    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);
    if(clusterListHandle->size() == 0) return;

    // get 3D vertices
    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);
    art::PtrVector<recob::Vertex> recoVtxList;
    double xyz[3] = {0,0,0};
    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii){
      art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      recoVtxList.push_back(vertex);
      vertex->XYZ(xyz);
//  mf::LogVerbatim("PFPAna")
//    <<"Reco Vtx "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2];
    }
    
    // get PFParticles
    art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
    evt.getByLabel(fPFParticleModuleLabel,PFPListHandle);
    art::PtrVector<recob::PFParticle> recoPFPList;
    for(unsigned int ii = 0; ii < PFPListHandle->size(); ++ii){
      art::Ptr<recob::PFParticle> pfp(PFPListHandle, ii);
      recoPFPList.push_back(pfp);
      mf::LogVerbatim("PFPAna")<<"PFParticle PDG "<<pfp->PdgCode();
    }
    
    // list of all true particles
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    sim::ParticleList const& plist = pi_serv->ParticleList();
    // list of all true particles that will be considered
    std::vector<const simb::MCParticle*> plist2;
    // true (reconstructed) hits for each particle in plist2
    std::vector<std::vector<art::Ptr<recob::Hit>>> hlist2;
    // index of cluster matched to each particle in plist2 in each plane
    std::vector<std::vector<short>> truToCl;
    // number of true hits in each plane and cluster
    std::vector<std::vector<unsigned short>> nTruHitInCl;
    //number of reconstructed hits in all clusters
    std::vector<unsigned short> nRecHitInCl;
    
    // calculate average EP2 for every event to facilitate code development
    // Beam Neutrinos - muons and not-muons

    float aveNuEP2mu = 0.;
    float numNuEP2mu = 0.;
    float aveNuEP2nm = 0.;
    float numNuEP2nm = 0.;
    // Cosmic Rays
    float aveCREP2 = 0.;
    float numCREP2 = 0.;

    // track ID of the neutrino
    int neutTrackID = -1;
    std::vector<int> tidlist;
    float neutEnergy = -1.;
    int neutIntType = -1;
    int neutCCNC = -1;
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      const simb::MCParticle* part = (*ipart).second;
      assert(part != 0);
      int pdg = abs(part->PdgCode());
      int trackID = part->TrackId();
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
      if(fSkipCosmics && theTruth->Origin() == simb::kCosmicRay) continue;

      if(fPrintLevel > 3) mf::LogVerbatim("PFPAna")
        <<"Pre-Cuts origin "<<theTruth->Origin()<<" trackID "<<trackID
        <<" PDG "<<part->PdgCode()
        <<" E "<<part->E()<<" mass "<<part->Mass()
        <<" Mother "<<part->Mother() + neutTrackID
        <<" Proc "<<part->Process();

      // Get the neutrino track ID. Assume that there is only one neutrino
      // interaction and it is first in the list of BeamNeutrino particles
      if(theTruth->Origin() == simb::kBeamNeutrino && neutTrackID < 0) {
        neutTrackID = trackID;
        simb::MCNeutrino theNeutrino = theTruth->GetNeutrino();
        neutEnergy = 1000. * theNeutrino.Nu().E();
        neutIntType = theNeutrino.InteractionType();
        neutCCNC = theNeutrino.CCNC();
//  mf::LogVerbatim("PFPAna")
//    <<"True Vtx "<<part->Vx()<<" "<<part->Vy()<<" "<<part->Vz();
        for(unsigned short iv = 0; iv < recoVtxList.size(); ++iv) {
          recoVtxList[iv]->XYZ(xyz);
          fNuVtx_dx->Fill(part->Vx() - xyz[0]);
          fNuVtx_dy->Fill(part->Vy() - xyz[1]);
          fNuVtx_dz->Fill(part->Vz() - xyz[2]);
        } // iv
      } // theTruth->Origin() == simb::kBeamNeutrino && neutTrackID < 

      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211)
         || (pdg == 321) || (pdg == 2212);


      if(!isCharged) continue;

      float KE = 1000 * (part->E() - part->Mass());
      // KE (MeV) cuts
      if(pdg ==   11) {
         if(fElecKERange[0] < 0) continue;
         // only allow primary electrons
         if(part->Process() != "primary") continue;
         if(KE < fElecKERange[0] || KE > fElecKERange[1]) continue;
      }
      if(pdg ==   13) {
         if(fMuonKERange[0] < 0) continue;
         if(KE < fMuonKERange[0] || KE > fMuonKERange[1]) continue;
      }
      if(pdg ==  211) {
         if(fPionKERange[0] < 0) continue;
         if(KE < fPionKERange[0] || KE > fPionKERange[1]) continue;
      }
      if(pdg ==  321) {
         if(fKaonKERange[0] < 0) continue;
         if(KE < fKaonKERange[0] || KE > fKaonKERange[1]) continue;
      }
      if(pdg == 2212) {
         if(fProtKERange[0] < 0) continue;
         if(KE < fProtKERange[0] || KE > fProtKERange[1]) continue;
      }
      // ignore secondaries from neutron interactions
      if(part->Process() == "NeutronInelastic") continue;
      plist2.push_back(part);
      tidlist.push_back(trackID);
      // initialize the true->(cluster,plane) association
      std::vector<short> temp {-1, -1, -1};
      truToCl.push_back(temp);
      // initialize the true hit count
      std::vector<unsigned short> temp2(3);
      nTruHitInCl.push_back(temp2);

      if(fPrintLevel > 2) mf::LogVerbatim("PFPAna")
        <<plist2.size() - 1
        <<" Origin "<<theTruth->Origin()<<" trackID "<<trackID
        <<" PDG "<<part->PdgCode()
        <<" KE "<<(int)KE
        <<" Mother "<<part->Mother() + neutTrackID
        <<" Proc "<<part->Process();
    }
    
    if(plist2.size() == 0) return;
    
    // get the hits (in all planes) that are matched to the true tracks
    hlist2 = bt_serv->TrackIdsToHits_Ps( tidlist, allhits);
    if(hlist2.size() != plist2.size()) {
      mf::LogError("PFPAna")
        <<"MC particle list size "<<plist2.size()
        <<" != size of MC particle true hits lists "<<hlist2.size();
      return;
    }
    tidlist.clear();

    // vector of (mother, daughter) pairs
    std::vector<std::pair<unsigned short, unsigned short>> moda;
    // Deal with mother-daughter tracks
    if(fMergeDaughters && neutTrackID >= 0) {
      // Assume that daughters appear later in the list. Step backwards
      // to accumulate all generations of daughters
      for(unsigned short dpl = plist2.size() - 1; dpl > 0; --dpl) {
        // no mother
        if(plist2[dpl]->Mother() == 0) continue;
        // electron
        if(abs(plist2[dpl]->PdgCode()) == 11) continue;
        // the actual mother trackID is offset from the neutrino trackID
        int motherID = neutTrackID + plist2[dpl]->Mother() - 1;
        // ensure that we are only looking at BeamNeutrino daughters
        if(motherID < 0) continue;
        // count the number of daughters
        int ndtr = 0;
        for(unsigned short kpl = 0; kpl < plist2.size(); ++kpl) {
          if(plist2[kpl]->Mother() == motherID) ++ndtr;
        }
        // require only one daughter
        if(ndtr > 1) continue;
        // find the mother in the list
        int mpl = -1;
        for(unsigned short jpl = dpl - 1; jpl > 0; --jpl) {
          if(plist2[jpl]->TrackId() == motherID) {
            mpl = jpl;
            break;
          }
        } // jpl
        // mother not found for some reason
        if(mpl < 0) continue;
        // ensure that PDG code for mother and daughter are the same
        if(plist2[dpl]->PdgCode() != plist2[mpl]->PdgCode()) continue;
        moda.push_back(std::make_pair(mpl, dpl));
      } //  dpl
    } // MergeDaughters
    
    // Now match reconstructed clusters to true particles.
    art::PtrVector<recob::Cluster> clusters;
    for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }
    
    fNClusters->Fill(clusterListHandle->size());
    nRecHitInCl.resize(clusters.size());
    
    // get the plane from the view. Perhaps there is a method that does
    // this somewhere...
    std::map< geo::View_t, unsigned int > ViewToPlane;
    for(unsigned int plane=0; plane < geom->Nplanes(); ++plane){
      geo::View_t view = geom->Plane(plane).View();
      ViewToPlane[view] = plane;
    }
    for(size_t icl = 0; icl < clusters.size(); ++icl){
      unsigned int plane = ViewToPlane[clusters[icl]->View()];
      std::vector< art::Ptr<recob::Hit> > cluhits = fmh.at(icl);
      fNHitInCluster->Fill(cluhits.size());
      nRecHitInCl[icl] = cluhits.size();
      // count the number of hits matched to each true particle in plist2
      std::vector<unsigned short> nHitInPl2(plist2.size());
      for(size_t iht = 0; iht < cluhits.size(); ++iht){
/*
  mf::LogVerbatim("PFPAna")
    <<"Clus Hit "<<cluhits[iht]->View()
    <<":"<<cluhits[iht]->WireID().Wire
    <<":"<<(int)cluhits[iht]->PeakTime();
*/
        // look for this hit in all of the truth hit lists
        short hitInPl2 = -1;
        for(unsigned short ipl = 0; ipl < plist2.size(); ++ipl) {
          unsigned short imat = 0;
          for(imat = 0; imat < hlist2[ipl].size(); ++imat) {
            if(cluhits[iht] == hlist2[ipl][imat]) break;
          } // imat
          if(imat < hlist2[ipl].size()) {
            hitInPl2 = ipl;
            break;
          }
        } // ipl
        if(hitInPl2 < 0) continue;
        // Assign the hit count to the mother if this is a daughter.
        // Mother-daughter pairs are entered in the moda vector in reverse
        // order, so assign daughter hits to the highest generation mother.
        for(unsigned short imd = 0; imd < moda.size(); ++imd) {
          if(moda[imd].second == hitInPl2) hitInPl2 = moda[imd].first;
        }
        // count
        ++nHitInPl2[hitInPl2];
      } // iht
      // Associate the cluster with the truth particle that has the highest
      // number of cluster hits
      unsigned short nhit = 0;
      short imtru = -1;
      for(unsigned int ipl = 0; ipl < nHitInPl2.size(); ++ipl) {
        if(nHitInPl2[ipl] > nhit) {
          nhit = nHitInPl2[ipl];
          imtru = ipl;
        }
      } // ipl
      // make the cluster->(true,plane) association and save the 
      // number of true hits in the cluster
      if(imtru != -1) {
        // clobber a previously made association?
        if(nhit > nTruHitInCl[imtru][plane]) {
          truToCl[imtru][plane] = icl;
          nTruHitInCl[imtru][plane] = nhit;
        }
      } // imtru != 1
    } // icl

    // ready to calculate Efficiency, Purity in each plane and EP2
    for(unsigned short ipl = 0; ipl < plist2.size(); ++ipl) {
      // ignore daughters
      bool skipit = false;
      for(unsigned short ii = 0; ii < moda.size(); ++ii) {
        if(moda[ii].second == ipl) {
          skipit = true;
          break;
        }
      } // ii
      if(skipit) continue;
      // ignore true particles with few true hits. Outside the detector
      // or not reconstructable?
      if(hlist2[ipl].size() < 3) continue;
      
      int trackID = plist2[ipl]->TrackId();
      art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
      bool isCosmic = (theTruth->Origin() == simb::kCosmicRay);
      float KE = 1000 * (plist2[ipl]->E() - plist2[ipl]->Mass());
      int PDG = abs(plist2[ipl]->PdgCode());
      
      std::vector<short> nTru(geom->Nplanes());
      std::vector<short> nRec(geom->Nplanes());
      std::vector<short> nTruRec(geom->Nplanes());
      std::vector<float> eff(geom->Nplanes());
      std::vector<float> pur(geom->Nplanes());
      std::vector<float> ep(geom->Nplanes());
      for(unsigned int plane = 0; plane < geom->Nplanes(); ++plane) {
        // count the number of true hits in this plane for the true particle.
        // First count the mother hits
        for(unsigned short ii = 0; ii < hlist2[ipl].size(); ++ii){
          if(ViewToPlane[hlist2[ipl][ii]->View()] == plane) ++nTru[plane];
        } // ii
//  mf::LogVerbatim("PFPAna")
//    <<"Chk mom "<<ipl<<" plane "<<plane<<" nTru "<<nTru[plane];
        // next look for daughters and count those hits in all generations
        unsigned short mom = ipl;
        std::vector<std::pair<unsigned short, unsigned short>>::reverse_iterator 
          rit = moda.rbegin();
        while(rit != moda.rend()) {
          if((*rit).first == mom) {
            unsigned short dau = (*rit).second;
            for(unsigned short jj = 0; jj < hlist2[dau].size(); ++jj) {
              if(ViewToPlane[hlist2[dau][jj]->View()] == plane) ++nTru[plane];
            } // jj
            // It is likely that one hit appears in the mother list
            // as well as the daughter list, so subtract one from the count
            --nTru[plane];
            mom = (*rit).second;
//  mf::LogVerbatim("PFPAna")<<"new mom "<<mom<<" nTru "<<nTru[plane];
          } // (*rit).first == mom
          ++rit;
        } // rit
//  mf::LogVerbatim("PFPAna")<<"Chk dau "<<nTru[plane];
        if(nTru[plane] == 0) {
//          mf::LogVerbatim("PFPAna")<<"No true hits in plane "<<plane
//            <<" for truth particle "<<ipl;
          continue;
        }
        short icl = truToCl[ipl][plane];
        nRec[plane] = nRecHitInCl[icl];
        nTruRec[plane] = nTruHitInCl[ipl][plane];
//  mf::LogVerbatim("PFPAna")<<"icl "<<icl<<" nRec "<<nRec[plane]
//    <<" nTruRec "<<nTruRec[plane];
        if(nTru[plane] > 0) 
          eff[plane] = (float)nTruRec[plane] / (float)nTru[plane];
        if(nRec[plane] > 0) 
          pur[plane] = (float)nTruRec[plane] / (float)nRec[plane];
        ep[plane] = eff[plane] * pur[plane];
      } // plane
      // sort the ep values in ascending order
      std::vector<float> temp;
      temp = ep;
      std::sort(temp.begin(), temp.end());
      // EP2 is the second highest value
      unsigned short ii = temp.size() - 2;
      float ep2 = temp[ii];
      // find the plane that defined EP2
      short ep2Plane = 0;
      short ep2Cluster = 0;
      for(unsigned short jj = 0; jj < temp.size(); ++jj) {
        if(ep[jj] == ep2) {
          ep2Plane = jj;
          ep2Cluster = truToCl[ipl][ep2Plane];
          break;
        }
      } // jj
      // find the US and DS ends of the cluster for printing
      std::array<double, 2> clBeg, clEnd;
      if(ep2Cluster >= 0) {
        clBeg[0] = clusters[ep2Cluster]->StartWire();
        clBeg[1] = clusters[ep2Cluster]->StartTick();
        clEnd[0] = clusters[ep2Cluster]->EndWire();
        clEnd[1] = clusters[ep2Cluster]->EndTick();
      }
      else {
        clBeg.fill(0.);
        clEnd.fill(0.);
      }
      // fill histograms
      if(isCosmic) {
        fCREP2->Fill(ep2);
        fCRE->Fill(eff[ep2Plane]);
        fCRP->Fill(pur[ep2Plane]);
        aveCREP2 += ep2;
        numCREP2 += 1.;
        if(fPrintLevel > 1) mf::LogVerbatim("PFPAna")
          <<">>>CREP2 "<<std::fixed<<std::setprecision(2)<<ep2
          <<" E "<<eff[ep2Plane]<<std::setprecision(2)<<" P "<<pur[ep2Plane]
          <<" P:W:T "<<ep2Plane<<":"<<(int)clBeg[0]<<":"<<(int)clBeg[1]
          <<"-"<<ep2Plane<<":"<<(int)clEnd[0]<<":"<<(int)clEnd[1]
          <<" PDG "<<PDG<<" KE "<<(int)KE<<" MeV";
      } // isCosmic
      else {
        float wght = 1.;
        if(fTrackWeightOption == 1) wght = KE;
        // accumulate statistics for muons and not-muons
        if(PDG == 13) {
          aveNuEP2mu += ep2 * wght;
          numNuEP2mu += wght;
        } else {
          aveNuEP2nm += ep2 * wght;
          numNuEP2nm += wght;
        }
        if(PDG == 11) {
          fNuKE_elec->Fill(KE, wght);
          fNuE_elec->Fill(eff[ep2Plane], wght);
          fNuP_elec->Fill(pur[ep2Plane], wght);
          fNuEP2_elec->Fill(ep2, wght);
          fNuEP2_KE_elec->Fill(KE, ep2, wght);
        } else if(PDG == 13) {
          fNuKE_muon->Fill(KE, wght);
          fNuE_muon->Fill(eff[ep2Plane], wght);
          fNuP_muon->Fill(pur[ep2Plane], wght);
          fNuEP2_muon->Fill(ep2, wght);
          fNuEP2_KE_muon->Fill(KE, ep2, wght);
        } else if(PDG == 211) {
          fNuKE_pion->Fill(KE, wght);
          fNuE_pion->Fill(eff[ep2Plane], wght);
          fNuP_pion->Fill(pur[ep2Plane], wght);
          fNuEP2_pion->Fill(ep2, wght);
          fNuEP2_KE_pion->Fill(KE, ep2, wght);
        } else if(PDG == 321) {
          fNuKE_kaon->Fill(KE, wght);
          fNuE_kaon->Fill(eff[ep2Plane], wght);
          fNuP_kaon->Fill(pur[ep2Plane], wght);
          fNuEP2_kaon->Fill(ep2, wght);
          fNuEP2_KE_kaon->Fill(KE, ep2, wght);
        } else if(PDG == 2212) {
          fNuKE_prot->Fill(KE, wght);
          fNuE_prot->Fill(eff[ep2Plane], wght);
          fNuP_prot->Fill(pur[ep2Plane], wght);
          fNuEP2_prot->Fill(ep2, wght);
          fNuEP2_KE_prot->Fill(KE, ep2, wght);
        }
        if(fPrintLevel > 1) mf::LogVerbatim("PFPAna")
          <<">>>NuEP2 "<<std::fixed<<std::setprecision(2)<<ep2
          <<" E "<<eff[ep2Plane]<<std::setprecision(2)<<" P "<<pur[ep2Plane]
          <<" P:W:T "<<ep2Plane<<":"<<(int)clBeg[0]<<":"<<(int)clBeg[1]
          <<"-"<<ep2Plane<<":"<<(int)clEnd[0]<<":"<<(int)clEnd[1]
          <<" PDG "<<PDG<<" KE "<<(int)KE<<" MeV ";
        if(fPrintLevel > 2) {
          // print out the begin/end true hits
          mf::LogVerbatim mfp("PFPAna");
          mfp<<" Truth P:W:T ";
          for(unsigned int plane = 0; plane < geom->Nplanes(); ++plane) {
            unsigned short loW = 9999;
            unsigned short loT = 0;
            unsigned short hiW = 0;
            unsigned short hiT = 0;
            for(unsigned short ii = 0; ii < hlist2[ipl].size(); ++ii){
              if(ViewToPlane[hlist2[ipl][ii]->View()] == plane) {
                art::Ptr<recob::Hit> theHit = hlist2[ipl][ii];
                if(theHit->WireID().Wire < loW) {
                  loW = theHit->WireID().Wire;
                  loT = theHit->PeakTime();
                }
                if(theHit->WireID().Wire > hiW) {
                  hiW = theHit->WireID().Wire;
                  hiT = theHit->PeakTime();
                }
              } // correct view
            } // ii
            mfp<<plane<<":"<<loW<<":"<<loT<<"-"<<plane<<":"<<hiW<<":"<<hiT<<" ";
          } // plane
        } // fPrintLevel > 2
      } // !isCosmic
    } // ipl

  float ave1 = -1.;
  if(numNuEP2mu > 0.) ave1 = aveNuEP2mu/numNuEP2mu;

  float ave2 = -1.;
  if(numNuEP2nm > 0.) ave2 = aveNuEP2nm/numNuEP2nm;

  float ave3 = -1.;
  if(numCREP2 > 0.) ave3 = aveCREP2/numCREP2;

  

    if(fPrintLevel > 0) {
      std::string nuType = "Other";
      if(neutCCNC == simb::kCC) {
        if(neutIntType == 1001) nuType = "CCQE";
        if(neutIntType == 1091) nuType = "DIS";
        if(neutIntType == 1097) nuType = "COH";
        if(neutIntType > 1002 && neutIntType < 1091) nuType = "RES";
      } else if(neutCCNC == simb::kNC) {
        nuType = "NC";
      } else {
        nuType = "Unknown";
      }
      mf::LogVerbatim("PFPAna")
      <<"EvtEP2 "<<evt.id().event()
      <<" NuType "<<nuType
      <<" Enu "<<std::fixed<<std::setprecision(0)<<neutEnergy
      <<std::right<<std::fixed<<std::setprecision(2)
      <<" NuMuons "<<ave1
      <<" NuPiKp "<<ave2
      <<" CosmicRays "<<ave3
      <<" CCNC "<<neutCCNC<<" IntType "<<neutIntType;
    }  
  } // analyze

  DEFINE_ART_MODULE(PFPAna)

  
} //end namespace

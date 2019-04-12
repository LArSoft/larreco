////////////////////////////////////////////////////////////////////////
//
// HitFinderAna class
//
// echurch@fnal.gov
//
//  This algorithm is designed to analyze hits on wires after deconvolution
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TFile.h>

// C++ includes
#include <algorithm>
#include <fstream>
#include <bitset>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/LArFFT.h"

#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"

#include "art/Framework/Core/EDAnalyzer.h"

#include <string>

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace hit {

  /// Base class for creation of raw signals on wires.
  class HitFinderAna : public art::EDAnalyzer {

  public:

    explicit HitFinderAna(fhicl::ParameterSet const& pset);
    virtual ~HitFinderAna();

    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::string            fFFTHitFinderModuleLabel;
    std::string            fLArG4ModuleLabel;

      TTree* fHTree;
      Int_t fRun;
      Int_t fEvt;
      Int_t fNp0;
      Int_t fNp1;
      Int_t fNp2;
      Int_t fN3p0;
      Int_t fN3p1;
      Int_t fN3p2;
      Float_t* fTimep0;
      Float_t* fTimep1;
      Float_t* fTimep2;
      Int_t* fWirep0;
      Int_t* fWirep1;
      Int_t* fWirep2;
      Float_t* fChgp0;
      Float_t* fChgp1;
      Float_t* fChgp2;
      Float_t* fXYZp0;
      Float_t* fXYZp1;
      Float_t* fXYZp2;

      Int_t*  fMCPdg0;
      Int_t*  fMCTId0;
      Float_t*  fMCE0;
      Int_t*  fMCPdg1;
      Int_t*  fMCTId1;
      Float_t*  fMCE1;
      Int_t*  fMCPdg2;
      Int_t*  fMCTId2;
      Float_t*  fMCE2;

  }; // class HitFinderAna

}

namespace hit{

  //-------------------------------------------------
  HitFinderAna::HitFinderAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  HitFinderAna::~HitFinderAna()
  {
  }

  void HitFinderAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fFFTHitFinderModuleLabel = p.get< std::string >("HitsModuleLabel");
    fLArG4ModuleLabel        = p.get< std::string >("LArGeantModuleLabel");
    return;
  }
  //-------------------------------------------------
  void HitFinderAna::beginJob()
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;
    fNp0 = 9000;
    fNp1 = 9000;
    fNp2 = 9000;

    fHTree = tfs->make<TTree>("HTree","HTree");
    fTimep0 = new Float_t[fNp0];
    fTimep1 = new Float_t[fNp1];
    fTimep2 = new Float_t[fNp2];
    fWirep0 = new Int_t[fNp0];
    fWirep1 = new Int_t[fNp1];
    fWirep2 = new Int_t[fNp2];
    fChgp0 = new Float_t[fNp0];
    fChgp1 = new Float_t[fNp1];
    fChgp2 = new Float_t[fNp2];
    fXYZp0 = new Float_t[fNp0*3];
    fXYZp1 = new Float_t[fNp1*3];
    fXYZp2 = new Float_t[fNp2*3];

    fMCPdg0 = new Int_t[fNp0];
    fMCPdg1 = new Int_t[fNp1];
    fMCPdg2 = new Int_t[fNp2];
    fMCTId0 = new Int_t[fNp0];
    fMCTId1 = new Int_t[fNp1];
    fMCTId2 = new Int_t[fNp2];
    fMCE0 = new Float_t[fNp0];
    fMCE1 = new Float_t[fNp1];
    fMCE2 = new Float_t[fNp2];

    fHTree->Branch("HEvt", &fEvt, "HEvt/I");
    fHTree->Branch("HRun", &fRun, "HRun/I");
    fHTree->Branch("HNp0", &fNp0, "HNp0/I");
    fHTree->Branch("HNp1", &fNp1, "HNp1/I");
    fHTree->Branch("HNp2", &fNp2, "HNp2/I");
    fHTree->Branch("HN3p0", &fN3p0, "HN3p0/I");
    fHTree->Branch("HN3p1", &fN3p1, "HN3p1/I");
    fHTree->Branch("HN3p2", &fN3p2, "HN3p2/I");
    fHTree->Branch("Htp0", fTimep0, "Htp0[HNp0]/F");
    fHTree->Branch("Htp1", fTimep1, "Htp1[HNp1]/F");
    fHTree->Branch("Htp2", fTimep2, "Htp2[HNp2]/F");
    fHTree->Branch("Hwp0", fWirep0, "Hwp0[HNp0]/I");
    fHTree->Branch("Hwp1", fWirep1, "Hwp1[HNp1]/I");
    fHTree->Branch("Hwp2", fWirep2, "Hwp2[HNp2]/I");
    fHTree->Branch("Hchgp0", fChgp0, "Hchgp0[HNp0]/F");
    fHTree->Branch("Hchgp1", fChgp1, "Hchgp1[HNp1]/F");
    fHTree->Branch("Hchgp2", fChgp2, "Hchgp2[HNp2]/F");
    fHTree->Branch("HMCXYZp0", fXYZp0, "HMCXYZp0[HN3p0]/F");
    fHTree->Branch("HMCXYZp1", fXYZp1, "HMCXYZp1[HN3p1]/F");
    fHTree->Branch("HMCXYZp2", fXYZp2, "HMCXYZp2[HN3p2]/F");
    fHTree->Branch("HMCPdgp0", fMCPdg0, "HMCPdgp0[HNp0]/I");
    fHTree->Branch("HMCPdgp1", fMCPdg1, "HMCPdgp1[HNp1]/I");
    fHTree->Branch("HMCPdgp2", fMCPdg2, "HMCPdgp2[HNp2]/I");
    fHTree->Branch("HMCTIdp0", fMCTId0, "HMCTIdp0[HNp0]/I");
    fHTree->Branch("HMCTIdp1", fMCTId1, "HMCTIdp1[HNp1]/I");
    fHTree->Branch("HMCTIdp2", fMCTId2, "HMCTIdp2[HNp2]/I");
    fHTree->Branch("HMCEp0", fMCE0, "HMCEp0[HNp0]/F");
    fHTree->Branch("HMCEp1", fMCE1, "HMCEp1[HNp1]/F");
    fHTree->Branch("HMCEp2", fMCE2, "HMCEp2[HNp2]/F");


    return;

  }

  //-------------------------------------------------
  void HitFinderAna::analyze(const art::Event& evt)
  {

    if (evt.isRealData()){
      throw cet::exception("HitFinderAna: ") << "Not for use on Data yet...\n";
    }

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fFFTHitFinderModuleLabel,hitHandle);

    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;

    sim::ParticleList const& _particleList = pi_serv->ParticleList();

    MF_LOG_VERBATIM("HitFinderAna") << _particleList;

    art::ServiceHandle<geo::Geometry const> geom;

    std::map<geo::PlaneID, std::vector< art::Ptr<recob::Hit> > > planeIDToHits;
    for(size_t h = 0; h < hitHandle->size(); ++h)
      planeIDToHits[hitHandle->at(h).WireID().planeID()].push_back(art::Ptr<recob::Hit>(hitHandle, h));


    for(auto mapitr : planeIDToHits){
      fNp0=0;       fN3p0=0;
      fNp1=0;       fN3p1=0;
      fNp2=0;       fN3p2=0;

      geo::PlaneID pid = mapitr.first;
      auto itr = mapitr.second.begin();
      while(itr != mapitr.second.end()) {

	fRun = evt.run();
	fEvt = evt.id().event();

	std::vector<sim::TrackIDE> trackides = bt_serv->HitToTrackIDEs(*itr);
	std::vector<sim::TrackIDE>::iterator idesitr = trackides.begin();
	std::vector<double> xyz = bt_serv->HitToXYZ(*itr);

	if (pid.Plane == 0 && fNp0 < 9000){
	  fTimep0[fNp0] = (*itr)->PeakTime();
	  fWirep0[fNp0] = (*itr)->WireID().Wire;
	  fChgp0[fNp0] = (*itr)->Integral();

	  for (unsigned int kk = 0; kk < 3; ++kk){
	    fXYZp0[fNp0*3+kk] = xyz[kk];
	  }


	  while( idesitr != trackides.end() ){
	    fMCTId0[fNp0] = (*idesitr).trackID;
	    if (_particleList.find((*idesitr).trackID) != _particleList.end()){
	      const simb::MCParticle* particle = _particleList.at( (*idesitr).trackID);
	      fMCPdg0[fNp0] = particle->PdgCode();
	      fMCE0[fNp0] = particle->E();
	    }
	    idesitr++;
	  }

	  ++fNp0;
	}

	else if (pid.Plane == 1 && fNp1 < 9000){
	  fTimep1[fNp1] = (*itr)->PeakTime();
	  fWirep1[fNp1] = (*itr)->WireID().Wire;
	  fChgp1[fNp1] = (*itr)->Integral();

	  for (unsigned int kk = 0; kk < 3; ++kk){
	    fXYZp1[fNp1*3+kk] = xyz[kk];
	  }

	  while( idesitr != trackides.end() ){
	    fMCTId1[fNp1] = (*idesitr).trackID;
	    if (_particleList.find((*idesitr).trackID) != _particleList.end()){
	      const simb::MCParticle* particle = _particleList.at( (*idesitr).trackID);
	      fMCPdg1[fNp1] = particle->PdgCode();
	      fMCE1[fNp1] = particle->E();
	    }
	    idesitr++;
	  }
	  ++fNp1;
	}

	else if (pid.Plane == 2  && fNp2 < 9000){
	  fTimep2[fNp2] = (*itr)->PeakTime();
	  fWirep2[fNp2] = (*itr)->WireID().Wire;
	  fChgp2[fNp2] = (*itr)->Integral();

	  for (unsigned int kk = 0; kk < 3; ++kk){
	    fXYZp2[fNp2*3+kk] = xyz[kk];
	  }

	  while( idesitr != trackides.end()){
	    fMCTId2[fNp2] = (*idesitr).trackID;
	    if (_particleList.find((*idesitr).trackID) != _particleList.end() ){
	      const simb::MCParticle* particle = _particleList.at( (*idesitr).trackID);
	      fMCPdg2[fNp2] = particle->PdgCode();
	      fMCE2[fNp2] = particle->E();
	    }
	    idesitr++;
	  }
	  ++fNp2;
	}

	fN3p0 = 3* fNp0;
	fN3p1 = 3* fNp1;
	fN3p2 = 3* fNp2;

	fHTree->Fill();
	itr++;
      } // loop on Hits
    } // loop on map

    return;
  }//end analyze method

}//end namespace

namespace hit{

  DEFINE_ART_MODULE(HitFinderAna)

} // end of hit namespace

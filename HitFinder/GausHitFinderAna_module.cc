#ifndef GAUSHITFINDERANA_H
#define GAUSHITFINDERANA_H
////////////////////////////////////////////////////////////////////////
//
// GausHitFinder class designed to analyze signal on a wire in the TPC
//
// jaasaadi@syr.edu
//
// Note: This Ana module has been based (stolen) from the FFTHitFinderAna thus
// there is still some unneeded hold overs that will get cleaned up later
////////////////////////////////////////////////////////////////////////


// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TFile.h>
#include "TComplex.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <bitset>
#include <vector>
#include <string>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"

namespace geo { class Geometry;   }
namespace sim { class SimChannel; }

namespace hit{

  /// Base class for creation of raw signals on wires. 
  class GausHitFinderAna : public art::EDAnalyzer {
    
  public:
        
    explicit GausHitFinderAna(fhicl::ParameterSet const& pset); 
    virtual ~GausHitFinderAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::string            fGausHitFinderModuleLabel;
    std::string            fLArG4ModuleLabel;
    std::string		   fHitCheaterModuleLabel;
    std::vector<const sim::SimChannel*>    fSimChannels;           ///< all the SimChannels for the event
      
    TH1F* fRun;
    TH1F* fEvt;
    TH1F* fWireNumbMulti1; //<---Wire number of hit with multiplicity of 1
    TH1F* fFitGoodnessMulti1; //<---Goodness of hit with mulitplicity of 1
    TH1F* fChargeMulti1; //<---Charge of hit with multiplicity of 1
      
    TH1F* fTruthPeakPosition;//<---Peak time of hit from backtracker with multiplicity of 1
    TH1F* fTruthPeakPositionPlane0Multi1;
    TH1F* fTruthPeakPositionPlane1Multi1;
    TH1F* fTruthPeakPositionPlane2Multi1;
      
    TH1F* fRecoPeakPositionMulti1;//<---Peak time of hit with multiplicity of 1
    TH1F* fRecoPeakPositionPlane0Multi1;//<---Peak time of hit with multiplicity of 1 in plane 0
    TH1F* fRecoPeakPositionPlane1Multi1;//<---Peak time of hit with multiplicity of 1 in plane 1
    TH1F* fRecoPeakPositionPlane2Multi1;//<---Peak time of hit with multiplicity of 1 in plane 2
      
    TH1F* fRecoPeakPositionUncertMulti1;//<---Peak time of hit with multiplicity of 1
    TH1F* fRecoPeakPositionUncertPlane0Multi1;//<---Peak time of hit with multiplicity of 1 in plane 0
    TH1F* fRecoPeakPositionUncertPlane1Multi1;//<---Peak time of hit with multiplicity of 1 in plane 1
    TH1F* fRecoPeakPositionUncertPlane2Multi1;//<---Peak time of hit with multiplicity of 1 in plane 2
      
    TH1F* fRecoPeakPositionUncertMultiGT1;//<---Peak time of hit with multiplicity of > 1
    TH1F* fRecoPeakPositionUncertPlane0MultiGT1;//<---Peak time of hit with multiplicity of > 1 in plane 0
    TH1F* fRecoPeakPositionUncertPlane1MultiGT1;//<---Peak time of hit with multiplicity of > 1 in plane 1
    TH1F* fRecoPeakPositionUncertPlane2MultiGT1;//<---Peak time of hit with multiplicity of > 1 in plane 2
      
    TH1F* fHitResidualAll;
    TH1F* fHitResidualMulti1;
    TH1F* fHitResidualPlane0Multi1;
    TH1F* fHitResidualPlane1Multi1;
    TH1F* fHitResidualPlane2Multi1;
      
    TTree* fHTree;
    //Int_t fRun;
    //Int_t fEvt;
    Int_t fnhits; //<---Number of Hits in the Event
    Int_t fnOnePulseHits; //<---Number of One Pulse hit events
    Int_t fmulitPulseHits; //<---Number of Multi pulse hit events
      
      
    Int_t fSingleHit; //<<---Set a indicator to know if this is a single pulse hit or multihit
    Int_t fMultiHit;  //<<---Set a indicator to know if this is a single pulse hit or multihit
      
    Float_t fgoodoffitn1; //<---Goodness of hit with mulitplicity of 1
    Float_t fChargen1; //<---Charge of hit with multiplicity of 1
    Float_t fSigmaChargen1; //<---Uncertainty of charge of hit with multiplicity of 1
    Float_t fWidthn1; //<---Width of the Hit (Endtime - Peaktime) with multiplicity of 1
    Float_t fPeakn1; //<---Peak time of hit with multiplicity of 1
    Float_t fPeakUncertn1; //<---Uncertainty in the peak position of the hit with multiplicity of 1
    Float_t fStartTimen1; //<---Start Position of the hit with multiplicity of 1
    Float_t fStartTimeUncertn1; //<---Start Time Uncertainty of the hit with multiplicity of 1
    Float_t fEndTimen1; //<---End Position of the hit with multiplicity of 1
    Float_t fEndTimeUncertn1; //<---End Time Uncertainty of the hit with multiplicity of 1
      
      
    Int_t fWirenGT1; //<---Wire number of hit with multiplicity of 1
    Float_t fWidthnGT1; //<---Width of the Hit (Endtime - Peaktime) with multiplicity > 1
    Float_t fPeaknGT1; //<---Peak time of hit with multiplicity > 1
    Float_t fPeakUncertnGT1; //<---Uncertainty in the peak position of the hit with multiplicity > 1
    Float_t fSigmaChargenGT1; //<---Uncertainty of charge of hit with multiplicity > 1
    Float_t fgoodoffitnGT1; //<---Goodness of hit with mulitplicity > 1
    Float_t fChargenGT1; //<---Charge of hit with multiplicity > 1
    Float_t fStartTimenGT1; //<---Start Position of the hit with multiplicity > 1
    Float_t fStartTimeUncertnGT1; //<---Start Time Uncertainty of the hit with multiplicity > 1
    Float_t fEndTimenGT1; //<---End Position of the hit with multiplicity > 1
    Float_t fEndTimeUncertnGT1; //<---End Time Uncertainty of the hit with multiplicity > 1
      
      

    Int_t fN3p0;
    Int_t fN3p1;
    Int_t fN3p2;
    Float_t* fPeakTime0;
    Float_t* fPeakTime1;
    Float_t* fPeakTime2;
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
      
  }; // class GausHitFinderAna


  //-------------------------------------------------
  GausHitFinderAna::GausHitFinderAna(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  GausHitFinderAna::~GausHitFinderAna()
  {
  }

  void GausHitFinderAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fGausHitFinderModuleLabel = p.get< std::string >("HitsModuleLabel");
    fLArG4ModuleLabel        = p.get< std::string >("LArGeantModuleLabel");
    return;
  }
  //-------------------------------------------------
  void GausHitFinderAna::beginJob() 
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    // ====================================
    // ==== Outputting TH1F Histograms ====
    // ====================================
    fRun 		= tfs->make<TH1F>("fRun", "Run Number", 1000, 0, 1000);
    fEvt 		= tfs->make<TH1F>("fEvt", "Event Number", 1000, 0, 1000);
    fWireNumbMulti1 	= tfs->make<TH1F>("fWireNumbMulti1", "Wire Number Multi 1", 4000, 0, 4000);
    fFitGoodnessMulti1	= tfs->make<TH1F>("fFitGoodnessMulti1", "Fit Goodness Multi 1", 200, 0, 100);
    fChargeMulti1	= tfs->make<TH1F>("fChargeMulti1", "Charge Multi 1", 5000, 0, 5000);
    
    // =================================================================
    // === Gaussian Hit Peak Position for hits with multiplicity = 1 ===
    fRecoPeakPositionMulti1	   = tfs->make<TH1F>("fRecoPeakPositionMulti1", "Reco Peak Position Multi1", 500, 0, 2000);
    fRecoPeakPositionPlane0Multi1  = tfs->make<TH1F>("fRecoPeakPositionPlane0Multi1", "Reco Peak Position Plane 0 Multi1", 500, 0, 2000);
    fRecoPeakPositionPlane1Multi1  = tfs->make<TH1F>("fRecoPeakPositionPlane1Multi1", "Reco Peak Position Plane 1 Multi1", 500, 0, 2000);
    fRecoPeakPositionPlane2Multi1  = tfs->make<TH1F>("fRecoPeakPositionPlane2Multi1", "Reco Peak Position Plane 2 Multi1", 500, 0, 2000);
    
    // =============================================================================
    // === Gaussian Hit Peak Position uncertainty for hits with multiplicity = 1 ===
    fRecoPeakPositionUncertMulti1	 = tfs->make<TH1F>("fRecoPeakPositionUncertMulti1", "Reco Peak Position Uncert Multi1", 100, 0, 1);
    fRecoPeakPositionUncertPlane0Multi1  = tfs->make<TH1F>("fRecoPeakPositionUncertPlane0Multi1", "Reco Peak Position Uncert Plane 0 Multi1", 100, 0, 1);
    fRecoPeakPositionUncertPlane1Multi1  = tfs->make<TH1F>("fRecoPeakPositionUncertPlane1Multi1", "Reco Peak Position Uncert Plane 1 Multi1", 100, 0, 1);
    fRecoPeakPositionUncertPlane2Multi1  = tfs->make<TH1F>("fRecoPeakPositionUncertPlane2Multi1", "Reco Peak Position Uncert Plane 2 Multi1", 100, 0, 1);
    
    // =============================================================================
    // === Gaussian Hit Peak Position uncertainty for hits with multiplicity > 1 ===
    fRecoPeakPositionUncertMultiGT1	 = tfs->make<TH1F>("fRecoPeakPositionUncertMultiGT1", "Reco Peak Position Uncert Multi1", 100, 0, 1);
    fRecoPeakPositionUncertPlane0MultiGT1  = tfs->make<TH1F>("fRecoPeakPositionUncertPlane0MultiGT1", "Reco Peak Position Uncert Plane 0 Multi1", 100, 0, 1);
    fRecoPeakPositionUncertPlane1MultiGT1  = tfs->make<TH1F>("fRecoPeakPositionUncertPlane1MultiGT1", "Reco Peak Position Uncert Plane 1 Multi1", 100, 0, 1);
    fRecoPeakPositionUncertPlane2MultiGT1  = tfs->make<TH1F>("fRecoPeakPositionUncertPlane2MultiGT1", "Reco Peak Position Uncert Plane 2 Multi1", 100, 0, 1);
    
    // ============================================
    // === Truth Peak position from BackTracker ===
    fTruthPeakPosition             = tfs->make<TH1F>("fTruthPeakPosition", "Truth Peak Position Multi1", 500, 0, 200);
    fTruthPeakPositionPlane0Multi1 = tfs->make<TH1F>("fTruthPeakPositionPlane0Multi1", "Truth Peak Position Plane 0 Multi1", 500, 0, 2000);
    fTruthPeakPositionPlane1Multi1 = tfs->make<TH1F>("fTruthPeakPositionPlane1Multi1", "Truth Peak Position Plane 1 Multi1", 500, 0, 2000);
    fTruthPeakPositionPlane2Multi1 = tfs->make<TH1F>("fTruthPeakPositionPlane2Multi1", "Truth Peak Position Plane 2 Multi1", 500, 0, 2000);
    
    
    fHitResidualAll	     = tfs->make<TH1F>("fHitResidualAll", "Hit Residual All", 1600, -400, 400);
    fHitResidualMulti1	     = tfs->make<TH1F>("fHitResidualMulti1", "Hit Residual Multi 1", 1600, -400, 400);
    fHitResidualPlane0Multi1 = tfs->make<TH1F>("fHitResidualPlane0Multi1", "Hit Residual Plane 0 Multi 1", 1600, -400, 400);
    fHitResidualPlane1Multi1 = tfs->make<TH1F>("fHitResidualPlane1Multi1", "Hit Residual Plane 1 Multi 1", 1600, -400, 400);
    fHitResidualPlane2Multi1 = tfs->make<TH1F>("fHitResidualPlane2Multi1", "Hit Residual Plane 2 Multi 1", 1600, -400, 400);
    
    fHTree = tfs->make<TTree>("HTree","HTree");


    //fHTree->Branch("HEvt", &fEvt, "HEvt/I");
    //fHTree->Branch("HRun", &fRun, "HRun/I");
    fHTree->Branch("NHits", &fnhits, "NHits/I");
    fHTree->Branch("SinglePulseEvent", &fSingleHit, "SingleHitEvent/I");
    fHTree->Branch("MulitPulseEvent", &fMultiHit, "MulitPulseEvent/I");

    
    

    fHTree->Branch("NOnePulseHit", &fnOnePulseHits, "NOnePulseHit/I");
    
    fHTree->Branch("ChargeMulti1", &fChargen1, "ChargeMulti1/F");
    fHTree->Branch("SigmaChargeMulti1", &fSigmaChargen1, "SigmaChargeMulti1/F");
    fHTree->Branch("WidthMulti1", &fWidthn1, "WidthMulti1/F");
    fHTree->Branch("PeakPosMulti1", &fPeakn1, "PeakPosMulti1/F");
    fHTree->Branch("PeakUncertMulti1", &fPeakUncertn1, "PeakUncertMulti1/F");
    fHTree->Branch("StartPosMulti1", &fStartTimen1, "StartPosMulti1/F");
    fHTree->Branch("StartPosUncertMulti1", &fStartTimeUncertn1, "StartPosUncertMulti1/F");
    fHTree->Branch("EndPosMulti1", &fEndTimen1, "EndPosMulti1/F");
    fHTree->Branch("EndPosUncertMulti1", &fEndTimeUncertn1, "EndPosUncertMulti1/F");
    
    
    fHTree->Branch("WireNumbernGT1", &fWirenGT1, "WireNumbernGT1/I");
    fHTree->Branch("NMultiPulseHit", &fmulitPulseHits, "NMultiPulseHit/I");
    fHTree->Branch("GOFMultiGT1", &fgoodoffitnGT1, "GOFMultiGT1/F");
    fHTree->Branch("ChargeMultiGT1", &fChargenGT1, "ChargeMultiGT1/F");
    fHTree->Branch("SigmaChargeMultiGT1", &fSigmaChargenGT1, "SigmaChargeMultiGT1/F");
    fHTree->Branch("WidthMultiGT1", &fWidthnGT1, "WidthMultiGT1/F");
    fHTree->Branch("PeakPosMultiGT1", &fPeaknGT1, "PeakPosMultiGT1/F");
    fHTree->Branch("PeakUncertMultiGT1", &fPeakUncertnGT1, "PeakUncertMultiGT1/F");
    fHTree->Branch("StartPosMultiGT1", &fStartTimenGT1, "StartPosMultiGT1/F");
    fHTree->Branch("StartPosUncertMultiGT1", &fStartTimeUncertnGT1, "StartPosUncertMultiGT1/F");
    fHTree->Branch("EndPosMultiGT1", &fEndTimenGT1, "EndPosMultiGT1/F");
    fHTree->Branch("EndPosUncertMultiGT1", &fEndTimeUncertnGT1, "EndPosUncertMultiGT1/F");
    
    
    /*fHTree->Branch("Hit0PeakTime", &fNp0, "Hit0PeakTime/I");
      fHTree->Branch("HNp1", &fNp1, "HNp1/I");
      fHTree->Branch("HNp2", &fNp2, "HNp2/I");
      fHTree->Branch("HN3p0", &fN3p0, "HN3p0/I");
      fHTree->Branch("HN3p1", &fN3p1, "HN3p1/I");
      fHTree->Branch("HN3p2", &fN3p2, "HN3p2/I");
      fHTree->Branch("Htp0", fPeakTime0, "Htp0[Hit0PeakTime]/F");
      fHTree->Branch("Htp1", fPeakTime1, "Htp1[HNp1]/F");
      fHTree->Branch("Htp2", fPeakTime2, "Htp2[HNp2]/F");
      fHTree->Branch("Hwp0", fWirep0, "Hwp0[Hit0PeakTime]/I");
      fHTree->Branch("Hwp1", fWirep1, "Hwp1[HNp1]/I");
      fHTree->Branch("Hwp2", fWirep2, "Hwp2[HNp2]/I");
      fHTree->Branch("Hchgp0", fChgp0, "Hchgp0[Hit0PeakTime]/F");
      fHTree->Branch("Hchgp1", fChgp1, "Hchgp1[HNp1]/F");
      fHTree->Branch("Hchgp2", fChgp2, "Hchgp2[HNp2]/F");
      fHTree->Branch("HMCXYZp0", fXYZp0, "HMCXYZp0[HN3p0]/F");
      fHTree->Branch("HMCXYZp1", fXYZp1, "HMCXYZp1[HN3p1]/F");
      fHTree->Branch("HMCXYZp2", fXYZp2, "HMCXYZp2[HN3p2]/F");
      fHTree->Branch("HMCPdgp0", fMCPdg0, "HMCPdgp0[Hit0PeakTime]/I");
      fHTree->Branch("HMCPdgp1", fMCPdg1, "HMCPdgp1[HNp1]/I");
      fHTree->Branch("HMCPdgp2", fMCPdg2, "HMCPdgp2[HNp2]/I");
      fHTree->Branch("HMCTIdp0", fMCTId0, "HMCTIdp0[Hit0PeakTime]/I");
      fHTree->Branch("HMCTIdp1", fMCTId1, "HMCTIdp1[HNp1]/I");
      fHTree->Branch("HMCTIdp2", fMCTId2, "HMCTIdp2[HNp2]/I");
      fHTree->Branch("HMCEp0", fMCE0, "HMCEp0[Hit0PeakTime]/F");
      fHTree->Branch("HMCEp1", fMCE1, "HMCEp1[HNp1]/F");
      fHTree->Branch("HMCEp2", fMCE2, "HMCEp2[HNp2]/F");*/

  
    return;

  }

  //-------------------------------------------------
  void GausHitFinderAna::analyze(const art::Event& evt)
  {

    // ##############################################
    // ### Outputting Run Number and Event Number ###
    // ##############################################
    std::cout << "run    : " << evt.run() <<" event  : "<<evt.id().event() << std::endl;
    fRun->Fill( evt.run() );
    fEvt->Fill( evt.id().event() );
    //fRun = evt.run();
    //fEvt = evt.id().event();
  
    // #########################################
    // ### Generic Variables used throughout ###
    // #########################################
    int NSinglePulseEvents = 0 , NMultiPulseEvents = 0;
    int SinglePulse = 0, Multipulse = 0;
    double TruthHitTime = 0 , TruthHitCalculated = 0;
  
    // ###############################################
    // ### Hit Variables to be saved to Histograms ###
    // ###############################################
    float GoodnessOfFit_multiplicity1 = 0; 	//<---Goodness of hit with mulitplicity of 1
    float Charge_multiplicity1 = 0;		//<---Charge of hit with multiplicity of 1
    float SigmaCharge_multiplicity1 = 0;		//<---Uncertainty of charge of hit with multiplicity of 1
    float PeakPositionmultiplicity1 = 0;
    float PeakPositionUncertmulitplicity1 = 0;
    float PeakPositionUncertmulitplicityGT1 = 0;
    // ####################################
    // ### Getting Geometry Information ###
    // ####################################
    art::ServiceHandle<geo::Geometry> geom;
  
    // #######################################
    // ### Getting Liquid Argon Properites ###
    // #######################################
    art::ServiceHandle<util::LArProperties> larp;
  
    // ###################################
    // ### Getting Detector Properties ###
    // ###################################
    art::ServiceHandle<util::DetectorProperties> detp;
  
  
    // ##################################################
    // ### Getting the Reconstructed Hits (hitHandle) ###
    // ##################################################
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fGausHitFinderModuleLabel,hitHandle);
  
    // #########################################
    // ### Putting Hits into a vector (hits) ###
    // #########################################
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitHandle);
  
    
    // ####################################
    // ### Using BackTracker HitCheater ###
    // ####################################
    art::ServiceHandle<cheat::BackTracker> bt;
  
  
    // ###############################################################
    // ### Integers used for setting Channel, TPC, Plane, and Wire ###
    // ###############################################################
    unsigned int p = 0, w = 0;
 
    fnhits = hitHandle->size();
  
    // ================================================
    // === Calculating Time Tick and Drift Velocity ===
    // ================================================
    double time_tick      = detp->SamplingRate()/1000.;
    double drift_velocity = larp->DriftVelocity(larp->Efield(),larp->Temperature());
  

    // #########################
    // ### Looping over Hits ###
    // #########################
    for(size_t nHits = 0; nHits < hitHandle->size(); ++nHits){
      
      // === Finding Channel associated with the hit ===
      art::Ptr<recob::Hit> hit(hitHandle, nHits);
      p = hit->WireID().Plane;
      w = hit->WireID().Wire;    	
      
      
      // ===================================================================
      // Using Track IDE's to locate the XYZ location from truth information
      // ===================================================================
      std::vector<cheat::TrackIDE> trackides = bt->HitToTrackID(hit);
      //std::vector<cheat::TrackIDE>::iterator idesitr = trackides.begin();
      std::vector<double> xyz = bt->HitToXYZ(hit);
      
      // ==============================================================
      // Calculating the truth tick position of the hit using 2 methods
      // Method 1: ConvertXtoTicks from the detector properties package
      // Method 2: Actually do the calculation myself to double check things
      // ==============================================================
      
      // ### Method 1 ###
      TruthHitTime = detp->ConvertXToTicks(xyz[0], p, hit->WireID().TPC, hit->WireID().Cryostat);
      
      // ### Method 2 ###
      // ================================================
      // Establishing the x-position of the current plane
      // ================================================ 
      const double origin[3] = {0.};
      double pos[3];
      geom->Plane(p).LocalToWorld(origin, pos);
      double planePos_timeCorr = (pos[0]/drift_velocity)*(1./time_tick)+60; 
      //<---x position of plane / drift velocity + 60 (Trigger offset)
      
      TruthHitCalculated = ( (xyz[0]) / (drift_velocity * time_tick) ) + planePos_timeCorr;
      
      //float TickstoX = detp->ConvertTicksToX(hit->StartTime(),p2,t2,c2);
      
      double hitresid = ( ( TruthHitTime - hit->PeakTime() ) / hit->SigmaPeakTime() );
      fHitResidualAll->Fill( hitresid );
      
      // ##################################################
      // ### Looking at "Hits" with a multiplicity == 1 ###
      // ##################################################
      if(hit->Multiplicity() == 1){
	
	NSinglePulseEvents++;
	SinglePulse = 1;
	
	fSingleHit         = SinglePulse;
	
	// *******************
	// *** Wire Number ***
	// *******************
	//---Filling Histogram---
	fWireNumbMulti1->Fill(w);
      	
	
	// ***********************
	// *** Goodness of Fit ***
	// ***********************
	GoodnessOfFit_multiplicity1	= hit->GoodnessOfFit();
	//---Filling Histogram---
	fFitGoodnessMulti1->Fill( GoodnessOfFit_multiplicity1);
	
	
	// **************
	// *** Charge ***
	// **************
	Charge_multiplicity1		= hit->Charge();
	//---Filling Histogram---
	fChargeMulti1->Fill( Charge_multiplicity1 );
	//FillHisto(Charge_multiplicity1, hname, 2000,0,1000);
	//---Filling TTree---
	fChargen1          = hit->Charge();
	
	// ***********************
	// *** Error on Charge ***
	// ***********************
	SigmaCharge_multiplicity1  = hit->SigmaCharge();
	//---Filling Histogram---
	
	//---Filling TTree---
	fSigmaChargen1     = SigmaCharge_multiplicity1;
	
	// *********************
	// *** Peak Position ***
	// *********************
	PeakPositionmultiplicity1 = hit->PeakTime();
	//---Filling Histogram---
	fRecoPeakPositionMulti1->Fill(PeakPositionmultiplicity1);
	fTruthPeakPosition->Fill(TruthHitTime);
	//---Filling TTree---
	fPeakn1            = PeakPositionmultiplicity1;
	
	// *********************************
	// *** Peak Position Uncertainty ***
	// *********************************
	PeakPositionUncertmulitplicity1 = hit->SigmaPeakTime();
	//---Filling Histogram---
	fRecoPeakPositionUncertMulti1->Fill(PeakPositionUncertmulitplicity1);
	//---Filling TTree---
	fPeakUncertn1      = PeakPositionUncertmulitplicity1;
	
	
	
	fWidthn1           = (hit->EndTime() - hit->PeakTime());
	
      	
      	
	fStartTimen1       = hit->StartTime();
	fStartTimeUncertn1 = hit->SigmaStartTime();
	fEndTimen1         = hit->EndTime();
	fEndTimeUncertn1   = hit->SigmaEndTime();
	
	
	double hitresid = ( ( TruthHitTime - hit->PeakTime() ) / hit->SigmaPeakTime() );
	fHitResidualMulti1->Fill( hitresid );
	
	if(p == 0){
	    
	  //---Filling Histograms---
	  fRecoPeakPositionPlane0Multi1->Fill( hit->PeakTime() );
	  fRecoPeakPositionUncertPlane0Multi1->Fill( hit->SigmaPeakTime() );
	  
	  fTruthPeakPositionPlane0Multi1->Fill( TruthHitCalculated );
	  fHitResidualPlane0Multi1->Fill( hitresid );
	  
	  
	}//<---End looking at plane 0
	
	if(p == 1){
	  //---Filling Histograms---
	  fRecoPeakPositionPlane1Multi1->Fill( hit->PeakTime() );
	  fRecoPeakPositionUncertPlane1Multi1->Fill( hit->SigmaPeakTime() );
	  
	  fTruthPeakPositionPlane1Multi1->Fill(TruthHitCalculated);
	  fHitResidualPlane1Multi1->Fill( hitresid );
	  
	}//<---End looking at Plane 1
	
	
	if(p == 2){
	  //---Filling Histograms---
	  fRecoPeakPositionPlane2Multi1->Fill( hit->PeakTime() );
	  fRecoPeakPositionUncertPlane2Multi1->Fill( hit->SigmaPeakTime() );
	  
	  fTruthPeakPositionPlane2Multi1->Fill(TruthHitCalculated);
	  fHitResidualPlane2Multi1->Fill( hitresid );
	  
	}//<---End looking at Plane 2
	
      }//<---End Hit Multiplicity == 1
      
      // ##################################################
      // ### Looking at "Hits" with a multiplicity == 1 ###
      // ##################################################
      if(hit->Multiplicity() > 1){
	Multipulse = 1;
	NMultiPulseEvents++;
	
	fMultiHit            = Multipulse;
	fWirenGT1            = w;
	fgoodoffitnGT1       = hit->GoodnessOfFit();
	fChargenGT1          = hit->Charge();
	fSigmaChargenGT1     = hit->SigmaCharge();
	fWidthnGT1           = (hit->EndTime() - hit->PeakTime());
	
	fPeaknGT1            = hit->PeakTime();
	fStartTimenGT1       = hit->StartTime();
	fStartTimeUncertnGT1 = hit->SigmaStartTime();
	fEndTimenGT1         = hit->EndTime();
	fEndTimeUncertnGT1   = hit->SigmaEndTime();
	
	// *********************************
	// *** Peak Position Uncertainty ***
	// *********************************
	PeakPositionUncertmulitplicityGT1 = hit->SigmaPeakTime();
	//---Filling Histogram---
	fRecoPeakPositionUncertMultiGT1->Fill(PeakPositionUncertmulitplicityGT1);
	//---Filling TTree---
	fPeakUncertnGT1      = PeakPositionUncertmulitplicityGT1;
	
	
	if(p == 0){
	  //fHitResidualnGT1Plane0     = ( ( TruthHitTime - hit->PeakTime() )/ hit->SigmaPeakTime() );
	  //fHitResidualnGT1Plane0     = ( ( TruthHitCalculated - hit->PeakTime() )/ hit->SigmaPeakTime() );
	  fRecoPeakPositionUncertPlane0MultiGT1->Fill( hit->SigmaPeakTime() );
	}
	
	if(p == 1){
	  //fHitResidualnGT1Plane1     = ( ( TruthHitTime - hit->PeakTime() )/ hit->SigmaPeakTime() );
	  //fHitResidualnGT1Plane1     = ( ( TruthHitCalculated - hit->PeakTime() )/ hit->SigmaPeakTime() );
	  fRecoPeakPositionUncertPlane1MultiGT1->Fill( hit->SigmaPeakTime() );
	}
	
	if(p == 2){
	  //fHitResidualnGT1Plane2     = ( ( TruthHitTime - hit->PeakTime() )/ hit->SigmaPeakTime() );
	  //fHitResidualnGT1Plane2     = ( ( TruthHitCalculated - hit->PeakTime() )/ hit->SigmaPeakTime() );
	  fRecoPeakPositionUncertPlane2MultiGT1->Fill( hit->SigmaPeakTime() );
	}		
	
      }//<---End Hit Multiplicity > 1
      
      
      /*std::cout<<"c = "<<c<<" t = "<<t<<" p = "<<p<<" w = "<<w<<std::endl;
	std::cout<<"Start Time        = "<<	hit->StartTime()	<<	std::endl;
	std::cout<<"Sigma Start Time  = "<<	hit->SigmaStartTime()	<<	std::endl;
	std::cout<<"End Time          = "<<	hit->EndTime()		<<	std::endl;
	std::cout<<"Sigma End Time    = "<<	hit->SigmaEndTime()	<<	std::endl;
	std::cout<<"Peak Time         = "<<	hit->PeakTime()		<<	std::endl;
	std::cout<<"Sigma Peak Time   = "<<	hit->SigmaPeakTime()	<<	std::endl;
	std::cout<<"Multiplicity      = "<<	hit->Multiplicity()	<<	std::endl;
	std::cout<<"Charge            = "<<	hit->Charge()		<<	std::endl;
	std::cout<<"Sigma Charge      = "<<	hit->SigmaCharge()	<<	std::endl;
	std::cout<<"Goodness Fit      = "<<	hit->GoodnessOfFit()	<<	std::endl;
	std::cout<<std::endl;*/
      
      fHTree->Fill();
      Multipulse = 0;
      SinglePulse = 0;
    }//<---End Loop over hits
    
    fnOnePulseHits = NSinglePulseEvents;
    fmulitPulseHits = NMultiPulseEvents;
    fHTree->Fill();
    return;
    
  }//end analyze method
   
  // --------------------------------------------------------
  DEFINE_ART_MODULE(GausHitFinderAna)

} // end of hit namespace



#endif // GAUSHITFINDERANA_H

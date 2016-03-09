//////////////////////////////////////////////////
//
// Calorimetry class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// ART port echurch@fnal.gov
//  This algorithm is designed to perform the calorimetric reconstruction 
//  of the 3D reconstructed tracks
////////////////////////////////////////////////////////////////////////
#ifndef CALO_H
#define CALO_H


extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "lardata/AnalysisAlg/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larcore/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/AnalysisBase/Calorimetry.h"
#include "lardata/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larreco/RecoAlg/TrackUtils.h" // lar::TrackPitchInView()
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>
#include <TVector3.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


///calorimetry
namespace calo {
   
  class Calorimetry : public art::EDProducer {
    
  public:
    
    explicit Calorimetry(fhicl::ParameterSet const& pset); 
    virtual ~Calorimetry();
    
    void beginJob(); 
    //    void endJob();

    void produce(art::Event& evt);

  private:
        
    void   ReadCaloTree();

    bool BeginsOnBoundary(art::Ptr<recob::Track> lar_track);
    bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);

    void GetPitch(art::Ptr<recob::Hit> hit, std::vector<double> trkx, std::vector<double> trky, std::vector<double> trkz, std::vector<double> trkw, std::vector<double> trkx0, double *xyz3d, double &pitch, double TickT0);

    std::string fTrackModuleLabel;
    std::string fSpacePointModuleLabel;
    std::string fT0ModuleLabel;
    bool fMakeTree;
    bool fUseArea;
    CalorimetryAlg caloAlg;
	
	
    TTree *ftree;
    int frun;           //Run 
    int fevent;         //Event
    int fntrack;        //track number
    int ftotTracks;        //track number
    double fTrkPitch;
    double fTrkPitchI;  //
    double fTrkPitchC;
    double fXStart;
    double fYStart;
    double fZStart;
    double fXEnd;
    double fYEnd;
    double fZEnd;

    int fnhits3D; 
    std::vector<double> fXHit;
    std::vector<double> fYHit;
    std::vector<double> fZHit;
    std::vector<double> fMIPs3D;   
 
    int fnhitsIND; 
    int fnspsIND;
    std::vector<int>    fwireIND;
    std::vector<double> ftimeIND;
    std::vector<double> fstimeIND;
    std::vector<double> fetimeIND;
    std::vector<double> fMIPsIND;
    std::vector<double> fdEdxIND;
    std::vector<double> fResRngIND;
	
    int fnhitsCOL;
    int fnspsCOL;
    std::vector<int>    fwireCOL;
    std::vector<double> ftimeCOL;
    std::vector<double> fstimeCOL;
    std::vector<double> fetimeCOL;
    std::vector<double> fMIPsCOL;
    std::vector<double> fdEdxCOL;
    std::vector<double> fResRngCOL;

    int fnhits;
    int fnsps;
    std::vector<int>    fwire;
    std::vector<double> ftime;
    std::vector<double> fstime;
    std::vector<double> fetime;
    std::vector<double> fMIPs;
    std::vector<double> fdQdx;
    std::vector<double> fdEdx;
    std::vector<double> fResRng;
    std::vector<double> fpitch;
    std::vector<TVector3> fXYZ;

    double TPCsize[3];

  protected: 
    
  
  }; // class Calorimetry

}

#endif // CALO_H

//-------------------------------------------------
calo::Calorimetry::Calorimetry(fhicl::ParameterSet const& pset)
  : fTrackModuleLabel(pset.get< std::string >("TrackModuleLabel")      ),
    fSpacePointModuleLabel (pset.get< std::string >("SpacePointModuleLabel")       ),
    fT0ModuleLabel (pset.get< std::string >("T0ModuleLabel") ),
    fMakeTree(pset.get< bool >("MakeTree") ),
    fUseArea(pset.get< bool >("UseArea") ),
    caloAlg(pset.get< fhicl::ParameterSet >("CaloAlg"))
{
  produces< std::vector<anab::Calorimetry>              >();
  produces< art::Assns<recob::Track, anab::Calorimetry> >();
}

//-------------------------------------------------
calo::Calorimetry::~Calorimetry()
{
  
}

//-------------------------------------------------
void calo::Calorimetry::beginJob()
{
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  /// ROOT tree for Calorimetry 
  if(fMakeTree) {
    ftree = tfs->make<TTree>("CaloTree","CaloTree");
    ftree->Branch("run",       &frun,       "run/I");
    ftree->Branch("event",     &fevent,     "event/I");
    ftree->Branch("itrack",    &fntrack,    "itrack/I");
    ftree->Branch("ntracks",   &ftotTracks, "ntracks/I");
    ftree->Branch("TrkPitchI", &fTrkPitchI, "TrkPitchI/F");
    ftree->Branch("TrkPitchC", &fTrkPitchC, "TrkPitchC/F");
    ftree->Branch("XStart",    &fXStart,    "XStart/F");
    ftree->Branch("YStart",    &fYStart,    "YStart/F");
    ftree->Branch("ZStart",    &fZStart,    "ZStart/F");  
    ftree->Branch("XEnd",      &fXEnd,      "XEnd/F");  
    ftree->Branch("YEnd",      &fYEnd,      "YEnd/F");  
    ftree->Branch("ZEnd",      &fZEnd,      "ZEnd/F");  
    ftree->Branch("nhits3D",   &fnhits3D,   "nhits3D/I"); 
    ftree->Branch("XHit",      &fXHit);		     
    ftree->Branch("YHit",      &fYHit);		     
    ftree->Branch("ZHit",      &fZHit);                 
    ftree->Branch("nhitsIND",  &fnhitsIND, "nhitsIND/I");	 
    ftree->Branch("nspsIND",   &fnspsIND,  "nspsIND/I");	 
    ftree->Branch("wireIND",   &fwireIND);		 
    ftree->Branch("timeIND",   &ftimeIND);		 
    ftree->Branch("stimeIND",  &fstimeIND);		 
    ftree->Branch("etimeIND",  &fetimeIND);		 
    ftree->Branch("MIPsIND",   &fMIPsIND);  
    ftree->Branch("dEdxIND",   &fdEdxIND);
    ftree->Branch("ResRngIND", &fResRngIND);
    ftree->Branch("nhitsCOL",  &fnhitsCOL, "nhitsCOL/I");
    ftree->Branch("nspsCOL",   &fnspsCOL,  "nspsCOL/I");
    ftree->Branch("wireCOL",   &fwireCOL);
    ftree->Branch("timeCOL",   &ftimeCOL);
    ftree->Branch("stimeCOL",  &fstimeCOL);
    ftree->Branch("etimeCOL",  &fetimeCOL);
    ftree->Branch("MIPsCOL",   &fMIPsCOL);
    ftree->Branch("dEdxCOL",   &fdEdxCOL);
    ftree->Branch("ResRngCOL", &fResRngCOL);
  }
   
  return;
}

//------------------------------------------------------------------------------------//
void calo::Calorimetry::produce(art::Event& evt)
{ 
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  // Get Geometry
  art::ServiceHandle<geo::Geometry> geom;
  
  //TPC dimensions  
  TPCsize[0] = (geom->DetHalfWidth())*2.;
  TPCsize[1] = (geom->DetHalfHeight())*2.;
  TPCsize[2] = (geom->DetLength());

  // channel quality
  lariov::ChannelStatusProvider const& channelStatus
    = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

  frun = evt.id().run();
  fevent = evt.id().event();

  fntrack = 0; 
  fnhits3D = 0; fnhitsIND = 0; fnhitsCOL = 0;

  size_t nplanes = geom->Nplanes();


  //create anab::Calorimetry objects and make association with recob::Track
  std::unique_ptr< std::vector<anab::Calorimetry> > calorimetrycol(new std::vector<anab::Calorimetry>);
  std::unique_ptr< art::Assns<recob::Track, anab::Calorimetry> > assn(new art::Assns<recob::Track, anab::Calorimetry>);

  //art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit>        fmht(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<anab::T0>          fmt0(trackListHandle, evt, fT0ModuleLabel);

  for(size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter){   
    //std::cout<<"!!!!!!NEW TRACK "<< int(trkIter) <<"!!!!!!!"<<std::endl;
    //fnhits3D = fmsp.at(trkIter).size();
    std::vector<double> larStart;
    std::vector<double> larEnd;
    //put xyz coordinates at begin/end of track into vectors(?)
    tracklist[trkIter]->Extent(larStart,larEnd);
    //store track directional cosines
    double trackCosStart[3]={0.,0.,0.};
    double trackCosEnd[3]={0.,0.,0.};
    tracklist[trkIter]->Direction(trackCosStart,trackCosEnd);
        
    bool startsonboundary = BeginsOnBoundary(tracklist[trkIter]);
    bool endsonboundary   = EndsOnBoundary(tracklist[trkIter]);

    int containment = 0;
    if(!startsonboundary && !endsonboundary ) containment = 0;  // contained track
    if(!startsonboundary &&  endsonboundary ) containment = 1;  // escaping track
    if( startsonboundary && !endsonboundary ) containment = 2;  // entering track
    if( startsonboundary &&  endsonboundary ) containment = 3;  // passing track 

    bool TrackStops = (containment == 0 || containment == 2);
    //      fMIPs3D = new double[fnhits3D];

    //recover the Induction (geo::kU) hits list
    
    //int npI = 0;
    
    // Some variables for the hit
    float time;          //hit time at maximum
    float stime;         //hit start time 
    float etime;         //hit end time 
    uint32_t     channel = 0;//channel number
    unsigned int cstat   = 0;    //hit cryostat number 
    unsigned int tpc     = 0;    //hit tpc number 
    unsigned int wire    = 0;   //hit wire number 
    unsigned int plane   = 0;  //hit plane number

    

    std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(trkIter);
    double T0 =0;
    double TickT0 =0;
    if ( fmt0.isValid() ) {
      std::vector< art::Ptr<anab::T0> > allT0 = fmt0.at(trkIter);
      if ( allT0.size() ) T0 = allT0[0]->Time();
      TickT0 = T0 / detprop->SamplingRate();    
    }
    
    std::vector< std::vector<unsigned int> > hits(nplanes);

    art::FindManyP<recob::SpacePoint> fmspts(allHits, evt, fSpacePointModuleLabel);
    for (size_t ah = 0; ah< allHits.size(); ++ah){
      hits[allHits[ah]->WireID().Plane].push_back(ah);
    }
    //get hits in each plane
    for (size_t ipl = 0; ipl < nplanes; ++ipl){//loop over all wire planes

      geo::PlaneID planeID;//(cstat,tpc,ipl);

      fwire.clear();
      ftime.clear();
      fstime.clear();
      fetime.clear();
      fMIPs.clear();
      fdQdx.clear();
      fdEdx.clear();
      fpitch.clear();
      fResRng.clear();
      fXYZ.clear();

      double Kin_En = 0.;
      double Trk_Length = 0.;
      std::vector<double> vdEdx;
      std::vector<double> vresRange;
      std::vector<double> vdQdx;
      std::vector<double> deadwire; //residual range for dead wires
      std::vector<TVector3> vXYZ;

      //range of wire signals
      unsigned int wire0 = 100000;
      unsigned int wire1 = 0;
      double PIDA = 0;
      int nPIDA = 0;

      // determine track direction. Fill residual range array
      bool GoingDS = true;
      // find the track direction by comparing US and DS charge BB
      double USChg = 0;
      double DSChg = 0;
      // temp array holding distance betweeen space points
      std::vector<double> spdelta;
      //int nht = 0; //number of hits
      fnsps = 0; //number of space points
      std::vector<double> ChargeBeg;
      std::stack<double> ChargeEnd;     

      // find the separation between all space points
      double xx = 0.,yy = 0.,zz = 0.;

      try{
	fTrkPitch = lar::TrackPitchInView(*tracklist[trkIter], geom->Plane(ipl).View());
      }
      catch( cet::exception &e){
	mf::LogWarning("Calorimetry") << "caught exception " 
				      << e << "\n setting pitch (C) to "
				      << util::kBogusD;
	fTrkPitch = 0;
      }

      //save track 3d points
      std::vector<double> trkx;
      std::vector<double> trky;
      std::vector<double> trkz;
      std::vector<double> trkw;
      std::vector<double> trkx0;
      for (size_t i = 0; i<hits[ipl].size(); ++i){
	//Get space points associated with the hit
	std::vector< art::Ptr<recob::SpacePoint> > sptv = fmspts.at(hits[ipl][i]);
	for (size_t j = 0; j < sptv.size(); ++j){
	  
	  double t = allHits[hits[ipl][i]]->PeakTime() - TickT0; // Want T0 here? Otherwise ticks to x is wrong?
	  double x = detprop->ConvertTicksToX(t, allHits[hits[ipl][i]]->WireID().Plane, allHits[hits[ipl][i]]->WireID().TPC, allHits[hits[ipl][i]]->WireID().Cryostat);
	  double w = allHits[hits[ipl][i]]->WireID().Wire;
	  trkx.push_back(sptv[j]->XYZ()[0]-detprop->ConvertTicksToX(TickT0, allHits[hits[ipl][i]]->WireID().Plane, allHits[hits[ipl][i]]->WireID().TPC, allHits[hits[ipl][i]]->WireID().Cryostat));
	  trky.push_back(sptv[j]->XYZ()[1]);
	  trkz.push_back(sptv[j]->XYZ()[2]);
	  trkw.push_back(w);
	  trkx0.push_back(x);
	}
      }
      for (size_t ihit = 0; ihit < hits[ipl].size(); ++ihit){//loop over all hits on each wire plane

	//std::cout<<ihit<<std::endl;
	
	if (!planeID.isValid){
	  plane = allHits[hits[ipl][ihit]]->WireID().Plane;
	  tpc   = allHits[hits[ipl][ihit]]->WireID().TPC;
	  cstat = allHits[hits[ipl][ihit]]->WireID().Cryostat;
	  planeID.Cryostat = cstat;
	  planeID.TPC = tpc;
	  planeID.Plane = plane;
	  planeID.isValid = true;
	}

	wire = allHits[hits[ipl][ihit]]->WireID().Wire;
	time = allHits[hits[ipl][ihit]]->PeakTime(); // What about here? T0 
	stime = allHits[hits[ipl][ihit]]->PeakTimeMinusRMS();
	etime = allHits[hits[ipl][ihit]]->PeakTimePlusRMS();
	
	double charge = allHits[hits[ipl][ihit]]->PeakAmplitude();
	if (fUseArea) charge = allHits[hits[ipl][ihit]]->Integral();
	//get 3d coordinate and track pitch for the current hit
	//not all hits are associated with space points, the method uses neighboring spacepts to interpolate
	double xyz3d[3];
	double pitch;
	GetPitch(allHits[hits[ipl][ihit]], trkx, trky, trkz, trkw, trkx0, xyz3d, pitch, TickT0);

	if (xyz3d[2]<-100) continue; //hit not on track
	if (pitch<=0) pitch = fTrkPitch;
	if (!pitch) continue;

        if(fnsps == 0) {
          xx = xyz3d[0];
          yy = xyz3d[1];
          zz = xyz3d[2];
          spdelta.push_back(0);
        } else {
          double dx = xyz3d[0] - xx;
          double dy = xyz3d[1] - yy;
          double dz = xyz3d[2] - zz;
          spdelta.push_back(sqrt(dx*dx + dy*dy + dz*dz));
          Trk_Length += spdelta.back();
          xx = xyz3d[0];
          yy = xyz3d[1];
          zz = xyz3d[2];
        }
	
	ChargeBeg.push_back(charge);
	ChargeEnd.push(charge);

	double MIPs = charge;
	double dQdx = MIPs/pitch;
	double dEdx = 0;
	if (fUseArea) dEdx = caloAlg.dEdx_AREA(allHits[hits[ipl][ihit]], pitch, T0);
	else dEdx = caloAlg.dEdx_AMP(allHits[hits[ipl][ihit]], pitch, T0);

	Kin_En = Kin_En + dEdx * pitch;	

	if (allHits[hits[ipl][ihit]]->WireID().Wire < wire0) wire0 = allHits[hits[ipl][ihit]]->WireID().Wire;
	if (allHits[hits[ipl][ihit]]->WireID().Wire > wire1) wire1 = allHits[hits[ipl][ihit]]->WireID().Wire;

	fMIPs.push_back(MIPs);
	fdEdx.push_back(dEdx);
	fdQdx.push_back(dQdx);
	fwire.push_back(wire);
	ftime.push_back(time);
	fstime.push_back(stime);
	fetime.push_back(etime);
	fpitch.push_back(pitch);
	TVector3 v(xyz3d[0],xyz3d[1],xyz3d[2]);
	//std::cout << "Adding these positions to v and then fXYZ " << xyz3d[0] << " " << xyz3d[1] << " " << xyz3d[2] << "\n" <<std::endl;
	fXYZ.push_back(v);
	++fnsps;
      }
      if (!fnsps){
	//std::cout << "Adding the aforementioned positions..." << std::endl;
	calorimetrycol->push_back(anab::Calorimetry(util::kBogusD,
						    vdEdx,
						    vdQdx,
						    vresRange,
						    deadwire,
						    util::kBogusD,
						    fpitch,
						    vXYZ,
						    planeID));
	util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);
	continue;
      }
      for (int isp = 0; isp<fnsps; ++isp){
	if (isp>3) break;
	USChg += ChargeBeg[isp];
      }
      int countsp = 0;
      while (!ChargeEnd.empty()){
	if (countsp>3) break;
	DSChg += ChargeEnd.top();
	ChargeEnd.pop();
	++countsp;
      }
      // Going DS if charge is higher at the end
      GoingDS = (DSChg > USChg);
      // determine the starting residual range and fill the array
      fResRng.resize(fnsps);
      if(GoingDS) {
        fResRng[fnsps - 1] = spdelta[fnsps - 1] / 2;
        for(int isp = fnsps - 2; isp > -1; isp--) {
          fResRng[isp] = fResRng[isp+1] + spdelta[isp+1];
        }
      } else {
        fResRng[0] = spdelta[1] / 2;
        for(int isp = 1; isp < fnsps; isp++) {
          fResRng[isp] = fResRng[isp-1] + spdelta[isp];
        }
      }
    
      LOG_DEBUG("CaloPrtTrk") << "Calorimetry Run/Evt: "
		 << frun   << " / " 
		 << fevent << " Track #"
		 << trkIter << " Plane #"
		 << ipl;
      switch (containment) {
      case 0:
	LOG_DEBUG("CaloPrtTrk") << ", contained.";
	break;
      case 1:
	LOG_DEBUG("CaloPrtTrk") << ", escaping.";
	break;
      case 2:
	LOG_DEBUG("CaloPrtTrk") << ", entering.";
	break;
      case 3:
	LOG_DEBUG("CaloPrtTrk") << ", passing.";
	break;
      default:
	LOG_DEBUG("CaloPrtTrk") << ", ??";
      }
      if(TrackStops) {
	if(GoingDS) {
	  LOG_DEBUG("CaloPrtTrk") <<" Going downstream";
	} 
	else {
	  LOG_DEBUG("CaloPrtTrk") <<" Going upstream";
	}
      }
      LOG_DEBUG("CaloPrtTrk")<<"\n";
      LOG_DEBUG("CaloPrtHit") << " pt wire  time  ResRng    MIPs   pitch   dE/dx    Ai X Y Z\n";

      if (ipl==0){
	fnspsIND = fnsps;
	fTrkPitchI = fTrkPitch;
      }
      else if (ipl == nplanes-1){
	fnspsCOL = fnsps;
	fTrkPitchC = fTrkPitch;
      }
      double Ai = -1;
      for (int i = 0; i < fnsps; ++i){//loop over all 3D points
	if (ipl==0){
	  fwireIND.push_back(fwire[i]);
	  ftimeIND.push_back(ftime[i]);
	  fstimeIND.push_back(fstime[i]);
	  fetimeIND.push_back(fetime[i]);
	  fMIPsIND.push_back(fMIPs[i]);
	  fdEdxIND.push_back(fdEdx[i]);
	  fResRngIND.push_back(fResRng[i]);
	}
	else if (ipl == nplanes-1){
	  fwireCOL.push_back(fwire[i]);
	  ftimeCOL.push_back(ftime[i]);
	  fstimeCOL.push_back(fstime[i]);
	  fetimeCOL.push_back(fetime[i]);
	  fMIPsCOL.push_back(fMIPs[i]);
	  fdEdxCOL.push_back(fdEdx[i]);
	  fResRngCOL.push_back(fResRng[i]);
	}
	if (i!=0 && i!= fnsps-1){//ignore the first and last point
	  vresRange.push_back(fResRng[i]);
	  vdEdx.push_back(fdEdx[i]);
	  vdQdx.push_back(fdQdx[i]);
	  vXYZ.push_back(fXYZ[i]);
	  // Calculate PIDA 
	  if(TrackStops){
	    Ai = fdEdx[i] * pow(fResRng[i],0.42);
	    nPIDA++;
	    PIDA += Ai;
	  }
	}
	LOG_DEBUG("CaloPrtHit") << std::setw(4) << i
	  //std::cout << std::setw(4) << i
		   <<std::setw(4)  << fwire[i]
		   << std::setw(6) << (int)ftime[i]
		   << std::setiosflags(std::ios::fixed | std::ios::showpoint)
		   << std::setprecision(2)
		   << std::setw(8) << fResRng[i]
		   << std::setprecision(1)
		   << std::setw(8) << fMIPs[i]
		   << std::setprecision(2)
		   << std::setw(8) << fpitch[i]
		   << std::setw(8) << fdEdx[i]
		   << std::setw(8) << Ai
		  << std::setw(8) << fXYZ[i].x()
		  << std::setw(8) << fXYZ[i].y()
		  << std::setw(8) << fXYZ[i].z()
		   << "\n";
      }//end looping over 3D points
      if(nPIDA > 0) {
	PIDA = PIDA / (double)nPIDA;
      } 
      else {
	PIDA = -1;
      }
      LOG_DEBUG("CaloPrtTrk") << "Plane # "<< ipl
		 << "TrkPitch= "
		 << std::setprecision(2) << fTrkPitch 
		 << " nhits= "        << fnsps
		 << "\n" 
		 << std::setiosflags(std::ios::fixed | std::ios::showpoint)
		 << "Trk Length= "       << std::setprecision(1)
		 << Trk_Length           << " cm,"
		 << " KE calo= "         << std::setprecision(1)
		 << Kin_En               << " MeV,"
		 << " PIDA= "            << PIDA
		 << "\n";
      
      // look for dead wires
      for (unsigned int iw = wire0; iw<wire1+1; ++iw){
	plane = allHits[hits[ipl][0]]->WireID().Plane;
	tpc   = allHits[hits[ipl][0]]->WireID().TPC;
	cstat = allHits[hits[ipl][0]]->WireID().Cryostat;
	channel = geom->PlaneWireToChannel(plane,iw,tpc,cstat);
	if (channelStatus.IsBad(channel)){
	  LOG_DEBUG("Calorimetry") << "Found dead wire at Plane = " << plane 
					 << " Wire =" << iw;
	  unsigned int closestwire = 0;
	  unsigned int endwire = 0;
	  int dwire = 100000;
	  double mindis = 100000;
	  double goodresrange = 0;
	  //hitCtr = 0;
	  for (size_t ihit = 0; ihit <hits[ipl].size(); ++ihit){
	    //	for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitsV.begin(); 
	    //	    hitIter != hitsV.end();  
	    //	    ++hitCtr, hitIter++){
	    channel = allHits[hits[ipl][ihit]]->Channel();
	    if (channelStatus.IsBad(channel)) continue;
	    // grab the space points associated with this hit
	    std::vector< art::Ptr<recob::SpacePoint> > sppv = fmspts.at(hits[ipl][ihit]);
	    if(sppv.size() < 1) continue;
	    // only use the first space point in the collection, really each hit should
	    // only map to 1 space point
	    const double xyz[3] = {sppv[0]->XYZ()[0],
				   sppv[0]->XYZ()[1],
				   sppv[0]->XYZ()[2]};
	    double dis1 = (larEnd[0]-xyz[0])*(larEnd[0]-xyz[0]) +
	      (larEnd[1]-xyz[1])*(larEnd[1]-xyz[1]) +
	      (larEnd[2]-xyz[2])*(larEnd[2]-xyz[2]);
	    if (dis1) dis1 = std::sqrt(dis1);
	    if (dis1 < mindis){
	      endwire = allHits[hits[ipl][ihit]]->WireID().Wire;
	      mindis = dis1;
	    }
	    if (std::abs(wire-iw) < dwire){
	      closestwire = allHits[hits[ipl][ihit]]->WireID().Wire;
	      dwire = abs(allHits[hits[ipl][ihit]]->WireID().Wire-iw);
	      goodresrange = dis1;
	    }
	  }
	  if (closestwire){
	    if (iw < endwire){
	      deadwire.push_back(goodresrange+(int(closestwire)-int(iw))*fTrkPitch);
	    }
	    else{
	      deadwire.push_back(goodresrange+(int(iw)-int(closestwire))*fTrkPitch);
	    }
	  }
	}
      }
      //std::cout << "Adding at the end but still same fXYZ" << std::endl;
      calorimetrycol->push_back(anab::Calorimetry(Kin_En,
						  vdEdx,
						  vdQdx,
						  vresRange,
						  deadwire,
						  Trk_Length,
						  fpitch,
						  vXYZ,
						  planeID));
      util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);
      
    }//end looping over planes
  }//end looping over tracks
  
  evt.put(std::move(calorimetrycol));
  evt.put(std::move(assn));

  return;
}

//--------------------------------------------------
bool calo::Calorimetry::BeginsOnBoundary(art::Ptr<recob::Track> lar_track)
{
  double fdBoundary = 1.5;
  std::vector<double> larStart, larEnd;
  //put xyz coordinates at begin/end of track into vectors(?)
  lar_track->Extent(larStart,larEnd);

  if(std::abs(larStart[0])            < fdBoundary ||
     std::abs(TPCsize[0]-larStart[0]) < fdBoundary ||
     std::abs(larStart[1]+TPCsize[1]/2) < fdBoundary ||
     std::abs(TPCsize[1]/2-larStart[1]) < fdBoundary ||
     std::abs(larStart[2])            < fdBoundary ||
     std::abs(TPCsize[2]-larStart[2]) < fdBoundary )   
    return true;  
  else return false;
}
    
//--------------------------------------------------
bool calo::Calorimetry::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
{
  double fdBoundary = 1.5;
  std::vector<double> larStart, larEnd;
  //put xyz coordinates at begin/end of track into vectors(?)
  lar_track->Extent(larStart, larEnd); 

  if(std::abs(larEnd[0])              < fdBoundary ||
     std::abs(larEnd[2])              < fdBoundary ||
     std::abs(TPCsize[0]-larEnd[0])   < fdBoundary ||
     std::abs(TPCsize[2]-larEnd[2])   < fdBoundary ||   
     std::abs(larEnd[1]+TPCsize[1]/2) < fdBoundary ||
     std::abs(TPCsize[1]/2-larEnd[1]) < fdBoundary )
    return true;  

  return false;
}

//------------------------------------------------------------------------------------//
void calo::Calorimetry::ReadCaloTree(){
  // This method is an example of how to read 
  // the dynamic vectors in the ROOT Tree

  ftree->Scan("run:event:TrkPitchI:TrkPitchC:XStart:nhits3D:nhitsIND:nhitsCOL");
  int nentries=(int)ftree->GetEntries();
  LOG_DEBUG("Calorimetry") << "nentries  " << nentries;

  std::vector<double> *MIPsCOL = 0;
  TBranch *bMIPsCOL = 0;
  ftree->SetBranchAddress("MIPsCOL",&MIPsCOL,&bMIPsCOL);

  std::vector<double> *MIPsIND = 0;
  TBranch *bMIPsIND = 0;
  ftree->SetBranchAddress("MIPsIND",&MIPsIND,&bMIPsIND);
  for (int i = 0; i < nentries; ++i) {
    LOG_DEBUG("Calorimetry") << " entry " << i;
    long int tentry = ftree->LoadTree(i);

    bMIPsCOL->GetEntry(tentry);
    LOG_DEBUG("Calorimetry") << "# of hits COLL " << MIPsCOL->size();
    for (size_t j = 0; j < MIPsCOL->size(); ++j) {
      LOG_DEBUG("Calorimetry") << " Coll " << MIPsCOL->at(j);
    }
    bMIPsIND->GetEntry(tentry);
    LOG_DEBUG("Calorimetry") << "# of hits IND " << MIPsIND->size();
    for (size_t j = 0; j < MIPsIND->size(); ++j) {
      LOG_DEBUG("Calorimetry") << " Ind " << MIPsIND->at(j);
    }
  }
}

void calo::Calorimetry::GetPitch(art::Ptr<recob::Hit> hit, std::vector<double> trkx, std::vector<double> trky, std::vector<double> trkz, std::vector<double> trkw, std::vector<double> trkx0, double *xyz3d, double &pitch, double TickT0){
  //Get 3d coordinates and track pitch for each hit
  //Find 5 nearest space points and determine xyz and curvature->track pitch
  
  //std::cout << "Start of get pitch" << std::endl;

  // Get services
  art::ServiceHandle<geo::Geometry> geom;
  auto const* dp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  
  //save distance to each spacepoint sorted by distance
  std::map<double,size_t> sptmap;
  //save the sign of distance
  std::map<size_t, int> sptsignmap;

  double wire_pitch = geom->WirePitch(0,1,0);

  double t0 = hit->PeakTime() - TickT0;
  double x0 = dp->ConvertTicksToX(t0, hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat);
  double w0 = hit->WireID().Wire;

  for (size_t i = 0; i<trkx.size(); ++i){
    double distance = pow((trkw[i]-w0)*wire_pitch,2)+pow(trkx0[i]-x0,2);
    if (distance>0) distance = sqrt(distance);
    //std::cout << "Dis " << distance << ", sqaured " << distance*distance << " = " << wire_pitch*wire_pitch <<"("<<trkw[i]<<"-"<<w0<<")^2 + ("<<trkx0[i]<<"-"<<x0<<")^2"<<std::endl;
    sptmap.insert(std::pair<double,size_t>(distance,i));
    if (w0-trkw[i]>0) sptsignmap.insert(std::pair<size_t,int>(i,1));
    else sptsignmap.insert(std::pair<size_t,int>(i,-1));
  }

  //x,y,z vs distance
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;
  std::vector<double> vs;

  double kx = 0, ky = 0, kz = 0;

  int np = 0;
  for (auto isp = sptmap.begin(); isp!=sptmap.end(); isp++){
//    const double *xyz = new double[3];
//    xyz = isp->second->XYZ();
    double xyz[3];
    xyz[0] = trkx[isp->second];
    xyz[1] = trky[isp->second];
    xyz[2] = trkz[isp->second];
        
    double distancesign = sptsignmap[isp->second];
    //std::cout<<np<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<(*isp).first<<std::endl;
    if (np==0&&isp->first>30){//hit not on track
      xyz3d[0] = -1000;
      xyz3d[1] = -1000;
      xyz3d[2] = -1000;
      pitch = -1;
      return;
    }
    //std::cout<<np<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<(*isp).first<<" Plane " << hit->WireID().Plane << " TPC " << hit->WireID().TPC << std::endl;
    if (np<5) {
      vx.push_back(xyz[0]);
      vy.push_back(xyz[1]);
      vz.push_back(xyz[2]);
      vs.push_back(isp->first*distancesign);
    }
    else {
      break;
    }
    np++;
    //delete [] xyz;
  }
  //std::cout<<"np="<<np<<std::endl;
  if (np>=2){//at least two points
    //std::cout << "At least 2 points.."<<std::endl;
    TGraph *xs = new TGraph(np,&vs[0],&vx[0]);
    //for (int i = 0; i<np; i++) std::cout<<i<<" "<<vs[i]<<" "<<vx[i]<<" "<<vy[i]<<" "<<vz[i]<<std::endl;
    try{
      if (np>2){
	xs->Fit("pol2","Q");
      }
      else{
	xs->Fit("pol1","Q");
      }
      TF1 *pol = 0;
      if (np>2) pol = (TF1*) xs->GetFunction("pol2");
      else pol = (TF1*) xs->GetFunction("pol1");
      xyz3d[0] = pol->Eval(0);
      kx = pol->GetParameter(1);
      //std::cout<<"X fit "<<xyz3d[0]<<" "<<kx<<std::endl;
    }
    catch(...){
      mf::LogWarning("Calorimetry::GetPitch") <<"Fitter failed";
      xyz3d[0] = vx[0];
    }
    delete xs;
    TGraph *ys = new TGraph(np,&vs[0],&vy[0]);
    try{
      if (np>2){
	ys->Fit("pol2","Q");
      }
      else{
	ys->Fit("pol1","Q");
      }
      TF1 *pol = 0;
      if (np>2) pol = (TF1*) ys->GetFunction("pol2");
      else pol = (TF1*) ys->GetFunction("pol1");
      xyz3d[1] = pol->Eval(0);
      ky = pol->GetParameter(1);
      //std::cout<<"Y fit "<<xyz3d[1]<<" "<<ky<<std::endl;
    }
    catch(...){
      mf::LogWarning("Calorimetry::GetPitch") <<"Fitter failed";
      xyz3d[1] = vy[0];
    }
    delete ys;
    TGraph *zs = new TGraph(np,&vs[0],&vz[0]);
    try{
      if (np>2){
	zs->Fit("pol2","Q");
      }
      else{
	zs->Fit("pol1","Q");
      }
      TF1 *pol = 0;
      if (np>2) pol = (TF1*) zs->GetFunction("pol2");
      else pol = (TF1*) zs->GetFunction("pol1");
      xyz3d[2] = pol->Eval(0);
      kz = pol->GetParameter(1);
      //std::cout<<"Z fit "<<xyz3d[2]<<" "<<kz<<std::endl;
    }
    catch(...){
      mf::LogWarning("Calorimetry::GetPitch") <<"Fitter failed";
      xyz3d[2] = vz[0];
    }
    delete zs;
  }
  else if (np){
    xyz3d[0] = vx[0];
    xyz3d[1] = vy[0];
    xyz3d[2] = vz[0];
  }
  else{
    xyz3d[0] = -1000;
    xyz3d[1] = -1000;
    xyz3d[2] = -1000;
    pitch = -1;
    return;
  }
  pitch = -1;
  if (kx*kx+ky*ky+kz*kz){
    double tot = sqrt(kx*kx+ky*ky+kz*kz);
    kx /= tot;
    ky /= tot;
    kz /= tot;
    //get pitch
    double wirePitch = geom->WirePitch(0,1,hit->WireID().Plane,hit->WireID().TPC,hit->WireID().Cryostat);
    double angleToVert = geom->Plane(hit->WireID().Plane,hit->WireID().TPC,hit->WireID().Cryostat).Wire(0).ThetaZ(false) - 0.5*TMath::Pi();
    double cosgamma = TMath::Abs(TMath::Sin(angleToVert)*ky+TMath::Cos(angleToVert)*kz);
    if (cosgamma>0) pitch = wirePitch/cosgamma;   

  }
  //std::cout << "At end of get pitch " << xyz3d[0] << " " << xyz3d[1] << " " << xyz3d[2] << " " << x0 << " " << std::endl;
}


namespace calo{

  DEFINE_ART_MODULE(Calorimetry)
  
} // end namespace 


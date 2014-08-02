////////////////////////////////////////////////////////////////////////
//
//  CosmicTracker
//
//  Tracker to reconstruct cosmic ray muons
// 
//  tjyang@fnal.gov
// 
//  This algorithm is based on orignal idea in Track3Dreco
//
//  ** Modified by Muhammad Elnimr to check multiple TPCs for the 35t prototype.
//   
//  mmelnimr@as.ua.edu
//  April 2014
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
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

// ROOT includes
#include "TVectorD.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1D.h"
#include "TVirtualFitter.h"

//2-D weighted fit of time vs wire
std::vector<double> vwire;
std::vector<double> vtime;
std::vector<double> vph;

struct CluLen{
  int index;
  double length;
};

bool myfunction (CluLen c1, CluLen c2) { return (c1.length>c2.length);}

void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  //minimisation function computing the sum of squares of residuals
  f = 0;
  for (size_t i = 0; i<vwire.size(); ++i){
    double x = vwire[i];
    double y = par[0]+par[1]*x+par[2]*x*x;
    //using ph^2 to suppress low ph hits
    f += vph[i]*vph[i]*(vtime[i]-y)*(vtime[i]-y);
  }
}

bool sp_sort_x0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[0] < xyz2[0];
}

bool sp_sort_x1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[0] > xyz2[0];
}

bool sp_sort_y0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[1] < xyz2[1];
}

bool sp_sort_y1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[1] > xyz2[1];
}

bool sp_sort_z0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] < xyz2[2];
}

bool sp_sort_z1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] > xyz2[2];
}

namespace trkf {
   
  struct SortByWire {
    bool operator() (art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2) const { 
      return 
	h1->Wire()->RawDigit()->Channel() < 
	h2->Wire()->RawDigit()->Channel() ;
    }
  };

  class CosmicTracker : public art::EDProducer {
    
  public:
    
    explicit CosmicTracker(fhicl::ParameterSet const& pset);
    ~CosmicTracker();
    
    //////////////////////////////////////////////////////////
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:

    double          fKScut;              ///< tolerance for cluster matching based on KS test.

    double          ftmatch;             ///< tolerance for time matching (in time samples) 
    
    double          fsmatch;             ///< tolerance for distance matching (in cm)

    double          ftoler1;             ///< tolerance for using hits in fit

    double          ftoler2;             ///< tolerance for using hits in fit

    std::string     fClusterModuleLabel; ///< label for input cluster collection
    
    bool             fdebug;

    bool            fEnableU;
    bool            fEnableV;
    bool            fEnableZ;

    int             fisohitcut;
    
    std::string     fSortDir;            ///< sort space points 

    bool            fCleanUpHits;        ///< flag to remove outlier hits

    bool            fDirSPS;             ///< calculate direction cosine for each space point
    //testing histograms
    std::vector<TH1D *> dtime;
    std::vector<TH1D *> testsig;
    std::vector<TH1D *> testpulse;
    TH1D *hks;
  
  }; // class CosmicTracker

}

namespace trkf {

  //-------------------------------------------------
  CosmicTracker::CosmicTracker(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<recob::Track>                        >();
    produces< std::vector<recob::SpacePoint>                   >();
    produces< art::Assns<recob::Track,      recob::Cluster>    >();
    produces< art::Assns<recob::Track,      recob::SpacePoint> >();
    produces< art::Assns<recob::SpacePoint, recob::Hit>        >();
    produces< art::Assns<recob::Track,      recob::Hit>        >();

    dtime  .resize(3);
    testsig.resize(3);
    testpulse.resize(3);
  }

  //-------------------------------------------------
  CosmicTracker::~CosmicTracker()
  {
  }

  //-------------------------------------------------
  void CosmicTracker::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
    fKScut                  = pset.get< double >("KScut");
    ftmatch                 = pset.get< double >("TMatch");
    fsmatch                 = pset.get< double >("SMatch");
    ftoler1                 = pset.get< double >("Toler1");
    ftoler2                 = pset.get< double >("Toler2");
    fdebug                  = pset.get< bool   >("Debug");
    fEnableU                = pset.get< bool   >("EnableU");
    fEnableV                = pset.get< bool   >("EnableV");
    fEnableZ                = pset.get< bool   >("EnableZ");
    fisohitcut              = pset.get< int    >("IsoHitCut");
    fSortDir                = pset.get< std::string >("SortDirection","+z");
    fCleanUpHits            = pset.get< bool   >("CleanUpHits");
    fDirSPS                 = pset.get< bool   >("DirSPS");
  }

  //-------------------------------------------------
  void CosmicTracker::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
 
    dtime[0] = tfs->make<TH1D>("dtime0","dtime0",100,-50,50);
    dtime[1] = tfs->make<TH1D>("dtime1","dtime1",100,-50,50);
    dtime[2] = tfs->make<TH1D>("dtime2","dtime2",100,-50,50);

    testsig[0] = tfs->make<TH1D>("testsig0","testsig0",4096,0,4096);
    testsig[1] = tfs->make<TH1D>("testsig1","testsig1",4096,0,4096);
    testsig[2] = tfs->make<TH1D>("testsig2","testsig2",4096,0,4096);

    testpulse[0] = tfs->make<TH1D>("testpulse0","testpulse0",4096,0,4096);
    testpulse[1] = tfs->make<TH1D>("testpulse1","testpulse1",4096,0,4096);
    testpulse[2] = tfs->make<TH1D>("testpulse2","testpulse2",4096,0,4096);


    hks = tfs->make<TH1D>("hks","hks",100,0,1);

    for (int i = 0; i<3; ++i) dtime[i]->Sumw2();


  }

  //-------------------------------------------------
  void CosmicTracker::endJob()
  {
  }

  //------------------------------------------------------------------------------------//
  void CosmicTracker::produce(art::Event& evt){
  
    // get services
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    std::unique_ptr<std::vector<recob::Track>      >              tcol (new std::vector<recob::Track>);	   
    std::unique_ptr<std::vector<recob::SpacePoint> > 	        spcol(new std::vector<recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::Cluster> >    tcassn(new art::Assns<recob::Track, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit> >        thassn(new art::Assns<recob::Track, recob::Hit>);
    std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> >   shassn(new art::Assns<recob::SpacePoint, recob::Hit>);

    double timetick = detprop->SamplingRate()*1e-3;    //time sample in us
    //double presamplings = detprop->TriggerOffset(); // presamplings in ticks  
    //double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
    double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
    double Efield_drift = larprop->Efield(0);  // Electric Field in the drift region in kV/cm
    double Temperature = larprop->Temperature();  // LAr Temperature in K

    double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
    double timepitch = driftvelocity*timetick;                         //time sample (cm) 

    int nts = detprop->NumberTimeSamples();
    int nplanes = geom->Nplanes();

    std::vector< std::vector<TH1D*> > signals(nplanes);
    std::vector< std::vector<TH1D*> > pulses(nplanes);

    // get input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    std::vector<art::Ptr<recob::Cluster> > clusterlist;
    if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);

    art::FindManyP<recob::Hit> fm(clusterListHandle, evt, fClusterModuleLabel);

    std::vector< std::vector<int> > Cls(nplanes);
    std::vector< std::vector<CluLen> > clulens(nplanes);
    for (size_t iclu = 0; iclu<clusterlist.size(); ++iclu){

      double w0 = clusterlist[iclu]->StartPos()[0];
      double w1 = clusterlist[iclu]->EndPos()[0];
      double t0 = clusterlist[iclu]->StartPos()[1];
      double t1 = clusterlist[iclu]->EndPos()[1];
      //      t0 -= detprop->GetXTicksOffset(clusterlist[iclu]->View(),0,0);
      //      t1 -= detprop->GetXTicksOffset(clusterlist[iclu]->View(),0,0);
 

      CluLen clulen;
      clulen.index = iclu;
      clulen.length = sqrt(pow((w0-w1)*wire_pitch,2)+pow(detprop->ConvertTicksToX(t0,clusterlist[iclu]->View(),1,0)-detprop->ConvertTicksToX(t1,clusterlist[iclu]->View(),1,0),2));


      std::vector< art::Ptr<recob::Hit> > hitlist = fm.at(iclu);
      std::sort(hitlist.begin(), hitlist.end(), trkf::SortByWire());
      auto theHit = hitlist.begin();
      int hcryo = (*theHit)->WireID().Cryostat;
      int htpc = (*theHit)->WireID().TPC;

      LOG_VERBATIM("CosmicTracker")  <<std::cout << "Cluster " << iclu << " view=" <<clusterlist[iclu]->View() <<  " cryostat=" << hcryo <<" tpc=" << htpc << " length = " << clulen.length <<" start wire=" << w0 << " end wire="<< w1 << " start time=" <<t0 << " end time=" << t1 << std::endl;

      switch(clusterlist[iclu]->View()){
      case geo::kU :
	if (fEnableU) clulens[0].push_back(clulen);
	break;
      case geo::kV :
	if (fEnableV) clulens[1].push_back(clulen);
	break;
      case geo::kZ :
	if (fEnableZ) clulens[2].push_back(clulen);
	break;
      default :
	break;
      }

    }

    //sort clusters based on 2D length
    for (size_t i = 0; i<clulens.size(); ++i){
      std::sort (clulens[i].begin(),clulens[i].end(), myfunction);
      for (size_t j = 0; j<clulens[i].size(); ++j){
	Cls[i].push_back(clulens[i][j].index);
      }
    }

    //calibrate drift times between wire planes using single muons
    std::vector< std::vector<double> > meantime(nplanes);
    for (int i = 0; i<nplanes; ++i){
      for (size_t ic = 0; ic < Cls[i].size(); ++ic){
	TH1D sig(Form("sig_%d_%d",i,int(ic)),Form("sig_%d_%d",i,int(ic)),nts,0,nts);
	TH1D sigint(Form("sigint_%d_%d",i,int(ic)),Form("sigint_%d_%d",i,int(ic)),nts,0,nts);    
	std::vector< art::Ptr<recob::Hit> > hitlist = fm.at(Cls[i][ic]);
	std::sort(hitlist.begin(), hitlist.end(), trkf::SortByWire());
	for(auto theHit = hitlist.begin(); theHit != hitlist.end();  theHit++){
	
	  double time = (*theHit)->PeakTime();
	  time -= detprop->GetXTicksOffset((*theHit)->WireID().Plane,
					   (*theHit)->WireID().TPC,
					   (*theHit)->WireID().Cryostat);

	  double charge = (*theHit)->Charge();
	  int bin = sig.FindBin(time);
	  sig.SetBinContent(bin,sig.GetBinContent(bin)+charge);
	  for (int j = bin; j<=sig.GetNbinsX(); ++j){
	    sigint.SetBinContent(j,sigint.GetBinContent(j)+charge);
	  }
	}
	if (sigint.Integral()) sigint.Scale(1./sigint.GetBinContent(sigint.GetNbinsX()));
	pulses[i].push_back(new TH1D(sig));
	signals[i].push_back(new TH1D(sigint));
	if (hitlist.size()>10){
	  meantime[i].push_back(sig.GetMean());
	}
      }
    }

    bool singletrack = true;
    for (int i = 0; i<nplanes; ++i){
      singletrack = singletrack && meantime[i].size()==1;
    }
    if (singletrack){
      for (int i = 0; i<nplanes; ++i){
	for (int j = i+1; j<nplanes; ++j){
	  dtime[i+j-1]->Fill(meantime[j][0]-meantime[i][0]);
	}
      }
      for (int i = 0; i<nplanes; ++i){
	for (size_t k = 0; k<signals[i].size(); ++k){
	  if (fm.at(Cls[i][k]).size()<10) continue;
	  for (int j = 0; j<signals[i][k]->GetNbinsX(); ++j){
	    double binc = signals[i][k]->GetBinContent(j+1);
	    testsig[i]->SetBinContent(j+1,binc);
	    testpulse[i]->SetBinContent(j+1,pulses[i][k]->GetBinContent(j+1));
	  }
	}
      }
    }


    //matching clusters between different views
    std::vector<int> matched(clusterlist.size());
    for (size_t i = 0; i<clusterlist.size(); ++i) matched[i] = 0;

    std::vector< std::vector<int> > matchedclusters;

    int tpcCheck1=-999;
    int tpcCheck2=-999;


    //  for (int i = 0; i<nplanes-1; ++i){
    //    for (int j = i+1; j<nplanes; ++j){
    for (int i = 0; i<nplanes; ++i){
      for (int j = 0; j<nplanes; ++j){
	for (size_t c1 = 0; c1<Cls[i].size(); ++c1){
	  //
	  //ADD check if the clusters are not in the same TPC
	  //
	  //
	  std::vector< art::Ptr<recob::Hit> > hitlist_TPCcheck1 = fm.at(Cls[i][c1]);
	  std::sort(hitlist_TPCcheck1.begin(), hitlist_TPCcheck1.end(), trkf::SortByWire());
	  for(auto theHit1 = hitlist_TPCcheck1.begin(); theHit1 != hitlist_TPCcheck1.end();  theHit1++){
	    tpcCheck1=(*theHit1)->WireID().TPC;
	    break;
	  }
	  for (size_t c2 = 0; c2<Cls[j].size(); ++c2){
	    std::vector< art::Ptr<recob::Hit> > hitlist_TPCcheck2 = fm.at(Cls[j][c2]);
	    std::sort(hitlist_TPCcheck2.begin(), hitlist_TPCcheck2.end(), trkf::SortByWire());
	    for(auto theHit2 = hitlist_TPCcheck2.begin(); theHit2 != hitlist_TPCcheck2.end();  theHit2++){
	      tpcCheck2=(*theHit2)->WireID().TPC;
	      break;
	    }
	    
	    // check if both are the same view
	    if (clusterlist[Cls[i][c1]]->View()==
		clusterlist[Cls[j][c2]]->View()) continue;
	    //check if both the same TPC
	    if(tpcCheck1!=tpcCheck2) continue;
	    // check if both are already in the matched list
	    if (matched[Cls[i][c1]]==1&&matched[Cls[j][c2]]==1) continue;
	    // KS test between two views in time
	    double ks = signals[i][c1]->KolmogorovTest(signals[j][c2]);
	    hks->Fill(ks);
	    int imatch = -1; //track candidate index
	    int iadd = -1; //cluster index to be inserted
	    if (ks>fKScut){//pass KS test
	      // check both clusters with all matched clusters
	      // if one is already matched, 
	      // check if need to add the other to the same track candidate
	      for (size_t l = 0; l<matchedclusters.size(); ++l){
		for (size_t m = 0; m<matchedclusters[l].size(); ++m){
		  if (matchedclusters[l][m]==Cls[i][c1]){
		    imatch = l; //track candidate
		    iadd = j; //consider the other cluster
		  }
		  else if (matchedclusters[l][m]==Cls[j][c2]){
		    imatch = l; //track candidate
		    iadd = i; //consider the other cluster
		  }
		}
	      }
	      if (imatch>=0){
		if (iadd == i){
		  bool matchview = false;
		  // check if one matched cluster has the same view
		  for (size_t ii = 0; ii<matchedclusters[imatch].size(); ++ii){
		    if (clusterlist[matchedclusters[imatch][ii]]->View()==
			clusterlist[Cls[i][c1]]->View()){
		      matchview = true;
		      //replace if the new cluster has more hits
		      if (fm.at(Cls[i][c1]).size()>fm.at(matchedclusters[imatch][ii]).size()){
			matched[matchedclusters[imatch][ii]] = 0;
			matchedclusters[imatch][ii] = Cls[i][c1];
			matched[Cls[i][c1]] = 1;
		      }
		    }
		  }
		  if (!matchview){//not matched view found, just add
		    matchedclusters[imatch].push_back(Cls[i][c1]);
		    matched[Cls[i][c1]] = 1;
		  }
		}
		else {
		  bool matchview = false;
		  for (size_t jj = 0; jj<matchedclusters[imatch].size(); ++jj){
		    if (clusterlist[matchedclusters[imatch][jj]]->View()==
			clusterlist[Cls[j][c2]]->View()){
		      matchview = true;
		      //replace if it has more hits
		      if (fm.at(Cls[j][c2]).size()>fm.at(matchedclusters[imatch][jj]).size()){
			matched[matchedclusters[imatch][jj]] = 0;
			matchedclusters[imatch][jj] = Cls[j][c2];
			matched[Cls[j][c2]] = 1;
		      }
		    }
		  }
		  if (!matchview){
		    matchedclusters[imatch].push_back(Cls[j][c2]);
		    matched[Cls[j][c2]] = 1;
		  }		
		}
	      }
	      else{
		std::vector<int> tmp;
		tmp.push_back(Cls[i][c1]);
		tmp.push_back(Cls[j][c2]);
		matchedclusters.push_back(tmp);
		matched[Cls[i][c1]]=1;
		matched[Cls[j][c2]]=1;
	      }
	    }//pass KS test
	  }//c2
	}//c1
      }//j
    }//i

    for (size_t i = 0; i<matchedclusters.size(); ++i){
      if (matchedclusters[i].size()) mf::LogVerbatim("CosmicTracker")<<"Track candidate "<<i<<":";
      for (size_t j = 0; j<matchedclusters[i].size(); ++j){
	mf::LogVerbatim("CosmicTracker")<<matchedclusters[i][j];
      }
    } 

    for (int i = 0; i<nplanes; ++i){
      for (size_t j = 0; j<signals[i].size(); ++j){
	delete signals[i][j];
	delete pulses[i][j];
      }
    }

    /////////////////////////////////////////////////////
    /////// 2D Track Matching and 3D Track Reconstruction
    /////////////////////////////////////////////////////
  
    ///Prepare fitter
    double arglist[10];
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    if (fCleanUpHits){
      fitter->SetFCN(myfcn);
      arglist[0] = -1;
      fitter->ExecuteCommand("SET PRIN",arglist,1);
    }
    //fit each cluster in 2D using pol2, iterate once to remove outliers
    for (size_t itrk = 0; itrk<matchedclusters.size(); ++itrk){//loop over tracks

      //all the clusters associated with the current track
      art::PtrVector<recob::Cluster> clustersPerTrack;
      for (size_t iclu = 0; iclu<matchedclusters[itrk].size(); ++iclu){
	art::Ptr <recob::Cluster> cluster(clusterListHandle,matchedclusters[itrk][iclu]);
	clustersPerTrack.push_back(cluster);
      }

      //save time/hit information along track trajectory
      std::vector<std::map<int,double> > vtimemap;
      std::vector<std::map<int,art::Ptr<recob::Hit> > > vhitmap;

      for (size_t iclu = 0; iclu<matchedclusters[itrk].size(); ++iclu){//loop over clusters

	vwire.clear();
	vtime.clear();
	vph.clear();
	//fit hits time vs wire with pol2
	std::vector< art::Ptr<recob::Hit> > hits = fm.at(matchedclusters[itrk][iclu]);
	std::sort(hits.begin(), hits.end(), trkf::SortByWire());
	if (fCleanUpHits){
	  double dtdw = 0;
	  if (clusterlist[matchedclusters[itrk][iclu]]->StartPos()[1]-
	      clusterlist[matchedclusters[itrk][iclu]]->EndPos()[1]){
	    dtdw = (clusterlist[matchedclusters[itrk][iclu]]->EndPos()[0]-
		    clusterlist[matchedclusters[itrk][iclu]]->StartPos()[0])/
	      (clusterlist[matchedclusters[itrk][iclu]]->EndPos()[1]-
	       clusterlist[matchedclusters[itrk][iclu]]->StartPos()[1]);
	  }
	  fitter->SetParameter(0,"p0",clusterlist[matchedclusters[itrk][iclu]]->StartPos()[1]-dtdw-detprop->GetXTicksOffset(clusterlist[matchedclusters[itrk][iclu]]->View(),0,1),0.1,0,0);
	  fitter->SetParameter(1,"p1",clusterlist[matchedclusters[itrk][iclu]]->dTdW(),0.1,0,0);
	  fitter->SetParameter(2,"p2",0,0.1,0,0);
	}
	for (size_t ihit = 0; ihit<hits.size(); ++ihit){//loop over hits
	  geo::WireID hitWireID = hits[ihit]->WireID();
	  unsigned int w = hitWireID.Wire;
	  vwire.push_back(w);
	  double time = hits[ihit]->PeakTime();
	  time -= detprop->GetXTicksOffset(hits[ihit]->WireID().Plane,
					   hits[ihit]->WireID().TPC,
					   hits[ihit]->WireID().Cryostat);
	  vtime.push_back(time);
	  vph.push_back(hits[ihit]->Charge());
	}

	if (fCleanUpHits){
	  arglist[0] = 0;
	  if (vwire.size()>2) fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  else{
	    fitter->SetParameter(0,"p0",vtime[0],0.1,0,0);
	    fitter->SetParameter(1,"p1",0,0.1,0,0);
	    fitter->SetParameter(2,"p2",0,0.1,0,0);
	  }
	  //remove outliers
	  for (auto iw = vwire.begin(), it = vtime.begin(), iph = vph.begin(); iw!=vwire.end(); ){
	    double y = fitter->GetParameter(0)+
	      fitter->GetParameter(1)*(*iw)+
	      fitter->GetParameter(2)*(*iw)*(*iw);
	    if (std::abs(*it-y)>ftoler1){
	      iw = vwire.erase(iw);
	      it = vtime.erase(it);
	      iph = vph.erase(iph);
	    }
	    else{
	      ++iw;
	      ++it;
	      ++iph;
	    }
	  }
	  
	  //refit
	  if (vwire.size()>2) fitter->ExecuteCommand("MIGRAD", arglist, 0);
	}
	
	std::map<int,double> timemap;
	std::map<int,double> phmap;
	std::map<int,art::Ptr<recob::Hit> > hitmap;

	//find hit on each wire along the fitted line
	for (size_t ihit = 0; ihit<hits.size(); ++ihit){//loop over hits
	  geo::WireID hitWireID = hits[ihit]->WireID();
	  unsigned int w = hitWireID.Wire;
	  vwire.push_back(w);
	  double time = hits[ihit]->PeakTime();
	  time -= detprop->GetXTicksOffset(hits[ihit]->WireID().Plane,
					   hits[ihit]->WireID().TPC,
					   hits[ihit]->WireID().Cryostat);
	  double ph = hits[ihit]->Charge();
	  if (fCleanUpHits){
	    if (ph>(phmap[w])){
	      double y = fitter->GetParameter(0)+
		fitter->GetParameter(1)*w+
		fitter->GetParameter(2)*w*w;
	      if (std::abs(time-y)<ftoler2){
		phmap[w] = ph;
		timemap[w] = time;
		hitmap[w] = hits[ihit];
	      }
	    }
	  }
	  else{
	    phmap[w] = ph;
	    timemap[w] = time;
	    hitmap[w] = hits[ihit];
	    //std::cout<<w<<" "<<time<<" "<<ph<<" "<<hits[ihit]->WireID().Plane<<std::endl;
	  }	  
	}//ihit
	vtimemap.push_back(timemap);
	vhitmap.push_back(hitmap);
      }//iclu
      
      
      // Remove isolated hits
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
	auto ihit = vhitmap[iclu].begin();
	for (auto itime = vtimemap[iclu].begin(); itime!=vtimemap[iclu].end();){
	  int diffw0 = 0;
	  int diffw1 = 0;
	  if (itime!=vtimemap[iclu].begin()){
	    auto itime0 = std::prev(itime,1);
	    diffw0 = std::abs(itime->first - itime0->first);
	  }
	  else diffw0 = 99999;
	  auto itime1 = std::next(itime,1);
	  if (itime1!=vtimemap[iclu].end()){
	    diffw1 = abs(itime->first - itime1->first);
	  }
	  else diffw1 = 99999;
	  if (diffw0>fisohitcut&&diffw1>fisohitcut){
	    vtimemap[iclu].erase(itime++);
	    vhitmap[iclu].erase(ihit++);
	  }
	  else{
	    ++itime;
	    ++ihit;
	  }
	}
      }

      // Find two clusters with the most numbers of hits, and time ranges
      int iclu1 = -1;
      int iclu2 = -1;
      int iclu3 = -1;
      unsigned maxnumhits0 = 0;
      unsigned maxnumhits1 = 0;
    
      std::vector<double> tmin(vtimemap.size());
      std::vector<double> tmax(vtimemap.size());
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
	tmin[iclu] = 1e9;
	tmax[iclu] = -1e9;
      }
    
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
	for (auto itime = vtimemap[iclu].begin(); itime!=vtimemap[iclu].end(); ++itime){
	  if (itime->second>tmax[iclu]){
	    tmax[iclu] = itime->second;
	  }
	  if (itime->second<tmin[iclu]){
	    tmin[iclu] = itime->second;
	  }
	}
	if (vtimemap[iclu].size()>maxnumhits0){
	  if (iclu1!=-1){
	    iclu2 = iclu1;
	    maxnumhits1 = maxnumhits0;
	  }
	  iclu1 = iclu;
	  maxnumhits0 = vtimemap[iclu].size();
	}
	else if (vtimemap[iclu].size()>maxnumhits1){
	  iclu2 = iclu;
	  maxnumhits1 = vtimemap[iclu].size();
	}
      }
    
      std::swap(iclu1,iclu2); //now iclu1 has fewer hits than iclu2

      for (int iclu = 0; iclu<(int)vtimemap.size(); ++iclu){
	if (iclu!=iclu1&&iclu!=iclu2) iclu3 = iclu;
      }
    
      if (iclu1!=-1&&iclu2!=-1){
	//select hits in a common time range
	auto ihit = vhitmap[iclu1].begin();
	auto itime = vtimemap[iclu1].begin();
	while (itime!=vtimemap[iclu1].end()){
	  if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
	      itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
	    vtimemap[iclu1].erase(itime++);
	    vhitmap[iclu1].erase(ihit++);
	  }
	  else{
	    ++itime;
	    ++ihit;
	  }
	}

	ihit = vhitmap[iclu2].begin();
	itime = vtimemap[iclu2].begin();
	while (itime!=vtimemap[iclu2].end()){
	  if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
	      itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
	    vtimemap[iclu2].erase(itime++);
	    vhitmap[iclu2].erase(ihit++);
	  }
	  else{
	    ++itime;
	    ++ihit;
	  }
	}
      
	//if one cluster is empty, replace it with iclu3
	if (!vtimemap[iclu1].size()){
	  if (iclu3!=-1){
	    std::swap(iclu3,iclu1);
	  }
	}
	if (!vtimemap[iclu2].size()){
	  if (iclu3!=-1){
	    std::swap(iclu3,iclu2);
	    std::swap(iclu1,iclu2);
	  }
	}
	if ((!vtimemap[iclu1].size())||(!vtimemap[iclu2].size())) continue;

	size_t spStart = spcol->size();
	std::vector<recob::SpacePoint> spacepoints;
	//TVector3 startpointVec,endpointVec, DirCos;

	bool rev = false;
	auto times1 = vtimemap[iclu1].begin();
	auto timee1 = vtimemap[iclu1].end();
	--timee1;
	auto times2 = vtimemap[iclu2].begin();
	auto timee2 = vtimemap[iclu2].end();
	--timee2;

	double ts1 = times1->second;
	double te1 = timee1->second;
	double ts2 = times2->second;
	double te2 = timee2->second;
            
	//find out if we need to flip ends
	if (std::abs(ts1-ts2)+std::abs(te1-te2)>std::abs(ts1-te2)+std::abs(te1-ts2)){
	  rev = true;
	}
//	std::cout<<times1->first<<" "<<times1->second<<std::endl;
//	std::cout<<timee1->first<<" "<<timee1->second<<std::endl;
//	std::cout<<times2->first<<" "<<times2->second<<std::endl;
//	std::cout<<timee2->first<<" "<<timee2->second<<std::endl;
	std::vector<double> vtracklength;
      
	for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
	
	  double tracklength = 0;
	  if (vtimemap[iclu].size()==1){
	    tracklength = wire_pitch;
	  }
	  else{
	    double t0 = 0., w0 = 0.;
	    for (auto iw = vtimemap[iclu].begin(); iw!=vtimemap[iclu].end(); ++iw){
	      if (iw==vtimemap[iclu].begin()){
		w0 = iw->first;
		t0 = iw->second;
	      }
	      else{
		tracklength += std::sqrt(std::pow((iw->first-w0)*wire_pitch,2)+std::pow((iw->second-t0)*timepitch,2));
		w0 = iw->first;
		t0 = iw->second;	     
	      }
	    }
	  }
	  vtracklength.push_back(tracklength);
	}
      
	std::map<int,int> maxhitsMatch;

	auto ihit1 = vhitmap[iclu1].begin();
	for (auto itime1 = vtimemap[iclu1].begin(); 
	     itime1 != vtimemap[iclu1].end(); 
	     ++itime1, ++ihit1){//loop over min-hits
	  art::PtrVector<recob::Hit> sp_hits;
	  sp_hits.push_back(ihit1->second);
	  double hitcoord[3];
	  double length1 = 0;
	  if (vtimemap[iclu1].size()==1){
	    length1 = wire_pitch;
	  }
	  else{
	    for (auto iw1 = vtimemap[iclu1].begin(); iw1!=itime1; ++iw1){
	      auto iw2 = iw1;
	      ++iw2;
	      length1 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)
				   +std::pow((iw1->second-iw2->second)*timepitch,2));
	    }
	  }
	  double difference = 1e10; //distance between two matched hits
	  auto matchedtime = vtimemap[iclu2].end();
	  auto matchedhit  = vhitmap[iclu2].end();
	
	  auto ihit2 = vhitmap[iclu2].begin();
	  for (auto itime2 = vtimemap[iclu2].begin(); 
	       itime2!=vtimemap[iclu2].end(); 
	       ++itime2, ++ihit2){//loop over max-hits
	    if (maxhitsMatch[itime2->first]) continue;
	    double length2 = 0;
	    if (vtimemap[iclu2].size()==1){
	      length2 = wire_pitch;
	    }
	    else{
	      for (auto iw1 = vtimemap[iclu2].begin(); iw1!=itime2; ++iw1){
		auto iw2 = iw1;
		++iw2;
		length2 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)+std::pow((iw1->second-iw2->second)*timepitch,2));
	      }
	    }
	    if (rev) length2 = vtracklength[iclu2] - length2;
	    length2 = vtracklength[iclu1]/vtracklength[iclu2]*length2;
	    bool timematch = std::abs(itime1->second-itime2->second)<ftmatch;
	    if (timematch &&std::abs(length2-length1)<difference){
	      difference = std::abs(length2-length1);
	      matchedtime = itime2;
	      matchedhit = ihit2;
	    }
	  }//loop over hits2
	  if (difference<fsmatch){
	    hitcoord[0] = matchedtime->second*detprop->GetXTicksCoefficient();
	    hitcoord[1] = -1e10;
	    hitcoord[2] = -1e10;
	    /*	    	    geom->ChannelsIntersect((ihit1->second)->Wire()->RawDigit()->Channel(),
				    (matchedhit->second)->Wire()->RawDigit()->Channel(),
				    hitcoord[1],hitcoord[2]);
	    */


	    //WireID is the exact segment of the wire where the hit is on (1 out of 3 for the 35t)
	    geo::WireID c1=(ihit1->second)->WireID();
	    geo::WireID c2=(matchedhit->second)->WireID();
	    //	    std::vector< geo::WireID > chan1wires, chan2wires; 
	    //	    chan1wires = geom->ChannelToWire(c1);
	    //	    chan2wires = geom->ChannelToWire(c2);
	    geo::WireIDIntersection tmpWIDI;
	    

	    //
	    //   Outputs for debugging purposes.
	    //
	    //
	    //
	    //
	    /*
	    double w1_Start[3] = {0.};
	    double w1_End[3]   = {0.};
	    double w2_Start[3] = {0.};
	    double w2_End[3]   = {0.};
	    // get the endpoints to see if i1 and i2 even intersect
	    geom->WireEndPoints(c1.Cryostat, c1.TPC, c1.Plane, c1.Wire, w1_Start, w1_End);
	    geom->WireEndPoints(c2.Cryostat, c2.TPC, c2.Plane, c2.Wire, w2_Start, w2_End);

	    mf::LogVerbatim("Summary") <<"TPC :c1 " << c1.TPC << "   c2 " << c2.TPC;
	    mf::LogVerbatim("Summary") <<"Cryo :c1 " << c1.Cryostat << "   c2 " << c2.Cryostat;
	    mf::LogVerbatim("Summary") <<"Plane:c1 " << c1.Plane << "   c2 " << c2.Plane;
	    mf::LogVerbatim("Summary") <<"Wire:c1 " << c1.Wire << "   c2 " << c2.Wire;

	    bool overlapY         = geom->ValueInRange( w1_Start[1], w2_Start[1], w2_End[1] ) ||
	                                  geom->ValueInRange( w1_End[1],   w2_Start[1], w2_End[1] );
	    bool overlapY_reverse = geom->ValueInRange( w2_Start[1], w1_Start[1], w1_End[1] ) ||
	                                   geom->ValueInRange( w2_End[1],   w1_Start[1], w1_End[1] );
	         
	    bool overlapZ         = geom->ValueInRange( w1_Start[2], w2_Start[2], w2_End[2] ) ||
	                                   geom->ValueInRange( w1_End[2],   w2_Start[2], w2_End[2] );
	    bool overlapZ_reverse = geom->ValueInRange( w2_Start[2], w1_Start[2], w1_End[2] ) ||
	                                   geom->ValueInRange( w2_End[2],   w1_Start[2], w1_End[2] );
	    if(std::abs(w2_Start[2] - w2_End[2]) < 0.01) overlapZ = overlapZ_reverse;
	    mf::LogVerbatim("Summary") << "overlapY:" << overlapY << "   " <<
	      "overlapY_reverse:" << overlapY_reverse <<"    "<<
	      "overlapZ:" << overlapZ <<"    "<<
	      "overlapZ_reverse:"<< overlapZ_reverse<<"     ";
	    */
	    //
	    //
	    //  End of outputs for debugging purposes
	    //
	    //
	   
		    bool sameTpcOrNot=geom->WireIDsIntersect(c1,c2, tmpWIDI);
		    
	    //   	    bool sameTpcOrNot=APAGeometryAlg::APAChannelsIntersect((ihit1->second)->Wire()->RawDigit()->Channel(),
	    //								   (matchedhit->second)->Wire()->RawDigit()->Channel(),
	    //		   						   tmpWIDI
	    //								   )
		    if(sameTpcOrNot)
		      {
			hitcoord[1]=tmpWIDI.y;
			hitcoord[2]=tmpWIDI.z;
		      }

	    if (hitcoord[1]>-1e9&&hitcoord[2]>-1e9){
	      maxhitsMatch[matchedtime->first] = 1;
	      sp_hits.push_back(matchedhit->second);
	    }
	  }
	  if (sp_hits.size()>1){
	    double err[6] = {util::kBogusD};
	    recob::SpacePoint mysp(hitcoord, 
				   err, 
				   util::kBogusD, 
				   spStart + spacepoints.size());//3d point at end of track
	    spacepoints.push_back(mysp);
	    spcol->push_back(mysp);	
	    util::CreateAssn(*this, evt, *spcol, sp_hits, *shassn);
	  }
	}//loop over hits1
      
	size_t spEnd = spcol->size();

	if (fSortDir == "+x"){
	  std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_x0);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_x0);
	}
	if (fSortDir == "-x"){
	  std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_x1);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_x1);
	}
	if (fSortDir == "+y"){
	  std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_y0);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_y0);
	}
	if (fSortDir == "-y"){
	  std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_y1);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_y1);
	}
	if (fSortDir == "+z"){
	  std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_z0);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_z0);
	}
	if (fSortDir == "-z"){
	  std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_z1);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_z1);
	}

	// Add the 3D track to the vector of the reconstructed tracks
	if(spacepoints.size()>0){

	  // make a vector of the trajectory points along the track
	  std::vector<TVector3> xyz(spacepoints.size());
	  for(size_t s = 0; s < spacepoints.size(); ++s){
	    xyz[s] = TVector3(spacepoints[s].XYZ());
	  }	
	  //Calculate track direction cosines 
	  TVector3 startpointVec,endpointVec, DirCos;
	  startpointVec = xyz[0];
	  endpointVec = xyz.back();
	  DirCos = endpointVec - startpointVec;
	  //SetMag casues a crash if the magnitude of the vector is zero
	  try
	    {
	      DirCos.SetMag(1.0);//normalize vector
	    }
	  catch(...){std::cout<<"The Spacepoint is infinitely small"<<std::endl;
	    continue;
	  }
	  //std::cout<<DirCos.x()<<" "<<DirCos.y()<<" "<<DirCos.z()<<std::endl;
	  std::vector<TVector3> dircos(spacepoints.size(), DirCos);	  
	  if (fDirSPS){//calculat direction for each spacepoint
	    for (int s = 0; s < int(xyz.size()); ++s){
	      int np = 0;
	      std::vector<double> vx;
	      std::vector<double> vy;
	      std::vector<double> vz;
	      std::vector<double> vs;
	      vx.push_back(xyz[s].x());
	      vy.push_back(xyz[s].y());
	      vz.push_back(xyz[s].z());
	      vs.push_back(0);
	      ++np;
	      for (int ip = 1; ip<int(xyz.size()); ++ip){
		if (s-ip>=0){
		  vx.push_back(xyz[s-ip].x());
		  vy.push_back(xyz[s-ip].y());
		  vz.push_back(xyz[s-ip].z());
		  double dis = 0;
		  for (int j = s-ip; j<s; ++j){
		    dis += -sqrt(pow(xyz[j].x()-xyz[j+1].x(),2)+
				 pow(xyz[j].y()-xyz[j+1].y(),2)+
				 pow(xyz[j].z()-xyz[j+1].z(),2));
		  }
		  vs.push_back(dis);
		  ++np;
		  if (np==5) break;
		}
		if (s+ip<int(xyz.size())){
		  vx.push_back(xyz[s+ip].x());
		  vy.push_back(xyz[s+ip].y());
		  vz.push_back(xyz[s+ip].z());
		  double dis = 0;
		  for (int j = s; j<s+ip; ++j){
		    dis += sqrt(pow(xyz[j].x()-xyz[j+1].x(),2)+
				pow(xyz[j].y()-xyz[j+1].y(),2)+
				pow(xyz[j].z()-xyz[j+1].z(),2));
		  }
		  vs.push_back(dis);
		  ++np;
		  if (np==5) break;
		}
	      }
	      double kx = 0, ky = 0, kz = 0;
	      if (np>=2){//at least two points
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
		  kx = pol->GetParameter(1);
		  //std::cout<<xyz3d[0]<<" "<<kx<<std::endl;
		}
		catch(...){
		  mf::LogWarning("CosmicTracker") <<"Fitter failed";
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
		  ky = pol->GetParameter(1);
		  //std::cout<<xyz3d[1]<<" "<<ky<<std::endl;
		}
		catch(...){
		  mf::LogWarning("CosmicTracker") <<"Fitter failed";
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
		  kz = pol->GetParameter(1);
		  //std::cout<<xyz3d[2]<<" "<<kz<<std::endl;
		}
		catch(...){
		  mf::LogWarning("CosmicTracker") <<"Fitter failed";
		}
		delete zs;
		if (kx||ky||kz){
		  double tot = sqrt(kx*kx+ky*ky+kz*kz);
		  kx /= tot;
		  ky /= tot;
		  kz /= tot;
		  dircos[s].SetXYZ(kx,ky,kz);
		  //std::cout<<s<<" "<<kx<<" "<<ky<<" "<<kz<<std::endl;
		}
	      }//np>=2
	    }//loop over space points
	  }
	  std::vector< std::vector<double> > dQdx;
	  std::vector<double> mom(2, util::kBogusD);
	  tcol->push_back(recob::Track(xyz, dircos, dQdx, mom, tcol->size()));
	
	  // make associations between the track and space points
	  util::CreateAssn(*this, evt, *tcol, *spcol, *tspassn, spStart, spEnd);
	
	  // now the track and clusters
	  util::CreateAssn(*this, evt, *tcol, clustersPerTrack, *tcassn);
	
	  // and the hits and track
	  if (!fdebug){
	    art::FindManyP<recob::Hit> fmh(clustersPerTrack, evt, fClusterModuleLabel);
	    for(size_t cpt = 0; cpt < clustersPerTrack.size(); ++cpt)
	      util::CreateAssn(*this, evt, *tcol, fmh.at(cpt), *thassn);
	  }
	  else{
	    std::vector<art::Ptr<recob::Hit> > trkhits;
	    for(auto ihit = vhitmap[iclu1].begin(); ihit!=vhitmap[iclu1].end();++ihit){
	      trkhits.push_back(ihit->second);
	    }
	    for(auto ihit = vhitmap[iclu2].begin(); ihit!=vhitmap[iclu2].end();++ihit){
	      trkhits.push_back(ihit->second);
	    }
	    util::CreateAssn(*this, evt, *tcol, trkhits, *thassn);
	  }
	}
      }//if iclu1&&iclu2
    }//itrk		 

    mf::LogVerbatim("Summary") << std::setfill('-') 
			       << std::setw(175) 
			       << "-" 
			       << std::setfill(' ');
    mf::LogVerbatim("Summary") << "CosmicTracker Summary:";
    for(unsigned int i = 0; i<tcol->size(); ++i) mf::LogVerbatim("Summary") << tcol->at(i) ;
    mf::LogVerbatim("Summary") << "CosmicTracker Summary End:";
    
    evt.put(std::move(tcol));
    evt.put(std::move(spcol));
    evt.put(std::move(tspassn));
    evt.put(std::move(tcassn));
    evt.put(std::move(thassn));
    evt.put(std::move(shassn));
    
    return;
  }
  
  DEFINE_ART_MODULE(CosmicTracker)

} // namespace

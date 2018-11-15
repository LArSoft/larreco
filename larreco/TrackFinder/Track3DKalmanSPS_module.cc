////////////////////////////////////////////////////////////////////////
//
// \file Track3DKalmanSPS.cxx
//
// \author echurch@fnal.gov
//
//  This algorithm is designed to reconstruct 3D tracks through  
//  GENFIT's Kalman filter.
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iterator> // std::distance()
#include <algorithm> // std::sort()

// ROOT includes
#include "TVectorD.h" // TVector3
#include "TFile.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TDatabasePDG.h"
#include "TTree.h"
#include "TMatrixT.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/exception.h" 

#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"

// GENFIT includes
#include "larreco/Genfit/GFAbsTrackRep.h"
#include "larreco/Genfit/RKTrackRep.h"
#include "larreco/Genfit/GFConstField.h"
#include "larreco/Genfit/GFFieldManager.h"
#include "larreco/Genfit/PointHit.h"
#include "larreco/Genfit/GFTrack.h"
#include "larreco/Genfit/GFKalman.h"
 
// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
//\todo Reconstruction Producers should never include SimulationBase objects
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/SpacePointAlg.h"





static bool sp_sort_3dz(const art::Ptr<recob::SpacePoint>& h1, const art::Ptr<recob::SpacePoint>& h2)
{
  const double* xyz1 = h1->XYZ();
  const double* xyz2 = h2->XYZ();
  return xyz1[2] < xyz2[2];
}
static bool sp_sort_3dy(const art::Ptr<recob::SpacePoint>& h1, const art::Ptr<recob::SpacePoint>& h2)
{
  const double* xyz1 = h1->XYZ();
  const double* xyz2 = h2->XYZ();
  return xyz1[1] < xyz2[1];
}
static bool sp_sort_3dx(const art::Ptr<recob::SpacePoint>& h1, const art::Ptr<recob::SpacePoint>& h2)
{
  const double* xyz1 = h1->XYZ();
  const double* xyz2 = h2->XYZ();
  return xyz1[0] < xyz2[0];
}

static bool sp_sort_nsppts(const art::PtrVector<recob::SpacePoint>& h1, const art::PtrVector<recob::SpacePoint>& h2)
{
  const unsigned int s1 = h1.size();
  const unsigned int s2 = h2.size();
  return s1 > s2;
}


namespace trkf {

  class Track3DKalmanSPS : public art::EDProducer {
    
  public:
    
    explicit Track3DKalmanSPS(fhicl::ParameterSet const& pset);
    virtual ~Track3DKalmanSPS();
    
    //////////////////////////////////////////////////////////
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    double energyLossBetheBloch(const double& mass,
				const double p
				);
  private:

    void rotationCov(TMatrixT<Double_t>  &cov, const TVector3 &u, const TVector3 &v);
    std::vector <double> dQdxCalc(const art::FindManyP<recob::Hit> &h, const art::PtrVector<recob::SpacePoint> &s, const TVector3 &p, const TVector3 &d );

    std::string     fClusterModuleLabel;// label for input collection
    std::string     fSpptModuleLabel;// label for input collection
    std::string     fGenieGenModuleLabel;// label for input MC single particle generator
    std::string     fG4ModuleLabel;// label for input MC single particle generator
    std::string     fSortDim; // direction in which to sort spacepoints

  //  TFile *fileGENFIT;
    TTree *tree;

    TMatrixT<Double_t> *stMCT;
    TMatrixT<Double_t> *covMCT;
    TMatrixT<Double_t> *stREC;
    TMatrixT<Double_t> *covREC;
    Float_t chi2;
    Float_t chi2ndf;
    int     fcont;

    Float_t *fpRECL;
    Float_t *fpREC;
    Float_t *fpMCMom;
    Float_t *fpMCPos;
    Float_t *fState0;
    Float_t *fCov0;
    int nfail;
    int ndf;
    int nchi2rePass;
    int ispptvec;
    int nspptvec;
    unsigned int evtt;
    unsigned int nTrks;
    unsigned int fptsNo;
    Float_t *fshx;
    Float_t *fshy;
    Float_t *fshz;
    Float_t *feshx;
    Float_t *feshy;
    Float_t *feshz;
    Float_t *feshyz;
    Float_t *fupdate;
    Float_t *fchi2hit;
    Float_t *fth;
    Float_t *feth;
    Float_t *fedudw;
    Float_t *fedvdw;
    Float_t *feu;
    Float_t *fev;
    Float_t *fsep;
    Float_t *fdQdx;
    unsigned int fDimSize; // if necessary will get this from pset in constructor.
    Float_t *fPCmeans;
    Float_t *fPCevals;
    Float_t *fPCsigmas;
    Float_t *fPC1;
    Float_t *fPC2;
    Float_t *fPC3;

    std::vector<double> fPosErr;
    std::vector<double> fMomErr;
    std::vector<double> fMomStart;
    double fPerpLim;
    bool fDoFit;
    int fNumIt;
    uint16_t fMinNumSppts;
    double fErrScaleS;
    double fErrScaleM;
    int fDecimate;
    double fMaxUpdate;
    int fDecimateU;
    double fDistanceU;
    double fMaxUpdateU;
    double fMomLow;
    double fMomHigh;
    int fPdg;
    double fChi2Thresh;
    int fMaxPass;

    genf::GFAbsTrackRep *repMC;
    genf::GFAbsTrackRep *rep;

  protected: 
    
  
  }; // class Track3DKalmanSPS

}

namespace trkf {

//-------------------------------------------------
  Track3DKalmanSPS::Track3DKalmanSPS(fhicl::ParameterSet const& pset) 
    : fDoFit(true)
    , fNumIt(5)
    , fMinNumSppts(5)
    , fErrScaleS(1.0)
    , fErrScaleM(1.0)
    , fDecimate(1)
    , fMaxUpdate(0.10)
    , fDecimateU(1)
    , fDistanceU(1.)
    , fMaxUpdateU(0.10)
    , fMomLow(0.001)
    , fMomHigh(100.)
    , fPdg(-13)
    , fChi2Thresh(12.0E12)
    , fMaxPass (1)
  {
    
    this->reconfigure(pset);
    
    produces< std::vector<recob::Track>                  >();
    produces<art::Assns<recob::Track, recob::Cluster>    >();
    produces<art::Assns<recob::Track, recob::SpacePoint> >();
    produces<art::Assns<recob::Track, recob::Hit>        >();

  }

//-------------------------------------------------
  void Track3DKalmanSPS::reconfigure(fhicl::ParameterSet const& pset) 
  {
    
    fClusterModuleLabel    = pset.get< std::string >("ClusterModuleLabel");
    fSpptModuleLabel       = pset.get< std::string >("SpptModuleLabel");
    fGenieGenModuleLabel   = pset.get< std::string >("GenieGenModuleLabel");
    fG4ModuleLabel         = pset.get< std::string >("G4ModuleLabel");
    fPosErr                = pset.get< std::vector < double >  >("PosErr3");   // resolution. cm
    fMomErr                = pset.get< std::vector < double >  >("MomErr3");   // GeV
    fMomStart              = pset.get< std::vector < double >  >("MomStart3"); // 
    fPerpLim               = pset.get< double  >("PerpLimit", 1.e6); // PCA cut.
    fDoFit                 = pset.get< bool  >("DoFit", true); // Der.
    fNumIt                 = pset.get< int  >("NumIt", 5); // Number x2 passes per fit.
    fMinNumSppts           = pset.get< int  >("MinNumSppts", 5); // Min number of sppts in vector to bother fitting
    fErrScaleS             = pset.get< double >("ErrScaleSim", 1.0); // error scale.
    fErrScaleM             = pset.get< double >("ErrScaleMeas", 1.0); // error scale.
    fDecimate              = pset.get< int  >("DecimateC", 40); // Sparsify data.
    fMaxUpdate             = pset.get< double >("MaxUpdateC", 0.1); // 0-out. 
    fDecimateU             = pset.get< int  >("DecimateU", 100);// Sparsify data.
    fDistanceU             = pset.get< double >("DistanceU", 10.0);// Require this separation on uncontained 2nd pass.
    fMaxUpdateU            = pset.get< double >("MaxUpdateU", 0.02); // 0-out. 
    fMomLow                = pset.get< double >("MomLow", 0.01); // Fit Range. 
    fMomHigh               = pset.get< double >("MomHigh", 20.); // Fit Range. 
    fPdg                   = pset.get< int  >("PdgCode", -13); // mu+ Hypothesis.
    fChi2Thresh            = pset.get< double >("Chi2HitThresh", 12.0E12); //For Re-pass.
    fSortDim               = pset.get< std::string> ("SortDirection", "z"); // case sensitive
    fMaxPass               = pset.get< int  >("MaxPass", 2); // mu+ Hypothesis.
    bool fGenfPRINT;
    if (pset.get_if_present("GenfPRINT", fGenfPRINT)) {
      LOG_WARNING("Track3DKalmanSPS_GenFit")
        << "Parameter 'GenfPRINT' has been deprecated.\n"
        "Please use the standard message facility to enable GenFit debug output.";
      // A way to enable debug output is all of the following:
      // - compile in debug mode (no optimization, no profiling)
      // - if that makes everything too noisy, add to have everything else quiet
      //   services.message.debugModules: [ "Track3DKalmanSPS" ]
      // - to print all the GenFit debug messages, set
      //   services.message.destinations.LogDebugFile.categories.Track3DKalmanSPS_GenFit.limit: -1
      //   (assuming there is a LogDebugFile destination already; for example
      //   see the settings in uboonecode/uboone/Utilities/services_microboone.fcl )
    }
  }

//-------------------------------------------------
  Track3DKalmanSPS::~Track3DKalmanSPS()
  {
  }

//-------------------------------------------------
// stolen, mostly, from GFMaterialEffects.
  double Track3DKalmanSPS::energyLossBetheBloch(const double& mass,
						const double p=1.5
						)
  {
    const double charge(1.0);
    const double mEE(188.); // eV 
    const double matZ(18.);
    const double matA(40.);
    const double matDensity(1.4);
    const double me(0.000511);
    
    double beta = p/std::sqrt(mass*mass+p*p);
    double gammaSquare = 1./(1.0 - beta*beta);
    // 4pi.r_e^2.N.me = 0.307075, I think.
    double dedx = 0.307075*matDensity*matZ/matA/(beta*beta)*charge*charge;
    double massRatio = me/mass;
    // me=0.000511 here is in GeV. So mEE comes in here in eV.
    double argument = gammaSquare*beta*beta*me*1.E3*2./((1.E-6*mEE) * std::sqrt(1+2*std::sqrt(gammaSquare)*massRatio + massRatio*massRatio));
    
    if (mass==0.0) return(0.0);
    if (argument <= exp(beta*beta))
      { 
	dedx = 0.;
      }
    else{
      dedx *= (log(argument)-beta*beta); // Bethe-Bloch [MeV/cm]
      dedx *= 1.E-3;  // in GeV/cm, hence 1.e-3
      if (dedx<0.) dedx = 0.;
    }
    return dedx;
  }

  void Track3DKalmanSPS::rotationCov(TMatrixT<Double_t> &cov, const TVector3 &u, const TVector3 &v)
  {
    TVector3 xhat(1.0,0.0,0.0);
    TVector3 yhat(0.0,1.0,0.0);
    TVector3 zhat(0.0,0.0,1.0);
    TVector3 w(u.Cross(v));
    TVector3 uprime(u);
    TVector3 vprime(w.Cross(xhat)); // vprime now lies in yz plane
    Double_t angle(v.Angle(vprime));/* This is the angle through which v 
				       must rotate. */
    uprime.Rotate(angle,w);// u now is rotated the same amount
    if (uprime*xhat<0)
      {
      	uprime.Rotate(TMath::Pi(),w);
	vprime.Rotate(TMath::Pi(),w);
	angle+=TMath::Pi();
      }
    // Build the block-diagonal 5x5 matrix 
    double c = TMath::Cos(angle), s = TMath::Sin(angle);
    TMatrixT<Double_t> rot(5,5);
    rot[0][0] = 1.0;
    rot[1][1] =  c;
    rot[1][2] =  s;
    rot[2][1] = -s;
    rot[2][2] =  c;
    rot[3][3] =  c;
    rot[3][4] =  s;
    rot[4][3] = -s;
    rot[4][4] =  c;
    
    cov=rot*cov;
  }  

   std::vector<double> Track3DKalmanSPS::dQdxCalc(const art::FindManyP<recob::Hit> &h, const art::PtrVector<recob::SpacePoint> &s, const TVector3 &dir, const TVector3 &loc )
     {
      // For now just Collection plane.
      // We should loop over all views, more generally.
      geo::SigType_t sig(geo::kCollection);
      art::ServiceHandle<geo::Geometry> geom;
      static art::PtrVector<recob::SpacePoint>::const_iterator sstart(s.begin());
      //      art::PtrVector<recob::SpacePoint>::const_iterator sppt = sstart;
      art::PtrVector<recob::SpacePoint>::const_iterator sppt = s.begin();
      std::vector <double> v;



      double mindist(100.0); // cm
      auto spptminIt(sppt);
      while (sppt != s.end())
	{
	  if (((**sppt).XYZ() - loc).Mag() < mindist)
	    {
	      double dist = ((**sppt).XYZ() - loc).Mag(); 
	      
	      // Jump out if we're as close as 0.1 mm away.
	      if (dist<mindist) 
		{
		  mindist = dist;  
		  spptminIt = sppt;
		  if (mindist < 0.01) break; 
		}
	    }
	  sppt++;
	}
      sstart = spptminIt; // for next time.
      unsigned int ind(std::distance(s.begin(),spptminIt));



      std::vector< art::Ptr<recob::Hit> > hitlist = h.at(ind);

      double wirePitch = 0.;
      double angleToVert = 0;
      //      unsigned int tpc1;
      unsigned int plane1;
      double charge = 0.;

      for(std::vector< art::Ptr<recob::Hit> >::const_iterator ihit = hitlist.begin();
	  ihit != hitlist.end(); ++ihit) 
	{
	  const recob::Hit& hit1 = **ihit;
	  //	  if (hit1.View() != view) continue;
	  if (hit1.SignalType() != sig) continue;
	  geo::WireID hit1WireID = hit1.WireID();
	  //	  tpc1 = hit1WireID.TPC;
	  plane1 = hit1WireID.Plane;
	  charge = hit1.Integral();
	  wirePitch = geom->WirePitch(plane1);
	  angleToVert = geom->Plane(plane1).Wire(0).ThetaZ(false) - 0.5*TMath::Pi();
	}
      
      double cosgamma = TMath::Abs(TMath::Sin(angleToVert)*dir.Y() +
				   TMath::Cos(angleToVert)*dir.Z());      
      // if(cosgamma < 1.e-5)
	//	throw cet::exception("Track") << "cosgamma is basically 0, that can't be right\n";

      v.push_back(charge/wirePitch/cosgamma);
      //      std::cout << " Track3DKalmanSPS::dQdxCalc() : For loc.XYZ() hit is ... " << ind << " and v is " << v.back() << std::endl;

      return v;

    }

//-------------------------------------------------
  void Track3DKalmanSPS::beginJob()
  {


    art::ServiceHandle<art::TFileService> tfs;
    
    stMCT  = new TMatrixT<Double_t>(5,1);
    covMCT = new TMatrixT<Double_t>(5,5);
    stREC  = new TMatrixT<Double_t>(5,1);
    covREC = new TMatrixT<Double_t>(5,5);
  
    fpMCMom = new Float_t[4];
    fpMCPos = new Float_t[4];
    fpREC = new Float_t[4];
    fpRECL = new Float_t[4];
    fState0 = new Float_t[5];
    fCov0 = new Float_t[25];
    fDimSize = 20000; // if necessary will get this from pset in constructor.
    
    fshx = new Float_t[fDimSize];
    fshy = new Float_t[fDimSize];
    fshz = new Float_t[fDimSize];
    feshx = new Float_t[fDimSize];
    feshy = new Float_t[fDimSize];
    feshz = new Float_t[fDimSize];
    feshyz = new Float_t[fDimSize];
    fupdate = new Float_t[fDimSize];
    fchi2hit = new Float_t[fDimSize];
    fth  = new Float_t[fDimSize];
    feth = new Float_t[fDimSize];
    fedudw = new Float_t[fDimSize];
    fedvdw = new Float_t[fDimSize];
    feu = new Float_t[fDimSize];
    fev = new Float_t[fDimSize];
    fsep = new Float_t[fDimSize];
    fdQdx = new Float_t[fDimSize];
  
    fPC1 = new Float_t[3];
    fPC2 = new Float_t[3];
    fPC3 = new Float_t[3];
    fPCmeans = new Float_t[3];
    fPCsigmas = new Float_t[3];
    fPCevals = new Float_t[3];
    
//   //TFile fileGENFIT("GENFITout.root","RECREATE");

  
    tree = tfs->make<TTree>("GENFITttree","GENFITttree");
    //tree->Branch("stMCT",&stMCT,"stMCT[5]/F"); // "TMatrixT<Double_t>"

  //tree->Branch("stMCT","TMatrixD",&stMCT,64000,0);
  //tree->Branch("covMCT",covMCT,"covMCT[25]/F");
    tree->Branch("covMCT","TMatrixD",&covMCT,64000,0);
    tree->Branch("stREC",fState0,"stREC[5]/F");
  //tree->Branch("stREC","TMatrixD",&stREC,64000,0);
    tree->Branch("covREC",fCov0,"covREC[25]/F");
  //tree->Branch("covREC","TMatrixD",&covREC,64000,0);
  
    tree->Branch("nchi2rePass",&nchi2rePass,"nchi2rePass/I");
    tree->Branch("ispptvec",&ispptvec,"ispptvec/I");
    tree->Branch("chi2",&chi2,"chi2/F");
    tree->Branch("nfail",&nfail,"nfail/I");
    tree->Branch("ndf",&ndf,"ndf/I");
    tree->Branch("evtNo",&evtt,"evtNo/I");
    tree->Branch("nspptvec",&nspptvec,"nspptvec/I");
    tree->Branch("chi2ndf",&chi2ndf,"chi2ndf/F");
  
    tree->Branch("trkNo",&nTrks,"trkNo/I");
    tree->Branch("ptsNo",&fptsNo,"ptsNo/I");
    tree->Branch("cont",&fcont,"cont/I"); //O? Yes, O. Not 0, not L, ...
    tree->Branch("shx",fshx,"shx[ptsNo]/F");
    tree->Branch("shy",fshy,"shy[ptsNo]/F");
    tree->Branch("shz",fshz,"shz[ptsNo]/F");
    tree->Branch("sep",fsep,"sep[ptsNo]/F");
    tree->Branch("dQdx",fdQdx,"dQdx[ptsNo]/F");
    tree->Branch("eshx",feshx,"eshx[ptsNo]/F");
    tree->Branch("eshy",feshy,"eshy[ptsNo]/F");
    tree->Branch("eshz",feshz,"eshz[ptsNo]/F");
    tree->Branch("eshyz",feshyz,"eshyz[ptsNo]/F");  
    tree->Branch("update",fupdate,"update[ptsNo]/F");
    tree->Branch("chi2hit",fchi2hit,"chi2hit[ptsNo]/F");
    tree->Branch("th",fth,"th[ptsNo]/F");  
    tree->Branch("eth",feth,"eth[ptsNo]/F");
    tree->Branch("edudw",fedudw,"edudw[ptsNo]/F");
    tree->Branch("edvdw",fedvdw,"edvdw[ptsNo]/F");
    tree->Branch("eu",feu,"eu[ptsNo]/F");
    tree->Branch("ev",fev,"ev[ptsNo]/F");


    tree->Branch("pcMeans", fPCmeans,"pcMeans[3]/F");
    tree->Branch("pcSigmas",fPCsigmas,"pcSigmas[3]/F");
    tree->Branch("pcEvals", fPCevals,"pcEvals[3]/F");
    tree->Branch("pcEvec1",fPC1,"pcEvec1[3]/F");
    tree->Branch("pcEvec2",fPC2,"pcEvec2[3]/F");
    tree->Branch("pcEvec3",fPC3,"pcEvec3[3]/F");


    tree->Branch("pMCMom",fpMCMom,"pMCMom[4]/F");
    tree->Branch("pMCPos",fpMCPos,"pMCPos[4]/F");
    tree->Branch("pRECKalF",fpREC,"pRECKalF[4]/F");
    tree->Branch("pRECKalL",fpRECL,"pRECKalL[4]/F");
  

  //TGeoManager* geomGENFIT = new TGeoManager("Geometry", "Geane geometry");
  //TGeoManager::Import("config/genfitGeom.root");
  //  gROOT->Macro("config/Geane.C"); 
 
  }

//-------------------------------------------------
  void Track3DKalmanSPS::endJob()
  {
    if (!rep) delete rep;
    if (!repMC) delete repMC;

  /*
  //  not sure why I can't do these, but at least some cause seg faults.
  delete[] stMCT;
  delete[] covMCT;
  delete[] stREC;
  delete[] covREC;
  */

    delete[] fpREC;
    delete[] fpRECL;
    delete[] fState0;
    delete[] fCov0;

    delete[] fshx;
    delete[] fshy;
    delete[] fshz;
    delete[] feshx;
    delete[] feshy;
    delete[] feshyz;
    delete[] feshz;
    delete[] fupdate;
    delete[] fchi2hit;
    delete[] fth;
    delete[] feth;
    delete[] fedudw;
    delete[] fedvdw;
    delete[] feu;
    delete[] fev;
    delete[] fsep;
    delete[] fdQdx;

    delete[] fPCmeans;
    delete[] fPCsigmas;
    delete[] fPCevals;
    delete[] fPC1;
    delete[] fPC2;
    delete[] fPC3;
    
  }


//------------------------------------------------------------------------------------//
void Track3DKalmanSPS::produce(art::Event& evt)
{ 

  rep=0;
  repMC=0;
  // get services
  art::ServiceHandle<geo::Geometry> geom;

  //////////////////////////////////////////////////////
  // Make a std::unique_ptr<> for the thing you want to put into the event
  // because that handles the memory management for you
  //////////////////////////////////////////////////////
  std::unique_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);
  std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>); 
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> > thassn(new art::Assns<recob::Track, recob::Hit>); 

  unsigned int tcnt = 0;

  // define TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();


  // get input Hit object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  art::Handle< std::vector< art::PtrVector < recob::SpacePoint > > > spptListHandle;
  evt.getByLabel(fSpptModuleLabel,spptListHandle);

  art::PtrVector<simb::MCTruth> mclist;

  /// \todo Should never test whether the event is real data in reconstruction algorithms
  /// \todo as that introduces potential data/MC differences that are very hard to track down
  /// \todo Remove this test as soon as possible please
  if (!evt.isRealData()){

      //      std::cout << "Track3DKalmanSPS: This is MC." << std::endl;
      // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;

      art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
      evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);

      for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
	{
	  art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
	  mclist.push_back(mctparticle);
	}

    }


  // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;

  // Put this back when Wes's reign of terror ends ...
  //  LOG_DEBUG("Track3DKalmanSPS") << "There are " <<  spptListHandle->size() << " Spacepoint PtrVectors (spacepoint clumps) in this event.";

  std::vector < art::PtrVector<recob::SpacePoint> > spptIn(spptListHandle->begin(),spptListHandle->end());
  // Get the spptvectors that are largest to be first, and smallest last.
  std::sort(spptIn.begin(), spptIn.end(), sp_sort_nsppts);
 
  std::vector < art::PtrVector<recob::SpacePoint> >::const_iterator sppt = spptIn.begin();
  auto spptB = sppt;


  TVector3 MCOrigin;
  TVector3 MCMomentum;
  // TVector3 posErr(.05,.05,.05); // resolution. 0.5mm
  // TVector3 momErr(.1,.1,0.2);   // GeV
  TVector3 posErr(fPosErr[0],fPosErr[1],fPosErr[2]); // resolution. 0.5mm
  TVector3 momErr(fMomErr[0],fMomErr[1],fMomErr[2]);   // GeV
  TVector3 momErrFit(fMomErr[0],fMomErr[1],fMomErr[2]);   // GeV

  // This is strictly for MC
  /// \todo Should never test whether the event is real data in reconstruction algorithms
  /// \todo as that introduces potential data/MC differences that are very hard to track down
  /// \todo Remove this test as soon as possible please
  if (!evt.isRealData())
    {
      // Below breaks are stupid, I realize. But rather than keep all the MC
      // particles I just take the first primary, e.g., muon and only keep its
      // info in the branches of the Ttree. I could generalize later, ...
      for( unsigned int ii = 0; ii < mclist.size(); ++ii )
	{
	    //art::Ptr<const simb::MCTruth> mc(mctruthListHandle,i);
	    art::Ptr<simb::MCTruth> mc(mclist[ii]);
	    for(int jj = 0; jj < mc->NParticles(); ++jj)
	      {
		simb::MCParticle part(mc->GetParticle(jj));
		MCOrigin.SetXYZ(part.Vx(),part.Vy(),part.Vz()); // V for Vertex
		MCMomentum.SetXYZ(part.Px(),part.Py(),part.Pz());
		LOG_DEBUG("Track3DKalmanSPS_GenFit")
		  << "FROM MC TRUTH, the particle's pdg code is: "<<part.PdgCode()<< " with energy = "<<part.E() <<", with energy = "<<part.E()
		  << "\n  vtx: " << genf::ROOTobjectToString(MCOrigin)
		  << "\n  momentum: " << genf::ROOTobjectToString(MCMomentum)
		  << "\n    (both in Global (not volTPC) coords)";
		
		repMC = new genf::RKTrackRep(MCOrigin,
					     MCMomentum,
					     posErr,
					     momErr,
					     part.PdgCode());  
		break;
	      }
	    break;
	  }
	  //for saving of MC truth    
	  stMCT->ResizeTo(repMC->getState());
	  *stMCT = repMC->getState();
	  covMCT-> ResizeTo(repMC->getCov());
	  *covMCT = repMC->getCov();
	  LOG_DEBUG("Track3DKalmanSPS_GenFit") << " repMC, covMC are ... \n"
	    << genf::ROOTobjectToString(repMC->getState())
	    << genf::ROOTobjectToString(repMC->getCov());

	} // !isRealData
      nTrks = 0;
      TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(fPdg);
      Double_t mass = part->Mass();


      size_t cntp(0);
      while (sppt!=spptIn.end()) 
	{

	  const art::PtrVector<recob::SpacePoint> & spacepoints = *sppt;

	  double fMaxUpdateHere(fMaxUpdateU);
	  int fDecimateHere(fDecimateU);
	  double fErrScaleSHere(fErrScaleS);
	  double fErrScaleMHere(fErrScaleM);
	  int rePass0(1);
	  unsigned int nTailPoints = 0; // 100;
	  if (spacepoints.size()<5) 
	    { sppt++; rePass0 = 3; continue;} // for now...
		  
	  LOG_DEBUG("Track3DKalmanSPS_GenFit")
	    <<"\n\t found "<<spacepoints.size()<<" 3D spacepoint(s) for this element of std::vector<art:PtrVector> spacepoints. \n";
	  
	  //const double resolution = posErr.Mag(); 
	  //	  

	  
	  // Let's find track's principle components.
	  // We will sort along that direction, rather than z.
	  // Further, we will skip outliers away from main axis.
	  
	  TPrincipal* principal = new TPrincipal(3,"ND");
	  
	  // I need to shuffle these around, so use copy constructor
	  // to make non-const version spacepointss.
	  art::PtrVector<recob::SpacePoint> spacepointss(spacepoints);

	  // What I need is a nearest neighbor sorting.
	  if (fSortDim.compare("y") && fSortDim.compare("x")) std::sort(spacepointss.begin(), spacepointss.end(), sp_sort_3dz);
	  if (!fSortDim.compare("y")) std::sort(spacepointss.begin(), spacepointss.end(), sp_sort_3dy);
	  if (!fSortDim.compare("x")) std::sort(spacepointss.begin(), spacepointss.end(), sp_sort_3dx);

	  for (unsigned int point=0;point<spacepointss.size();++point)
	    {
	      //	      std::cout << "Spacepoint " << point << " added:" << spacepointss[point]->XYZ()[0]<< ", " << spacepointss[point]->XYZ()[1]<< ", " << spacepointss[point]->XYZ()[2]<< ". " << std::endl;
	      if (point<(spacepointss.size()-nTailPoints))
		{
		  principal->AddRow(spacepointss[point]->XYZ());
		}
	    }
	  principal->MakePrincipals();
	  /*
	    principal->Test();
	    principal->MakeHistograms();
	    principal->Print("MSEV");
	  */
	  const TVectorD* evals = principal->GetEigenValues(); 
	  const TMatrixD* evecs = principal->GetEigenVectors();
	  const TVectorD* means = principal->GetMeanValues();
	  const TVectorD* sigmas = principal->GetSigmas();
	/*
	  std::vector<TVector3*> pcs;
	  Double_t* pc = new Double_t[3];
	  for (unsigned int point=0;point<spacepointss.size();++point)
	    {
	      principal->X2P((Double_t *)(spacepointss[point]->XYZ()),pc);
	      pcs.push_back((TVector3 *)pc); // !!!
	    }
	  delete [] pc;
	*/
	  Double_t tmp[3], tmp2[3];
	  principal->X2P((Double_t *)(means->GetMatrixArray()),tmp);
	  principal->X2P((Double_t *)(sigmas->GetMatrixArray()),tmp2);
	  for (unsigned int ii=0;ii<3;++ii)
	    {
	      fPCmeans[ii] = (Float_t )(tmp[ii]);
	      fPCsigmas[ii] = (Float_t )(tmp2[ii]);
	      fPCevals[ii] = (Float_t )(evals->GetMatrixArray())[ii];
	      // This method requires apparently pulling all 9
	      // elements. Maybe 3 works. 
	      // Certainly, w can't be a scalar, I discovered.
	      double w[9];
	      evecs->ExtractRow(ii,0,w);
	      fPC1[ii] = w[0];
	      fPC2[ii] = w[1];
	      fPC3[ii] = w[2];
	    }
	  Double_t tmp3[3], tmp4[3], tmp5[3];
	  principal->X2P((Double_t *)fPC1,tmp3);
	  principal->X2P((Double_t *)fPC2,tmp4);
	  principal->X2P((Double_t *)fPC3,tmp5);


	  // Use a mip approximation assuming straight lines
	  // and a small angle wrt beam. 
	  fMomStart[0] = spacepointss[spacepointss.size()-1]->XYZ()[0] - spacepointss[0]->XYZ()[0];
	  fMomStart[1] = spacepointss[spacepointss.size()-1]->XYZ()[1] - spacepointss[0]->XYZ()[1];
	  fMomStart[2] = spacepointss[spacepointss.size()-1]->XYZ()[2] - spacepointss[0]->XYZ()[2];
	  // This presumes a 0.8 GeV/c particle
	  double dEdx = energyLossBetheBloch(mass, 1.0);
	  // mom is really KE. 
	  TVector3 mom(dEdx*fMomStart[0],dEdx*fMomStart[1],dEdx*fMomStart[2]);
	  double pmag2 = pow(mom.Mag()+mass, 2. - mass*mass);
	  mom.SetMag(std::sqrt(pmag2));
	  // Over-estimate by just enough for contained particles (5%).
	  mom.SetMag(1.0 * mom.Mag()); 
	  // My true 0.5 GeV/c muons need a yet bigger over-estimate.
	  //if (mom.Mag()<0.7) mom.SetMag(1.2*mom.Mag());  
	  //	  if (mom.Mag()>2.0) mom.SetMag(10.0*mom.Mag());  
	  //	  mom.SetMag(3*mom.Mag()); // EC, 15-Feb-2012. TEMPORARY!!!
	  // If 1st/last point is close to edge of TPC, this track is 
	  // uncontained.Give higher momentum starting value in 
	  // that case.
	  bool uncontained(false);
	  double close(5.); // cm. 
	  double epsMag(0.001);// cm. 
	  double epsX(250.0);  // cm. 
	  double epsZ(0.001);  // cm. 

	  if (
	      spacepointss[spacepointss.size()-1]->XYZ()[0] > (2.*geom->DetHalfWidth(0,0)-close) || spacepointss[spacepointss.size()-1]->XYZ()[0] < close ||
	      spacepointss[0]->XYZ()[0] > (2.*geom->DetHalfWidth(0,0)-close) || spacepointss[0]->XYZ()[0] < close ||
	      spacepointss[spacepointss.size()-1]->XYZ()[1] > (1.*geom->DetHalfHeight(0,0)-close) || (spacepointss[spacepointss.size()-1]->XYZ()[1] < -1.*geom->DetHalfHeight(0,0)+close) ||
	      spacepointss[0]->XYZ()[1] > (1.*geom->DetHalfHeight(0,0)-close) || spacepointss[0]->XYZ()[1] < (-1.*geom->DetHalfHeight(0,0)+close) ||
	      spacepointss[spacepointss.size()-1]->XYZ()[2] > (geom->DetLength(0,0)-close) || spacepointss[spacepointss.size()-1]->XYZ()[2] < close ||
	      spacepointss[0]->XYZ()[2] > (geom->DetLength(0,0)-close) || spacepointss[0]->XYZ()[2] < close
	      )
	    uncontained = true; 
	  fErrScaleSHere = fErrScaleS;
	  fErrScaleMHere = fErrScaleM;
	  
	  if (uncontained) 
	    {		      
	      // Big enough to not run out of gas right at end of
	      // track and give large angular deviations which
	      // will kill the fit.
	      mom.SetMag(2.0 * mom.Mag()); 
	      LOG_DEBUG("Track3DKalmanSPS_GenFit")<<"Uncontained track ... ";
	      fDecimateHere = fDecimateU;
	      fMaxUpdateHere = fMaxUpdateU;
	    }
	  else
	    {
	      LOG_DEBUG("Track3DKalmanSPS_GenFit")<<"Contained track ... Run "<<evt.run()<<" Event "<<evt.id().event();
	      // Don't decimate contained tracks as drastically, 
	      // and omit only very large corrections ...
	      // which hurt only high momentum tracks.
	      fDecimateHere = fDecimate;
	      fMaxUpdateHere = fMaxUpdate;
	    }
	  fcont = (int) (!uncontained);

	  // This seems like best place to jump back to for a re-pass.
	  unsigned short rePass = rePass0; // 1 by default; 
	  unsigned short maxPass(fMaxPass);
	  unsigned short tcnt1(0);
	  while (rePass<=maxPass)
	    {

	      TVector3 momM(mom);
	      TVector3 momErrFit(momM[0]/3.0,
				 momM[1]/3.0,
				 momM[2]/3.0);   // GeV
	  
	      genf::GFFieldManager::getInstance()->init(new genf::GFConstField(0.0,0.0,0.0));
	      genf::GFDetPlane planeG((TVector3)(spacepointss[0]->XYZ()),momM);
	  

	      //      std::cout<<"Track3DKalmanSPS about to do GAbsTrackRep."<<std::endl;
	      // Initialize with 1st spacepoint location and ...
	      rep = new genf::RKTrackRep(//posM-.5/momM.Mag()*momM,
					 (TVector3)(spacepointss[0]->XYZ()),
					 momM,
					 posErr,
					 momErrFit,
					 fPdg);  // mu+ hypothesis
	      //      std::cout<<"Track3DKalmanSPS: about to do GFTrack. repDim is " << rep->getDim() <<std::endl;
	  
	  
	      genf::GFTrack fitTrack(rep);//initialized with smeared rep
	      fitTrack.setPDG(fPdg);
	      // Gonna sort in z cuz I want to essentially transform here to volTPC coords.
	      // volTPC coords, cuz that's what the Geant3/Geane stepper wants, as that's its understanding
	      // from the Geant4 geometry, which it'll use. EC, 7-Jan-2011.
	      int ihit = 0;
	      fptsNo = 0;
	      std::vector <unsigned int> spptSurvivedIndex; 
	      std::vector <unsigned int> spptSkippedIndex; 
	      unsigned int ppoint(0);
	      for (unsigned int point=0;point<spacepointss.size();++point)
		{
		  double sep;
		  // Calculate the distance in 2nd and 3rd PCs and
		  // reject spt if it's too far out. Remember, the 
		  // sigmas are std::sqrt(eigenvals).
		  double tmp[3];
		  principal->X2P((Double_t *)(spacepointss[point]->XYZ()),tmp);
		  sep = std::sqrt(tmp[1]*tmp[1]/fPCevals[1]+tmp[2]*tmp[2]/fPCevals[2]);
		  if ((std::abs(sep) > fPerpLim) && (point<(spacepointss.size()-nTailPoints)) && rePass<=1)
		    {
		      //		      std::cout << "Spacepoint " << point << " DROPPED, cuz it's sufficiently far from the PCA major axis!!!:" << spacepointss[point]->XYZ()[0]<< ", " << spacepointss[point]->XYZ()[1]<< ", " << spacepointss[point]->XYZ()[2]<< ". " << std::endl;
		      spptSkippedIndex.push_back(point);
		      continue;
		    }
		  // If point is too close in Mag or Z or too far in X from last kept point drop it.
		  // I think this is largely redundant with PCA cut.
		  TVector3 one(spacepointss[point]->XYZ());
		  TVector3 two(spacepointss[ppoint]->XYZ());
		  if (rePass==2 && uncontained) 
		    {
		      epsMag = fDistanceU; // cm
		      fNumIt = 4;
		      fErrScaleMHere = 0.1; 
		      // Above allows us to pretend as though measurements 
		      // are perfect, which we can ostensibly do now with 
		      // clean set of sppts. This creates larger gains, bigger
		      // updates: bigger sensitivity to multiple scattering.

		      //		      std::cout << "Spacepoint " << point << " ?DROPPED? magnitude and TV3 diff to ppoint is :" << (((TVector3)(spacepointss[point]->XYZ()-spacepointss[ppoint]->XYZ())).Mag()) << " and " << one[0] << ", " << one[1] << ", " << one[2] << two[0] << ", " << two[1] << ", " << two[2] << ". " << std::endl;
		    }
		  else if (rePass==2 && !uncontained) 
		    {

		      //		      fNumIt = 2;
		      //		      std::cout << "Spacepoint " << point << " ?DROPPED? magnitude and TV3 diff to ppoint is :" << (((TVector3)(spacepointss[point]->XYZ()-spacepointss[ppoint]->XYZ())).Mag()) << " and " << one[0] << ", " << one[1] << ", " << one[2] << two[0] << ", " << two[1] << ", " << two[2] << ". " << std::endl;
		    }
		  if (point>0 && 
		      (
		       (one-two).Mag()<epsMag || // too close
		       ((one-two).Mag()>8.0&&rePass==1) || // too far
		       std::abs(spacepointss[point]->XYZ()[2]-spacepointss[ppoint]->XYZ()[2])<epsZ || 
		       std::abs(spacepointss[point]->XYZ()[0]-spacepointss[ppoint]->XYZ()[0])>epsX  
		       )
		      )
		    {
		      //		      std::cout << "Spacepoint " << point << " DROPPED, cuz it's too far in x or too close in magnitude or z to previous used spacepoint!!!:" << spacepointss[point]->XYZ()[0]<< ", " << spacepointss[point]->XYZ()[1]<< ", " << spacepointss[point]->XYZ()[2]<< ". " << std::endl;
		      //		      std::cout << "Prev used Spacepoint " << spacepointss[ppoint]->XYZ()[0]<< ", " << spacepointss[ppoint]->XYZ()[1]<< ", " << spacepointss[ppoint]->XYZ()[2]<< ". " << std::endl;
		      spptSkippedIndex.push_back(point);
		      continue;
		    }
		  
		  if (point%fDecimateHere && rePass<=1) // Jump out of loop except on every fDecimate^th pt. fDecimate==1 never sees continue. 
		    {
		      /* Replace continue with a counter that will be used
			 to index into vector of GFKalman fits.
		      */
		      // spptSkippedIndex.push_back(point);
		      continue;
		    }

		  ppoint=point;
		  TVector3 spt3 = (TVector3)(spacepointss[point]->XYZ());
		  std::vector <double> err3;
		  err3.push_back(spacepointss[point]->ErrXYZ()[0]);
		  err3.push_back(spacepointss[point]->ErrXYZ()[2]);
		  err3.push_back(spacepointss[point]->ErrXYZ()[4]);
		  err3.push_back(spacepointss[point]->ErrXYZ()[5]); // lower triangle diags.
		  if (fptsNo<fDimSize)
		    {
		      fshx[fptsNo] = spt3[0];
		      fshy[fptsNo] = spt3[1];
		      fshz[fptsNo] = spt3[2];
		      feshx[fptsNo] = err3[0];
		      feshy[fptsNo] = err3[1];
		      feshz[fptsNo] = err3[3];
		      feshyz[fptsNo] = err3[2];
		      fsep[fptsNo] = sep;

		      if (fptsNo>1)
			{
			  TVector3 pointer(fshx[fptsNo]-fshx[fptsNo-1],fshy[fptsNo]-fshy[fptsNo-1],fshz[fptsNo]-fshz[fptsNo-1]);
			  TVector3 pointerPrev(fshx[fptsNo-1]-fshx[fptsNo-2],fshy[fptsNo-1]-fshy[fptsNo-2],fshz[fptsNo-1]-fshz[fptsNo-2]);
			  fth[fptsNo] = (pointer.Unit()).Angle(pointerPrev.Unit());
			}
		      feth[fptsNo] = 0.0;
		      fedudw[fptsNo] = 0.0;
		      fedvdw[fptsNo] = 0.0;
		      feu[fptsNo] = 0.0;
		      fev[fptsNo] = 0.0;
		      fupdate[fptsNo] = 0.0;
		    }
	      
	      
		  LOG_DEBUG("Track3DKalmanSPS_GenFit") << "ihit xyz..." << spt3[0]<<","<< spt3[1]<<","<< spt3[2];

		  fitTrack.addHit(new genf::PointHit(spt3,err3),
				  1,//dummy detector id
				  ihit++
				  );
		  spptSurvivedIndex.push_back(point);
		  fptsNo++;
		} // end loop over spacepoints.
	  
	      if (fptsNo<=fMinNumSppts) // Cuz 1st 2 in each direction don't count. Should have, say, 3 more.
		{ 
		  LOG_DEBUG("Track3DKalmanSPS_GenFit") << "Bailing cuz only " << fptsNo << " spacepoints.";
		  rePass++;
		  continue;
		} 
	      LOG_DEBUG("Track3DKalmanSPS_GenFit") << "Fitting on " << fptsNo << " spacepoints.";
	      //      std::cout<<"Track3DKalmanSPS about to do GFKalman."<<std::endl;
	      genf::GFKalman k;
	      k.setBlowUpFactor(5); // 500 out of box. EC, 6-Jan-2011.
	      k.setMomHigh(fMomHigh); // Don't fit above this many GeV.
	      k.setMomLow(fMomLow);   // Don't fit below this many GeV.
	  
	      k.setInitialDirection(+1); // Instead of 1 out of box. EC, 6-Jan-2011.
	      k.setNumIterations(fNumIt);
	      k.setMaxUpdate(fMaxUpdateHere); // 0 out abs(update) bigger than this.		  
	      k.setErrorScaleSTh(fErrScaleSHere);
	      k.setErrorScaleMTh(fErrScaleMHere);
	  
	      bool skipFill = false;
	      //      std::cout<<"Track3DKalmanSPS back from setNumIterations."<<std::endl;
	      std::vector < TMatrixT<double> > hitMeasCov;
	      std::vector < TMatrixT<double> > hitUpdate;
	      std::vector < TMatrixT<double> > hitCov;
	      std::vector < TMatrixT<double> > hitCov7x7;
	      std::vector < TMatrixT<double> > hitState;
	      std::vector < double >           hitChi2;
	      std::vector <TVector3> hitPlaneXYZ;
	      std::vector <TVector3> hitPlaneUxUyUz;
	      std::vector <TVector3> hitPlaneU;
	      std::vector <TVector3> hitPlaneV;
	  
	      try{
		//	std::cout<<"Track3DKalmanSPS about to processTrack."<<std::endl;
		if (fDoFit) k.processTrack(&fitTrack);
		//std::cout<<"Track3DKalmanSPS back from processTrack."<<std::endl;
	      }
	      //catch(GFException& e){
	      catch(cet::exception &e){
		LOG_ERROR("Track3DKalmanSPS") << "just caught a cet::exception: " << e.what()
		  << "\nExceptions won't be further handled; skip filling big chunks of the TTree.";
		skipFill = true;
		//	exit(1);
	      }
	  
	      if(rep->getStatusFlag()==0) // 0 is successful completion
		{
		  if(mf::isDebugEnabled()) {
		    
		    std::ostringstream dbgmsg;
		    dbgmsg << "Original plane:";
		    planeG.Print(dbgmsg);
		    
		    dbgmsg << "Current (fit) reference Plane:";
		    rep->getReferencePlane().Print(dbgmsg);
		    
		    dbgmsg << "Last reference Plane:";
		    rep->getLastPlane().Print(dbgmsg);
		    
		    if (planeG != rep->getReferencePlane())
		      dbgmsg <<"  => original hit plane (not surprisingly) not current reference Plane!";
		    
		    LOG_DEBUG("Track3DKalmanSPS_GenFit") << dbgmsg.str();
		  }
		  if (!skipFill)
		    {
		      hitMeasCov = fitTrack.getHitMeasuredCov();
		      hitUpdate = fitTrack.getHitUpdate();
		      hitCov = fitTrack.getHitCov();
		      hitCov7x7 = fitTrack.getHitCov7x7();
		      hitState = fitTrack.getHitState();
		      hitChi2 = fitTrack.getHitChi2();
		      hitPlaneXYZ = fitTrack.getHitPlaneXYZ();
		      hitPlaneUxUyUz = fitTrack.getHitPlaneUxUyUz();
		      hitPlaneU = fitTrack.getHitPlaneU();
		      hitPlaneV = fitTrack.getHitPlaneV();
		      unsigned int totHits = hitState.size(); 
		  
		      //		  for (unsigned int ihit=0; ihit<fptsNo; ihit++)
		      // Pick up info from last fwd Kalman pass.
		      unsigned int jhit=0;
		      for (unsigned int ihit=totHits-2*totHits/(2*fNumIt); ihit<(totHits-totHits/(2*fNumIt)); ihit++) // was ihit<ihit<(totHits-fptsNo)<7
			{
			  feth[jhit] = (Float_t ) (hitMeasCov.at(ihit)[0][0]); // eth
			  fedudw[jhit] = (Float_t ) (hitMeasCov.at(ihit)[1][1]); 
			  fedvdw[jhit] = (Float_t ) (hitMeasCov.at(ihit)[2][2]); 
			  feu[jhit] = (Float_t ) (hitMeasCov.at(ihit)[3][3]); 
			  fev[jhit] = (Float_t ) (hitMeasCov.at(ihit)[4][4]);
			  fupdate[jhit] = (Float_t ) (hitUpdate.at(ihit)[0][0]);
			  fchi2hit[jhit] = (Float_t ) (hitChi2.at(ihit));
			  jhit++;
			}

		      stREC->ResizeTo(rep->getState());
		      *stREC = rep->getState();
		      covREC->ResizeTo(rep->getCov());
		      *covREC = rep->getCov();
		      double dum[5];
		      double dum2[5];
		      for (unsigned int ii=0;ii<5;ii++)
			{
			  stREC->ExtractRow(ii,0,dum);
			  fState0[ii] = dum[0];
			  covREC->ExtractRow(ii,0,dum2);
			  for (unsigned int jj=0;jj<5;jj++)
			    {
			      fCov0[ii*5+jj] = dum2[jj];
			    }
			}
		      LOG_DEBUG("Track3DKalmanSPS_GenFit")
		        << " First State and Cov:" << genf::ROOTobjectToString(*stREC)
		        << genf::ROOTobjectToString(*covREC);
		      chi2 = (Float_t)(rep->getChiSqu());
		      ndf = rep->getNDF();
		      nfail = fitTrack.getFailedHits();
		      nchi2rePass = (int)rePass;
		      ispptvec=1+std::distance(spptB, sppt);
		      chi2ndf = (Float_t)(chi2/ndf);
		  
		      nTrks++;
		      LOG_DEBUG("Track3DKalmanSPS_GenFit") << "Track3DKalmanSPS about to do tree->Fill(). Chi2/ndf is " << chi2/ndf << ".";
		      fpMCMom[3] = MCMomentum.Mag();
		      for (int ii=0;ii<3;++ii)
			{
			  fpMCMom[ii] = MCMomentum[ii];
			  fpMCPos[ii] = MCOrigin[ii];
			  fpREC[ii]   = hitPlaneUxUyUz.at(totHits-2*totHits/(2*fNumIt))[ii];
			  fpRECL[ii]  = hitPlaneUxUyUz.at(totHits-totHits/(2*fNumIt)-1)[ii];
			}
		  	      
		      evtt = (unsigned int) evt.id().event();
		      nspptvec = (unsigned int)  spptListHandle->size();

		      cntp++;
		      std::vector < std::vector <double> > dQdx;
		      // Calculate LastFwdPass quantities.
		      std::vector < TMatrixT<double> > hitCovLFP;
		      std::vector <TVector3> hitPlaneXYZLFP;
		      std::vector <TVector3> hitPlaneUxUyUzLFP;
		      std::vector <TVector3> hitPlanePxPyPzLFP;
		      std::vector <TVector3> hitPlaneULFP;
		      std::vector <TVector3> hitPlaneVLFP;
		      std::vector <double> pLFP;
		      std::vector < TMatrixT<double> > c7x7LFP;

		      art::FindManyP<recob::Hit> hitAssns(spacepointss, evt, fSpptModuleLabel);
		      for (unsigned int ii=0; ii<totHits/(2*fNumIt); ii++)
			{
			  pLFP.push_back(1./hitState.at(totHits-2*totHits/(2*fNumIt)+ii)[0][0]);
			  // hitCov -> hitCov7x7 !! EC, 11-May-2012.
			  c7x7LFP.push_back(hitCov7x7.at(totHits-2*totHits/(2*fNumIt)+ii));
			  hitCovLFP.push_back(hitCov.at(totHits-2*totHits/(2*fNumIt)+ii));
			  hitPlaneXYZLFP.push_back(hitPlaneXYZ.at(totHits-2*totHits/(2*fNumIt)+ii));
			  hitPlaneUxUyUzLFP.push_back(hitPlaneUxUyUz.at(totHits-2*totHits/(2*fNumIt)+ii));
			  hitPlanePxPyPzLFP.push_back(hitPlaneUxUyUzLFP.back()*pLFP.back());
			  hitPlaneULFP.push_back(hitPlaneU.at(totHits-2*totHits/(2*fNumIt)+ii));
			  hitPlaneVLFP.push_back(hitPlaneV.at(totHits-2*totHits/(2*fNumIt)+ii));
			  // Transform cov appropriate for track rotated 
			  // about w, forcing  
			  // v to be in y-z plane and u pointing in 
			  // +-ive x direction, per TrackAna convention.

			  rotationCov(hitCovLFP.back(),
				      hitPlaneULFP.back(),
				      hitPlaneVLFP.back()
				      );
			  dQdx.push_back(dQdxCalc(hitAssns,
						  spacepointss,
						  hitPlaneUxUyUzLFP.back(),
						  hitPlaneXYZLFP.back()
						  )
					 );
			  fdQdx[ii] = dQdx.back().back();

			}
		      fpREC[3]  = rep->getMom(rep->getReferencePlane()).Mag();
		      fpRECL[3] = pLFP[1];

		      tree->Fill();
		      
		      // Put newest track on stack for this set of sppts,
		      // remove previous one.
		      recob::tracking::SMatrixSym55 covVtx, covEnd;
		      for (unsigned int i=0; i<5; i++) {
			for (unsigned int j=i; j<5; j++) {
			  covVtx(i,j) = hitCovLFP.front()(i,j);
			  covEnd(i,j) = hitCovLFP.back()(i,j);
			}
		      }
		      recob::Track  the3DTrack(recob::TrackTrajectory(recob::tracking::convertCollToPoint(hitPlaneXYZLFP),
								      recob::tracking::convertCollToVector(hitPlanePxPyPzLFP),
								      recob::Track::Flags_t(hitPlaneXYZLFP.size()), true),
					       0, -1., 0, covVtx, covEnd, tcnt++);
		      if (rePass==1) tcnt1++; // won't get here if Trackfit failed.
		      if (rePass!=1 && tcnt1) tcol->pop_back();
		      tcol->push_back(the3DTrack); 
		      util::CreateAssn(*this, evt, *tcol, spacepointss, *tspassn);
		      art::PtrVector<recob::Hit> hits;// = hitAssns;
		      for (unsigned int ii=0; ii < spacepointss.size(); ++ii)
			{
			  //			  for (unsigned int jj=0; jj < hitAssns.at(ii).size(); ++jj)
			    //			    hits.push_back(hitAssns.at(ii).at(jj));
			  hits.insert(hits.end(),hitAssns.at(ii).begin(),hitAssns.at(ii).end());
			}
		      util::CreateAssn(*this, evt, *tcol, hits, *thassn, tcol->size()-1);
		    } // end !skipFill
		} // getStatusFlag
	  

	      if (!rep) delete rep;
	      rePass++;
	      // need to first excise bad spacepoints. 
	      // Grab up large Chi2hits first.
	      art::PtrVector<recob::SpacePoint> spacepointssExcise;
	      for (unsigned int ind=0;ind<spptSurvivedIndex.size();++ind)
		{
		  // Stricter to chuck sppts from uncontained then contained trks.
		  if ((uncontained&&fchi2hit[ind] >fChi2Thresh)     || 
		      (!uncontained&&fchi2hit[ind]>1.e9) || 
		      fchi2hit[ind]<0.0 
		      // =0 eliminates ruled-out large updates. Not obviously
		      // helpful.
		      // add a restriction on dQdx here ...
		      ) 
		    {
		      art::PtrVector<recob::SpacePoint>::iterator spptIt = spacepointss.begin()+spptSurvivedIndex[ind];
		      spacepointssExcise.push_back(*spptIt);
		    }
		}
	      // Now grab up those sppts which we skipped and don't want
	      // to reconsider cuz they're too close to each other or
	      // cuz they're too far in x, e.g.
	      for (unsigned int ind=0;ind<spptSkippedIndex.size();++ind)
		{
		  art::PtrVector<recob::SpacePoint>::iterator spptIt = spacepointss.begin()+spptSkippedIndex[ind];
		  spacepointssExcise.push_back(*spptIt);

		}
	      // Get rid of redundantly Excised sppts before proceeding.
	      std::stable_sort(spacepointss.begin(),spacepointss.end());
	      std::stable_sort(spacepointssExcise.begin(),spacepointssExcise.end());
	      std::set_union(spacepointssExcise.begin(),spacepointssExcise.end(),
			     spacepointssExcise.begin(),spacepointssExcise.end(),
			     spacepointssExcise.begin()
			     );
	      // Now excise. New spacepointss will be smaller for second pass.
	      art::PtrVector<recob::SpacePoint>::iterator diffSpptIt =
	      std::set_difference(spacepointss.begin(),spacepointss.end(),
				  spacepointssExcise.begin(),spacepointssExcise.end(),
				  spacepointss.begin()
				  );
	      spacepointss.erase(diffSpptIt,spacepointss.end());

	      // calculate new seed momentum, and errors as merited
	      if (rePass==2/* && uncontained */)
		{
		  if (fpREC[3]<fMomHigh && fpREC[3]>fMomLow)
		    {
		      double kick(0.9); //Try to get away with a smaller start
		      // for contained tracks. While for uncontained tracks
		      // let's start up at a higher momentum and come down.
		      // Though, 2 (1) GeV/c tracks are too low (high), so
		      // instead let's actually lower starting value on
		      // this second pass. -- EC 7-Mar-2013
		      if  (uncontained) kick = 0.5;
		      for (int ii=0;ii<3;++ii)
			{
			  //mom[ii] = fpREC[ii]*fpREC[3]*kick;
			  mom[ii] = momM[ii]*kick;
			}
		    }
		  else if (uncontained)
		    {
		      double unstick(1.0);
		      if  (fpREC[3]>=fMomHigh) unstick = 0.3;
		      for (int ii=0;ii<3;++ii)
			{
			  mom[ii] = momM[ii]*unstick;
			}
		    }
		  else 
		      for (int ii=0;ii<3;++ii)
			{
			  mom[ii] = 1.1*momM[ii];
			}
		    
		}

	    } // end while rePass<=maxPass

	  sppt++;
	  
	} // end while on elements of std::vector of art::PtrVectors of SpPts.
      
      if (!repMC) delete repMC;
      
      evt.put(std::move(tcol)); 
      // and now the spacepoints
      evt.put(std::move(tspassn));
      // and the hits. Note that these are all the hits from all the spacepoints considered,
      // even though they're not all contributing to the tracks.
      evt.put(std::move(thassn));
}

  DEFINE_ART_MODULE(Track3DKalmanSPS)

}  // end namespace

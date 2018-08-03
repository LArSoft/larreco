// -------------------------------------------------
// makes shower profile templates
//
// Author: Rory Fitzpatrick (roryfitz@umich.edu)
// Created: 8/3/18
// -------------------------------------------------

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/MCBase/MCDataHolder.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TROOT.h>
#include <TStyle.h>

namespace shower {

  class TCShowerElectronLikelihood : public art::EDAnalyzer {

  public:

    explicit TCShowerElectronLikelihood(fhicl::ParameterSet const& pset);
    virtual ~TCShowerElectronLikelihood();

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void analyze(const art::Event& evt);
   
    void getShowerProfile(std::vector< art::Ptr<recob::Hit> > showerhits, TVector3 shwvtx, TVector3 shwdir);
    int findEnergyBin();
 
  protected:

  private: 

    std::string fTemplateFile;
    std::string fROOTfile;

    TH3F* longTemplate;
    TH3F* tranTemplate;

    //TTree* fTree;

    TH1F* longProfile;
    TH1F* tranProfile;

    const int LBINS = 20;
    const int LMIN = 0;
    const int LMAX = 5;

    const int TBINS = 20;
    const int TMIN = -5;
    const int TMAX = 5;

    const int EBINS = 10;
    const double EMIN = 0.5;
    const double EMAX = 10.5;

    const double X0 = 14;

    std::string fHitModuleLabel;
    std::string fShowerModuleLabel;
    std::string fTemplateModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fDigitModuleLabel;

    calo::CalorimetryAlg fCalorimetryAlg;

  }; // class TCShowerElectronLikelihood

} // shower

// -------------------------------------------------

shower::TCShowerElectronLikelihood::TCShowerElectronLikelihood(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fHitModuleLabel           (pset.get< std::string >("HitModuleLabel", "trajcluster" ) ),
  fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel", "tcshower" ) ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel", "generator") ),
  fCalorimetryAlg           (pset.get< fhicl::ParameterSet >("CalorimetryAlg") ) {
  this->reconfigure(pset);
} // TCShowerElectronLikelihood

// -------------------------------------------------

shower::TCShowerElectronLikelihood::~TCShowerElectronLikelihood() {
} // ~TCShowerTemplateMaker

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::reconfigure(fhicl::ParameterSet const& pset) {
  fTemplateFile           = pset.get< std::string >("TemplateFile");
  cet::search_path sp("FW_SEARCH_PATH");
  if( !sp.find_file(fTemplateFile, fROOTfile) )
    throw cet::exception("TCShowerElectronLikelihood") << "cannot find the root template file: \n" 
						       << fTemplateFile
						       << "\n bail ungracefully.\n";

  TFile *file = TFile::Open(fROOTfile.c_str());

  longTemplate = (TH3F*)file->Get("tcshowertemplate/fLongitudinal");
  tranTemplate = (TH3F*)file->Get("tcshowertemplate/fTransverse");

} // reconfigure

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::beginJob() {

  art::ServiceHandle<art::TFileService> tfs;
  //fTree = tfs->make<TTree>("tcshowerana", "tcshowerana");

  longProfile = tfs->make<TH1F>("longProfile", "longitudinal shower profile;t;Q", LBINS, LMIN, LMAX);
  tranProfile = tfs->make<TH1F>("tranProfile", "transverse shower profile;dist (cm);Q", TBINS, TMIN, TMAX);

} // beginJob

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::analyze(const art::Event& evt) {

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if (evt.getByLabel(fShowerModuleLabel,showerListHandle))
    art::fill_ptr_vector(showerlist, showerListHandle);

  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  art::FindManyP<recob::Hit> shwfm(showerListHandle, evt, fShowerModuleLabel);

  if (showerlist.size()) {
    std::vector< art::Ptr<recob::Hit> > showerhits = shwfm.at(0);
    getShowerProfile(showerhits, showerlist[0]->ShowerStart(), showerlist[0]->Direction());
    findEnergyBin();

    // check true shower energy
    if (mclist.size()) {
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
      if (mctruth->NeutrinoSet()) {
	if (std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 12 && mctruth->GetNeutrino().CCNC() == 0) {
	  double elep =  mctruth->GetNeutrino().Lepton().E();
	  std::cout << "true shower energy: " << elep << std::endl;
	}
      }
    }
  }

  //fTree->Fill();

} // analyze

// -------------------------------------------------

void shower::TCShowerElectronLikelihood::getShowerProfile(std::vector< art::Ptr<recob::Hit> > showerhits, TVector3 shwvtx, TVector3 shwdir) {

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;

  auto collectionPlane = geo::PlaneID(0, 0, 1);

  double shwVtxTime = detprop->ConvertXToTicks(shwvtx[0], collectionPlane);
  double shwVtxWire = geom->WireCoordinate(shwvtx[1], shwvtx[2], collectionPlane);

  double shwTwoTime = detprop->ConvertXToTicks(shwvtx[0]+shwdir[0], collectionPlane);
  double shwTwoWire = geom->WireCoordinate(shwvtx[1]+shwdir[1], shwvtx[2]+shwdir[2], collectionPlane);

  for (size_t i = 0; i < showerhits.size(); ++i) {
    if (showerhits[i]->WireID().Plane != collectionPlane.Plane) continue;

    double wirePitch = geom->WirePitch(showerhits[i]->WireID());
    double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns

    double xvtx = shwVtxTime * tickToDist;
    double yvtx = shwVtxWire * wirePitch;

    double xtwo = shwTwoTime * tickToDist;
    double ytwo = shwTwoWire * wirePitch;

    double xtwoorth = (ytwo - yvtx) + xvtx;
    double ytwoorth = -(xtwo - xvtx) + yvtx;

    double xhit = showerhits[i]->PeakTime() * tickToDist;
    double yhit = showerhits[i]->WireID().Wire * wirePitch;

    double ldist = std::abs((ytwoorth-yvtx)*xhit - (xtwoorth-xvtx)*yhit + xtwoorth*yvtx - ytwoorth*xvtx)/std::sqrt( pow((ytwoorth-yvtx), 2) + pow((xtwoorth-xvtx), 2) );
    double tdist = ((ytwo-yvtx)*xhit - (xtwo-xvtx)*yhit + xtwo*yvtx - ytwo*xvtx)/std::sqrt( pow((ytwo-yvtx), 2) + pow((xtwo-xvtx), 2) );

    double to3D = 1. / sqrt( pow(xvtx-xtwo,2) + pow(yvtx-ytwo,2) ) ; // distance between two points in 3D space is one 
    ldist *= to3D;
    tdist *= to3D;

    double Q = showerhits[i]->Integral() * fCalorimetryAlg.LifetimeCorrection(showerhits[i]->PeakTime());

    longProfile->Fill(ldist/X0, Q);
    tranProfile->Fill(tdist, Q);

  } // loop through showerhits

  return;

} // getShowerProfile

// -------------------------------------------------

int shower::TCShowerElectronLikelihood::findEnergyBin() {

  // confirm equal number of distance bins in templates and candidates

  // make TProfile for each energy bin from fTemplate histograms

  // calculate chi2 between fProfile and fTemplate for each energy bin (add longitudinal + transverse, normalize by dof)

  // pick energy bin with lowest chi2

  // return bin number

  return -1;

} // findEnergyBin

// -------------------------------------------------

DEFINE_ART_MODULE(shower::TCShowerElectronLikelihood)

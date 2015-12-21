////////////////////////////////////////////////////////////////////////
//
// Chi2PIDAlg class
//
// tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

//#include "RecoBase/Track.h"
#include "lardata/AnalysisBase/Calorimetry.h"
#include "lardata/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/AnalysisAlg/Chi2PIDAlg.h"

// ROOT includes
#include "TFile.h"
#include "TProfile.h"

// Framework includes
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

//------------------------------------------------------------------------------
pid::Chi2PIDAlg::Chi2PIDAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

//------------------------------------------------------------------------------
pid::Chi2PIDAlg::~Chi2PIDAlg()
{
}

//------------------------------------------------------------------------------
void pid::Chi2PIDAlg::reconfigure(fhicl::ParameterSet const& pset)
{
  fTemplateFile           = pset.get< std::string >("TemplateFile");
  //fCalorimetryModuleLabel = pset.get< std::string >("CalorimetryModuleLabel");

  cet::search_path sp("FW_SEARCH_PATH");
 
  if( !sp.find_file(fTemplateFile, fROOTfile) )
    throw cet::exception("Chi2ParticleID") << "cannot find the root template file: \n" 
					   << fTemplateFile
					   << "\n bail ungracefully.\n";
  TFile *file = TFile::Open(fROOTfile.c_str());
  dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
  dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
  dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
  dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

  return;
}


//------------------------------------------------------------------------------
void pid::Chi2PIDAlg::DoParticleID(art::Ptr<anab::Calorimetry> calo,
				  anab::ParticleID &pidOut){
  int npt = 0;
  double chi2pro = 0;
  double chi2ka = 0;
  double chi2pi = 0;
  double chi2mu = 0;
  double trkpitchc = calo->TrkPitchC();
  double avgdedx = 0;
  double PIDA = 0; //by Bruce Baller
  std::vector<double> trkdedx = calo->dEdx();
  std::vector<double> trkres = calo->ResidualRange();
  std::vector<double> deadwireresrc = calo->DeadWireResRC();
  pidOut.fPlaneID = calo->PlaneID();

  int used_trkres = 0;
  for (unsigned i = 0; i<trkdedx.size(); ++i){//hits
    avgdedx += trkdedx[i];
    if(trkres[i] < 30) {
      PIDA += trkdedx[i]*pow(trkres[i],0.42);
      used_trkres++;
    }
    if (trkdedx[i]>1000) continue; //protect against large pulse height
    int bin = dedx_range_pro->FindBin(trkres[i]);
    if (bin>=1&&bin<=dedx_range_pro->GetNbinsX()){
      double bincpro = dedx_range_pro->GetBinContent(bin);
      if (bincpro<1e-6){//for 0 bin content, using neighboring bins
	bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2;
      }
      double bincka = dedx_range_ka->GetBinContent(bin);
      if (bincka<1e-6){
	bincka = (dedx_range_ka->GetBinContent(bin-1)+dedx_range_ka->GetBinContent(bin+1))/2;
      }
      double bincpi = dedx_range_pi->GetBinContent(bin);
      if (bincpi<1e-6){
	bincpi = (dedx_range_pi->GetBinContent(bin-1)+dedx_range_pi->GetBinContent(bin+1))/2;
      }
      double bincmu = dedx_range_mu->GetBinContent(bin);
      if (bincmu<1e-6){
	bincmu = (dedx_range_mu->GetBinContent(bin-1)+dedx_range_mu->GetBinContent(bin+1))/2;
      }
      double binepro = dedx_range_pro->GetBinError(bin);
      if (binepro<1e-6){
	binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2;
      }
      double bineka = dedx_range_ka->GetBinError(bin);
      if (bineka<1e-6){
	bineka = (dedx_range_ka->GetBinError(bin-1)+dedx_range_ka->GetBinError(bin+1))/2;
      }
      double binepi = dedx_range_pi->GetBinError(bin);
      if (binepi<1e-6){
	binepi = (dedx_range_pi->GetBinError(bin-1)+dedx_range_pi->GetBinError(bin+1))/2;
      }
      double binemu = dedx_range_mu->GetBinError(bin);
      if (binemu<1e-6){
	binemu = (dedx_range_mu->GetBinError(bin-1)+dedx_range_mu->GetBinError(bin+1))/2;
      }
      //double errke = 0.05*trkdedx[i];   //5% KE resolution
      double errdedx = 0.04231+0.0001783*trkdedx[i]*trkdedx[i]; //resolution on dE/dx
      errdedx *= trkdedx[i];
      chi2pro += pow((trkdedx[i]-bincpro)/std::sqrt(pow(binepro,2)+pow(errdedx,2)),2);
      chi2ka += pow((trkdedx[i]-bincka)/std::sqrt(pow(bineka,2)+pow(errdedx,2)),2);
      chi2pi += pow((trkdedx[i]-bincpi)/std::sqrt(pow(binepi,2)+pow(errdedx,2)),2);
      chi2mu += pow((trkdedx[i]-bincmu)/std::sqrt(pow(binemu,2)+pow(errdedx,2)),2);
      //std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
      ++npt;
    }
  }

  //anab::ParticleID pidOut;
  if (npt){
    pidOut.fNdf        = npt;
    pidOut.fChi2Proton = chi2pro/npt;
    pidOut.fChi2Kaon   = chi2ka/npt;
    pidOut.fChi2Pion   = chi2pi/npt;
    pidOut.fChi2Muon   = chi2mu/npt;
    double chi2[4] = {chi2pro/npt,chi2ka/npt,chi2pi/npt,chi2mu/npt};
    double pdg[4] = {2212,321,211,13};
    double chi2min = 1e20;
    int imin = -1;
    double chi2min2 = 1e20;
    //int imin2;
    // find the minimal chi2 and next-to-minimal chi2
    for (int ichi2 = 0; ichi2<4; ++ichi2){
      if (chi2[ichi2]<chi2min){
	imin = ichi2;
	chi2min2 = chi2min;
	chi2min = chi2[ichi2];
      }
      else if (chi2[ichi2]<chi2min2){
	//imin2 = ichi2;
	chi2min2 = chi2[ichi2];
      }
    }
    if (imin>-1){
      pidOut.fPdg = pdg[imin];
      pidOut.fMinChi2 = chi2min;
      pidOut.fDeltaChi2 = chi2min2 - chi2min;
    }
  }
  //if (trkdedx.size()) pidOut.fPIDA = PIDA/trkdedx.size();
  if(used_trkres > 0) pidOut.fPIDA = PIDA/used_trkres;
  double missinge = 0;
  double missingeavg = 0;
  for (unsigned i = 0; i<deadwireresrc.size(); ++i){
    int bin = dedx_range_pro->FindBin(deadwireresrc[i]);
    //std::cout<<i<<" "<<deadwireresrc[i]<<" "<<bin<<std::endl;
    if (bin<1) continue;
    if (bin>dedx_range_pro->GetNbinsX()) bin = dedx_range_pro->GetNbinsX();
    if (pidOut.fPdg==2212){
      missinge += dedx_range_pro->GetBinContent(bin)*trkpitchc;
    }
    else if (pidOut.fPdg==321){
      missinge += dedx_range_ka->GetBinContent(bin)*trkpitchc;
    }
    else if (pidOut.fPdg==211){
      missinge += dedx_range_pi->GetBinContent(bin)*trkpitchc;
    }
    else if (pidOut.fPdg==13){
      missinge += dedx_range_mu->GetBinContent(bin)*trkpitchc;
    }
    //std::cout<<bin<<" "<<dedx_range_pro->GetBinContent(bin)*trkpitchc<<std::endl;
  }
  if (trkdedx.size()) missingeavg = avgdedx/trkdedx.size()*trkpitchc*deadwireresrc.size();
  //std::cout<<trkIter<<" "<<pid<<std::endl;
  pidOut.fMissingE = missinge;
  pidOut.fMissingEavg = missingeavg;
}


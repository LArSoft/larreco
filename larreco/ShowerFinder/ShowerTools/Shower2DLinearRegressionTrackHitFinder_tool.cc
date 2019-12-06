//############################################################################
//### Name:        Shower2DLinearRegressionTrackHitFinder                  ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk)         ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the initial shower track using a rms   ###
//###              based method to define when the shower starts to        ###
//###              shower. This methd is derived from the EMShower_module  ###
//############################################################################
#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/TRACSAlg.h"

//C++ Includes
#include <iostream>
#include <math.h>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools{

  class Shower2DLinearRegressionTrackHitFinder:IShowerTool {
  public:
    
    Shower2DLinearRegressionTrackHitFinder(const fhicl::ParameterSet& pset);
    
    ~Shower2DLinearRegressionTrackHitFinder();
    
    //Calculate the 2D initial track hits
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;
    
  private:
    
    //Function to find the 
    std::vector<art::Ptr<recob::Hit> > FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits);
    
    //Function to perform a weighted regression fit.
    Int_t WeightedFit(const Int_t n, const Double_t *x, const Double_t *y,
		      const Double_t *w,  Double_t *parm);
    
    //fcl parameters
    unsigned int               fNfitpass;           //Number of time to fit the straight
                                                    //line the hits.
    std::vector<unsigned int>  fNfithits;           //Max number of hits to fit to.
    std::vector<double>        fToler;              //Tolerance or each interaction. 
                                                    //Defined as the perpendicualar 
                                                    //distance from the best fit line.
    bool                       fApplyChargeWeight;  //Apply charge weighting to the fit.
    art::InputTag              fPFParticleModuleLabel;
    art::InputTag              fHitsModuleLabel;
    std::string                fShowerStartPositionInputLabel;
    std::string                fShowerDirectionInputLabel;
    std::string                fInitialTrackHitsOutputLabel;
    std::string                fInitialTrackSpacePointsOutputLabel;
  };
  

  Shower2DLinearRegressionTrackHitFinder::Shower2DLinearRegressionTrackHitFinder(
      const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fNfitpass(pset.get<unsigned int>("Nfitpass")),
    fNfithits(pset.get<std::vector<unsigned int> >("Nfithits")),
    fToler(pset.get<std::vector<double> >("Toler")),
    fApplyChargeWeight(pset.get<bool>("ApplyChargeWeight")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel")),
    fHitsModuleLabel(pset.get<art::InputTag>("HitsModuleLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel")),
    fInitialTrackHitsOutputLabel(pset.get<std::string>("InitialTrackHitsOutputLabel")),
    fInitialTrackSpacePointsOutputLabel(pset.get<std::string>("InitialTrackSpacePointsOutputLabel"))
  {
    if (fNfitpass!=fNfithits.size() ||
        fNfitpass!=fToler.size()) {
      throw art::Exception(art::errors::Configuration)
        << "Shower2DLinearRegressionTrackHitFinderEMShower: fNfithits and fToler need to have size fNfitpass";
    }
  }

  Shower2DLinearRegressionTrackHitFinder::~Shower2DLinearRegressionTrackHitFinder()
  {
  }

  int Shower2DLinearRegressionTrackHitFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      mf::LogError("Shower2DLinearRegressionTrackHitFinder")
        << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)){
      mf::LogError("Shower2DLinearRegressionTrackHitFinder")
        << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("Shower2DLinearRegressionTrackHitFinder") 
	<< "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("Shower2DLinearRegressionTrackHitFinder") 
	<< "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

    if(clusters.size()<2){
      mf::LogError("Shower2DLinearRegressionTrackHitFinder")
        << "Not enough clusters: "<<clusters.size() << std::endl;
      return 1;
    }

    //Get the hit association
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);
    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > plane_clusters;
    //Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){

      //Get the hits
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());

      for (auto hit : hits) {
        geo::WireID wire = hit->WireID();
        geo::PlaneID plane = wire.asPlaneID();
        plane_clusters[plane].push_back(hit);
      }

      // Was having issues with clusters having hits in multiple planes breaking PMA
      // So switched to the method above. May want to switch back when using PandoraTrack
      //plane_clusters[plane].insert(plane_clusters[plane].end(),hits.begin(),hits.end());
    }

    std::vector<art::Ptr<recob::Hit> > InitialTrackHits;
    //Loop over the clusters and order the hits and get the initial track hits in that plane
    for(auto const& cluster: plane_clusters){

      //Get the hits
      std::vector<art::Ptr<recob::Hit> > hits = cluster.second;

      //Order the hits
      IShowerTool::GetTRACSAlg().OrderShowerHits(hits,ShowerStartPosition,ShowerDirection);

      //Find the initial track hits
      std::vector<art::Ptr<recob::Hit> > trackhits = FindInitialTrackHits(hits);

      InitialTrackHits.insert(InitialTrackHits.end(),trackhits.begin(),trackhits.end());
    }

    //Holders for the initial track values.
    ShowerEleHolder.SetElement(InitialTrackHits, fInitialTrackHitsOutputLabel);

    //Get the associated spacepoints
    //Get the hits
    art::Handle<std::vector<recob::Hit> > hitHandle;
    if (!Event.getByLabel(fHitsModuleLabel, hitHandle)){
      throw cet::exception("Shower2DLinearRegressionTrackHitFinder")
        << "Could not get the hits." << std::endl;
      return 1;
    }

    //get the sp<->hit association
    art::FindManyP<recob::SpacePoint> fmsp(hitHandle,Event,fPFParticleModuleLabel);
    if(!fmsp.isValid()){
      throw cet::exception("Shower2DLinearRegressionTrackHitFinder")
        << "Spacepoint and hit association not valid. Stopping." << std::endl;
      return 1;
    }

    //Get the spacepoints associated to the track hit
    std::vector<art::Ptr<recob::SpacePoint > > intitaltrack_sp;
    for(auto const& hit: InitialTrackHits){
      std::vector<art::Ptr<recob::SpacePoint > > sps = fmsp.at(hit.key());
      for(auto const sp: sps){
        intitaltrack_sp.push_back(sp);
      }
    }
    ShowerEleHolder.SetElement(intitaltrack_sp, fInitialTrackSpacePointsOutputLabel);
    return 0;
  }

  //Function to calculate the what are the initial tracks hits. Adapted from EMShower FindInitialTrackHits
  std::vector<art::Ptr<recob::Hit> > Shower2DLinearRegressionTrackHitFinder::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits){

    std::vector<art::Ptr<recob::Hit> > trackHits;

    double parm[2];
    int fitok = 0;
    std::vector<double> wfit;
    std::vector<double> tfit;
    std::vector<double> cfit;

    for (size_t i = 0; i<fNfitpass; ++i){

      // Fit a straight line through hits
      unsigned int nhits = 0;
      for (auto &hit: hits){

        //Not sure I am a fan of doing things in wire tick space. What if id doesn't not iterate properly or the
        //two planes in each TPC are not symmetric.
        TVector2 coord = IShowerTool::GetTRACSAlg().HitCoordinates(hit);

        if (i==0||(std::abs((coord.Y()-(parm[0]+coord.X()*parm[1]))*cos(atan(parm[1])))<fToler[i-1])||fitok==1){
          ++nhits;
          if (nhits==fNfithits[i]+1) break;
          wfit.push_back(coord.X());
          tfit.push_back(coord.Y());

          if(fApplyChargeWeight) { cfit.push_back(hit->Integral());
          } else { cfit.push_back(1.); };
          if (i==fNfitpass-1) {
            trackHits.push_back(hit);
          }
        }
      }

      if (i<fNfitpass-1&&wfit.size()){
        fitok = WeightedFit(wfit.size(), &wfit[0], &tfit[0], &cfit[0], &parm[0]);
      }

      wfit.clear();
      tfit.clear();
      cfit.clear();
    }
    return trackHits;
  }

  //Stolen from EMShowerAlg, a linear regression fitting function
  Int_t Shower2DLinearRegressionTrackHitFinder::WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm){

    Double_t sumx=0.;
    Double_t sumx2=0.;
    Double_t sumy=0.;
    Double_t sumy2=0.;
    Double_t sumxy=0.;
    Double_t sumw=0.;
    Double_t eparm[2];

    parm[0]  = 0.;
    parm[1]  = 0.;
    eparm[0] = 0.;
    eparm[1] = 0.;

    for (Int_t i=0; i<n; i++) {
      sumx += x[i]*w[i];
      sumx2 += x[i]*x[i]*w[i];
      sumy += y[i]*w[i];
      sumy2 += y[i]*y[i]*w[i];
      sumxy += x[i]*y[i]*w[i];
      sumw += w[i];
    }

    if (sumx2*sumw-sumx*sumx==0.) return 1;
    if (sumx2-sumx*sumx/sumw==0.) return 1;

    parm[0] = (sumy*sumx2-sumx*sumxy)/(sumx2*sumw-sumx*sumx);
    parm[1] = (sumxy-sumx*sumy/sumw)/(sumx2-sumx*sumx/sumw);

    eparm[0] = sumx2*(sumx2*sumw-sumx*sumx);
    eparm[1] = (sumx2-sumx*sumx/sumw);

    if (eparm[0]<0. || eparm[1]<0.) return 1;

    eparm[0] = sqrt(eparm[0])/(sumx2*sumw-sumx*sumx);
    eparm[1] = sqrt(eparm[1])/(sumx2-sumx*sumx/sumw);

    return 0;

  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::Shower2DLinearRegressionTrackHitFinder)

////////////////////////////////////////////////////////////////////////
// Class:       EMShowerPDF
// Module Type: analyser
// File:        EMShowerPDF_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), July 2015
//
// Analyser module to calculate the PDFs needed to form the 3D merging
// likelihoods for EMShower reconstruction. 
////////////////////////////////////////////////////////////////////////

// Framework includes:
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
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FindManyP.h"
#include "MCCheater/BackTracker.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Shower.h"
#include "Utilities/AssociationUtil.h"

// ROOT includes
#include "TPrincipal.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TF2.h"
#include "TTree.h"
#include "TSpline.h"

namespace shower {
  class EMShowerPDF;
}

class shower::EMShowerPDF : public art::EDAnalyzer {
public:

  explicit EMShowerPDF(fhicl::ParameterSet const &pset);
  virtual ~EMShowerPDF();

  void analyze(art::Event const &evt);
  void reconfigure(fhicl::ParameterSet const &p);
  void beginJob();
  void endJob();

private:

  art::ServiceHandle<cheat::BackTracker> backtracker;
  std::map<int,int> trueTrackMap;

  art::ServiceHandle<art::TFileService> tfs;

  TH1D *hAngleMerge;
  TH1D *hEigenvalueMerge;
  TH1D *hAngEigenMerge;
  TH1D *hAngleNoMerge;
  TH1D *hEigenvalueNoMerge;
  TH1D *hAngEigenNoMerge;
  TH2D *hAngleTracksMerge;
  TH2D *hAngleTracksNoMerge;
  TH2D *hAngleSeparationMerge;
  TH2D *hAngleSeparationNoMerge;
  TH1D *hTubeDiameterMerge;
  TH1D *hTubeDiameterNoMerge;
  TH1D *hTubeDiameterNormMerge;
  TH1D *hTubeDiameterNormNoMerge;
  TH1D *hDOCAMerge;
  TH1D *hDOCANoMerge;
  TH1D *hLikelihoodMerge;
  TH1D *hLikelihoodNoMerge;

  TTree *fTree;
  double fAngle;
  double fLongTrackAngle;
  double fShortTrackAngle;
  double fLongTrackLength;
  double fShortTrackLength;
  double fSeparation;
  double fTubeDiameter;
  double fTubeDiameterNorm;
  double fTotalLength;
  double fDOCA;
  double fAngleProb;
  double fTubeProb;
  double fDOCAProb;
  double fTrueMerge;

  // PDFs
  TH1D *hAnglePDF;
  TH1D *hTubePDF;
  TH1D *hDOCAPDF;

  // Splines
  bool fMakeSplines;
  TString fSplineFileName = "EMSplines.root";
  TFile *fSplineFile;
  TSpline3 *fAngleSpline, *fTubeSpline, *fDOCASpline;

};

shower::EMShowerPDF::EMShowerPDF(fhicl::ParameterSet const &pset) : EDAnalyzer(pset) {
  this->reconfigure(pset);
}

shower::EMShowerPDF::~EMShowerPDF() { }

void shower::EMShowerPDF::beginJob() {

  gStyle->SetOptStat(0);

  hAngleMerge = new TH1D("AngleMerge","Angle",20,0,1.6);
  hEigenvalueMerge = new TH1D("EigenvalueMerge","Eigenvalue",100,0.99,1);
  hAngEigenMerge = new TH1D("AngEigenMerge","AngEigen",50,0,1.6);
  hAngleNoMerge = new TH1D("AngleNoMerge","Angle",20,0,1.6);
  hEigenvalueNoMerge = new TH1D("EigenvalueNoMerge","Eigenvalue",100,0.99,1);
  hAngEigenNoMerge = new TH1D("AngEigenNoMerge","AngEigen",50,0,1.6);
  hAngleTracksMerge = new TH2D("AngleTracksMerge",";Long track angle;Short track angle;",50,0,2,50,0,2);
  hAngleTracksNoMerge = new TH2D("AngleTracksNoMerge",";Long track angle;Short track angle;",50,0,2,50,0,2);
  hAngleSeparationMerge = new TH2D("AngleSeparationMerge",";Separation;Angle difference;",100,0,100,50,0,2);
  hAngleSeparationNoMerge = new TH2D("AngleSeparationNoMerge",";Separation;Angle difference;",100,0,100,50,0,2);
  hTubeDiameterMerge = new TH1D("TubeDiameterMerge",";Diameter;",10,0,10);
  hTubeDiameterNoMerge = new TH1D("TubeDiameterNoMerge",";Diameter;",10,0,10);
  hTubeDiameterNormMerge = new TH1D("TubeDiameterNormMerge",";Norm Diameter;",20,0,0.5);
  hTubeDiameterNormNoMerge = new TH1D("TubeDiameterNormNoMerge",";Norm Diameter;",20,0,0.5);
  hDOCAMerge = new TH1D("DOCAMerge",";DOCA;",40,0,100);
  hDOCANoMerge = new TH1D("DOCANoMerge",";DOCA;",40,0,100);
  hLikelihoodMerge = new TH1D("LikelihoodMerge",";Likelihood;",50,0,0.5);
  hLikelihoodNoMerge = new TH1D("LikelihoodNoMerge",";Likelihood;",50,0,0.5);

  fTree = tfs->make<TTree>("MatchingVariables","MatchingVariables");
  fTree->Branch("Angle",&fAngle);
  fTree->Branch("LongTrackAngle",&fLongTrackAngle);
  fTree->Branch("ShortTrackAngle",&fShortTrackAngle);
  fTree->Branch("LongTrackLength",&fLongTrackLength);
  fTree->Branch("ShortTrackLength",&fShortTrackLength);
  fTree->Branch("Separation",&fSeparation);
  fTree->Branch("TubeDiameter",&fTubeDiameter);
  fTree->Branch("TubeDiameterNorm",&fTubeDiameterNorm);
  fTree->Branch("TotalLength",&fTotalLength);
  fTree->Branch("DOCA",&fDOCA);
  fTree->Branch("AngleProb",&fAngleProb);
  fTree->Branch("TubeProb",&fTubeProb);
  fTree->Branch("DOCAProb",&fDOCAProb);
  fTree->Branch("TrueMerge",&fTrueMerge);

  // PDFs
  hAnglePDF = new TH1D("AnglePDF","Angle PDF",20,0,1.6);
  hTubePDF = new TH1D("TubePDF","Tube PDF",10,0,10);
  hDOCAPDF = new TH1D("DOCAPDF","DOCA PDF",40,0,100);

  // Splines
  if (!fMakeSplines) {
    fSplineFile = TFile::Open(fSplineFileName);
    fAngleSpline = (TSpline3*)fSplineFile->Get("AngleSpline");
    fTubeSpline = (TSpline3*)fSplineFile->Get("TubeSpline");
    fDOCASpline = (TSpline3*)fSplineFile->Get("DOCASpline");
  }

}

void shower::EMShowerPDF::reconfigure(fhicl::ParameterSet const &p) {
  fMakeSplines = false;
}

void shower::EMShowerPDF::analyze(art::Event const &evt) {

  // Output -- showers and associations with hits and clusters
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Cluster> > clusterAssociations(new art::Assns<recob::Shower, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitAssociations(new art::Assns<recob::Shower, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Shower, recob::SpacePoint> > spacePointAssociations(new art::Assns<recob::Shower, recob::SpacePoint>); 
  std::unique_ptr<art::Assns<recob::Shower, recob::Track> > trackAssociations(new art::Assns<recob::Shower, recob::Track>);

  // Get the tracks from the event record
  art::Handle<std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if (evt.getByLabel("emshower3d",trackHandle))
    art::fill_ptr_vector(tracks, trackHandle);

  // Get the cluster handle to find all associated hits
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  evt.getByLabel("blurredclustering",clusterHandle);

  // Find the associations for the tracks
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, "blurredclustering");
  art::FindManyP<recob::SpacePoint> fmsp(trackHandle, evt, "emshower3d");
  art::FindManyP<recob::Cluster> fmc(trackHandle, evt, "emshower3d");
  art::FindManyP<recob::Vertex> fmv(trackHandle, evt, "emshower3d");

  std::cout << std::endl << "Event " << evt.event() << " has " << tracks.size() << " tracks" << std::endl;

  for (unsigned int trackIt = 0; trackIt < tracks.size(); ++trackIt) {

    std::vector<art::Ptr<recob::SpacePoint> > points = fmsp.at(trackIt);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(trackIt);
    std::vector<art::Ptr<recob::Vertex> > vertices = fmv.at(trackIt);

    art::Ptr<recob::Track> track = tracks.at(trackIt);

    double XYZ[3];
    if (vertices.size() > 0) vertices.at(0)->XYZ(XYZ);

    // Use truth info to find which particle this track is associated with
    std::map<int,double> trackMap;

    for (auto &trackCluster : clusters) {

      std::vector<art::Ptr<recob::Hit> > clusterHits = fmh.at(trackCluster->ID());

      for (auto &clusterHit : clusterHits) {

	std::vector<sim::TrackIDE> ides = backtracker->HitToTrackID(clusterHit);

	for (auto &ide : ides)
	  trackMap[ide.trackID] += ide.energy;

      }
    
    }

    // Find the true particle associated with this track
    double highEnergy = 0;
    int bestTrack = 0;
    for (auto &track : trackMap) {
      if (track.second > highEnergy) {
	highEnergy = track.second;
	bestTrack = track.first;
      }
    }
    trueTrackMap[track->ID()] = bestTrack;

    // Dump some stuff
    std::cout << "Track " << track->ID() << " has " << track->NumberTrajectoryPoints() << " trajectory points and has start position ("
  	      << track->Vertex().X() << "," << track->Vertex().Y() << "," << track->Vertex().Z()
  	      << ") and direction ("
  	      << track->VertexDirection().X() << "," << track->VertexDirection().Y() << "," << track->VertexDirection().Z()
  	      << "), has " << vertices.size() << " vert(ex)(ices) and is associated with particle " << bestTrack << std::endl;
  }

  // Calculate the PDFs
  for (unsigned int track1It = 0; track1It < tracks.size(); ++track1It) {

    for (unsigned int track2It = track1It+1; track2It < tracks.size(); ++track2It) {

      // Get LArSoft objects
      art::Ptr<recob::Track> track1 = tracks.at(track1It);
      art::Ptr<recob::Track> track2 = tracks.at(track2It);

      // Find if these tracks correspond to same true particle
      if (trueTrackMap[track1->ID()] == trueTrackMap[track2->ID()])
	fTrueMerge = true;
      else fTrueMerge = false;

      // Find total hits for each track
      std::vector<art::Ptr<recob::Cluster> > clusters1 = fmc.at(track1It);
      std::vector<art::Ptr<recob::Cluster> > clusters2 = fmc.at(track2It);
      int nhits1 = 0, nhits2 = 0;
      for (auto const& cluster1 : clusters1) {
	std::vector<art::Ptr<recob::Hit> > hits1 = fmh.at(cluster1->ID());
        nhits1 += hits1.size();
      }
      for (auto const& cluster2 : clusters2) {
	std::vector<art::Ptr<recob::Hit> > hits2 = fmh.at(cluster2->ID());
        nhits2 += hits2.size();
      }

      // Angle between segments
      fAngle = track1->VertexDirection().Angle(track2->VertexDirection());
      if (fAngle > 1.57) fAngle = 3.14159 - fAngle;

      // Find the length of the tracks
      TVector3 lengthTrack1 = track1->Vertex() - track1->End();
      TVector3 lengthTrack2 = track2->Vertex() - track2->End();
      fLongTrackLength = std::max(lengthTrack1.Mag(), lengthTrack2.Mag());
      fShortTrackLength = std::min(lengthTrack1.Mag(), lengthTrack2.Mag());      
      fTotalLength = fLongTrackLength + fShortTrackLength;

      // Find the half point along each track
      TVector3 halfPointTrack1 = lengthTrack1 * 0.5;
      TVector3 halfPointTrack2 = lengthTrack2 * 0.5;
      TVector3 direction = halfPointTrack1 - halfPointTrack2;

      // Find info about this track
      fSeparation = direction.Mag();
      double angle1 = track1->VertexDirection().Angle(direction);
      double angle2 = track2->VertexDirection().Angle(direction);
      if (angle1 > 1.57) angle1 = 3.14159 - angle1;
      if (angle2 > 1.57) angle2 = 3.14159 - angle2;

      if (nhits1 > nhits2) {
	fLongTrackAngle = angle1;
	fShortTrackAngle = angle2;
      }
      else {
	fLongTrackAngle = angle2;
	fShortTrackAngle = angle1;
      }

      // Project onto plane
      TVector3 projection1 = lengthTrack1 - ((lengthTrack1.Dot(direction.Unit())) * direction.Unit());
      TVector3 projection2 = lengthTrack2 - ((lengthTrack2.Dot(direction.Unit())) * direction.Unit());

      // Find magnitude of this projection
      double magProjection1 = projection1.Mag();
      double magProjection2 = projection2.Mag();
      fTubeDiameter = magProjection1 > magProjection2 ? magProjection1 : magProjection2;
      fTubeDiameterNorm = (double)fTubeDiameter/(double)fTotalLength;

      // POCA/DOCA
      TVector3 direction1 = (track1->VertexDirection()).Unit();
      TVector3 direction2 = (track2->VertexDirection()).Unit();
      TVector3 position1 = track1->Vertex();
      TVector3 position2 = track2->Vertex();

      // The following follows from some vector algebra;
      // I've convinced myself it's correct
      double D = direction1.Dot(direction2);
      TVector3 P = position2 - position1;

      double parameter1 = (((D*direction1)-direction2).Dot(P)) * 1./(1 - D*D);
      double parameter2 = (parameter1 - direction1.Dot(P)) * (1./D);

      TVector3 closestPoint1 = position1 + parameter1 * direction1;
      TVector3 closestPoint2 = position2 + parameter2 * direction2;

      fDOCA = (closestPoint1 - closestPoint2).Mag();

      // Save stuff
      std::cout << "Tracks " << track1It << " and " << track2It << " Merge? " << fTrueMerge << ", angles " << angle1 << " and " << angle2 << ", and nhits " << nhits1 << " and " << nhits2 << std::endl;

      if (!fMakeSplines) {
	fAngleProb = fAngleSpline->Eval(fAngle);
	fTubeProb = fTubeSpline->Eval(fTubeDiameter);
	fDOCAProb = fDOCASpline->Eval(fDOCA);
	double likelihood = fAngleProb * fDOCAProb * fTubeProb;
	if (fTrueMerge) hLikelihoodMerge->Fill(likelihood);
	else hLikelihoodNoMerge->Fill(likelihood);
      }

      if (fTrueMerge) {
	hAngleMerge->Fill(fAngle);
	hAngleTracksMerge->Fill(fLongTrackAngle, fShortTrackAngle);
	hAngleSeparationMerge->Fill(fSeparation, std::abs(angle1 - angle2));
	hTubeDiameterMerge->Fill(fTubeDiameter);
	hTubeDiameterNormMerge->Fill((double)fTubeDiameter/(double)fTotalLength);
	hDOCAMerge->Fill(fDOCA);
      }

      else {
	hAngleNoMerge->Fill(fAngle);
	hAngleTracksNoMerge->Fill(fLongTrackAngle, fShortTrackAngle);
	hAngleSeparationNoMerge->Fill(fSeparation, std::abs(angle1 - angle2));
	hTubeDiameterNoMerge->Fill(fTubeDiameter);
	hTubeDiameterNormNoMerge->Fill((double)fTubeDiameter/(double)fTotalLength);
	hDOCANoMerge->Fill(fDOCA);
      }

      fTree->Fill();

    }
  }

  std::cout << std::endl;

}

void shower::EMShowerPDF::endJob() {

  TCanvas *cAngle = new TCanvas("cAngle","",800,600);
  //if (hAngleMerge->Integral() != 0) hAngleMerge->Scale(1/hAngleMerge->Integral());
  hAngleMerge->Draw();
  //if (hAngleNoMerge->Integral() != 0) hAngleNoMerge->Scale(1/hAngleNoMerge->Integral());
  hAngleNoMerge->SetLineColor(2);
  hAngleNoMerge->Draw("same");
  cAngle->SaveAs("Angle.png");

  TCanvas *cEigenvalue = new TCanvas("cEigenvalue","",800,600);
  if (hEigenvalueMerge->Integral() != 0) hEigenvalueMerge->Scale(1/hEigenvalueMerge->Integral());
  hEigenvalueMerge->Draw();
  if (hEigenvalueNoMerge->Integral() != 0) hEigenvalueNoMerge->Scale(1/hEigenvalueNoMerge->Integral());
  hEigenvalueNoMerge->SetLineColor(2);
  hEigenvalueNoMerge->Draw("same");
  cEigenvalue->SaveAs("Eigenvalue.png");

  TCanvas *cAngEigen = new TCanvas("cAngEigen","",800,600);
  if (hAngEigenMerge->Integral() != 0) hAngEigenMerge->Scale(1/hAngEigenMerge->Integral());
  hAngEigenMerge->Draw();
  if (hAngEigenNoMerge->Integral() != 0) hAngEigenNoMerge->Scale(1/hAngEigenNoMerge->Integral());
  hAngEigenNoMerge->SetLineColor(2);
  hAngEigenNoMerge->Draw("same");
  cAngEigen->SaveAs("AngEigen.png");

  TCanvas *cAngleTracks = new TCanvas("cAngleTracks","",800,600);
  hAngleTracksMerge->SetMarkerStyle(8);
  hAngleTracksMerge->SetMarkerSize(1.5);
  hAngleTracksMerge->Draw();
  hAngleTracksNoMerge->SetMarkerStyle(8);
  hAngleTracksNoMerge->SetMarkerSize(1.5);
  hAngleTracksNoMerge->SetMarkerColor(2);
  hAngleTracksNoMerge->Draw("same");
  cAngleTracks->SaveAs("AngleTracks.png");

  TCanvas *cAngleSeparation = new TCanvas("cAngleSeparation","",800,600);
  hAngleSeparationMerge->SetMarkerStyle(8);
  hAngleSeparationMerge->SetMarkerSize(1.5);
  hAngleSeparationMerge->Draw();
  hAngleSeparationNoMerge->SetMarkerColor(2);
  hAngleSeparationNoMerge->SetMarkerStyle(8);
  hAngleSeparationNoMerge->SetMarkerSize(1.5);
  hAngleSeparationNoMerge->Draw("same");
  cAngleSeparation->SaveAs("AngleSeparation.png");

  TCanvas *cTubeDiameter = new TCanvas("cTubeDiameter","",800,600);
  if (hTubeDiameterMerge->Integral() != 0) hTubeDiameterMerge->Scale(1/hTubeDiameterMerge->Integral());
  hTubeDiameterMerge->Draw();
  if (hTubeDiameterNoMerge->Integral() != 0) hTubeDiameterNoMerge->Scale(1/hTubeDiameterNoMerge->Integral());
  hTubeDiameterNoMerge->SetLineColor(2);
  hTubeDiameterNoMerge->Draw("same");
  cTubeDiameter->SaveAs("TubeDiameter.png");

  TCanvas *cTubeDiameterNorm = new TCanvas("cTubeDiameterNorm","",800,600);
  if (hTubeDiameterNormMerge->Integral() != 0) hTubeDiameterNormMerge->Scale(1/hTubeDiameterNormMerge->Integral());
  hTubeDiameterNormMerge->Draw();
  if (hTubeDiameterNormNoMerge->Integral() != 0) hTubeDiameterNormNoMerge->Scale(1/hTubeDiameterNormNoMerge->Integral());
  hTubeDiameterNormNoMerge->SetLineColor(2);
  hTubeDiameterNormNoMerge->Draw("same");
  cTubeDiameterNorm->SaveAs("TubeDiameterNorm.png");

  TCanvas *cDOCA = new TCanvas("cDOCA","",800,600);
  hDOCAMerge->Draw();
  hDOCANoMerge->SetLineColor(2);
  hDOCANoMerge->Draw("same");
  cDOCA->SaveAs("DOCA.png");

  // Make the PDF hists
  TH1D *hAngleSum = (TH1D*)hAngleMerge->Clone();
  hAngleSum->Add(hAngleNoMerge);
  hAnglePDF = (TH1D*)hAngleMerge->Clone();
  hAnglePDF->Divide(hAngleSum);

  TH1D *hTubeSum = (TH1D*)hTubeDiameterMerge->Clone();
  hTubeSum->Add(hTubeDiameterNoMerge);
  hTubePDF = (TH1D*)hTubeDiameterMerge->Clone();
  hTubePDF->Divide(hTubeSum);

  TH1D *hDOCASum = (TH1D*)hDOCAMerge->Clone();
  hDOCASum->Add(hDOCANoMerge);
  hDOCAPDF = (TH1D*)hDOCAMerge->Clone();
  hDOCAPDF->Divide(hDOCASum);

  // Fit splines
  TSpline3 sAnglePDF = TSpline3(hAnglePDF);
  TSpline3 sTubePDF = TSpline3(hTubePDF);
  TSpline3 sDOCAPDF = TSpline3(hDOCAPDF);
  if (fMakeSplines) {
    TFile *splineFile = TFile::Open(fSplineFileName,"RECREATE");
    sAnglePDF.SetName("AngleSpline");
    sTubePDF.SetName("TubeSpline");
    sDOCAPDF.SetName("DOCASpline");
    sAnglePDF.Write();
    sTubePDF.Write();
    sDOCAPDF.Write();
    splineFile->Close();
  }
  else {
    TCanvas *cLikelihood = new TCanvas("Likelihood","",800,600);
    hLikelihoodMerge->Draw();
    hLikelihoodNoMerge->SetLineColor(2);
    hLikelihoodNoMerge->Draw("same");
    cLikelihood->SaveAs("Likelihood.png");
  }

  TCanvas *cAnglePDF = new TCanvas("AnglePDF","",800,600);
  hAnglePDF->Draw();
  sAnglePDF.SetLineColor(2);
  sAnglePDF.Draw("same");
  cAnglePDF->SaveAs("AnglePDF.png");

  TCanvas *cTubePDF = new TCanvas("TubePDF","",800,600);
  hTubePDF->Draw();
  sTubePDF.SetLineColor(2);
  sTubePDF.Draw("same");
  cTubePDF->SaveAs("TubePDF.png");

  TCanvas *cDOCAPDF = new TCanvas("DOCAPDF","",800,600);
  hDOCAPDF->Draw();
  sDOCAPDF.SetLineColor(2);
  sDOCAPDF.Draw("same");
  cDOCAPDF->SaveAs("DOCAPDF.png");

  if (!fMakeSplines)
    fSplineFile->Close();

}


DEFINE_ART_MODULE(shower::EMShowerPDF)

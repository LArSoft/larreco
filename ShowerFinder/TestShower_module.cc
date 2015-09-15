////////////////////////////////////////////////////////////////////////
// Class:       TestShower
// Module Type: producer
// File:        TestShower_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), July 2015
//
// 
// 
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
#include "art/Framework/Core/EDProducer.h"
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

namespace shower {
  class TestShower;
}

class shower::TestShower : public art::EDProducer {
public:

  explicit TestShower(fhicl::ParameterSet const &pset);
  virtual ~TestShower();

  void produce(art::Event &evt);
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
  TH1D *hPOCAMerge;
  TH1D *hPOCANoMerge;

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
  double fPOCA;
  double fTrueMerge;

};

shower::TestShower::TestShower(fhicl::ParameterSet const &pset) {
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
  produces<art::Assns<recob::Shower, recob::Track> >();
}

shower::TestShower::~TestShower() { }

void shower::TestShower::beginJob() {

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
  hTubeDiameterMerge = new TH1D("TubeDiameterMerge",";Diameter;",20,0,10);
  hTubeDiameterNoMerge = new TH1D("TubeDiameterNoMerge",";Diameter;",20,0,10);
  hTubeDiameterNormMerge = new TH1D("TubeDiameterNormMerge",";Norm Diameter;",20,0,0.5);
  hTubeDiameterNormNoMerge = new TH1D("TubeDiameterNormNoMerge",";Norm Diameter;",20,0,0.5);
  hPOCAMerge = new TH1D("POCAMerge",";POCA;",100,0,100);
  hPOCANoMerge = new TH1D("POCANoMerge",";POCA;",100,0,100);

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
  fTree->Branch("POCA",&fPOCA);
  fTree->Branch("TrueMerge",&fTrueMerge);

}

void shower::TestShower::reconfigure(fhicl::ParameterSet const &p) {
}

void shower::TestShower::produce(art::Event &evt) {

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

    double highEnergy = 0;
    int bestTrack = 0;
    for (auto &track : trackMap) {
      if (track.second > highEnergy) {
	highEnergy = track.second;
	bestTrack = track.first;
      }
    }

    trueTrackMap[track->ID()] = bestTrack;

    std::cout << "Track " << track->ID() << " has " << track->NumberTrajectoryPoints() << " trajectory points and has start position ("
  	      << track->Vertex().X() << "," << track->Vertex().Y() << "," << track->Vertex().Z()
  	      << ") and direction ("
  	      << track->VertexDirection().X() << "," << track->VertexDirection().Y() << "," << track->VertexDirection().Z()
  	      << "), has " << vertices.size() << " vert(ex)(ices) and is associated with particle " << bestTrack << std::endl;
  }

  // std::vector<unsigned int> mergedTracks;

  // // Sort the tracks by size
  // std::sort(tracks.begin(), tracks.end(), [](const art::Ptr<recob::Track> &a, const art::Ptr<recob::Track> &b) {return a->NumberTrajectoryPoints() > b->NumberTrajectoryPoints();} );

  // // Until all tracks are merged, create new showers
  // bool mergedAllTracks = false;
  // if (!tracks.size()) mergedAllTracks = true;
  // while (!mergedAllTracks) {

  //   // New shower and associations
  //   art::PtrVector<recob::Hit> showerHits;
  //   art::PtrVector<recob::Cluster> showerClusters;
  //   art::PtrVector<recob::Track> showerTracks;
  //   art::PtrVector<recob::SpacePoint> showerSpacePoints;

  //   // Put the largest unmerged track into this shower
  //   for (unsigned int initTrack = 0; initTrack < tracks.size(); ++initTrack) {
  //     if (std::find(mergedTracks.begin(), mergedTracks.end(), initTrack) != mergedTracks.end()) continue;
  //     mergedTracks.push_back(initTrack);
  //     showerTracks.push_back(tracks.at(initTrack));
  //     for (auto &showerCluster : fmc.at(initTrack)) {
  // 	showerClusters.push_back(showerCluster);
  // 	for (auto &showerHit : fmh.at(showerCluster->ID()))
  // 	  showerHits.push_back(showerHit);
  //     }
  //     for (auto &showerSpacePoint : fmsp.at(initTrack))
  // 	showerSpacePoints.push_back(showerSpacePoint);
  //     break;
  //   }

  //   // Find all aligned tracks to this
  //   bool mergedAllToThisTrack = false;
  //   while (!mergedAllToThisTrack) {

  //     // Look at all tracks
  //     int nadded = 0;
  //     for (unsigned int trialTrack = 0; trialTrack < tracks.size(); ++trialTrack) {

  // 	if (std::find(mergedTracks.begin(), mergedTracks.end(), trialTrack) != mergedTracks.end()) continue;

  // 	// Calculate the PCA for each
  // 	TPrincipal *pca = new TPrincipal(3,"");

  // 	for (auto &showerSpacePoint : showerSpacePoints)
  // 	  pca->AddRow(showerSpacePoint->XYZ());

  // 	for (auto &trialTrackSpacePoint : fmsp.at(trialTrack))
  // 	  pca->AddRow(trialTrackSpacePoint->XYZ());

  // 	pca->MakePrincipals();

  // 	std::cout << "Eigenvalue for trial track " << tracks.at(trialTrack)->ID() << " is " << (*pca->GetEigenValues())[0] << std::endl;

  // 	if ((*pca->GetEigenValues())[0] > 0.9) {

  // 	  showerTracks.push_back(tracks.at(trialTrack));
  // 	  for (auto &showerCluster : fmc.at(trialTrack)) {
  // 	    showerClusters.push_back(showerCluster);
  // 	    for (auto &showerHit : fmh.at(showerCluster->ID()))
  // 	      showerHits.push_back(showerHit);
  // 	  }
  // 	  for (auto &showerSpacePoint : fmsp.at(trialTrack))
  // 	    showerSpacePoints.push_back(showerSpacePoint);	  

  // 	  mergedTracks.push_back(trialTrack);
  // 	  ++nadded;

  // 	}

  // 	delete pca;

  //     } // loop over tracks to add

  //     if (nadded == 0) mergedAllToThisTrack = true;

  //   } // while loop

  //   // Total charge
  //   std::map<int,double> charge;
  //   for (auto &hit : showerHits)
  //     charge[hit->View()] += (hit->SummedADC() * std::exp((hit->PeakTime() * 500)/(3000000)));

  //   std::vector<double> v = {0.0};
  //   std::vector<double> energy = {charge[0],charge[1],charge[2]};
  //   TVector3 v3 = TVector3(0,0,0);
  //   TVector3 VertexDirection = showerTracks.at(0)->VertexDirection();
  //   TVector3 Vertex = showerTracks.at(0)->Vertex();
  //   showers->emplace_back(VertexDirection, v3, Vertex, v3, energy, v, v, v, 0, showers->size());
  //   if (mergedTracks.size() == tracks.size()) mergedAllTracks = true;

  //   util::CreateAssn(*this, evt, *(showers.get()), showerHits,        *(hitAssociations.get()));
  //   util::CreateAssn(*this, evt, *(showers.get()), showerClusters,    *(clusterAssociations.get()));
  //   util::CreateAssn(*this, evt, *(showers.get()), showerTracks,      *(trackAssociations.get()));
  //   util::CreateAssn(*this, evt, *(showers.get()), showerSpacePoints, *(spacePointAssociations.get()));

  // } // end creating showers

  std::vector<int> mergedTracks;

  // Sort the tracks by size
  std::sort(tracks.begin(), tracks.end(), [](const art::Ptr<recob::Track> &a, const art::Ptr<recob::Track> &b) {return a->NumberTrajectoryPoints() > b->NumberTrajectoryPoints();} );

  for (unsigned int track = 0; track < tracks.size(); ++track) {

    if (std::find(mergedTracks.begin(), mergedTracks.end(), tracks.at(track)->ID()) != mergedTracks.end())
      continue;

    mergedTracks.push_back(tracks.at(track)->ID());

    for (unsigned int newTrack = 0; newTrack < tracks.size(); ++newTrack) {

      if (std::find(mergedTracks.begin(), mergedTracks.end(), tracks.at(newTrack)->ID()) != mergedTracks.end())
	continue;

      if (trueTrackMap[tracks.at(track)->ID()] == trueTrackMap[tracks.at(newTrack)->ID()])
	fTrueMerge = true;
      else fTrueMerge = false;

      TPrincipal *pca = new TPrincipal(3,"");
      TGraph2D *graph = new TGraph2D();

      for (auto &trackSpacePoint : fmsp.at(track)) {
	pca->AddRow(trackSpacePoint->XYZ());
	graph->SetPoint(graph->GetN(), trackSpacePoint->XYZ()[0], trackSpacePoint->XYZ()[1], trackSpacePoint->XYZ()[2]);
	//std::cout << "Points... " << trackSpacePoint->XYZ()[0] << " " << trackSpacePoint->XYZ()[1] << " " << trackSpacePoint->XYZ()[2] << std::endl;
      }

      for (auto &newTrackSpacePoint : fmsp.at(newTrack)) {
	pca->AddRow(newTrackSpacePoint->XYZ());
	graph->SetPoint(graph->GetN(), newTrackSpacePoint->XYZ()[0], newTrackSpacePoint->XYZ()[1], newTrackSpacePoint->XYZ()[2]);
	//std::cout << "More points... " << newTrackSpacePoint->XYZ()[0] << " " << newTrackSpacePoint->XYZ()[1] << " " << newTrackSpacePoint->XYZ()[2] << std::endl;
      }

      pca->MakePrincipals();

      // graph->Fit("pol1","QC0");
      // //TF1 *fit = graph->GetFunction("pol1");
      // //double gradient = std::abs(fit->GetParameter(1));
      // TCanvas *c = new TCanvas("c","",800,600);
      // graph->Draw();
      // c->SaveAs("test.png");
      // graph->SaveAs("TestGraph");
      // delete c;
      // std::cin.get();
      // //delete fit;
      // //delete graph;

      double eigenvalue = (*pca->GetEigenValues())[0];
      double angle = tracks.at(track)->VertexDirection().Angle(tracks.at(newTrack)->VertexDirection());
      if (angle > 1.57) angle = 3.14159 - angle;

      std::cout << "Tracks " << tracks.at(track)->ID() << " and " << tracks.at(newTrack)->ID() << ": angle is " << angle << " and eigenvalue is " << eigenvalue << ". Should they be merged? " << fTrueMerge << std::endl;

      if (fTrueMerge) {
	//hAngleMerge->Fill(angle);
	hEigenvalueMerge->Fill(eigenvalue);
	hAngEigenMerge->Fill(angle * eigenvalue);
      }
      else {
	//hAngleNoMerge->Fill(angle);
	hEigenvalueNoMerge->Fill(eigenvalue);
	hAngEigenNoMerge->Fill(angle * eigenvalue);
      }

    }

  }

  // Look at making new angle plots
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

      double D = direction1.Dot(direction2);
      TVector3 P = position2 - position1;

      double parameter1 = (((D*direction1)-direction2).Dot(P)) * 1./(1 - D*D);
      double parameter2 = (parameter1 - direction1.Dot(P)) * (1./D);

      TVector3 closestPoint1 = position1 + parameter1 * direction1;
      TVector3 closestPoint2 = position2 + parameter2 * direction2;

      fPOCA = (closestPoint1 - closestPoint2).Mag();

      // Save stuff
      std::cout << "Tracks " << track1It << " and " << track2It << " Merge? " << fTrueMerge << ", angles " << angle1 << " and " << angle2 << ", and nhits " << nhits1 << " and " << nhits2 << std::endl;

      fTree->Fill();

      if (fTrueMerge) {
	hAngleMerge->Fill(fAngle);
	hAngleTracksMerge->Fill(fLongTrackAngle, fShortTrackAngle);
	hAngleSeparationMerge->Fill(fSeparation, std::abs(angle1 - angle2));
	hTubeDiameterMerge->Fill(fTubeDiameter);
	hTubeDiameterNormMerge->Fill((double)fTubeDiameter/(double)fTotalLength);
	hPOCAMerge->Fill(fPOCA);
      }

      else {
	hAngleNoMerge->Fill(fAngle);
	hAngleTracksNoMerge->Fill(fLongTrackAngle, fShortTrackAngle);
	hAngleSeparationNoMerge->Fill(fSeparation, std::abs(angle1 - angle2));
	hTubeDiameterNoMerge->Fill(fTubeDiameter);
	hTubeDiameterNormNoMerge->Fill((double)fTubeDiameter/(double)fTotalLength);
	hPOCANoMerge->Fill(fPOCA);
      }

    }
  }

  // Put in event

  std::cout << "Before merging, there were " << tracks.size() << " tracks; there are " << showers->size() << " showers" << std::endl;

  evt.put(std::move(showers));
  evt.put(std::move(hitAssociations));
  evt.put(std::move(clusterAssociations));
  evt.put(std::move(trackAssociations));
  evt.put(std::move(spacePointAssociations));

  std::cout << std::endl;

}

void shower::TestShower::endJob() {

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

  TCanvas *cPOCA = new TCanvas("cPOCA","",800,600);
  hPOCAMerge->Draw();
  hPOCANoMerge->SetLineColor(2);
  hPOCANoMerge->Draw("same");
  cPOCA->SaveAs("POCA.png");

}


DEFINE_ART_MODULE(shower::TestShower)

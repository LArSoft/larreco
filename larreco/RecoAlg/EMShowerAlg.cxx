///////////////////////////////////////////////////////////////////
// Implementation of the EMShower algorithm
//
// Forms EM showers from clusters and associated tracks.
// Also provides methods for finding the vertex and further
// properties of the shower.
//
// Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
////////////////////////////////////////////////////////////////////

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcorealg/Geometry/Exceptions.h"
#include "larreco/RecoAlg/EMShowerAlg.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
#include "TMath.h"
#include "TMathBase.h"
#include "TMultiGraph.h"
#include "TProfile.h"
#include "TFile.h"

shower::EMShowerAlg::EMShowerAlg(fhicl::ParameterSet const& pset) : fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
								    fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
								    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
								    fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")) {

  fMinTrackLength  = pset.get<double>("MinTrackLength");
  fdEdxTrackLength = pset.get<double>("dEdxTrackLength");
  fSpacePointSize  = pset.get<double>("SpacePointSize");
  fNfitpass        = pset.get<unsigned int>("Nfitpass");
  fNfithits        = pset.get<std::vector<unsigned int> >("Nfithits");
  fToler           = pset.get<std::vector<double> >("Toler");
  if (fNfitpass!=fNfithits.size() ||
      fNfitpass!=fToler.size()) {
    throw art::Exception(art::errors::Configuration)
      << "EMShowerAlg: fNfithits and fToler need to have size fNfitpass";
  }
  fDetector = pset.get<std::string>("Detector","dune35t");

  hTrueDirection = tfs->make<TH1I>("trueDir","",2,0,2);

  //this->MakePicture();
  // tmp
  fMakeGradientPlot = pset.get<bool>("MakeGradientPlot",false);
  fMakeRMSGradientPlot = pset.get<bool>("MakeRMSGradientPlot",false);
  fNumShowerSegments = pset.get<int>("NumShowerSegments",4);
  hNumHitsInSegment = tfs->make<TProfile>("NumHitsInSegment","",100,0,100);
  hNumSegments = tfs->make<TProfile>("NumSegments","",100,0,100);

}

void shower::EMShowerAlg::AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
						     art::FindManyP<recob::Hit> const& fmh,
						     art::FindManyP<recob::Track> const& fmt,
						     std::map<int,std::vector<int> >& clusterToTracks,
                                                     std::map<int,std::vector<int> >& trackToClusters) const {

  std::vector<int> clustersToIgnore = {-999};
  this->AssociateClustersAndTracks(clusters, fmh, fmt, clustersToIgnore, clusterToTracks, trackToClusters);

  return;

}

void shower::EMShowerAlg::AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
						     art::FindManyP<recob::Hit> const& fmh,
						     art::FindManyP<recob::Track> const& fmt,
						     std::vector<int> const& clustersToIgnore,
						     std::map<int,std::vector<int> >& clusterToTracks,
                                                     std::map<int,std::vector<int> >& trackToClusters) const {

  // Look through all the clusters
  for (std::vector<art::Ptr<recob::Cluster> >::const_iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt) {

    // Get the hits in this cluster
    std::vector<art::Ptr<recob::Hit> > clusterHits = fmh.at(clusterIt->key());

    // Look at all these hits and find the associated tracks
    for (std::vector<art::Ptr<recob::Hit> >::iterator clusterHitIt = clusterHits.begin(); clusterHitIt != clusterHits.end(); ++clusterHitIt) {

      // Get the tracks associated with this hit
      std::vector<art::Ptr<recob::Track> > clusterHitTracks = fmt.at(clusterHitIt->key());
      if (clusterHitTracks.size() > 1) { std::cout << "More than one track associated with this hit!" << std::endl; continue; }
      if (clusterHitTracks.size() < 1) continue;
      if (clusterHitTracks.at(0)->Length() < fMinTrackLength) {
	if (fDebug > 2)
	  std::cout << "Track " << clusterHitTracks.at(0)->ID() << " is too short! (" << clusterHitTracks.at(0)->Length() << ")" << std::endl;
	continue;
      }

      // Add this cluster to the track map
      int track = clusterHitTracks.at(0).key();
      //int trackID = clusterHitTracks.at(0)->ID();
      int cluster = (*clusterIt).key();
      if (std::find(clustersToIgnore.begin(), clustersToIgnore.end(), cluster) != clustersToIgnore.end())
	continue;
      if (std::find(trackToClusters[track].begin(), trackToClusters[track].end(), cluster) == trackToClusters[track].end())
	trackToClusters[track].push_back(cluster);
      if (std::find(clusterToTracks[cluster].begin(), clusterToTracks[cluster].end(), track) == clusterToTracks[cluster].end())
	clusterToTracks[cluster].push_back(track);

    }

  }

  return;

}

void shower::EMShowerAlg::CheckIsolatedHits(std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const {

  std::map<int,std::vector<int> > firstTPC;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    firstTPC[showerHitsIt->second.at(0)->WireID().TPC].push_back(showerHitsIt->first);

  // If all in the same TPC then that's great!
  if (firstTPC.size() == 1)
    return;

  // If they are in more than two TPCs, not much we can do
  else if (firstTPC.size() > 2)
    return;

  // If we get to this point, there should be something we can do!

  // Find the problem plane
  int problemPlane = -1;
  for (std::map<int,std::vector<int> >::iterator firstTPCIt = firstTPC.begin(); firstTPCIt != firstTPC.end(); ++firstTPCIt)
    if (firstTPCIt->second.size() == 1)
      problemPlane = firstTPCIt->second.at(0);

  // Require three hits
  if (showerHitsMap.at(problemPlane).size() < 3)
    return;

  // and get the other planes with at least three hits
  std::vector<int> otherPlanes;
  for (int plane = 0; plane < (int)fGeom->MaxPlanes(); ++plane)
    if (plane != problemPlane and
	showerHitsMap.count(plane) and
	showerHitsMap.at(plane).size() >= 3)
      otherPlanes.push_back(plane);

  if (otherPlanes.size() == 0)
    return;

  // Look at the hits after the first one
  unsigned int wrongTPC = showerHitsMap.at(problemPlane).at(0)->WireID().TPC;
  unsigned int nHits = 0;
  for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = showerHitsMap.at(problemPlane).begin();
       hitIt != showerHitsMap.at(problemPlane).end() and (*hitIt)->WireID().TPC == wrongTPC;
       ++hitIt)
    ++nHits;

  // If there are more than two hits in the 'wrong TPC', we can't be sure it is indeed wrong
  if (nHits > 2)
    return;

  // See if at least the next four times as many hits are in a different TPC
  std::map<int,int> otherTPCs;
  for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = std::next(showerHitsMap.at(problemPlane).begin(),nHits);
       hitIt != showerHitsMap.at(problemPlane).end() and std::distance(std::next(showerHitsMap.at(problemPlane).begin(),nHits),hitIt) < 4*nHits;
       ++hitIt)
    ++otherTPCs[(*hitIt)->WireID().TPC];

  if (otherTPCs.size() > 1)
    return;

  // If we get this far, we can move the problem hits from the front of the shower to the back
  std::map<int,int> tpcCount;
  for (std::vector<int>::iterator otherPlaneIt = otherPlanes.begin(); otherPlaneIt != otherPlanes.end(); ++otherPlaneIt)
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = std::next(showerHitsMap.at(*otherPlaneIt).begin());
	 hitIt != showerHitsMap.at(*otherPlaneIt).end() and hitIt != std::next(showerHitsMap.at(*otherPlaneIt).begin(),2);
	 ++hitIt)
      ++tpcCount[(*hitIt)->WireID().TPC];

  // Remove the first hit if it is in the wrong TPC
  if (tpcCount.size() == 1 and tpcCount.begin()->first == (int)(*std::next(showerHitsMap.at(problemPlane).begin(),nHits))->WireID().TPC) {
    std::vector<art::Ptr<recob::Hit> > naughty_hits;
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = showerHitsMap.at(problemPlane).begin(); hitIt != std::next(showerHitsMap.at(problemPlane).begin(),nHits); ++hitIt) {
      naughty_hits.push_back(*hitIt);
      showerHitsMap.at(problemPlane).erase(hitIt);
    }
    for (std::vector<art::Ptr<recob::Hit> >::iterator naughty_hitIt = naughty_hits.begin(); naughty_hitIt != naughty_hits.end(); ++naughty_hitIt)
      showerHitsMap.at(problemPlane).push_back(*naughty_hitIt);
  }

  return;

}

bool shower::EMShowerAlg::CheckShowerHits(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) const {

  bool consistencyCheck = true;

  if (showerHitsMap.size() < 2)
    consistencyCheck = true;

  else if (showerHitsMap.size() == 2) {

    // With two views, we can check:
    //  -- timing between views is consistent
    //  -- the 3D start point makes sense when projected back onto the individual planes

    std::vector<art::Ptr<recob::Hit> > startHits;
    std::vector<int> planes;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      startHits.push_back(showerHitsIt->second.front());
      planes.push_back(showerHitsIt->first);
    }

    TVector3 showerStartPos = Construct3DPoint(startHits.at(0), startHits.at(1));
    TVector2 proj1 = Project3DPointOntoPlane(showerStartPos, planes.at(0));
    TVector2 proj2 = Project3DPointOntoPlane(showerStartPos, planes.at(1));

    double timingDifference = TMath::Abs( startHits.at(0)->PeakTime() - startHits.at(1)->PeakTime() );
    double projectionDifference = ( (HitPosition(startHits.at(0)) - proj1).Mod() + (HitPosition(startHits.at(1)) - proj2).Mod() ) / (double)2;

    if (timingDifference > 40 or
	projectionDifference > 5 or
	showerStartPos.X() == -9999 or showerStartPos.Y() == -9999 or showerStartPos.Z() == -9999)
      consistencyCheck = false;

    if (fDebug > 1)
      std::cout << "Timing difference is " << timingDifference << " and projection distance is " << projectionDifference << " (start is (" << showerStartPos.X() << ", " << showerStartPos.Y() << ", " << showerStartPos.Z() << ")" << std::endl;

  }

  else if (showerHitsMap.size() == 3) {

    // With three views, we can check:
    //  -- the timing between views is consistent
    //  -- the 3D start point formed by two views and projected back into the third is close to the start point in that view

    std::map<int,art::Ptr<recob::Hit> > start2DMap;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitIt = showerHitsMap.begin(); showerHitIt != showerHitsMap.end(); ++showerHitIt)
      start2DMap[showerHitIt->first] = showerHitIt->second.front();

    std::map<int,double> projDiff;
    std::map<int,double> timingDiff;

    for (int plane = 0; plane < 3; ++plane) {

      std::vector<int> otherPlanes;
      for (int otherPlane = 0; otherPlane < 3; ++otherPlane)
	if (otherPlane != plane)
	  otherPlanes.push_back(otherPlane);

      TVector3 showerStartPos = Construct3DPoint(start2DMap.at(otherPlanes.at(0)), start2DMap.at(otherPlanes.at(1)));
      TVector2 showerStartProj = Project3DPointOntoPlane(showerStartPos, plane);

      if (fDebug > 2) {
	std::cout << "Plane... " << plane << std::endl;
	std::cout << "Start position in this plane is " << HitPosition(start2DMap.at(plane)).X() << ", " << HitPosition(start2DMap.at(plane)).Y() << ")" << std::endl;
	std::cout << "Shower start from other two planes is (" << showerStartPos.X() << ", " << showerStartPos.Y() << ", " << showerStartPos.Z() << ")" << std::endl;
	std::cout << "Projecting the other two planes gives position (" << showerStartProj.X() << ", " << showerStartProj.Y() << ")" << std::endl;
      }

      double projDiff = TMath::Abs((showerStartProj-HitPosition(start2DMap.at(plane))).Mod());
      double timeDiff = TMath::Max(TMath::Abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(0))->PeakTime()),
				   TMath::Abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(1))->PeakTime()));

      if (fDebug > 1)
	std::cout << "Plane " << plane << " has projDiff " << projDiff << " and timeDiff " << timeDiff << std::endl;

      if (projDiff > 5 or timeDiff > 40)
	consistencyCheck = false;

    }

  }

  if (fDebug > 1)
    std::cout << "Consistency check is " << consistencyCheck << std::endl;

  return consistencyCheck;

}

std::vector<int> shower::EMShowerAlg::CheckShowerPlanes(std::vector<std::vector<int> > const& initialShowers,
							std::vector<art::Ptr<recob::Cluster> > const& clusters,
                                                        art::FindManyP<recob::Hit> const& fmh) const {

  std::vector<int> clustersToIgnore;

  // // Look at each shower
  // for (std::vector<std::vector<int> >::const_iterator initialShowerIt = initialShowers.begin(); initialShowerIt != initialShowers.end(); ++initialShowerIt) {

  //   // Map the clusters and cluster hits to each view
  //   std::map<int,std::vector<art::Ptr<recob::Cluster> > > planeClusters;
  //   std::map<int,std::vector<art::Ptr<recob::Hit> > > planeClusterHits;
  //   for (std::vector<int>::const_iterator clusterIt = initialShowerIt->begin(); clusterIt != initialShowerIt->end(); ++clusterIt) {
  //     art::Ptr<recob::Cluster> cluster = clusters.at(*clusterIt);
  //     std::vector<art::Ptr<recob::Hit> > hits = fmh.at(cluster.key());
  //     planeClusters[cluster->Plane().Plane].push_back(cluster);
  //     for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
  // 	planeClusterHits[cluster.key()].push_back(*hitIt);
  //   }

  //   // Can't do much with fewer than three views
  //   if (planeClusters.size() < 3)
  //     continue;

  //   // Look at the average RMS of clusters in each view
  //   std::map<int,double> avRMS;
  //   for (std::map<int,std::vector<art::Ptr<recob::Cluster> > >::iterator planeClusterIt = planeClusters.begin(); planeClusterIt != planeClusters.end(); ++planeClusterIt) {
  //     std::cout << "Plane " << planeClusterIt->first << std::endl;
  //     for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeClusterIt->second.begin(); clusterIt != planeClusterIt->second.end(); ++clusterIt) {
  // 	double rms = ShowerHitRMS(planeClusterHits.at(clusterIt->key()));
  // 	std::cout << "Cluster " << clusterIt->key() << " has RMS " << rms << std::endl;
  //     }
  //   }
  // }


  // Look at each shower
  for (std::vector<std::vector<int> >::const_iterator initialShowerIt = initialShowers.begin(); initialShowerIt != initialShowers.end(); ++initialShowerIt) {

    if (std::distance(initialShowers.begin(),initialShowerIt) > 0)
      continue;

    // Map the clusters and cluster hits to each view
    std::map<int,std::vector<art::Ptr<recob::Cluster> > > planeClusters;
    std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHits;
    for (std::vector<int>::const_iterator clusterIt = initialShowerIt->begin(); clusterIt != initialShowerIt->end(); ++clusterIt) {
      art::Ptr<recob::Cluster> cluster = clusters.at(*clusterIt);
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(cluster.key());
      planeClusters[cluster->Plane().Plane].push_back(cluster);
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
  	planeHits[(*hitIt)->WireID().Plane].push_back(*hitIt);
    }

    TFile* outFile = new TFile("chargeDistributions.root","RECREATE");
    std::map<int,TH1D*> chargeDist;
    for (std::map<int,std::vector<art::Ptr<recob::Cluster> > >::iterator planeIt = planeClusters.begin(); planeIt != planeClusters.end(); ++planeIt) {
      for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeIt->second.begin(); clusterIt != planeIt->second.end(); ++clusterIt) {
	chargeDist[planeIt->first] = new TH1D(std::string("chargeDist_Plane"+std::to_string(planeIt->first)+"_Cluster"+std::to_string(clusterIt->key())).c_str(),"",150,0,1000);
	std::vector<art::Ptr<recob::Hit> > hits = fmh.at(clusterIt->key());
	for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
	  chargeDist[planeIt->first]->Fill((*hitIt)->Integral());
	outFile->cd();
	chargeDist[planeIt->first]->Write();
      }
    }
    outFile->Close();
    delete outFile;

    // Can't do much with fewer than three views
    if (planeClusters.size() < 3)
      continue;

    // Look at how many clusters each plane has, and the proportion of hits each one uses
    std::map<int,std::vector<double> > planeClusterSizes;
    for (std::map<int,std::vector<art::Ptr<recob::Cluster> > >::iterator planeClustersIt = planeClusters.begin(); planeClustersIt != planeClusters.end(); ++planeClustersIt) {
      for (std::vector<art::Ptr<recob::Cluster> >::iterator planeClusterIt = planeClustersIt->second.begin(); planeClusterIt != planeClustersIt->second.end(); ++planeClusterIt) {
  	std::vector<art::Ptr<recob::Hit> > hits = fmh.at(planeClusterIt->key());
        planeClusterSizes[planeClustersIt->first].push_back((double)hits.size()/(double)planeHits.at(planeClustersIt->first).size());
      }
    }

    // Find the average hit fraction across all clusters in the plane
    std::map<int,double> planeClustersAvSizes;
    for (std::map<int,std::vector<double> >::iterator planeClusterSizesIt = planeClusterSizes.begin(); planeClusterSizesIt != planeClusterSizes.end(); ++planeClusterSizesIt) {
      double average = 0;
      for (std::vector<double>::iterator planeClusterSizeIt = planeClusterSizesIt->second.begin(); planeClusterSizeIt != planeClusterSizesIt->second.end(); ++planeClusterSizeIt)
  	average += *planeClusterSizeIt;
      average /= planeClusterSizesIt->second.size();
      planeClustersAvSizes[planeClusterSizesIt->first] = average;
    }

    // Now decide if there is one plane which is ruining the reconstruction
    // If two planes have a low average cluster fraction and one high, this plane likely merges two particle deposits together
    int badPlane = -1;
    for (std::map<int,double>::iterator clusterAvSizeIt = planeClustersAvSizes.begin(); clusterAvSizeIt != planeClustersAvSizes.end(); ++clusterAvSizeIt) {

      // Get averages from other planes and add in quadrature
      std::vector<double> otherAverages;
      for (std::map<int,double>::iterator otherClustersAvSizeIt = planeClustersAvSizes.begin(); otherClustersAvSizeIt != planeClustersAvSizes.end(); ++otherClustersAvSizeIt)
    	if (clusterAvSizeIt->first != otherClustersAvSizeIt->first)
    	  otherAverages.push_back(otherClustersAvSizeIt->second);
      double quadOtherPlanes = 0;
      for (std::vector<double>::iterator otherAvIt = otherAverages.begin(); otherAvIt != otherAverages.end(); ++otherAvIt)
    	quadOtherPlanes += TMath::Power(*otherAvIt,2);
      quadOtherPlanes = TMath::Sqrt(quadOtherPlanes);

      // If this plane has an average higher than the quadratic sum of the others, it may be bad
      if (clusterAvSizeIt->second >= quadOtherPlanes)
    	badPlane = clusterAvSizeIt->first;

    }

    if (badPlane != -1) {
      if (fDebug > 1)
  	std::cout << "Bad plane is " << badPlane << std::endl;
      for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeClusters.at(badPlane).begin(); clusterIt != planeClusters.at(badPlane).end(); ++clusterIt)
  	clustersToIgnore.push_back(clusterIt->key());
    }

  }

  return clustersToIgnore;

}

TVector3 shower::EMShowerAlg::Construct3DPoint(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2) const {

  // x is average of the two x's
  double x = (fDetProp->ConvertTicksToX(hit1->PeakTime(), hit1->WireID().planeID()) + fDetProp->ConvertTicksToX(hit2->PeakTime(), hit2->WireID().planeID())) / (double)2;

  // y and z got from the wire interections
  geo::WireIDIntersection intersection;
  fGeom->WireIDsIntersect(hit1->WireID(), hit2->WireID(), intersection);

  return TVector3(x, intersection.y, intersection.z);

}

std::unique_ptr<recob::Track> shower::EMShowerAlg::ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& hits1,
								  std::vector<art::Ptr<recob::Hit> > const& hits2,
                                                                  std::map<int,TVector2> const& showerCentreMap) const {

  std::unique_ptr<recob::Track> track;

  std::vector<art::Ptr<recob::Hit> > track1, track2;

  // Check the TPCs
  if ((*hits1.begin())->WireID().TPC != (*hits2.begin())->WireID().TPC) {
    mf::LogWarning("EMShowerAlg") << "Warning: attempting to construct a track from two different TPCs.  Returning a null track.";
    return track;
  }

  // Check for tracks crossing TPC boundaries
  std::map<int,int> tpcMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits1.begin(); hitIt != hits1.end(); ++hitIt)
    ++tpcMap[(*hitIt)->WireID().TPC];
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits2.begin(); hitIt != hits2.end(); ++hitIt)
    ++tpcMap[(*hitIt)->WireID().TPC];
  if (tpcMap.size() > 1) {
    mf::LogWarning("EMShowerAlg") << "Warning: attempting to construct a track which crosses more than one TPC -- PMTrack can't handle this right now.  Returning a track made just from hits in the first TPC it traverses.";
    unsigned int firstTPC1 = hits1.at(0)->WireID().TPC, firstTPC2 = hits2.at(0)->WireID().TPC;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits1.begin(); hitIt != hits1.end(); ++hitIt)
      if ((*hitIt)->WireID().TPC == firstTPC1) track1.push_back(*hitIt);
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits2.begin(); hitIt != hits2.end(); ++hitIt)
      if ((*hitIt)->WireID().TPC == firstTPC2) track2.push_back(*hitIt);
  }
  else {
    track1 = hits1;
    track2 = hits2;
  }

  if (fDebug > 1) {
    std::cout << "About to make a track from these hits:" << std::endl;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit1 = track1.begin(); hit1 != track1.end(); ++hit1)
      std::cout << "Hit (" << HitCoordinates(*hit1).X() << ", " << HitCoordinates(*hit1).Y() << ") (real wire " << (*hit1)->WireID().Wire << ") in TPC " << (*hit1)->WireID().TPC << std::endl;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit2 = track2.begin(); hit2 != track2.end(); ++hit2)
      std::cout << "Hit (" << HitCoordinates(*hit2).X() << ", " << HitCoordinates(*hit2).Y() << ") (real wire " << (*hit2)->WireID().Wire << ") in TPC " << (*hit2)->WireID().TPC << std::endl;
  }

  TVector3 trackStart = Construct3DPoint(track1.at(0), track2.at(0));
  pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(track1, track2, trackStart);

  if (!pmatrack) {
    mf::LogInfo("EMShowerAlg") << "Skipping this event because not enough hits in two views";
    return track;
  }

  std::vector<TVector3> xyz, dircos;
  // std::vector<std::vector<double> > dEdx; // Right now, not finding the dE/dx for these tracks.  Can extend if needed.

  for (unsigned int i = 0; i < pmatrack->size(); i++) {

    xyz.push_back((*pmatrack)[i]->Point3D());

    if (i < pmatrack->size()-1) {
      size_t j = i+1;
      double mag = 0.0;
      TVector3 dc(0., 0., 0.);
      while ((mag == 0.0) and (j < pmatrack->size())) {
      	dc = (*pmatrack)[j]->Point3D();
      	dc -= (*pmatrack)[i]->Point3D();
      	mag = dc.Mag();
      	++j;
      }
      if (mag > 0.0) dc *= 1.0 / mag;
      else if (!dircos.empty()) dc = dircos.back();
      // TVector3 dc((*pmatrack)[i+1]->Point3D());
      // dc -= (*pmatrack)[i]->Point3D();
      // dc *= 1.0 / dc.Mag();
      dircos.push_back(dc);
    }
    else dircos.push_back(dircos.back());

  }

  // Orient the track correctly
  std::map<int,double> distanceToVertex, distanceToEnd;
  TVector3 vertex = *xyz.begin(), end = *xyz.rbegin();

  // Loop over all the planes and find the distance from the vertex and end projections to the centre in each plane
  for (std::map<int,TVector2>::const_iterator showerCentreIt = showerCentreMap.begin(); showerCentreIt != showerCentreMap.end(); ++showerCentreIt) {

    // Project the vertex and the end point onto this plane
    TVector2 vertexProj = Project3DPointOntoPlane(vertex, showerCentreIt->first);
    TVector2 endProj    = Project3DPointOntoPlane(end, showerCentreIt->first);

    // Find the distance of each to the centre of the cluster
    distanceToVertex[showerCentreIt->first] = (vertexProj - showerCentreIt->second).Mod();
    distanceToEnd[showerCentreIt->first] = (endProj - showerCentreIt->second).Mod();

  }

  // Find the average distance to the vertex and the end across the planes
  double avDistanceToVertex = 0, avDistanceToEnd = 0;
  for (std::map<int,double>::iterator distanceToVertexIt = distanceToVertex.begin(); distanceToVertexIt != distanceToVertex.end(); ++distanceToVertexIt)
    avDistanceToVertex += distanceToVertexIt->second;
  avDistanceToVertex /= distanceToVertex.size();

  for (std::map<int,double>::iterator distanceToEndIt = distanceToEnd.begin(); distanceToEndIt != distanceToEnd.end(); ++distanceToEndIt)
    avDistanceToEnd += distanceToEndIt->second;
  if (distanceToEnd.size() != 0)
    avDistanceToEnd /= distanceToEnd.size();

  if (fDebug > 2)
    std::cout << "Distance to vertex is " << avDistanceToVertex << " and distance to end is " << avDistanceToEnd << std::endl;

  // Change order if necessary
  if (avDistanceToEnd > avDistanceToVertex) {
    std::reverse(xyz.begin(), xyz.end());
    std::transform(dircos.begin(), dircos.end(), dircos.begin(), [](TVector3 const& vec){return -1*vec;});
  }

  if (xyz.size() != dircos.size())
    mf::LogError("EMShowerAlg") << "Problem converting pma::Track3D to recob::Track";

  track = std::make_unique<recob::Track>(recob::TrackTrajectory(recob::tracking::convertCollToPoint(xyz),
								recob::tracking::convertCollToVector(dircos),
								recob::Track::Flags_t(xyz.size()), false),
					 0, -1., 0, recob::tracking::SMatrixSym55(), recob::tracking::SMatrixSym55(), -1);

  return track;

}

std::unique_ptr<recob::Track> shower::EMShowerAlg::ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
                                                                  std::vector<art::Ptr<recob::Hit> > const& track2) const {

  std::map<int,TVector2> showerCentreMap;

  return this->ConstructTrack(track1, track2, showerCentreMap);

}

double shower::EMShowerAlg::FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, std::unique_ptr<recob::Track> const& track) const {

  double totalCharge = 0, totalDistance = 0, avHitTime = 0;
  unsigned int nHits = 0;

  if (!track)
    return -999;

  // size_t trajectory_point = 0;
  // double wirePitch   = fGeom->WirePitch(trackHits.at(0)->View(), 1, 0);
  // double angleToVert = fGeom->WireAngleToVertical(trackHits.at(0)->View(), 1, 0) - 0.5*::util::pi<>();
  // const TVector3& dir = track->DirectionAtPoint(trajectory_point);
  // //(sin(angleToVert),cos(angleToVert)) is the direction perpendicular to wire
  // double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() +
  // 			     std::cos(angleToVert)*dir.Z());

  // if(cosgamma < 1.e-5)
  //   throw cet::exception("Track") << "cosgamma is basically 0, that can't be right\n";
  // double pitch =  wirePitch/cosgamma;

  // Get the pitch
  double pitch = 0;
  try { pitch = lar::util::TrackPitchInView(*track, trackHits.at(0)->View()); }
  catch(...) { pitch = 0; }

  // Deal with large pitches
  if (pitch > fdEdxTrackLength) {
    double dEdx = fCalorimetryAlg.dEdx_AREA(*trackHits.begin(), pitch);
    return dEdx;
  }

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt) {
    if (totalDistance + pitch < fdEdxTrackLength) {
      totalDistance += pitch;
      totalCharge += (*trackHitIt)->Integral();
      avHitTime += (*trackHitIt)->PeakTime();
      ++nHits;
    }
  }

  avHitTime /= (double)nHits;

  double dQdx = totalCharge / totalDistance;
  double dEdx = fCalorimetryAlg.dEdx_AREA(dQdx, avHitTime, trackHits.at(0)->WireID().Plane);

  return dEdx;

}

void shower::EMShowerAlg::FindInitialTrack(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap,
					   std::unique_ptr<recob::Track>& initialTrack,
                                           std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits, int plane) const {

  /// Finding the initial track requires three stages:
  ///  -- put the hits in the correct order in each view
  ///  -- find the initial track-like hits in each view
  ///  -- use these to construct a track

  // // First, order the hits into the correct shower order in each plane
  // if (fDebug > 1)
  //   std::cout << " ------------------ Ordering shower hits -------------------- " << std::endl;
  // std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHitsMap = OrderShowerHits(hits, plane);
  // if (fDebug > 1)
  //   std::cout << " ------------------ End ordering shower hits -------------------- " << std::endl;

  // Now find the hits belonging to the track
  if (fDebug > 1)
    std::cout << " ------------------ Finding initial track hits -------------------- " << std::endl;
  initialTrackHits = FindShowerStart(showerHitsMap);
  if (fDebug > 1) {
    std::cout << "Here are the initial shower hits... " << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator initialTrackHitsIt = initialTrackHits.begin(); initialTrackHitsIt != initialTrackHits.end(); ++initialTrackHitsIt) {
      std::cout << "  Plane " << initialTrackHitsIt->first << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator initialTrackHitIt = initialTrackHitsIt->second.begin(); initialTrackHitIt != initialTrackHitsIt->second.end(); ++initialTrackHitIt)
	std::cout << "    Hit is (" << HitCoordinates(*initialTrackHitIt).X() << " (real hit " << (*initialTrackHitIt)->WireID() << "), " << HitCoordinates(*initialTrackHitIt).Y() << ")" << std::endl;
    }
  }
  if (fDebug > 1)
    std::cout << " ------------------ End finding initial track hits -------------------- " << std::endl;

  // Now we have the track hits -- can make a track!
  if (fDebug > 1)
    std::cout << " ------------------ Finding initial track -------------------- " << std::endl;
  initialTrack = MakeInitialTrack(initialTrackHits, showerHitsMap);
  if (initialTrack and fDebug > 1) {
    std::cout << "The track start is (" << initialTrack->Vertex().X() << ", " << initialTrack->Vertex().Y() << ", " << initialTrack->Vertex().Z() << ")" << std::endl;
    std::cout << "The track direction is (" << initialTrack->VertexDirection().X() << ", " << initialTrack->VertexDirection().Y() << ", " << initialTrack->VertexDirection().Z() << ")" << std::endl;
  }
  if (fDebug > 1)
    std::cout << " ------------------ End finding initial track -------------------- " << std::endl;

  // // Fill correct or incorrect direction histogram
  // std::map<int,int> trackHits;
  // for (art::PtrVector<recob::Hit>::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
  //   ++trackHits[FindTrackID(*hitIt)];
  // int trueTrack = -9999;
  // for (std::map<int,int>::iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt)
  //   if (trackHitIt->second/(double)hits.size() > 0.8)
  //     trueTrack = trackHitIt->first;
  // if (trueTrack != -9999) {
  //   const simb::MCParticle* trueParticle = backtracker->TrackIDToParticle(trueTrack);
  //   TVector3 trueStart = trueParticle->Position().Vect();
  //   if (initialTrack) {
  //     TVector3 recoStart = initialTrack->Vertex();
  //     if ((recoStart-trueStart).Mag() < 5)
  // 	hTrueDirection->Fill(1);
  //     else
  // 	hTrueDirection->Fill(0);
  //   }
  // }
  // else
  //   hTrueDirection->Fill(0);

  return;

}

std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::FindOrderOfHits(std::vector<art::Ptr<recob::Hit> > const& hits, bool perpendicular) const {

  // Find the charge-weighted centre (in [cm]) of this shower
  TVector2 centre = ShowerCentre(hits);

  // Find a rough shower 'direction'
  TVector2 direction = ShowerDirection(hits);

  if (perpendicular)
    direction = direction.Rotate(TMath::Pi()/2);

  // Find how far each hit (projected onto this axis) is from the centre
  TVector2 pos;
  std::map<double,art::Ptr<recob::Hit> > hitProjection;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
    pos = HitPosition(*hit) - centre;
    hitProjection[direction*pos] = *hit;
  }

  // Get a vector of hits in order of the shower
  std::vector<art::Ptr<recob::Hit> > showerHits;
  std::transform(hitProjection.begin(), hitProjection.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });

  // Make gradient plot
  if (fMakeGradientPlot) {
    std::map<int,TGraph*> graphs;
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = showerHits.begin(); hitIt != showerHits.end(); ++hitIt) {
      int tpc = (*hitIt)->WireID().TPC;
      if (graphs[tpc] == nullptr)
	graphs[tpc] = new TGraph();
      graphs[tpc]->SetPoint(graphs[tpc]->GetN(), HitPosition(*hitIt).X(), HitPosition(*hitIt).Y());
      //graphs[tpc]->SetPoint(graphs[tpc]->GetN(), HitCoordinates(*hitIt).X(), HitCoordinates(*hitIt).Y());
    }
    TMultiGraph* multigraph = new TMultiGraph();
    for (std::map<int,TGraph*>::iterator graphIt = graphs.begin(); graphIt != graphs.end(); ++graphIt) {
      graphIt->second->SetMarkerColor(graphIt->first);
      graphIt->second->SetMarkerStyle(8);
      graphIt->second->SetMarkerSize(2);
      multigraph->Add(graphIt->second);
    }
    TCanvas* canvas = new TCanvas();
    multigraph->Draw("AP");
    TLine line;
    line.SetLineColor(2);
    line.DrawLine(centre.X()-1000*direction.X(),centre.Y()-1000*direction.Y(),centre.X()+1000*direction.X(),centre.Y()+1000*direction.Y());
    canvas->SaveAs("Gradient.png");
    delete canvas; delete multigraph;
  }

  return showerHits;

}

std::vector<std::vector<int> > shower::EMShowerAlg::FindShowers(std::map<int,std::vector<int> > const& trackToClusters) const {

  // Showers are vectors of clusters
  std::vector<std::vector<int> > showers;

  // Loop over all tracks
  for (std::map<int,std::vector<int> >::const_iterator trackToClusterIt = trackToClusters.begin(); trackToClusterIt != trackToClusters.end(); ++ trackToClusterIt) {

    // Find which showers already made are associated with this track
    std::vector<int> matchingShowers;
    for (unsigned int shower = 0; shower < showers.size(); ++shower)
      for (std::vector<int>::const_iterator cluster = trackToClusterIt->second.begin(); cluster != trackToClusterIt->second.end(); ++cluster)
	if ( (std::find(showers.at(shower).begin(), showers.at(shower).end(), *cluster) != showers.at(shower).end()) and
	     (std::find(matchingShowers.begin(), matchingShowers.end(), shower)) == matchingShowers.end() )
	  matchingShowers.push_back(shower);

    // THINK THERE PROBABLY CAN BE MORE THAN ONE!
    // IN FACT, THIS WOULD BE A SUCCESS OF THE SHOWERING METHOD!
    // // Shouldn't be more than one
    // if (matchingShowers.size() > 1)
    //   mf::LogInfo("EMShowerAlg") << "Number of showers this track matches is " << matchingShowers.size() << std::endl;

    // New shower
    if (matchingShowers.size() < 1)
      showers.push_back(trackToClusterIt->second);

    // Add to existing shower
    else {
      for (std::vector<int>::const_iterator cluster = trackToClusterIt->second.begin(); cluster != trackToClusterIt->second.end(); ++cluster)
	if (std::find(showers.at(matchingShowers.at(0)).begin(), showers.at(matchingShowers.at(0)).end(), *cluster) == showers.at(matchingShowers.at(0)).end())
	  showers.at(matchingShowers.at(0)).push_back(*cluster);
    }
  }

  return showers;

}

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::FindShowerStart(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap) const {

  std::map<int,std::vector<art::Ptr<recob::Hit> > > initialHitsMap;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator orderedShowerIt = orderedShowerMap.begin(); orderedShowerIt != orderedShowerMap.end(); ++orderedShowerIt) {

    std::vector<art::Ptr<recob::Hit> > initialHits;
    const std::vector<art::Ptr<recob::Hit> > orderedShower = orderedShowerIt->second;

    // Find if the shower is travelling along ticks or wires
    bool wireDirection = true;
    std::vector<int> wires;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt)
      wires.push_back(std::round(HitCoordinates(*hitIt).X()));
    std::sort(wires.begin(), wires.end());
    if (TMath::Abs(*wires.begin()-std::round(HitCoordinates(*orderedShower.begin()).X())) > 3 and
    	TMath::Abs(*wires.rbegin()-std::round(HitCoordinates(*orderedShower.begin()).X())) > 3)
      wireDirection = false;

    // Deal with showers travelling along wires
    if (wireDirection) {
      bool increasing = HitCoordinates(*orderedShower.rbegin()).X() > HitCoordinates(*orderedShower.begin()).X();
      std::map<int,std::vector<art::Ptr<recob::Hit> > > wireHitMap;
      int multipleWires = 0;
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt)
	wireHitMap[std::round(HitCoordinates(*hitIt).X())].push_back(*hitIt);
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt) {
	int wire = std::round(HitCoordinates(*hitIt).X());
	if (wireHitMap[wire].size() > 1) {
	  ++multipleWires;
	  if (multipleWires > 5) break;
	  continue;
	}
	else if ( (increasing and wireHitMap[wire+1].size() > 1 and (wireHitMap[wire+2].size() > 1 or wireHitMap[wire+3].size() > 1)) or
		  (!increasing and wireHitMap[wire-1].size() > 1 and (wireHitMap[wire-2].size() > 1 or wireHitMap[wire-3].size() > 1)) ) {
	  if ( (increasing and (wireHitMap[wire+4].size() < 2 and wireHitMap[wire+5].size() < 2 and wireHitMap[wire+6].size() < 2 and wireHitMap[wire+9].size() > 1)) or
	       (!increasing and (wireHitMap[wire-4].size() < 2 and wireHitMap[wire-5].size() < 2 and wireHitMap[wire-6].size() < 2) and wireHitMap[wire-9].size() > 1) )
	    initialHits.push_back(*hitIt);
	  else
	    break;
	}
	else
	  initialHits.push_back(*hitIt);
      }
      if (!initialHits.size()) initialHits.push_back(*orderedShower.begin());
    }

    // Deal with showers travelling along ticks
    else {
      bool increasing = HitCoordinates(*orderedShower.rbegin()).Y() > HitCoordinates(*orderedShower.begin()).Y();
      std::map<int,std::vector<art::Ptr<recob::Hit> > > tickHitMap;
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt)
	tickHitMap[std::round(HitCoordinates(*hitIt).Y())].push_back(*hitIt);
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt) {
	int tick = std::round(HitCoordinates(*hitIt).Y());
	if ( (increasing and (tickHitMap[tick+1].size() or tickHitMap[tick+2].size() or tickHitMap[tick+3].size() or tickHitMap[tick+4].size() or tickHitMap[tick+5].size())) or
	     (!increasing and (tickHitMap[tick-1].size() or tickHitMap[tick-2].size() or tickHitMap[tick-3].size() or tickHitMap[tick-4].size() or tickHitMap[tick-5].size())) )
	  break;
	else
	  initialHits.push_back(*hitIt);
      }
      if (!initialHits.size()) initialHits.push_back(*orderedShower.begin());
    }

    // Need at least two hits to make a track
    if (initialHits.size() == 1 and orderedShower.size() > 2)
      initialHits.push_back(orderedShower.at(1));

    // Quality check -- make sure there isn't a single hit in a different TPC (artefact of tracking failure)
    std::vector<art::Ptr<recob::Hit> > newInitialHits;
    std::map<int,int> tpcHitMap;
    std::vector<int> singleHitTPCs;
    for (std::vector<art::Ptr<recob::Hit> >::iterator initialHitIt = initialHits.begin(); initialHitIt != initialHits.end(); ++initialHitIt)
      ++tpcHitMap[(*initialHitIt)->WireID().TPC];
    for (std::map<int,int>::iterator tpcIt = tpcHitMap.begin(); tpcIt != tpcHitMap.end(); ++tpcIt)
      if (tpcIt->second == 1) singleHitTPCs.push_back(tpcIt->first);
    if (singleHitTPCs.size()) {
      if (fDebug > 2)
	for (std::vector<int>::iterator tpcIt = singleHitTPCs.begin(); tpcIt != singleHitTPCs.end(); ++tpcIt)
	  std::cout << "Removed hits in TPC " << *tpcIt << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator initialHitIt = initialHits.begin(); initialHitIt != initialHits.end(); ++initialHitIt)
	if (std::find(singleHitTPCs.begin(), singleHitTPCs.end(), (*initialHitIt)->WireID().TPC) == singleHitTPCs.end())
	  newInitialHits.push_back(*initialHitIt);
      if (!newInitialHits.size()) newInitialHits.push_back(*orderedShower.begin());
    }
    else
      newInitialHits = initialHits;

    initialHitsMap[orderedShowerIt->first] = newInitialHits;

  }

  return initialHitsMap;

}

std::map<int,std::map<int,bool> > shower::EMShowerAlg::GetPlanePermutations(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const {

  // The map to return
  std::map<int,std::map<int,bool> > permutations;

  // Get the properties of the shower hits across the planes which will be used to
  // determine the likelihood of a particular reorientation permutation
  //   -- relative width in the wire direction (if showers travel along the wire
  //      direction in a particular plane)
  //   -- the RMS gradients (how likely it is the RMS of the hit positions from
  //      central axis increases along a particular orientation)

  // Find the RMS, RMS gradient and wire widths
  std::map<int,double> planeRMSGradients, planeRMS;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    planeRMS[showerHitsIt->first] = ShowerHitRMS(showerHitsIt->second);
    planeRMSGradients[showerHitsIt->first] = ShowerHitRMSGradient(showerHitsIt->second);
  }

  // Order these backwards so they can be used to discriminate between planes
  std::map<double,int> gradientMap;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    gradientMap[TMath::Abs(planeRMSGradients.at(showerHitsIt->first))] = showerHitsIt->first;

  std::map<double,int> wireWidthMap = RelativeWireWidth(showerHitsMap);

  if (fDebug > 1)
    for (std::map<double,int>::const_iterator wireWidthIt = wireWidthMap.begin(); wireWidthIt != wireWidthMap.end(); ++wireWidthIt)
      std::cout << "Plane " << wireWidthIt->second << " has relative width in wire of " << wireWidthIt->first << "; and an RMS gradient of " << planeRMSGradients[wireWidthIt->second] << std::endl;

  // Find the correct permutations
  int perm = 0;
  std::vector<std::map<int,bool> > usedPermutations;

  // Most likely is to not change anything
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    permutations[perm][showerHitsIt->first] = 0;
  ++perm;

  // Use properties of the shower to determine the middle cases
  for (std::map<double,int>::iterator wireWidthIt = wireWidthMap.begin(); wireWidthIt != wireWidthMap.end(); ++wireWidthIt) {
    std::map<int,bool> permutation;
    permutation[wireWidthIt->second] = 1;
    for (std::map<double,int>::iterator wireWidth2It = wireWidthMap.begin(); wireWidth2It != wireWidthMap.end(); ++wireWidth2It)
      if (wireWidthIt->second != wireWidth2It->second)
	permutation[wireWidth2It->second] = 0;
    if (std::find(usedPermutations.begin(), usedPermutations.end(), permutation) == usedPermutations.end()) {
      permutations[perm] = permutation;
      usedPermutations.push_back(permutation);
      ++perm;
    }
  }
  for (std::map<double,int>::reverse_iterator wireWidthIt = wireWidthMap.rbegin(); wireWidthIt != wireWidthMap.rend(); ++wireWidthIt) {
    std::map<int,bool> permutation;
    permutation[wireWidthIt->second] = 0;
    for (std::map<double,int>::iterator wireWidth2It = wireWidthMap.begin(); wireWidth2It != wireWidthMap.end(); ++wireWidth2It)
      if (wireWidthIt->second != wireWidth2It->second)
	permutation[wireWidth2It->second] = 1;
    if (std::find(usedPermutations.begin(), usedPermutations.end(), permutation) == usedPermutations.end()) {
      permutations[perm] = permutation;
      usedPermutations.push_back(permutation);
      ++perm;
    }
  }

  // Least likely is to change everything
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    permutations[perm][showerHitsIt->first] = 1;
  ++perm;

  if (fDebug > 1) {
    std::cout << "Here are the permutations!" << std::endl;
    for (std::map<int,std::map<int,bool> >::iterator permIt = permutations.begin(); permIt != permutations.end(); ++permIt) {
      std::cout << "Permutation " << permIt->first << std::endl;
      for (std::map<int,bool>::iterator planePermIt = permIt->second.begin(); planePermIt != permIt->second.end(); ++planePermIt)
	std::cout << "  Plane " << planePermIt->first << " has value " << planePermIt->second << std::endl;
    }
  }

  return permutations;

}

std::unique_ptr<recob::Track> shower::EMShowerAlg::MakeInitialTrack(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap,
                                                                    std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) const {

  // Can't do much with just one view
  if (initialHitsMap.size() < 2) {
    mf::LogInfo("EMShowerAlg") << "Only one useful view for this shower; nothing can be done" << std::endl;
    return std::unique_ptr<recob::Track>();
  }

  std::vector<std::pair<int,int> > initialHitsSize;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator initialHitIt = initialHitsMap.begin(); initialHitIt != initialHitsMap.end(); ++initialHitIt)
    initialHitsSize.push_back(std::make_pair(initialHitIt->first, initialHitIt->second.size()));

  // Sort the planes by number of hits
  std::sort(initialHitsSize.begin(), initialHitsSize.end(), [](std::pair<int,int> const& size1, std::pair<int,int> const& size2) { return size1.second > size2.second; });

  // Pick the two planes to use to make the track
  //   -- if more than two planes, can choose the two views
  //      more accurately
  //   -- if not, just use the two views available

  std::vector<int> planes;

  if (initialHitsSize.size() > 2 and !CheckShowerHits(showerHitsMap)) {
    int planeToIgnore = WorstPlane(showerHitsMap);
    if (fDebug > 1)
      std::cout << "Igoring plane " << planeToIgnore << " in creation of initial track" << std::endl;
    for (std::vector<std::pair<int,int> >::iterator hitsSizeIt = initialHitsSize.begin(); hitsSizeIt != initialHitsSize.end() and planes.size() < 2; ++hitsSizeIt) {
      if (hitsSizeIt->first == planeToIgnore)
	continue;
      planes.push_back(hitsSizeIt->first);
    }
  }
  else
    planes = {initialHitsSize.at(0).first, initialHitsSize.at(1).first};

  return ConstructTrack(initialHitsMap.at(planes.at(0)), initialHitsMap.at(planes.at(1)));

}

recob::Shower shower::EMShowerAlg::MakeShower(art::PtrVector<recob::Hit> const& hits,
					      std::unique_ptr<recob::Track> const& initialTrack,
                                              std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap) const {

  //return recob::Shower();

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  int bestPlane = -1;
  unsigned int highestNumberOfHits = 0;
  std::vector<double> totalEnergy, totalEnergyError, dEdx, dEdxError;

  // Look at each of the planes
  for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {

    // If there's hits on this plane...
    if (planeHitsMap.count(plane)) {
      dEdx.push_back(FinddEdx(initialHitsMap.at(plane), initialTrack));
      totalEnergy.push_back(fShowerEnergyAlg.ShowerEnergy(planeHitsMap.at(plane), plane));
      if (planeHitsMap.at(plane).size() > highestNumberOfHits and initialHitsMap.count(plane)) {
	bestPlane = plane;
	highestNumberOfHits = planeHitsMap.at(plane).size();
      }
    }

    // If not...
    else {
      dEdx.push_back(0);
      totalEnergy.push_back(0);
    }

  }

  TVector3 direction, directionError, showerStart, showerStartError;
  if (initialTrack) {
    direction = initialTrack->VertexDirection<TVector3>();
    showerStart = initialTrack->Vertex<TVector3>();
  }

  if (fDebug > 0) {
    std::cout << "Best plane is " << bestPlane << std::endl;
    std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
    std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
    std::cout << "The shower start is (" << showerStart.X() << ", " << showerStart.Y() << ", " << showerStart.Z() << ")" << std::endl;
    std::cout << "The shower direction is (" << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << std::endl;
  }

  return recob::Shower(direction, directionError, showerStart, showerStartError, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane);

}

recob::Shower shower::EMShowerAlg::MakeShower(art::PtrVector<recob::Hit> const& hits,
					      art::Ptr<recob::Vertex> const& vertex,
                                              int & iok) const {

  iok = 1;

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->WireID().Plane].push_back(*hit);

  std::vector<std::vector<art::Ptr<recob::Hit> > > initialTrackHits(3);

  int pl0 = -1;
  int pl1 = -1;
  unsigned maxhits0 = 0;
  unsigned maxhits1 = 0;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {

    std::vector<art::Ptr<recob::Hit> > showerHits;
    OrderShowerHits(planeHits->second, showerHits, vertex);
    //if (!isCleanShower(showerHits)) continue;
    FindInitialTrackHits(showerHits, vertex, initialTrackHits[planeHits->first]);
    if ((planeHits->second).size()>maxhits0){
      if (pl0!=-1){
	maxhits1 = maxhits0;
	pl1 = pl0;
      }
      pl0 = planeHits->first;
      maxhits0 = (planeHits->second).size();
    }
    else if ((planeHits->second).size()>maxhits1){
      pl1 = planeHits->first;
      maxhits1 = (planeHits->second).size();
    }

  }
  //std::cout<<pl0<<" "<<pl1<<std::endl;
//  if (pl0!=-1&&pl1!=-1) {
//    pl0 = 1;
//    pl1 = 2;
//  }
  if (pl0!=-1&&pl1!=-1
      &&initialTrackHits[pl0].size()>=2
      &&initialTrackHits[pl1].size()>=2
      &&initialTrackHits[pl0][0]->WireID().TPC==
      initialTrackHits[pl1][0]->WireID().TPC){
    double xyz[3];
    vertex->XYZ(xyz);
    TVector3 vtx(xyz);
//    std::vector<art::Ptr<recob::Hit>> alltrackhits;
//    for (size_t i = 0; i<3; ++i){
//      for (auto const&hit : initialTrackHits[i]){
//	alltrackhits.push_back(hit);
//      }
//    }
    //std::cout<<"vertex "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;
    //for (auto const&hit : initialTrackHits[pl0]) std::cout<<*hit<<std::endl;
    //for (auto const&hit : initialTrackHits[pl1]) std::cout<<*hit<<std::endl;
    pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(initialTrackHits[pl0], initialTrackHits[pl1]);
    //std::cout<<pmatrack->size()<<std::endl;
    //pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(alltrackhits);
    std::vector<TVector3> spts;

//  points are shifted now by tracking code, commented out shifts here
//    double xshift = pmatrack->GetXShift();
//    bool has_shift = (xshift != 0.0);
    for (size_t i = 0; i<pmatrack->size(); ++i){
      if ((*pmatrack)[i]->IsEnabled()){
	TVector3 p3d = (*pmatrack)[i]->Point3D();
//	if (has_shift) p3d.SetX(p3d.X() + xshift);
	//std::cout<<p3d.X()<<" "<<p3d.Y()<<" "<<p3d.Z()<<std::endl;
	spts.push_back(p3d);
      }
    }
    if (spts.size()>=2){ //at least two space points
      TVector3 shwxyz, shwxyzerr;
      TVector3 shwdir, shwdirerr;
      std::vector<double> totalEnergy, totalEnergyError, dEdx, dEdxError;
      int bestPlane = pl0;
      double minpitch = 1000;
      std::vector<TVector3> dirs;
      if ((spts[0]-vtx).Mag()<(spts.back()-vtx).Mag()){
	shwxyz = spts[0];
	size_t i = 5;
	if (spts.size()-1<5) i = spts.size()-1;
	shwdir = spts[i] - spts[0];
	shwdir = shwdir.Unit();
      }
      else{
	shwxyz = spts.back();
	size_t i = 0;
	if (spts.size()>6) i = spts.size() - 6;
	shwdir = spts[i] - spts[spts.size()-1];
	shwdir = shwdir.Unit();
      }
      //std::cout<<shwxyz.X()<<" "<<shwxyz.Y()<<" "<<shwxyz.Z()<<std::endl;
      //std::cout<<shwdir.X()<<" "<<shwdir.Y()<<" "<<shwdir.Z()<<std::endl;
      for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {
	if (planeHitsMap.find(plane)!=planeHitsMap.end()){
	  totalEnergy.push_back(fShowerEnergyAlg.ShowerEnergy(planeHitsMap[plane], plane));
	}
	else{
	  totalEnergy.push_back(0);
	}
	if (initialTrackHits[plane].size()){
	  double fdEdx = 0;
	  double totQ = 0;
	  double avgT = 0;
	  double pitch = 0;
	  double wirepitch = fGeom->WirePitch(initialTrackHits[plane][0]->WireID().planeID());
	  double angleToVert = fGeom->WireAngleToVertical(fGeom->Plane(plane).View(),initialTrackHits[plane][0]->WireID().planeID()) - 0.5*TMath::Pi();
	  double cosgamma = std::abs(sin(angleToVert)*shwdir.Y()+
				     cos(angleToVert)*shwdir.Z());
	  if (cosgamma>0) pitch = wirepitch/cosgamma;
	  if (pitch){
	    if (pitch<minpitch){
	      minpitch = pitch;
	      bestPlane = plane;
	    }
	    int nhits = 0;
	    //std::cout<<"pitch = "<<pitch<<std::endl;
            std::vector<float> vQ;
	    for (auto const& hit: initialTrackHits[plane]){
	      //std::cout<<hit->WireID()<<" "<<hit->PeakTime()<<" "<<std::abs((hit->WireID().Wire-initialTrackHits[plane][0]->WireID().Wire)*pitch)<<" "<<fdEdxTrackLength<<std::endl;
	      int w1 = hit->WireID().Wire;
	      int w0 = initialTrackHits[plane][0]->WireID().Wire;
	      if (std::abs((w1-w0)*pitch)<fdEdxTrackLength){
                vQ.push_back(hit->Integral());
		totQ += hit->Integral();
		avgT+= hit->PeakTime();
		++nhits;
		//std::cout<<hit->WireID()<<" "<<hit->PeakTime()<<" "<<hit->Integral()<<" "<<totQ<<" "<<avgT<<std::endl;
	      }
	    }
	    if (totQ) {
	      //double dQdx = totQ/(nhits*pitch);
              //std::sort(vQ.begin(), vQ.end());
              //double dQdx = vQ[vQ.size()/2]/pitch;
              double dQdx = TMath::Median(vQ.size(), &vQ[0])/pitch;
	      fdEdx = fCalorimetryAlg.dEdx_AREA(dQdx, avgT/nhits, initialTrackHits[plane][0]->WireID().Plane);
	    }
	  }
	  dEdx.push_back(fdEdx);
	}
	else{
	  dEdx.push_back(0);
	}
      }
      iok = 0;
      if (fDebug > 1) {
	std::cout << "Best plane is " << bestPlane << std::endl;
	std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
	std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
	std::cout << "The shower start is (" << shwxyz.X() << ", " << shwxyz.Y() << ", " << shwxyz.Z() << ")" << std::endl;
	shwxyz.Print();
      }

      return recob::Shower(shwdir, shwdirerr, shwxyz, shwxyzerr, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane);
    }
  }
  return recob::Shower();
}

std::vector<recob::SpacePoint> shower::EMShowerAlg::MakeSpacePoints(std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHits,
                                                                    std::vector<std::vector<art::Ptr<recob::Hit> > >& hitAssns) const {

  // Space points to return
  std::vector<recob::SpacePoint> spacePoints;

  // // Order by charge
  // for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = showerHits.begin(); planeIt != showerHits.end(); ++planeIt)
  //   std::sort(planeIt->second.begin(), planeIt->second.end(), [](const art::Ptr<recob::Hit>& h1, const art::Ptr<recob::Hit>& h2) { return h1->Integral() > h2->Integral(); });

  // Make space points
  // Use the following procedure:
  //  -- Consider hits plane by plane
  //  -- For each hit on the first plane, consider the 3D point made by combining with each hit from the second plane
  //  -- Project this 3D point back into the two planes
  //  -- Determine how close to a the original hits this point lies
  //  -- If close enough, make a 3D space point from this point
  // //  -- Discard these used hits in future iterations, along with hits in the
  // //       third plane (if exists) close to the projection of the point into this plane

  // Container to hold used hits
  std::vector<int> usedHits;

  // Look through plane by plane
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitIt = showerHits.begin(); showerHitIt != showerHits.end(); ++showerHitIt) {

    // Find the other planes with hits
    std::vector<int> otherPlanes;
    for (unsigned int otherPlane = 0; otherPlane < fGeom->MaxPlanes(); ++otherPlane)
      if ((int)otherPlane != showerHitIt->first and showerHits.count(otherPlane))
	otherPlanes.push_back(otherPlane);

    // Can't make space points if we only have one view
    if (otherPlanes.size() == 0)
      return spacePoints;

    // Look at all hits on this plane
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator planeHitIt = showerHitIt->second.begin(); planeHitIt != showerHitIt->second.end(); ++planeHitIt) {

      if (std::find(usedHits.begin(), usedHits.end(), planeHitIt->key()) != usedHits.end())
	continue;

      // Make a 3D point with every hit on the second plane
      const std::vector<art::Ptr<recob::Hit> > otherPlaneHits = showerHits.at(otherPlanes.at(0));
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator otherPlaneHitIt = otherPlaneHits.begin();
	   otherPlaneHitIt != otherPlaneHits.end() and std::find(usedHits.begin(), usedHits.end(), planeHitIt->key()) == usedHits.end();
	   ++otherPlaneHitIt) {

	if ((*otherPlaneHitIt)->WireID().TPC != (*planeHitIt)->WireID().TPC or
	    std::find(usedHits.begin(), usedHits.end(), otherPlaneHitIt->key()) != usedHits.end())
	  continue;

	TVector3 point = Construct3DPoint(*planeHitIt, *otherPlaneHitIt);
	std::vector<art::Ptr<recob::Hit> > pointHits;
	bool truePoint = false;

	if (otherPlanes.size() > 1) {

	  TVector2 projThirdPlane = Project3DPointOntoPlane(point, otherPlanes.at(1));
	  const std::vector<art::Ptr<recob::Hit> > otherOtherPlaneHits = showerHits.at(otherPlanes.at(1));

	  for (std::vector<art::Ptr<recob::Hit> >::const_iterator otherOtherPlaneHitIt = otherOtherPlaneHits.begin();
	       otherOtherPlaneHitIt != otherOtherPlaneHits.end() and !truePoint;
	       ++otherOtherPlaneHitIt) {

	    if ((*otherOtherPlaneHitIt)->WireID().TPC == (*planeHitIt)->WireID().TPC and
		(projThirdPlane-HitPosition(*otherOtherPlaneHitIt)).Mod() < fSpacePointSize) {

	      truePoint = true;

	      // Remove hits used to make the point
	      usedHits.push_back(planeHitIt->key());
	      usedHits.push_back(otherPlaneHitIt->key());
	      usedHits.push_back(otherOtherPlaneHitIt->key());

	      pointHits.push_back(*planeHitIt);
	      pointHits.push_back(*otherPlaneHitIt);
	      pointHits.push_back(*otherOtherPlaneHitIt);

	    }
	  }
	}

	else if ((Project3DPointOntoPlane(point, (*planeHitIt)->WireID().Plane) - HitPosition(*planeHitIt)).Mod() < fSpacePointSize and
		 (Project3DPointOntoPlane(point, (*otherPlaneHitIt)->WireID().Plane) - HitPosition(*otherPlaneHitIt)).Mod() < fSpacePointSize) {

	  truePoint = true;

	  usedHits.push_back(planeHitIt->key());
	  usedHits.push_back(otherPlaneHitIt->key());

	  pointHits.push_back(*planeHitIt);
	  pointHits.push_back(*otherPlaneHitIt);

	}

	// Make space point
	if (truePoint) {
	  double xyz[3] = {point.X(), point.Y(), point.Z()};
	  double xyzerr[6] = {fSpacePointSize, fSpacePointSize, fSpacePointSize, fSpacePointSize, fSpacePointSize, fSpacePointSize};
	  double chisq = 0.;
	  spacePoints.emplace_back(xyz, xyzerr, chisq);
	  hitAssns.push_back(pointHits);
	}

      } // end loop over second plane hits

    } // end loop over first plane hits

  } // end loop over planes

  if (fDebug > 0) {
    std::cout << "-------------------- Space points -------------------" << std::endl;
    std::cout << "There are " << spacePoints.size() << " space points:" << std::endl;
    if (fDebug > 1)
      for (std::vector<recob::SpacePoint>::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {
	const double* xyz = spacePointIt->XYZ();
	std::cout << "  Space point (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")" << std::endl;
      }
  }

  return spacePoints;

}

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::OrderShowerHits(art::PtrVector<recob::Hit> const& shower, int plane) const {

  /// Ordering the shower hits requires three stages:
  ///  -- putting all the hits in a given plane in some kind of order
  ///  -- use the properties of the hits in all three planes to check this order
  ///  -- orient the hits correctly using properties of the shower

  // ------------- Put hits in order ------------

  // // Find the true start for reconstruction validation purposes
  // std::vector<art::Ptr<recob::Hit> > showerHitsVector;
  // std::map<int,int> trueParticleHits;
  // for (art::PtrVector<recob::Hit>::const_iterator showerHitIt = shower.begin(); showerHitIt != shower.end(); ++showerHitIt) {
  //   showerHitsVector.push_back(*showerHitIt);
  //   ++trueParticleHits[FindParticleID(*showerHitIt)];
  // }
  // int trueParticleID = FindTrueParticle(showerHitsVector);
  // const simb::MCParticle* trueParticle = bt_serv->TrackIDToParticle(trueParticleID);
  // TVector3 trueStart3D = trueParticle->Position().Vect();
  // double purity = trueParticleHits[trueParticleID]/(double)shower.size();

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
    showerHitsMap[(*hit)->WireID().Plane].push_back(*hit);

  // if (purity < 0.9)
  //   return showerHitsMap;

  // Order the hits, get the RMS and the RMS gradient for the hits in this plane
  std::map<int,double> planeRMSGradients, planeRMS;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    if (plane != showerHitsIt->first and plane != -1)
      continue;
    std::vector<art::Ptr<recob::Hit> > orderedHits = FindOrderOfHits(showerHitsIt->second);
    planeRMS[showerHitsIt->first] = ShowerHitRMS(orderedHits);
    //TVector2 trueStart2D = Project3DPointOntoPlane(trueStart3D, showerHitsIt->first);
    planeRMSGradients[showerHitsIt->first] = ShowerHitRMSGradient(orderedHits);
    showerHitsMap[showerHitsIt->first] = orderedHits;
  }

  if (fDebug > 1)
    for (std::map<int,double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end(); ++planeRMSIt)
      std::cout << "Plane " << planeRMSIt->first << " has RMS " << planeRMSIt->second << " and RMS gradient " << planeRMSGradients.at(planeRMSIt->first) << std::endl;

  if (fDebug > 2) {
    std::cout << std::endl << "Hits in order; after ordering, before reversing..." << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      std::cout << "  Plane " << showerHitsIt->first << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = showerHitsIt->second.begin(); hitIt != showerHitsIt->second.end(); ++hitIt)
	std::cout << "    Hit at (" << HitCoordinates(*hitIt).X() << ", " << HitCoordinates(*hitIt).Y() << ") -- real wire " << (*hitIt)->WireID() << ", hit position (" << HitPosition(*hitIt).X() << ", " << HitPosition(*hitIt).Y() << ")" << std::endl;
    }
  }

  // ------------- Check between the views to ensure consistency of ordering -------------

  // Check between the views to make sure there isn't a poorly formed shower in just one view
  // First, determine the average RMS and RMS gradient across the other planes
  std::map<int,double> planeOtherRMS, planeOtherRMSGradients;
  for (std::map<int,double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end(); ++planeRMSIt) {
    planeOtherRMS[planeRMSIt->first] = 0;
    planeOtherRMSGradients[planeRMSIt->first] = 0;
    int nOtherPlanes = 0;
    for (int plane = 0; plane < (int)fGeom->MaxPlanes(); ++plane) {
      if (plane != planeRMSIt->first and planeRMS.count(plane)) {
	planeOtherRMS[planeRMSIt->first] += planeRMS.at(plane);
	planeOtherRMSGradients[planeRMSIt->first] += planeRMSGradients.at(plane);
	++nOtherPlanes;
      }
    }
    planeOtherRMS[planeRMSIt->first] /= (double)nOtherPlanes;
    planeOtherRMSGradients[planeRMSIt->first] /= (double)nOtherPlanes;
  }

  // Look to see if one plane has a particularly high RMS (compared to the others) whilst having a similar gradient
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    if (planeRMS.at(showerHitsIt->first) > planeOtherRMS.at(showerHitsIt->first) * 2) {
      //and TMath::Abs(planeRMSGradients.at(showerHitsIt->first) / planeOtherRMSGradients.at(showerHitsIt->first)) < 0.1) {
      if (fDebug > 1)
	std::cout << "Plane " << showerHitsIt->first << " was perpendicular... recalculating" << std::endl;
      std::vector<art::Ptr<recob::Hit> > orderedHits = this->FindOrderOfHits(showerHitsIt->second, true);
      showerHitsMap[showerHitsIt->first] = orderedHits;
      planeRMSGradients[showerHitsIt->first] = this->ShowerHitRMSGradient(orderedHits);
    }
  }

  // ------------- Orient the shower correctly ---------------

  if (fDebug > 1) {
    std::cout << "Before reversing, here are the start and end points..." << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
      std::cout << "  Plane " << showerHitsIt->first << " has start (" << HitCoordinates(showerHitsIt->second.front()).X() << ", " << HitCoordinates(showerHitsIt->second.front()).Y() << ") (real wire " << showerHitsIt->second.front()->WireID() << ") and end (" << HitCoordinates(showerHitsIt->second.back()).X() << ", " << HitCoordinates(showerHitsIt->second.back()).Y() << ") (real wire " << showerHitsIt->second.back()->WireID() << ")" << std::endl;
  }

  // Use the RMS gradient information to get an initial ordering
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    double gradient = planeRMSGradients.at(showerHitsIt->first);
    if (gradient < 0)
      std::reverse(showerHitsIt->second.begin(), showerHitsIt->second.end());
  }

  if (fDebug > 2) {
    std::cout << std::endl << "Hits in order; after reversing, before checking isolated hits..." << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      std::cout << "  Plane " << showerHitsIt->first << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = showerHitsIt->second.begin(); hitIt != showerHitsIt->second.end(); ++hitIt)
	std::cout << "    Hit at (" << HitCoordinates(*hitIt).X() << ", " << HitCoordinates(*hitIt).Y() << ") -- real wire " << (*hitIt)->WireID() << ", hit position (" << HitPosition(*hitIt).X() << ", " << HitPosition(*hitIt).Y() << ")" << std::endl;
    }
  }

  CheckIsolatedHits(showerHitsMap);

  if (fDebug > 2) {
    std::cout << std::endl << "Hits in order; after checking isolated hits..." << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      std::cout << "  Plane " << showerHitsIt->first << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = showerHitsIt->second.begin(); hitIt != showerHitsIt->second.end(); ++hitIt)
	std::cout << "    Hit at (" << HitCoordinates(*hitIt).X() << ", " << HitCoordinates(*hitIt).Y() << ") -- real wire " << (*hitIt)->WireID() << ", hit position (" << HitPosition(*hitIt).X() << ", " << HitPosition(*hitIt).Y() << ")" << std::endl;
    }
  }

  if (fDebug > 1) {
    std::cout << "After reversing and checking isolated hits, here are the start and end points..." << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
      std::cout << "  Plane " << showerHitsIt->first << " has start (" << HitCoordinates(showerHitsIt->second.front()).X() << ", " << HitCoordinates(showerHitsIt->second.front()).Y() << ") (real wire " << showerHitsIt->second.front()->WireID() << ") and end (" << HitCoordinates(showerHitsIt->second.back()).X() << ", " << HitCoordinates(showerHitsIt->second.back()).Y() << ")" << std::endl;
  }

  // Check for views in which the shower travels almost along the wire planes
  // (shown by a small relative wire width)
  std::map<double,int> wireWidths = RelativeWireWidth(showerHitsMap);
  std::vector<int> badPlanes;
  if (fDebug > 1)
    std::cout << "Here are the relative wire widths... " << std::endl;
  for (std::map<double,int>::iterator wireWidthIt = wireWidths.begin(); wireWidthIt != wireWidths.end(); ++wireWidthIt) {
    if (fDebug > 1)
      std::cout << "Plane " << wireWidthIt->second << " has relative wire width " << wireWidthIt->first << std::endl;
    if (wireWidthIt->first < 1/(double)wireWidths.size())
      badPlanes.push_back(wireWidthIt->second);
  }

  std::map<int,std::vector<art::Ptr<recob::Hit> > > ignoredPlanes;
  if (badPlanes.size() == 1) {
    int badPlane = badPlanes.at(0);
    if (fDebug > 1)
      std::cout << "Ignoring plane " << badPlane << " when orientating" << std::endl;
    ignoredPlanes[badPlane] = showerHitsMap.at(badPlane);
    showerHitsMap.erase(showerHitsMap.find(badPlane));
  }

  // Consider all possible permutations of planes (0/1, oriented correctly/incorrectly)
  std::map<int,std::map<int,bool> > permutations = GetPlanePermutations(showerHitsMap);

  // Go through all permutations until we have a satisfactory orientation
  std::map<int,std::vector<art::Ptr<recob::Hit> > > originalShowerHitsMap = showerHitsMap;
  int n = 0;
  while (!CheckShowerHits(showerHitsMap) and n < TMath::Power(2,(int)showerHitsMap.size())) {
    if (fDebug > 1)
      std::cout << "Permutation " << n << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      std::vector<art::Ptr<recob::Hit> > hits = originalShowerHitsMap.at(showerHitsIt->first);
      if (permutations[n][showerHitsIt->first] == 1) {
  	std::reverse(hits.begin(), hits.end());
  	showerHitsMap[showerHitsIt->first] = hits;
      }
      else
  	showerHitsMap[showerHitsIt->first] = hits;
    }
    ++n;
  }

  // Go back to original if still no consistency
  if (!CheckShowerHits(showerHitsMap))
    showerHitsMap = originalShowerHitsMap;

  if (fDebug > 2) {
    std::cout << "End of OrderShowerHits: here are the order of hits:" << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHitsIt = showerHitsMap.begin(); planeHitsIt != showerHitsMap.end(); ++planeHitsIt) {
      std::cout << "  Plane " << planeHitsIt->first << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = planeHitsIt->second.begin(); hitIt != planeHitsIt->second.end(); ++hitIt)
  	std::cout << "    Hit (" << HitCoordinates(*hitIt).X() << " (real wire " << (*hitIt)->WireID() << "), " << HitCoordinates(*hitIt).Y() << ") -- pos (" << HitPosition(*hitIt).X() << ", " << HitPosition(*hitIt).Y() << ")" << std::endl;
    }
  }

  if (ignoredPlanes.size())
    showerHitsMap[ignoredPlanes.begin()->first] = ignoredPlanes.begin()->second;

  return showerHitsMap;

}

void shower::EMShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
					    std::vector<art::Ptr<recob::Hit> >& showerHits,
                                            art::Ptr<recob::Vertex> const& vertex) const {

  showerHits = FindOrderOfHits(shower);

  // Find TPC for the vertex
  double xyz[3];
  vertex->XYZ(xyz);
  geo::TPCID tpc = fGeom->FindTPCAtPosition(xyz);
  if (!tpc.isValid&&showerHits.size()) tpc = geo::TPCID(showerHits[0]->WireID());
  //std::cout<<tpc<<std::endl;
  // Find hits in the same TPC
  art::Ptr<recob::Hit> hit0, hit1;
  for (auto &hit: showerHits){
    if (hit->WireID().TPC==tpc.TPC){
      if (hit0.isNull()){
	hit0 = hit;
      }
      hit1 = hit;
    }
  }
  if (hit0.isNull()||hit1.isNull()) return;
  TVector2 coord0 = TVector2(hit0->WireID().Wire, hit0->PeakTime());
  TVector2 coord1 = TVector2(hit1->WireID().Wire, hit1->PeakTime());
  TVector2 coordvtx = TVector2(fGeom->WireCoordinate(xyz[1], xyz[2], hit0->WireID().planeID()),
			       fDetProp->ConvertXToTicks(xyz[0],  hit0->WireID().planeID()));
//  std::cout<<coord0.X()<<" "<<coord0.Y()<<std::endl;
//  std::cout<<coord1.X()<<" "<<coord1.Y()<<std::endl;
//  std::cout<<coordvtx.X()<<" "<<coordvtx.Y()<<std::endl;
//  std::cout<<hit0->WireID()<<" "<<hit1->WireID()<<std::endl;
  if ((coord1-coordvtx).Mod()<(coord0-coordvtx).Mod()){
    std::reverse(showerHits.begin(), showerHits.end());
  }
  //std::cout<<showerHits[0]->WireID()<<" "<<showerHits.back()->WireID()<<std::endl;
}

void shower::EMShowerAlg::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >const& showerHits,
					       art::Ptr<recob::Vertex> const& vertex,
                                               std::vector<art::Ptr<recob::Hit> >& trackHits) const {

  // Find TPC for the vertex
  //std::cout<<"here"<<std::endl;
  double xyz[3];
  vertex->XYZ(xyz);
  //std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;
  geo::TPCID tpc = fGeom->FindTPCAtPosition(xyz);
  //std::cout<<tpc<<std::endl;
  //vertex cannot be projected into a TPC, find the TPC that has the most hits
  if (!tpc.isValid){
    std::map<geo::TPCID, unsigned int> tpcmap;
    unsigned maxhits = 0;
    for (auto const&hit : showerHits){
      ++tpcmap[geo::TPCID(hit->WireID())];
    }
    for (auto const&t : tpcmap){
      if (t.second > maxhits){
	maxhits = t.second;
	tpc = t.first;
      }
    }
  }
  //std::cout<<tpc<<std::endl;
    //if (!tpc.isValid&&showerHits.size()) tpc = geo::TPCID(showerHits[0]->WireID());
  if (!tpc.isValid) return;
  //std::cout<<"here 1"<<std::endl;

  double parm[2];
  int fitok = 0;
  std::vector<double> wfit;
  std::vector<double> tfit;
  std::vector<double> cfit;

  for (size_t i = 0; i<fNfitpass; ++i){

    // Fit a straight line through hits
    unsigned int nhits = 0;
    for (auto &hit: showerHits){
      //std::cout<<i<<" "<<hit->WireID()<<" "<<tpc<<std::endl;
      if (hit->WireID().TPC==tpc.TPC){
	TVector2 coord = HitCoordinates(hit);
	//std::cout<<i<<" "<<hit->WireID()<<" "<<hit->PeakTime()<<std::endl;
	if (i==0||(std::abs((coord.Y()-(parm[0]+coord.X()*parm[1]))*cos(atan(parm[1])))<fToler[i-1])||fitok==1){
	  ++nhits;
	  if (nhits==fNfithits[i]+1) break;
	  wfit.push_back(coord.X());
	  tfit.push_back(coord.Y());
	//cfit.push_back(hit->Integral());
	  cfit.push_back(1.);
	  if (i==fNfitpass-1) {
	    trackHits.push_back(hit);
	  }
	//std::cout<<*hit<<std::endl;
//
//<<hit->PeakTime()<<" "<<std::abs((coord.Y()-(parm[0]+coord.X()*parm[1]))*cos(atan(parm[1])))<<std::endl;
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

}


TVector2 shower::EMShowerAlg::HitCoordinates(art::Ptr<recob::Hit> const& hit) const {

  return TVector2(GlobalWire(hit->WireID()), hit->PeakTime());

}

TVector2 shower::EMShowerAlg::HitPosition(art::Ptr<recob::Hit> const& hit) const {

  geo::PlaneID planeID = hit->WireID().planeID();

  return HitPosition(HitCoordinates(hit), planeID);

}

TVector2 shower::EMShowerAlg::HitPosition(TVector2 const& pos, geo::PlaneID planeID) const {

  return TVector2(pos.X() * fGeom->WirePitch(planeID),
		  fDetProp->ConvertTicksToX(pos.Y(), planeID));

}

double shower::EMShowerAlg::GlobalWire(const geo::WireID& wireID) const {

  double globalWire = -999;

  // Induction
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    double wireCentre[3];
    fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }

  // Collection
  else {
    // FOR COLLECTION WIRES, HARD CODE THE GEOMETRY FOR GIVEN DETECTORS
    // THIS _SHOULD_ BE TEMPORARY. GLOBAL WIRE SUPPORT IS BEING ADDED TO THE LARSOFT GEOMETRY AND SHOULD BE AVAILABLE SOON
    if (fDetector == "dune35t") {
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      if (wireID.TPC == 0 or wireID.TPC == 1) globalWire = wireID.Wire;
      else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5) globalWire = nwires + wireID.Wire;
      else if (wireID.TPC == 6 or wireID.TPC == 7) globalWire = (2*nwires) + wireID.Wire;
      else mf::LogError("BlurredClusterAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC << " (geometry" << fDetector << ")";
    }
    else if (fDetector == "dune10kt") {
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      // Detector geometry has four TPCs, two on top of each other, repeated along z...
      int block = wireID.TPC / 4;
      globalWire = (nwires*block) + wireID.Wire;
    }
    else {
      double wireCentre[3];
      fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);
      if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
      else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
    }
  }

  return globalWire;

}

std::map<double,int> shower::EMShowerAlg::RelativeWireWidth(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const {

  // Get the wire widths
  std::map<int,int> planeWireLength;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    planeWireLength[showerHitsIt->first] = TMath::Abs(HitCoordinates(showerHitsIt->second.front()).X() - HitCoordinates(showerHitsIt->second.back()).X());

  // Find the relative wire width for each plane with respect to the others
  std::map<int,double> planeOtherWireLengths;
  for (std::map<int,int>::iterator planeWireLengthIt = planeWireLength.begin(); planeWireLengthIt != planeWireLength.end(); ++planeWireLengthIt) {
    double quad = 0.;
    for (int plane = 0; plane < (int)fGeom->MaxPlanes(); ++plane)
      if (plane != planeWireLengthIt->first and planeWireLength.count(plane))
	quad += TMath::Power(planeWireLength[plane],2);
    quad = TMath::Sqrt(quad);
    planeOtherWireLengths[planeWireLengthIt->first] = planeWireLength[planeWireLengthIt->first] / (double)quad;
  }

  // Order these backwards
  std::map<double,int> wireWidthMap;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    wireWidthMap[planeOtherWireLengths.at(showerHitsIt->first)] = showerHitsIt->first;

  return wireWidthMap;

}

TVector2 shower::EMShowerAlg::ShowerDirection(const std::vector<art::Ptr<recob::Hit> >& showerHits) const {

  TVector2 pos;
  //double weight;
  double weight = 1;
  //int nhits = 0;
  double sumx=0., sumy=0., sumx2=0., sumxy=0., sumweight = 0.;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = showerHits.begin(); hit != showerHits.end(); ++hit) {
    //++nhits;
    pos = HitPosition(*hit);
    weight = TMath::Power((*hit)->Integral(),2);
    sumweight += weight;
    sumx += weight * pos.X();
    sumy += weight * pos.Y();
    sumx2 += weight * pos.X() * pos.X();
    sumxy += weight * pos.X() * pos.Y();
  }
  //double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
  double gradient = (sumweight * sumxy - sumx * sumy) / (sumweight * sumx2 - sumx * sumx);
  TVector2 direction = TVector2(1,gradient).Unit();

  return direction;

}

TVector2 shower::EMShowerAlg::ShowerCentre(const std::vector<art::Ptr<recob::Hit> >& showerHits) const {

  TVector2 pos, chargePoint = TVector2(0,0);
  double totalCharge = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = showerHits.begin(); hit != showerHits.end(); ++hit) {
    pos = HitPosition(*hit);
    chargePoint += (*hit)->Integral() * pos;
    totalCharge += (*hit)->Integral();
  }
  TVector2 centre = chargePoint / totalCharge;

  return centre;

}

double shower::EMShowerAlg::ShowerHitRMS(const std::vector<art::Ptr<recob::Hit> >& showerHits) const {

  TVector2 direction = ShowerDirection(showerHits);
  TVector2 centre = ShowerCentre(showerHits);

  std::vector<double> distanceToAxis;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitsIt = showerHits.begin(); showerHitsIt != showerHits.end(); ++showerHitsIt) {
    TVector2 proj = (HitPosition(*showerHitsIt) - centre).Proj(direction) + centre;
    distanceToAxis.push_back((HitPosition(*showerHitsIt) - proj).Mod());
  }
  double RMS = TMath::RMS(distanceToAxis.begin(), distanceToAxis.end());

  return RMS;

}

double shower::EMShowerAlg::ShowerHitRMSGradient(const std::vector<art::Ptr<recob::Hit> >& showerHits, TVector2 trueStart) const {

  // Don't forget to clean up the header file!

  // Find a rough shower 'direction' and centre
  TVector2 direction = ShowerDirection(showerHits);

  //std::cout << std::endl << "Number of hits in this shower: " << showerHits.size() << std::endl;

  // for (int i = 0; i < 100; ++i) {

  //   // Bin the hits into discreet chunks
  //   //int nSegmentHits = fNumSegmentHits;
  //   int nSegmentHits = i;
  //   int nShowerSegments = showerHits.size() / (double)nSegmentHits;
  //   //int nShowerSegments = fNumShowerSegments;
  //   double lengthOfShower = (HitPosition(showerHits.back()) - HitPosition(showerHits.front())).Mod();
  //   double lengthOfSegment = lengthOfShower / (double)nShowerSegments;
  //   std::map<int,std::vector<art::Ptr<recob::Hit> > > showerSegments;
  //   std::map<int,double> segmentCharge;
  //   for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitsIt = showerHits.begin(); showerHitsIt != showerHits.end(); ++showerHitsIt) {
  //     showerSegments[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment].push_back(*showerHitsIt);
  //     segmentCharge[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment] += (*showerHitsIt)->Integral();
  //   }

  //   TGraph* graph = new TGraph();
  //   std::vector<std::pair<int,double> > binVsRMS;

  //   // Loop over the bins to find the distribution of hits as the shower progresses
  //   for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerSegmentIt = showerSegments.begin(); showerSegmentIt != showerSegments.end(); ++showerSegmentIt) {

  //     // Get the mean position of the hits in this bin
  //     TVector2 meanPosition(0,0);
  //     for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt)
  // 	meanPosition += HitPosition(*hitInSegmentIt);
  //     meanPosition /= (double)showerSegmentIt->second.size();

  //     // Get the RMS of this bin
  //     std::vector<double> distanceToAxisBin;
  //     for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt) {
  // 	TVector2 proj = (HitPosition(*hitInSegmentIt) - meanPosition).Proj(direction) + meanPosition;
  // 	distanceToAxisBin.push_back((HitPosition(*hitInSegmentIt) - proj).Mod());
  //     }

  //     double RMSBin = TMath::RMS(distanceToAxisBin.begin(), distanceToAxisBin.end());
  //     binVsRMS.push_back(std::make_pair(showerSegmentIt->first, RMSBin));
  //     if (fMakeRMSGradientPlot)
  // 	graph->SetPoint(graph->GetN(), showerSegmentIt->first, RMSBin);

  //   }

  //   // Get the gradient of the RMS-bin plot

  //   // OLD
  //   int nhits1 = 0;
  //   double sumx1=0., sumy1=0., sumx21=0., sumxy1=0.;
  //   for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
  //     ++nhits1;
  //     sumx1 += binVsRMSIt->first;
  //     sumy1 += binVsRMSIt->second;
  //     sumx21 += binVsRMSIt->first * binVsRMSIt->first;
  //     sumxy1 += binVsRMSIt->first * binVsRMSIt->second;
  //   }
  //   double RMSgradient1 = (nhits1 * sumxy1 - sumx1 * sumy1) / (nhits1 * sumx21 - sumx1 * sumx1);

  //   // NEW
  //   double sumx=0., sumy=0., sumx2=0., sumxy=0., sumweight = 0.;
  //   for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
  //     //double weight = showerSegments.at(binVsRMSIt->first).size();
  //     double weight = TMath::Power(segmentCharge.at(binVsRMSIt->first),2);
  //     sumweight += weight;
  //     sumx += weight * binVsRMSIt->first;
  //     sumy += weight * binVsRMSIt->second;
  //     sumx2 += weight * binVsRMSIt->first * binVsRMSIt->first;
  //     sumxy += weight * binVsRMSIt->first * binVsRMSIt->second;
  //   }
  //   double RMSgradient = (sumweight * sumxy - sumx * sumy) / (sumweight * sumx2 - sumx * sumx);

  //   // PCA
  //   TPrincipal* pca = new TPrincipal(2,"");
  //   double point[2];
  //   for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
  //     if (fDebug > 2)
  // 	std::cout << "Number of hits in segment " << binVsRMSIt->first << " is " << showerSegments.at(binVsRMSIt->first).size() << std::endl;
  //     for (int nhit = 0; nhit < (int)showerSegments.at(binVsRMSIt->first).size(); ++nhit) {
  // 	point[0] = binVsRMSIt->first;
  // 	point[1] = binVsRMSIt->second;
  // 	pca->AddRow(point);
  //     }
  //   }
  //   pca->MakePrincipals();
  //   TVector2 RMSgradientPCA = TVector2( (*pca->GetEigenVectors())[0][0], (*pca->GetEigenVectors())[1][0] ).Unit();

  //   if (fDebug > 1)
  //     std::cout << "RMS gradient without weighting: " << RMSgradient1 << ", and with weighting: " << RMSgradient << "; PCA with weighting: " << RMSgradientPCA.Y()/RMSgradientPCA.X() << std::endl;

  //   if (fMakeRMSGradientPlot) {
  //     TVector2 direction = TVector2(1,RMSgradient).Unit();
  //     TVector2 direction1 = TVector2(1,RMSgradient1).Unit();
  //     TCanvas* canv = new TCanvas();
  //     graph->Draw();
  //     TVector2 centre = TVector2(graph->GetMean(1), graph->GetMean(2));
  //     TLine line, line1, line2;
  //     line.SetLineColor(3);
  //     line1.SetLineColor(4);
  //     line2.SetLineColor(5);
  //     line.DrawLine(centre.X()-1000*direction.X(),centre.Y()-1000*direction.Y(),centre.X()+1000*direction.X(),centre.Y()+1000*direction.Y());
  //     line1.DrawLine(centre.X()-1000*direction1.X(),centre.Y()-1000*direction1.Y(),centre.X()+1000*direction1.X(),centre.Y()+1000*direction1.Y());
  //     line2.DrawLine(centre.X()-1000*RMSgradientPCA.X(),centre.Y()-1000*RMSgradientPCA.Y(),centre.X()+1000*RMSgradientPCA.X(),centre.Y()+1000*RMSgradientPCA.Y());
  //     canv->SaveAs("RMSGradient.png");
  //     delete canv;
  //   }
  //   delete graph; delete pca;

  //   //std::cout << "Number of hits per bin " << i << " gives gradient of " << RMSgradient << std::endl;

  //   double distanceFromStart = 0, distanceFromEnd = 0;
  //   if (RMSgradient > 0) {
  //     distanceFromStart = (trueStart - HitPosition(showerHits.front())).Mod();
  //     distanceFromEnd = (trueStart - HitPosition(showerHits.back())).Mod();
  //   }
  //   if (RMSgradient < 0) {
  //    distanceFromStart = (trueStart - HitPosition(showerHits.back())).Mod();
  //    distanceFromEnd = (trueStart - HitPosition(showerHits.front())).Mod();
  //   }
  //   bool correct = distanceFromStart < distanceFromEnd;
  //   hNumHitsInSegment->Fill(i,correct);

  // }




  // for (int i = 0; i < 100; ++i) {

  //   // Bin the hits into discreet chunks
  //   int nShowerSegments = i;
  //   double lengthOfShower = (HitPosition(showerHits.back()) - HitPosition(showerHits.front())).Mod();
  //   double lengthOfSegment = lengthOfShower / (double)nShowerSegments;
  //   std::map<int,std::vector<art::Ptr<recob::Hit> > > showerSegments;
  //   std::map<int,double> segmentCharge;
  //   for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitsIt = showerHits.begin(); showerHitsIt != showerHits.end(); ++showerHitsIt) {
  //     showerSegments[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment].push_back(*showerHitsIt);
  //     segmentCharge[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment] += (*showerHitsIt)->Integral();
  //   }

  //   TGraph* graph = new TGraph();
  //   std::vector<std::pair<int,double> > binVsRMS;

  //   // Loop over the bins to find the distribution of hits as the shower progresses
  //   for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerSegmentIt = showerSegments.begin(); showerSegmentIt != showerSegments.end(); ++showerSegmentIt) {

  //     // Get the mean position of the hits in this bin
  //     TVector2 meanPosition(0,0);
  //     for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt)
  // 	meanPosition += HitPosition(*hitInSegmentIt);
  //     meanPosition /= (double)showerSegmentIt->second.size();

  //     // Get the RMS of this bin
  //     std::vector<double> distanceToAxisBin;
  //     for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt) {
  // 	TVector2 proj = (HitPosition(*hitInSegmentIt) - meanPosition).Proj(direction) + meanPosition;
  // 	distanceToAxisBin.push_back((HitPosition(*hitInSegmentIt) - proj).Mod());
  //     }

  //     double RMSBin = TMath::RMS(distanceToAxisBin.begin(), distanceToAxisBin.end());
  //     binVsRMS.push_back(std::make_pair(showerSegmentIt->first, RMSBin));
  //     if (fMakeRMSGradientPlot)
  // 	graph->SetPoint(graph->GetN(), showerSegmentIt->first, RMSBin);

  //   }

  //   // Get the gradient of the RMS-bin plot

  //   // OLD
  //   int nhits1 = 0;
  //   double sumx1=0., sumy1=0., sumx21=0., sumxy1=0.;
  //   for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
  //     ++nhits1;
  //     sumx1 += binVsRMSIt->first;
  //     sumy1 += binVsRMSIt->second;
  //     sumx21 += binVsRMSIt->first * binVsRMSIt->first;
  //     sumxy1 += binVsRMSIt->first * binVsRMSIt->second;
  //   }
  //   double RMSgradient1 = (nhits1 * sumxy1 - sumx1 * sumy1) / (nhits1 * sumx21 - sumx1 * sumx1);

  //   // NEW
  //   double sumx=0., sumy=0., sumx2=0., sumxy=0., sumweight = 0.;
  //   for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
  //     //double weight = showerSegments.at(binVsRMSIt->first).size();
  //     double weight = TMath::Power(segmentCharge.at(binVsRMSIt->first),2);
  //     sumweight += weight;
  //     sumx += weight * binVsRMSIt->first;
  //     sumy += weight * binVsRMSIt->second;
  //     sumx2 += weight * binVsRMSIt->first * binVsRMSIt->first;
  //     sumxy += weight * binVsRMSIt->first * binVsRMSIt->second;
  //   }
  //   double RMSgradient = (sumweight * sumxy - sumx * sumy) / (sumweight * sumx2 - sumx * sumx);

  //   // PCA
  //   TPrincipal* pca = new TPrincipal(2,"");
  //   double point[2];
  //   for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
  //     if (fDebug > 2)
  // 	std::cout << "Number of hits in segment " << binVsRMSIt->first << " is " << showerSegments.at(binVsRMSIt->first).size() << std::endl;
  //     for (int nhit = 0; nhit < (int)showerSegments.at(binVsRMSIt->first).size(); ++nhit) {
  // 	point[0] = binVsRMSIt->first;
  // 	point[1] = binVsRMSIt->second;
  // 	pca->AddRow(point);
  //     }
  //   }
  //   pca->MakePrincipals();
  //   TVector2 RMSgradientPCA = TVector2( (*pca->GetEigenVectors())[0][0], (*pca->GetEigenVectors())[1][0] ).Unit();

  //   if (fDebug > 1)
  //     std::cout << "RMS gradient without weighting: " << RMSgradient1 << ", and with weighting: " << RMSgradient << "; PCA with weighting: " << RMSgradientPCA.Y()/RMSgradientPCA.X() << std::endl;

  //   if (fMakeRMSGradientPlot) {
  //     TVector2 direction = TVector2(1,RMSgradient).Unit();
  //     TVector2 direction1 = TVector2(1,RMSgradient1).Unit();
  //     TCanvas* canv = new TCanvas();
  //     graph->Draw();
  //     TVector2 centre = TVector2(graph->GetMean(1), graph->GetMean(2));
  //     TLine line, line1, line2;
  //     line.SetLineColor(3);
  //     line1.SetLineColor(4);
  //     line2.SetLineColor(5);
  //     line.DrawLine(centre.X()-1000*direction.X(),centre.Y()-1000*direction.Y(),centre.X()+1000*direction.X(),centre.Y()+1000*direction.Y());
  //     line1.DrawLine(centre.X()-1000*direction1.X(),centre.Y()-1000*direction1.Y(),centre.X()+1000*direction1.X(),centre.Y()+1000*direction1.Y());
  //     line2.DrawLine(centre.X()-1000*RMSgradientPCA.X(),centre.Y()-1000*RMSgradientPCA.Y(),centre.X()+1000*RMSgradientPCA.X(),centre.Y()+1000*RMSgradientPCA.Y());
  //     canv->SaveAs("RMSGradient.png");
  //     delete canv;
  //   }
  //   delete graph; delete pca;

  //   //std::cout << "Number of hits per bin " << i << " gives gradient of " << RMSgradient << std::endl;

  //   double distanceFromStart = 0, distanceFromEnd = 0;
  //   if (RMSgradient > 0) {
  //     distanceFromStart = (trueStart - HitPosition(showerHits.front())).Mod();
  //     distanceFromEnd = (trueStart - HitPosition(showerHits.back())).Mod();
  //   }
  //   if (RMSgradient < 0) {
  //    distanceFromStart = (trueStart - HitPosition(showerHits.back())).Mod();
  //    distanceFromEnd = (trueStart - HitPosition(showerHits.front())).Mod();
  //   }
  //   bool correct = distanceFromStart < distanceFromEnd;
  //   hNumSegments->Fill(i,correct);

  // }
















  // Bin the hits into discreet chunks
  int nShowerSegments = fNumShowerSegments;
  double lengthOfShower = (HitPosition(showerHits.back()) - HitPosition(showerHits.front())).Mod();
  double lengthOfSegment = lengthOfShower / (double)nShowerSegments;
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerSegments;
  std::map<int,double> segmentCharge;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitsIt = showerHits.begin(); showerHitsIt != showerHits.end(); ++showerHitsIt) {
    showerSegments[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment].push_back(*showerHitsIt);
    segmentCharge[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment] += (*showerHitsIt)->Integral();
  }

  TGraph* graph = new TGraph();
  std::vector<std::pair<int,double> > binVsRMS;

  // Loop over the bins to find the distribution of hits as the shower progresses
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerSegmentIt = showerSegments.begin(); showerSegmentIt != showerSegments.end(); ++showerSegmentIt) {

    // Get the mean position of the hits in this bin
    TVector2 meanPosition(0,0);
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt)
      meanPosition += HitPosition(*hitInSegmentIt);
    meanPosition /= (double)showerSegmentIt->second.size();

    // Get the RMS of this bin
    std::vector<double> distanceToAxisBin;
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt) {
      TVector2 proj = (HitPosition(*hitInSegmentIt) - meanPosition).Proj(direction) + meanPosition;
      distanceToAxisBin.push_back((HitPosition(*hitInSegmentIt) - proj).Mod());
    }

    double RMSBin = TMath::RMS(distanceToAxisBin.begin(), distanceToAxisBin.end());
    binVsRMS.push_back(std::make_pair(showerSegmentIt->first, RMSBin));
    if (fMakeRMSGradientPlot)
      graph->SetPoint(graph->GetN(), showerSegmentIt->first, RMSBin);

  }

  // Get the gradient of the RMS-bin plot
  double sumx=0., sumy=0., sumx2=0., sumxy=0., sumweight = 0.;
  for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
    //double weight = showerSegments.at(binVsRMSIt->first).size();
    //double weight = 1;
    double weight = segmentCharge.at(binVsRMSIt->first);
    sumweight += weight;
    sumx += weight * binVsRMSIt->first;
    sumy += weight * binVsRMSIt->second;
    sumx2 += weight * binVsRMSIt->first * binVsRMSIt->first;
    sumxy += weight * binVsRMSIt->first * binVsRMSIt->second;
  }
  double RMSgradient = (sumweight * sumxy - sumx * sumy) / (sumweight * sumx2 - sumx * sumx);

  if (fMakeRMSGradientPlot) {
    TVector2 direction = TVector2(1,RMSgradient).Unit();
    TCanvas* canv = new TCanvas();
    graph->Draw();
    graph->GetXaxis()->SetTitle("Shower segment");
    graph->GetYaxis()->SetTitle("RMS of hit distribution");
    TVector2 centre = TVector2(graph->GetMean(1), graph->GetMean(2));
    TLine line;
    line.SetLineColor(2);
    line.DrawLine(centre.X()-1000*direction.X(),centre.Y()-1000*direction.Y(),centre.X()+1000*direction.X(),centre.Y()+1000*direction.Y());
    canv->SaveAs("RMSGradient.png");
    delete canv;
  }
  delete graph;

  return RMSgradient;

}

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane(TVector3 const& point, int plane, int cryostat) const {

  TVector2 wireTickPos = TVector2(-999., -999.);

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  geo::TPCID tpcID = fGeom->FindTPCAtPosition(pointPosition);
  int tpc = 0;
  if (tpcID.isValid)
    tpc = tpcID.TPC;
  else
    return wireTickPos;

  // Construct wire ID for this point projected onto the plane
  geo::PlaneID planeID = geo::PlaneID(cryostat, tpc, plane);
  geo::WireID wireID;
  try{
    wireID = fGeom->NearestWireID(point, planeID);
  }
  catch(geo::InvalidWireError const& e) {
    wireID = e.suggestedWireID(); // pick the closest valid wire
  }

  wireTickPos = TVector2(GlobalWire(wireID),
                         fDetProp->ConvertXToTicks(point.X(), planeID));

  // wireTickPos = TVector2(fGeom->WireCoordinate(point.Y(), point.Z(), planeID.Plane, tpc % 2, planeID.Cryostat),
  // 			 fDetProp->ConvertXToTicks(point.X(), planeID.Plane, tpc % 2, planeID.Cryostat));

  //return wireTickPos;
  return HitPosition(wireTickPos, planeID);

}

int shower::EMShowerAlg::WorstPlane(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) const {

  // Get the width of the shower in wire coordinate
  std::map<int,int> planeWireLength;
  std::map<int,double> planeOtherWireLengths;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    planeWireLength[showerHitsIt->first] = TMath::Abs(HitCoordinates(showerHitsIt->second.front()).X() - HitCoordinates(showerHitsIt->second.back()).X());
  for (std::map<int,int>::iterator planeWireLengthIt = planeWireLength.begin(); planeWireLengthIt != planeWireLength.end(); ++planeWireLengthIt) {
    double quad = 0.;
    for (int plane = 0; plane < (int)fGeom->MaxPlanes(); ++plane)
      if (plane != planeWireLengthIt->first and planeWireLength.count(plane))
	quad += TMath::Power(planeWireLength[plane],2);
    quad = TMath::Sqrt(quad);
    planeOtherWireLengths[planeWireLengthIt->first] = planeWireLength[planeWireLengthIt->first] / (double)quad;
  }

  if (fDebug > 1)
    for (std::map<int,double>::const_iterator planeOtherWireLengthIt = planeOtherWireLengths.begin(); planeOtherWireLengthIt != planeOtherWireLengths.end(); ++planeOtherWireLengthIt)
      std::cout << "Plane " << planeOtherWireLengthIt->first << " has " << planeWireLength[planeOtherWireLengthIt->first] << " wire width and therefore has relative width in wire of " << planeOtherWireLengthIt->second << std::endl;

  std::map<double,int> wireWidthMap;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    double wireWidth = planeWireLength.at(showerHitsIt->first);
    wireWidthMap[wireWidth] = showerHitsIt->first;
  }

  return wireWidthMap.begin()->second;

}

Int_t shower::EMShowerAlg::WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm) const {

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

bool shower::EMShowerAlg::isCleanShower(std::vector<art::Ptr<recob::Hit> > const& hits) const {

  if (!hits.size()) return false;
  if (hits.size()>2000) return true;
  if (hits.size()<20) return true;
  std::map<int, int> hitmap;
  unsigned nhits = 0;
  for (auto const&hit : hits){
    ++nhits;
    if (nhits>2)
      ++hitmap[hit->WireID().Wire];
    if (nhits==20) break;
  }
  //std::cout<<hits.size()<<" "<<float(nhits-2)/hitmap.size()<<std::endl;
  if (float(nhits-2)/hitmap.size()>1.4) return false;
  else return true;
}

shower::HitPosition::HitPosition()
  : fGeom(lar::providerFrom<geo::Geometry>())
  , fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {}


// // Code to make wire maps showing the global wire coordinates
// struct TPCWire {
//   int fPlane, fTPC, fWire, fGlobalWire;
//   TVector3 fStart, fEnd;
//   TPCWire(int plane, int tpc, int wire, int globalWire, TVector3 start, TVector3 end) {
//     fPlane = plane;
//     fTPC = tpc;
//     fWire = wire;
//     fGlobalWire = globalWire;
//     fStart = start;
//     fEnd = end;
//   }
//   TPCWire() { }
//   void SetProps(int plane, int tpc, int wire, int globalWire, TVector3 start, TVector3 end) {
//     fPlane = plane;
//     fTPC = tpc;
//     fWire = wire;
//     fGlobalWire = globalWire;
//     fStart = start;
//     fEnd = end;
//   }
// };

// void shower::EMShowerAlg::MakePicture() {

//   std::vector<TPCWire> allWires;

//   for (geo::WireID const& wireID : fGeom->IterateWireIDs()) {

//     if (wireID.TPC % 2 == 0)
//       continue;

//     double xyzStart[3], xyzEnd[3];
//     fGeom->WireEndPoints(wireID, xyzStart, xyzEnd);
//     int globalWire = GlobalWire(wireID);

//     allWires.emplace_back(wireID.Plane, wireID.TPC, wireID.Wire, globalWire, TVector3(xyzStart[0],xyzStart[1],xyzStart[2]), TVector3(xyzEnd[0],xyzEnd[1],xyzEnd[2]));

//   } // for all wires

//   TCanvas* uplane = new TCanvas();
//   TCanvas* vplane = new TCanvas();
//   TCanvas* zplane = new TCanvas();
//   TText* number = new TText();
//   number->SetTextSize(0.002);
//   TLine* line = new TLine(0,100,0,100);
//   line->SetLineWidth(0.5);
//   uplane->Range(-5,-100,160,140);
//   vplane->Range(-5,-100,160,140);
//   zplane->Range(-5,-100,160,140);

//   for (std::vector<TPCWire>::iterator wireIt = allWires.begin(); wireIt != allWires.end(); ++wireIt) {
//     if (wireIt->fPlane == 0)
//       uplane->cd();
//     else if (wireIt->fPlane == 1)
//       vplane->cd();
//     else if (wireIt->fPlane == 2)
//       zplane->cd();
//     line->DrawLine(wireIt->fStart.Z(), wireIt->fStart.Y(), wireIt->fEnd.Z(), wireIt->fEnd.Y());
//     number->DrawText(wireIt->fStart.Z(), wireIt->fStart.Y(), (std::to_string(wireIt->fGlobalWire)+" ("+std::to_string(wireIt->fWire)+")").c_str());
//     number->DrawText(wireIt->fEnd.Z(), wireIt->fEnd.Y(), (std::to_string(wireIt->fGlobalWire)+" ("+std::to_string(wireIt->fWire)+")").c_str());
//   }

//   uplane->SaveAs("UPlane.pdf");
//   vplane->SaveAs("VPlane.pdf");
//   zplane->SaveAs("ZPlane.pdf");

// }

//
// // Timing code...
//
// auto start_time = std::chrono::high_resolution_clock::now();
// // Put stuff here!
// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count();
// std::cout << "Duration is " << duration/1000000.0 << " s " << std::endl;;
//

int shower::EMShowerAlg::FindTrueParticle(const std::vector<art::Ptr<recob::Hit> >& showerHits) const {

  /// Returns the true particle most likely associated with this shower

  // Make a map of the tracks which are associated with this shower and the charge each contributes
  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitIt = showerHits.begin(); showerHitIt != showerHits.end(); ++showerHitIt) {
    art::Ptr<recob::Hit> hit = *showerHitIt;
    int trackID = FindParticleID(hit);
    trackMap[trackID] += hit->Integral();
  }

  // Pick the track with the highest charge as the 'true track'
  double highestCharge = 0;
  int showerTrack = 0;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      showerTrack  = trackIt->first;
    }
  }

  return showerTrack;

}

int shower::EMShowerAlg::FindParticleID(const art::Ptr<recob::Hit>& hit) const {

  /// Returns the true track ID associated with this hit (if more than one, returns the one with highest energy)

  double particleEnergy = 0;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }

  return likelyTrackID;

}

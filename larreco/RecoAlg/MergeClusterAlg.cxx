////////////////////////////////////////////////////////////////////
// Merge Cluster algorithm
//
// Runs on the output of previous clustering algorithms to merge
// clusters together which lie on a straight line and are within
// some separation threshold.
// Runs recursively over all clusters, including new ones formed
// in the algorithm.
//
// M Wallbank (m.wallbank@sheffield.ac.uk), July 2015
////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/MergeClusterAlg.h"

#include "TPrincipal.h"
#include "TTree.h"

cluster::MergeClusterAlg::MergeClusterAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  fTree = tfs->make<TTree>("MatchingVariables", "MatchingVariables");
  fTree->Branch("Angle", &fAngle);
  fTree->Branch("Eigenvalue", &fEigenvalue);
  fTree->Branch("Cluster1Size", &fCluster1Size);
  fTree->Branch("Cluster2Size", &fCluster2Size);
  fTree->Branch("Length1", &fLength1);
  fTree->Branch("Length2", &fLength2);
  fTree->Branch("Separation", &fSeparation);
  fTree->Branch("CrossingDistance", &fCrossingDistance);
  fTree->Branch("ProjectedWidth", &fProjectedWidth);
  fTree->Branch("Overlap", &fOverlap);
  fTree->Branch("TrueMerge", &fTrueMerge);
}

void cluster::MergeClusterAlg::FindClusterEndPoints(art::PtrVector<recob::Hit> const& cluster,
                                                    TVector2 const& centre,
                                                    TVector2 const& direction,
                                                    TVector2& start,
                                                    TVector2& end) const
{

  /// Find estimates of cluster start/end points

  TVector2 pos;
  std::map<double, TVector2> hitProjection;

  // Project all hits onto line to determine end points
  for (auto& hit : cluster) {
    pos = HitCoordinates(hit) - centre;
    hitProjection[direction * pos] = pos;
  }

  // Project end points onto line which passes through centre of cluster
  start = hitProjection.begin()->second.Proj(direction) + centre;
  end = hitProjection.rbegin()->second.Proj(direction) + centre;

  return;
}

double cluster::MergeClusterAlg::FindClusterOverlap(TVector2 const& direction,
                                                    TVector2 const& centre,
                                                    TVector2 const& start1,
                                                    TVector2 const& end1,
                                                    TVector2 const& start2,
                                                    TVector2 const& end2) const
{

  /// Calculates the overlap of the clusters on the line projected between them

  double clusterOverlap = 0;

  // Project onto the average direction through both clusters
  double s1 = (start1 - centre) * direction;
  double e1 = (end1 - centre) * direction;
  double s2 = (start2 - centre) * direction;
  double e2 = (end2 - centre) * direction;

  // Make sure end > start
  if (s1 > e1) {
    std::cout << "s1>e1: " << s1 << " and " << e1 << std::endl;
    double tmp = e1;
    e1 = s1;
    s1 = tmp;
  }
  if (s2 > e2) {
    std::cout << "s1>e1: " << s1 << " and " << e1 << std::endl;
    double tmp = e2;
    e2 = s2;
    s2 = tmp;
  }

  // Find the overlap of the clusters on the centre line
  if ((e1 > s2) && (e2 > s1)) clusterOverlap = std::min((e1 - s2), (e2 - s1));

  return clusterOverlap;
}

double cluster::MergeClusterAlg::FindCrossingDistance(TVector2 const& direction1,
                                                      TVector2 const& centre1,
                                                      TVector2 const& direction2,
                                                      TVector2 const& centre2) const
{

  /// Finds the distance between the crossing point of the lines and the closest line centre

  // Find intersection point of two lines drawn through the centre of the clusters
  double dcross = (direction1.X() * direction2.Y()) - (direction1.Y() * direction2.X());
  TVector2 p = centre2 - centre1;
  double pcrossd = (p.X() * direction2.Y()) - (p.Y() * direction2.X());
  TVector2 crossing = centre1 + ((pcrossd / dcross) * direction1);

  // Get distance from this point to the clusters
  double crossingDistance = std::min((centre1 - crossing).Mod(), (centre2 - crossing).Mod());

  return crossingDistance;
}

double cluster::MergeClusterAlg::FindMinSeparation(art::PtrVector<recob::Hit> const& cluster1,
                                                   art::PtrVector<recob::Hit> const& cluster2) const
{

  /// Calculates the minimum separation between two clusters

  double minDistance = 99999.;

  // Loop over the two clusters to find the smallest distance
  for (auto const& hit1 : cluster1) {
    for (auto const& hit2 : cluster2) {

      TVector2 pos1 = HitCoordinates(hit1);
      TVector2 pos2 = HitCoordinates(hit2);

      double distance = (pos1 - pos2).Mod();

      if (distance < minDistance) minDistance = distance;
    }
  }

  return minDistance;
}

double cluster::MergeClusterAlg::FindProjectedWidth(TVector2 const& centre1,
                                                    TVector2 const& start1,
                                                    TVector2 const& end1,
                                                    TVector2 const& centre2,
                                                    TVector2 const& start2,
                                                    TVector2 const& end2) const
{

  /// Projects clusters parallel to the line which runs through their centres and finds the minimum containing width

  // Get the line running through the centre of the two clusters
  TVector2 parallel = (centre2 - centre1).Unit();
  TVector2 perpendicular = parallel.Rotate(TMath::Pi() / 2);

  // Project the cluster vector onto this perpendicular line
  double s1 = (start1 - centre1) * perpendicular;
  double e1 = (end1 - centre1) * perpendicular;
  double s2 = (start2 - centre2) * perpendicular;
  double e2 = (end2 - centre2) * perpendicular;

  // Find the width in each direction
  double projectionStart = std::max(TMath::Abs(s1), TMath::Abs(s2));
  double projectionEnd = std::max(TMath::Abs(e1), TMath::Abs(e2));

  double projectionWidth = projectionStart + projectionEnd;

  return projectionWidth;
}

double cluster::MergeClusterAlg::GlobalWire(geo::WireID const& wireID) const
{

  /// Find the global wire position

  auto const wireCenter = fGeom->WireIDToWireGeo(wireID).GetCenter<geo::Point_t>();

  if (fGeom->SignalType(wireID) == geo::kInduction) {
    return fGeom->WireCoordinate(wireCenter,
                                 geo::PlaneID{wireID.Cryostat, wireID.TPC % 2, wireID.Plane});
  }

  double globalWire;
  if (wireID.TPC % 2 == 0)
    globalWire = wireID.Wire + ((wireID.TPC / 2) * fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat));
  else
    globalWire =
      wireID.Wire + ((int)(wireID.TPC / 2) * fGeom->Nwires(wireID.Plane, 1, wireID.Cryostat));

  return globalWire;
}

TVector2 cluster::MergeClusterAlg::HitCoordinates(art::Ptr<recob::Hit> const& hit) const
{

  /// Return the coordinates of this hit in global wire/tick space

  return TVector2(GlobalWire(hit->WireID()), hit->PeakTime());
}

int cluster::MergeClusterAlg::MergeClusters(
  std::vector<art::PtrVector<recob::Hit>> const& planeClusters,
  std::vector<art::PtrVector<recob::Hit>>& clusters) const
{

  /// Merges clusters which lie along a straight line

  // // ////// MAKE SOME MESSY CODE! CHECK FEATURES OF THESE CLUSTERS

  // // Truth matching
  // for (unsigned int cluster = 0; cluster < planeClusters.size(); ++cluster) {
  //   std::map<int,double> trackMap;
  //   art::PtrVector<recob::Hit> hits = planeClusters.at(cluster);
  //   for (auto &hit : hits) {
  //     std::vector<sim::TrackIDE> ides = backtracker->HitToTrackID(hit);
  //     for (auto &ide : ides)
  //    trackMap[ide.trackID] += ide.energy;
  //   }
  //   // Find the true particle associated with this track
  //   double highEnergy = 0;
  //   int bestTrack = 0;
  //   for (auto &track : trackMap) {
  //     if (track.second > highEnergy) {
  //    highEnergy = track.second;
  //    bestTrack = track.first;
  //     }
  //   }
  //   trueClusterMap[cluster] = bestTrack;
  // }

  // for (unsigned int cluster1It = 0; cluster1It < planeClusters.size(); ++cluster1It) {
  //   for (unsigned int cluster2It = cluster1It+1; cluster2It < planeClusters.size(); ++cluster2It) {

  //     const art::PtrVector<recob::Hit> cluster1 = planeClusters.at(cluster1It);
  //     const art::PtrVector<recob::Hit> cluster2 = planeClusters.at(cluster2It);

  //     // true merge
  //     if (trueClusterMap[cluster1It] == trueClusterMap[cluster2It])
  //    fTrueMerge = true;
  //     else fTrueMerge = false;

  //     // geometry
  //     fCluster1Size = cluster1.size();
  //     fCluster2Size = cluster2.size();
  //     fSeparation = this->FindMinSeparation(cluster1, cluster2);

  //     // PCA
  //     TPrincipal *pca = new TPrincipal(2,"");
  //     TPrincipal *pca1 = new TPrincipal(2,"");
  //     TPrincipal *pca2 = new TPrincipal(2,"");
  //     double hits[2];
  //     TVector2 pos;

  //     // Cluster centre
  //     TVector2 chargePoint1 = TVector2(0,0), chargePoint2 = TVector2(0,0);
  //     double totalCharge1 = 0, totalCharge2 = 0;

  //     for (auto &hit1 : cluster1) {
  //    pos = HitCoordinates(hit1);
  //    hits[0] = pos.X();
  //    hits[1] = pos.Y();
  //    pca->AddRow(hits);
  //    pca1->AddRow(hits);
  //    chargePoint1 += hit1->Integral() * pos;
  //    totalCharge1 += hit1->Integral();
  //     }
  //     for (auto &hit2 : cluster2) {
  //    pos = HitCoordinates(hit2);
  //    hits[0] = pos.X();
  //    hits[1] = pos.Y();
  //    pca->AddRow(hits);
  //    pca2->AddRow(hits);
  //    chargePoint2 += hit2->Integral() * pos;
  //    totalCharge2 += hit2->Integral();
  //     }

  //     pca->MakePrincipals();
  //     pca1->MakePrincipals();
  //     pca2->MakePrincipals();

  //     // Properties of these clusters
  //     TVector2 direction1 = TVector2( (*pca1->GetEigenVectors())[0][0], (*pca1->GetEigenVectors())[1][0] ).Unit();
  //     TVector2 direction2 = TVector2( (*pca2->GetEigenVectors())[0][0], (*pca2->GetEigenVectors())[1][0] ).Unit();
  //     TVector2 directionAv = ((direction1+direction2)/2).Unit();
  //     TVector2 centre1 = chargePoint1 / totalCharge1;
  //     TVector2 centre2 = chargePoint2 / totalCharge2;
  //     TVector2 centre = (centre1+centre2)/2;
  //     TVector2 start1, end1;
  //     TVector2 start2, end2;
  //     FindClusterEndPoints(cluster1, centre1, direction1, start1, end1);
  //     FindClusterEndPoints(cluster2, centre2, direction2, start2, end2);
  //     fLength1 = (end1-start1).Mod();
  //     fLength2 = (end2-start2).Mod();

  //     // Properties of the pair of clusters
  //     fCrossingDistance = FindCrossingDistance(direction1, centre1, direction2, centre2);
  //     fProjectedWidth = FindProjectedWidth(centre1, start1, end1, centre2, start2, end2);
  //     fAngle = direction1.DeltaPhi(direction2);
  //     if (fAngle > 1.57) fAngle = 3.14159 - fAngle;
  //     fOverlap = FindClusterOverlap(directionAv, centre, start1, end1, start2, end2);
  //     fSeparation = FindMinSeparation(cluster1, cluster2);
  //     fEigenvalue = (*pca->GetEigenValues())[0];

  //     fTree->Fill();

  //     // std::cout << std::endl << "Plane " << fPlane << ": Clusters " << cluster1It << " and " << cluster2It << " have overlap " << fOverlap << " and start and end ... " << std::endl;
  //     // start1.Print();
  //     // end1.Print();
  //     // start2.Print();
  //     // end2.Print();

  //     // // Find if this is merged!
  //     // if (fCrossingDistance < 6 + (5 / (fAngle - 0.05)))
  //     //     fMerge = true;
  //     // else fMerge = false;

  //     // if (fCluster1Size >= 10 && fCluster2Size >= 10) std::cout << "Merge " << fMerge << " and true merge " << fTrueMerge << std::endl;

  //   }
  // }

  // ----------------------------- END OF MESSY CODE! --------------------------------------------------------------------------------------------------------------

  std::vector<unsigned int> mergedClusters;

  std::vector<art::PtrVector<recob::Hit>> oldClusters = planeClusters;

  // Sort the clusters by size
  std::sort(oldClusters.begin(),
            oldClusters.end(),
            [](const art::PtrVector<recob::Hit>& a, const art::PtrVector<recob::Hit>& b) {
              return a.size() > b.size();
            });

  // Find the numbers of clusters above size threshold
  unsigned int nclusters = 0;
  for (auto& cluster : oldClusters)
    if (cluster.size() >= fMinMergeClusterSize) ++nclusters;

  // Until all clusters are merged, create new clusters
  bool mergedAllClusters = false;
  while (!mergedAllClusters) {

    // New cluster
    art::PtrVector<recob::Hit> cluster;

    // Put the largest unmerged cluster in this new cluster
    for (unsigned int initCluster = 0; initCluster < oldClusters.size(); ++initCluster) {
      if (oldClusters.at(initCluster).size() < fMinMergeClusterSize or
          std::find(mergedClusters.begin(), mergedClusters.end(), initCluster) !=
            mergedClusters.end())
        continue;
      cluster = oldClusters.at(initCluster);
      mergedClusters.push_back(initCluster);
      break;
    }

    // Merge all aligned clusters to this
    bool mergedAllToThisCluster = false;
    while (!mergedAllToThisCluster) {

      // Look at all clusters and merge
      int nadded = 0;
      for (unsigned int trialCluster = 0; trialCluster < oldClusters.size(); ++trialCluster) {

        if (oldClusters.at(trialCluster).size() < fMinMergeClusterSize or
            std::find(mergedClusters.begin(), mergedClusters.end(), trialCluster) !=
              mergedClusters.end())
          continue;

        // Calculate the PCA for each
        TPrincipal *pca1 = new TPrincipal(2, ""), *pca2 = new TPrincipal(2, "");
        double hits[2];
        TVector2 pos;

        // Cluster centre
        TVector2 chargePoint1 = TVector2(0, 0), chargePoint2 = TVector2(0, 0);
        double totalCharge1 = 0, totalCharge2 = 0;

        for (auto& hit1 : cluster) {
          pos = HitCoordinates(hit1);
          hits[0] = pos.X();
          hits[1] = pos.Y();
          pca1->AddRow(hits);
          chargePoint1 += hit1->Integral() * pos;
          totalCharge1 += hit1->Integral();
        }
        for (auto& hit2 : oldClusters.at(trialCluster)) {
          pos = HitCoordinates(hit2);
          hits[0] = pos.X();
          hits[1] = pos.Y();
          pca2->AddRow(hits);
          chargePoint2 += hit2->Integral() * pos;
          totalCharge2 += hit2->Integral();
        }

        pca1->MakePrincipals();
        pca2->MakePrincipals();

        // Properties of these clusters
        TVector2 direction1 =
          TVector2((*pca1->GetEigenVectors())[0][0], (*pca1->GetEigenVectors())[1][0]).Unit();
        TVector2 direction2 =
          TVector2((*pca2->GetEigenVectors())[0][0], (*pca2->GetEigenVectors())[1][0]).Unit();
        TVector2 direction = ((direction1 + direction2) / 2).Unit();
        TVector2 centre1 = chargePoint1 / totalCharge1;
        TVector2 centre2 = chargePoint2 / totalCharge2;
        TVector2 centre = (centre1 + centre2) / 2;
        TVector2 start1, end1;
        TVector2 start2, end2;
        FindClusterEndPoints(cluster, centre1, direction1, start1, end1);
        FindClusterEndPoints(oldClusters.at(trialCluster), centre2, direction2, start2, end2);
        double length1 = (end1 - start1).Mod();
        double length2 = (end2 - start2).Mod();

        // Properties of the pair of clusters
        double crossingDistance = FindCrossingDistance(direction1, centre1, direction2, centre2);
        double projectedWidth = FindProjectedWidth(centre1, start1, end1, centre2, start2, end2);
        double angle = direction1.DeltaPhi(direction2);
        if (angle > 1.57) angle = 3.14159 - angle;
        double overlap = FindClusterOverlap(direction, centre, start1, end1, start2, end2);
        double separation = FindMinSeparation(cluster, oldClusters.at(trialCluster));

        if (separation > fMaxMergeSeparation) continue;
        if (PassCuts(angle,
                     crossingDistance,
                     projectedWidth,
                     separation,
                     overlap,
                     TMath::Max(length1, length2))) {

          for (auto& hit : oldClusters.at(trialCluster))
            cluster.push_back(hit);

          mergedClusters.push_back(trialCluster);
          ++nadded;
        }

        delete pca1;
        delete pca2;

      } // loop over clusters to add

      if (nadded == 0) mergedAllToThisCluster = true;

    } // while loop

    clusters.push_back(cluster);
    if (mergedClusters.size() == nclusters) mergedAllClusters = true;
  }

  return clusters.size();
}

bool cluster::MergeClusterAlg::PassCuts(double const& angle,
                                        double const& crossingDistance,
                                        double const& projectedWidth,
                                        double const& separation,
                                        double const& overlap,
                                        double const& longLength) const
{

  /// Boolean function which decides whether or not two clusters should be merged, depending on their properties

  bool passCrossingDistanceAngle = false;
  if (crossingDistance < (-2 + (5 / (1 * TMath::Abs(angle)) - 0))) passCrossingDistanceAngle = true;

  bool passSeparationAngle = false;
  if (separation < (200 * TMath::Abs(angle) + 40)) passSeparationAngle = true;

  bool passProjectedWidth = false;
  if (((double)projectedWidth / (double)longLength) < fProjWidthThreshold)
    passProjectedWidth = true;

  return passCrossingDistanceAngle and passSeparationAngle and passProjectedWidth;
}

void cluster::MergeClusterAlg::reconfigure(fhicl::ParameterSet const& p)
{
  fMinMergeClusterSize = p.get<int>("MinMergeClusterSize");
  fMaxMergeSeparation = p.get<double>("MaxMergeSeparation");
  fProjWidthThreshold = p.get<double>("ProjWidthThreshold");
}

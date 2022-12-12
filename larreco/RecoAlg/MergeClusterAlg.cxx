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
#include "larcore/Geometry/ExptGeoHelperInterface.h"

#include "TPrincipal.h"
#include "TTree.h"

cluster::MergeClusterAlg::MergeClusterAlg(fhicl::ParameterSet const& pset)
  : fChannelMapAlg{art::ServiceHandle<geo::ExptGeoHelperInterface const>()->ChannelMapAlgPtr()}
{
  fMinMergeClusterSize = pset.get<int>("MinMergeClusterSize");
  fMaxMergeSeparation = pset.get<double>("MaxMergeSeparation");
  fProjWidthThreshold = pset.get<double>("ProjWidthThreshold");
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

  std::map<double, TVector2> hitProjection;

  // Project all hits onto line to determine end points
  for (auto& hit : cluster) {
    TVector2 pos = HitCoordinates(hit) - centre;
    hitProjection[direction * pos] = pos;
  }

  // Project end points onto line which passes through centre of cluster
  start = hitProjection.begin()->second.Proj(direction) + centre;
  end = hitProjection.rbegin()->second.Proj(direction) + centre;
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
  return std::min((centre1 - crossing).Mod(), (centre2 - crossing).Mod());
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

  return projectionStart + projectionEnd; // FIXME (KJK): Really?  The width is the start + the end?
}

double cluster::MergeClusterAlg::GlobalWire(geo::WireID const& wireID) const
{
  /// Find the global wire position

  auto const wireCenter = fGeom->Wire(wireID).GetCenter();
  geo::PlaneID const planeID{wireID.Cryostat, wireID.TPC % 2, wireID.Plane};

  if (fChannelMapAlg->SignalType(wireID) == geo::kInduction) {
    return fGeom->Plane(planeID).WireCoordinate(wireCenter);
  }
  return wireID.Wire + ((wireID.TPC / 2) * fGeom->Nwires(planeID));
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

bool cluster::MergeClusterAlg::PassCuts(double const angle,
                                        double const crossingDistance,
                                        double const projectedWidth,
                                        double const separation,
                                        double const overlap,
                                        double const longLength) const
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

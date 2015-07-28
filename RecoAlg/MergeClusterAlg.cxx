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

#include "RecoAlg/MergeClusterAlg.h"

cluster::MergeClusterAlg::MergeClusterAlg(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
}

void cluster::MergeClusterAlg::reconfigure(fhicl::ParameterSet const& p) {
  fMinMergeClusterSize = p.get<int>   ("MinMergeClusterSize");
  fMaxMergeSeparation  = p.get<double>("MaxMergeSeparation");
  fMergingThreshold    = p.get<double>("MergingThreshold");
}

TVector2 cluster::MergeClusterAlg::ConvertWireDriftToCm(unsigned int wire, float drift, unsigned int plane, unsigned int tpc, unsigned int cryo) {

  /// Convert the wire/tick coordinates roughly into cm
  /// Taken from RecoAlg/PMAlg/Utilities.cxx (written by D.Stefan & R.Sulej)

  return TVector2(fGeom->TPC(tpc, cryo).Plane(plane).WirePitch() * wire,
		  fDetProp->ConvertTicksToX(drift, plane, tpc, cryo)
		  );
}

int cluster::MergeClusterAlg::MergeClusters(std::vector<art::PtrVector<recob::Hit> > *planeClusters, std::vector<art::PtrVector<recob::Hit> > &clusters, unsigned int plane, unsigned int tpc, unsigned int cryo) {

  /// Merges clusters which lie along a straight line

  std::vector<unsigned int> mergedClusters;

  // Sort the clusters by size
  std::sort(planeClusters->begin(), planeClusters->end(), [](const art::PtrVector<recob::Hit> &a, const art::PtrVector<recob::Hit> &b) {return a.size() > b.size();} );

  // Find the numbers of clusters above size threshold
  unsigned int nclusters = 0;
  for (auto &cluster : *planeClusters)
    if (cluster.size() >= fMinMergeClusterSize) ++nclusters;

  // Until all clusters are merged, create new clusters
  bool mergedAllClusters = false;
  while (!mergedAllClusters) {

    // New cluster
    art::PtrVector<recob::Hit> cluster;

    // Put the largest unmerged cluster in this new cluster
    for (unsigned int initCluster = 0; initCluster < planeClusters->size(); ++initCluster) {
      if (planeClusters->at(initCluster).size() < fMinMergeClusterSize or std::find(mergedClusters.begin(), mergedClusters.end(), initCluster) != mergedClusters.end()) continue;
      cluster = planeClusters->at(initCluster);
      mergedClusters.push_back(initCluster);
      break;
    }
    
    // Merge all aligned clusters to this
    bool mergedAllToThisCluster = false;
    while (!mergedAllToThisCluster) {

      // Look at all clusters and merge
      int nadded = 0;
      for (unsigned int trialCluster = 0; trialCluster < planeClusters->size(); ++trialCluster) {

  	if (planeClusters->at(trialCluster).size() < fMinMergeClusterSize or std::find(mergedClusters.begin(), mergedClusters.end(), trialCluster) != mergedClusters.end()) continue;

	// Calculate the PCA for each
	TPrincipal *pca = new TPrincipal(2,"");
	double hits[2];
	TVector2 pos;
	std::pair<double,double> mergedPos(1000,0), trialPos(1000,0);

	for (auto &mergedClusterHits : cluster) {
	  pos = ConvertWireDriftToCm(mergedClusterHits->WireID().Wire, (int)mergedClusterHits->PeakTime(), plane, tpc, cryo);
	  if (pos.Mod() < mergedPos.first) mergedPos.first = pos.Mod();
	  else if (pos.Mod() > mergedPos.second) mergedPos.second = pos.Mod();
	  hits[0] = pos.X();
	  hits[1] = pos.Y();
	  pca->AddRow(hits);
	}
	for (auto &trialClusterHits : planeClusters->at(trialCluster)) {
	  pos = ConvertWireDriftToCm(trialClusterHits->WireID().Wire, (int)trialClusterHits->PeakTime(), plane, tpc, cryo);
	  if (pos.Mod() < trialPos.first) trialPos.first = pos.Mod();
	  else if (pos.Mod() > trialPos.second) trialPos.second = pos.Mod();
	  hits[0] = pos.X();
	  hits[1] = pos.Y();
	  pca->AddRow(hits);
	}

	pca->MakePrincipals();

	// Merge these clusters if they are part of the same straight line
	bool passParallelCut = (*pca->GetEigenValues())[0] > fMergingThreshold;
	bool passProximityCut = (std::abs(trialPos.first - mergedPos.second) < fMaxMergeSeparation) or (std::abs(trialPos.second - mergedPos.first) < fMaxMergeSeparation);

	//std::cout << "Event " << fEvent << ", tpc " << fTPC << " and plane " << fPlane << ". Clusters have eigenvalue of " << (*pca->GetEigenValues())[0] << " and are separated by a distance of around " << std::abs(trialPos.second - mergedPos.first) << " but do they pass the proximity cut? " << passProximityCut << std::endl;

	if (passParallelCut and passProximityCut) {

	  for (auto &hit : planeClusters->at(trialCluster))
	    cluster.push_back(hit);

	  mergedClusters.push_back(trialCluster);
	  ++nadded;

	}

	delete pca;

      } // loop over clusters to add

      if (nadded == 0) mergedAllToThisCluster = true;

    } // while loop

    clusters.push_back(cluster);
    if (mergedClusters.size() == nclusters) mergedAllClusters = true;

  }

  return clusters.size();

}

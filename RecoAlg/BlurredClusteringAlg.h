////////////////////////////////////////////////////////////////////
// Implementation of the Blurred Clustering algorithm
//
// Converts a hit map into a 2D image of the hits before convoling
// with a Gaussian function to introduce a weighted blurring.
// Clustering proceeds on this blurred image to create more
// complete clusters.
//
// M Wallbank (m.wallbank@sheffield.ac.uk), May 2015
////////////////////////////////////////////////////////////////////

#ifndef BlurredClustering_h
#define BlurredClustering_h

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "DataProviders/DetectorProperties.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "MCCheater/BackTracker.h"

// ROOT & C++
#include <TTree.h>
#include <TH2F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TString.h>
#include <TMarker.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLine.h>
#include <TPrincipal.h>
#include <TMath.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TVector2.h>

#include <string>
#include <vector>
#include <map>
#include <sstream>


namespace cluster {
  class BlurredClusteringAlg;
}

class cluster::BlurredClusteringAlg {
public:

  BlurredClusteringAlg(fhicl::ParameterSet const& pset);
  virtual ~BlurredClusteringAlg();

  void reconfigure(fhicl::ParameterSet const&p);
  art::PtrVector<recob::Hit> ConvertBinsToRecobHits(TH2F* image, std::vector<int> const& bins);
  art::Ptr<recob::Hit> ConvertBinToRecobHit(TH2F* image, int const& bin);
  void ConvertBinsToClusters(TH2F *image, std::vector<std::vector<int> > const& allClusterBins, std::vector<art::PtrVector<recob::Hit> >& clusters);
  void CreateDebugPDF(int run, int subrun, int event);
  TH2F ConvertRecobHitsToTH2(std::vector<art::Ptr<recob::Hit> > const& hits);
  TH2F* Convolve(TH2F* image, std::map<int,double> const& kernel, int const& width, int const& height, const char *new_name = 0);
  void FindBlurringParameters(int& blurwire, int& blurtick, int& sigmawire, int& sigmatick);
  int FindClusters(TH2F* image, std::vector<std::vector<int> >& allcluster);
  int FindGlobalWire(geo::WireID const& wireID);
  TH2F* GaussianBlur(TH2F* image);
  unsigned int GetMinSize() { return fMinSize; }
  double GetTimeOfBin(TH2F* image, int const& bin);
  unsigned int NumNeighbours(int const& nx, std::vector<bool> const& used, int const& bin);
  bool PassesTimeCut(std::vector<double> const& times, double const& time);
  void RemoveTrackHits(std::vector<art::Ptr<recob::Hit> > const& ihits, std::vector<art::Ptr<recob::Track> > const& tracks, std::vector<art::Ptr<recob::SpacePoint> > const& spacePoints, art::FindManyP<recob::Track> const& fmth, art::FindManyP<recob::Track> const& fmtsp, art::FindManyP<recob::Hit> const& fmh, std::vector<art::Ptr<recob::Hit> >& hits, int Event, int Run);
  void SaveImage(TH2F* image, std::vector<art::PtrVector<recob::Hit> > const& allClusters, int pad, int tpc, int plane);
  void SaveImage(TH2F* image, int pad, int tpc, int plane);
  void SaveImage(TH2F* image, std::vector<std::vector<int> > const& allClusterBins, int pad, int tpc, int plane);

  std::map<int,std::map<int,art::Ptr<recob::Hit> > > fHitMap;

private:

  // Parameters used in the Blurred Clustering algorithm
  int          fBlurWire;                 // blur radius for Gauss kernel in the wire direction
  int          fBlurTick;                 // blur radius for Gauss kernel in the tick direction
  double       fBlurSigma;                // sigma for Gaussian kernel
  int          fClusterWireDistance;      // how far to cluster from seed in wire direction
  int          fClusterTickDistance;      // how far to cluster from seed in tick direction
  unsigned int fMinMergeClusterSize;      // minimum size of a cluster to consider merging it to another
  double       fMergingThreshold;        // the PCA eigenvalue needed to consider two clusters a merge
  unsigned int fNeighboursThreshold;      // min. number of neighbors to add to cluster
  int          fMinNeighbours;            // minumum number of neighbors to keep in the cluster
  unsigned int fMinSize;                  // minimum size for cluster
  double       fMinSeed;                  // minimum seed after blurring needed before clustering proceeds
  double       fTimeThreshold;            // time threshold for clustering
  double       fChargeThreshold;          // charge threshold for clustering

  // Wire and tick information for histograms
  int fLowerHistTick, fUpperHistTick;
  int fLowerHistWire, fUpperHistWire;

  // Blurring stuff
  int fLastBlurWire;
  int fLastBlurTick;
  double fLastSigma;
  std::map<int,double> fLastKernel;

  // For the debug pdf
  TCanvas *fDebugCanvas;
  std::string fDebugPDFName;

  // art service handles
  art::ServiceHandle<geo::Geometry> fGeom;
  dataprov::DetectorProperties const* fDetProp;

};

#endif

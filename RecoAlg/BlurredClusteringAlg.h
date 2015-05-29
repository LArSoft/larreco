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
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Hit.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/Geometry.h"

// ROOT & C++
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
  art::PtrVector<recob::Hit> ConvertBinsToRecobHits(TH2F *image, std::vector<int> &bins);
  art::Ptr<recob::Hit> ConvertBinToRecobHit(TH2F * image, int bin);
  std::vector<art::PtrVector<recob::Hit> > ConvertBinsToClusters(TH2F *image, std::vector<art::Ptr<recob::Hit> > *allHits, std::vector<std::vector<int> > &allClusterBins);
  void CreateDebugPDF(int fEvent, int fRun, int fSubrun, bool debug);
  TH2F ConvertRecobHitsToTH2(std::vector<art::Ptr<recob::Hit> > *hits);
  TH2* Convolve(TH2 *image, std::map<int,double> kernel, int width, int height, const char *new_name = 0);
  int FindClusters(TH2F *image, std::vector<std::vector<int> > &allcluster);
  TH2* GaussianBlur(TH2 *image);
  unsigned int GetMinSize() { return fMinSize; }
  double GetTimeOfBin(TH2F *image, int bin);
  unsigned int NumNeighbours(int nx, std::vector<bool> *used, int bin);
  bool PassesTimeCut(std::vector<double> &times, double time);
  void SaveImage(TH2F *image, std::vector<art::PtrVector<recob::Hit> > &allClusters, int pad);
  void SaveImage(TH2F *image, int pad);
  void SaveImage(TH2F *image, std::vector<std::vector<int> > &allClusterBins, int pad);

  unsigned int fPlane;
  unsigned int fTPC;
  std::map<int,std::map<int,art::Ptr<recob::Hit> > > fHitMap;

private:

  unsigned int fNWires, fNTicks;
  int fLowerHistTick, fUpperHistTick;
  int fClusterCount; // not being used right now

  // For the debug pdf
  TCanvas *fDebugCanvas;
  std::string fDebugPDFName;
  bool fCreateDebugPDF;

  // Blurring stuff
  int fLastBlurN;
  double fLastSigma;
  std::map<int,double> fLastKernel;

  /// Parameters used in the Blurred Clustering algorithm
  int          fBlurN;
  double       fBlurSigma;             // parameters for gaussian blur
  unsigned int fNeighboursThreshold;   // min. number of neighbors to add to cluster
  int          fMinNeighbours;         // minumum number of neighbors to keep in the cluster
  unsigned int fMinSize;               // minimum size for cluster
  double       fMinSeed;               // minimum seed after blurring needed before clustering proceeds
  double       fTimeThreshold;         // time threshold for clustering
  double       fChargeThreshold;       // charge threshold for clustering

  // Create geometry and detector property handle
  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

};

#endif

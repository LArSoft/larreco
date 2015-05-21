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

#include "RecoAlg/BlurredClusteringAlg.h"

cluster::BlurredClusteringAlg::BlurredClusteringAlg(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset); 

  fClusterCount = 0;
  fNWires = -1;
  fNTicks = -1;
  fLastKernel.clear();
  fLastBlurN = -1000;
  fLastSigma = -1000;
  fMinHits = 3;
  fMinSeed = 0.1;
  fBlurredMin = 0.07;

  // For the debug PDF
  fDebugCanvas = NULL;
  fDebugPDFName = "";

  std::cout << "Beginning blurred" << std::endl;
  mf::LogInfo("Blurred Clustering") << "*** Parameters for Blurred Clustering: ***";
  mf::LogInfo("Blurred Clustering") << "BlurN:               " << fBlurN;
  mf::LogInfo("Blurred Clustering") << "BlurSigma:           " << fBlurSigma;
  mf::LogInfo("Blurred Clustering") << "NeighboursThreshold: " << fNeighboursThreshold;
  mf::LogInfo("Blurred Clustering") << "MinNeighbours:       " << fMinNeighbours;
  mf::LogInfo("Blurred Clustering") << "MinSize:             " << fMinSize;
  mf::LogInfo("Blurred Clustering") << "TimeThreshold:       " << fTimeThreshold;
}

cluster::BlurredClusteringAlg::~BlurredClusteringAlg() {
  if (fDebugCanvas) {
    std::string closeName = fDebugPDFName;
    closeName.append("]");
    fDebugCanvas->Print(closeName.c_str());
    delete fDebugCanvas;
  }
  if (!fLastKernel.empty())
    fLastKernel.clear();
}

void cluster::BlurredClusteringAlg::reconfigure(fhicl::ParameterSet const& p) {
  fBlurN               = p.get<int>   ("BlurN");
  fBlurSigma           = p.get<double>("BlurSigma");
  fNeighboursThreshold = p.get<int>   ("NeighboursThreshold");
  fMinNeighbours       = p.get<int>   ("MinNeighbours");
  fMinSize             = p.get<int>   ("MinSize");
  fTimeThreshold       = p.get<double>("TimeThreshold");
}


// Create the PDF to save debug images
void cluster::BlurredClusteringAlg::CreateDebugPDF(int event, int run, int subrun, bool debug) {

  fCreateDebugPDF = debug;

  if (!fDebugCanvas) {

    // Create the grayscale palette for the Z axis
    Double_t Red[2] = { 1.00, 0.00 };
    Double_t Green[2] = { 1.00, 0.00 };
    Double_t Blue[2] = { 1.00, 0.00 };
    Double_t Length[2] = { 0.00, 1.00 };
    TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, 1000);
    gStyle->SetOptStat(110000);

    // Decide what to call this PDF
    std::ostringstream oss;
    oss << "BlurredImages_Run" << run << "_Subrun" << subrun;
    fDebugPDFName = oss.str();
    fDebugCanvas = new TCanvas(fDebugPDFName.c_str(), "Image canvas", 1000, 500);
    fDebugPDFName.append(".pdf");

    std::string openName = fDebugPDFName;
    openName.append("[");
    fDebugCanvas->Print(openName.c_str());
    fDebugCanvas->Divide(2, 2);
    fDebugCanvas->SetGrid();
  }

  // Clear the pads on the canvas
  for (int i = 1; i <= 4; i++) {
    fDebugCanvas->GetPad(i)->Clear();
  }

  std::ostringstream oss;
  oss << "Event " << event;
  fDebugCanvas->cd(1);
  TLatex l;
  l.SetTextSize(0.15);
  l.DrawLatex(0.1, 0.1, oss.str().c_str());
  fDebugCanvas->Print(fDebugPDFName.c_str());

}


// Converts a vector of bins into a hit selection - not all the hits in the bins vector are real hits
art::PtrVector<recob::Hit> cluster::BlurredClusteringAlg::ConvertBinsToRecobHits(TH2F *image, std::vector<int> &bins) {

  // Create the vector of hits to output
  art::PtrVector<recob::Hit> hits;

  // Look through the hits in the cluster
  for (std::vector<int>::iterator it = bins.begin(); it != bins.end(); it++) {

    // Take each hit and convert it to a recob::Hit
    int bin = *it;
    art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);

    // If this hit was a real hit put it in the hit selection
    if (!hit.isNull())  /// CHECK THIS!!!!!!!!
      hits.push_back(hit);
  }

  // Return the vector of hits to make cluster
  return hits;
}


// Converts a bin into a recob::Hit (not all of these bins correspond to recob::Hits - some are fake hits created by the blurring)
art::Ptr<recob::Hit> cluster::BlurredClusteringAlg::ConvertBinToRecobHit(TH2F *image, int bin) {

  // Take the hit
  art::Ptr<recob::Hit> hit;

  // Find the point in wire/tick space which corresponds to this bin
  int wire, tick, z;
  image->GetBinXYZ(bin,wire,tick,z);

  // If this particular bin corresponds to an entry in the hit map then this was a real hit
  if (fHitMap.find(wire-1) != fHitMap.end()) {
    if (fHitMap[wire-1].find(tick+fLowerHistTick-1) != fHitMap[wire-1].end()) {
      hit = fHitMap[wire-1][tick+fLowerHistTick-1];
    }
  }

  // Return this hit place in art::PtrVector<recob::Hit>
  return hit;
}


// Takes a vector of clusters (itself a vector of hits) and turns them into clusters using the initial hit selection (hits)
std::vector<art::PtrVector<recob::Hit> > cluster::BlurredClusteringAlg::ConvertBinsToClusters(TH2F *image, std::vector<art::Ptr<recob::Hit> > *hits, std::vector<std::vector<int> > &allClusterBins) {

  // Output is a vector of clusters
  std::vector<art::PtrVector<recob::Hit> > clusters;

  // Loop through the clusters (each a vector of bins)
  for (std::vector<std::vector<int> >::iterator clustIt = allClusterBins.begin(); clustIt != allClusterBins.end(); clustIt++) {
    std::vector<int> bins = *clustIt;

    // Convert the clusters (vectors of bins) to hits in a vector of recob::Hits
    art::PtrVector<recob::Hit> clusHits = ConvertBinsToRecobHits(image, bins);

    mf::LogInfo("BlurredClustering") << "Cluster made from " << bins.size() << " bins, of which " << clusHits.size() << " were real hits";

    clusters.push_back(clusHits);
  }

  SaveImage(image, clusters, 4);

  return clusters;
}


/// Takes hit map and returns a TH2 histogram of bar vs layer, with charge on z-axis
TH2F cluster::BlurredClusteringAlg::ConvertRecobHitsToTH2(std::vector<art::Ptr<recob::Hit> >* hits) {

  // Define the size of this particular plane
  art::Ptr<recob::Hit> firstHit = *(hits->begin());
  fNWires = fGeom->Nwires(firstHit->WireID().Plane,firstHit->WireID().TPC,firstHit->WireID().Cryostat);
  //fNTicks = fDetProp->ReadOutWindowSize();
  // Try looping over all hits so the histogram isn't so huge!
  int lowerTick = fDetProp->ReadOutWindowSize(), upperTick = 0;
  for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits->begin(); hitIt != hits->end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    if (hit->PeakTime() < lowerTick) lowerTick = hit->PeakTime();
    if (hit->PeakTime() > upperTick) upperTick = hit->PeakTime();
  }
  fLowerHistTick = lowerTick-20;
  fUpperHistTick = upperTick+20;

  // Use a map to keep a track of the real hits and their wire/ticks
  fHitMap.clear();

  std::stringstream planeImage;
  planeImage << "blurred_plane" << firstHit->WireID().Plane << "_image";

  // Create a TH2 histogram
  //TH2F image(planeImage.str().c_str(), planeImage.str().c_str(), fNWires, -0.5, fNWires - 0.5, fNTicks, -0.5, fNTicks - 0.5);
  TH2F image(planeImage.str().c_str(), planeImage.str().c_str(), fNWires, -0.5, fNWires - 0.5, (fUpperHistTick-fLowerHistTick), fLowerHistTick-0.5, fUpperHistTick-0.5);
  image.Clear();
  image.SetXTitle("Wire number");
  image.SetYTitle("Tick number");
  image.SetZTitle("Charge");

  // Look through the hits
  for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits->begin(); hitIt != hits->end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    unsigned int wire   = hit->WireID().Wire;
    int          tick   = (int)hit->PeakTime();
    float        charge = hit->SummedADC();

    // Fill hit map and keep a note of all real hits for later
    if (charge > image.GetBinContent(wire+1,(tick-fLowerHistTick)+1)) {
      image.Fill(wire,tick,charge);
      fHitMap[wire][tick] = hit;
    }
  }

  return image;
}


// Convolves the Gaussian kernel with the image to blurrr
TH2 *cluster::BlurredClusteringAlg::Convolve(TH2* image, std::map<int,double> kernel, int width, int height, const char *new_name) {

  // Get the magnitude of the bins in the kernel
  double mag = 0;
  for(std::map<int,double>::iterator it = kernel.begin(); it != kernel.end(); ++it)
    mag += it->second;

  // Copy the histogram to blur
  TString name = new_name == 0 ? TString::Format("%s_convolved", image->GetName()) : TString(new_name);
  TH2 *copy = (TH2*)image->Clone(name);

  int nbinsx = copy->GetNbinsX();
  int nbinsy = copy->GetNbinsY();

  // Loop through all the bins in the histogram to blur
  for (int x = 1; x <= nbinsx; ++x) {
    for (int y = 1; y <= nbinsy; ++y) {

      // Create new variables for each bin
      double newval = 0;
      double norml = 0;

      // Loop over a blurring region (nothing to do with current bin - will be used as a distance from current bin to blur later on)
      for (int blurx = -width / 2; blurx < (width + 1) / 2; blurx++) {
	for (int blury = -height / 2; blury < (height + 1) / 2; blury++) {

	  // Form a weight for each bin based on the Gaussian kernel
	  // NB/ This is simply a measure of how large the Gaussian weight is this far from the seed - nothing to do with current bins
	  double weight = kernel[width * (height / 2 + blury) + (width / 2 + blurx)];

	  // Look at the blurred region - if still within the hit map, increase the weight according to the value of the kernel at this blurring distance
	  if (x + blurx >= 1 && x + blurx <= nbinsx && y + blury >= 1 && y + blury <= nbinsy) {
	    newval += weight * image->GetBinContent(x + blurx, y + blury);
	    norml += weight;
	  }

	}
      } // End loop over blurring region

      // Set the new value for this bin in the blurred copy of the image
      newval /= (norml / mag);
      copy->SetBinContent(x, y, newval);

    }
  } // End loop over bins in histogram

  // Return the blurred histogram
  return copy; 
}


/// Find clusters in the histogram
int cluster::BlurredClusteringAlg::FindClusters(TH2F *image, std::vector<std::vector<int> > &allcluster) {

  // Save images to debug pdf if fCreateDebugPDF is true
  SaveImage(image, 1);

  // Blur the image using a Gaussian blur
  TH2F *blurred = (TH2F*) GaussianBlur(image);

  SaveImage(blurred, 2);

  // Vectors to hold cluster information
  std::vector<int> cluster;
  std::vector<double> times;

  // Size of image in x and y
  const int nbinsx = image->GetNbinsX() + 2;
  const int nbinsy = image->GetNbinsY() + 2;
  const int nbins = nbinsx * nbinsy;

  // Vectors to hold hit information
  std::vector<bool> used(nbins);
  std::vector<std::pair<double, int> > values(nbins);

  // Place the bin number and contents as a pair in the values vector
  for (int xbin = 1; xbin <= nbinsx; xbin++) {
    for (int ybin = 1; ybin <= nbinsy; ybin++) {
      int bin = image->GetBin(xbin, ybin);
      values.push_back(std::pair<double, int>(blurred->GetArray()[bin], bin));
    }
  }

  // Sort the values into charge order
  std::sort(values.rbegin(), values.rend());

  // Count the number of iterations of the cluster forming loop (== number of clusters)
  int niter = 0;


  // Clustering loops
  // First loop - considers highest charge hits in decreasing order, and puts them in a new cluster if they aren't already clustered (makes new cluster every iteration)
  // Second loop - looks at the direct neighbours of this seed and clusters to this if above charge/time thresholds. Runs recursively over all hits in cluster (inc. new ones)
  while (true) {

    // Start a new cluster each time loop is executed
    cluster.clear();
    times.clear();

    // Get the highest charge bin (go no further if below seed threshold) [currently 0.1]
    double blurred_binval = values[niter].first;
    if (blurred_binval < fMinSeed)
      break;

    // Iterate through the bins from highest charge down
    int bin = values[niter++].second;

    // Put this bin in used if not already there
    if (used[bin])
      continue;
    used[bin] = true;

    // Start a new cluster
    cluster.push_back(bin);

    // Get the time of this hit
    double time = GetTimeOfBin(image, bin);
    if (time > 0)
      times.push_back(time);


    // Now cluster neighbouring hits to this seed
    while (true) {

      int nadded = 0;

      for (unsigned int clusBin = 0; clusBin < cluster.size(); ++clusBin) {

	// Get x and y values for bin (c++ returns a%b = a if a<b)
        int binx, biny;
        binx = cluster[clusBin] % nbinsx;
        biny = ((cluster[clusBin] - binx) / nbinsx) % nbinsy;

	// Look for hits in the directly neighbouring x/y bins
        for (int x = binx - 1; x <= binx + 1; x++) {
          for (int y = biny - 1; y <= biny + 1; y++) {
            if (x == binx && y == biny)
              continue;

	    // Get this bin
            bin = image->GetBin(x, y);
            if (used[bin])
              continue;

	    // Get the blurred value and time for this bin
            blurred_binval = blurred->GetArray()[bin];
            time = GetTimeOfBin(image, bin); /// NB for 'fake' hits, time is defaulted to -10000

	    // Check real hits pass time cut (ignores fake hits)
            if (time > 0 && times.size() > 0 && ! PassesTimeCut(times, time))
	      continue;

	    // Add to cluster if bin value is above threshold
            if (blurred_binval > fBlurredMin) {
              used[bin] = true;
              cluster.push_back(bin);
              nadded++;

              if (time > 0) {
                times.push_back(time);
              }

            } // End of adding blurred bin to cluster

          }
       } // End of looking at directly neighbouring bins

      } // End of looping over bins already in this cluster

      if (nadded == 0)
        break;

    } // End of adding hits to this cluster

    // Check this cluster is above minimum size
    if (cluster.size() < fMinSize) {
      for (unsigned int i = 0; i < cluster.size(); i++)
        used[cluster[i]] = false;
      continue;
    }

    // Fill in holes in the cluster
    for (unsigned int clusBin = 0; clusBin < cluster.size(); clusBin++) {

      // Looks at directly neighbouring bins (and not itself)
      for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
          if (!x && !y) continue;

	  // Look at neighbouring bins to the clustered bin which are inside the cluster
          int neighbouringBin = cluster[clusBin] + x + (y * nbinsx);
          if (neighbouringBin < nbinsx || neighbouringBin % nbinsx == 0 || neighbouringBin % nbinsx == nbinsx - 1 || neighbouringBin >= nbinsx * (nbinsy - 1))
            continue;

          double time = GetTimeOfBin(image, neighbouringBin);

	  // If not already clustered and passes neighbour/time thresholds, add to cluster
          if ( !used[neighbouringBin] && (NumNeighbours(nbinsx, &used, neighbouringBin) > fNeighboursThreshold) && PassesTimeCut(times, time) ) {
            used[neighbouringBin] = true;
            cluster.push_back(neighbouringBin);

            if (time > 0) {
              times.push_back(time);
            }
          } // End of clustering neighbouring bin

        }
      } // End of looping over neighbouring bins

    } // End of looping over bins already in cluster

    mf::LogVerbatim("Blurred Clustering") << "Size of cluster after filling in holes: " << cluster.size();


    // Remove peninsulas - hits which have too few neighbouring hits in the cluster (defined by fMinNeighbours)
    while (true) {
      int nremoved = 0;

      // Loop over all the bins in the cluster
      for (int clusBin = cluster.size() - 1; clusBin >= 0; clusBin--) {
        bin = cluster[clusBin];

	// If bin is in cluster ignore
        if (bin < nbinsx || bin % nbinsx == 0 || bin % nbinsx == nbinsx - 1 || bin >= nbinsx * (nbinsy - 1)) continue;

	// Remove hit if it has too few neighbouring hits
        if ((int) NumNeighbours(nbinsx, &used, bin) < fMinNeighbours) {
          used[bin] = false;
          nremoved++;
          cluster.erase(cluster.begin() + clusBin);
          blurred_binval = blurred->GetArray()[bin];
        }
      }

      if (!nremoved)
        break;
    }

    mf::LogVerbatim("Blurred Clustering") << "Size of cluster after removing peninsulas: " << cluster.size();


    // Disregard cluster if not of minimum size
    if (cluster.size() < fMinSize) {
      for (unsigned int i = 0; i < cluster.size(); i++)
        used[cluster[i]] = false;
      continue;
    }

    // Put this cluster in the vector of clusters
    allcluster.push_back(cluster);

  } // End loop over this cluster


  SaveImage(blurred, allcluster, 3);

  blurred->Delete();

  /// Return the number of clusters found in this hit map
  return allcluster.size();
}


// Applies Gaussian blur to image
TH2 *cluster::BlurredClusteringAlg::GaussianBlur(TH2 *image) {

  // Create Gaussian kernel
  std::map<int,double> kernel;
  int width = 2 * fBlurN + 1;
  int height = 2 * fBlurN + 1;

  if (fBlurSigma == 0)
    return (TH2F*) image->Clone(image->GetName() + TString("_blur"));

  // If the parameters match the last parameters, used the same last kernel
  if (fLastBlurN == fBlurN && fLastSigma == fBlurSigma && !fLastKernel.empty())
    kernel = fLastKernel;

  // Otherwise, allocate and compute a new kernel, freeing the fLastKernel if it has been allocated.
  else {

    if (!fLastKernel.empty())
      fLastKernel.clear();

    // Smear out according to the blur radius fBlurN
    for (int i = -fBlurN; i <= fBlurN; i++) {
      for (int j = -fBlurN; j <= fBlurN; j++) {
        double sigmai = fBlurSigma / 1.5;
        double sigmaj = fBlurSigma;

	// Fill kernel
	double sig2i = 2. * sigmai * sigmai;
	double sig2j = 2. * sigmaj * sigmaj;

	int key = (width * (j + fBlurN)) + (i + fBlurN);
	double value = 1. / sqrt(sig2i * M_PI) * exp(-i * i / sig2i) * 1. / sqrt(sig2j * M_PI) * exp(-j * j / sig2j);
	kernel[key] = value;

      }
    } // End loop over blurring region

    fLastKernel = kernel;
    fLastBlurN  = fBlurN;
    fLastSigma  = fBlurSigma;
  }

  // Return a convolution of this Gaussian kernel with the unblurred hit map
  return (TH2F*) Convolve(image, kernel, width, height, TString::Format("%s_blur", image->GetName()));

}


// Returns the hit time of a hit in a particular bin
double cluster::BlurredClusteringAlg::GetTimeOfBin(TH2F *image, int bin) {

  double time = -10000;

  art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);
  if (!hit.isNull()) // CHECK THIS!!!!!!!!!!
    time = hit->PeakTime();

  return time;

}


// Determines the number of clustered neighbours of a hit
unsigned int cluster::BlurredClusteringAlg::NumNeighbours(int nbinsx, std::vector<bool> *used, int bin) {
  unsigned int neighbours = 0; /// Everybody wants good neighbours

  // Loop over all directly neighbouring hits (not itself)
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      if (!x && !y) continue;

      // Determine bin
      int neighbouringBin = bin + x + (y * nbinsx); /// 2D hists can be considered a string of bins - the equation to convert between them is [bin = x + (nbinsx * y)]

      // If this bin is in the cluster, increase the neighbouring bin counter
      if (used->at(neighbouringBin))
	neighbours++;
    }
  }

  // Return the number of neighbours in the cluster of a particular hit
  return neighbours;
}


// Determine if a hit is within a time threshold of any other hits in a cluster
bool cluster::BlurredClusteringAlg::PassesTimeCut(std::vector<double> &times, double time) {

  for (std::vector<double>::iterator timeIt = times.begin(); timeIt != times.end(); timeIt++) {
    if (std::abs(time - *timeIt) < fTimeThreshold) return true;
  }

  return false;
}


// Save the images for debugging
void cluster::BlurredClusteringAlg::SaveImage(TH2F *image, std::vector<art::PtrVector<recob::Hit> > &allClusters, int pad) {
  if (!fCreateDebugPDF) {
    return;
  }

  /// Make a vector of clusters
  std::vector<std::vector<int> > allClusterBins;

  for (std::vector<art::PtrVector<recob::Hit> >::iterator clusterIt = allClusters.begin(); clusterIt != allClusters.end(); clusterIt++) {
    art::PtrVector<recob::Hit> cluster = *clusterIt;

    if (!cluster.size())
      continue;

    std::vector<int> clusterBins;

    for (art::PtrVector<recob::Hit>::iterator hitIt = cluster.begin(); hitIt != cluster.end(); hitIt++) {
      art::Ptr<recob::Hit> hit = *hitIt;
      unsigned int wire = hit->WireID().Wire;
      float tick = hit->PeakTime();
      int bin = image->GetBin(wire+1,(tick-fLowerHistTick)+1);
      if (cluster.size() < fMinHits)
        bin *= -1;

      clusterBins.push_back(bin);
    }

    allClusterBins.push_back(clusterBins);
  }

  SaveImage(image, allClusterBins, pad);
}

void cluster::BlurredClusteringAlg::SaveImage(TH2F *image, int pad) {
  std::vector<std::vector<int> > allClusterBins;
  SaveImage(image, allClusterBins, pad);
}

void cluster::BlurredClusteringAlg::SaveImage(TH2F *image, std::vector<std::vector<int> > &allClusterBins, int pad) {
  if (!fCreateDebugPDF) {
    return;
  }

  /// Save the different stages on different pads
  fDebugCanvas->cd(pad);
  std::string stage;

  switch (pad) {
    case 1:
      stage = "Stage 1: Unblurred";
      break;
    case 2:
      stage = "Stage 2: Blurred";
      break;
    case 3:
      stage = "Stage 3: Blurred with clusters overlaid";
      break;
    case 4:
      stage = "Stage 4: Output clusters";
      break;
    default:
      stage = "Unknown stage";
      break;
  }

  std::stringstream title;
  title << stage << " - TPC " << fTPC << ", Plane " << fPlane;

  image->SetName(title.str().c_str());
  image->SetTitle(title.str().c_str());
  image->DrawClone("COLZ");

  int clusterNum = 2;

  for (std::vector<std::vector<int> >::iterator it = allClusterBins.begin(); it != allClusterBins.end(); it++, clusterNum++) {
    std::vector<int> bins = *it;
    TMarker mark(0, 0, 20);
    mark.SetMarkerColor(clusterNum);
    mark.SetMarkerSize(0.1);

    for (std::vector<int>::iterator binIt = bins.begin(); binIt != bins.end(); binIt++) {
      int bin = *binIt;
      int wire, tick, z;

      /// Hit from a cluster that we aren't going to save
      if (bin < 0) {
        bin *= -1;
        mark.SetMarkerStyle(24);
      }

      image->GetBinXYZ(bin,wire,tick,z);
      mark.DrawMarker(wire-1, tick+fLowerHistTick-1);
      mark.SetMarkerStyle(20);
    }
  }

  if (pad == 4)
    fDebugCanvas->Print(fDebugPDFName.c_str());
}

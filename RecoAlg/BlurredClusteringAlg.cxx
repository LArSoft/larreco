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

  fLastKernel.clear();
  fLastBlurWire = -1000;
  fLastBlurTick = -1000;
  fLastSigma = -1000;

  // For the debug PDF
  fDebugCanvas = NULL;
  fDebugPDFName = "";

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
  fBlurWire            = p.get<int>   ("BlurWire");
  fBlurTick            = p.get<int>   ("BlurTick");
  fBlurSigma           = p.get<double>("BlurSigma");
  fClusterWireDistance = p.get<int>   ("ClusterWireDistance");
  fClusterTickDistance = p.get<int>   ("ClusterTickDistance");
  fNeighboursThreshold = p.get<int>   ("NeighboursThreshold");
  fMinNeighbours       = p.get<int>   ("MinNeighbours");
  fMinSize             = p.get<int>   ("MinSize");
  fMinSeed             = p.get<double>("MinSeed");
  fTimeThreshold       = p.get<double>("TimeThreshold");
  fChargeThreshold     = p.get<double>("ChargeThreshold");
}


void cluster::BlurredClusteringAlg::CreateDebugPDF(int run, int subrun, int event) {

  /// Create the PDF to save debug images

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


art::PtrVector<recob::Hit> cluster::BlurredClusteringAlg::ConvertBinsToRecobHits(std::vector<std::vector<double> > const& image, std::vector<int> const& bins) {

  /// Converts a vector of bins into a hit selection - not all the hits in the bins vector are real hits

  // Create the vector of hits to output
  art::PtrVector<recob::Hit> hits;

  // Look through the hits in the cluster
  for (std::vector<int>::const_iterator binIt = bins.begin(); binIt != bins.end(); binIt++) {

    // Take each hit and convert it to a recob::Hit
    art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, *binIt);

    // If this hit was a real hit put it in the hit selection
    if (!hit.isNull())
      hits.push_back(hit);
  }

  // Return the vector of hits to make cluster
  return hits;
}


art::Ptr<recob::Hit> cluster::BlurredClusteringAlg::ConvertBinToRecobHit(std::vector<std::vector<double> > const& image, int bin) {

  /// Converts a bin into a recob::Hit (not all of these bins correspond to recob::Hits - some are fake hits created by the blurring)

  // Take the hit
  art::Ptr<recob::Hit> hit;

  // Find the point in wire/tick space which corresponds to this bin
  int xbin = bin % image.size();
  int ybin = bin / image.size();
  int wire = xbin + fLowerHistWire;
  int tick = ybin + fLowerHistWire;

  // If this particular bin corresponds to an entry in the hit map then this was a real hit
  if (fHitMap.find(wire) != fHitMap.end()) {
    if (fHitMap[wire].find(tick) != fHitMap[wire].end()) {
      hit = fHitMap[wire][tick];
    }
  }

  // Return this hit place in art::PtrVector<recob::Hit>
  return hit;

}

//   // Find the point in wire/tick space which corresponds to this bin
//   int wire, tick, z;
//   image->GetBinXYZ(bin,wire,tick,z);

//   // If this particular bin corresponds to an entry in the hit map then this was a real hit
//   if (fHitMap.find(wire+fLowerHistWire-1) != fHitMap.end()) {
//     if (fHitMap[wire+fLowerHistWire-1].find(tick+fLowerHistTick-1) != fHitMap[wire+fLowerHistWire-1].end()) {
//       hit = fHitMap[wire+fLowerHistWire-1][tick+fLowerHistTick-1];
//     }
//   }

//   // Return this hit place in art::PtrVector<recob::Hit>
//   return hit;
// }

void cluster::BlurredClusteringAlg::ConvertBinsToClusters(std::vector<std::vector<double> > const& image,
							  std::vector<std::vector<int> > const& allClusterBins,
							  std::vector<art::PtrVector<recob::Hit> >& clusters) {

  /// Takes a vector of clusters (itself a vector of hits) and turns them into clusters using the initial hit selection

  // Loop through the clusters (each a vector of bins)
  for (std::vector<std::vector<int> >::const_iterator clustIt = allClusterBins.begin(); clustIt != allClusterBins.end(); clustIt++) {
    std::vector<int> bins = *clustIt;

    // Convert the clusters (vectors of bins) to hits in a vector of recob::Hits
    art::PtrVector<recob::Hit> clusHits = ConvertBinsToRecobHits(image, bins);

    mf::LogInfo("BlurredClustering") << "Cluster made from " << bins.size() << " bins, of which " << clusHits.size() << " were real hits";

    // Make sure the clusters are above the minimum cluster size
    if (clusHits.size() < fMinSize) {
      mf::LogVerbatim("BlurredClustering") << "Cluster of size " << clusHits.size() << " not saved since it is smaller than the minimum cluster size, set to " << fMinSize;
      continue;
    }

    clusters.push_back(clusHits);
  }

  return;

}


std::vector<std::vector<double> > cluster::BlurredClusteringAlg::ConvertRecobHitsToTH2(std::vector<art::Ptr<recob::Hit> > const& hits) {

  /// Takes hit map and returns a TH2 histogram of bar vs layer, with charge on z-axis

  // Define the size of this particular plane -- dynamically to avoid huge histograms
  int lowerTick = fDetProp->ReadOutWindowSize(), upperTick = 0, lowerWire = fGeom->MaxWires(), upperWire = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    int histWire = GlobalWire((*hitIt)->WireID());
    if ((*hitIt)->PeakTime() < lowerTick) lowerTick = (*hitIt)->PeakTime();
    if ((*hitIt)->PeakTime() > upperTick) upperTick = (*hitIt)->PeakTime();
    if (histWire < lowerWire) lowerWire = histWire;
    if (histWire > upperWire) upperWire = histWire;
  }
  fLowerHistTick = lowerTick-20;
  fUpperHistTick = upperTick+20;
  fLowerHistWire = lowerWire-20;
  fUpperHistWire = upperWire+20;

  // Use a map to keep a track of the real hits and their wire/ticks
  fHitMap.clear();

  std::stringstream planeImage;
  planeImage << "blurred_image";

  // // Create a TH2 histogram
  // TH2F image(planeImage.str().c_str(), planeImage.str().c_str(), (fUpperHistWire-fLowerHistWire), fLowerHistWire-0.5, fUpperHistWire-0.5, (fUpperHistTick-fLowerHistTick), fLowerHistTick-0.5, fUpperHistTick-0.5);
  // image.Clear();
  // image.SetXTitle("Wire number");
  // image.SetYTitle("Tick number");
  // image.SetZTitle("Charge");

  // Create a 2D vector
  std::vector<std::vector<double> > image(fUpperHistWire-fLowerHistWire, std::vector<double>(fUpperHistTick-fLowerHistTick, 0));

  // Look through the hits
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    unsigned int wire = GlobalWire((*hitIt)->WireID());
    int   tick   = (int)(*hitIt)->PeakTime();
    float charge = (*hitIt)->SummedADC();

    // Fill hit map and keep a note of all real hits for later
    if (charge > image.at(wire-fLowerHistWire).at(tick-fLowerHistTick)) {
      //image.Fill(wire,tick,charge);
      image.at(wire-fLowerHistWire).at(tick-fLowerHistTick) = charge;
      fHitMap[wire][tick] = (*hitIt);
    }
  }

  return image;

}


int cluster::BlurredClusteringAlg::ConvertWireTickToBin(std::vector<std::vector<double> > const& image, int xbin, int ybin) {

  /// Converts an xbin and a ybin to a global bin number                                                                                                                        
  return (ybin * image.size()) + xbin;

}


int cluster::BlurredClusteringAlg::ConvertBinToCharge(std::vector<std::vector<double> > const& image, int bin) {

  /// Returns the charge stored in the global bin value                                                                                                                         
  int x = bin % image.size();
  int y = bin / image.size();

  return image.at(x).at(y);

}


std::vector<std::vector<double> > cluster::BlurredClusteringAlg::Convolve(std::vector<std::vector<double> > const& image,
									  std::vector<double> const& kernel, 
									  int width, int height) {//, const char *new_name) {

  /// Convolves the Gaussian kernel with the image to blurrr

  // Get the magnitude of the bins in the kernel
  double mag = 0;
  for(std::vector<double>::const_iterator kernelIt = kernel.begin(); kernelIt != kernel.end(); ++kernelIt)
    mag += *kernelIt;

  // Copy the histogram to blur
  // TString name = new_name == 0 ? TString::Format("%s_convolved", image->GetName()) : TString(new_name);
  // TH2F *copy = (TH2F*)image->Clone(name);
  std::vector<std::vector<double> > copy = image;

  int nbinsx = copy.size();
  int nbinsy = copy.at(0).size();

  // Loop through all the bins in the histogram to blur
  for (int x = 0; x < nbinsx; ++x) {
    for (int y = 0; y < nbinsy; ++y) {

      // Create new variables for each bin
      double newval = 0;
      double norml = 0;

      // Loop over a blurring region (nothing to do with current bin - will be used as a distance from current bin to blur later on)
      for (int blurx = -width / 2; blurx < (width + 1) / 2; ++blurx) {
	for (int blury = -height / 2; blury < (height + 1) / 2; ++blury) {

	  // Form a weight for each bin based on the Gaussian kernel
	  // NB/ This is simply a measure of how large the Gaussian weight is this far from the seed - nothing to do with current bins
	  double weight = kernel.at(width * (height / 2 + blury) + (width / 2 + blurx));

	  // Look at the blurred region - if still within the hit map, increase the weight according to the value of the kernel at this blurring distance
	  // if (x + blurx >= 1 and x + blurx <= nbinsx and y + blury >= 1 and y + blury <= nbinsy) {
	  //   newval += weight * image->GetBinContent(x + blurx, y + blury);
	  //   norml += weight;
	  // }
	  if (x + blurx >= 0 and x + blurx < nbinsx and y + blury >= 0 and y + blury < nbinsy) {
	    newval += weight * image.at(x+blurx).at(y+blury);
	    norml += weight;
	  }

	}
      } // End loop over blurring region

      // Set the new value for this bin in the blurred copy of the image
      newval /= (norml / mag);
      //copy->SetBinContent(x, y, newval);
      copy.at(x).at(y) = newval;

    }
  } // End loop over bins in histogram

  // Return the blurred histogram
  return copy; 

}


void cluster::BlurredClusteringAlg::FindBlurringParameters(int& blurwire, int& blurtick, int& sigmawire, int& sigmatick) {

  /// Dynamically find the blurring radii in each direction

  // Calculate least squares slope
  int x, y;
  double nhits=0, sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (const auto &wireIt : fHitMap) {
    for (const auto &tickIt : wireIt.second) {
      ++nhits;
      x = wireIt.first;
      y = tickIt.first;
      sumx += x;
      sumy += y;
      sumx2 += x*x;
      sumxy += x*y;
    }
  }
  double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);

  // Get the rough unit vector for the trajectories
  TVector2 unit = TVector2(1,gradient).Unit();

  // Use this direction to scale the blurring radii and Gaussian sigma
  blurwire = std::max(std::abs(std::round(fBlurWire * unit.X())),1.);
  blurtick = std::max(std::abs(std::round(fBlurTick * unit.Y())),1.);

  sigmawire = std::max(std::abs(std::round(fBlurSigma * unit.X())),1.);
  sigmatick = std::max(std::abs(std::round(fBlurSigma * unit.Y())),1.);

  // std::cout << "Blurring: wire " << blurwire << " and tick " << blurtick << "; sigma: wire " << sigmawire << " and tick " << sigmatick << std::endl;

  return;

}


int cluster::BlurredClusteringAlg::FindClusters(std::vector<std::vector<double> > const& blurred, std::vector<std::vector<int> >& allcluster) {

  /// Find clusters in the histogram

  // Vectors to hold cluster information
  std::vector<int> cluster;
  std::vector<double> times;

  // Size of image in x and y
  // const int nbinsx = blurred->GetNbinsX() + 2;
  // const int nbinsy = blurred->GetNbinsY() + 2;
  const int nbinsx = blurred.size();// + 2;
  const int nbinsy = blurred.at(0).size();// + 2;
  const int nbins = nbinsx * nbinsy;

  // Vectors to hold hit information
  std::vector<bool> used(nbins);
  std::vector<std::pair<double, int> > values(nbins);

  // Place the bin number and contents as a pair in the values vector
  for (int xbin = 0; xbin < nbinsx; ++xbin) {
    for (int ybin = 0; ybin < nbinsy; ++ybin) {
      //int bin = blurred->GetBin(xbin, ybin);
      //values.push_back(std::pair<double, int>(blurred->GetArray()[bin], bin));
      int bin = ConvertWireTickToBin(blurred, xbin, ybin);
      values.push_back(std::make_pair(ConvertBinToCharge(blurred, bin), bin));
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

    // Get the highest charge bin (go no further if below seed threshold)
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
    double time = GetTimeOfBin(blurred, bin);
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

	// Look for hits in the neighbouring x/y bins
        for (int x = binx - fClusterWireDistance; x <= binx + fClusterWireDistance; x++) {
          for (int y = biny - fClusterTickDistance; y <= biny + fClusterTickDistance; y++) {
            if (x == binx && y == biny)
              continue;

	    // Get this bin
            //bin = blurred->GetBin(x, y);
	    bin = ConvertWireTickToBin(blurred, x, y); // NOT SURE!
            if (used[bin])
              continue;

	    // Get the blurred value and time for this bin
            //blurred_binval = blurred->GetArray()[bin];
	    blurred_binval = ConvertBinToCharge(blurred, bin);
            time = GetTimeOfBin(blurred, bin); // NB for 'fake' hits, time is defaulted to -10000

	    // Check real hits pass time cut (ignores fake hits)
            if (time > 0 && times.size() > 0 && ! PassesTimeCut(times, time))
	      continue;

	    // Add to cluster if bin value is above threshold
            if (blurred_binval > fChargeThreshold) {
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

          double time = GetTimeOfBin(blurred, neighbouringBin);

	  // If not already clustered and passes neighbour/time thresholds, add to cluster
          if ( !used[neighbouringBin] && (NumNeighbours(nbinsx, used, neighbouringBin) > fNeighboursThreshold) && PassesTimeCut(times, time) ) {
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
        if ((int) NumNeighbours(nbinsx, used, bin) < fMinNeighbours) {
          used[bin] = false;
          nremoved++;
          cluster.erase(cluster.begin() + clusBin);
          //blurred_binval = blurred->GetArray()[bin];
	  blurred_binval = ConvertBinToCharge(blurred, bin);
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

  // Return the number of clusters found in this hit map
  return allcluster.size();

}


int cluster::BlurredClusteringAlg::GlobalWire(geo::WireID const& wireID) {

  /// Find the global wire position

  double wireCentre[3];
  fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);

  double globalWire = -999;
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }
  else {
    unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
    if (wireID.TPC == 0 or wireID.TPC == 1) globalWire = wireID.Wire;
    else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5) globalWire = nwires + wireID.Wire;
    else if (wireID.TPC == 6 or wireID.TPC == 7) globalWire = (2*nwires) + wireID.Wire;
    else mf::LogError("BlurredClusterAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC;
  }

  return globalWire;

}


std::vector<std::vector<double> > cluster::BlurredClusteringAlg::GaussianBlur(std::vector<std::vector<double> > const& image) {

  /// Applies Gaussian blur to image

  if (fBlurSigma == 0) {
    //return (TH2F*) image->Clone(image->GetName() + TString("_blur"));
    return image;
  }

  // Find the blurring parameters
  int blurwire, blurtick, sigmawire, sigmatick;
  FindBlurringParameters(blurwire, blurtick, sigmawire, sigmatick);

  // Create Gaussian kernel
  //std::map<int,double> kernel;
  std::vector<double> kernel((2*blurwire+1)*(2*blurtick+1),0);
  int width = 2 * blurwire + 1;
  int height = 2 * blurtick + 1;

  // If the parameters match the last parameters, used the same last kernel
  if (fLastBlurWire == blurwire && fLastBlurTick == blurtick && fLastSigma == fBlurSigma && !fLastKernel.empty())
    kernel = fLastKernel;

  // Otherwise, compute a new kernel
  else {

    if (!fLastKernel.empty())
      fLastKernel.clear();

    // Smear out according to the blur radii in each direction
    for (int i = -blurwire; i <= blurwire; i++) {
      for (int j = -blurtick; j <= blurtick; j++) {

	// Fill kernel
	double sig2i = 2. * sigmawire * sigmawire;
	double sig2j = 2. * sigmatick * sigmatick;

	int key = (width * (j + blurtick)) + (i + blurwire);
	double value = 1. / sqrt(sig2i * M_PI) * exp(-i * i / sig2i) * 1. / sqrt(sig2j * M_PI) * exp(-j * j / sig2j);
	kernel.at(key) = value;

      }
    } // End loop over blurring region

    fLastKernel   = kernel;
    fLastBlurWire = blurwire;
    fLastBlurTick = blurtick;
    fLastSigma    = fBlurSigma;
  }

  // Return a convolution of this Gaussian kernel with the unblurred hit map
  return Convolve(image, kernel, width, height);//, TString::Format("%s_blur", image->GetName()));

}


double cluster::BlurredClusteringAlg::GetTimeOfBin(std::vector<std::vector<double> > const& image, int bin) {

  /// Returns the hit time of a hit in a particular bin

  double time = -10000;

  art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);
  if (!hit.isNull())
    time = hit->PeakTime();

  return time;

}


unsigned int cluster::BlurredClusteringAlg::NumNeighbours(int nbinsx, std::vector<bool> const& used, int bin) {

  /// Determines the number of clustered neighbours of a hit

  unsigned int neighbours = 0;

  // Loop over all directly neighbouring hits (not itself)
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      if (!x && !y) continue;

      // Determine bin
      int neighbouringBin = bin + x + (y * nbinsx); /// 2D hists can be considered a string of bins - the equation to convert between them is [bin = x + (nbinsx * y)]

      // If this bin is in the cluster, increase the neighbouring bin counter
      if (used.at(neighbouringBin))
	neighbours++;
    }
  }

  // Return the number of neighbours in the cluster of a particular hit
  return neighbours;
}


bool cluster::BlurredClusteringAlg::PassesTimeCut(std::vector<double> const& times, double time) {

  /// Determine if a hit is within a time threshold of any other hits in a cluster

  for (std::vector<double>::const_iterator timeIt = times.begin(); timeIt != times.end(); timeIt++) {
    if (std::abs(time - *timeIt) < fTimeThreshold) return true;
  }

  return false;
}


void cluster::BlurredClusteringAlg::SaveImage(TH2F* image, std::vector<art::PtrVector<recob::Hit> > const& allClusters, int pad, int tpc, int plane) {

  /// Save the images for debugging

  // Make a vector of clusters
  std::vector<std::vector<int> > allClusterBins;

  for (std::vector<art::PtrVector<recob::Hit> >::const_iterator clusterIt = allClusters.begin(); clusterIt != allClusters.end(); clusterIt++) {
    art::PtrVector<recob::Hit> cluster = *clusterIt;

    if (!cluster.size())
      continue;

    std::vector<int> clusterBins;

    for (art::PtrVector<recob::Hit>::iterator hitIt = cluster.begin(); hitIt != cluster.end(); hitIt++) {
      art::Ptr<recob::Hit> hit = *hitIt;
      unsigned int wire = GlobalWire(hit->WireID());
      float tick = hit->PeakTime();
      int bin = image->GetBin((wire-fLowerHistWire)+1,(tick-fLowerHistTick)+1);
      if (cluster.size() < fMinSize)
        bin *= -1;

      clusterBins.push_back(bin);
    }

    allClusterBins.push_back(clusterBins);
  }

  SaveImage(image, allClusterBins, pad, tpc, plane);
}

void cluster::BlurredClusteringAlg::SaveImage(TH2F* image, int pad, int tpc, int plane) {
  std::vector<std::vector<int> > allClusterBins;
  SaveImage(image, allClusterBins, pad, tpc, plane);
}

void cluster::BlurredClusteringAlg::SaveImage(TH2F* image, std::vector<std::vector<int> > const& allClusterBins, int pad, int tpc, int plane) {

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
  title << stage << " -- TPC " << tpc << ", Plane " << plane;// << " (Event " << fEvent << ")";

  image->SetName(title.str().c_str());
  image->SetTitle(title.str().c_str());
  image->DrawCopy("colz");

  // Draw the clustered hits on the histograms
  int clusterNum = 2;
  for (std::vector<std::vector<int> >::const_iterator it = allClusterBins.begin(); it != allClusterBins.end(); it++, clusterNum++) {
    std::vector<int> bins = *it;
    TMarker mark(0, 0, 20);
    mark.SetMarkerColor(clusterNum);
    mark.SetMarkerSize(0.1);

    for (std::vector<int>::iterator binIt = bins.begin(); binIt != bins.end(); binIt++) {
      int bin = *binIt;
      int wire, tick, z;

      // Hit from a cluster that we aren't going to save
      if (bin < 0) {
        bin *= -1;
        mark.SetMarkerStyle(24);
      }

      image->GetBinXYZ(bin,wire,tick,z);
      mark.DrawMarker(wire+fLowerHistWire-1, tick+fLowerHistTick-1);
      mark.SetMarkerStyle(20);
    }
  }

  if (pad == 4) {
    fDebugCanvas->Print(fDebugPDFName.c_str());
    fDebugCanvas->Clear("D");
  }

}

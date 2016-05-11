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

#include "larreco/RecoAlg/BlurredClusteringAlg.h"

cluster::BlurredClusteringAlg::BlurredClusteringAlg(fhicl::ParameterSet const& pset) {

  this->reconfigure(pset);
  this->MakeKernels();

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
}

void cluster::BlurredClusteringAlg::reconfigure(fhicl::ParameterSet const& p) {
  fBlurWire            = p.get<int>   ("BlurWire");
  fBlurTick            = p.get<int>   ("BlurTick");
  fSigmaWire           = p.get<double>("SigmaWire");
  fSigmaTick           = p.get<double>("SigmaTick");
  fMaxTickWidthBlur    = p.get<int>   ("MaxTickWidthBlur");
  fClusterWireDistance = p.get<int>   ("ClusterWireDistance");
  fClusterTickDistance = p.get<int>   ("ClusterTickDistance");
  fNeighboursThreshold = p.get<int>   ("NeighboursThreshold");
  fMinNeighbours       = p.get<int>   ("MinNeighbours");
  fMinSize             = p.get<int>   ("MinSize");
  fMinSeed             = p.get<double>("MinSeed");
  fTimeThreshold       = p.get<double>("TimeThreshold");
  fChargeThreshold     = p.get<double>("ChargeThreshold");
  fDebug               = p.get<bool>  ("Debug",false);
  fDetector            = p.get<std::string>("Detector","dune35t");

  fKernelWidth = 2 * fBlurWire + 1;
  fKernelHeight = 2 * fBlurTick*fMaxTickWidthBlur + 1;

  fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

void cluster::BlurredClusteringAlg::CreateDebugPDF(int run, int subrun, int event) {

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

  int wire = bin % image.size();
  int tick = bin / image.size();

  return fHitMap[wire][tick];

}

void cluster::BlurredClusteringAlg::ConvertBinsToClusters(std::vector<std::vector<double> > const& image,
							  std::vector<std::vector<int> > const& allClusterBins,
							  std::vector<art::PtrVector<recob::Hit> >& clusters) {

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

std::vector<std::vector<double> > cluster::BlurredClusteringAlg::ConvertRecobHitsToVector(std::vector<art::Ptr<recob::Hit> > const& hits) {

  // Define the size of this particular plane -- dynamically to avoid huge histograms
  int lowerTick = fDetProp->ReadOutWindowSize(), upperTick = 0, lowerWire = fGeom->MaxWires(), upperWire = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    int histWire = GlobalWire((*hitIt)->WireID());
    if ((*hitIt)->PeakTime() < lowerTick) lowerTick = (*hitIt)->PeakTime();
    if ((*hitIt)->PeakTime() > upperTick) upperTick = (*hitIt)->PeakTime();
    if (histWire < lowerWire) lowerWire = histWire;
    if (histWire > upperWire) upperWire = histWire;
    //std::cout << "Hit at " << (*hitIt)->WireID() << " has peak time " << (*hitIt)->PeakTime() << " and RMS " << (*hitIt)->RMS() << ": lower and upper tick is " << (*hitIt)->StartTick() << " and " << (*hitIt)->EndTick() << ")" << std::endl;
  }
  fLowerTick = lowerTick-20;
  fUpperTick = upperTick+20;
  fLowerWire = lowerWire-20;
  fUpperWire = upperWire+20;

  // Use a map to keep a track of the real hits and their wire/ticks
  fHitMap.clear();
  fHitMap.resize(fUpperWire-fLowerWire, std::vector<art::Ptr<recob::Hit> >(fUpperTick-fLowerTick, art::Ptr<recob::Hit>()));

  // Create a 2D vector
  std::vector<std::vector<double> > image(fUpperWire-fLowerWire, std::vector<double>(fUpperTick-fLowerTick, 0));

  // Look through the hits
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    int wire = GlobalWire((*hitIt)->WireID());
    int tick = (int)(*hitIt)->PeakTime();
    float charge = (*hitIt)->Integral();

    // Fill hit map and keep a note of all real hits for later
    if (charge > image.at(wire-fLowerWire).at(tick-fLowerTick)) {
      image.at(wire-fLowerWire).at(tick-fLowerTick) = charge;
      fHitMap[wire-fLowerWire][tick-fLowerTick] = (*hitIt);
    }
  }

  // Keep a note of dead wires
  fDeadWires.clear();
  fDeadWires.resize(fUpperWire-fLowerWire,false);
  geo::PlaneID planeID = hits.front()->WireID().planeID();

  for (int wire = fLowerWire; wire < fUpperWire; ++wire) {
    raw::ChannelID_t channel = fGeom->PlaneWireToChannel(planeID.Plane,wire,planeID.TPC,planeID.Cryostat);
    fDeadWires[wire-fLowerWire] = !fChanStatus.IsGood(channel);
  }

  return image;

}

int cluster::BlurredClusteringAlg::ConvertWireTickToBin(std::vector<std::vector<double> > const& image, int xbin, int ybin) {

  return (ybin * image.size()) + xbin;

}

double cluster::BlurredClusteringAlg::ConvertBinToCharge(std::vector<std::vector<double> > const& image, int bin) {

  int x = bin % image.size();
  int y = bin / image.size();

  return image.at(x).at(y);

}

std::pair<int,int> cluster::BlurredClusteringAlg::DeadWireCount(int wire_bin, int width) {

  std::pair<int,int> deadWires = std::make_pair<int,int>(0,0);

  int lower_bin = width / 2;
  int upper_bin = (width+1) / 2;

  for (int wire = TMath::Max(wire_bin + fLowerWire - lower_bin, fLowerWire); wire < TMath::Min(wire_bin + fLowerWire + upper_bin, fUpperWire); ++wire) {
    if (fDeadWires[wire-fLowerWire]) {
      if (wire < wire_bin + fLowerWire)
	++deadWires.first;
      else if (wire > wire_bin + fLowerWire)
	++deadWires.second;
    }
  }

  return deadWires;

}

void cluster::BlurredClusteringAlg::FindBlurringParameters(int& blur_wire, int& blur_tick, int& sigma_wire, int& sigma_tick) {

  // Calculate least squares slope
  int x, y;
  double nhits=0, sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (unsigned int wireIt = 0; wireIt < fHitMap.size(); ++wireIt) {
    for (unsigned int tickIt = 0; tickIt < fHitMap.at(wireIt).size(); ++tickIt) {
      if (fHitMap[wireIt][tickIt].isNull())
	continue;
      ++nhits;
      x = wireIt + fLowerWire;
      y = tickIt + fLowerTick;
      sumx += x;
      sumy += y;
      sumx2 += x*x;
      sumxy += x*y;
    }
  }
  double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);

  // Get the rough unit vector for the trajectories
  TVector2 unit = TVector2(1,gradient).Unit();

  // Catch vertical gradients
  if (std::isnan(gradient))
    unit = TVector2(0,1);

  // Use this direction to scale the blurring radii and Gaussian sigma
  blur_wire = std::max(std::abs(std::round(fBlurWire * unit.X())),1.);
  blur_tick = std::max(std::abs(std::round(fBlurTick * unit.Y())),1.);

  sigma_wire = std::max(std::abs(std::round(fSigmaWire * unit.X())),1.);
  sigma_tick = std::max(std::abs(std::round(fSigmaTick * unit.Y())),1.);

  // std::cout << "Gradient is " << gradient << ", giving x and y components " << unit.X() << " and " << unit.Y() << std::endl;
  // std::cout << "Blurring: wire " << blurwire << " and tick " << blurtick << "; sigma: wire " << sigmawire << " and tick " << sigmatick << std::endl;

  return;

}

int cluster::BlurredClusteringAlg::FindClusters(std::vector<std::vector<double> > const& blurred, std::vector<std::vector<int> >& allcluster) {

  // Vectors to hold cluster information
  std::vector<int> cluster;
  std::vector<double> times;

  // Size of image in x and y
  const int nbinsx = blurred.size();
  const int nbinsy = blurred.at(0).size();
  const int nbins = nbinsx * nbinsy;

  // Vectors to hold hit information
  std::vector<bool> used(nbins);
  std::vector<std::pair<double, int> > values;//(nbins);

  // Place the bin number and contents as a pair in the values vector
  for (int xbin = 0; xbin < nbinsx; ++xbin) {
    for (int ybin = 0; ybin < nbinsy; ++ybin) {
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
            if ( (x == binx and y == biny) or (x >= nbinsx or y >= nbinsy) or (x < 0 or y < 0) )
              continue;

	    // Get this bin
	    bin = ConvertWireTickToBin(blurred, x, y);
	    if (bin >= nbinsx * nbinsy or bin < 0)
	      continue;
            if (used[bin])
              continue;

	    // Get the blurred value and time for this bin
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

  double wireCentre[3];
  fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);

  double globalWire = -999;
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }
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
    else if (fDetector == "dunefd") {
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      // Detector geometry has four TPCs, two on top of each other, repeated along z...
      int block = wireID.TPC / 4;
      globalWire = (nwires*block) + wireID.Wire;
    }
    else {
      if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
      else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
    }
  }

  return std::round(globalWire);

}

std::vector<std::vector<double> > cluster::BlurredClusteringAlg::GaussianBlur(std::vector<std::vector<double> > const& image) {

  if (fSigmaWire == 0 and fSigmaTick == 0)
    return image;

  // Find the blurring parameters
  int blur_wire, blur_tick, sigma_wire, sigma_tick;
  FindBlurringParameters(blur_wire, blur_tick, sigma_wire, sigma_tick);

  // Convolve the Gaussian
  int width = 2 * blur_wire + 1;
  int height = 2 * blur_tick + 1;
  int nbinsx = image.size();
  int nbinsy = image.at(0).size();

  // Blurred histogram and normalisation for each bin
  std::vector<std::vector<double> > copy(nbinsx, std::vector<double>(nbinsy, 0));

  // Loop through all the bins in the histogram to blur
  for (int x = 0; x < nbinsx; ++x) {
    for (int y = 0; y < nbinsy; ++y) {

      if (image[x][y] == 0)
      	continue;

      // Scale the tick blurring based on the width of the hit
      int tick_scale = TMath::Sqrt(TMath::Power(fHitMap[x][y]->RMS(),2) + TMath::Power(sigma_tick,2)) / (double)sigma_tick;
      tick_scale = TMath::Max(TMath::Min(tick_scale,fMaxTickWidthBlur),1);
      std::vector<double> correct_kernel = fAllKernels[sigma_wire][sigma_tick*tick_scale];

      // Find any dead wires in the potential blurring region
      std::pair<int,int> num_deadwires = DeadWireCount(x, width);
      //num_deadwires = std::make_pair<int,int>(0,0);

      // Note of how many dead wires we have passed whilst blurring in the wire direction
      // If blurring below the seed hit, need to keep a note of how many dead wires to come
      // If blurring above, need to keep a note of how many dead wires have passed
      int dead_wires_passed = num_deadwires.first;

      // bool dead = false;
      // if (num_deadwires.first != 0 or num_deadwires.second != 0)
      // 	dead = true;

      // if (dead) {
      // 	std::cout << "Wire is " << x+fLowerWire << std::endl;
      // 	std::cout << "Width is " << width << std::endl;
      // }

      // Loop over the blurring region around this hit
      for (int blurx = -(width/2+num_deadwires.first); blurx < (width+1)/2+num_deadwires.second; ++blurx) {
  	for (int blury = -height/2*tick_scale; blury < ((((height+1)/2)-1)*tick_scale)+1; ++blury) {

	  // if (dead)
	  //   std::cout << "Start... dead_wires_passed is " << dead_wires_passed << " and blurx is " << blurx << std::endl;

	  if (blurx < 0 and fDeadWires[x+blurx])
	    dead_wires_passed -= 1;

  	  // Smear the charge of this hit
	  double weight = correct_kernel[fKernelWidth * (fKernelHeight / 2 + blury) + (fKernelWidth / 2 + (blurx - dead_wires_passed))];
  	  if (x + blurx >= 0 and x + blurx < nbinsx and y + blury >= 0 and y + blury < nbinsy)
  	    copy[x+blurx][y+blury] += weight * image[x][y];

	  if (blurx > 0 and fDeadWires[x+blurx])
	    dead_wires_passed += 1;

	  // if (dead)
	  //   std::cout << "Start... dead_wires_passed is " << dead_wires_passed << " and blurx is " << blurx << std::endl;

  	}
      } // blurring region

    }
  } // hits to blur

  // HAVE REMOVED NOMALISATION CODE
  // WHEN USING DIFFERENT KERNELS, THERE'S NO EASY WAY OF DOING THIS...
  // RECONSIDER...

  // Return the blurred histogram
  return copy;

}

double cluster::BlurredClusteringAlg::GetTimeOfBin(std::vector<std::vector<double> > const& image, int bin) {

  double time = -10000;

  art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);
  if (!hit.isNull())
    time = hit->PeakTime();

  return time;

}

void cluster::BlurredClusteringAlg::MakeKernels() {

  // Kernel size is the largest possible given the hit width rescaling
  fAllKernels.clear();
  fAllKernels.resize(fSigmaWire+1,std::vector<std::vector<double> >(fSigmaTick*fMaxTickWidthBlur+1,std::vector<double>(fKernelWidth*fKernelHeight,0.)));

  // Ranges of kernels to make
  // Complete range of sigmas possible after dynamic fixing and hit width convolution
  for (int sigma_wire = 1; sigma_wire <= fSigmaWire; ++sigma_wire) {
    for (int sigma_tick = 1; sigma_tick <= fSigmaTick*fMaxTickWidthBlur; ++sigma_tick) {

      // New kernel
      std::vector<double> kernel(fKernelWidth*fKernelHeight,0);

      // Smear out according to the blur radii in each direction
      for (int i = -fBlurWire; i <= fBlurWire; i++) {
	for (int j = -fBlurTick*fMaxTickWidthBlur; j <= fBlurTick*fMaxTickWidthBlur; j++) {

	  // Fill kernel
	  double sig2i = 2. * sigma_wire * sigma_wire;
	  double sig2j = 2. * sigma_tick * sigma_tick;

	  int key = (fKernelWidth * (j + fBlurTick*fMaxTickWidthBlur)) + (i + fBlurWire);
	  double value = 1. / sqrt(sig2i * M_PI) * exp(-i * i / sig2i) * 1. / sqrt(sig2j * M_PI) * exp(-j * j / sig2j);
	  kernel.at(key) = value;

	}
      } // End loop over blurring region

      fAllKernels[sigma_wire][sigma_tick] = kernel;

    }
  }

  return;

}

TH2F* cluster::BlurredClusteringAlg::MakeHistogram(std::vector<std::vector<double> > const& image, TString name) {

  TH2F* hist = new TH2F(name,name,fUpperWire-fLowerWire,fLowerWire-0.5,fUpperWire-0.5,fUpperTick-fLowerTick,fLowerTick-0.5,fUpperTick-0.5);
  hist->Clear();
  hist->SetXTitle("Wire number");
  hist->SetYTitle("Tick number");
  hist->SetZTitle("Charge");

  for (unsigned int imageWireIt = 0; imageWireIt < image.size(); ++imageWireIt) {
    int wire = imageWireIt + fLowerWire;
    for (unsigned int imageTickIt = 0; imageTickIt < image.at(imageWireIt).size(); ++imageTickIt) {
      int tick = imageTickIt + fLowerTick;
      hist->Fill(wire, tick, image.at(imageWireIt).at(imageTickIt));
    }
  }

  return hist;

}

unsigned int cluster::BlurredClusteringAlg::NumNeighbours(int nbinsx, std::vector<bool> const& used, int bin) {

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

  for (std::vector<double>::const_iterator timeIt = times.begin(); timeIt != times.end(); timeIt++) {
    if (std::abs(time - *timeIt) < fTimeThreshold) return true;
  }

  return false;
}

void cluster::BlurredClusteringAlg::SaveImage(TH2F* image, std::vector<art::PtrVector<recob::Hit> > const& allClusters, int pad, int tpc, int plane) {

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
      int bin = image->GetBin((wire-fLowerWire)+1,(tick-fLowerTick)+1);
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
      mark.DrawMarker(wire+fLowerWire-1, tick+fLowerTick-1);
      mark.SetMarkerStyle(20);
    }
  }

  if (pad == 4) {
    fDebugCanvas->Print(fDebugPDFName.c_str());
    fDebugCanvas->Clear("D");
  }

}

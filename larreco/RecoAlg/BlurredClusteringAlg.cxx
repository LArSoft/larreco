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

#include "cetlib/pow.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/RecoAlg/BlurredClusteringAlg.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
  // Since we perform a static_cast of an unsigned int to an int, we
  // must make sure that the configuration does not produce data loss
  // during the cast.
  unsigned int verify_safe_cast(unsigned int const value,
                                std::string const& parameter_name)
  {
    if (value > std::numeric_limits<int>::max()) {
      throw art::Exception{art::errors::Configuration}
      << "The " << parameter_name << " value (" << value << ") exceeds the maximum value for an integer:\n"
      << std::numeric_limits<int>::max() << '\n';
    }
    return value;
  }
}

cluster::BlurredClusteringAlg::BlurredClusteringAlg(fhicl::ParameterSet const& pset)
  : fDebug{pset.get<bool>("Debug", false)}
  , fDetector{pset.get<std::string>("Detector", "dune35t")}
  , fBlurWire{verify_safe_cast(pset.get<unsigned int>("BlurWire"), "BlurWire")}
  , fBlurTick{pset.get<unsigned int>("BlurTick")}
  , fSigmaWire{pset.get<double>("SigmaWire")}
  , fSigmaTick{pset.get<double>("SigmaTick")}
  , fMaxTickWidthBlur{pset.get<unsigned int>("MaxTickWidthBlur")}
  , fClusterWireDistance{pset.get<unsigned int>("ClusterWireDistance")}
  , fClusterTickDistance{pset.get<unsigned int>("ClusterTickDistance")}
  , fNeighboursThreshold{pset.get<unsigned int>("NeighboursThreshold")}
  , fMinNeighbours{pset.get<unsigned int>("MinNeighbours")}
  , fMinSize{pset.get<unsigned int>("MinSize")}
  , fMinSeed{pset.get<double>("MinSeed")}
  , fTimeThreshold{pset.get<double>("TimeThreshold")}
  , fChargeThreshold{pset.get<double>("ChargeThreshold")}
  , fKernelWidth{2 * fBlurWire + 1}
  , fKernelHeight{2 * fBlurTick*fMaxTickWidthBlur + 1}
  , fDetProp{lar::providerFrom<detinfo::DetectorPropertiesService>()}
{
  verify_safe_cast(fBlurTick*fMaxTickWidthBlur, "fBlurTick*fMaxTickWidthBlur");
  this->MakeKernels();
}

cluster::BlurredClusteringAlg::~BlurredClusteringAlg()
{
  if (fDebugCanvas) {
    std::string closeName = fDebugPDFName;
    closeName.append("]");
    fDebugCanvas->Print(closeName.c_str());
    delete fDebugCanvas;
  }
}

void
cluster::BlurredClusteringAlg::CreateDebugPDF(int run, int subrun, int event)
{
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

art::PtrVector<recob::Hit>
cluster::BlurredClusteringAlg::ConvertBinsToRecobHits(std::vector<std::vector<double>> const& image,
                                                      std::vector<unsigned int> const& bins) const
{
  art::PtrVector<recob::Hit> hits;

  for (auto const bin : bins) {
    art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);

    if (!hit.isNull())
      hits.push_back(hit);
  }

  // Return the vector of hits to make cluster
  return hits;
}

art::Ptr<recob::Hit>
cluster::BlurredClusteringAlg::ConvertBinToRecobHit(std::vector<std::vector<double>> const& image,
                                                    unsigned int const bin) const
{
  unsigned int const wire = bin % image.size();
  unsigned int const tick = bin / image.size();

  return fHitMap[wire][tick];
}

void
cluster::BlurredClusteringAlg::ConvertBinsToClusters(std::vector<std::vector<double>> const& image,
                                                     std::vector<std::vector<unsigned int>> const& allClusterBins,
                                                     std::vector<art::PtrVector<recob::Hit>>& clusters) const
{
  // Loop through the clusters (each a vector of bins)
  for (auto const& bins : allClusterBins) {
    // Convert the clusters (vectors of bins) to hits in a vector of recob::Hits
    art::PtrVector<recob::Hit> clusHits = ConvertBinsToRecobHits(image, bins);

    mf::LogInfo("BlurredClustering") << "Cluster made from " << bins.size() << " bins, of which " << clusHits.size() << " were real hits";

    // Make sure the clusters are above the minimum cluster size
    if (clusHits.size() < fMinSize) {
      mf::LogVerbatim("BlurredClustering") << "Cluster of size " << clusHits.size() << " not saved since it is smaller than the minimum cluster size, set to " << fMinSize;
      continue;
    }

    clusters.push_back(std::move(clusHits));
  }
}

std::vector<std::vector<double>>
cluster::BlurredClusteringAlg::ConvertRecobHitsToVector(std::vector<art::Ptr<recob::Hit>> const& hits) {
  // Define the size of this particular plane -- dynamically to avoid huge histograms
  unsigned int lowerTick = fDetProp->ReadOutWindowSize(), upperTick = 0, lowerWire = fGeom->MaxWires(), upperWire = 0;
  for (auto const& hit : hits) {
    unsigned int histWire = GlobalWire(hit->WireID());
    if (hit->PeakTime() < lowerTick) lowerTick = hit->PeakTime();
    if (hit->PeakTime() > upperTick) upperTick = hit->PeakTime();
    if (histWire < lowerWire) lowerWire = histWire;
    if (histWire > upperWire) upperWire = histWire;
  }
  fLowerTick = lowerTick-20;
  fUpperTick = upperTick+20;
  fLowerWire = lowerWire-20;
  fUpperWire = upperWire+20;

  // Use a map to keep a track of the real hits and their wire/ticks
  fHitMap.clear();
  fHitMap.resize(fUpperWire-fLowerWire, std::vector<art::Ptr<recob::Hit>>(fUpperTick-fLowerTick, art::Ptr<recob::Hit>()));

  // Create a 2D vector
  std::vector<std::vector<double>> image(fUpperWire-fLowerWire, std::vector<double>(fUpperTick-fLowerTick, 0));

  // Look through the hits
  for (auto const& hit : hits) {
    int wire = GlobalWire(hit->WireID());
    int tick = (int)hit->PeakTime();
    float charge = hit->Integral();

    // Fill hit map and keep a note of all real hits for later
    if (charge > image.at(wire-fLowerWire).at(tick-fLowerTick)) {
      image.at(wire-fLowerWire).at(tick-fLowerTick) = charge;
      fHitMap[wire-fLowerWire][tick-fLowerTick] = hit;
    }
  }

  // Keep a note of dead wires
  fDeadWires.clear();
  fDeadWires.resize(fUpperWire-fLowerWire,false);
  geo::PlaneID planeID = hits.front()->WireID().planeID();

  for (unsigned wire = fLowerWire; wire < fUpperWire; ++wire) {
    raw::ChannelID_t channel = fGeom->PlaneWireToChannel(planeID.Plane,wire,planeID.TPC,planeID.Cryostat);
    // 'wire' guaranteed to be greater than or equal to fLowerWire.
    fDeadWires[wire-fLowerWire] = !fChanStatus.IsGood(channel);
  }

  return image;
}

unsigned int
cluster::BlurredClusteringAlg::ConvertWireTickToBin(std::vector<std::vector<double>> const& image,
                                                    unsigned int const xbin,
                                                    unsigned int const ybin) const
{
  return (ybin * image.size()) + xbin;
}

double
cluster::BlurredClusteringAlg::ConvertBinToCharge(std::vector<std::vector<double>> const& image,
                                                  unsigned int const bin) const
{
  unsigned int const x = bin % image.size();
  unsigned int const y = bin / image.size();

  return image.at(x).at(y);
}

std::pair<unsigned int, unsigned int>
cluster::BlurredClusteringAlg::DeadWireCount(unsigned int const wire_bin,
                                             unsigned int const width) const
{
  auto deadWires = std::make_pair(0u, 0u);

  unsigned int const lower_bin = width / 2;
  unsigned int const upper_bin = (width+1) / 2;

  auto const offset = wire_bin + fLowerWire;
  for (auto wire = std::max(offset - lower_bin, fLowerWire); wire < std::min(offset + upper_bin, fUpperWire); ++wire) {
    // 'wire is guaranteed to be greater than or equal to fLowerWire
    if (fDeadWires[wire-fLowerWire]) {
      if (wire < offset)
        ++deadWires.first;
      else if (wire > offset)
        ++deadWires.second;
    }
  }

  return deadWires;
}

std::array<unsigned int, 4>
cluster::BlurredClusteringAlg::FindBlurringParameters()
{
  // Calculate least squares slope
  double nhits=0, sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (unsigned int wireIt = 0; wireIt < fHitMap.size(); ++wireIt) {
    for (unsigned int tickIt = 0; tickIt < fHitMap.at(wireIt).size(); ++tickIt) {
      if (fHitMap[wireIt][tickIt].isNull())
        continue;
      ++nhits;
      unsigned int const x = wireIt + fLowerWire;
      unsigned int const y = tickIt + fLowerTick;
      sumx += x;
      sumy += y;
      sumx2 += x*x;
      sumxy += x*y;
    }
  }
  double const gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);

  // Get the rough unit vector for the trajectories
  //   - catch vertical gradients
  TVector2 const unit = std::isnan(gradient) ? TVector2(0, 1) : TVector2(1,gradient).Unit();

  unsigned int const blur_wire = std::max(std::abs(std::round(fBlurWire * unit.X())), 1.);
  unsigned int const blur_tick = std::max(std::abs(std::round(fBlurTick * unit.Y())), 1.);
  unsigned int const sigma_wire = std::max(std::abs(std::round(fSigmaWire * unit.X())), 1.);
  unsigned int const sigma_tick = std::max(std::abs(std::round(fSigmaTick * unit.Y())), 1.);

  // Use this direction to scale the blurring radii and Gaussian sigma
  return {{blur_wire, blur_tick, sigma_wire, sigma_tick}};
}

unsigned int
cluster::BlurredClusteringAlg::FindClusters(std::vector<std::vector<double>> const& blurred,
                                            std::vector<std::vector<unsigned int>>& allcluster)
{
  // Size of image in x and y
  unsigned const nbinsx = blurred.size();
  unsigned const nbinsy = blurred.at(0).size();
  unsigned const nbins = nbinsx * nbinsy;

  // Vectors to hold hit information
  std::vector<bool> used(nbins);
  std::vector<std::pair<double, int>> values;//(nbins);

  // Place the bin number and contents as a pair in the values vector
  for (unsigned int xbin = 0; xbin < nbinsx; ++xbin) {
    for (unsigned int ybin = 0; ybin < nbinsy; ++ybin) {
      unsigned int const bin = ConvertWireTickToBin(blurred, xbin, ybin);
      values.emplace_back(ConvertBinToCharge(blurred, bin), bin);
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
    std::vector<unsigned int> cluster;
    std::vector<double> times;

    // Get the highest charge bin (go no further if below seed threshold)
    double blurred_binval = values[niter].first;
    if (blurred_binval < fMinSeed)
      break;

    // Iterate through the bins from highest charge down
    int const bin = values[niter++].second;

    // Put this bin in used if not already there
    if (used[bin])
      continue;
    used[bin] = true;

    // Start a new cluster
    cluster.push_back(bin);

    // Get the time of this hit
    double const time = GetTimeOfBin(blurred, bin);
    if (time > 0)
      times.push_back(time);


    // Now cluster neighbouring hits to this seed
    while (true) {

      bool added_bin{false};

      // N.B. We cannot use a range-for because the cluster container
      // can be extended during this loop.
      for (unsigned int clusBin = 0; clusBin < cluster.size(); ++clusBin) {

        // Get x and y values for bin (c++ returns a%b = a if a<b)
        unsigned int const binx = cluster[clusBin] % nbinsx;
        unsigned int const biny = ((cluster[clusBin] - binx) / nbinsx) % nbinsy;

        // Look for hits in the neighbouring x/y bins
        auto const x_start = binx < fClusterWireDistance ? 0 : binx - fClusterWireDistance;
        auto const y_start = biny < fClusterTickDistance ? 0 : biny - fClusterTickDistance;
        for (auto x = x_start; x <= binx + fClusterWireDistance; ++x) {
          for (auto y = y_start; y <= biny + fClusterTickDistance; ++y) {
            if ((x == binx and y == biny) or (x >= nbinsx or y >= nbinsy))
              continue;

            // Get this bin
            unsigned int const bin = ConvertWireTickToBin(blurred, x, y);
            if (bin >= nbinsx * nbinsy or bin < 0)
              continue;
            if (used[bin])
              continue;

            // Get the blurred value and time for this bin
            int const blurred_binval = ConvertBinToCharge(blurred, bin);
            int const time = GetTimeOfBin(blurred, bin); // NB for 'fake' hits, time is defaulted to -10000

            // Check real hits pass time cut (ignores fake hits)
            if (time > 0 && times.size() > 0 && ! PassesTimeCut(times, time))
              continue;

            // Add to cluster if bin value is above threshold
            if (blurred_binval > fChargeThreshold) {
              used[bin] = true;
              cluster.push_back(bin);
              added_bin = true;
              if (time > 0) {
                times.push_back(time);
              }
            } // End of adding blurred bin to cluster

          }
        } // End of looking at directly neighbouring bins

      } // End of looping over bins already in this cluster

      if (!added_bin)
        break;

    } // End of adding hits to this cluster

    // Check this cluster is above minimum size
    if (cluster.size() < fMinSize) {
      for (unsigned int i = 0; i < cluster.size(); i++)
        used[cluster[i]] = false;
      continue;
    }

    // Fill in holes in the cluster
    // N.B. We cannot use a range-for because the cluster container
    // can be extended during this loop.
    for (unsigned int clusBin = 0; clusBin < cluster.size(); clusBin++) {

      // Looks at directly neighbouring bins (and not itself)
      for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
          if (x == 0 && y == 0) continue;

          // Look at neighbouring bins to the clustered bin which are inside the cluster
          unsigned int const neighbouringBin = cluster[clusBin] + x + (y * nbinsx);
          if (neighbouringBin < nbinsx || neighbouringBin % nbinsx == 0 || neighbouringBin % nbinsx == nbinsx - 1 || neighbouringBin >= nbinsx * (nbinsy - 1))
            continue;

          double const time = GetTimeOfBin(blurred, neighbouringBin);

          // If not already clustered and passes neighbour/time thresholds, add to cluster
          if (!used[neighbouringBin] && (NumNeighbours(nbinsx, used, neighbouringBin) > fNeighboursThreshold) && PassesTimeCut(times, time)) {
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
      bool removed_bins{false};
      // Loop over all the bins in the cluster
      for (int clusBin = cluster.size() - 1; clusBin >= 0; clusBin--) {
        auto const bin = cluster[clusBin];

        // If bin is in cluster ignore
        if (bin < nbinsx || bin % nbinsx == 0 || bin % nbinsx == nbinsx - 1 || bin >= nbinsx * (nbinsy - 1)) continue;

        // Remove hit if it has too few neighbouring hits
        if (NumNeighbours(nbinsx, used, bin) < fMinNeighbours) {
          used[bin] = false;
          removed_bins = true;
          cluster.erase(cluster.begin() + clusBin);
          blurred_binval = ConvertBinToCharge(blurred, bin);
        }
      }

      if (!removed_bins) break;
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

  return allcluster.size();
}

int
cluster::BlurredClusteringAlg::GlobalWire(const geo::WireID& wireID)
{
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
      else mf::LogError("BlurredClusterAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC << " (geometry " << fDetector << ")";
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

  return std::round(globalWire);
}

std::vector<std::vector<double>>
cluster::BlurredClusteringAlg::GaussianBlur(std::vector<std::vector<double>> const& image)
{
  if (fSigmaWire == 0 and fSigmaTick == 0)
    return image;

  // Find the blurring parameters
  auto const [blur_wire, blur_tick, sigma_wire, sigma_tick] = FindBlurringParameters();

  // Convolve the Gaussian
  unsigned int const width = 2 * blur_wire + 1;
  unsigned int const height = 2 * blur_tick + 1;
  unsigned int const nbinsx = image.size();
  unsigned int const nbinsy = image.at(0).size();

  // Blurred histogram and normalisation for each bin
  std::vector<std::vector<double>> copy(nbinsx, std::vector<double>(nbinsy, 0));

  // Loop through all the bins in the histogram to blur
  for (unsigned int x = 0; x < nbinsx; ++x) {
    for (unsigned int y = 0; y < nbinsy; ++y) {
      if (image[x][y] == 0)
        continue;

      // Scale the tick blurring based on the width of the hit
      unsigned int tick_scale = std::sqrt(cet::square(fHitMap[x][y]->RMS()) + cet::square(sigma_tick)) / (double)sigma_tick;
      tick_scale = std::max(std::min(tick_scale, fMaxTickWidthBlur), 1u);
      std::vector<double> correct_kernel = fAllKernels[sigma_wire][sigma_tick*tick_scale];

      // Find any dead wires in the potential blurring region
      auto const [num_under_offset, num_over_offset] = DeadWireCount(x, width);

      // Note of how many dead wires we have passed whilst blurring in the wire direction
      // If blurring below the seed hit, need to keep a note of how many dead wires to come
      // If blurring above, need to keep a note of how many dead wires have passed
      int dead_wires_passed = num_under_offset;

      // Loop over the blurring region around this hit
      auto const blurx_end = static_cast<int>((width+1)/2+num_over_offset);
      auto const blury_end = static_cast<int>(((((height+1)/2)-1)*tick_scale)+1);
      for (int blurx = -(width/2+num_under_offset); blurx < blurx_end; ++blurx) {
        if (x + blurx < 0 or x + blurx >= nbinsx) continue;
        for (int blury = -height/2*tick_scale; blury < blury_end; ++blury) {
          if (y + blury < 0 or y + blury >= nbinsy) continue;

          if (blurx < 0 and fDeadWires[x+blurx])
            dead_wires_passed -= 1;

          // Smear the charge of this hit
          double weight = correct_kernel[fKernelWidth * (fKernelHeight / 2 + blury) + (fKernelWidth / 2 + (blurx - dead_wires_passed))];
          copy[x+blurx][y+blury] += weight * image[x][y];

          if (blurx > 0 and fDeadWires[x+blurx])
            dead_wires_passed += 1;

        }
      } // blurring region
    }
  } // hits to blur

  return copy;
}

double
cluster::BlurredClusteringAlg::GetTimeOfBin(std::vector<std::vector<double>> const& image,
                                            int const bin)
{
  art::Ptr<recob::Hit> const hit = ConvertBinToRecobHit(image, bin);
  return hit.isNull() ? -10000 : hit->PeakTime();
}

void
cluster::BlurredClusteringAlg::MakeKernels()
{
  // Kernel size is the largest possible given the hit width rescaling
  fAllKernels.resize(fSigmaWire+1, std::vector<std::vector<double>>(fSigmaTick*fMaxTickWidthBlur+1,std::vector<double>(fKernelWidth*fKernelHeight,0.)));

  // Ranges of kernels to make
  // Complete range of sigmas possible after dynamic fixing and hit width convolution
  for (unsigned int sigma_wire = 1; sigma_wire <= fSigmaWire; ++sigma_wire) {
    for (unsigned int sigma_tick = 1; sigma_tick <= fSigmaTick*fMaxTickWidthBlur; ++sigma_tick) {

      // New kernel
      std::vector<double> kernel(fKernelWidth*fKernelHeight,0);

      // Smear out according to the blur radii in each direction

      // We perform a static_cast to int because it is unlikely that a
      // user would ever specify a value of fBlurWire, fBlurTick, and
      // fMaxTickWidthBlur that would exceed result in a value greater
      // than std::numeric_limits<int>::max().
      for (int i = -fBlurWire; i <= static_cast<int>(fBlurWire); ++i) {
        for (int j = -fBlurTick*fMaxTickWidthBlur; j <= static_cast<int>(fBlurTick*fMaxTickWidthBlur); ++j) {

          // Fill kernel
          double const sig2i = 2. * sigma_wire * sigma_wire;
          double const sig2j = 2. * sigma_tick * sigma_tick;

          int const key = (fKernelWidth * (j + fBlurTick*fMaxTickWidthBlur)) + (i + fBlurWire);
          double const value = 1. / sqrt(sig2i * M_PI) * exp(-i * i / sig2i) * 1. / sqrt(sig2j * M_PI) * exp(-j * j / sig2j);
          kernel.at(key) = value;
        }
      } // End loop over blurring region

      fAllKernels[sigma_wire][sigma_tick] = kernel;
    }
  }
}

TH2F*
cluster::BlurredClusteringAlg::MakeHistogram(std::vector<std::vector<double>> const& image,
                                             TString const name)
{
  auto hist = new TH2F(name, name,
                       fUpperWire-fLowerWire, fLowerWire-0.5, fUpperWire-0.5,
                       fUpperTick-fLowerTick, fLowerTick-0.5, fUpperTick-0.5);
  hist->SetXTitle("Wire number");
  hist->SetYTitle("Tick number");
  hist->SetZTitle("Charge");

  for (unsigned int imageWireIt = 0; imageWireIt < image.size(); ++imageWireIt) {
    int const wire = imageWireIt + fLowerWire;
    for (unsigned int imageTickIt = 0; imageTickIt < image.at(imageWireIt).size(); ++imageTickIt) {
      int const tick = imageTickIt + fLowerTick;
      hist->Fill(wire, tick, image.at(imageWireIt).at(imageTickIt));
    }
  }

  return hist;
}

unsigned int
cluster::BlurredClusteringAlg::NumNeighbours(int const nbinsx,
                                             std::vector<bool> const& used,
                                             int const bin)
{
  unsigned int neighbours = 0;

  // Loop over all directly neighbouring hits (not itself)
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      if (x == 0 && y == 0) continue;

      // Determine bin
      int neighbouringBin = bin + x + (y * nbinsx); /// 2D hists can be considered a string of bins - the equation to convert between them is [bin = x + (nbinsx * y)]

      // If this bin is in the cluster, increase the neighbouring bin counter
      if (used.at(neighbouringBin))
        ++neighbours;
    }
  }

  // Return the number of neighbours in the cluster of a particular hit
  return neighbours;
}

bool
cluster::BlurredClusteringAlg::PassesTimeCut(std::vector<double> const& times,
                                             double const time)
{
  for (auto const t : times) {
    if (std::abs(time - t) < fTimeThreshold) return true;
  }
  return false;
}

void
cluster::BlurredClusteringAlg::SaveImage(TH2F* image,
                                         std::vector<art::PtrVector<recob::Hit>> const& allClusters,
                                         int const pad,
                                         int const tpc,
                                         int const plane)
{
  // Make a vector of clusters
  std::vector<std::vector<unsigned int>> allClusterBins;

  for (auto const& cluster : allClusters) {
    if (cluster.empty())
      continue;

    std::vector<unsigned int> clusterBins;
    for (auto const& hit : cluster) {
      unsigned int wire = GlobalWire(hit->WireID());
      float tick = hit->PeakTime();
      // ROOT will return a non-negative value for the bin number.
      unsigned int const bin = image->GetBin(wire-fLowerWire+1, tick-fLowerTick+1);
      clusterBins.push_back(bin);
    }

    allClusterBins.push_back(clusterBins);
  }

  SaveImage(image, allClusterBins, pad, tpc, plane);
}

void
cluster::BlurredClusteringAlg::SaveImage(TH2F* image,
                                         int const pad,
                                         int const tpc,
                                         int const plane)
{
  std::vector<std::vector<unsigned int>> allClusterBins;
  SaveImage(image, allClusterBins, pad, tpc, plane);
}

void
cluster::BlurredClusteringAlg::SaveImage(TH2F* image,
                                         std::vector<std::vector<unsigned int>> const& allClusterBins,
                                         int const pad,
                                         int const tpc,
                                         int const plane)
{
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
  for (auto it = cbegin(allClusterBins), e = cend(allClusterBins); it != e; ++it, ++clusterNum) {
    auto const& bins = *it;

    // Distinguish between clusters that are not saved (i.e. failing
    // the minimimum size) and those that are.
    auto const markerStyle = bins.size() < fMinSize ? 24 : 20;

    TMarker mark(0, 0, markerStyle);
    mark.SetMarkerColor(clusterNum);
    mark.SetMarkerSize(0.1);

    for (auto const bin : bins) {
      int wire, tick, z;
      image->GetBinXYZ(bin, wire, tick, z);
      mark.DrawMarker(wire+fLowerWire-1, tick+fLowerTick-1);
    }
  }

  if (pad == 4) {
    fDebugCanvas->Print(fDebugPDFName.c_str());
    fDebugCanvas->Clear("D");
  }
}

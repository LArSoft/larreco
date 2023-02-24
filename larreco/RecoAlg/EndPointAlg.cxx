////////////////////////////////////////////////////////////////////////
//
// EndPointAlg class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find (weak) vertices from hits after deconvolution and hit finding.
//  A weak end point is a end point that has been found using a dedicated end point finding algorithm only. A
//  strong end point is a end point that has been found using a dedicated end point finding algorithm and matched
//  to a crossing of two or more HoughLineFinder lines. The End PointMatch module finds strong vertices.
////////////////////////////////////////////////////////////////////////
/// The algorithm is based on:
///C. Harris and M. Stephens (1988). "A combined corner and edge detector". Proceedings of the 4th Alvey
///Vision Conference. pp. 147-151.
///B. Morgan (2010). "Interest Point Detection for Reconstruction in High Granularity Tracking Detectors".
///arXiv:1006.3012v1 [physics.ins-det]
//Thanks to B. Morgan of U. of Warwick for comments and suggestions

#include <algorithm>
#include <fstream>
#include <math.h>
#include <vector>

#include "TMath.h"
#include "TString.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/EndPointAlg.h"

//-----------------------------------------------------------------------------
cluster::EndPointAlg::EndPointAlg(fhicl::ParameterSet const& p)
{
  fTimeBins = p.get<int>("TimeBins");
  fMaxCorners = p.get<int>("MaxCorners");
  fGsigma = p.get<double>("Gsigma");
  fWindow = p.get<int>("Window");
  fThreshold = p.get<double>("Threshold");
  fSaveVertexMap = p.get<int>("SaveVertexMap");
}

//-----------------------------------------------------------------------------
double cluster::EndPointAlg::Gaussian(int x, int y, double sigma) const
{
  double const Norm = 1. / std::sqrt(2 * TMath::Pi() * cet::square(sigma));
  return Norm * exp(-cet::sum_of_squares(x, y) / (2 * cet::square(sigma)));
}

//-----------------------------------------------------------------------------
double cluster::EndPointAlg::GaussianDerivativeX(int x, int y) const
{
  double const Norm = 1. / (std::sqrt(2 * TMath::Pi()) * cet::cube(fGsigma));
  return Norm * (-x) * exp(-cet::sum_of_squares(x, y) / (2 * cet::square(fGsigma)));
}

//-----------------------------------------------------------------------------
double cluster::EndPointAlg::GaussianDerivativeY(int x, int y) const
{
  double const Norm = 1. / (std::sqrt(2 * TMath::Pi()) * cet::cube(fGsigma));
  return Norm * (-y) * exp(-cet::sum_of_squares(x, y) / (2 * cet::square(fGsigma)));
}

//-----------------------------------------------------------------------------
//this method saves a BMP image of the vertex map space, which can be viewed with gimp
void cluster::EndPointAlg::VSSaveBMPFile(const char* fileName,
                                         unsigned char* pix,
                                         int dx,
                                         int dy) const
{
  std::ofstream bmpFile(fileName, std::ios::binary);
  bmpFile.write("B", 1);
  bmpFile.write("M", 1);
  int bitsOffset = 54 + 256 * 4;
  int size = bitsOffset + dx * dy; //header plus 256 entry LUT plus pixels
  bmpFile.write((const char*)&size, 4);
  int reserved = 0;
  bmpFile.write((const char*)&reserved, 4);
  bmpFile.write((const char*)&bitsOffset, 4);
  int bmiSize = 40;
  bmpFile.write((const char*)&bmiSize, 4);
  bmpFile.write((const char*)&dx, 4);
  bmpFile.write((const char*)&dy, 4);
  short planes = 1;
  bmpFile.write((const char*)&planes, 2);
  short bitCount = 8;
  bmpFile.write((const char*)&bitCount, 2);
  int i, temp = 0;
  for (i = 0; i < 6; i++)
    bmpFile.write((const char*)&temp, 4); // zero out optional color info
  // write a linear LUT
  char lutEntry[4]; // blue,green,red
  lutEntry[3] = 0;  // reserved part
  for (i = 0; i < 256; i++) {
    lutEntry[0] = i;
    lutEntry[1] = i + 1;
    lutEntry[2] = i + 2;
    bmpFile.write(lutEntry, sizeof lutEntry);
  }
  // write the actual pixels
  bmpFile.write((const char*)pix, dx * dy);
}

//......................................................
size_t cluster::EndPointAlg::EndPoint(const art::PtrVector<recob::Cluster>& clusIn,
                                      std::vector<recob::EndPoint2D>& vtxcol,
                                      std::vector<art::PtrVector<recob::Hit>>& vtxHitsOut,
                                      art::Event const& evt,
                                      std::string const& label) const
{
  auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const det_prop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);

  //Point to a collection of vertices to output.
  std::vector<art::Ptr<recob::Hit>> hit;

  art::FindManyP<recob::Hit> fmh(clusIn, evt, label);

  int flag = 0;
  int windex = 0; // the wire index to make sure the end point finder does not
                  // fall off the edge of the hit map
  int tindex = 0; // the time index to make sure the end point finder does not
                  // fall off the edge of the hit map
  int n = 0;      // index of window cell. There are 49 cells in the 7X7 Gaussian and
                  // Gaussian derivative windows
  unsigned int numberwires = 0;
  const unsigned int numbertimesamples = det_prop.ReadOutWindowSize();
  const float TicksPerBin = numbertimesamples / fTimeBins;

  double MatrixAAsum = 0.;
  double MatrixBBsum = 0.;
  double MatrixCCsum = 0.;
  std::vector<double> Cornerness2;
  //gaussian window definitions. The cell weights are calculated here to help the algorithm's speed
  double w[49] = {0.};
  double wx[49] = {0.};
  double wy[49] = {0.};
  int ctr = 0;
  for (int i = -3; i < 4; ++i) {
    for (int j = 3; j > -4; --j) {
      w[ctr] = Gaussian(i, j, fGsigma);
      wx[ctr] = GaussianDerivativeX(i, j);
      wy[ctr] = GaussianDerivativeY(i, j);
      ++ctr;
    }
  }

  unsigned int wire = 0;
  unsigned int wire2 = 0;
  for (auto view : wireReadoutGeom.Views()) {
    art::PtrVector<recob::Hit> vHits;
    art::PtrVector<recob::Cluster>::const_iterator clusterIter = clusIn.begin();
    hit.clear();
    size_t cinctr = 0;
    while (clusterIter != clusIn.end()) {
      if ((*clusterIter)->View() == view) {
        hit = fmh.at(cinctr);
      } //end if cluster is in the correct view
      clusterIter++;
      ++cinctr;
    }

    if (hit.size() == 0) continue;

    geo::WireID wid = hit[0]->WireID();
    numberwires = wireReadoutGeom.Plane(wid.asPlaneID()).Nwires();
    mf::LogInfo("EndPointAlg") << " --- endpoints check " << numberwires << " " << numbertimesamples
                               << " " << fTimeBins;

    std::vector<std::vector<double>> MatrixAsum(numberwires);
    std::vector<std::vector<double>> MatrixBsum(numberwires);
    std::vector<std::vector<double>> hit_map(numberwires); //the map of hits

    //the index of the hit that corresponds to the potential corner
    std::vector<std::vector<int>> hit_loc(numberwires);

    std::vector<std::vector<double>> Cornerness(numberwires); //the "weight" of a corner

    for (unsigned int wi = 0; wi < numberwires; ++wi) {
      hit_map[wi].resize(fTimeBins, 0);
      hit_loc[wi].resize(fTimeBins, -1);
      Cornerness[wi].resize(fTimeBins, 0);
      MatrixAsum[wi].resize(fTimeBins, 0);
      MatrixBsum[wi].resize(fTimeBins, 0);
    }
    for (unsigned int i = 0; i < hit.size(); ++i) {
      wire = hit[i]->WireID().Wire;
      const float center = hit[i]->PeakTime(), sigma = hit[i]->RMS();
      const int iFirstBin = int((center - 3 * sigma) / TicksPerBin),
                iLastBin = int((center + 3 * sigma) / TicksPerBin);
      for (int iBin = iFirstBin; iBin <= iLastBin; ++iBin) {
        const float bin_center = iBin * TicksPerBin;
        hit_map[wire][iBin] += Gaussian(bin_center, center, sigma);
      }
    }

    // Gaussian derivative convolution
    for (unsigned int wire = 1; wire < numberwires - 1; ++wire) {

      for (int timebin = 1; timebin < fTimeBins - 1; ++timebin) {
        MatrixAsum[wire][timebin] = 0.;
        MatrixBsum[wire][timebin] = 0.;
        n = 0;
        for (int i = -3; i <= 3; ++i) {
          windex = wire + i;
          if (windex < 0) windex = 0;
          // this is ok, because the line before makes sure it's not negative
          else if ((unsigned int)windex >= numberwires)
            windex = numberwires - 1;

          for (int j = -3; j <= 3; ++j) {
            tindex = timebin + j;
            if (tindex < 0)
              tindex = 0;
            else if (tindex >= fTimeBins)
              tindex = fTimeBins - 1;

            MatrixAsum[wire][timebin] += wx[n] * hit_map[windex][tindex];
            MatrixBsum[wire][timebin] += wy[n] * hit_map[windex][tindex];
            ++n;
          } // end loop over j
        }   // end loop over i
      }     // end loop over time bins
    }       // end loop over wires

    //calculate the cornerness of each pixel while making sure not to fall off the hit map.
    for (unsigned int wire = 1; wire < numberwires - 1; ++wire) {

      for (int timebin = 1; timebin < fTimeBins - 1; ++timebin) {
        MatrixAAsum = 0.;
        MatrixBBsum = 0.;
        MatrixCCsum = 0.;
        //Gaussian smoothing convolution
        n = 0;
        for (int i = -3; i <= 3; ++i) {
          windex = wire + i;
          if (windex < 0) windex = 0;
          // this is ok, because the line before makes sure it's not negative
          else if ((unsigned int)windex >= numberwires)
            windex = numberwires - 1;

          for (int j = -3; j <= 3; ++j) {
            tindex = timebin + j;
            if (tindex < 0)
              tindex = 0;
            else if (tindex >= fTimeBins)
              tindex = fTimeBins - 1;

            MatrixAAsum += w[n] * pow(MatrixAsum[windex][tindex], 2);
            MatrixBBsum += w[n] * pow(MatrixBsum[windex][tindex], 2);
            MatrixCCsum += w[n] * MatrixAsum[windex][tindex] * MatrixBsum[windex][tindex];
            ++n;
          } // end loop over j
        }   // end loop over i

        if ((MatrixAAsum + MatrixBBsum) > 0)
          Cornerness[wire][timebin] =
            (MatrixAAsum * MatrixBBsum - pow(MatrixCCsum, 2)) / (MatrixAAsum + MatrixBBsum);
        else
          Cornerness[wire][timebin] = 0;

        if (Cornerness[wire][timebin] > 0) {
          for (unsigned int i = 0; i < hit.size(); ++i) {
            wire2 = hit[i]->WireID().Wire;
            //make sure the end point candidate coincides with an actual hit.
            if (wire == wire2 && std::abs(hit[i]->TimeDistanceAsRMS(
                                   timebin * (numbertimesamples / fTimeBins))) < 1.) {
              //this index keeps track of the hit number
              hit_loc[wire][timebin] = i;
              Cornerness2.push_back(Cornerness[wire][timebin]);
              break;
            }
          } // end loop over hits
        }   // end if cornerness > 0
      }     // end loop over time bins
    }       // end wire loop

    std::sort(Cornerness2.rbegin(), Cornerness2.rend());

    auto const& plane = wireReadoutGeom.Plane({0, 0, 0});
    for (int vertexnum = 0; vertexnum < fMaxCorners; ++vertexnum) {
      flag = 0;
      for (unsigned int wire = 0; wire < numberwires && flag == 0; ++wire) {
        for (int timebin = 0; timebin < fTimeBins && flag == 0; ++timebin) {
          if (Cornerness2.size() > (unsigned int)vertexnum)
            if (Cornerness[wire][timebin] == Cornerness2[vertexnum] &&
                Cornerness[wire][timebin] > 0. && hit_loc[wire][timebin] > -1) {
              ++flag;

              //thresholding
              if (Cornerness2.size())
                if (Cornerness[wire][timebin] < (fThreshold * Cornerness2[0]))
                  vertexnum = fMaxCorners;
              vHits.push_back(hit[hit_loc[wire][timebin]]);

              // get the total charge from the associated hits
              double totalQ = 0.;
              for (size_t vh = 0; vh < vHits.size(); ++vh)
                totalQ += vHits[vh]->Integral();

              recob::EndPoint2D endpoint(hit[hit_loc[wire][timebin]]->PeakTime(),
                                         hit[hit_loc[wire][timebin]]->WireID(),
                                         Cornerness[wire][timebin],
                                         vtxcol.size(),
                                         view,
                                         totalQ);
              vtxcol.push_back(endpoint);
              vtxHitsOut.push_back(vHits);
              vHits.clear();

              // non-maximal suppression on a square window. The wire coordinate units are
              // converted to time ticks so that the window is truly square.
              // Note that there are 1/0.0743=13.46 time samples per 4.0 mm (wire pitch in ArgoNeuT),
              // assuming a 1.5 mm/us drift velocity for a 500 V/cm E-field

              double drifttick = det_prop.DriftVelocity(det_prop.Efield(), det_prop.Temperature());
              drifttick *= sampling_rate(clock_data) * 1.e-3;
              double wirepitch = plane.WirePitch();
              double corrfactor = drifttick / wirepitch;

              for (int wireout =
                     (int)wire -
                     (int)((fWindow * (numbertimesamples / fTimeBins) * corrfactor) + .5);
                   wireout <=
                   (int)wire + (int)((fWindow * (numbertimesamples / fTimeBins) * corrfactor) + .5);
                   ++wireout) {
                for (int timebinout = timebin - fWindow; timebinout <= timebin + fWindow;
                     timebinout++) {
                  if (std::sqrt(pow(wire - wireout, 2) + pow(timebin - timebinout, 2)) <
                      fWindow) //circular window
                    Cornerness[wireout][timebinout] = 0;
                }
              }
            }
        }
      }
    }
    Cornerness2.clear();
    hit.clear();
    if (clusterIter != clusIn.end()) clusterIter++;

    if ((int)wid.Plane == fSaveVertexMap) {
      unsigned char* outPix = new unsigned char[fTimeBins * numberwires];
      //finds the maximum cell in the map for image scaling
      int cell = 0;
      int pix = 0;
      int maxCell = 0;
      //int xmaxx, ymaxx;
      for (int y = 0; y < fTimeBins; ++y) {
        for (unsigned int x = 0; x < numberwires; ++x) {
          cell = (int)(hit_map[x][y] * 1000);
          if (cell > maxCell) { maxCell = cell; }
        }
      }
      for (int y = 0; y < fTimeBins; ++y) {
        for (unsigned int x = 0; x < numberwires; ++x) {
          //scales the pixel weights based on the maximum cell value
          if (maxCell > 0) pix = (int)((1500000 * hit_map[x][y]) / maxCell);
          outPix[y * numberwires + x] = pix;
        }
      }

      // add 3x3 pixel squares to with the harris vertex finders to the .bmp file

      for (unsigned int ii = 0; ii < vtxcol.size(); ++ii) {
        if (vtxcol[ii].View() == (unsigned int)view) {
          pix = (int)(255);
          outPix[(int)(vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) * numberwires +
                 vtxcol[ii].WireID().Wire] = pix;
          outPix[(int)(vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) * numberwires +
                 vtxcol[ii].WireID().Wire - 1] = pix;
          outPix[(int)(vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) * numberwires +
                 vtxcol[ii].WireID().Wire + 1] = pix;
          outPix[(int)((vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) - 1) *
                   numberwires +
                 vtxcol[ii].WireID().Wire] = pix;
          outPix[(int)((vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) + 1) *
                   numberwires +
                 vtxcol[ii].WireID().Wire] = pix;
          outPix[(int)((vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) - 1) *
                   numberwires +
                 vtxcol[ii].WireID().Wire - 1] = pix;
          outPix[(int)((vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) - 1) *
                   numberwires +
                 vtxcol[ii].WireID().Wire - 2] = pix;
          outPix[(int)((vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) + 1) *
                   numberwires +
                 vtxcol[ii].WireID().Wire + 1] = pix;
          outPix[(int)((vtxcol[ii].DriftTime() * (fTimeBins / numbertimesamples)) + 1) *
                   numberwires +
                 vtxcol[ii].WireID().Wire + 2] = pix;
        }
      }

      VSSaveBMPFile(Form("harrisvertexmap_%d_%d.bmp", (*clusIn.begin())->ID(), wid.Plane),
                    outPix,
                    numberwires,
                    fTimeBins);
      delete[] outPix;
    }
  } // end loop over views

  return vtxcol.size();
}

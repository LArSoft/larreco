/*!
 * Title:   RFFHitFinderAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Class that runs the RFF HitFinder. Implements an RFFHitFitter, and takes
 * the result and stores it in recob::Hit objects.
 *
 * Input:  recob::Wire
 * Output: recob::Hit
*/

#include "fhiclcpp/ParameterSet.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larreco/HitFinder/RFFHitFitter.h"

#include "RFFHitFinderAlg.h"

#include <numeric>

hit::RFFHitFinderAlg::RFFHitFinderAlg(fhicl::ParameterSet const& p)
{
  fMatchThresholdVec = p.get<std::vector<float>>("MeanMatchThreshold");
  fMergeMultiplicityVec = p.get<std::vector<unsigned int>>("MinMergeMultiplicity");
  fAmpThresholdVec = p.get<std::vector<float>>("AmplitudeThreshold", std::vector<float>(1, 0.0));
}

void hit::RFFHitFinderAlg::SetFitterParamsVectors(unsigned int const n_planes)
{
  //If size zero, throw. If size one, assume same for all planes.
  //If size > 1 but < n_planes, throw. If size = n_plane, good.

  if (fMatchThresholdVec.size() == 0 || fMergeMultiplicityVec.size() == 0 ||
      fAmpThresholdVec.size() == 0)
    throw std::runtime_error("Error in RFFHitFinderAlg: Configured with zero planes.");

  if ((fMatchThresholdVec.size() > 1 && fMatchThresholdVec.size() < n_planes) ||
      (fMergeMultiplicityVec.size() > 1 && fMergeMultiplicityVec.size() < n_planes) ||
      (fAmpThresholdVec.size() > 1 && fAmpThresholdVec.size() < n_planes))
    throw std::runtime_error("Error in RFFHitFinderAlg: Configured with incorrect n_planes.");

  if (fMatchThresholdVec.size() == 1) fMatchThresholdVec.resize(n_planes, fMatchThresholdVec[0]);

  if (fMergeMultiplicityVec.size() == 1)
    fMergeMultiplicityVec.resize(n_planes, fMergeMultiplicityVec[0]);

  if (fAmpThresholdVec.size() == 1) fAmpThresholdVec.resize(n_planes, fAmpThresholdVec[0]);
}

void hit::RFFHitFinderAlg::SetFitterParams(unsigned int p)
{
  fFitter.SetFitterParams(fMatchThresholdVec[p], fMergeMultiplicityVec[p], fAmpThresholdVec[p]);
}

void hit::RFFHitFinderAlg::Run(std::vector<recob::Wire> const& wireVector,
                               std::vector<recob::Hit>& hitVector,
                               geo::WireReadoutGeom const& wireReadoutGeom)
{
  hitVector.reserve(wireVector.size());
  for (auto const& wire : wireVector) {
    geo::SigType_t const& sigtype = wireReadoutGeom.SignalType(wire.Channel());
    geo::WireID const& wireID = wireReadoutGeom.ChannelToWire(wire.Channel()).at(0);

    SetFitterParams(wire.View());

    for (auto const& roi : wire.SignalROI().get_ranges()) {
      fFitter.RunFitter(roi.data());

      const float summedADCTotal = std::accumulate(roi.data().begin(), roi.data().end(), 0.0);
      const raw::TDCtick_t startTick = roi.begin_index();
      const raw::TDCtick_t endTick = roi.begin_index() + roi.size();

      EmplaceHit(hitVector, wire, summedADCTotal, startTick, endTick, sigtype, wireID);
    } //end loop over ROIs on wire

  } //end loop over wires
}

void hit::RFFHitFinderAlg::EmplaceHit(std::vector<recob::Hit>& hitVector,
                                      recob::Wire const& wire,
                                      float const& summedADCTotal,
                                      raw::TDCtick_t const& startTick,
                                      raw::TDCtick_t const& endTick,
                                      geo::SigType_t const& sigtype,
                                      geo::WireID const& wireID)
{

  float totalArea = 0.0;
  std::vector<float> areaVector(fFitter.NHits());
  std::vector<float> areaErrorVector(fFitter.NHits());
  std::vector<float> areaFracVector(fFitter.NHits());

  for (size_t ihit = 0; ihit < fFitter.NHits(); ihit++) {
    areaVector[ihit] = fFitter.AmplitudeVector()[ihit] * fFitter.SigmaVector()[ihit] * SQRT_TWO_PI;
    areaErrorVector[ihit] =
      SQRT_TWO_PI * std::sqrt(fFitter.AmplitudeVector()[ihit] * fFitter.SigmaErrorVector()[ihit] *
                                fFitter.AmplitudeVector()[ihit] * fFitter.SigmaErrorVector()[ihit] +
                              fFitter.AmplitudeErrorVector()[ihit] * fFitter.SigmaVector()[ihit] *
                                fFitter.AmplitudeErrorVector()[ihit] * fFitter.SigmaVector()[ihit]);
    totalArea += areaVector[ihit];
  }

  for (size_t ihit = 0; ihit < fFitter.NHits(); ihit++) {
    areaFracVector[ihit] = areaVector[ihit] / totalArea;

    hitVector.emplace_back(wire.Channel(),
                           startTick,
                           endTick,
                           fFitter.MeanVector()[ihit] + (float)startTick,
                           fFitter.MeanErrorVector()[ihit],
                           fFitter.SigmaVector()[ihit],
                           fFitter.AmplitudeVector()[ihit],
                           fFitter.AmplitudeErrorVector()[ihit],
                           summedADCTotal * areaFracVector[ihit],
                           areaVector[ihit],
                           areaErrorVector[ihit],
                           fFitter.NHits(),
                           ihit,
                           -999.,
                           -999,
                           wire.View(),
                           sigtype,
                           wireID);
  }
}

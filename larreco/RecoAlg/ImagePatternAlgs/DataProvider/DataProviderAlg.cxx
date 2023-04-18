////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan, May 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ImagePatternAlgs/DataProvider/DataProviderAlg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "lardataobj/RecoBase/Wire.h"
#include "CLHEP/Random/RandGauss.h"

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

namespace detinfo {
  class DetectorClocksData;
}

img::DataProviderAlg::DataProviderAlg(const Config& config)
  : fAlgView{}
  , fDownscaleMode(img::DataProviderAlg::kMax)
  , fDriftWindow(10)
  , fCalorimetryAlg(config.CalorimetryAlg())
  , fGeometry(art::ServiceHandle<geo::Geometry const>().get())
  , fAdcSumOverThr(0)
  , fAdcSumThr(10)
  , // set fixed threshold of 10 ADC counts for counting the sum
  fAdcAreaOverThr(0)
  , fNoiseSigma(0)
  , fCoherentSigma(0)
{
  fCalibrateLifetime = config.CalibrateLifetime();
  fCalibrateAmpl = config.CalibrateAmpl();

  fAmplCalibConst.resize(fGeometry->MaxPlanes());
  if (fCalibrateAmpl) {
    mf::LogInfo("DataProviderAlg") << "Using calibration constants:";
    for (size_t p = 0; p < fAmplCalibConst.size(); ++p) {
      try {
        fAmplCalibConst[p] = 1.2e-3 * fCalorimetryAlg.ElectronsFromADCPeak(1.0, p);
        mf::LogInfo("DataProviderAlg") << "   plane:" << p << " const:" << 1.0 / fAmplCalibConst[p];
      }
      catch (...) {
        fAmplCalibConst[p] = 1.0;
      }
    }
  }
  else {
    mf::LogInfo("DataProviderAlg") << "No plane-to-plane calibration.";
    for (size_t p = 0; p < fAmplCalibConst.size(); ++p) {
      fAmplCalibConst[p] = 1.0;
    }
  }

  fDriftWindow = config.DriftWindow();
  fDownscaleFullView = config.DownscaleFullView();
  fDriftWindowInv = 1.0 / fDriftWindow;

  std::string mode_str = config.DownscaleFn();
  mf::LogVerbatim("DataProviderAlg") << "Downscale mode is: " << mode_str;
  if (mode_str == "maxpool") {
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
    fDownscaleMode = img::DataProviderAlg::kMax;
  }
  else if (mode_str == "maxmean") {
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMaxMean(dst, adc, tick0); };
    fDownscaleMode = img::DataProviderAlg::kMaxMean;
  }
  else if (mode_str == "mean") {
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMean(dst, adc, tick0); };
    fDownscaleMode = img::DataProviderAlg::kMean;
  }
  else {
    mf::LogError("DataProviderAlg") << "Downscale mode string not recognized, set to max pooling.";
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
    fDownscaleMode = img::DataProviderAlg::kMax;
  }

  fAdcMax = config.AdcMax();
  fAdcMin = config.AdcMin();
  fAdcOffset = config.OutMin();
  fAdcScale = (config.OutMax() - fAdcOffset) / (fAdcMax - fAdcMin);
  fAdcZero = fAdcOffset + fAdcScale * (0 - fAdcMin); // level of zero ADC after scaling

  if (fAdcMax <= fAdcMin) {
    throw cet::exception("img::DataProviderAlg") << "Misconfigured: AdcMax <= AdcMin" << std::endl;
  }
  if (fAdcScale == 0) {
    throw cet::exception("img::DataProviderAlg") << "Misconfigured: OutMax == OutMin" << std::endl;
  }

  fBlurKernel = config.BlurKernel();
  fNoiseSigma = config.NoiseSigma();
  fCoherentSigma = config.CoherentSigma();
}
// ------------------------------------------------------

img::DataProviderAlg::~DataProviderAlg() = default;
// ------------------------------------------------------

img::DataProviderAlgView img::DataProviderAlg::resizeView(
  detinfo::DetectorClocksData const& clock_data,
  detinfo::DetectorPropertiesData const& det_prop,
  size_t wires,
  size_t drifts)
{
  img::DataProviderAlgView result;
  result.fNWires = wires;
  result.fNDrifts = drifts;
  result.fNScaledDrifts = drifts / fDriftWindow;
  result.fNCachedDrifts = fDownscaleFullView ? result.fNScaledDrifts : drifts;

  result.fWireChannels.resize(wires, raw::InvalidChannelID);

  result.fWireDriftData.resize(wires, std::vector<float>(result.fNCachedDrifts, fAdcZero));

  result.fLifetimeCorrFactors.resize(drifts);
  if (fCalibrateLifetime) {
    for (size_t t = 0; t < drifts; ++t) {
      result.fLifetimeCorrFactors[t] = fCalorimetryAlg.LifetimeCorrection(clock_data, det_prop, t);
    }
  }
  else {
    for (size_t t = 0; t < drifts; ++t) {
      result.fLifetimeCorrFactors[t] = 1.0;
    }
  }
  return result;
}
// ------------------------------------------------------

float img::DataProviderAlg::poolMax(int wire, int drift, size_t r) const
{
  size_t rw = r, rd = r;
  if (!fDownscaleFullView) { rd *= fDriftWindow; }

  size_t didx = getDriftIndex(drift);
  int d0 = didx - rd;
  if (d0 < 0) { d0 = 0; }
  int d1 = didx + rd;
  if (d1 >= (int)fAlgView.fNCachedDrifts) { d1 = fAlgView.fNCachedDrifts - 1; }

  int w0 = wire - rw;
  if (w0 < 0) { w0 = 0; }
  int w1 = wire + rw;
  if (w1 >= (int)fAlgView.fNWires) { w1 = fAlgView.fNWires - 1; }

  float adc, max_adc = 0;
  for (int w = w0; w <= w1; ++w) {
    auto const* col = fAlgView.fWireDriftData[w].data();
    for (int d = d0; d <= d1; ++d) {
      adc = col[d];
      if (adc > max_adc) { max_adc = adc; }
    }
  }

  return max_adc;
}
// ------------------------------------------------------

//float img::DataProviderAlg::poolSum(int wire, int drift, size_t r) const
//{
//    size_t rw = r, rd = r;
//    if (!fDownscaleFullView) { rd *= fDriftWindow; }
//
//    size_t didx = getDriftIndex(drift);
//    int d0 = didx - rd; if (d0 < 0) { d0 = 0; }
//    int d1 = didx + rd; if (d1 >= (int)fNCachedDrifts) { d1 = fNCachedDrifts - 1; }
//
//    int w0 = wire - rw; if (w0 < 0) { w0 = 0; }
//    int w1 = wire + rw; if (w1 >= (int)fNWires) { w1 = fNWires - 1; }
//
//    float sum = 0;
//    for (int w = w0; w <= w1; ++w)
//    {
//        auto const * col = fWireDriftData[w].data();
//        for (int d = d0; d <= d1; ++d) { sum += col[d]; }
//    }
//
//    return sum;
//}
// ------------------------------------------------------
std::vector<float> img::DataProviderAlg::downscaleMax(std::size_t dst_size,
                                                      std::vector<float> const& adc,
                                                      size_t tick0) const
{
  size_t kStop = dst_size;
  std::vector<float> result(dst_size);
  if (adc.size() < kStop) { kStop = adc.size(); }
  for (size_t i = 0, k0 = 0; i < kStop; ++i, k0 += fDriftWindow) {
    size_t k1 = k0 + fDriftWindow;

    float max_adc = adc[k0] * fAlgView.fLifetimeCorrFactors[k0 + tick0];
    for (size_t k = k0 + 1; k < k1; ++k) {
      float ak = adc[k] * fAlgView.fLifetimeCorrFactors[k + tick0];
      if (ak > max_adc) max_adc = ak;
    }
    result[i] = max_adc;
  }
  scaleAdcSamples(result);
  return result;
}

std::vector<float> img::DataProviderAlg::downscaleMaxMean(std::size_t dst_size,
                                                          std::vector<float> const& adc,
                                                          size_t tick0) const
{
  size_t kStop = dst_size;
  std::vector<float> result(dst_size);
  if (adc.size() < kStop) { kStop = adc.size(); }
  for (size_t i = 0, k0 = 0; i < kStop; ++i, k0 += fDriftWindow) {
    size_t k1 = k0 + fDriftWindow;
    size_t max_idx = k0;
    float max_adc = adc[k0] * fAlgView.fLifetimeCorrFactors[k0 + tick0];
    for (size_t k = k0 + 1; k < k1; ++k) {
      float ak = adc[k] * fAlgView.fLifetimeCorrFactors[k + tick0];
      if (ak > max_adc) {
        max_adc = ak;
        max_idx = k;
      }
    }

    size_t n = 1;
    if (max_idx > 0) {
      max_adc += adc[max_idx - 1] * fAlgView.fLifetimeCorrFactors[max_idx - 1 + tick0];
      n++;
    }
    if (max_idx + 1 < adc.size()) {
      max_adc += adc[max_idx + 1] * fAlgView.fLifetimeCorrFactors[max_idx + 1 + tick0];
      n++;
    }

    result[i] = max_adc / n;
  }
  scaleAdcSamples(result);
  return result;
}

std::vector<float> img::DataProviderAlg::downscaleMean(std::size_t dst_size,
                                                       std::vector<float> const& adc,
                                                       size_t tick0) const
{
  size_t kStop = dst_size;
  std::vector<float> result(dst_size);
  if (adc.size() < kStop) { kStop = adc.size(); }
  for (size_t i = 0, k0 = 0; i < kStop; ++i, k0 += fDriftWindow) {
    size_t k1 = k0 + fDriftWindow;

    float sum_adc = 0;
    for (size_t k = k0; k < k1; ++k) {
      if (k + tick0 < fAlgView.fLifetimeCorrFactors.size())
        sum_adc += adc[k] * fAlgView.fLifetimeCorrFactors[k + tick0];
    }
    result[i] = sum_adc * fDriftWindowInv;
  }
  scaleAdcSamples(result);
  return result;
}

std::optional<std::vector<float>> img::DataProviderAlg::setWireData(std::vector<float> const& adc,
                                                                    size_t wireIdx) const
{
  if (wireIdx >= fAlgView.fWireDriftData.size()) { return std::nullopt; }
  auto& wData = fAlgView.fWireDriftData[wireIdx];

  if (fDownscaleFullView) {
    if (!adc.empty()) { return downscale(wData.size(), adc, 0); }
    return std::nullopt;
  }
  if (adc.empty()) { return std::nullopt; }
  if (adc.size() <= wData.size()) { return adc; }
  return std::vector<float>(adc.begin(), adc.begin() + wData.size());
}
// ------------------------------------------------------

bool img::DataProviderAlg::setWireDriftData(detinfo::DetectorClocksData const& clock_data,
                                            detinfo::DetectorPropertiesData const& det_prop,
                                            const std::vector<recob::Wire>& wires,
                                            unsigned int plane,
                                            unsigned int tpc,
                                            unsigned int cryo,
                                            art::Timestamp t)
{
  mf::LogInfo("DataProviderAlg") << "Create image for cryo:" << cryo << " tpc:" << tpc
                                 << " plane:" << plane;

  fCryo = cryo;
  fTPC = tpc;
  fPlane = plane;

  fAdcSumOverThr = 0;
  fAdcAreaOverThr = 0;

  size_t nwires = fGeometry->Nwires({cryo, tpc, plane});
  size_t ndrifts = det_prop.NumberTimeSamples();

  fAlgView = resizeView(clock_data, det_prop, nwires, ndrifts);

  auto const& channelStatus =
    art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

  bool allWrong = true;
  for (auto const& wire : wires) {
    auto wireChannelNumber = wire.Channel();
    if (!channelStatus.IsGood(t.value(), wireChannelNumber)) { continue; }

    size_t w_idx = 0;
    for (auto const& id : fGeometry->ChannelToWire(wireChannelNumber)) {
      if ((id.Plane == plane) && (id.TPC == tpc) && (id.Cryostat == cryo)) {
        w_idx = id.Wire;

        auto adc = wire.Signal();
        if (adc.size() < ndrifts) {
          mf::LogWarning("DataProviderAlg") << "Wire ADC vector size lower than NumberTimeSamples.";
          continue; // not critical, maybe other wires are OK, so continue
        }
        auto wire_data = setWireData(adc, w_idx);
        if (!wire_data) {
          mf::LogWarning("DataProviderAlg") << "Wire data not set.";
          continue; // also not critical, try to set other wires
        }
        fAlgView.fWireDriftData[w_idx] = *wire_data;
        for (auto v : adc) {
          if (v >= fAdcSumThr) {
            fAdcSumOverThr += v;
            fAdcAreaOverThr++;
          }
        }

        fAlgView.fWireChannels[w_idx] = wireChannelNumber;
        allWrong = false;
      }
    }
  }
  if (allWrong) {
    mf::LogError("DataProviderAlg")
      << "Wires data not set in the cryo:" << cryo << " tpc:" << tpc << " plane:" << plane;
    return false;
  }

  applyBlur();
  addWhiteNoise();
  addCoherentNoise();

  return true;
}
// ------------------------------------------------------

float img::DataProviderAlg::scaleAdcSample(float val) const
{
  val *= fAmplCalibConst[fPlane]; // prescale by plane-to-plane calibration factors

  if (val < fAdcMin) { val = fAdcMin; } // saturate
  else if (val > fAdcMax) {
    val = fAdcMax;
  }

  return fAdcOffset +
         fAdcScale *
           (val - fAdcMin); // shift and scale to the output range, shift to the output min
}
// ------------------------------------------------------
void img::DataProviderAlg::scaleAdcSamples(std::vector<float>& values) const
{
  float calib = fAmplCalibConst[fPlane];
  auto* data = values.data();

  size_t k = 0, size4 = values.size() >> 2, size = values.size();
  for (size_t i = 0; i < size4; ++i) // vectorize if you can
  {
    data[k] *= calib; // prescale by plane-to-plane calibration factors
    data[k + 1] *= calib;
    data[k + 2] *= calib;
    data[k + 3] *= calib;

    if (data[k] < fAdcMin) { data[k] = fAdcMin; } // saturate min
    if (data[k + 1] < fAdcMin) { data[k + 1] = fAdcMin; }
    if (data[k + 2] < fAdcMin) { data[k + 2] = fAdcMin; }
    if (data[k + 3] < fAdcMin) { data[k + 3] = fAdcMin; }

    if (data[k] > fAdcMax) { data[k] = fAdcMax; } // saturate max
    if (data[k + 1] > fAdcMax) { data[k + 1] = fAdcMax; }
    if (data[k + 2] > fAdcMax) { data[k + 2] = fAdcMax; }
    if (data[k + 3] > fAdcMax) { data[k + 3] = fAdcMax; }

    data[k] = fAdcOffset +
              fAdcScale *
                (data[k] - fAdcMin); // shift and scale to the output range, shift to the output min
    data[k + 1] = fAdcOffset + fAdcScale * (data[k + 1] - fAdcMin);
    data[k + 2] = fAdcOffset + fAdcScale * (data[k + 2] - fAdcMin);
    data[k + 3] = fAdcOffset + fAdcScale * (data[k + 3] - fAdcMin);

    k += 4;
  }
  while (k < size) {
    data[k] = scaleAdcSample(data[k]);
    ++k;
  } // do the tail
}
// ------------------------------------------------------

void img::DataProviderAlg::applyBlur()
{
  if (fBlurKernel.size() < 2) return;

  size_t margin_left = (fBlurKernel.size() - 1) >> 1,
         margin_right = fBlurKernel.size() - margin_left - 1;

  std::vector<std::vector<float>> src(fAlgView.fWireDriftData.size());
  for (size_t w = 0; w < fAlgView.fWireDriftData.size(); ++w) {
    src[w] = fAlgView.fWireDriftData[w];
  }

  for (size_t w = margin_left; w < fAlgView.fWireDriftData.size() - margin_right; ++w) {
    for (size_t d = 0; d < fAlgView.fWireDriftData[w].size(); ++d) {
      float sum = 0;
      for (size_t i = 0; i < fBlurKernel.size(); ++i) {
        sum += fBlurKernel[i] * src[w + i - margin_left][d];
      }
      fAlgView.fWireDriftData[w][d] = sum;
    }
  }
}
// ------------------------------------------------------

// MUST give the same result as get_patch() in scripts/utils.py
bool img::DataProviderAlg::patchFromDownsampledView(size_t wire,
                                                    float drift,
                                                    size_t size_w,
                                                    size_t size_d,
                                                    std::vector<std::vector<float>>& patch) const
{
  int halfSizeW = size_w / 2;
  int halfSizeD = size_d / 2;

  int w0 = wire - halfSizeW;
  int w1 = wire + halfSizeW;

  size_t sd = (size_t)(drift / fDriftWindow);
  int d0 = sd - halfSizeD;
  int d1 = sd + halfSizeD;

  int wsize = fAlgView.fWireDriftData.size();
  for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch) {
    auto& dst = patch[wpatch];
    if ((w >= 0) && (w < wsize)) {
      auto& src = fAlgView.fWireDriftData[w];
      int dsize = src.size();
      for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch) {
        if ((d >= 0) && (d < dsize)) { dst[dpatch] = src[d]; }
        else {
          dst[dpatch] = fAdcZero;
        }
      }
    }
    else {
      std::fill(dst.begin(), dst.end(), fAdcZero);
    }
  }

  return true;
}

bool img::DataProviderAlg::patchFromOriginalView(size_t wire,
                                                 float drift,
                                                 size_t size_w,
                                                 size_t size_d,
                                                 std::vector<std::vector<float>>& patch) const
{
  int dsize = fDriftWindow * size_d;
  int halfSizeW = size_w / 2;
  int halfSizeD = dsize / 2;

  int w0 = wire - halfSizeW;
  int w1 = wire + halfSizeW;

  int d0 = int(drift) - halfSizeD;
  int d1 = int(drift) + halfSizeD;

  if (d0 < 0) d0 = 0;

  std::vector<float> tmp(dsize);
  int wsize = fAlgView.fWireDriftData.size();
  for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch) {
    if ((w >= 0) && (w < wsize)) {
      auto& src = fAlgView.fWireDriftData[w];
      int src_size = src.size();
      for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch) {
        if ((d >= 0) && (d < src_size)) { tmp[dpatch] = src[d]; }
        else {
          tmp[dpatch] = fAdcZero;
        }
      }
    }
    else {
      std::fill(tmp.begin(), tmp.end(), fAdcZero);
    }
    patch[wpatch] = downscale(patch[wpatch].size(), tmp, d0);
  }

  return true;
}
// ------------------------------------------------------

void img::DataProviderAlg::addWhiteNoise()
{
  if (fNoiseSigma == 0) return;

  double effectiveSigma = scaleAdcSample(fNoiseSigma);
  if (fDownscaleFullView) effectiveSigma /= fDriftWindow;

  CLHEP::RandGauss gauss(fRndEngine);
  std::vector<double> noise(fAlgView.fNCachedDrifts);
  for (auto& wire : fAlgView.fWireDriftData) {
    gauss.fireArray(fAlgView.fNCachedDrifts, noise.data(), 0., effectiveSigma);
    for (size_t d = 0; d < wire.size(); ++d) {
      wire[d] += noise[d];
    }
  }
}
// ------------------------------------------------------

void img::DataProviderAlg::addCoherentNoise()
{
  if (fCoherentSigma == 0) return;

  double effectiveSigma = scaleAdcSample(fCoherentSigma);
  if (fDownscaleFullView) effectiveSigma /= fDriftWindow;

  CLHEP::RandGauss gauss(fRndEngine);
  std::vector<double> amps1(fAlgView.fWireDriftData.size());
  std::vector<double> amps2(1 + (fAlgView.fWireDriftData.size() / 32));
  gauss.fireArray(amps1.size(), amps1.data(), 1., 0.1); // 10% wire-wire ampl. variation
  gauss.fireArray(amps2.size(), amps2.data(), 1., 0.1); // 10% group-group ampl. variation

  double group_amp = 1.0;
  std::vector<double> noise(fAlgView.fNCachedDrifts);
  for (size_t w = 0; w < fAlgView.fWireDriftData.size(); ++w) {
    if ((w & 31) == 0) {
      group_amp = amps2[w >> 5]; // div by 32
      gauss.fireArray(fAlgView.fNCachedDrifts, noise.data(), 0., effectiveSigma);
    } // every 32 wires

    auto& wire = fAlgView.fWireDriftData[w];
    for (size_t d = 0; d < wire.size(); ++d) {
      wire[d] += group_amp * amps1[w] * noise[d];
    }
  }
}
// ------------------------------------------------------

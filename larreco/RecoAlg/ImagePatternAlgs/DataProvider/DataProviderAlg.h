////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Authors:     D.Stefan (Dorota.Stefan@ncbj.gov.pl),         from DUNE, CERN/NCBJ, since May 2016
//              R.Sulej (Robert.Sulej@cern.ch),               from DUNE, FNAL/NCBJ, since May 2016
//              P.Plonski,                                    from DUNE, WUT,       since May 2016
//
//
// Algorithm for making 2D image-like data from recob::Wire's. Used by CNN codes for training data
// preparation and application of trained models to new data. Also used by PMA to keep images of
// 2D projections used for the track validation.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef DataProviderAlg_h
#define DataProviderAlg_h

// Framework includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// LArSoft includes
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "CLHEP/Random/JamesRandom.h" // for testing on noise, not used by any reco

// ROOT & C++
#include <memory>

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace img {
  class DataProviderAlg;
  struct DataProviderAlgView {
    unsigned int fNWires;
    unsigned int fNDrifts;
    unsigned int fNScaledDrifts;
    unsigned int fNCachedDrifts;
    std::vector<raw::ChannelID_t> fWireChannels;
    std::vector<std::vector<float>> fWireDriftData;
    std::vector<float> fLifetimeCorrFactors;
  };
}

/// Base class providing data for training / running image based classifiers. It can be used
/// also for any other algorithms where 2D projection image is useful. Currently the image
/// is 32-bit fp / pixel, as sson as have time will template it so e.g. byte pixels would
/// be possible.
class img::DataProviderAlg {
public:
  enum EDownscaleMode { kMax = 1, kMaxMean = 2, kMean = 3 };

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Table<calo::CalorimetryAlg::Config> CalorimetryAlg{
      Name("CalorimetryAlg"),
      Comment("Used to eliminate amplitude variation due to electron lifetime.")};

    fhicl::Atom<float> AdcMax{Name("AdcMax"), Comment("Saturation max value")};
    fhicl::Atom<float> AdcMin{Name("AdcMin"), Comment("Saturation min value")};
    fhicl::Atom<float> OutMax{Name("OutMax"), Comment("Output max value")};
    fhicl::Atom<float> OutMin{Name("OutMin"), Comment("Output min value")};

    fhicl::Atom<bool> CalibrateAmpl{Name("CalibrateAmpl"),
                                    Comment("Calibrate ADC values with CalAmpConstants")};

    fhicl::Atom<bool> CalibrateLifetime{Name("CalibrateLifetime"),
                                        Comment("Calibrate ADC values with the electron lifetime")};

    fhicl::Atom<unsigned int> DriftWindow{Name("DriftWindow"),
                                          Comment("Downsampling window (in drift ticks).")};

    fhicl::Atom<std::string> DownscaleFn{Name("DownscaleFn"), Comment("Downsampling function")};

    fhicl::Atom<bool> DownscaleFullView{
      Name("DownscaleFullView"),
      Comment("Downsample full view (faster / lower location precision)")};

    fhicl::Sequence<float> BlurKernel{Name("BlurKernel"), Comment("Blur kernel in wire direction")};

    fhicl::Atom<float> NoiseSigma{Name("NoiseSigma"), Comment("White noise sigma")};

    fhicl::Atom<float> CoherentSigma{Name("CoherentSigma"), Comment("Coherent noise sigma")};
  };

  DataProviderAlg(const fhicl::ParameterSet& pset)
    : DataProviderAlg(fhicl::Table<Config>(pset, {})())
  {}

  DataProviderAlg(const Config& config);

  virtual ~DataProviderAlg();


  bool setWireDriftData(const detinfo::DetectorClocksData& clock_data,
                        const detinfo::DetectorPropertiesData& det_prop,
                        const std::vector<recob::Wire>&
                          wires, // once per plane: setup ADC buffer, collect & downscale ADC's
                        unsigned int plane,
                        unsigned int tpc,
                        unsigned int cryo);

  std::vector<float> const&
  wireData(size_t widx) const
  {
    return fAlgView.fWireDriftData[widx];
  }

  /// Return patch of data centered on the wire and drift, witht the size in (downscaled) pixels givent
  /// with patchSizeW and patchSizeD.  Pad with the zero-level calue if patch extends beyond the event
  /// projection.
  std::vector<std::vector<float>>
  getPatch(size_t wire, float drift, size_t patchSizeW, size_t patchSizeD) const
  {
    bool ok = false;
    std::vector<std::vector<float>> patch;
    if (fDownscaleFullView) {
      ok = patchFromDownsampledView(wire, drift, patchSizeW, patchSizeD, patch);
    }
    else {
      ok = patchFromOriginalView(wire, drift, patchSizeW, patchSizeD, patch);
    }

    if (ok)
      return patch;
    throw cet::exception("img::DataProviderAlg") << "Patch filling failed." << std::endl;
  }

  /// Return value from the ADC buffer, or zero if coordinates are out of the view;
  /// will scale the drift according to the downscale settings.
  float
  getPixelOrZero(int wire, int drift) const
  {
    size_t didx = getDriftIndex(drift), widx = (size_t)wire;

    if ((widx < fAlgView.fWireDriftData.size()) && (didx < fAlgView.fNCachedDrifts)) {
      return fAlgView.fWireDriftData[widx][didx];
    }
      return 0;
  }

  double
  getAdcSum() const
  {
    return fAdcSumOverThr;
  }
  size_t
  getAdcArea() const
  {
    return fAdcAreaOverThr;
  }

  /// Pool max value in a patch around the wire/drift pixel.
  float poolMax(int wire, int drift, size_t r = 0) const;

  /// Pool sum of pixels in a patch around the wire/drift pixel.
  //float poolSum(int wire, int drift, size_t r = 0) const;

  unsigned int
  Cryo() const
  {
    return fCryo;
  }
  unsigned int
  TPC() const
  {
    return fTPC;
  }
  unsigned int
  Plane() const
  {
    return fPlane;
  }

  unsigned int
  NWires() const
  {
    return fAlgView.fNWires;
  }
  unsigned int
  NScaledDrifts() const
  {
    return fAlgView.fNScaledDrifts;
  }
  unsigned int
  NCachedDrifts() const
  {
    return fAlgView.fNCachedDrifts;
  }
  unsigned int
  DriftWindow() const
  {
    return fDriftWindow;
  }

  /// Level of zero ADC after scaling.
  float
  ZeroLevel() const
  {
    return fAdcZero;
  }

  double
  LifetimeCorrection(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     double tick) const
  {
    return fCalorimetryAlg.LifetimeCorrection(clock_data, det_prop, tick);
  }

protected:
  DataProviderAlgView fAlgView;
  EDownscaleMode fDownscaleMode;
  //std::function<void (std::vector<float> &, std::vector<float> const &, size_t)> fnDownscale;

  size_t fDriftWindow;
  bool fDownscaleFullView;
  float fDriftWindowInv;

  std::vector<float> downscaleMax(std::size_t dst_size,
                                  std::vector<float> const& adc,
                                  size_t tick0) const;
  std::vector<float> downscaleMaxMean(std::size_t dst_size,
                                      std::vector<float> const& adc,
                                      size_t tick0) const;
  std::vector<float> downscaleMean(std::size_t dst_size,
                                   std::vector<float> const& adc,
                                   size_t tick0) const;
  std::vector<float>
  downscale(std::size_t dst_size, std::vector<float> const& adc, size_t tick0) const
  {
    switch (fDownscaleMode) {
    case img::DataProviderAlg::kMean: return downscaleMean(dst_size, adc, tick0);
    case img::DataProviderAlg::kMaxMean: return downscaleMaxMean(dst_size, adc, tick0);
    case img::DataProviderAlg::kMax: return downscaleMax(dst_size, adc, tick0);
    }
    throw cet::exception("img::DataProviderAlg") << "Downscale mode not supported." << std::endl;
  }

  size_t
  getDriftIndex(float drift) const
  {
    if (fDownscaleFullView)
      return (size_t)(drift * fDriftWindowInv);
    else
      return (size_t)drift;
  }

  std::optional<std::vector<float>> setWireData(std::vector<float> const& adc,
                                                size_t wireIdx) const;

  bool patchFromDownsampledView(size_t wire,
                                float drift,
                                size_t size_w,
                                size_t size_d,
                                std::vector<std::vector<float>>& patch) const;
  bool patchFromOriginalView(size_t wire,
                             float drift,
                             size_t size_w,
                             size_t size_d,
                             std::vector<std::vector<float>>& patch) const;

  virtual DataProviderAlgView resizeView(detinfo::DetectorClocksData const& clock_data,
                          detinfo::DetectorPropertiesData const& det_prop,
                          size_t wires,
                          size_t drifts);

  // Calorimetry needed to equalize ADC amplitude along drift:
  calo::CalorimetryAlg fCalorimetryAlg;

  // Geometry and detector properties:
  geo::GeometryCore const* fGeometry;

private:
  float scaleAdcSample(float val) const;
  void scaleAdcSamples(std::vector<float>& values) const;
  std::vector<float> fAmplCalibConst;
  bool fCalibrateAmpl, fCalibrateLifetime;
  unsigned int fCryo = 9999, fTPC = 9999, fPlane = 9999;
  float fAdcMax, fAdcMin, fAdcScale, fAdcOffset, fAdcZero;
  double fAdcSumOverThr, fAdcSumThr;
  size_t fAdcAreaOverThr;

  CLHEP::HepJamesRandom fRndEngine;

  void applyBlur();
  std::vector<float> fBlurKernel; // blur not applied if empty

  void addWhiteNoise();
  float fNoiseSigma; // noise not added if sigma=0

  void addCoherentNoise();
  float fCoherentSigma; // noise not added if sigma=0
};
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

#endif

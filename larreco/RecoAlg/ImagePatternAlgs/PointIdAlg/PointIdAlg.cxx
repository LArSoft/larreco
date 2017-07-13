////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan, May 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/Geometry/ChannelMapAlg.h" // geo::InvalidWireIDError
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGauss.h"

#include <sys/stat.h>

nnet::DataProviderAlg::DataProviderAlg(const Config& config) :
	fCryo(9999), fTPC(9999), fView(9999),
	fNWires(0), fNDrifts(0), fNScaledDrifts(0), fNCachedDrifts(0),
	fDownscaleMode(nnet::DataProviderAlg::kMax), fDriftWindow(10),
	fCalorimetryAlg(config.CalorimetryAlg()),
	fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
	fNoiseSigma(0), fCoherentSigma(0)
{
	fGeometry = &*(art::ServiceHandle<geo::Geometry>());

	this->reconfigure(config); 
}
// ------------------------------------------------------

nnet::DataProviderAlg::~DataProviderAlg(void)
{
}
// ------------------------------------------------------

void nnet::DataProviderAlg::reconfigure(const Config& config)
{
	fCalorimetryAlg.reconfigure(config.CalorimetryAlg());
	fCalibrateAmpl = config.CalibrateAmpl();
	if (fCalibrateAmpl)
	{
	    fAmplCalibConst.resize(fGeometry->MaxPlanes());
	    mf::LogInfo("DataProviderAlg") << "Using calibration constants:";
	    for (size_t p = 0; p < fAmplCalibConst.size(); ++p)
	    {
	        try
	        {
	            fAmplCalibConst[p] = 1.2e-3 * fCalorimetryAlg.ElectronsFromADCPeak(1.0, p);
    	        mf::LogInfo("DataProviderAlg") << "   plane:" << p << " const:" << 1.0 / fAmplCalibConst[p];
    	    }
    	    catch (...) { fAmplCalibConst[p] = 1.0; }
	    }
	}

	fDriftWindow = config.DriftWindow();
	fDownscaleFullView = config.DownscaleFullView();
	fDriftWindowInv = 1.0 / fDriftWindow;

	std::string mode_str = config.DownscaleFn();
	if (mode_str == "maxpool")      fDownscaleMode = nnet::DataProviderAlg::kMax;
	else if (mode_str == "maxmean") fDownscaleMode = nnet::DataProviderAlg::kMaxMean;
	else if (mode_str == "mean")    fDownscaleMode = nnet::DataProviderAlg::kMean;
	else
	{
		mf::LogError("DataProviderAlg") << "Downscale mode string not recognized, set to max pooling.";
		fDownscaleMode = nnet::DataProviderAlg::kMax;
	}

    fBlurKernel = config.BlurKernel();
    fNoiseSigma = config.NoiseSigma();
    fCoherentSigma = config.CoherentSigma();
}
// ------------------------------------------------------

void nnet::DataProviderAlg::resizeView(size_t wires, size_t drifts)
{
    fNWires = wires; fNDrifts = drifts;
    fNScaledDrifts = drifts / fDriftWindow;

    if (fDownscaleFullView) { fNCachedDrifts = fNScaledDrifts; }
    else { fNCachedDrifts = fNDrifts; }

    fWireChannels.resize(wires);
    std::fill(fWireChannels.begin(), fWireChannels.end(), raw::InvalidChannelID);

    fWireDriftData.resize(wires);
    for (auto & w : fWireDriftData)
    {
    	w.resize(fNCachedDrifts);
    	std::fill(w.begin(), w.end(), 0.0F);
    }

    fLifetimeCorrFactors.resize(fNDrifts);
    for (size_t t = 0; t < fNDrifts; ++t)
    {
        fLifetimeCorrFactors[t] = fCalorimetryAlg.LifetimeCorrection(t);
    }
}
// ------------------------------------------------------

float nnet::DataProviderAlg::poolMax(int wire, int drift, size_t r) const
{
    size_t rw = r, rd = r;
    if (!fDownscaleFullView) { rd *= fDriftWindow; }

    size_t didx = getDriftIndex(drift);
    int d0 = didx - rd; if (d0 < 0) { d0 = 0; }
    int d1 = didx + rd; if (d1 >= (int)fNCachedDrifts) { d1 = fNCachedDrifts - 1; }

    int w0 = wire - rw; if (w0 < 0) { w0 = 0; }
    int w1 = wire + rw; if (w1 >= (int)fNWires) { w1 = fNWires - 1; }

    float adc, max_adc = 0;
    for (int w = w0; w <= w1; ++w)
    {
        auto const * col = fWireDriftData[w].data();
        for (int d = d0; d <= d1; ++d)
        {
            adc = col[d]; if (adc > max_adc) { max_adc = adc; }
        }
    }

    return max_adc;
}
// ------------------------------------------------------

float nnet::DataProviderAlg::poolSum(int wire, int drift, size_t r) const
{
    size_t rw = r, rd = r;
    if (!fDownscaleFullView) { rd *= fDriftWindow; }

    size_t didx = getDriftIndex(drift);
    int d0 = didx - rd; if (d0 < 0) { d0 = 0; }
    int d1 = didx + rd; if (d1 >= (int)fNCachedDrifts) { d1 = fNCachedDrifts - 1; }

    int w0 = wire - rw; if (w0 < 0) { w0 = 0; }
    int w1 = wire + rw; if (w1 >= (int)fNWires) { w1 = fNWires - 1; }

    float sum = 0;
    for (int w = w0; w <= w1; ++w)
    {
        auto const * col = fWireDriftData[w].data();
        for (int d = d0; d <= d1; ++d) { sum += col[d]; }
    }

    return sum;
}
// ------------------------------------------------------

void nnet::DataProviderAlg::downscaleMax(std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) const
{
	for (size_t i = 0; i < dst.size(); ++i)
	{
		size_t k0 = i * fDriftWindow;
		size_t k1 = (i + 1) * fDriftWindow;

		float max_adc = adc[k0] * fLifetimeCorrFactors[k0 + tick0];
		for (size_t k = k0 + 1; k < k1; ++k)
		{
			float ak = adc[k] * fLifetimeCorrFactors[k + tick0];
			if (ak > max_adc) max_adc = ak;
		}

		dst[i] = scaleAdcSample(max_adc);
	}
}

void nnet::DataProviderAlg::downscaleMaxMean(std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) const
{
	for (size_t i = 0; i < dst.size(); ++i)
	{
		size_t k0 = i * fDriftWindow;
		size_t k1 = (i + 1) * fDriftWindow;

		size_t max_idx = k0;
		float max_adc = adc[k0] * fLifetimeCorrFactors[k0 + tick0];
		for (size_t k = k0 + 1; k < k1; ++k)
		{
			float ak = adc[k] * fLifetimeCorrFactors[k + tick0];
			if (ak > max_adc) { max_adc = ak; max_idx = k; }
		}

		size_t n = 1;
		if (max_idx > 0) { max_adc += adc[max_idx - 1] * fLifetimeCorrFactors[max_idx - 1 + tick0]; n++; }
		if (max_idx + 1 < adc.size()) { max_adc += adc[max_idx + 1] * fLifetimeCorrFactors[max_idx + 1 + tick0]; n++; }

		dst[i] = scaleAdcSample(max_adc / n);
	}
}

void nnet::DataProviderAlg::downscaleMean(std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) const
{
	for (size_t i = 0; i < dst.size(); ++i)
	{
		size_t k0 = i * fDriftWindow;
		size_t k1 = (i + 1) * fDriftWindow;

		float sum_adc = 0;
		for (size_t k = k0; k < k1; ++k)
		{
			sum_adc += adc[k] * fLifetimeCorrFactors[k + tick0];
		}

		if (sum_adc != 0) { dst[i] = scaleAdcSample(sum_adc * fDriftWindowInv); }
		else { dst[i] = 0; }
	}
}

bool nnet::DataProviderAlg::setWireData(std::vector<float> const & adc, size_t wireIdx)
{
   	if (wireIdx >= fWireDriftData.size()) return false;
   	auto & wData = fWireDriftData[wireIdx];

    if (fDownscaleFullView)
    {
        if (adc.size() / fDriftWindow <= fNCachedDrifts) { return downscale(wData, adc); }
        else { return false; }
    }
    else
    {
        if (adc.size() <= fNCachedDrifts) // copy ADC's, no downsampling nor scaling
        {
            for (size_t i = 0; i < adc.size(); ++i) { wData[i] = adc[i]; }
        }
        else { return false; }
    }
    return true;
}
// ------------------------------------------------------

bool nnet::DataProviderAlg::setWireDriftData(const std::vector<recob::Wire> & wires,
	unsigned int view, unsigned int tpc, unsigned int cryo)
{
	fCryo = cryo; fTPC = tpc; fView = view;

	size_t nwires = fGeometry->Nwires(view, tpc, cryo);
	size_t ndrifts = fDetProp->NumberTimeSamples();

	resizeView(nwires, ndrifts);

    bool allWrong = true;
    for (auto const & wire : wires)
	{
		auto wireChannelNumber = wire.Channel();

		size_t w_idx = 0;
		for (auto const& id : fGeometry->ChannelToWire(wireChannelNumber))
		{
			if ((id.Cryostat == cryo) && (id.TPC == tpc) && (id.Plane == view))
			{
			    w_idx = id.Wire;

			    auto adc = wire.Signal();
			    if (adc.size() < ndrifts)
			    {
			    	mf::LogWarning("DataProviderAlg") << "Wire ADC vector size lower than NumberTimeSamples.";
			    	continue; // not critical, maybe other wires are OK, so continue
			    }

			    if (!setWireData(adc, w_idx))
			    {
			    	mf::LogWarning("DataProviderAlg") << "Wire data not set.";
			    	continue; // also not critical, try to set other wires
			    }

			    fWireChannels[w_idx] = wireChannelNumber;
			    allWrong = false;
			}
		}
	}
	if (allWrong)
	{
	    mf::LogError("DataProviderAlg") << "Wires data not set in the cryo:"
	        << cryo << " tpc:" << tpc << " plane:" << view << " (skip this plane)";
	    return false;
	}
	
    applyBlur();
    addWhiteNoise();
    addCoherentNoise();
	
	return true;
}
// ------------------------------------------------------

float nnet::DataProviderAlg::scaleAdcSample(float val) const
{
    if (val < -50.) val = -50.;
    if (val > 150.) val = 150.;

    if (fCalibrateAmpl) { val *= fAmplCalibConst[fView]; }

    return 0.1 * val;
}
// ------------------------------------------------------

void nnet::DataProviderAlg::applyBlur()
{
    if (fBlurKernel.size() < 2) return;

    size_t margin_left = (fBlurKernel.size()-1) >> 1, margin_right = fBlurKernel.size() - margin_left - 1;

    std::vector< std::vector<float> > src(fWireDriftData.size());
    for (size_t w = 0; w < fWireDriftData.size(); ++w) { src[w] = fWireDriftData[w]; }

    for (size_t w = margin_left; w < fWireDriftData.size() - margin_right; ++w)
    {
        for (size_t d = 0; d < fWireDriftData[w].size(); ++d)
        {
            float sum = 0;
            for (size_t i = 0; i < fBlurKernel.size(); ++i)
            {
                sum += fBlurKernel[i] * src[w + i - margin_left][d];
            }
            fWireDriftData[w][d] = sum;
        }
    }
}
// ------------------------------------------------------

void nnet::DataProviderAlg::addWhiteNoise()
{
    if (fNoiseSigma == 0) return;

    double effectiveSigma = scaleAdcSample(fNoiseSigma);
    if (fDownscaleFullView) effectiveSigma /= fDriftWindow;

    CLHEP::RandGauss gauss(fRndEngine);
    std::vector<double> noise(fNCachedDrifts);
    for (auto & wire : fWireDriftData)
    {
        gauss.fireArray(fNCachedDrifts, noise.data(), 0., effectiveSigma);
        for (size_t d = 0; d < wire.size(); ++d)
        {
            wire[d] += noise[d];
        }
    }
}
// ------------------------------------------------------

void nnet::DataProviderAlg::addCoherentNoise()
{
    if (fCoherentSigma == 0) return;

    double effectiveSigma = scaleAdcSample(fCoherentSigma);
    if (fDownscaleFullView) effectiveSigma /= fDriftWindow;

    CLHEP::RandGauss gauss(fRndEngine);
    std::vector<double> amps1(fWireDriftData.size());
    std::vector<double> amps2(1 + (fWireDriftData.size() / 32));
    gauss.fireArray(amps1.size(), amps1.data(), 1., 0.1); // 10% wire-wire ampl. variation
    gauss.fireArray(amps2.size(), amps2.data(), 1., 0.1); // 10% group-group ampl. variation

    double group_amp = 1.0;
    std::vector<double> noise(fNCachedDrifts);
    for (size_t w = 0; w < fWireDriftData.size(); ++w)
    {
        if ((w & 31) == 0)
        {
            group_amp = amps2[w >> 5]; // div by 32
            gauss.fireArray(fNCachedDrifts, noise.data(), 0., effectiveSigma);
        } // every 32 wires

        auto & wire = fWireDriftData[w];
        for (size_t d = 0; d < wire.size(); ++d)
        {
            wire[d] += group_amp * amps1[w] * noise[d];
        }
    }
}
// ------------------------------------------------------


// ------------------------------------------------------
// -------------------ModelInterface---------------------
// ------------------------------------------------------

std::string nnet::ModelInterface::findFile(const char* fileName) const
{
    std::string fname_out;
    cet::search_path sp("FW_SEARCH_PATH");
    if (!sp.find_file(fileName, fname_out))
    {
        struct stat buffer;
        if (stat(fileName, &buffer) == 0) { fname_out = fileName; }
        else
        {
            throw art::Exception(art::errors::NotFound)
                << "Could not find the model file " << fileName;
        }
    }
    return fname_out;
}

// ------------------------------------------------------
// -----------------MlpModelInterface--------------------
// ------------------------------------------------------

nnet::MlpModelInterface::MlpModelInterface(const char* xmlFileName) :
	m(nnet::ModelInterface::findFile(xmlFileName).c_str())
{
	mf::LogInfo("MlpModelInterface") << "MLP model loaded.";
}
// ------------------------------------------------------

bool nnet::MlpModelInterface::Run(std::vector< std::vector<float> > const & inp2d)
{
	auto input = nnet::PointIdAlg::flattenData2D(inp2d);
	if (input.size() == m.GetInputLength())
	{
		m.Run(input);
		return true;
	}

	mf::LogError("MlpModelInterface") << "Flattened patch size does not match MLP model.";
	return false;
}
// ------------------------------------------------------

std::vector<float> nnet::MlpModelInterface::GetAllOutputs(void) const
{
	std::vector<float> result(m.GetOutputLength(), 0);
	for (size_t o = 0; o < result.size(); ++o)
	{
		result[o] = m.GetOneOutput(o);
	}
	return result;
}
// ------------------------------------------------------

float nnet::MlpModelInterface::GetOneOutput(int neuronIndex) const
{
	if ((int)neuronIndex < m.GetOutputLength())
	{
		return m.GetOneOutput(neuronIndex);
	}

	mf::LogError("MlpModelInterface") << "Output index does not match MLP model.";
	return 0.;
}
// ------------------------------------------------------


// ------------------------------------------------------
// ----------------KerasModelInterface-------------------
// ------------------------------------------------------

nnet::KerasModelInterface::KerasModelInterface(const char* modelFileName) :
	m(nnet::ModelInterface::findFile(modelFileName).c_str())
{
	mf::LogInfo("KerasModelInterface") << "Keras model loaded.";
}
// ------------------------------------------------------

bool nnet::KerasModelInterface::Run(std::vector< std::vector<float> > const & inp2d)
{
	std::vector< std::vector< std::vector<float> > > inp3d;
	inp3d.push_back(inp2d); // lots of copy, should add 2D to keras...

	keras::DataChunk *sample = new keras::DataChunk2D();
	sample->set_data(inp3d); // and more copy...
	fOutput = m.compute_output(sample);
	delete sample;

	return true;
}
// ------------------------------------------------------

std::vector<float> nnet::KerasModelInterface::GetAllOutputs(void) const
{
	return fOutput;
}
// ------------------------------------------------------

float nnet::KerasModelInterface::GetOneOutput(int neuronIndex) const
{
	if (neuronIndex < (int)fOutput.size()) return fOutput[neuronIndex];

	mf::LogError("KerasModelInterface") << "Output index does not match Keras model.";
	return 0.;
}
// ------------------------------------------------------


// ------------------------------------------------------
// --------------------PointIdAlg------------------------
// ------------------------------------------------------

nnet::PointIdAlg::PointIdAlg(const Config& config) : nnet::DataProviderAlg(config),
	fNNet(0),
	fPatchSizeW(32), fPatchSizeD(32),
	fCurrentWireIdx(99999), fCurrentScaledDrift(99999)
{
	this->reconfigure(config); 
}
// ------------------------------------------------------

nnet::PointIdAlg::~PointIdAlg(void)
{
	deleteNNet();
}
// ------------------------------------------------------

void nnet::PointIdAlg::reconfigure(const Config& config)
{
	fNNetModelFilePath = config.NNetModelFile();

	fPatchSizeW = config.PatchSizeW();
	fPatchSizeD = config.PatchSizeD();

	deleteNNet();

	if ((fNNetModelFilePath.length() > 4) &&
	    (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 4, 4, ".xml") == 0))
	{
		fNNet = new nnet::MlpModelInterface(fNNetModelFilePath.c_str());
	}
	else if ((fNNetModelFilePath.length() > 5) &&
	    (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 5, 5, ".nnet") == 0))
	{
		fNNet = new nnet::KerasModelInterface(fNNetModelFilePath.c_str());
	}
	else
	{
		mf::LogError("PointIdAlg") << "Loading model from file failed.";
		return;
	}

    resizePatch();
}
// ------------------------------------------------------

void nnet::PointIdAlg::resizePatch(void)
{
	fWireDriftPatch.resize(fPatchSizeW);
	for (auto & r : fWireDriftPatch) r.resize(fPatchSizeD);
}
// ------------------------------------------------------

size_t nnet::PointIdAlg::NClasses(void) const
{
	if (fNNet) return fNNet->GetOutputLength();
	else return 0;
}
// ------------------------------------------------------

float nnet::PointIdAlg::predictIdValue(unsigned int wire, float drift, size_t outIdx) const
{
	float result = 0.;

	if (!bufferPatch(wire, drift))
	{
		mf::LogError("PointIdAlg") << "Patch buffering failed.";
		return result;
	}

	if (fNNet)
	{
		if (fNNet->Run(fWireDriftPatch))
		{
			result = fNNet->GetOneOutput(outIdx);
		}
		else mf::LogError("PointIdAlg") << "Problem with applying model to input.";
	}

	return result;
}
// ------------------------------------------------------

std::vector<float> nnet::PointIdAlg::predictIdVector(unsigned int wire, float drift) const
{
	std::vector<float> result(NClasses(), 0);
	if (result.empty()) return result;

	if (!bufferPatch(wire, drift))
	{
		mf::LogError("PointIdAlg") << "Patch buffering failed.";
		return result;
	}

	if (fNNet)
	{
		if (fNNet->Run(fWireDriftPatch))
		{
			for (size_t o = 0; o < result.size(); ++o)
			{
				result[o] = fNNet->GetOneOutput(o);
			}
		}
		else mf::LogError("PointIdAlg") << "Problem with applying model to input.";
	}

	return result;
}
// ------------------------------------------------------

// MUST give the same result as get_patch() in scripts/utils.py
bool nnet::PointIdAlg::patchFromDownsampledView(size_t wire, float drift) const
{
	size_t sd = (size_t)(drift / fDriftWindow);
	if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == sd))
		return true; // still within the current position

	fCurrentWireIdx = wire;
	fCurrentScaledDrift = sd;

	int halfSizeW = fPatchSizeW / 2;
	int halfSizeD = fPatchSizeD / 2;

	int w0 = fCurrentWireIdx - halfSizeW;
	int w1 = fCurrentWireIdx + halfSizeW;

	int d0 = fCurrentScaledDrift - halfSizeD;
	int d1 = fCurrentScaledDrift + halfSizeD;

    int wsize = fWireDriftData.size();
	for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch)
	{
		auto & dst = fWireDriftPatch[wpatch];
		if ((w >= 0) && (w < wsize))
		{
			auto & src = fWireDriftData[w];
			int dsize = src.size();
			for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch)
			{
				if ((d >= 0) && (d < dsize))
				{
					dst[dpatch] = src[d];
				}
				else
				{
					dst[dpatch] = 0;
				}
			}
		}
		else
		{
			std::fill(dst.begin(), dst.end(), 0);
		}
	}

	return true;
}

bool nnet::PointIdAlg::patchFromOriginalView(size_t wire, float drift) const
{
	fCurrentWireIdx = wire;
	fCurrentScaledDrift = drift;

    int dsize = fDriftWindow * fPatchSizeD;
	int halfSizeW = fPatchSizeW / 2;
	int halfSizeD = dsize / 2;

	int w0 = fCurrentWireIdx - halfSizeW;
	int w1 = fCurrentWireIdx + halfSizeW;

	int d0 = fCurrentScaledDrift - halfSizeD;
	int d1 = fCurrentScaledDrift + halfSizeD;

    std::vector<float> tmp(dsize);
    int wsize = fWireDriftData.size();
	for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch)
	{
		if ((w >= 0) && (w < wsize))
		{
			auto & src = fWireDriftData[w];
			int src_size = src.size();
			for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch)
			{
				if ((d >= 0) && (d < src_size))
				{
					tmp[dpatch] = src[d];
				}
				else
				{
					tmp[dpatch] = 0;
				}
			}
		}
		else
		{
			std::fill(tmp.begin(), tmp.end(), 0);
		}

        downscale(fWireDriftPatch[wpatch], tmp, d0);
	}

	return true;
}
// ------------------------------------------------------

std::vector<float> nnet::PointIdAlg::flattenData2D(std::vector< std::vector<float> > const & patch)
{
	std::vector<float> flat;
	if (patch.empty() || patch.front().empty())
	{
		mf::LogError("DataProviderAlg") << "Patch is empty.";
		return flat;
	}

	flat.resize(patch.size() * patch.front().size());

	for (size_t w = 0, i = 0; w < patch.size(); ++w)
	{
		auto const & wire = patch[w];
		for (size_t d = 0; d < wire.size(); ++d, ++i)
		{
			flat[i] = wire[d];
		}
	}

	return flat;
}
// ------------------------------------------------------

bool nnet::PointIdAlg::isInsideFiducialRegion(unsigned int wire, float drift) const
{
	size_t marginW = fPatchSizeW / 8; // fPatchSizeX/2 will make patch always completely filled
	size_t marginD = fPatchSizeD / 8;

	size_t scaledDrift = (size_t)(drift / fDriftWindow);
	if ((wire >= marginW) && (wire < fNWires - marginW) &&
	    (scaledDrift >= marginD) && (scaledDrift < fNScaledDrifts - marginD)) return true;
	else return false;
}
// ------------------------------------------------------

// ------------------------------------------------------
// ------------------TrainingDataAlg---------------------
// ------------------------------------------------------

nnet::TrainingDataAlg::TrainingDataAlg(const Config& config) : nnet::DataProviderAlg(config)
{
	this->reconfigure(config); 
}
// ------------------------------------------------------

nnet::TrainingDataAlg::~TrainingDataAlg(void)
{
}
// ------------------------------------------------------

void nnet::TrainingDataAlg::reconfigure(const Config& config)
{
	fWireProducerLabel = config.WireLabel();
	fHitProducerLabel = config.HitLabel();
	fTrackModuleLabel = config.TrackLabel();
	fSimulationProducerLabel = config.SimulationLabel();
	fSaveVtxFlags = config.SaveVtxFlags();

    fAdcDelay = config.AdcDelayTicks();

	for(int x = 0; x < 100; x++) {
	  events_per_bin.push_back(0);
	}


}
// ------------------------------------------------------

void nnet::TrainingDataAlg::resizeView(size_t wires, size_t drifts)
{
	nnet::DataProviderAlg::resizeView(wires, drifts);

	fWireDriftEdep.resize(wires);
	for (auto & w : fWireDriftEdep)
	{
		w.resize(fNCachedDrifts);
		std::fill(w.begin(), w.end(), 0.0F);
	}

	fWireDriftPdg.resize(wires);
	for (auto & w : fWireDriftPdg)
	{
		w.resize(fNCachedDrifts);
		std::fill(w.begin(), w.end(), 0);
	}
}
// ------------------------------------------------------

bool nnet::TrainingDataAlg::setWireEdepsAndLabels(
	std::vector<float> const & edeps, std::vector<int> const & pdgs, size_t wireIdx)
{
	if ((wireIdx >= fWireDriftEdep.size()) || (edeps.size() != pdgs.size())) { return false; }

	size_t dstep = 1;
	if (fDownscaleFullView) { dstep = fDriftWindow; }

	if (edeps.size() / dstep > fNCachedDrifts) { return false; }

	auto & wEdep = fWireDriftEdep[wireIdx];
	auto & wPdg = fWireDriftPdg[wireIdx];

	for (size_t i = 0; i < fNCachedDrifts; ++i)
	{
		size_t i0 = i * dstep;
		size_t i1 = (i + 1) * dstep;

		int best_pdg = pdgs[i0] & nnet::TrainingDataAlg::kPdgMask;
		int vtx_flags = pdgs[i0] & nnet::TrainingDataAlg::kVtxMask;
		int type_flags = pdgs[i0] & nnet::TrainingDataAlg::kTypeMask;
		float max_edep = edeps[i0];
		for (size_t k = i0 + 1; k < i1; ++k)
		{
			float ek = edeps[k];
			if (ek > max_edep)
			{
				max_edep = ek;
				best_pdg = pdgs[k] & nnet::TrainingDataAlg::kPdgMask; // remember best matching pdg
			}
			type_flags |= pdgs[k] & nnet::TrainingDataAlg::kTypeMask; // accumulate track type flags
			vtx_flags |= pdgs[k] & nnet::TrainingDataAlg::kVtxMask;   // accumulate all vtx flags
		}

		wEdep[i] = max_edep;

        best_pdg |= type_flags;
		if (fSaveVtxFlags) best_pdg |= vtx_flags;
		wPdg[i] = best_pdg;
	}

	return true;
}
// ------------------------------------------------------

nnet::TrainingDataAlg::WireDrift nnet::TrainingDataAlg::getProjection(const TLorentzVector& tvec, unsigned int view) const
{
	auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	nnet::TrainingDataAlg::WireDrift wd;
	wd.Wire = 0; wd.Drift = 0; wd.TPC = -1; wd.Cryo = -1;

    try
    {
    	double vtx[3] = {tvec.X(), tvec.Y(), tvec.Z()};
	    if (fGeometry->FindTPCAtPosition(vtx).isValid)
	    {
	    	geo::TPCID tpcid = fGeometry->FindTPCAtPosition(vtx);
	    	unsigned int tpc = tpcid.TPC, cryo = tpcid.Cryostat;

	    	// correct for the time offset
	    	float dx = tvec.T() * 1.e-3 * detprop->DriftVelocity();
	    	int driftDir = fGeometry->TPC(tpcid).DetectDriftDirection();
			if (driftDir == 1) { dx *= -1; }
			else if (driftDir != -1)
			{
			    throw cet::exception("nnet::TrainingDataAlg") << "drift direction is not X." << std::endl;
			}
			vtx[0] = tvec.X() + dx;
				
	    	wd.Wire = fGeometry->NearestWire(vtx, view, tpc, cryo);
	    	wd.Drift = fDetProp->ConvertXToTicks(vtx[0], view, tpc, cryo);
	    	wd.TPC = tpc; wd.Cryo = cryo;
	    }
	}
	catch (const geo::InvalidWireIDError & e)
	{
	    mf::LogWarning("TrainingDataAlg") << "Vertex projection out of wire planes, just skipping this vertex.";
	}
	catch (...)
	{
	    mf::LogWarning("TrainingDataAlg") << "Vertex projection out of wire planes, skip MC vertex.";
	}
	return wd;
}
// ------------------------------------------------------

bool nnet::TrainingDataAlg::isElectronEnd(const simb::MCParticle & particle,
    const std::unordered_map< int, const simb::MCParticle* > & particleMap) const
{
    const float minElectronLength2 = 2.5*2.5;
    const float maxDeltaLength2 = 0.15*0.15;

    int pdg = abs(particle.PdgCode());
	if (pdg != 11) return false; // should be applied only to electrons

	size_t nSec = particle.NumberDaughters();
	for (size_t d = 0; d < nSec; ++d)
	{
		auto d_search = particleMap.find(particle.Daughter(d));
		if (d_search != particleMap.end())
		{
			auto const & daughter = *((*d_search).second);
			int d_pdg = abs(daughter.PdgCode());
			if (d_pdg != 22) { return false; } // not the end of the shower
		}
	}

    float trkLength2 = 0;
	auto const * p = &particle;
	bool branching = false;
	while (!branching)
	{
        trkLength2 += particleRange2(*p);
        auto m_search = particleMap.find(p->Mother());
		if (m_search != particleMap.end())
        {
			p = (*m_search).second;
			int m_pdg = abs(p->PdgCode());
			if (m_pdg == 11)
			{
			    nSec = p->NumberDaughters();
			    size_t ne = 0;
			    for (size_t d = 0; d < nSec; ++d)
			    {
			        auto d_search = particleMap.find(p->Daughter(d));
			        if (d_search != particleMap.end())
			        {
			            auto const & daughter = *((*d_search).second);
			            int d_pdg = abs(daughter.PdgCode());
			            if (d_pdg == 11)
			            {
                        	if (particleRange2(daughter) > maxDeltaLength2) { ne++; }
			            }
			        }
			    }
			    if (ne > 1) { branching = true; }
			}
			else break;
        }
        else break;
    }

    return (trkLength2 > minElectronLength2);
}

bool nnet::TrainingDataAlg::isMuonDecaying(const simb::MCParticle & particle,
    const std::unordered_map< int, const simb::MCParticle* > & particleMap) const
{
    bool hasElectron = false, hasNuMu = false, hasNuE = false;

    int pdg = abs(particle.PdgCode());
	if ((pdg == 13) && (particle.EndProcess() == "FastScintillation")) // potential muon decay at rest
	{
		unsigned int nSec = particle.NumberDaughters();
		for (size_t d = 0; d < nSec; ++d)
		{
			auto d_search = particleMap.find(particle.Daughter(d));
			if (d_search != particleMap.end())
			{
				auto const & daughter = *((*d_search).second);
				int d_pdg = abs(daughter.PdgCode());
				if (d_pdg == 11) hasElectron = true;
				else if (d_pdg == 14) hasNuMu = true;
				else if (d_pdg == 12) hasNuE = true;
			}
		}
	}

	return (hasElectron && hasNuMu && hasNuE);
}

void nnet::TrainingDataAlg::collectVtxFlags(
	std::unordered_map< size_t, std::unordered_map< int, int > > & wireToDriftToVtxFlags,
	const std::unordered_map< int, const simb::MCParticle* > & particleMap,
	unsigned int view) const
{
	std::cout << "collectVtxFlags" << std::endl;
	for (auto const & p : particleMap)
	{
		auto const & particle = *p.second;

		double ekStart = 1000. * (particle.E() - particle.Mass());
		double ekEnd = 1000. * (particle.EndE() - particle.Mass());

		int pdg = abs(particle.PdgCode());
		int flagsStart = nnet::TrainingDataAlg::kNone;
		int flagsEnd = nnet::TrainingDataAlg::kNone;
	
		switch (pdg)
		{
			case 22:   // gamma
				if ((particle.EndProcess() == "conv") &&
				    (ekStart > 40.0)) // conversion, gamma > 40MeV
				{
					//std::cout << "---> gamma conversion at " << ekStart << std::endl;
					flagsEnd = nnet::TrainingDataAlg::kConv;
				}
				break;

			case 11:   // e+/-
			    if (isElectronEnd(particle, particleMap))
			    {
			        flagsEnd = nnet::TrainingDataAlg::kElectronEnd;
			    }
				break;


			case 13:   // mu+/-
			    if (isMuonDecaying(particle, particleMap))
			    {
			        //std::cout << "---> mu decay to electron" << std::endl;
			        flagsEnd = nnet::TrainingDataAlg::kDecay;
			    }
				break;

			case 111:  // pi0
				//std::cout << "---> pi0" << std::endl;
				flagsStart = nnet::TrainingDataAlg::kPi0;
				break;

			case 321:  // K+/-
			case 211:  // pi+/-
			case 2212: // proton
				if (ekStart > 50.0)
				{
					if (particle.Mother() != 0)
					{
						auto search = particleMap.find(particle.Mother());
						if (search != particleMap.end())
						{
							auto const & mother = *((*search).second);
							int m_pdg = abs(mother.PdgCode());
							unsigned int nSec = mother.NumberDaughters();
							unsigned int nVisible = 0;
							if (nSec > 1)
							{
								for (size_t d = 0; d < nSec; ++d)
								{
									auto d_search = particleMap.find(mother.Daughter(d));
									if (d_search != particleMap.end())
									{
										auto const & daughter = *((*d_search).second);
										int d_pdg = abs(daughter.PdgCode());
										if (((d_pdg == 2212) || (d_pdg == 211) || (d_pdg == 321)) &&
										    (1000. * (daughter.E() - daughter.Mass()) > 50.0))
										{
											++nVisible;
										}
									}
								}
							}
							// hadron with Ek > 50MeV (so well visible) and
							// produced by another hadron (but not neutron, so not single track from nothing) or
							// at least secondary hadrons with Ek > 50MeV (so this is a good kink or V-like)
							if (((m_pdg != pdg) && (m_pdg != 2112)) || ((m_pdg != 2112) && (nVisible > 0)) || ((m_pdg == 2112) && (nVisible > 1)))
							{
								// std::cout << "---> hadron at " << ekStart
								//	<< ", pdg: " << pdg << ", mother pdg: " << m_pdg
								//	<< ", vis.daughters: " << nVisible << std::endl;
								flagsStart = nnet::TrainingDataAlg::kHadr;
							} 
						}
						// else std::cout << "---> mother not found for tid: " << particle.Mother() << std::endl;
					}

					if (particle.EndProcess() == "FastScintillation") // potential decay at rest
					{
						unsigned int nSec = particle.NumberDaughters();
						for (size_t d = 0; d < nSec; ++d)
						{
							auto d_search = particleMap.find(particle.Daughter(d));
							if (d_search != particleMap.end())
							{
								auto const & daughter = *((*d_search).second);
								int d_pdg = abs(daughter.PdgCode());
								if ((pdg == 321) && (d_pdg == 13))
								{
									//std::cout << "---> K decay to mu" << std::endl;
									flagsEnd = nnet::TrainingDataAlg::kDecay;
									break;
								}
								if ((pdg == 211) && (d_pdg == 13))
								{
									//std::cout << "---> pi decay to mu" << std::endl;
									flagsEnd = nnet::TrainingDataAlg::kDecay;
									break;
								}
							}
						}
					}

					if ((particle.EndProcess() == "Decay") && (ekEnd > 200.0)) // decay in flight
					{
						unsigned int nSec = particle.NumberDaughters();
						for (size_t d = 0; d < nSec; ++d)
						{
							auto d_search = particleMap.find(particle.Daughter(d));
							if (d_search != particleMap.end())
							{
								auto const & daughter = *((*d_search).second);
								int d_pdg = abs(daughter.PdgCode());
								if ((pdg == 321) && (d_pdg == 13))
								{
									//std::cout << "---> in-flight K decay to mu" << std::endl;
									flagsEnd = nnet::TrainingDataAlg::kHadr;
									break;
								}
								if ((pdg == 211) && (d_pdg == 13))
								{
									//std::cout << "---> in-flight pi decay to mu" << std::endl;
									flagsEnd = nnet::TrainingDataAlg::kHadr;
									break;
								}
							}
						}
					}
				}
				break;

			default: continue;
		}
		
		if (particle.Process() == "primary")
		{
			flagsStart |= nnet::TrainingDataAlg::kNuPri;
		}
		
		
		if (flagsStart != nnet::TrainingDataAlg::kNone)
		{
			auto wd = getProjection(particle.Position(), view);
			
			if ((wd.TPC == (int)fTPC) && (wd.Cryo == (int)fCryo))
			{
				wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= flagsStart;
				// std::cout << "---> flagsStart:" << flagsStart << " view:" << view << " wire:" << wd.Wire << " drift:" << wd.Drift << std::endl;
			}
			// else std::cout << "---> not in current TPC" << std::endl;
		}
		if (flagsEnd != nnet::TrainingDataAlg::kNone)
		{
			auto wd = getProjection(particle.EndPosition(), view);
			if ((wd.TPC == (int)fTPC) && (wd.Cryo == (int)fCryo))
			{
			    //if (flagsEnd == nnet::TrainingDataAlg::kElectronEnd) { std::cout << "---> clear electron endpoint" << std::endl; }
				wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= flagsEnd;
				//if (flagsEnd == nnet::TrainingDataAlg::kElectronEnd)
				//    std::cout << "---> flagsEnd:" << flagsEnd << " view:" << view << " wire:" << wd.Wire << " drift:" << wd.Drift << std::endl;
			}
			// else std::cout << "---> not in current TPC" << std::endl;
		}

		//if (ekStart > 30.0)
		//{
		//	std::cout << particle.PdgCode() << ", " << ekStart << ": "
		//		<< particle.Process() << " --> " << particle.EndProcess()
		//		<< " " << ekEnd	<< std::endl;
		//}
	}
}
// ------------------------------------------------------

bool nnet::TrainingDataAlg::setDataEventData(const art::Event& event,
	unsigned int view, unsigned int tpc, unsigned int cryo)
{

  art::Handle< std::vector<recob::Wire> > wireHandle;
  std::vector< art::Ptr<recob::Wire> > Wirelist;

  if(event.getByLabel(fWireProducerLabel, wireHandle))
    art::fill_ptr_vector(Wirelist, wireHandle);

  if(!setWireDriftData(*wireHandle, view, tpc, cryo)) {
    mf::LogError("TrainingDataAlg") << "Wire data not set.";
    return false;
  }

  // Hit info
  art::Handle< std::vector<recob::Hit> > HitHandle;
  std::vector< art::Ptr<recob::Hit> > Hitlist;

  if(event.getByLabel(fHitProducerLabel, HitHandle))
    art::fill_ptr_vector(Hitlist, HitHandle);

  // Track info
  art::Handle< std::vector<recob::Track> > TrackHandle;
  std::vector< art::Ptr<recob::Track> > Tracklist;

  if(event.getByLabel(fTrackModuleLabel, TrackHandle))
    art::fill_ptr_vector(Tracklist, TrackHandle);

  art::FindManyP<recob::Track> ass_trk_hits(HitHandle,   event, fTrackModuleLabel);

  // Loop over wires (sorry about hard coded value) to fill in 1) pdg and 2) charge depo
  for (size_t widx = 0; widx < 240; ++widx) {
    
    std::vector< float > labels_deposit(fNDrifts, 0);  // full-drift-length buffers
    std::vector< int > labels_pdg(fNDrifts, 0);

    // First, the charge depo
    for(size_t subwidx = 0; subwidx < Wirelist.size(); ++subwidx) {
      if(widx+240 == Wirelist[subwidx]->Channel()) {	
	labels_deposit = Wirelist[subwidx]->Signal();
	break;
      }
    }
   
    // Second, the pdg code
    // This code finds the angle of the track and records
    //  events based on its angle to try to get an isometric sample
    //  instead of just a bunch of straight tracks

    // Meta code:
    // For each hit:
    //  find farthest hit from point
    //  then find farthest hit from THAT one
    //  should be start and end of track, then just use trig

    for(size_t iHit = 0; iHit < Hitlist.size(); ++iHit) {

      if(Hitlist[iHit]->Channel() != widx+240) { continue; }
      if(Hitlist[iHit]->View() != 1) { continue; }

      // Make sure there is a track association
      if(ass_trk_hits.at(iHit).size() == 0) { continue; }
      
      // Not sure about this
      // Cutting on length to not just get a bunch of shower stubs
      // Might add a lot of bias though
      if(ass_trk_hits.at(iHit)[0]->Length() < 5) { continue; }

      // Search for farest hit from this one
      int far_index = 0;
      double far_dist = 0;

      for(size_t jHit = 0; jHit < Hitlist.size(); ++jHit) {
	if(jHit == iHit) { continue; }
	if(Hitlist[jHit]->View() != 1) { continue; }

	if(ass_trk_hits.at(jHit).size() == 0) { continue; }
	if(ass_trk_hits.at(jHit)[0]->ID() != 
	   ass_trk_hits.at(iHit)[0]->ID()) { continue; }

	double dist = sqrt((Hitlist[iHit]->Channel()-Hitlist[jHit]->Channel()) *
			   (Hitlist[iHit]->Channel()-Hitlist[jHit]->Channel()) +
			   (Hitlist[iHit]->PeakTime()-Hitlist[jHit]->PeakTime()) * 
			   (Hitlist[iHit]->PeakTime()-Hitlist[jHit]->PeakTime()));

	if(far_dist < dist){
	  far_dist = dist;
	  far_index = jHit;
	}
      }

      // Search for the other end of the track
      int other_end = 0;
      int other_dist = 0;

      for(size_t jHit = 0; jHit < Hitlist.size(); ++jHit) {
	if(jHit == iHit or int(jHit) == far_index) { continue; }
	if(Hitlist[jHit]->View() != 1) { continue; }

	if(ass_trk_hits.at(jHit).size() == 0) { continue; }
	if(ass_trk_hits.at(jHit)[0]->ID() != 
	   ass_trk_hits.at(iHit)[0]->ID()) { continue; }

	double dist = sqrt((Hitlist[far_index]->Channel()-Hitlist[jHit]->Channel()) *
			   (Hitlist[far_index]->Channel()-Hitlist[jHit]->Channel()) +
			   (Hitlist[far_index]->PeakTime()-Hitlist[jHit]->PeakTime()) * 
			   (Hitlist[far_index]->PeakTime()-Hitlist[jHit]->PeakTime()));

	if(other_dist < dist){
	  other_dist = dist;
	  other_end = jHit;
	}
      }

      // We have the end points now
      double del_wire = double(Hitlist[other_end]->Channel() - Hitlist[far_index]->Channel());
      double del_time = double(Hitlist[other_end]->PeakTime() - Hitlist[far_index]->PeakTime());
      double hypo = sqrt(del_wire * del_wire + del_time * del_time);      

      if(hypo == 0) { continue; } // Should never happen, but doing it anyway

      double cosser = TMath::Abs(del_wire / hypo);
      double norm_ang = TMath::ACos(cosser) * 2 / TMath::Pi();

      // Using events_per_bin to keep track of number of hits per angle (normalized to 0 to 1)

      int binner = int(norm_ang * 100);      
      if(binner == 100) { binner = 99; } // Dealing with rounding errors

      // So we should get a total of 5000 * 100 = 50,000 if we use the whole set
      if(events_per_bin.at(binner) > 5000) { continue; }
      events_per_bin.at(binner) += 1;
      
      // If survives everything, saves the pdg
      labels_pdg[Hitlist[iHit]->PeakTime()] = 211; // Same as pion for now

    }

    setWireEdepsAndLabels(labels_deposit, labels_pdg, widx);

  } // for each Wire

  /*
  for(size_t i = 0; i < events_per_bin.size(); i ++) {
    std::cout << i << ") " << events_per_bin[i] << " - ";
  }
  */

  return true;

}

bool nnet::TrainingDataAlg::setEventData(const art::Event& event,
	unsigned int view, unsigned int tpc, unsigned int cryo)
{
	art::ValidHandle< std::vector<recob::Wire> > wireHandle
		= event.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);

	if (!setWireDriftData(*wireHandle, view, tpc, cryo))
	{
		mf::LogError("TrainingDataAlg") << "Wire data not set.";
		return false;
	}

	art::ServiceHandle<sim::LArG4Parameters> larParameters;
	double electronsToGeV = 1. / larParameters->GeVToElectrons();

	auto particleHandle = event.getValidHandle< std::vector<simb::MCParticle> >(fSimulationProducerLabel);

	auto simChannelHandle = event.getValidHandle< std::vector<sim::SimChannel> >(fSimulationProducerLabel);

    std::unordered_map< int, const simb::MCParticle* > particleMap;
	for (auto const & particle : *particleHandle)
    {
		particleMap[particle.TrackId()] = &particle;
	}

	std::unordered_map< size_t, std::unordered_map< int, int > > wireToDriftToVtxFlags;
	if (fSaveVtxFlags) collectVtxFlags(wireToDriftToVtxFlags, particleMap, view);

	std::map< int, int > trackToPDG;
    for (size_t widx = 0; widx < fNWires; ++widx)
	{
		auto wireChannelNumber = fWireChannels[widx];
		if (wireChannelNumber == raw::InvalidChannelID) continue;

		std::vector< float > labels_deposit(fNDrifts, 0);         // full-drift-length buffers,
		std::vector< int > labels_pdg(labels_deposit.size(), 0);  // both of the same size,
		int labels_size = labels_deposit.size();                  // cached as int for comparisons below

		std::map< int, std::map< int, double > > timeToTrackToCharge;
		for (auto const & channel : *simChannelHandle)
		{
			if (channel.Channel() != wireChannelNumber) continue;		

			auto const & timeSlices = channel.TDCIDEMap();
			for (auto const & timeSlice : timeSlices)
			{
				int time = timeSlice.first;

				auto const & energyDeposits = timeSlice.second;
				for (auto const & energyDeposit : energyDeposits)
				{
					int pdg = 0;
					int tid = energyDeposit.trackID;
					if (tid < 0) // negative tid means it is EM activity, and -tid is the mother
					{
					    pdg = 11; tid = -tid;
					    
						auto search = particleMap.find(tid);
					    if (search == particleMap.end())
					    {
						    mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND";
						    continue;
					    }
					    auto const & mother = *((*search).second); // mother particle of this EM
    					int mPdg = abs(mother.PdgCode());
                        if ((mPdg == 13) || (mPdg == 211) || (mPdg == 2212))
                        {
                            if (energyDeposit.numElectrons > 10) pdg |= nnet::TrainingDataAlg::kDelta; // tag delta ray
                        }
					}
					else
					{
						auto search = particleMap.find(tid);
					    if (search == particleMap.end())
					    {
						    mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND";
						    continue;
					    }
					    auto const & particle = *((*search).second);
					    pdg = abs(particle.PdgCode());

                        auto msearch = particleMap.find(particle.Mother());
	    				if (msearch != particleMap.end())
	    				{
	    				    auto const & mother = *((*msearch).second);
                            if (pdg == 11) // electron, check if it is Michel or primary electron
                            {
	    		                if (nnet::TrainingDataAlg::isMuonDecaying(mother, particleMap))
	    		                {
                			        pdg |= nnet::TrainingDataAlg::kMichel; // tag Michel
	    		                }
	    		                else if (mother.Mother() < 0)
	    		                {
	    		                    pdg |= nnet::TrainingDataAlg::kPriEl; // tag primary
	    		                }
                            }
                            else if (pdg == 13) // muon, check if primary
                            {
                                if (mother.Mother() < 0)
                                {
                                    pdg |= nnet::TrainingDataAlg::kPriMu; // tag primary
                                }
                            }
                        }
					}

					trackToPDG[energyDeposit.trackID] = pdg;

					double energy = energyDeposit.numElectrons * electronsToGeV;
					timeToTrackToCharge[time][energyDeposit.trackID] += energy;

	      		} // loop over energy deposits
      		} // loop over time slices
      	} // for each SimChannel

        int type_pdg_mask = nnet::TrainingDataAlg::kTypeMask | nnet::TrainingDataAlg::kPdgMask;
		for (auto const & ttc : timeToTrackToCharge)
		{
			float max_deposit = 0.0;
			int max_pdg = 0;
			for (auto const & tc : ttc.second) {

				if( tc.second > max_deposit ) 
				{
					max_deposit = tc.second;
					max_pdg = trackToPDG[tc.first];
				}			
			}

			if (ttc.first < labels_size)
			{
			    int tick_idx = ttc.first + fAdcDelay;
				if (tick_idx < labels_size)
				{
				    labels_deposit[tick_idx] = max_deposit;
				    labels_pdg[tick_idx] = max_pdg & type_pdg_mask;
				}
			}
		}

		for (auto const & drift_flags : wireToDriftToVtxFlags[widx])
		{
			int drift = drift_flags.first, flags = drift_flags.second;
			if ((drift >= 0) && (drift < labels_size))
			{
				labels_pdg[drift] |= flags;
			}
		}

		setWireEdepsAndLabels(labels_deposit, labels_pdg, widx);

	} // for each Wire

	return true;
}
// ------------------------------------------------------

bool nnet::TrainingDataAlg::findCrop(float max_e_cut, unsigned int & w0, unsigned int & w1, unsigned int & d0, unsigned int & d1) const
{
    if (fWireDriftEdep.empty() || fWireDriftEdep.front().empty()) return false;

    float max_cut = 0.25 * max_e_cut;

    w0 = 0;
    float cut = 0;
    while (w0 < fWireDriftEdep.size())
    {
        for (auto const d : fWireDriftEdep[w0]) cut += d;
        if (cut < max_cut) w0++;
        else break;
    }
    w1 = fWireDriftEdep.size() - 1;
    cut = 0;
    while (w1 > w0)
    {
        for (auto const d : fWireDriftEdep[w1]) cut += d;
        if (cut < max_cut) w1--;
        else break;
    }
    w1++;

    d0 = 0;
    cut = 0;
    while (d0 < fWireDriftEdep.front().size())
    {
        for (size_t i = w0; i < w1; ++i) cut += fWireDriftEdep[i][d0];
        if (cut < max_cut) d0++;
        else break;
    }
    d1 = fWireDriftEdep.front().size() - 1;
    cut = 0;
    while (d1 > d0)
    {
        for (size_t i = w0; i < w1; ++i) cut += fWireDriftEdep[i][d1];
        if (cut < max_cut) d1--;
        else break;
    }
    d1++;

    unsigned int margin = 20;
    if ((w1 - w0 > 8) && (d1 - d0 > 8))
    {
        if (w0 < margin) w0 = 0;
        else w0 -= margin;

        if (w1 > fWireDriftEdep.size() - margin) w1 = fWireDriftEdep.size();
        else w1 += margin;
        
        if (d0 < margin) d0 = 0;
        else d0 -= margin;
        
        if (d1 > fWireDriftEdep.front().size() - margin) d1 = fWireDriftEdep.front().size();
        else d1 += margin;
        
        return true;
    }
    else return false;
}
// ------------------------------------------------------


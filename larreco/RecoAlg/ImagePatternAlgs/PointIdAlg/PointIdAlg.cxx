////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan, May 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

nnet::DataProviderAlg::DataProviderAlg(const fhicl::ParameterSet& pset) :
	fCryo(9999), fTPC(9999), fView(9999),
	fNWires(0), fNDrifts(0), fNScaledDrifts(0),
	fDriftWindow(10), fPatchSizeW(32), fPatchSizeD(32),
	fDownscaleMode(nnet::DataProviderAlg::kMax),
	fCurrentWireIdx(99999), fCurrentScaledDrift(99999),
	fCalorimetryAlg(pset.get< fhicl::ParameterSet >("CalorimetryAlg")),
	fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
	fGeometry = &*(art::ServiceHandle<geo::Geometry>());

	this->reconfigure(pset); 
}
// ------------------------------------------------------

nnet::DataProviderAlg::~DataProviderAlg(void)
{
}
// ------------------------------------------------------

void nnet::DataProviderAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fCalorimetryAlg.reconfigure(p.get< fhicl::ParameterSet >("CalorimetryAlg"));

	fDriftWindow = p.get< unsigned int >("DriftWindow");
	fPatchSizeW = p.get< unsigned int >("PatchSizeW");
	fPatchSizeD = p.get< unsigned int >("PatchSizeD");

	std::string mode_str = p.get< std::string >("DownscaleFn");
	if (mode_str == "maxpool")      fDownscaleMode = nnet::DataProviderAlg::kMax;
	else if (mode_str == "maxmean") fDownscaleMode = nnet::DataProviderAlg::kMaxMean;
	else if (mode_str == "mean")    fDownscaleMode = nnet::DataProviderAlg::kMean;
	else
	{
		mf::LogError("DataProviderAlg") << "Downscale mode string not recognized, set to max pooling.";
		fDownscaleMode = nnet::DataProviderAlg::kMax;
	}

	resizePatch();
}
// ------------------------------------------------------

void nnet::DataProviderAlg::resizePatch(void)
{
	fWireDriftPatch.resize(fPatchSizeW);
	for (auto & r : fWireDriftPatch) r.resize(fPatchSizeD);
}
// ------------------------------------------------------

void nnet::DataProviderAlg::resizeView(size_t wires, size_t drifts)
{
	fNWires = wires; fNDrifts = drifts;
	fNScaledDrifts = drifts / fDriftWindow;

	fWireChannels.resize(wires);
	std::fill(fWireChannels.begin(), fWireChannels.end(), raw::InvalidChannelID);

	fWireDriftData.resize(wires);
	for (auto & w : fWireDriftData)
	{
		w.resize(fNScaledDrifts);
		std::fill(w.begin(), w.end(), 0.0F);
	}
}
// ------------------------------------------------------

float nnet::DataProviderAlg::scaleAdcSample(float val) const
{
	if (val < -50.) val = -50.;
	if (val > 150.) val = 150.;
	return 0.1 * val;
}

void nnet::DataProviderAlg::downscaleMax(std::vector<float> & dst, std::vector<float> const & adc) const
{
	for (size_t i = 0; i < dst.size(); ++i)
	{
		size_t i0 = i * fDriftWindow;
		size_t i1 = (i + 1) * fDriftWindow;

		float max_adc = scaleAdcSample(adc[i0] * fCalorimetryAlg.LifetimeCorrection(i0));
		for (size_t k = i0 + 1; k < i1; ++k)
		{
			float ak = scaleAdcSample(adc[k] * fCalorimetryAlg.LifetimeCorrection(k));
			if (ak > max_adc) max_adc = ak;
		}

		dst[i] = max_adc;
	}
}

void nnet::DataProviderAlg::downscaleMaxMean(std::vector<float> & dst, std::vector<float> const & adc) const
{
	for (size_t i = 0; i < dst.size(); ++i)
	{
		size_t i0 = i * fDriftWindow;
		size_t i1 = (i + 1) * fDriftWindow;

		size_t max_idx = i0;
		float max_adc = scaleAdcSample(adc[i0] * fCalorimetryAlg.LifetimeCorrection(i0));
		for (size_t k = i0 + 1; k < i1; ++k)
		{
			float ak = scaleAdcSample(adc[k] * fCalorimetryAlg.LifetimeCorrection(k));
			if (ak > max_adc) { max_adc = ak; max_idx = k; }
		}

		size_t n = 1;
		if (max_idx > 0) { max_adc += scaleAdcSample(adc[max_idx - 1] * fCalorimetryAlg.LifetimeCorrection(max_idx - 1)); n++; }
		if (max_idx + 1 < adc.size()) { max_adc += scaleAdcSample(adc[max_idx + 1] * fCalorimetryAlg.LifetimeCorrection(max_idx + 1)); n++; }

		dst[i] = max_adc / n;
	}
}

void nnet::DataProviderAlg::downscaleMean(std::vector<float> & dst, std::vector<float> const & adc) const
{
	for (size_t i = 0; i < dst.size(); ++i)
	{
		size_t i0 = i * fDriftWindow;
		size_t i1 = (i + 1) * fDriftWindow;

		float sum_adc = 0;
		for (size_t k = i0; k < i1; ++k)
		{
			sum_adc += scaleAdcSample(adc[k] * fCalorimetryAlg.LifetimeCorrection(k));
		}

		dst[i] = sum_adc / fDriftWindow;
	}
}

bool nnet::DataProviderAlg::setWireData(std::vector<float> const & adc, size_t wireIdx)
{
	if ((wireIdx >= fWireDriftData.size()) ||
	    (adc.size() / fDriftWindow > fNScaledDrifts)) return false;

	auto & wData = fWireDriftData[wireIdx];

	switch (fDownscaleMode)
	{
		case nnet::DataProviderAlg::kMax:     downscaleMax(wData, adc);     break;
		case nnet::DataProviderAlg::kMaxMean: downscaleMaxMean(wData, adc); break;
		case nnet::DataProviderAlg::kMean:    downscaleMean(wData, adc);    break;
		default: return false;
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

    for (auto const & wire : wires)
	{
		auto wireChannelNumber = wire.Channel();

		size_t w_idx = 0;
		bool right_plane = false;
		for (auto const& id : fGeometry->ChannelToWire(wireChannelNumber))
		{
			w_idx = id.Wire;

			if ((id.Cryostat == cryo) && (id.TPC == tpc) && (id.Plane == view))
			{
				right_plane = true; break;
			}
		}
		if (right_plane)
		{
			auto adc = wire.Signal();
			if (adc.size() != ndrifts)
			{
				mf::LogError("DataProviderAlg") << "ADC vector sizes not match.";
				return false;
			}

			if (!setWireData(adc, w_idx))
			{
				mf::LogError("DataProviderAlg") << "Wire data not set.";
				return false;
			}

			fWireChannels[w_idx] = wireChannelNumber;
		}
	}
	return true;
}
// ------------------------------------------------------

bool nnet::DataProviderAlg::bufferPatch(size_t wire, float drift) const
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

	for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch)
	{
		auto & dst = fWireDriftPatch[wpatch];
		if ((w >= 0) && (w < (int)fWireDriftData.size()))
		{
			auto & src = fWireDriftData[w];
			for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch)
			{
				if ((d >= 0) && (d < (int)src.size()))
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
// ------------------------------------------------------

std::vector<float> nnet::DataProviderAlg::flattenData2D(std::vector< std::vector<float> > const & patch)
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

bool nnet::DataProviderAlg::isInsideFiducialRegion(unsigned int wire, float drift) const
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
// -----------------MlpModelInterface--------------------
// ------------------------------------------------------

nnet::MlpModelInterface::MlpModelInterface(const char* xmlFileName) :
	m(xmlFileName)
{
	mf::LogInfo("MlpModelInterface") << "MLP model loaded.";
}
// ------------------------------------------------------

bool nnet::MlpModelInterface::Run(std::vector< std::vector<float> > const & inp2d)
{
	auto input = nnet::DataProviderAlg::flattenData2D(inp2d);
	if (input.size() == m.GetInputLength())
	{
		m.Run(input);
		return true;
	}
	else
	{
		mf::LogError("MlpModelInterface") << "Flattened patch size does not match MLP model.";
		return false;
	}
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
	else
	{
		mf::LogError("MlpModelInterface") << "Output index does not match MLP model.";
		return 0.;
	}
}
// ------------------------------------------------------


// ------------------------------------------------------
// ----------------KerasModelInterface-------------------
// ------------------------------------------------------

nnet::KerasModelInterface::KerasModelInterface(const char* modelFileName) :
	m(modelFileName)
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
	fOutput = m.compute_output(sample); // add using reference to input so no need for new/delete

	// anyway time is spent in 2D convolutions, not much to be improved in this simple approach...

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
	else
	{
		mf::LogError("KerasModelInterface") << "Output index does not match Keras model.";
		return 0.;
	}
}
// ------------------------------------------------------


// ------------------------------------------------------
// --------------------PointIdAlg------------------------
// ------------------------------------------------------

nnet::PointIdAlg::PointIdAlg(const fhicl::ParameterSet& pset) : nnet::DataProviderAlg(pset),
	fNNet(0)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

nnet::PointIdAlg::~PointIdAlg(void)
{
	deleteNNet();
}
// ------------------------------------------------------

void nnet::PointIdAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fNNetModelFilePath = p.get< std::string >("NNetModelFile");

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

float nnet::PointIdAlg::predictIdValue(std::vector< art::Ptr<recob::Hit> > const & hits, size_t outIdx) const
{
	double pmin = 1.0e-6, pmax = 1.0 - pmin;
	double log_pmin = log(pmin), log_pmax = log(pmax);
	double totarea = 0.0;
	size_t nhits = 0;

	double resultA = 0.0, resultB = 0.0;

	for (auto const & h : hits)
	{
		unsigned int wire = h->WireID().Wire;
		float drift = h->PeakTime();

		if (!isInsideFiducialRegion(wire, drift)) continue;

		double area = 1.0; //h->SummedADC();

		double voutA = predictIdValue(wire, drift);
		double voutB = 1.0 - voutA;

		if (voutA < pmin) { voutA = log_pmin; voutB = log_pmax; }
		else if (voutA > pmax) { voutA = log_pmax; voutB = log_pmin; }
		else { voutA = log(voutA); voutB = log(voutB); }

		resultA += area * voutA;
		resultB += area * voutB;

		totarea += area;
		nhits++;
	}

	if (nhits)
	{
		resultA = exp(resultA / totarea);
		resultB = exp(resultB / totarea);

		resultA = resultA / (resultA + resultB);
	}
	else resultA = 0.5;

	return (float)resultA;
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

std::vector<float> nnet::PointIdAlg::predictIdVector(std::vector< art::Ptr<recob::Hit> > const & hits) const
{
	std::vector<float> result(NClasses(), 0);
	if (result.empty()) return result;

	double pmin = 1.0e-6, pmax = 1.0 - pmin;
	double log_pmin = log(pmin), log_pmax = log(pmax);
	double totarea = 0.0;
	size_t nhits = 0;

	for (auto const & h : hits)
	{
		unsigned int wire = h->WireID().Wire;
		float drift = h->PeakTime();

		if (!isInsideFiducialRegion(wire, drift)) continue;

		double area = 1.0; // h->SummedADC();

		auto vout = predictIdVector(wire, drift);
		for (size_t i = 0; i < vout.size(); ++i)
		{
			if (vout[i] < pmin) vout[i] = log_pmin;
			else if (vout[i] > pmax) vout[i] = log_pmax;
			else vout[i] = log(vout[i]);

			result[i] += area * vout[i];
		}
		totarea += area;
		nhits++;
	}

	if (nhits)
	{
		double totp = 0.0;
		for (size_t i = 0; i < result.size(); ++i)
		{
			result[i] = exp(result[i] / totarea);
			totp += result[i];
		}
		for (size_t i = 0; i < result.size(); ++i)
		{
			result[i] /= totp;
		}
	}
	else std::fill(result.begin(), result.end(), 1.0 / result.size());

	return result;
}
// ------------------------------------------------------

// ------------------------------------------------------
// ------------------TrainingDataAlg---------------------
// ------------------------------------------------------

nnet::TrainingDataAlg::TrainingDataAlg(const fhicl::ParameterSet& pset) : nnet::DataProviderAlg(pset)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

nnet::TrainingDataAlg::~TrainingDataAlg(void)
{
}
// ------------------------------------------------------

void nnet::TrainingDataAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fWireProducerLabel = p.get< std::string >("WireLabel");
	fSimulationProducerLabel = p.get< std::string >("SimulationLabel");
	fSaveVtxFlags = p.get< bool >("SaveVtxFlags");
}
// ------------------------------------------------------

void nnet::TrainingDataAlg::resizeView(size_t wires, size_t drifts)
{
	nnet::DataProviderAlg::resizeView(wires, drifts);

	fWireDriftEdep.resize(wires);
	for (auto & w : fWireDriftEdep)
	{
		w.resize(fNScaledDrifts);
		std::fill(w.begin(), w.end(), 0.0F);
	}

	fWireDriftPdg.resize(wires);
	for (auto & w : fWireDriftPdg)
	{
		w.resize(fNScaledDrifts);
		std::fill(w.begin(), w.end(), 0);
	}
}
// ------------------------------------------------------

bool nnet::TrainingDataAlg::setWireEdepsAndLabels(
	std::vector<float> const & edeps, std::vector<int> const & pdgs, size_t wireIdx)
{
	if ((wireIdx >= fWireDriftEdep.size()) ||
	    (edeps.size() != pdgs.size()) ||
	    (edeps.size() / fDriftWindow > fNScaledDrifts)) return false;

	auto & wEdep = fWireDriftEdep[wireIdx];
	auto & wPdg = fWireDriftPdg[wireIdx];

	for (size_t i = 0; i < fNScaledDrifts; ++i)
	{
		size_t i0 = i * fDriftWindow;
		size_t i1 = (i + 1) * fDriftWindow;

		int best_pdg = pdgs[i0] & nnet::TrainingDataAlg::kPdgMask;
		int vtx_flags = pdgs[i0] & nnet::TrainingDataAlg::kVtxMask;
		float max_edep = edeps[i0];
		for (size_t k = i0 + 1; k < i1; ++k)
		{
			float ek = edeps[k];
			if (ek > max_edep)
			{
				max_edep = ek;
				best_pdg = pdgs[k] & nnet::TrainingDataAlg::kPdgMask; // remember best matching pdg
			}
			vtx_flags |= pdgs[k] & nnet::TrainingDataAlg::kVtxMask;   // accumulate all vtx flags
		}

		wEdep[i] = max_edep;

		if (fSaveVtxFlags) best_pdg |= vtx_flags;
		wPdg[i] = best_pdg;
	}

	return true;
}
// ------------------------------------------------------

nnet::TrainingDataAlg::WireDrift nnet::TrainingDataAlg::getProjection(double x, double y, double z, unsigned int view) const
{
	nnet::TrainingDataAlg::WireDrift wd;
	wd.Wire = 0; wd.Drift = 0; wd.TPC = -1;

	double vtx[3] = {x, y, z};
	if (fGeometry->FindTPCAtPosition(vtx).isValid)
	{
		unsigned int cryo = fGeometry->FindCryostatAtPosition(vtx);
		unsigned int tpc = fGeometry->FindTPCAtPosition(vtx).TPC;

		wd.Wire = fGeometry->NearestWire(vtx, view, tpc, cryo);
		wd.Drift = fDetProp->ConvertXToTicks(x, view, tpc, cryo);
		wd.TPC = tpc;
	}
	return wd;
}
// ------------------------------------------------------

void nnet::TrainingDataAlg::collectVtxFlags(
	std::map< size_t, std::map< int, int > > & wireToDriftToVtxFlags,
	const std::map< int, const simb::MCParticle* > & particleMap,
	unsigned int view) const
{
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

			case 13:   // mu+/-
				if ((particle.EndProcess() == "FastScintillation")) // potential decay at rest
				{
					unsigned int nSec = particle.NumberDaughters();
					
					for (size_t d = 0; d < nSec; ++d)
					{
						auto d_search = particleMap.find(particle.Daughter(d));
						if (d_search != particleMap.end())
						{
							auto const & daughter = *((*d_search).second);
							int d_pdg = abs(daughter.PdgCode());
							if (d_pdg == 11)
							{
								//std::cout << "---> mu decay to electron" << std::endl;
								flagsEnd = nnet::TrainingDataAlg::kDecay;
								break;
							}
						}
					}
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
		if (flagsStart != nnet::TrainingDataAlg::kNone)
		{
			auto wd = getProjection(particle.Vx(), particle.Vy(), particle.Vz(), view);
			
			if (wd.TPC == (int)fTPC)
			{
				wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= flagsStart;
				// std::cout << "---> flagsStart:" << flagsStart << " view:" << view << " wire:" << wd.Wire << " drift:" << wd.Drift << std::endl;
			}
			// else std::cout << "---> not in current TPC" << std::endl;
		}
		if (flagsEnd != nnet::TrainingDataAlg::kNone)
		{
			auto wd = getProjection(particle.EndX(), particle.EndY(), particle.EndZ(), view);
			if (wd.TPC == (int)fTPC)
			{
				wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= flagsEnd;
				// std::cout << "---> flagsEnd:" << flagsEnd << " view:" << view << " wire:" << wd.Wire << " drift:" << wd.Drift << std::endl;
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

	art::ValidHandle< std::vector<simb::MCParticle> > particleHandle
		= event.getValidHandle< std::vector<simb::MCParticle> >(fSimulationProducerLabel);

	art::ValidHandle< std::vector<sim::SimChannel> > simChannelHandle
		= event.getValidHandle< std::vector<sim::SimChannel> >(fSimulationProducerLabel);

    std::map< int, const simb::MCParticle* > particleMap;
	for (auto const & particle : (*particleHandle))
    {
		particleMap[particle.TrackId()] = &particle;
	}

	std::map< size_t, std::map< int, int > > wireToDriftToVtxFlags;
	if (fSaveVtxFlags) collectVtxFlags(wireToDriftToVtxFlags, particleMap, view);

	std::map< int, int > trackToPDG;
    for (size_t widx = 0; widx < fNWires; ++widx)
	{
		auto wireChannelNumber = fWireChannels[widx];
		if (wireChannelNumber == raw::InvalidChannelID) continue;

		std::vector< float > labels_deposit(fNDrifts, 0);  // full-drift-length buffers
		std::vector< int > labels_pdg(fNDrifts, 0);

		std::map< int, std::map< int, double > > timeToTrackToCharge;
		
		for (auto const & channel : (*simChannelHandle))
		{
			auto simChannelNumber = channel.Channel();

			if (simChannelNumber != wireChannelNumber) continue;		

			auto const & timeSlices = channel.TDCIDEMap();
			for (auto const & timeSlice : timeSlices)
			{
				int time = timeSlice.first;

				auto const & energyDeposits = timeSlice.second;
				for (auto const & energyDeposit : energyDeposits)
				{
					int pdg = 0;
					int tid = energyDeposit.trackID;
					if (tid < 0) { pdg = 11; tid = -tid; } // negative tid means it is EM activity

					auto search = particleMap.find(tid);
					if (search == particleMap.end())
					{
						mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND";
						continue;
					}

					auto const & particle = *((*search).second);
					if (!pdg) pdg = particle.PdgCode(); // not EM activity so read what PDG it is

					trackToPDG[energyDeposit.trackID] = abs(pdg);

					double energy = energyDeposit.numElectrons * electronsToGeV;
					timeToTrackToCharge[time][energyDeposit.trackID] += energy;

	      		} // loop over energy deposits
      		} // loop over time slices
      	} // for each SimChannel

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
			
			if (ttc.first < (int)labels_deposit.size())
			{
				labels_deposit[ttc.first] = max_deposit;
				labels_pdg[ttc.first]     = max_pdg & 0xFFFF;
			}
		}

		for (auto const & drift_flags : wireToDriftToVtxFlags[widx])
		{
			int drift = drift_flags.first, flags = drift_flags.second;
			if ((drift >= 0) && (drift < (int)labels_pdg.size()))
			{
				labels_pdg[drift] |= flags;
			}
		}

		setWireEdepsAndLabels(labels_deposit, labels_pdg, widx);

	} // for each Wire

	return true;
}
// ------------------------------------------------------


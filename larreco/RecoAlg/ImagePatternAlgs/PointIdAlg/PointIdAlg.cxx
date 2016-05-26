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

#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "SimulationBase/MCParticle.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardata/RecoBase/Wire.h"


#include "messagefacility/MessageLogger/MessageLogger.h"

nnet::DataProviderAlg::DataProviderAlg(const fhicl::ParameterSet& pset) :
	fCryo(9999), fTPC(9999), fView(9999),
	fNWires(0), fNDrifts(0), fNScaledDrifts(0),
	fDriftWindow(10), fPatchSize(32),
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
	fWireProducerLabel = p.get< std::string >("WireLabel");

	resizePatch();
}
// ------------------------------------------------------

void nnet::DataProviderAlg::resizePatch(void)
{
	fWireDriftPatch.resize(fPatchSize);
	for (auto & r : fWireDriftPatch) r.resize(fPatchSize);
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
bool nnet::DataProviderAlg::setWireData(std::vector<float> const & adc, size_t wireIdx)
{
	if ((wireIdx >= fWireDriftData.size()) ||
	    (adc.size() / fDriftWindow > fNScaledDrifts)) return false;

	auto & wData = fWireDriftData[wireIdx];

	for (size_t i = 0; i < fNScaledDrifts; ++i)
	{
		size_t i0 = i * fDriftWindow;
		size_t i1 = (i + 1) * fDriftWindow;

		float max_adc = scaleAdcSample(adc[i0] * fCalorimetryAlg.LifetimeCorrection(i0));
		for (size_t k = i0 + 1; k < i1; ++k)
		{
			float ak = scaleAdcSample(adc[k] * fCalorimetryAlg.LifetimeCorrection(k));
			if (ak > max_adc) max_adc = ak;
		}

		wData[i] = max_adc;
	}

	return true;
}
// ------------------------------------------------------

bool nnet::DataProviderAlg::setWireDriftData(const art::Event& event,
	unsigned int view, unsigned int tpc, unsigned int cryo)
{
	fCryo = cryo; fTPC = tpc; fView = view;

	art::ValidHandle< std::vector<recob::Wire> > wireHandle
		= event.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);

	size_t nwires = fGeometry->Nwires(view, tpc, cryo);
	size_t ndrifts = fDetProp->NumberTimeSamples();

	resizeView(nwires, ndrifts);

    for (auto const & wire : *wireHandle)
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

	int halfSize = fPatchSize / 2;

	int w0 = fCurrentWireIdx - halfSize;
	int w1 = fCurrentWireIdx + halfSize;

	int d0 = fCurrentScaledDrift - halfSize;
	int d1 = fCurrentScaledDrift + halfSize;

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

std::vector<float> nnet::DataProviderAlg::patchData1D(void) const
{
	std::vector<float> flat;
	flat.resize(fPatchSize * fPatchSize);

	for (size_t w = 0, i = 0; w < fWireDriftPatch.size(); ++w)
	{
		auto const & wire = fWireDriftPatch[w];
		for (size_t d = 0; d < wire.size(); ++d, ++i)
			flat[i] = wire[d];
	}

	return flat;
}
// ------------------------------------------------------

bool nnet::DataProviderAlg::isInsideFiducialRegion(unsigned int wire, float drift) const
{
	size_t halfPatch = fPatchSize / 8; // fPatchSize/2 will make patch always completely filled
	size_t scaledDrift = (size_t)(drift / fDriftWindow);
	if ((wire >= halfPatch) && (wire < fNWires - halfPatch) &&
	    (scaledDrift >= halfPatch) && (scaledDrift < fNScaledDrifts - halfPatch)) return true;
	else return false;
}
// ------------------------------------------------------

// ------------------------------------------------------
// -----------------MlpModelInterface--------------------
// ------------------------------------------------------

nnet::MlpModelInterface::MlpModelInterface(const char* xmlFileName) :
	m(xmlFileName)
{
}
// ------------------------------------------------------


// ------------------------------------------------------
// ----------------KerasModelInterface-------------------
// ------------------------------------------------------

nnet::KerasModelInterface::KerasModelInterface(const char* modelFileName) :
	m(modelFileName)
{
}
// ------------------------------------------------------

// ------------------------------------------------------
// --------------------PointIdAlg------------------------
// ------------------------------------------------------

nnet::PointIdAlg::PointIdAlg(const fhicl::ParameterSet& pset) : nnet::DataProviderAlg(pset),
	fMLP(0), fCNN(0)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

nnet::PointIdAlg::~PointIdAlg(void)
{
	deleteMLP();
	deleteCNN();
}
// ------------------------------------------------------

void nnet::PointIdAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fCalorimetryAlg.reconfigure(p.get< fhicl::ParameterSet >("CalorimetryAlg"));
	fWireProducerLabel = p.get< std::string >("WireLabel");
	fNNetModelFilePath = p.get< std::string >("NNetModelFile");

	deleteMLP();
	deleteCNN();

	// decide if mlp or cnn...
	// read nnet model and weights, set patch and scalling sizes...
	// ... ...

	fMLP = new nnet::NNReader(fNNetModelFilePath.c_str());
	fDriftWindow = 10; // should be in nnet xml?
	fPatchSize = 32; // derive from nnet input size?


	resizePatch();

}
// ------------------------------------------------------

size_t nnet::PointIdAlg::NClasses(void) const
{
	if (fMLP) return fMLP->GetOutputLength();
	else if (fCNN)
	{
		return 0;
	}
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

	if (fMLP)
	{
		auto input = patchData1D();
		if ((input.size() == fMLP->GetInputLength()) &&
		    ((int)outIdx < fMLP->GetOutputLength()))
		{
			fMLP->Run(input);
			result = fMLP->GetOneOutput(outIdx);
		}
		else mf::LogError("PointIdAlg") << "Patch size or output index does not match MLP model.";
	}
	else if (fCNN)
	{
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

	if (fMLP)
	{
		auto input = patchData1D();
		if (input.size() == fMLP->GetInputLength())
		{
			fMLP->Run(input);
			for (size_t o = 0; o < result.size(); ++o)
				result[o] = fMLP->GetOneOutput(o);
		}
		else mf::LogError("PointIdAlg") << "Patch size does not match MLP model.";
	}
	else if (fCNN)
	{
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
	fCalorimetryAlg.reconfigure(p.get< fhicl::ParameterSet >("CalorimetryAlg"));
	fWireProducerLabel = p.get< std::string >("WireLabel");
	fSimulationProducerLabel = p.get< std::string >("SimulationLabel");

	fDriftWindow = p.get< unsigned int >("DriftWindow");
	fPatchSize = p.get< unsigned int >("PatchSize");

	resizePatch();
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

		int best_pdg = 0;
		float max_edep = edeps[i0];
		for (size_t k = i0 + 1; k < i1; ++k)
		{
			float ek = edeps[k];
			if (ek > max_edep)
			{
				max_edep = ek;
				best_pdg = pdgs[k];
			}
		}

		wEdep[i] = max_edep;
		wPdg[i] = best_pdg;
	}

	return true;
}
// ------------------------------------------------------

bool nnet::TrainingDataAlg::setEventData(const art::Event& event,
	unsigned int view, unsigned int tpc, unsigned int cryo)
{
	if (!setWireDriftData(event, view, tpc, cryo))
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

					const simb::MCParticle& particle = *((*search).second);
					if (!pdg) pdg = particle.PdgCode(); // not EM activity so read what PDG it is

					trackToPDG[energyDeposit.trackID] = abs(pdg);

					double energy = energyDeposit.numElectrons * electronsToGeV;
					timeToTrackToCharge[time][energyDeposit.trackID] += energy;

	      		} // loop over energy deposits
      		} // loop over time slices
      	} // for each SimChannel

		for (auto const & tttc : timeToTrackToCharge)
		{
			float max_deposit = 0.0;
			int max_pdg = 0;		
			for (auto const & tc : tttc.second) {

				if( tc.second > max_deposit ) 
				{
					max_deposit = tc.second;
					max_pdg = trackToPDG[tc.first];
				}			
			}
			labels_deposit[tttc.first] = max_deposit;
			labels_pdg[tttc.first] 	   = max_pdg;
		}

		setWireEdepsAndLabels(labels_deposit, labels_pdg, widx);

	} // for each Wire

	return true;
}
// ------------------------------------------------------


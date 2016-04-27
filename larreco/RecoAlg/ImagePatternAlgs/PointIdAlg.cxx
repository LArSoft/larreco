////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski (pplonski86@gmail.com) and R.Sulej (Robert.Sulej@cern.ch), May 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardata/RecoBase/Wire.h"


#include "messagefacility/MessageLogger/MessageLogger.h"

nnet::PointIdAlg::PointIdAlg(const fhicl::ParameterSet& pset) :
	fCryo(9999), fTPC(9999), fView(9999),
	fNWires(0), fNScaledDrifts(0),
	fDriftWindow(10), fPatchSize(32),
	fCurrentWireIdx(99999), fCurrentScaledDrift(99999),
	fCalorimetryAlg(pset.get< fhicl::ParameterSet >("CalorimetryAlg")),
	fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
	fGeometry = &*(art::ServiceHandle<geo::Geometry>());

	this->reconfigure(pset); 
}
// ------------------------------------------------------

nnet::PointIdAlg::~PointIdAlg(void)
{
}
// ------------------------------------------------------

void nnet::PointIdAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fCalorimetryAlg.reconfigure(p.get< fhicl::ParameterSet >("CalorimetryAlg"));
	fWireProducerLabel = p.get< std::string >("WireLabel");
	fNNetModelFilePath = p.get< std::string >("NNetModelFile");

	if (fNNetModelFilePath == "training")
	{
		fDriftWindow = p.get< unsigned int >("DriftWindow");
		fPatchSize = p.get< unsigned int >("PatchSize");
	}
	else // read nnet model and weights
	{
		// ... ...
	}

	resizePatch();

}
// ------------------------------------------------------

void nnet::PointIdAlg::resizePatch(void)
{
	fWireDriftPatch.resize(fPatchSize);
	for (auto & r : fWireDriftData) r.resize(fPatchSize);
}
// ------------------------------------------------------

void nnet::PointIdAlg::resizeView(size_t wires, size_t drifts)
{
	fNWires = wires; fNScaledDrifts = drifts / fDriftWindow;

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

bool nnet::PointIdAlg::setWireData(std::vector<float> const & adc, size_t wireIdx)
{
	if ((wireIdx >= fWireDriftData.size()) ||
	    (adc.size() / fDriftWindow > fNScaledDrifts)) return false;

	auto & wData = fWireDriftData[wireIdx];

	for (size_t i = 0; i < fNScaledDrifts; ++i)
	{
		size_t i0 = i * fDriftWindow;
		size_t i1 = (i + 1) * fDriftWindow;

		float max_adc = adc[i0] * fCalorimetryAlg.LifetimeCorrection(i0);
		for (size_t k = i0 + 1; k < i1; ++k)
		{
			float ak = adc[k] * fCalorimetryAlg.LifetimeCorrection(k);
			if (ak > max_adc) max_adc = ak;
		}

		wData[i] = max_adc;
	}

	return true;
}
// ------------------------------------------------------

bool nnet::PointIdAlg::setWireDriftData(const art::Event& event,
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
				mf::LogError("PointIdAlg") << "ADC vector sizes not match.";
				return false;
			}

			if (!setWireData(adc, w_idx))
			{
				mf::LogError("PointIdAlg") << "Wire data not set.";
				return false;
			}

			fWireChannels[w_idx] = wireChannelNumber;
		}
	}
	return true;
}
// ------------------------------------------------------

bool nnet::PointIdAlg::bufferPatch(size_t wire, size_t drift) const
{
	if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == drift / fDriftWindow))
		return true; // still within the current position

	fCurrentWireIdx = wire;
	fCurrentScaledDrift = drift / fDriftWindow;

	int halfSize = fDriftWindow / 2;

	int w0 = fCurrentWireIdx - halfSize;
	if (w0 < 0) w0 = 0;

	int w1 = fCurrentWireIdx + halfSize;
	if (w1 > (int)fNWires) w1 = fNWires;

	int d0 = fCurrentScaledDrift - halfSize;
	if (d0 < 0) d0 = 0;

	int d1 = fCurrentScaledDrift + halfSize;
	if (d1 > (int)fNScaledDrifts) d1 = fNScaledDrifts;

	for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch)
	{
		auto & src = fWireDriftData[w];
		auto & dst = fWireDriftPatch[wpatch];

		std::fill(dst.begin(), dst.end(), 0);

		for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch)
		{
			dst[dpatch] = src[d];
		}
	}

	return true;
}
// ------------------------------------------------------

std::vector<float> nnet::PointIdAlg::patchData1D(void) const
{
	std::vector<float> flat;
	flat.reserve(fPatchSize * fPatchSize);

	for (size_t w = 0, i = 0; w < fWireDriftPatch.size(); ++w)
	{
		auto const & wire = fWireDriftPatch[w];
		for (size_t d = 0; d < wire.size(); ++d, ++i)
			flat[i] = wire[d];
	}

	return flat;
}
// ------------------------------------------------------

float nnet::PointIdAlg::predictIdValue(unsigned int wire, float drift) const
{
	float result = 0.;

	return result;
}
// ------------------------------------------------------

std::vector<float> nnet::PointIdAlg::predictIdVector(unsigned int wire, float drift) const
{
	std::vector<float> result;

	return result;
}
// ------------------------------------------------------


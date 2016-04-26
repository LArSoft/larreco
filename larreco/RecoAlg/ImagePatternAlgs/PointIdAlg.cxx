////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski (pplonski..@gmail.com) and R.Sulej (Robert.Sulej@cern.ch), May 2016
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

nnet::PointIdAlg::PointIdAlg(const fhicl::ParameterSet& pset)
  : fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

nnet::PointIdAlg::~PointIdAlg(void)
{
}
// ------------------------------------------------------

void nnet::PointIdAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fDriftWindow = 10;        // <----- read it from the nnet model file

	size_t patchSize = 32;    // <----- read it from the nnet model file
	resizePatch(patchSize);

}
// ------------------------------------------------------

void nnet::PointIdAlg::resizeView(size_t wires, size_t drifts)
{
	fWireDriftData.resize(wires);
	for (auto & w : fWireDriftData) w.resize(drifts);
}
// ------------------------------------------------------

void nnet::PointIdAlg::resizePatch(size_t size)
{
	fWireDriftPatch.resize(size);
	for (auto & row : fWireDriftData) row.resize(size);
}

void nnet::PointIdAlg::setWireDriftData(unsigned int view, unsigned int tpc, unsigned int cryo)
{
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


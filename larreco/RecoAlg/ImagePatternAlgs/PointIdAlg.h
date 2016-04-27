////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski (pplonski..@gmail.com) and R.Sulej (Robert.Sulej@cern.ch), May 2016
//
// Point Identification Algorithm
//
//      Run CNN or MLP trained to classify a point in 2D projection. Various features can be
//      recognized, depending on the net model/weights used.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PointIdAlg_h
#define PointIdAlg_h

// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/AnalysisAlg/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT & C++
#include <memory>

namespace nnet
{
	class PointIdAlg;
}

class nnet::PointIdAlg
{
public:

	PointIdAlg(const fhicl::ParameterSet& pset);
	virtual ~PointIdAlg(void);

	void reconfigure(const fhicl::ParameterSet& p);  // read-in nnet, setup patch buffer, ...

	bool setWireDriftData(const art::Event& event,   // once per view: setup ADC buffer, collect & downscale ADC's
		unsigned int view, unsigned int tpc, unsigned int cryo);

	std::vector<float> const & wireData(size_t widx) const { return fWireDriftData[widx]; }

	std::vector< std::vector<float> > const & patchData2D(void) const { return fWireDriftPatch; }
	std::vector<float> patchData1D(void) const;  // flat vector made of the patch data, wire after wire

	unsigned int Cryo(void) const { return fCryo; }
	unsigned int TPC(void) const { return fTPC; }
	unsigned int View(void) const { return fView; }

	unsigned int NWires(void) const { return fNWires; }
	unsigned int NScaledDrifts(void) const { return fNScaledDrifts; }

	float predictIdValue(unsigned int wire, float drift) const;  // calculate single-value prediction (2-class probability) for [wire, drift] point
	std::vector<float> predictIdVector(unsigned int wire, float drift) const;  // calculate multi-class probabilities for [wire, drift] point

private:
	unsigned int fCryo, fTPC, fView;
	unsigned int fNWires, fNScaledDrifts;

	std::vector< raw::ChannelID_t > fWireChannels;              // wire channels (may need this connection...), InvalidChannelID if not used
	std::vector< std::vector<float> > fWireDriftData;           // 2D data for entire projection, drifts scaled down
	mutable std::vector< std::vector<float> > fWireDriftPatch;  // placeholder for patch around identified point

	size_t fDriftWindow, fPatchSize;

	mutable size_t fCurrentWireIdx, fCurrentScaledDrift;
	bool bufferPatch(size_t wire, size_t drift) const;

	bool setWireData(std::vector<float> const & adc, size_t wireIdx);

	void resizeView(size_t wires, size_t drifts);
	void resizePatch(void);

	std::string fWireProducerLabel;
	std::string fNNetModelFilePath;

	// Calorimetry needed to equalize ADC amplitude along drift:
	calo::CalorimetryAlg  fCalorimetryAlg;

	// Geometry and detector properties:
	geo::GeometryCore const* fGeometry;
	detinfo::DetectorProperties const* fDetProp;
};

#endif

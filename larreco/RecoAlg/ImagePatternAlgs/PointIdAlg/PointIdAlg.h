////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski (pplonski86@gmail.com) and R.Sulej (Robert.Sulej@cern.ch), May 2016
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

#include "larreco/RecoAlg/ImagePatternAlgs/MLP/NNReader.h"

// ROOT & C++
#include <memory>

namespace nnet
{
	class DataProviderAlg;
	class PointIdAlg;
	class TrainingDataAlg;
}

/// Base class providing data for training / running classifiers.
class nnet::DataProviderAlg
{
public:

	DataProviderAlg(const fhicl::ParameterSet& pset);
	virtual ~DataProviderAlg(void);

	virtual void reconfigure(const fhicl::ParameterSet& p);  // setup patch buffer, ...

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

protected:
	unsigned int fCryo, fTPC, fView;
	unsigned int fNWires, fNDrifts, fNScaledDrifts;

	std::vector< raw::ChannelID_t > fWireChannels;              // wire channels (may need this connection...), InvalidChannelID if not used
	std::vector< std::vector<float> > fWireDriftData;           // 2D data for entire projection, drifts scaled down
	mutable std::vector< std::vector<float> > fWireDriftPatch;  // placeholder for patch around identified point

	size_t fDriftWindow, fPatchSize;

	mutable size_t fCurrentWireIdx, fCurrentScaledDrift;
	bool bufferPatch(size_t wire, size_t drift) const;

	bool setWireData(std::vector<float> const & adc, size_t wireIdx);

	virtual void resizeView(size_t wires, size_t drifts);
	void resizePatch(void);

	std::string fWireProducerLabel;

	// Calorimetry needed to equalize ADC amplitude along drift:
	calo::CalorimetryAlg  fCalorimetryAlg;

	// Geometry and detector properties:
	geo::GeometryCore const* fGeometry;
	detinfo::DetectorProperties const* fDetProp;
};

class nnet::PointIdAlg : public nnet::DataProviderAlg
{
public:

	PointIdAlg(const fhicl::ParameterSet& pset);
	virtual ~PointIdAlg(void);

	virtual void reconfigure(const fhicl::ParameterSet& p) override;  // read-in nnet

	float predictIdValue(unsigned int wire, float drift) const;  // calculate single-value prediction (2-class probability) for [wire, drift] point
	std::vector<float> predictIdVector(unsigned int wire, float drift) const;  // calculate multi-class probabilities for [wire, drift] point

private:
	std::string fNNetModelFilePath;

	// CNN and MLP models:
	nnet::NNReader* fMLP;
	void deleteMLP(void) { if (fMLP) delete fMLP; fMLP = 0; }

	void* fCNN; // to be defined..
	void deleteCNN(void) { fCNN = 0; } // { if (fCNN) delete fCNN; fCNN = 0; }
};

class nnet::TrainingDataAlg : public nnet::DataProviderAlg
{
public:

	TrainingDataAlg(const fhicl::ParameterSet& pset);
	virtual ~TrainingDataAlg(void);

	virtual void reconfigure(const fhicl::ParameterSet& p) override;

	bool setEventData(const art::Event& event,   // collect & downscale ADC's, charge deposits, pdg labels
		unsigned int view, unsigned int tpc, unsigned int cryo);

	std::vector<float> const & wireEdep(size_t widx) const { return fWireDriftEdep[widx]; }
	std::vector<int> const & wirePdg(size_t widx) const { return fWireDriftPdg[widx]; }

protected:

	virtual void resizeView(size_t wires, size_t drifts) override;

private:

	bool setWireEdepsAndLabels(
		std::vector<float> const & edeps,
		std::vector<int> const & pdgs,
		size_t wireIdx);

	std::vector< std::vector<float> > fWireDriftEdep;
	std::vector< std::vector<int> > fWireDriftPdg;

	std::string fSimulationProducerLabel;
};

#endif

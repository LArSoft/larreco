////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan, May 2016
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
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/AnalysisAlg/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "SimulationBase/MCParticle.h"

#include "larreco/RecoAlg/ImagePatternAlgs/MLP/NNReader.h"
#include "larreco/RecoAlg/ImagePatternAlgs/Keras/keras_model.h"

// ROOT & C++
#include <memory>

namespace nnet
{
	class DataProviderAlg;
	class ModelInterface;
	class MlpModelInterface;
	class KerasModelInterface;
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

	static std::vector<float> flattenData2D(std::vector< std::vector<float> > const & patch);

	std::vector< std::vector<float> > const & patchData2D(void) const { return fWireDriftPatch; }
	std::vector<float> patchData1D(void) const { return flattenData2D(fWireDriftPatch); }  // flat vector made of the patch data, wire after wire

	unsigned int Cryo(void) const { return fCryo; }
	unsigned int TPC(void) const { return fTPC; }
	unsigned int View(void) const { return fView; }

	unsigned int NWires(void) const { return fNWires; }
	unsigned int NScaledDrifts(void) const { return fNScaledDrifts; }

	bool isInsideFiducialRegion(unsigned int wire, float drift) const;

protected:
	unsigned int fCryo, fTPC, fView;
	unsigned int fNWires, fNDrifts, fNScaledDrifts;

	std::vector< raw::ChannelID_t > fWireChannels;              // wire channels (may need this connection...), InvalidChannelID if not used
	std::vector< std::vector<float> > fWireDriftData;           // 2D data for entire projection, drifts scaled down
	mutable std::vector< std::vector<float> > fWireDriftPatch;  // placeholder for patch around identified point

	size_t fDriftWindow, fPatchSize;

	mutable size_t fCurrentWireIdx, fCurrentScaledDrift;
	bool bufferPatch(size_t wire, float drift) const;

	bool setWireData(std::vector<float> const & adc, size_t wireIdx);
	float scaleAdcSample(float val) const;

	virtual void resizeView(size_t wires, size_t drifts);
	void resizePatch(void);

	std::string fWireProducerLabel;

	// Calorimetry needed to equalize ADC amplitude along drift:
	calo::CalorimetryAlg  fCalorimetryAlg;

	// Geometry and detector properties:
	geo::GeometryCore const* fGeometry;
	detinfo::DetectorProperties const* fDetProp;
};
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class nnet::ModelInterface
{
public:
	virtual ~ModelInterface(void) { }

	unsigned int GetInputLength(void) const { return GetInputCols() * GetInputRows(); }
	virtual unsigned int GetInputCols(void) const = 0;
	virtual unsigned int GetInputRows(void) const = 0;
	virtual int GetOutputLength(void) const = 0;

	virtual bool Run(std::vector< std::vector<float> > const & inp2d) = 0;
	virtual std::vector<float> GetAllOutputs(void) const = 0;
	virtual float GetOneOutput(int neuronIndex) const = 0;

protected:
	ModelInterface(void) { }

private:
};
// ------------------------------------------------------

class nnet::MlpModelInterface : public nnet::ModelInterface
{
public:
	MlpModelInterface(const char* xmlFileName);

	virtual unsigned int GetInputRows(void) const { return m.GetInputLength(); }
	virtual unsigned int GetInputCols(void) const { return 1; }
	virtual int GetOutputLength(void) const { return m.GetOutputLength(); }

	virtual bool Run(std::vector< std::vector<float> > const & inp2d);
	virtual float GetOneOutput(int neuronIndex) const;
	virtual std::vector<float> GetAllOutputs(void) const;

private:
	nnet::NNReader m;
};
// ------------------------------------------------------

class nnet::KerasModelInterface : public nnet::ModelInterface
{
public:
	KerasModelInterface(const char* modelFileName);

	virtual unsigned int GetInputRows(void) const { return m.get_input_rows(); }
	virtual unsigned int GetInputCols(void) const { return m.get_input_cols(); }
	virtual int GetOutputLength(void) const { return m.get_output_length(); }

	virtual bool Run(std::vector< std::vector<float> > const & inp2d);
	virtual float GetOneOutput(int neuronIndex) const;
	virtual std::vector<float> GetAllOutputs(void) const;

private:
	std::vector<float> fOutput; // buffer for output values
	keras::KerasModel m; // network model
};
// ------------------------------------------------------

class nnet::PointIdAlg : public nnet::DataProviderAlg
{
public:

	PointIdAlg(const fhicl::ParameterSet& pset);
	virtual ~PointIdAlg(void);

	virtual void reconfigure(const fhicl::ParameterSet& p) override;  // read-in nnet

	size_t NClasses(void) const;

	// calculate single-value prediction (2-class probability) for [wire, drift] point
	float predictIdValue(unsigned int wire, float drift, size_t outIdx = 0) const;
	float predictIdValue(std::vector< art::Ptr<recob::Hit> > const & hits, size_t outIdx = 0) const;

	// calculate multi-class probabilities for [wire, drift] point
	std::vector<float> predictIdVector(unsigned int wire, float drift) const;
	std::vector<float> predictIdVector(std::vector< art::Ptr<recob::Hit> > const & hits) const;

private:
	std::string fNNetModelFilePath;
	nnet::ModelInterface* fNNet;

	void deleteNNet(void) { if (fNNet) delete fNNet; fNNet = 0; }
};
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class nnet::TrainingDataAlg : public nnet::DataProviderAlg
{
public:

	enum EVtxId
	{
		kNone  = 0,
		kNuNC  = 0x0010000, kNuCC = 0x0020000,                      // nu interaction type
		kNuE   = 0x0100000, kNuMu = 0x0200000, kNuTau = 0x0400000,  // nu flavor
		kHadr  = 0x1000000,  // hadronic inelastic scattering
		kPi0   = 0x2000000,  // pi0 produced in this vertex
		kDecay = 0x4000000,  // point of particle decay
		kConv  = 0x8000000   // gamma conversion
	};

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

	struct WireDrift
	{
		size_t Wire;
		int Drift;
		int TPC;
	};

	WireDrift getProjection(double x, double y, double z, unsigned int view) const;

	bool setWireEdepsAndLabels(
		std::vector<float> const & edeps,
		std::vector<int> const & pdgs,
		size_t wireIdx);

	void collectVtxFlags(
		std::map< size_t, std::map< int, int > > & wireToDriftToVtxFlags,
		const std::map< int, const simb::MCParticle* > & particleMap,
		unsigned int view) const;

	std::vector< std::vector<float> > fWireDriftEdep;
	std::vector< std::vector<int> > fWireDriftPdg;

	std::string fSimulationProducerLabel;
	bool fSaveVtxFlags;
};
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

#endif

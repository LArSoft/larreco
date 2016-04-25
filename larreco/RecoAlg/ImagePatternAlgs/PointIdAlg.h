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

	void reconfigure(const fhicl::ParameterSet& p);

	void setWireDriftData(unsigned int view, unsigned int tpc, unsigned int cryo);

	float predictIdValue(unsigned int wire, float drift) const;
	std::vector<float> predictIdVector(unsigned int wire, float drift) const;

private:

	std::vector< std::vector<float> > fWireDriftData;           // 2D data for entire projection, drifts scaled down
	mutable std::vector< std::vector<float> > fWireDriftPatch;  // placeholder for patch around identified point

	void resize(size_t wires, size_t drifts);

	// Geometry and detector properties
	art::ServiceHandle<geo::Geometry> fGeom;
	detinfo::DetectorProperties const* fDetProp;
};

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Projection Matching Algorithm
// see RecoAlg/PMAlg/PmaTrack3D.h for more details.
//
//      Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
//      Based on the "Precise 3D track reco..." AHEP (2013) 260820, with all the tricks that we
//      developed later and with the work for the full-event topology optimization that is still
//      under construction (and porting to LArSoft implementation).
//
//      The algorithm class provides functionality to build the track from selected hits. These
//      can be detailed tracks or just simple segments. The parameters of optimization algorithm,
//      fixed nodes and 3D reference points can be configured here.
//      Please, check the track finding modules to find a way of selecting appropriate clusteres:
//        - PMAlgTrackMaker_module.cc
//
//      Note: not all parameters of the track optimization are available through .fcl, soon there
//      will be more.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ProjectionMatchingAlg_h
#define ProjectionMatchingAlg_h

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "RecoAlg/PMAlg/PmaTrack3D.h"

// ROOT & C++
#include <memory>

namespace pma
{
	class ProjectionMatchingAlg;
}

class pma::ProjectionMatchingAlg
{
public:

	ProjectionMatchingAlg(const fhicl::ParameterSet& pset);
	virtual ~ProjectionMatchingAlg(void);

	void reconfigure(const fhicl::ParameterSet& p);

	double validate(const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		unsigned int testView) const;

	unsigned int testHits(const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits) const
	{ return trk.TestHits(hits, fHitTestingDist2D); }

	pma::Track3D* buildTrack(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2,
		unsigned int testView) const;

	pma::Track3D* extendTrack(
		const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		bool add_nodes) const;

private:

	// Parameters used in the algorithm

	double fOptimizationEps;       // relative change in the obj.function that ends optimization,
	                               // then next nodes are added or track building is finished

	double fFineTuningEps;         // relative change in the obj.function that ends final tuning

	double fTrkValidationDist2D;   // max. distance [cm] used in the track validation in the "third" plane
	double fHitTestingDist2D;      // max. distance [cm] used in testing comp. of hits with the track

	// Geometry and detector properties
	art::ServiceHandle<geo::Geometry> fGeom;
	art::ServiceHandle<util::DetectorProperties> fDetProp;


	static size_t getSegCount(size_t trk_size);
};

#endif

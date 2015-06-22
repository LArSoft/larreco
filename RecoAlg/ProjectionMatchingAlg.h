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
//      under construction.
//
//      The algorithm class provides functionality to build a track from selected hits. These
//      can be detailed tracks or just simple segments (if the number of nodes to add is set to 0).
//      The parameters of optimization algorithm, fixed nodes and 3D reference points can be configured here.
//      Please, check the track finding module to find a way of selecting appropriate clusteres:
//        PMAlgTrackMaker_module.cc
//
//      Note: not all parameters of the track optimization are available through .fcl, soon more
//      will be added, thanks for patience...
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

	/// Calculate the fraction of the track that is closer than fTrkValidationDist2D
	/// to any hit from hits in the testView (a view that was not used to build the track).
	double validate(const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		unsigned int testView) const;

	/// Count the number of hits that are closer than fHitTestingDist2D to the track 2D projection.
	unsigned int testHits(const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits) const
	{ return trk.TestHits(hits, fHitTestingDist2D); }

	/// Build a track from two sets of hits (they should origin from two wire planes).
	/// Number of segments used to create the track depends on the number of hits.
	pma::Track3D* buildTrack(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2) const;

	/// Build a straight segment from two sets of hits (they should origin from two wire planes).
	/// Method is intendet for short tracks or shower initial parts, where only a few hits
	/// per plane are available and there is no chance to see a curvature or any other features.
	pma::Track3D* buildSegment(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2) const;

	/// Build a straight segment from two sets of hits (they should origin from two wire planes),
	/// starting from a given point (like vertex known from another algorithm).
	/// Method is intendet for short tracks or shower initial parts, where only a few hits
	/// per plane are available and there is no chance to see a curvature or any other features.
	pma::Track3D* buildSegment(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2,
		const TVector3& point) const;

	/// Add more hits to an existing track, reoptimize, optionally add more nodes.
	pma::Track3D* extendTrack(
		const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		bool add_nodes) const;

	/// Try to correct track direction of the stopping particle:
	///   dir: kForward  - particle stop is at the end of the track
	///        kBackward - particle stop is at the begining of the track
	/// compares dQ/dx of n hits at each end of the track (default is based on the track length).
	void autoFlip(pma::Track3D& trk,
		pma::Track3D::EDirection dir = Track3D::kForward,
		unsigned int n = 0) const;

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

	// Calculate good number of segments depending on the number of hits.
	static size_t getSegCount(size_t trk_size);
};

#endif

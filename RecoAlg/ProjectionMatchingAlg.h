////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Projection Matching Algorithm
// see RecoAlg/PMAlg/PmaTrack3D.h for more details.
//
//      Build 3D segments and whole tracks by matching an object 2D projections to hits, simultaneously
//      in multiple wire planes. Based on the algorithm first presented in "Precise 3D track reco..."
//      AHEP (2013) 260820, with all the tricks that we developed later and with the work for the full-event
//      topology optimization that is still under construction now (2015).
//
//      The algorithm class provides functionality to build a track from selected hits. These
//      can be detailed tracks or just simple segments (if the number of nodes to add is set to 0).
//      The parameters of optimization algorithm, fixed nodes and 3D reference points can be configured here.
//      Please, check the track making module to find a way of selecting appropriate clusteres:
//        PMAlgTrackMaker_module.cc
//
//      Note: not all parameters of the track optimization are available through .fcl, more are
//      being added still, thanks for patience...
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

	/// Calculate the fraction of trajectory seen by two 2D projections at least; even a
	/// prfect track starts/stops with the hit from one 2D view, then hits from other views
	/// come, which results with the fraction value high, but always < 1.0; wrong cluster
	/// matchings give significantly lower values.
	double twoViewFraction(pma::Track3D& trk) const;

	/// Count the number of hits that are closer than fHitTestingDist2D to the track 2D projection.
	unsigned int testHits(const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits) const
	{ return trk.TestHits(hits, fHitTestingDist2D); }

	/// Build a track from two sets of hits (they should origin from two wire planes);
	/// number of segments used to create the track depends on the number of hits;
	/// optional vmin is the minimum fraction of hits seen from two views.
	pma::Track3D* buildTrack(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2) const;

	/// Build a straight segment from two sets of hits (they should origin from two wire planes);
	/// method is intendet for short tracks or shower initial parts, where only a few hits
	/// per plane are available and there is no chance to see a curvature or any other features.
	pma::Track3D* buildSegment(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2) const;

	/// Build a straight segment from two sets of hits (they should origin from two wire planes),
	/// starting from a given point (like vertex known from another algorithm); method is intendet
	/// for short tracks or shower initial parts, where only a few hits per plane are available
	/// and there is no chance to see a curvature or any other features.
	pma::Track3D* buildSegment(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2,
		const TVector3& point) const;

	/// Add more hits to an existing track, reoptimize, optionally add more nodes.
	pma::Track3D* extendTrack(
		const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		bool add_nodes) const;

	/// Add src to dst as it was its continuation; nodes of src are added to dst after
	/// its own nodes, hits of src are added to hits of dst, then dst is reoptimized.
	void mergeTracks(pma::Track3D& dst, const pma::Track3D& src, bool reopt) const;

	/// Try to correct track direction of the stopping particle:
	///   dir: kForward  - particle stop is at the end of the track;
	///        kBackward - particle stop is at the beginning of the track;
	/// dQ/dx difference has to be above thr to actually flip the track;
	/// compares dQ/dx of n hits at each end of the track (default is based on the track length).
	void autoFlip(pma::Track3D& trk,
		pma::Track3D::EDirection dir = Track3D::kForward,
		 double thr = 0.0, unsigned int n = 0) const;

	/// Intendet to calculate dQ/dx in the initial part of EM cascade; collection
	/// view is used by default, but it works also with other projections.
	double selectInitialHits(pma::Track3D& trk, unsigned int view = geo::kZ) const;

private:

	// Parameters used in the algorithm

	double fOptimizationEps;       // relative change in the obj.function that ends optimization,
	                               // then next nodes are added or track building is finished

	double fFineTuningEps;         // relative change in the obj.function that ends final tuning

	double fTrkValidationDist2D;   // max. distance [cm] used in the track validation in the "third" plane
	double fHitTestingDist2D;      // max. distance [cm] used in testing comp. of hits with the track

	double fMinTwoViewFraction;    // min. length fraction covered with multiple 2D view hits intertwinded with each other

	// Geometry and detector properties
	art::ServiceHandle<geo::Geometry> fGeom;
	art::ServiceHandle<util::DetectorProperties> fDetProp;

	// Calculate good number of segments depending on the number of hits.
	static size_t getSegCount(size_t trk_size);
};

#endif

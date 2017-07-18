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
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ProjectionMatchingAlg_h
#define ProjectionMatchingAlg_h

// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"

#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

// ROOT & C++
#include "TH1F.h"
#include <memory>

namespace pma
{
	class ProjectionMatchingAlg;
}

class pma::ProjectionMatchingAlg
{
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Atom<double> OptimizationEps {
			Name("OptimizationEps"),
			Comment("relative change of the obj.fn which stops optimization after adding a node")
		};

		fhicl::Atom<double> FineTuningEps {
			Name("FineTuningEps"),
			Comment("relative change of the obj.fn which stops fine-tuning of optimized track")
		};

		fhicl::Atom<double> TrkValidationDist2D {
			Name("TrkValidationDist2D"),
			Comment("max. distance [cm] used in the track validation in the third plane")
		};

		fhicl::Atom<double> HitTestingDist2D {
			Name("HitTestingDist2D"),
			Comment("max. distance [cm] used in testing compatibility of hits with the track")
		};

		fhicl::Atom<double> MinTwoViewFraction {
			Name("MinTwoViewFraction"),
			Comment("min. fraction of track length covered with hits from many 2D views intertwinted with each other")
		};

		fhicl::Atom<double> NodeMargin3D {
			Name("NodeMargin3D"),
			Comment("margin in [cm] around TPC for allowed track node positions")
		};

		fhicl::Atom<double> HitWeightU {
			Name("HitWeightU"),
			Comment("weights used for hits in U plane")
		};

		fhicl::Atom<double> HitWeightV {
			Name("HitWeightV"),
			Comment("weights used for hits in V plane")
		};

		fhicl::Atom<double> HitWeightZ {
			Name("HitWeightZ"),
			Comment("weights used for hits in Z plane")
		};
    };

	ProjectionMatchingAlg(const Config& config);
	void reconfigure(const Config& config);

	ProjectionMatchingAlg(const fhicl::ParameterSet& pset) :
		ProjectionMatchingAlg(fhicl::Table<Config>(pset, {})())
	{}

	virtual ~ProjectionMatchingAlg(void) {}

	/// Calculate the fraction of the track that is close to non-empty pixel (above thr value)
	/// in the ADC image of the testView (a view that was not used to build the track).
	double validate_on_adc(const pma::Track3D& trk,
		const nnet::DataProviderAlg & adcImage, float thr) const;

	/// Calculate the fraction of the track that is closer than fTrkValidationDist2D
	/// to any hit from hits in the testView (a view that was not used to build the track).
	/// Creates also histograms of values in pixels for the passing and rejected points on
	/// the track, so the threshold value for the ADC-based calibration can be estimated.
	double validate_on_adc_test(const pma::Track3D& trk,
		const nnet::DataProviderAlg & adcImage,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		TH1F * histoPassing, TH1F * histoRejected) const;

	/// Calculate the fraction of the track that is closer than fTrkValidationDist2D
	/// to any hit from hits in the testView (a view that was not used to build the track).
	double validate(const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		unsigned int testView) const;

	/// Calculate the fraction of the 3D segment that is closer than fTrkValidationDist2D
	/// to any hit from hits in the testView of TPC/Cryo.
	double validate(const TVector3& p0, const TVector3& p1,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		unsigned int testView, unsigned int tpc, unsigned int cryo) const;

	/// Calculate the fraction of trajectory seen by two 2D projections at least; even a
	/// prfect track starts/stops with the hit from one 2D view, then hits from other views
	/// come, which results with the fraction value high, but always < 1.0; wrong cluster
	/// matchings or incomplete tracks give significantly lower values.
	double twoViewFraction(pma::Track3D& trk) const;

	/// Count the number of hits that are closer than eps * fHitTestingDist2D to the track 2D projection.
	unsigned int testHits(const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		double eps = 1.0) const
	{ return trk.TestHits(hits, eps * fHitTestingDist2D); }

	/// Test if hits at the track endpoinds do not stick out of TPC which they belong to.
	/// Here one can implement some configurable margin if needed for real data imeprfections.
	bool isContained(const pma::Track3D& trk, float margin = 0.0F) const
	{
		return (trk.FirstElement()->SameTPC(trk.front()->Point3D(), margin) &&
			trk.LastElement()->SameTPC(trk.back()->Point3D(), margin));
	}

	/// Build a track from two sets of hits from single TPC, hits should origin from at least two
    /// wire planes; number of segments used to create the track depends on the number of hits.
	pma::Track3D* buildTrack(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2 = std::vector< art::Ptr<recob::Hit> >()) const;

	/// Build a track from sets of hits, multiple TPCs are OK (like taken from PFParticles),
    /// as far as hits origin from at least two wire planes.
	pma::Track3D* buildMultiTPCTrack(const std::vector< art::Ptr<recob::Hit> >& hits) const;

	/// Build a shower segment from sets of hits and attached to the provided vertex.
	pma::Track3D* buildShowerSeg(
		const std::vector< art::Ptr<recob::Hit> >& hits, 
		const pma::Vector3D & vtx) const;

	/// Build a straight segment from two sets of hits (they should origin from two wire planes);
	/// method is intendet for short tracks or shower initial parts, where only a few hits
	/// per plane are available and there is no chance to see a curvature or any other features.
	pma::Track3D* buildSegment(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2 = std::vector< art::Ptr<recob::Hit> >()) const;

	/// Build a straight segment from two sets of hits (they should origin from two wire planes),
	/// starting from a given point (like vertex known from another algorithm); method is intendet
	/// for short tracks or shower initial parts, where only a few hits per plane are available
	/// and there is no chance to see a curvature or any other features.
	pma::Track3D* buildSegment(
		const std::vector< art::Ptr<recob::Hit> >& hits_1,
		const std::vector< art::Ptr<recob::Hit> >& hits_2,
		const TVector3& point) const;

	/// Build a straight segment from set of hits (they should origin from two wire planes at least),
	/// starting from a given point.
	pma::Track3D* buildSegment(
		const std::vector< art::Ptr<recob::Hit> >& hits,
		const TVector3& point) const;

	/// Get rid of small groups of hits around cascades; used to calculate cascade starting direction
	/// using the compact core cluster.
	void FilterOutSmallParts(
		double r2d,
		const std::vector< art::Ptr<recob::Hit> >& hits_in,
		std::vector< art::Ptr<recob::Hit> >& hits_out,
		const TVector2& vtx2d) const;

	void RemoveNotEnabledHits(pma::Track3D& trk) const;

	/// Add more hits to an existing track, reoptimize, optionally add more nodes.
	pma::Track3D* extendTrack(
		const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits,
		bool add_nodes) const;

	/// Add 3D reference points to clean endpoints of a track (both need to be in the same TPC).
	void guideEndpoints(pma::Track3D& trk,
		const std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > >& hits) const;

	/// Add 3D reference points to clean endpoint of a track.
	void guideEndpoints(pma::Track3D& trk, pma::Track3D::ETrackEnd endpoint,
		const std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > >& hits) const;

	std::vector< pma::Hit3D* > trimTrackToVolume(pma::Track3D& trk, TVector3 p0, TVector3 p1) const;

	/// Flip tracks to get second as a continuation of first; returns false if not
	/// possible (tracks in reversed order).
	bool alignTracks(pma::Track3D& first, pma::Track3D& second) const;

	/// Add src to dst as it was its continuation; nodes of src are added to dst after
	/// its own nodes, hits of src are added to hits of dst, then dst is reoptimized.
	void mergeTracks(pma::Track3D& dst, pma::Track3D& src, bool reopt) const;

	/// Try to correct track direction of the stopping particle:
	///   dir: kForward  - particle stop is at the end of the track;
	///        kBackward - particle stop is at the beginning of the track;
	/// dQ/dx difference has to be above thr to actually flip the track;
	/// compares dQ/dx of n hits at each end of the track (default is based on the track length).
	void autoFlip(pma::Track3D& trk,
		pma::Track3D::EDirection dir = Track3D::kForward,
		double thr = 0.0, unsigned int n = 0) const { trk.AutoFlip(dir, thr, n); };

	/// Intendet to calculate dQ/dx in the initial part of EM cascade; collection
	/// view is used by default, but it works also with other projections.
	double selectInitialHits(pma::Track3D& trk, unsigned int view = geo::kZ, unsigned int* nused = 0) const;

private:

	// Helpers for guideEndpoints
	bool chkEndpointHits(int wire, int wdir, double drift_x, int view,
		unsigned int tpc, unsigned int cryo,
		const pma::Track3D& trk,
		const std::vector< art::Ptr<recob::Hit> >& hits) const;
	bool addEndpointRef(pma::Track3D& trk,
		const std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > >& hits,
		std::pair<int, int> const * wires, double const * xPos,
		unsigned int tpc, unsigned int cryo) const;

	// Helpers for FilterOutSmallParts 
	bool GetCloseHits(
		double r2d, 
		const std::vector< art::Ptr<recob::Hit> >& hits_in, 
		std::vector<size_t>& used,
		std::vector< art::Ptr<recob::Hit> >& hits_out) const;

	bool Has(const std::vector<size_t>& v, size_t idx) const;

	// Make segment shorter depending on mse 
	void ShortenSeg(pma::Track3D& trk, const geo::TPCGeo& tpcgeom) const;

	// Control length of the track and number of hits which are still enabled
	bool TestTrk(pma::Track3D& trk, const geo::TPCGeo& tpcgeom) const;

	// Calculate good number of segments depending on the number of hits.
	static size_t getSegCount(size_t trk_size);


	// Parameters used in the algorithm

	double fOptimizationEps;       // relative change in the obj.function that ends optimization,
	                               // then next nodes are added or track building is finished

	double fFineTuningEps;         // relative change in the obj.function that ends final tuning

	double fTrkValidationDist2D;   // max. distance [cm] used in the track validation in the "third" plane
	double fHitTestingDist2D;      // max. distance [cm] used in testing comp. of hits with the track

	double fMinTwoViewFraction;    // min. length fraction covered with multiple 2D view hits intertwinted with each other

	// Geometry and detector properties
	art::ServiceHandle<geo::Geometry> fGeom;
	detinfo::DetectorProperties const* fDetProp;
};

#endif

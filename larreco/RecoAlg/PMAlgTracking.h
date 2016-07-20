///////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTracking
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), June 2016
//
// Single track reconstruction toolkit based on Projection Matching Algorithm. Uses cluster collections
// to find single tracks (modularized version of our original code).
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgTracking_h
#define PMAlgTracking_h


// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlgVertexing.h"

// ROOT & C++
#include <memory>

namespace pma
{
	recob::Track convertFrom(const pma::Track3D& src, unsigned int tidx);

	class PMAlgTrackingBase;
	class PMAlgFitter;
	class PMAlgTracker;
}

class pma::PMAlgTrackingBase
{
public:

	const pma::TrkCandidateColl & Result(void) { return fResult; }

	std::vector< std::pair< TVector3, std::vector< std::pair< size_t, bool > > > >
	getVertices(bool onlyBranching = false) const
	{ return fPMAlgVertexing.getVertices(fResult, onlyBranching); }

	std::vector< std::pair< TVector3, size_t > > getKinks(void) const
	{ return fPMAlgVertexing.getKinks(fResult); }

protected:

	PMAlgTrackingBase(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
		const pma::ProjectionMatchingAlg::Config& pmalgConfig,
		const pma::PMAlgVertexing::Config& pmvtxConfig);
	~PMAlgTrackingBase(void);

	bool has(const std::vector<int> & v, int i) const
	{
		for (auto c : v) { if (c == i) return true; }
		return false;
	}

	void guideEndpoints(void);

	pma::cryo_tpc_view_hitmap fHitMap;

	pma::ProjectionMatchingAlg fProjectionMatchingAlg;
	pma::PMAlgVertexing fPMAlgVertexing;

	pma::TrkCandidateColl fResult;
};


class pma::PMAlgFitter : public pma::PMAlgTrackingBase
{
public:

	PMAlgFitter(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
		const std::vector< recob::Cluster > & clusters,
		const std::vector< recob::PFParticle > & pfparticles,
		const art::FindManyP< recob::Hit > & hitsFromClusters,
		const art::FindManyP< recob::Cluster > & clusFromPfps,
		const art::FindManyP< recob::Vertex > & vtxFromPfps,
		const pma::ProjectionMatchingAlg::Config& pmalgConfig,
		const pma::PMAlgVertexing::Config& pmvtxConfig);

	int build(bool runVertexing,
		const std::vector<int> & trackingOnlyPdg = std::vector<int>(),
		const std::vector<int> & trackingSkipPdg = std::vector<int>());

private:

	void buildTracks(const std::vector<int> & trackingOnlyPdg, const std::vector<int> & trackingSkipPdg);
	void buildShowers(const std::vector<int> & trackingOnlyPdg, const std::vector<int> & trackingSkipPdg);

	std::vector< std::vector< art::Ptr<recob::Hit> > > fCluHits;
	std::map< int, std::vector< art::Ptr<recob::Cluster> > > fPfpClusters;
	std::map< int, pma::Vector3D > fPfpVtx;
	std::map< int, int > fPfpPdgCodes;
};

class pma::PMAlgTracker : public pma::PMAlgTrackingBase
{
public:

	PMAlgTracker(const std::vector< art::Ptr<recob::Hit> > & allhitlist,
		const pma::ProjectionMatchingAlg::Config& pmalgConfig,
		const pma::PMAlgVertexing::Config& pmvtxConfig) : PMAlgTrackingBase(allhitlist, pmalgConfig, pmvtxConfig)
	{}

private:
};

#endif

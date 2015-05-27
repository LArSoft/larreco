////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks with the Projection Matching Algorithm.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
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

#include <memory>

namespace trkf {

class PMAlgTrackMaker : public art::EDProducer {
public:
  explicit PMAlgTrackMaker(fhicl::ParameterSet const & p);

  PMAlgTrackMaker(PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker(PMAlgTrackMaker &&) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker &&) = delete;

  void reconfigure(fhicl::ParameterSet const& p);

  void produce(art::Event & e) override;


private:
  // *** various methods to create tracks from clusters ***
  int fromMaxCluster(const art::Event& evt, std::vector< pma::Track3D* >& result);
  int fromExistingAssocs(const art::Event& evt, std::vector< pma::Track3D* >& result);
  // ------------------------------------------------------


  // ************* some common functionality **************
  int maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindMany< recob::Hit >& fb,
	geo::View_t view, int* tpc);

  recob::Track convertFrom(const pma::Track3D& src);
  size_t collectHits(const pma::Track3D& src, std::vector< recob::SpacePoint > spoints);
  // ------------------------------------------------------


  std::string fCluModuleLabel; // label for input cluster collection
  int fCluMatchingAlg;         // which algorithm for cluster association

};
// ------------------------------------------------------

PMAlgTrackMaker::PMAlgTrackMaker(fhicl::ParameterSet const & p)
{
	this->reconfigure(p);
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
}
// ------------------------------------------------------

void PMAlgTrackMaker::reconfigure(fhicl::ParameterSet const& pset)
{
	fCluModuleLabel = pset.get< std::string >("ClusterModuleLabel");
	fCluMatchingAlg = pset.get< int >("CluMatchingAlg");
}
// ------------------------------------------------------

recob::Track PMAlgTrackMaker::convertFrom(const pma::Track3D& src)
{
	std::vector< TVector3 > xyz, dircos;

	for (size_t i = 0; i < src.Nodes().size() - 1; i++)
		if (src.Nodes()[i]->Point3D() != src.Nodes()[i + 1]->Point3D())
		{
			xyz.push_back(src.Nodes()[i]->Point3D());

			TVector3 dc(src.Nodes()[i + 1]->Point3D());
			dc -= src.Nodes()[i]->Point3D();
			dc *= 1.0 / dc.Mag();
			dircos.push_back(dc);
		}
	xyz.push_back(src.Nodes().back()->Point3D());
	dircos.push_back(dircos.back());

	return recob::Track(xyz, dircos);
}
// ------------------------------------------------------

size_t PMAlgTrackMaker::collectHits(
	const pma::Track3D& src, std::vector< recob::SpacePoint > spoints)
{
	return 0;
}
// ------------------------------------------------------

void PMAlgTrackMaker::produce(art::Event & evt)
{
	std::vector< pma::Track3D* > result;

	int retCode = 0;
	switch (fCluMatchingAlg)
	{
		case 1: retCode = fromMaxCluster(evt, result); break;
		case 2: retCode = fromExistingAssocs(evt, result); break;
	}
	switch (retCode)
	{
		case -2: mf::LogError("Summary") << "problem"; break;
		case -1: mf::LogWarning("Summary") << "no input"; break;
		case  0: mf::LogVerbatim("Summary") << "no tracks done"; break;
		default:
			if (retCode < 0) mf::LogVerbatim("Summary") << "unknown result";
			else mf::LogVerbatim("Summary") << retCode << " tracks ready";
			break;
	}

	if (result.size()) // ok, there is something to save
	{
		std::unique_ptr< std::vector< recob::Track > > tracks(new std::vector< recob::Track >);
		std::unique_ptr< std::vector< recob::SpacePoint > > allsp(new std::vector< recob::SpacePoint >);

		std::unique_ptr< art::Assns< recob::Track, recob::SpacePoint > > trk2sp(new art::Assns< recob::Track, recob::SpacePoint >);
		std::unique_ptr< art::Assns< recob::SpacePoint, recob::Hit > > sp2hit(new art::Assns< recob::SpacePoint, recob::Hit >);

		size_t nHits, spStart = 0, spEnd = 0;
		for (size_t i = 0; i < result.size(); i++)
		{
			tracks->push_back(convertFrom(*(result[i])));

			spStart = allsp->size();
			nHits = collectHits(*(result[i]), *allsp);
			spEnd = allsp->size();

			if (nHits) util::CreateAssn(*this, evt, *tracks, *allsp, *trk2sp, spStart, spEnd);
		}

		evt.put(std::move(tracks));
		evt.put(std::move(allsp));
		evt.put(std::move(trk2sp));
		evt.put(std::move(sp2hit));
	}
}
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromMaxCluster(const art::Event& evt, std::vector< pma::Track3D* >& result)
{
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	if (evt.getByLabel(fCluModuleLabel, cluListHandle))
	{
		art::FindMany< recob::Hit > fb(cluListHandle, evt, fCluModuleLabel);

		int max_clu_tpc = -1;
		int max_clu_idx = maxCluster(cluListHandle, fb, geo::kZ, &max_clu_tpc);
		if (max_clu_idx >= 0)
		{
			//recob::Cluster const& clu = (*cluListHandle)[max_clu_idx];

			std::vector< recob::Hit const* > v_coll, v_ind2, v_ind1;
			fb.get(max_clu_idx, v_coll);

			max_clu_idx = maxCluster(cluListHandle, fb, geo::kV, &max_clu_tpc);
			if (max_clu_idx >= 0) fb.get(max_clu_idx, v_ind2);

			max_clu_idx = maxCluster(cluListHandle, fb, geo::kU, &max_clu_tpc);
			if (max_clu_idx >= 0) fb.get(max_clu_idx, v_ind1);

			pma::Track3D* trk = new pma::Track3D();
			if (v_coll.size() > 10) { trk->AddHits(v_coll); mf::LogVerbatim("PMAlgTrackMaker") << "add coll:" << v_coll.size(); }
			if (v_ind2.size() > 10) { trk->AddHits(v_ind2); mf::LogVerbatim("PMAlgTrackMaker") << "add ind2:" << v_ind2.size(); }
			if (v_ind1.size() > 10) { trk->AddHits(v_ind1); mf::LogVerbatim("PMAlgTrackMaker") << "add ind1:" << v_ind1.size(); }

			mf::LogVerbatim("PMAlgTrackMaker") << "track size: " << trk->size();
			std::vector< int > tpcs = trk->TPCs();
			for (size_t t = 0; t < tpcs.size(); ++t)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "  tpc:" << tpcs[t];
			}
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "  #coll:" << trk->NHits(geo::kZ)
				<< "  #ind2:" << trk->NHits(geo::kV)
				<< "  #ind1:" << trk->NHits(geo::kU);

			trk->Initialize();

			result.push_back(trk);
		}
		else
		{
			mf::LogWarning("PMAlgTrackMaker") << "small clusters only";
			return 0;
		}
	}
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "no clusters";
		return -1;
	}

	return result.size();
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindMany< recob::Hit >& fb,
	geo::View_t view, int* tpc)
{
	int idx = -1, tpc_result = -1;
	size_t min_clu_size = 10, s_max = 0, s;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		std::vector< recob::Hit const* > v;
		fb.get(i, v);

		if ((v.front()->WireID().Plane == view) &&
		    ((*tpc == -1) || ((int)(v.front()->WireID().TPC) == *tpc)))
		{
			s = v.size();
			if ((s > min_clu_size) && (s > s_max))
			{
				tpc_result = v.front()->WireID().TPC;
				s_max = s; idx = i;
			}
		}
	}
	if (*tpc == -1) *tpc = tpc_result;
	return idx;
}
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromExistingAssocs(const art::Event& evt, std::vector< pma::Track3D* >& result)
{
	return 0;
}
// ------------------------------------------------------
// ------------------------------------------------------

DEFINE_ART_MODULE(PMAlgTrackMaker)

} // namespace trkf


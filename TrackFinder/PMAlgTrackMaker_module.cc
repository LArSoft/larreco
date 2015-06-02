////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks using Projection Matching Algorithm,
// see RecoAlg/PMAlg/PmaTrack3D.h for details.
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
    const art::FindManyP< recob::Hit >& fb,
    geo::View_t view, int* tpc);
  int maxCluster(
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fb,
    float tmin, float tmax,
    geo::View_t view, int tpc);

  recob::Track convertFrom(const pma::Track3D& src);
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
	produces< art::Assns<recob::Track, recob::Hit> >();
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
	std::vector< std::vector<double> > dst_dQdx;

	std::map< size_t, std::vector<double> > src_dQdx;
	src.GetRawdEdxSequence(src_dQdx, geo::kZ);

	std::cout << "***************" << std::endl;
	for (auto const& entry : src_dQdx)
	{
		std::cout << entry.first << " " << entry.second[5] << " " << entry.second[6] << std::endl;
	}

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

	return recob::Track(xyz, dircos, dst_dQdx);
}
// ------------------------------------------------------

void PMAlgTrackMaker::produce(art::Event& evt)
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
			else if (retCode == 1) mf::LogVerbatim("Summary") << retCode << " track ready";
			else mf::LogVerbatim("Summary") << retCode << " tracks ready";
			break;
	}

	if (result.size()) // ok, there is something to save
	{
		std::unique_ptr< std::vector< recob::Track > > tracks(new std::vector< recob::Track >);
		std::unique_ptr< std::vector< recob::SpacePoint > > allsp(new std::vector< recob::SpacePoint >);

		std::unique_ptr< art::Assns< recob::Track, recob::Hit > > trk2hit(new art::Assns< recob::Track, recob::Hit >);
		std::unique_ptr< art::Assns< recob::Track, recob::SpacePoint > > trk2sp(new art::Assns< recob::Track, recob::SpacePoint >);
		std::unique_ptr< art::Assns< recob::SpacePoint, recob::Hit > > sp2hit(new art::Assns< recob::SpacePoint, recob::Hit >);

		size_t spStart = 0, spEnd = 0;
		double sp_pos[3], sp_err[6];
		for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

		for (size_t t = 0; t < result.size(); t++)
		{
			tracks->push_back(convertFrom(*(result[t])));

			std::vector< art::Ptr< recob::Hit > > hits2d;
			art::PtrVector< recob::Hit > sp_hits;

			spStart = allsp->size();
			for (size_t h = 0; h < result[t]->size(); h++)
			{
				pma::Hit3D* h3d = (*result[t])[h];
				hits2d.push_back(h3d->Hit2DPtr());

				if ((h == 0) ||
				      (sp_pos[0] != h3d->Point3D().X()) ||
				      (sp_pos[1] != h3d->Point3D().Y()) ||
				      (sp_pos[2] != h3d->Point3D().Z()))
				{
					if (sp_hits.size()) // hits assigned to the previous sp
					{
						util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
						sp_hits.clear();
					}
					sp_pos[0] = h3d->Point3D().X();
					sp_pos[1] = h3d->Point3D().Y();
					sp_pos[2] = h3d->Point3D().Z();
					allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
				}
				sp_hits.push_back(h3d->Hit2DPtr());
			}
			if (sp_hits.size()) // hits assigned to the last sp
			{
				util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
			}
			spEnd = allsp->size();

			if (hits2d.size())
			{
				util::CreateAssn(*this, evt, *tracks, hits2d, *trk2hit);
				util::CreateAssn(*this, evt, *tracks, *allsp, *trk2sp, spStart, spEnd);
			}
		}

		// data prods done, delete all pma::Track3D's
		for (size_t t = 0; t < result.size(); t++) delete result[t];

		evt.put(std::move(tracks));
		evt.put(std::move(allsp));

		evt.put(std::move(trk2hit));
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
		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);

		int max_clu_tpc = -1;
		int max_clu_idx = maxCluster(cluListHandle, fbp, geo::kZ, &max_clu_tpc);
		if (max_clu_idx >= 0)
		{	
			std::vector< art::Ptr<recob::Hit> > v_coll = fbp.at(max_clu_idx);

			float tmax = v_coll.front()->PeakTime(), tmin = v_coll.front()->PeakTime(), t;
			for (size_t j = 0; j < v_coll.size(); ++j)
			{
				t = v_coll[j]->PeakTime();
				if (t > tmax) { tmax = t; }
				if (t < tmin) { tmin = t; }
			}

			std::vector< art::Ptr<recob::Hit> > v_ind2;
			max_clu_idx = maxCluster(cluListHandle, fbp, tmin, tmax, geo::kV, max_clu_tpc);
			if (max_clu_idx >= 0) v_ind2 = fbp.at(max_clu_idx);

			std::vector< art::Ptr<recob::Hit> > v_ind1;
			max_clu_idx = maxCluster(cluListHandle, fbp, tmin, tmax, geo::kU, max_clu_tpc);
			if (max_clu_idx >= 0) v_ind1 = fbp.at(max_clu_idx);

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

			mf::LogVerbatim("PMAlgTrackMaker") << "  initialize trk";
			trk->Initialize();
			mf::LogVerbatim("PMAlgTrackMaker") << "  optimize trk";
			trk->Optimize(2);
			mf::LogVerbatim("PMAlgTrackMaker") << "  sort trk";
			trk->SortHits();

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
	const art::FindManyP< recob::Hit >& fbp,
	geo::View_t view, int* tpc)
{
	int idx = -1, tpc_result = -1;
	size_t min_clu_size = 10, s_max = 0, s;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);

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

int PMAlgTrackMaker::maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	float tmin, float tmax,
	geo::View_t view, int tpc)
{
	int idx = -1;
	size_t min_clu_size = 10, s_max = 0, s;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);

		if ((v.front()->WireID().Plane == view) &&
		    ((int)(v.front()->WireID().TPC) == tpc))
		{
			s = 0;
			for (size_t j = 0; j < v.size(); ++j)
				if ((v[j]->PeakTime() >= tmin) && (v[j]->PeakTime() <= tmax))
					s++;

			if ((s > min_clu_size) && (s > s_max))
			{
				s_max = s; idx = i;
			}
		}
	}
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


////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks using Projection Matching Algorithm,
// see RecoAlg/PMAlg/PmaTrack3D.h for details.
//
// Progress:
//    May-June 2015:  track finding and validating, no attempts to build
//                    complex structures yet
//
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

#include "RecoAlg/PMAlg/Utilities.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"

#include <memory>

namespace trkf {

typedef std::map< size_t, std::vector<double> > dedx_map;
typedef std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > view_hitmap;
typedef std::map< unsigned int, view_hitmap > tpc_view_hitmap;
typedef std::map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

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
  cryo_tpc_view_hitmap c_t_v_hits;
  bool splitHits(const art::Event& evt);

  std::vector<size_t> used_clusters, tried_clusters;
  bool has(const std::vector<size_t>& v, size_t idx)
  {
  	for (auto c : v) if (c == idx) return true;
  	return false;
  }

  int maxCluster(
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    geo::View_t view, int* tpc);
  int maxCluster(
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    float tmin, float tmax, size_t min_clu_size,
    geo::View_t view, int tpc);

  bool validate(const pma::Track3D& trk, unsigned int testView);
  recob::Track convertFrom(const pma::Track3D& src);
  // ------------------------------------------------------

  art::ServiceHandle< geo::Geometry > fGeom;

  std::string fHitModuleLabel; // label for hits collection (used for trk validation)
  std::string fCluModuleLabel; // label for input cluster collection
  int fCluMatchingAlg;         // which algorithm for cluster association
  bool fDebugMode;             // for debugging purposes, off by default

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
	fHitModuleLabel = pset.get< std::string >("HitModuleLabel");
	fCluModuleLabel = pset.get< std::string >("ClusterModuleLabel");
	fCluMatchingAlg = pset.get< int >("CluMatchingAlg");
	fDebugMode = pset.get< bool >("DebugMode");
}
// ------------------------------------------------------

recob::Track PMAlgTrackMaker::convertFrom(const pma::Track3D& src)
{
	std::vector< TVector3 > xyz, dircos;
	std::vector< std::vector<double> > dst_dQdx; // [view][dQ/dx]
	dst_dQdx.push_back(std::vector<double>()); // kU
	dst_dQdx.push_back(std::vector<double>()); // kV
	dst_dQdx.push_back(std::vector<double>()); // kZ

	unsigned int cryo = (unsigned int)src.Cryos().front();
	unsigned int tpc = (unsigned int)src.TPCs().front();

	std::map< unsigned int, dedx_map > src_dQdx;
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kU))
	{
		src_dQdx[geo::kU] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kU], geo::kU);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kV))
	{
		src_dQdx[geo::kV] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kV], geo::kV);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kZ))
	{
		src_dQdx[geo::kZ] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kZ], geo::kZ);
	}

	if (fDebugMode) // trajectory from nodes (for debuging, dQ/dx not saved)
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "DEBUG: save only nodes.";
		for (size_t i = 0; i < src.Nodes().size(); i++)
		{
			xyz.push_back(src.Nodes()[i]->Point3D());

			if (i < src.Nodes().size() - 1)
			{
				TVector3 dc(src.Nodes()[i + 1]->Point3D());
				dc -= src.Nodes()[i]->Point3D();
				dc *= 1.0 / dc.Mag();
				dircos.push_back(dc);
			}
			else dircos.push_back(dircos.back());
		}
		return recob::Track(xyz, dircos);
	}
	else // trajectory from hits (use this!)
	{
		for (size_t i = 0; i < src.size(); i++)
		{
			xyz.push_back(src[i]->Point3D());

			if (i < src.size() - 1)
			{
				TVector3 dc(src[i + 1]->Point3D());
				dc -= src[i]->Point3D();
				dc *= 1.0 / dc.Mag();
				dircos.push_back(dc);
			}
			else dircos.push_back(dircos.back());

			double dQ = 0., dx = 0.;
			dst_dQdx[geo::kU].push_back(0.);
			dst_dQdx[geo::kV].push_back(0.);
			dst_dQdx[geo::kZ].push_back(0.);

			for (auto const& m : src_dQdx)
			{
				auto it = m.second.find(i);
				if (it != m.second.end())
				{
					dQ = it->second[5];
					dx = it->second[6];
					if (dx > 0.) dst_dQdx[m.first][i] = dQ/dx;
					break;
				}
			}
		}
		if (xyz.size() != dircos.size())
		{
			mf::LogError("PMAlgTrackMaker") << "pma::Track3D to recob::Track conversion problem.";
		}
		return recob::Track(xyz, dircos, dst_dQdx);
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::validate(const pma::Track3D& trk, unsigned int testView)
{
	double step = 0.3;
	double max_d = 1.0;
	double d2, max_d2 = max_d * max_d;
	unsigned int nAll = 0, nPassed = 0;

	TVector3 p(trk.front()->Point3D());
	for (size_t i = 0; i < trk.Nodes().size() - 1; i++)
	{
		unsigned int tpc = trk.Nodes()[i]->TPC();
		unsigned int cryo = trk.Nodes()[i]->Cryo();

		std::vector< art::Ptr<recob::Hit> >& hits = c_t_v_hits[cryo][tpc][testView];

		TVector3 vNext(trk.Nodes()[i + 1]->Point3D());
		TVector3 vThis(trk.Nodes()[i]->Point3D());
		TVector3 dc(vNext); dc -= vThis;
		dc *= step / dc.Mag();

		double f = pma::GetSegmentProjVector(p, vThis, vNext);
		while (f < 1.0)
		{
			TVector2 p2d = pma::GetProjectionToPlane(p, testView, tpc, cryo);

			for (const auto& h : hits)
			{
				d2 = pma::Dist2(p2d, pma::WireDriftToCm(h->WireID().Wire, h->PeakTime(), testView, tpc, cryo));
				if (d2 < max_d2) { nPassed++; break; }
			}
			nAll++;

			p += dc; f = pma::GetSegmentProjVector(p, vThis, vNext);
		}
		p = vNext;
	}
	double v = nPassed / (double)nAll;
	mf::LogError("PMAlgTrackMaker") << "validate trk: " << v;

	if (v > 0.8) return true;
	else return false;
}
// ------------------------------------------------------

bool PMAlgTrackMaker::splitHits(const art::Event& evt)
{
	art::Handle< std::vector<recob::Hit> > hitListHandle;
	std::vector< art::Ptr<recob::Hit> > hitlist;
	if (evt.getByLabel(fHitModuleLabel, hitListHandle))
	{
		art::fill_ptr_vector(hitlist, hitListHandle);
		
		unsigned int cryo, tpc, view;
		for (auto const& h : hitlist)
		{
			cryo = h->WireID().Cryostat;
			tpc = h->WireID().TPC;
			view = h->WireID().Plane;

			c_t_v_hits[cryo][tpc][view].push_back(h);
		}
		return true;
	}
	else return false;
}

void PMAlgTrackMaker::produce(art::Event& evt)
{
	std::vector< pma::Track3D* > result;

	if (!splitHits(evt))
	{
		mf::LogError("PMAlgTrackMaker") << "Hits not found in the event.";
		return;
	}

	int retCode = 0;
	switch (fCluMatchingAlg)
	{
		default:
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

		for (auto const& trk : result)
		{
			tracks->push_back(convertFrom(*trk));

			std::vector< art::Ptr< recob::Hit > > hits2d;
			art::PtrVector< recob::Hit > sp_hits;

			spStart = allsp->size();
			for (size_t h = 0; h < trk->size(); h++)
			{
				pma::Hit3D* h3d = (*trk)[h];
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
	size_t minBuildSize = 4;

	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	if (evt.getByLabel(fCluModuleLabel, cluListHandle))
	{
		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);

		tried_clusters.clear();
		used_clusters.clear();  // start matching

		int max_clu_tpc = -1;
		int max_Coll_idx = maxCluster(cluListHandle, fbp, geo::kZ, &max_clu_tpc);
		unsigned int nCollHits = (*cluListHandle)[max_Coll_idx].NHits();

		if ((max_Coll_idx >= 0) && (nCollHits > minBuildSize))
		{
			tried_clusters.clear(); // tried for selected coll cluster

			int idx, max_Ind2_idx, max_Ind1_idx;
			std::vector< int > trk_clusters;

			std::vector< art::Ptr<recob::Hit> > v_ind, v_coll = fbp.at(max_Coll_idx);

			float tmax = v_coll.front()->PeakTime(), tmin = v_coll.front()->PeakTime(), t;
			for (size_t j = 0; j < v_coll.size(); ++j)
			{
				t = v_coll[j]->PeakTime();
				if (t > tmax) { tmax = t; }
				if (t < tmin) { tmin = t; }
			}

			pma::Track3D* trk = new pma::Track3D(); // track candidate
			tried_clusters.push_back(max_Coll_idx);
			trk_clusters.push_back(max_Coll_idx);

			mf::LogVerbatim("PMAlgTrackMaker") << "add coll:" << v_coll.size();
			trk->AddHits(v_coll);

			// start matching here

			max_Ind2_idx = maxCluster(cluListHandle, fbp, tmin, tmax, minBuildSize, geo::kV, max_clu_tpc);
			max_Ind1_idx = maxCluster(cluListHandle, fbp, tmin, tmax, minBuildSize, geo::kU, max_clu_tpc);

			unsigned int nInd2Hits = 0, nInd1Hits = 0;
			if (max_Ind2_idx >= 0) nInd2Hits = (*cluListHandle)[max_Ind2_idx].NHits();
			if (max_Ind1_idx >= 0) nInd1Hits = (*cluListHandle)[max_Ind1_idx].NHits();

			bool try_build = true;
			unsigned int testView = geo::kUnknown;
			if ((nInd2Hits > nInd1Hits) && (nInd2Hits > minBuildSize))
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "add ind2" << nInd2Hits;
				idx = max_Ind2_idx;
				testView = geo::kU;
			}
			else if (nInd1Hits > minBuildSize)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "add ind1:" << nInd1Hits;
				idx = max_Ind1_idx;
				testView = geo::kV;
			}
			else try_build = false;

			if (try_build)
			{
				v_ind = fbp.at(idx);
				tried_clusters.push_back(idx);
				trk_clusters.push_back(idx);
				trk->AddHits(v_ind);

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


				int nNodes, nSegments = nCollHits / (2 * (int)sqrt(nCollHits));
				if (nSegments > 1) nNodes = nSegments - 1;
				else nNodes = 0;

				mf::LogVerbatim("PMAlgTrackMaker") << "  initialize trk";
				trk->Initialize();
				mf::LogVerbatim("PMAlgTrackMaker") << "  optimize trk (" << nSegments << " seg)";
				if (nNodes) trk->Optimize(nNodes, 0.01F);   // build nodes
				double g = trk->Optimize(0, 0.0001F);       // final tuning
				mf::LogVerbatim("PMAlgTrackMaker") << "  done, g = " << g;
				mf::LogVerbatim("PMAlgTrackMaker") << "  sort trk";
				trk->SortHits();

				if (validate(*trk, testView))
				{
					for (auto c : trk_clusters) used_clusters.push_back(c);
					result.push_back(trk);
				}
				else
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "track rejected";
					delete trk;
				}
				// go to next matching
			}
			else
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "no matching clusters";
				delete trk;
			}
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
		if (has(used_clusters, i) || has(tried_clusters, i)) continue;

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
	float tmin, float tmax, size_t min_clu_size,
	geo::View_t view, int tpc)
{
	int idx = -1;
	size_t s_max = 0, s;
	double fraction = 0.0;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		if (has(used_clusters, i) || has(tried_clusters, i)) continue;

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
				s_max = s; idx = i; fraction = s / (double)v.size();
			}
		}
	}
	if (fraction > 0.5) return idx;
	else return -1;
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


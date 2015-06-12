////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks using Projection Matching Algorithm,
// see RecoAlg/ProjectionMatchingAlg.h for details.
//
// The module is intended to loop over clusters wiht various logics. The ProjectionMatchingAlg
// class is used to build / modify / verify tracks using settings that are configured for it.
//
// Progress:
//    May-June 2015:  track finding and validation, quite a conservative iterative merging
//                    of matching clusters and growing tracks, no attempts to consciously
//                    build multi-track structures yet, however:
//                    cosmic tracking works fine since they are sets of independent tracks
//
////////////////////////////////////////////////////////////////////////////////////////////////////

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

#include "RecoAlg/ProjectionMatchingAlg.h"

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
  // *** methods to create tracks from clusters ***

  // loop over all clusters and build as much as possible to find
  int fromMaxCluster(const art::Event& evt, std::vector< pma::Track3D* >& result);

  void fromMaxCluster_tpc(
    std::vector< pma::Track3D* >& result,
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    size_t minBuildSize, unsigned int tpc, unsigned int cryo);

  pma::Track3D* extendTrack(
	const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	double gmax, double vmin, unsigned int testView,
	bool add_nodes);

  int matchCluster(
    const pma::Track3D& trk,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minSize, unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo,
	double fraction);
  // ------------------------------------------------------

  // build tracks from clusters associated by any other module - not yet implemented
  int fromExistingAssocs(const art::Event& evt, std::vector< pma::Track3D* >& result);
  // ------------------------------------------------------


  // ************* some common functionality **************
  cryo_tpc_view_hitmap c_t_v_hits;
  bool sortHits(const art::Event& evt);

  std::vector<size_t> used_clusters;
  std::map< unsigned int, std::vector<size_t> > tried_clusters;
  bool has(const std::vector<size_t>& v, size_t idx)
  {
  	for (auto c : v) if (c == idx) return true;
  	return false;
  }

  int maxCluster(
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);
  int maxCluster(
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    float tmin, float tmax, size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);

  double validate(const pma::Track3D& trk, unsigned int testView);
  recob::Track convertFrom(const pma::Track3D& src);
  // ------------------------------------------------------

  art::ServiceHandle< geo::Geometry > fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

  // ******************** parameters **********************
  std::string fHitModuleLabel; // label for hits collection (used for trk validation)
  std::string fCluModuleLabel; // label for input cluster collection
  int fCluMatchingAlg;         // which algorithm for cluster association
  bool fDebugMode;             // for debugging purposes, off by default

  pma::ProjectionMatchingAlg fProjectionMatchingAlg;
};
// ------------------------------------------------------

PMAlgTrackMaker::PMAlgTrackMaker(fhicl::ParameterSet const & p) :
	fProjectionMatchingAlg(p.get< fhicl::ParameterSet >("ProjectionMatchingAlg"))
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

	fProjectionMatchingAlg.reconfigure(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg"));
}
// ------------------------------------------------------

recob::Track PMAlgTrackMaker::convertFrom(const pma::Track3D& src)
{
	std::vector< TVector3 > xyz, dircos;
	std::vector< std::vector<double> > dst_dQdx; // [view][dQ/dx]
	dst_dQdx.push_back(std::vector<double>()); // kU
	dst_dQdx.push_back(std::vector<double>()); // kV
	dst_dQdx.push_back(std::vector<double>()); // kZ

	unsigned int cryo = src.FrontCryo();
	unsigned int tpc = src.FrontTPC();

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

double PMAlgTrackMaker::validate(const pma::Track3D& trk, unsigned int testView)
{
	if (testView == geo::kUnknown)
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "Skip validation.";
		return 1.0;
	}
	else mf::LogVerbatim("PMAlgTrackMaker") << "Validation in plane: " << testView;

	std::vector< art::Ptr<recob::Hit> >& hits = c_t_v_hits[trk.FrontCryo()][trk.FrontTPC()][testView];
	if (hits.size() > 10)
	{
		return fProjectionMatchingAlg.validate(trk, hits, testView);
	}
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "  too few hits for validation: " << hits.size();
		return 1.0;
	}
}
// ------------------------------------------------------

pma::Track3D* PMAlgTrackMaker::extendTrack(const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	double gmax, double vmin, unsigned int testView,
	bool add_nodes)
{
	pma::Track3D* copy = fProjectionMatchingAlg.extendTrack(trk, hits, add_nodes);
	double g1 = copy->GetObjFunction();
	double v1 = validate(*copy, testView);

	if ((g1 < 0.15) && (g1 <= gmax) && (v1 >= vmin))
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "  track EXTENDED";
		copy->SortHits();
		return copy;
	}
	else
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "  track NOT extended";
		delete copy;
		return 0;
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::sortHits(const art::Event& evt)
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
// ------------------------------------------------------

void PMAlgTrackMaker::produce(art::Event& evt)
{
	std::vector< pma::Track3D* > result;

	if (!sortHits(evt))
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
// ------------------------------------------------------

int PMAlgTrackMaker::fromMaxCluster(const art::Event& evt, std::vector< pma::Track3D* >& result)
{
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	if (evt.getByLabel(fCluModuleLabel, cluListHandle))
	{
		tried_clusters.clear(); used_clusters.clear();

		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);

		size_t minBuildSizeLarge = 20;
		size_t minBuildSizeSmall = 5;

		// find reasonably large parts
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(result, cluListHandle, fbp, minBuildSizeLarge, tpc_iter->TPC, tpc_iter->Cryostat);
		}

		// loop again to find small things
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(result, cluListHandle, fbp, minBuildSizeSmall, tpc_iter->TPC, tpc_iter->Cryostat);
		}
	}
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "no clusters";
		return -1;
	}

	return result.size();
}

int PMAlgTrackMaker::matchCluster(const pma::Track3D& trk,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minSize, unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo,
	double fraction)
{
	double f, fmax = 0.0;
	unsigned int n, max = 0;
	int idx = -1;
	for (size_t i = 0; i < clusters->size(); ++i)
	{
		unsigned int view = (*clusters)[i].View();
		unsigned int nhits = (*clusters)[i].NHits();

		if (has(used_clusters, i) ||                             // don't try already used clusters
		    (view == testView) ||                                // don't use clusters from validation view
		    ((preferedView != geo::kUnknown)&&(view != preferedView)) || // only prefered view if specified
		    (nhits < minSize))                                   // skip small clusters
		    continue;

		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);
		n = fProjectionMatchingAlg.testHits(trk, v);
		f = n / (double)v.size();
		if ((f > fraction) && (n > max))
		{
			max = n; fmax = f; idx = i;
		}
	}

	if (idx >= 0) mf::LogVerbatim("PMAlgTrackMaker") << "max matching hits: " << max << " (" << fmax << ")";
	else mf::LogVerbatim("PMAlgTrackMaker") << "no clusters to extend the track";

	return idx;
}


void PMAlgTrackMaker::fromMaxCluster_tpc(
	std::vector< pma::Track3D* >& result,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minBuildSize, unsigned int tpc, unsigned int cryo)
{
	tried_clusters[geo::kZ].clear();

	size_t minSizeCompl = minBuildSize / 5; // smaller minimum required in complementary views

	int max_Coll_idx = 0;
	while (max_Coll_idx >= 0) // loop over coll
	{
		max_Coll_idx = maxCluster(clusters, fbp, minBuildSize, geo::kZ, tpc, cryo);
		tried_clusters[geo::kZ].push_back(max_Coll_idx);

		if (max_Coll_idx >= 0)
		{
			unsigned int nCollHits = (*clusters)[max_Coll_idx].NHits();
			mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** Z ***  size: " << nCollHits;

			std::vector< art::Ptr<recob::Hit> > v_coll = fbp.at(max_Coll_idx);

			float tmax = fDetProp->ConvertTicksToX(v_coll.front()->PeakTime(), geo::kZ, tpc, cryo);
			float t, tmin = tmax;
			for (size_t j = 1; j < v_coll.size(); ++j)
			{
				t = fDetProp->ConvertTicksToX(v_coll[j]->PeakTime(), geo::kZ, tpc, cryo);
				if (t > tmax) { tmax = t; }
				if (t < tmin) { tmin = t; }
			}

			tried_clusters[geo::kU].clear();
			tried_clusters[geo::kV].clear();

			bool try_build = true;
			while (try_build) // loop over ind
			{
				int idx, max_Ind2_idx, max_Ind1_idx;
				max_Ind2_idx = maxCluster(clusters, fbp, tmin, tmax, minSizeCompl, geo::kV, tpc, cryo);
				max_Ind1_idx = maxCluster(clusters, fbp, tmin, tmax, minSizeCompl, geo::kU, tpc, cryo);

				unsigned int nInd2Hits = 0, nInd1Hits = 0;
				if (max_Ind2_idx >= 0) nInd2Hits = (*clusters)[max_Ind2_idx].NHits();
				if (max_Ind1_idx >= 0) nInd1Hits = (*clusters)[max_Ind1_idx].NHits();

				unsigned int testView = geo::kUnknown;
				if ((nInd2Hits > nInd1Hits) && (nInd2Hits > minSizeCompl))
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** V ***  size: " << nInd2Hits;
					tried_clusters[geo::kV].push_back(max_Ind2_idx);
					idx = max_Ind2_idx;
					testView = geo::kU;
				}
				else if (nInd1Hits > minSizeCompl)
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** U ***  size: " << nInd1Hits;
					tried_clusters[geo::kU].push_back(max_Ind1_idx);
					idx = max_Ind1_idx;
					testView = geo::kV;
				}
				else try_build = false;

				if (try_build)
				{
					if (!fGeom->TPC(tpc, cryo).HasPlane(testView)) testView = geo::kUnknown;

					pma::Track3D* trk = fProjectionMatchingAlg.buildTrack(v_coll, fbp.at(idx), testView);

					double g0 = trk->GetObjFunction();
					double v0 = validate(*trk, testView);
					if ((g0 < 0.1) && (v0 > 0.7)) // good track, try to extend it and then save to results
					{
						used_clusters.push_back(max_Coll_idx);
						used_clusters.push_back(idx);

						size_t minSize = 5;    // min size for clusters matching
						double fraction = 0.4; // min fraction of close hits
						
						idx = 0;
						while (idx >= 0) // try to collect matching clusters, use **any** plane except validation
						{
							idx = matchCluster(*trk, clusters, fbp, minSize, geo::kUnknown, testView, tpc, cryo, fraction);
							if (idx >= 0)
							{
								// try building extended copy:    src,    hits,    max.obj.fn, min.valid, valid.plane, add nodes
								pma::Track3D* copy = extendTrack(*trk, fbp.at(idx), 2.0 * g0,  0.98 * v0,  testView,    true);
								if (copy)
								{
									used_clusters.push_back(idx);
									delete trk; trk = copy;
								}
								else idx = -1;
							}
						}

						mf::LogVerbatim("PMAlgTrackMaker") << "merge clusters from the validation plane";
						fraction = 0.7; // only well matching the existing track

						idx = 0;
						while (idx >= 0)
						{	//                     match clusters from the plane used previously for the validation
							idx = matchCluster(*trk, clusters, fbp, minSize, testView, geo::kUnknown, tpc, cryo, fraction);
							if (idx >= 0)
							{
								// only loose validation using kZ, no new nodes:
								pma::Track3D* copy = extendTrack(*trk, fbp.at(idx), 2.0 * g0, 0.7, geo::kZ, false);
								if (copy)
								{
									used_clusters.push_back(idx);
									delete trk; trk = copy;
								}
								else idx = -1;
							}
						}

						result.push_back(trk);

						break; // found match, end loop over ind

					}
					else
					{
						mf::LogVerbatim("PMAlgTrackMaker") << "track REJECTED: g = " << g0 << "; v = " << v0;
						delete trk;
					}
				}
				else
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "no matching clusters";
				}
			} // end loop over ind
		}
		else
		{
			mf::LogWarning("PMAlgTrackMaker") << "small clusters only";
		}
	} // end loop over coll
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		if (has(used_clusters, i) || has(tried_clusters[view], i)) continue;

		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);

		if ((v.front()->WireID().Plane == view) &&
		    (v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = v.size();
			if ((s > min_clu_size) && (s > s_max))
			{
				s_max = s; idx = i;
			}
		}
	}
	return idx;
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	float tmin, float tmax, size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;
	double fraction = 0.0;
	float x;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		if (has(used_clusters, i) ||
		    has(tried_clusters[view], i) ||
		    ((*clusters)[i].NHits() <  min_clu_size) ||
		    ((*clusters)[i].View() != view)) continue;
		    
		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = 0;
			for (size_t j = 0; j < v.size(); ++j)
			{
				x = fDetProp->ConvertTicksToX(v[j]->PeakTime(), view, tpc, cryo);
				if ((x >= tmin) && (x <= tmax)) s++;
			}

			if (s > s_max)
			{
				s_max = s; idx = i; fraction = s / (double)v.size();
			}
		}
	}
	if (fraction > 0.4) return idx;
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


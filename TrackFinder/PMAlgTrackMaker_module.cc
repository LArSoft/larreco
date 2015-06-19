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
//                    build multi-track structures yet, however cosmic tracking works fine
//                    as they are sets of independent tracks
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

#include "MCCheater/BackTracker.h"
#include "Simulation/ParticleList.h"
#include "SimulationBase/MCParticle.h"

#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"

#include "TTree.h"

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

  void beginJob();

  void reconfigure(fhicl::ParameterSet const& p);

  void produce(art::Event & e) override;


private:
  // *** methods to create tracks from clusters ***

  // loop over all clusters and build as much as possible to find
  // the logic implemented here for sure is not exhaustive, it was
  // checked on long cosmic tracks and low energy stopping tracks
  // and seems to be a good example how to use the algorithm
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

  // display what was used and what is left, just for development and debugging
  void listUsedClusters(art::Handle< std::vector<recob::Cluster> > clusters) const;
  // ------------------------------------------------------

  // build tracks from clusters associated by any other module - not yet implemented
  int fromExistingAssocs(const art::Event& evt, std::vector< pma::Track3D* >& result);
  // ------------------------------------------------------


  // ************* some common functionality **************
  cryo_tpc_view_hitmap c_t_v_hits;
  bool sortHits(const art::Event& evt);

  std::vector< size_t > used_clusters, checked_clusters;
  std::map< unsigned int, std::vector<size_t> > tried_clusters;
  bool has(const std::vector<size_t>& v, size_t idx) const
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

  bool isMcStopping(void) const; // to be moved to the testing module
  // ------------------------------------------------------

  art::ServiceHandle< geo::Geometry > fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

  // ******************* tree output **********************
  int fEvNumber;        // event number
  int fTrkIndex;        // track index in the event, same for all dQ/dx points of the track
  int fPlaneIndex;      // wire plane index of the dQ/dx data point
  int fIsStopping;      // tag tracks of stopping particles
  double fdQdx;         // dQ/dx data point stored for each plane
  double fRange;        // residual range at dQ/dx data point, tracks are auto-flipped
  TTree* fTree;

  // ******************** parameters **********************
  std::string fHitModuleLabel; // label for hits collection (used for trk validation)
  std::string fCluModuleLabel; // label for input cluster collection
  int fCluMatchingAlg;         // which algorithm for cluster association
  bool fAutoFlip_dQdx;         // set the track direction to increasing dQ/dx
  bool fSave_dQdx;             // for debugging purposes, off by default

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

void PMAlgTrackMaker::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	fTree = tfs->make<TTree>("PMAlgTrackMaker", "tracks info");
	fTree->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fTree->Branch("fTrkIndex", &fTrkIndex, "fTrkIndex/I");
	fTree->Branch("fPlaneIndex", &fPlaneIndex, "fPlaneIndex/I");
	fTree->Branch("fIsStopping", &fIsStopping, "fIsStopping/I");
	fTree->Branch("fdQdx", &fdQdx, "fdQdx/D");
	fTree->Branch("fRange", &fRange, "fRange/D");
}

void PMAlgTrackMaker::reconfigure(fhicl::ParameterSet const& pset)
{
	fHitModuleLabel = pset.get< std::string >("HitModuleLabel");
	fCluModuleLabel = pset.get< std::string >("ClusterModuleLabel");
	fCluMatchingAlg = pset.get< int >("CluMatchingAlg");
	fAutoFlip_dQdx = pset.get< bool >("AutoFlip_dQdx");
	fSave_dQdx = pset.get< bool >("Save_dQdx");

	fProjectionMatchingAlg.reconfigure(pset.get< fhicl::ParameterSet >("ProjectionMatchingAlg"));
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

	// trajectory from nodes (for debuging only, dQ/dx not saved)
/*	{
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
	} */

	// trajectory from hits (use this!)
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

			double dQdx;
			for (auto const& m : src_dQdx)
			{
				auto it = m.second.find(i);
				if (it != m.second.end())
				{
					dQ = it->second[5];
					dx = it->second[6];
					if (dx > 0.) dQdx = dQ/dx;
					else dQdx = 0.;

					dst_dQdx[m.first][i] = dQdx;

					break;
				}
			}
		}

		if (fSave_dQdx)
		{
			fIsStopping = (int)isMcStopping();
			for (unsigned int view = 0; view < fGeom->Nviews(); view++)
				if (fGeom->TPC(tpc, cryo).HasPlane(view))
				{
					fPlaneIndex = view;
					for (auto const& data : src_dQdx[view])
						for (size_t i = 0; i < data.second.size() - 1; i++)
						{
							double dQ = data.second[5];
							double dx = data.second[6];
							fRange = data.second[7];
							if (dx > 0.)
							{
								fdQdx = dQ/dx;
								fTree->Fill();
							}
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

bool PMAlgTrackMaker::isMcStopping(void) const
{
	art::ServiceHandle< cheat::BackTracker > bt;
	const sim::ParticleList& plist = bt->ParticleList();
	const simb::MCParticle* particle = plist.Primary(0);

	if (particle)
	{
//		std::cout << "...:SIM:... " << particle->EndProcess()
//			<< " n:" << particle->NumberDaughters()
//			<< " m:" << particle->Mass()
//			<< " E:" << particle->EndE() << std::endl;
		return (particle->NumberDaughters() == 0);
		//return (particle->EndE() - particle->Mass() < 0.001);
	}
	else return false;
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
	fEvNumber = evt.id().event();

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

		fTrkIndex = 0;
		for (auto const& trk : result)
		{
			// flip the track by dQ/dx
			if (fAutoFlip_dQdx) fProjectionMatchingAlg.autoFlip(*trk, pma::Track3D::kBackward);

			tracks->push_back(convertFrom(*trk));
			fTrkIndex++;

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
		tried_clusters.clear();
		checked_clusters.clear();
		used_clusters.clear();

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

		// used for development
		//listUsedClusters(cluListHandle);
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
	checked_clusters.clear();

	size_t minSizeCompl = minBuildSize / 5; // smaller minimum required in complementary views

	int max_first_idx = 0;
	while (max_first_idx >= 0) // loop over clusters, any view, from the largest
	{
		max_first_idx = maxCluster(clusters, fbp, minBuildSize, geo::kUnknown, tpc, cryo); // any view
		if (max_first_idx >= 0)
		{
			const recob::Cluster& clu_first = (*clusters)[max_first_idx];

			geo::View_t first_view = clu_first.View();
			geo::View_t sec_view_a, sec_view_b;
			switch (first_view)
			{
				case geo::kU: sec_view_a = geo::kZ; sec_view_b = geo::kV; break;
				case geo::kV: sec_view_a = geo::kZ; sec_view_b = geo::kU; break;
				case geo::kZ: sec_view_a = geo::kV; sec_view_b = geo::kU; break;
				default:
					mf::LogError("PMAlgTrackMaker") << "Not a 2D view.";
					return;
			}

			tried_clusters[geo::kU].clear();
			tried_clusters[geo::kV].clear();
			tried_clusters[geo::kZ].clear();

			//tried_clusters[first_view].push_back(max_first_idx);
			checked_clusters.push_back(max_first_idx);

			unsigned int nFirstHits = clu_first.NHits();
			mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << first_view << " ***  size: " << nFirstHits;

			std::vector< art::Ptr<recob::Hit> > v_first = fbp.at(max_first_idx);

			float tmax = fDetProp->ConvertTicksToX(v_first.front()->PeakTime(), first_view, tpc, cryo);
			float t, tmin = tmax;
			for (size_t j = 1; j < v_first.size(); ++j)
			{
				t = fDetProp->ConvertTicksToX(v_first[j]->PeakTime(), first_view, tpc, cryo);
				if (t > tmax) { tmax = t; }
				if (t < tmin) { tmin = t; }
			}

			bool try_build = true;
			while (try_build) // loop over ind
			{
				int idx, max_sec_a_idx, max_sec_b_idx;
				max_sec_a_idx = maxCluster(clusters, fbp, tmin, tmax, minSizeCompl, sec_view_a, tpc, cryo);
				max_sec_b_idx = maxCluster(clusters, fbp, tmin, tmax, minSizeCompl, sec_view_b, tpc, cryo);

				unsigned int nSecHitsA = 0, nSecHitsB = 0;
				if (max_sec_a_idx >= 0) nSecHitsA = (*clusters)[max_sec_a_idx].NHits();
				if (max_sec_b_idx >= 0) nSecHitsB = (*clusters)[max_sec_b_idx].NHits();

				unsigned int testView = geo::kUnknown;
				if ((nSecHitsA > nSecHitsB) && (nSecHitsA > minSizeCompl))
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_a << " ***  size: " << nSecHitsA;
					tried_clusters[sec_view_a].push_back(max_sec_a_idx);
					idx = max_sec_a_idx;
					testView = sec_view_b;
				}
				else if (nSecHitsB > minSizeCompl)
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_b << " ***  size: " << nSecHitsB;
					tried_clusters[sec_view_b].push_back(max_sec_b_idx);
					idx = max_sec_b_idx;
					testView = sec_view_a;
				}
				else try_build = false;

				if (try_build)
				{
					if (!fGeom->TPC(tpc, cryo).HasPlane(testView)) testView = geo::kUnknown;

					pma::Track3D* trk = fProjectionMatchingAlg.buildTrack(v_first, fbp.at(idx));

					double g0 = trk->GetObjFunction();
					double v0 = validate(*trk, testView);
					if ((g0 < 0.1) && (v0 > 0.7)) // good track, try to extend it and then save to results
					{
						used_clusters.push_back(max_first_idx);
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
						while ((idx >= 0) && (testView != geo::kUnknown))
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
	} // end loop over clusters, any view, from the largest
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
		if (has(used_clusters, i) ||
		    has(checked_clusters, i) ||
		    has(tried_clusters[view], i) ||
		   ((view != geo::kUnknown) && ((*clusters)[i].View() != view)))
		continue;

		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);

		if ((v.front()->WireID().TPC == tpc) &&
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
		    has(checked_clusters, i) ||
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

void PMAlgTrackMaker::listUsedClusters(art::Handle< std::vector<recob::Cluster> > clusters) const
{
	mf::LogVerbatim("PMAlgTrackMaker") << std::endl << "----------- matched clusters: -----------";
	for (size_t i = 0; i < clusters->size(); ++i)
		if (has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    view: " << (*clusters)[i].View()
				<< "; size: " << (*clusters)[i].NHits();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "--------- not matched clusters: ---------";
	for (size_t i = 0; i < clusters->size(); ++i)
		if (!has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    view: " << (*clusters)[i].View()
				<< "; size: " << (*clusters)[i].NHits();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "-----------------------------------------";
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


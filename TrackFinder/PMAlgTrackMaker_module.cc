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
#include "RecoAlg/PMAlg/Utilities.h"

#include "TTree.h"

#include <memory>

namespace trkf {

typedef std::map< size_t, std::vector<double> > dedx_map;
typedef std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > view_hitmap;
typedef std::map< unsigned int, view_hitmap > tpc_view_hitmap;
typedef std::map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

struct TrkCandidate {
	pma::Track3D* Track;
	std::vector< size_t > Clusters;
	double Mse, Validation;
	bool Good;
};

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

  bool extendTrack(
	TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes);

  int matchCluster(
    const TrkCandidate& trk,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo);

  // display what was used and what is left, just for development and debugging
  void listUsedClusters(art::Handle< std::vector<recob::Cluster> > clusters) const;
  // ------------------------------------------------------

  // build tracks from clusters associated by any other module - not yet implemented
  int fromExistingAssocs(const art::Event& evt, std::vector< pma::Track3D* >& result);
  // ------------------------------------------------------


  // ************* some common functionality **************
  cryo_tpc_view_hitmap c_t_v_hits;
  bool sortHits(const art::Event& evt);

  std::vector< size_t > used_clusters, initial_clusters;
  std::map< unsigned int, std::vector<size_t> > tried_clusters;
  bool has(const std::vector<size_t>& v, size_t idx) const
  {
  	for (auto c : v) if (c == idx) return true;
  	return false;
  }

  std::vector< TrkCandidate > fCandidates;

  int maxCluster(
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);
  int maxCluster(
	size_t first_idx,
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
  bool fFlipToBeam;            // set the track direction to increasing Z values
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
	fFlipToBeam = pset.get< bool >("FlipToBeam");
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
	return recob::Track(xyz, dircos, dst_dQdx, std::vector< double >(2, util::kBogusD), fTrkIndex);
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
		return (particle->EndE() - particle->Mass() < 0.001);
		//return (particle->NumberDaughters() == 0);
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

bool PMAlgTrackMaker::extendTrack(TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes)
{
	double m_max = 2.0 * candidate.Mse; // max acceptable MSE value
	if (m_max < 0.05) m_max = 0.05;     // this is still good, low MSE value

	double v_min1 = 0.98 * candidate.Validation;
	double v_min2 = 0.9 * candidate.Validation;

	pma::Track3D* copy = fProjectionMatchingAlg.extendTrack(*(candidate.Track), hits, add_nodes);
	double m1 = copy->GetMse();
	double v1 = validate(*copy, testView);

	if (((m1 < candidate.Mse) && (v1 >= v_min2)) ||
	    ((m1 < 0.5) && (m1 <= m_max) && (v1 >= v_min1)))
	{
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track EXTENDED, MSE = " << m1 << ", v = " << v1;
		delete candidate.Track;    // delete previous version of the track
		candidate.Track = copy;    // replace with the new track
		copy->SortHits();          // sort hits in the new track

		candidate.Mse = m1;        // save info
		candidate.Validation = v1;

		return true;
	}
	else
	{
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track NOT extended, MSE = " << m1 << ", v = " << v1;
		delete copy;
		return false;
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

		double dQdxFlipThr = 0.0;
		if (fFlipToBeam) dQdxFlipThr = 0.25;

		fTrkIndex = 0;
		for (auto const& trk : result)
		{
			// flip the track by dQ/dx
			if (fFlipToBeam)
			{
				double z0 = trk->front()->Point3D().Z();
				double z1 = trk->back()->Point3D().Z();
				if (z0 < z1) trk->Flip(); // track points stored backwards, so start shold be higher z...
			}
			if (fAutoFlip_dQdx) fProjectionMatchingAlg.autoFlip(*trk, pma::Track3D::kBackward, dQdxFlipThr);

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
		initial_clusters.clear();
		tried_clusters.clear();
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
		listUsedClusters(cluListHandle);
	}
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "no clusters";
		return -1;
	}

	return result.size();
}

int PMAlgTrackMaker::matchCluster(const TrkCandidate& trk,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo)
{
	double f, fmax = 0.0;
	unsigned int n, max = 0;
	int idx = -1;
	for (size_t i = 0; i < clusters->size(); ++i)
	{
		unsigned int view = (*clusters)[i].View();
		unsigned int nhits = (*clusters)[i].NHits();

		if (has(used_clusters, i) ||                             // don't try already used clusters
			has(trk.Clusters, i) ||                              // don't try clusters from this candidate
		    (view == testView) ||                                // don't use clusters from validation view
		    ((preferedView != geo::kUnknown)&&(view != preferedView)) || // only prefered view if specified
		    (nhits < minSize))                                   // skip small clusters
		    continue;

		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);
		n = fProjectionMatchingAlg.testHits(*(trk.Track), v);
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
	initial_clusters.clear();

	size_t minSizeCompl = minBuildSize / 5; // smaller minimum required in complementary views

	int max_first_idx = 0;
	while (max_first_idx >= 0) // loop over clusters, any view, starting from the largest
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

			tried_clusters[first_view].push_back(max_first_idx);
			initial_clusters.push_back(max_first_idx);

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


			fCandidates.clear(); // possible solutions of the selected cluster and clusters in complementary views

			bool try_build = true;
			while (try_build) // loop over complementary views
			{
				TrkCandidate candidate;
				candidate.Clusters.push_back(max_first_idx);

				int idx, max_sec_a_idx, max_sec_b_idx;
				max_sec_a_idx = maxCluster(max_first_idx, clusters, fbp, tmin, tmax, minSizeCompl, sec_view_a, tpc, cryo);
				max_sec_b_idx = maxCluster(max_first_idx, clusters, fbp, tmin, tmax, minSizeCompl, sec_view_b, tpc, cryo);

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

					candidate.Clusters.push_back(idx);
					candidate.Track = fProjectionMatchingAlg.buildTrack(v_first, fbp.at(idx));

					double m0 = candidate.Track->GetMse();
					double v0 = validate(*(candidate.Track), testView);
					if ((m0 < 0.15) && (v0 > 0.7)) // good candidate, try to extend it
					{
						mf::LogVerbatim("PMAlgTrackMaker")
							<< "  track candidate, MSE = " << m0 << ", v = " << v0;

						candidate.Mse = m0;
						candidate.Validation = v0;
						candidate.Good = true;

						size_t minSize = 5;      // min size for clusters matching
						double fraction = 0.4;   // min fraction of close hits

						idx = 0;
						while (idx >= 0) // try to collect matching clusters, use **any** plane except validation
						{
							idx = matchCluster(candidate, clusters, fbp, minSize, fraction, geo::kUnknown, testView, tpc, cryo);
							if (idx >= 0)
							{
								// try building extended copy:
								//                src,        hits,    valid.plane, add nodes
								if (extendTrack(candidate, fbp.at(idx),  testView,    true))
								{
									candidate.Clusters.push_back(idx);
								}
								else idx = -1;
							}
						}

						mf::LogVerbatim("PMAlgTrackMaker") << "merge clusters from the validation plane";
						fraction = 0.7; // only well matching the existing track

						idx = 0;
						while ((idx >= 0) && (testView != geo::kUnknown))
						{	//                     match clusters from the plane used previously for the validation
							idx = matchCluster(candidate, clusters, fbp, minSize, fraction, testView, geo::kUnknown, tpc, cryo);
							if (idx >= 0)
							{
								// no validation, no new nodes:
								if (extendTrack(candidate, fbp.at(idx), geo::kUnknown, false))
								{
									candidate.Clusters.push_back(idx);
								}
								else idx = -1;
							}
						}

						candidate.Validation = validate(*(candidate.Track), testView);
					}
					else
					{
						mf::LogVerbatim("PMAlgTrackMaker") << "track REJECTED, MSE = " << m0 << "; v = " << v0;
						candidate.Good = false; // save also bad matches to avoid trying again the same pair of clusters
					}
					fCandidates.push_back(candidate);
				}
				else
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "no matching clusters";
				}
			} // end loop over complementary views

			if (fCandidates.size()) // save best candidate, release other tracks and clusters
			{                       // now "best" = min. MSE, but other condition might be better (validation, length, etc)

				std::cout << " *********** " << fCandidates.size() << " candidates" << std::endl;

				size_t best_trk = 0;
				double min_mse = 10.;
				for (size_t t = 0; t < fCandidates.size(); t++)
					if (fCandidates[t].Good && (fCandidates[t].Mse < min_mse))
					{
						min_mse = fCandidates[t].Mse;
						best_trk = t;
					}

				for (size_t t = 0; t < fCandidates.size(); t++)
				{
					if (fCandidates[t].Good)
					{
						std::cout << "             "
							<< fCandidates[t].Mse << " "
							<< fCandidates[t].Track->size() << " "
							<< fCandidates[t].Track->Length() << std::endl;
					}

					if (fCandidates[t].Good && (t == best_trk))
					{
						result.push_back(fCandidates[best_trk].Track);
						for (auto c : fCandidates[best_trk].Clusters)
							used_clusters.push_back(c);
					}
					else fCandidates[t].Track;
				}
			}
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
		    has(initial_clusters, i) ||
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

int PMAlgTrackMaker::maxCluster(size_t first_idx,
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
		    has(initial_clusters, i) ||
		    has(tried_clusters[view], i) ||
		    ((*clusters)[i].NHits() <  min_clu_size) ||
		    ((*clusters)[i].View() != view)) continue;

		bool pair_checked = false;
		for (auto const & c : fCandidates)
			if (has(c.Clusters, first_idx) && has(c.Clusters, i))
			{
				pair_checked = true; break;
			}
		if (pair_checked) continue;
		    
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


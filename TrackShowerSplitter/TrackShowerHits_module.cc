////////////////////////////////////////////////////////////////////////
// Class:       TrackShowerHits
// Module Type: producer
// File:        TrackShowerHits_module.cc
// Authors: dorota.stefan@cern.ch robert.sulej@cern.ch
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

#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoAlg/PMAlg/Utilities.h"
#include "Utilities/AssociationUtil.h"

#include "TrackShowerSplitter/Segmentation2D/TssHit2D.h"
#include "TrackShowerSplitter/Segmentation2D/SimpleClustering.h"

#include <memory>

namespace tss {

typedef std::map< unsigned int, std::vector< tss::Hit2D > > view_hitmap;
typedef std::map< unsigned int, view_hitmap > tpc_view_hitmap;
typedef std::map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

class TrackShowerHits : public art::EDProducer {
public:
	explicit TrackShowerHits(fhicl::ParameterSet const & p);

	TrackShowerHits(TrackShowerHits const &) = delete;
	TrackShowerHits(TrackShowerHits &&) = delete;
	TrackShowerHits & operator = (TrackShowerHits const &) = delete;
	TrackShowerHits & operator = (TrackShowerHits &&) = delete;

	void reconfigure(fhicl::ParameterSet const& p);

	void produce(art::Event & e) override;

private:

	cryo_tpc_view_hitmap c_t_v_hits;
	bool sortHits(const art::Event& evt);

	art::ServiceHandle<geo::Geometry> fGeom;

	tss::SimpleClustering simpleClustering;

	std::string fHitModuleLabel;

};
// ------------------------------------------------------

TrackShowerHits::TrackShowerHits(fhicl::ParameterSet const & p)
{
	this->reconfigure(p);

	produces< std::vector<recob::Hit> >();
	produces< std::vector<recob::Cluster> >();
	produces< art::Assns<recob::Cluster, recob::Hit> >();
}
// ------------------------------------------------------

void TrackShowerHits::reconfigure(fhicl::ParameterSet const& pset)
{
        fHitModuleLabel = pset.get< std::string >("HitModuleLabel");
}
// ------------------------------------------------------

bool TrackShowerHits::sortHits(const art::Event& evt)
{
	c_t_v_hits.clear();

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

			c_t_v_hits[cryo][tpc][view].emplace_back(tss::Hit2D(h));
		}
		return true;
	}
	else return false;
}
// ------------------------------------------------------

void TrackShowerHits::produce(art::Event & evt)
{
	std::unique_ptr< std::vector< recob::Hit > > track_hits(new std::vector< recob::Hit >);

	std::unique_ptr< std::vector< recob::Cluster > > clusters(new std::vector< recob::Cluster >);
	std::unique_ptr< art::Assns< recob::Cluster, recob::Hit > > clu2hit(new art::Assns< recob::Cluster, recob::Hit >);

	if (sortHits(evt))
	{
		unsigned int cidx = 0;
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			for (const auto & v : c_t_v_hits[tpc_iter->Cryostat][tpc_iter->TPC])
			{
				auto cls = simpleClustering.run(v.second);
				for (const auto & c : cls)
				{
					if (c.hits().size() < 2) continue;

					//const tss::Hit2D* ho = c.outermost();
					//std::cout << "outermost: " << ho->Point2D().X() << " " << ho->Point2D().Y() << std::endl;

					clusters->emplace_back(
						recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
							c.hits().size(), 0.0F, 0.0F, cidx,  (geo::View_t)c.hits().front()->View(), c.hits().front()->Hit2DPtr()->WireID().planeID()));

					std::vector< art::Ptr< recob::Hit > > hits2d;
					hits2d.reserve(c.hits().size());
					for (auto h2d : c.hits())
					{
						hits2d.push_back(h2d->Hit2DPtr());
					}
					if (hits2d.size())
					{
						util::CreateAssn(*this, evt, *clusters, hits2d, *clu2hit);
					}

					++cidx;
				}
			}
		}

	}
	else mf::LogWarning("TrackShowerHits") << "Hits not found in the event.";

	evt.put(std::move(track_hits));
	evt.put(std::move(clusters));
	evt.put(std::move(clu2hit));
}
// ------------------------------------------------------

DEFINE_ART_MODULE(TrackShowerHits)

}


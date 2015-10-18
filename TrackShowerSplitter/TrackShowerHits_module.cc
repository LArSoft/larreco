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
#include "Utilities/AssociationUtil.h"

#include <memory>

namespace tss {

typedef std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > view_hitmap;
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

	bool hitsTouching(const recob::Hit & h1, const recob::Hit & h2) const;
	bool hitsTouching(const std::vector< art::Ptr<recob::Hit> > & c1, const recob::Hit & h2) const;
	bool hitsTouching(const std::vector< art::Ptr<recob::Hit> > & c1, const std::vector< art::Ptr<recob::Hit> > & c2) const;

	void simpleClustering(const std::vector< art::Ptr<recob::Hit> > & inp);

	art::ServiceHandle<geo::Geometry> fGeom;

	std::string fHitModuleLabel;

};
// ------------------------------------------------------

TrackShowerHits::TrackShowerHits(fhicl::ParameterSet const & p)
{
	this->reconfigure(p);
	produces< std::vector<recob::Hit> >();
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

			c_t_v_hits[cryo][tpc][view].push_back(h);
		}
		return true;
	}
	else return false;
}
// ------------------------------------------------------

bool TrackShowerHits::hitsTouching(const recob::Hit & h1, const recob::Hit & h2) const
{
	if ((h1.WireID().Wire == h2.WireID().Wire) &&
	    (h1.PeakTime() == h2.PeakTime())) return false;

	bool touches = false;
	if ((h1.WireID().Wire >= h2.WireID().Wire - 1) &&
	    (h1.WireID().Wire <= h2.WireID().Wire + 1))
	{
		if (((h2.StartTick() <= h1.StartTick()) && (h1.StartTick() <= h2.EndTick() + 1)) ||
			((h2.StartTick() <= h1.EndTick() + 1) && (h1.EndTick() <= h2.EndTick())) ||
			((h2.StartTick() >= h1.StartTick()) && (h1.EndTick() >= h2.EndTick())))
		{
			touches = true;
		}
	}
	return touches;
}
// ------------------------------------------------------

bool TrackShowerHits::hitsTouching(const std::vector< art::Ptr<recob::Hit> > & c1, const recob::Hit & h2) const
{
	for (size_t i = 0; i < c1.size(); i++)
	{
		if (hitsTouching(*(c1[i]), h2)) return true;	
	}
	return false;
}
// ------------------------------------------------------

bool TrackShowerHits::hitsTouching(const std::vector< art::Ptr<recob::Hit> > & c1, const std::vector< art::Ptr<recob::Hit> > & c2) const
{
	for (unsigned int i = 0; i < c1.size(); i++)
	{
		if (hitsTouching(c2, *(c1[i]))) return true;
	}
	return false;
}
// ------------------------------------------------------

void TrackShowerHits::simpleClustering(const std::vector< art::Ptr<recob::Hit> > & inp)
{
	std::vector< std::vector< art::Ptr<recob::Hit> > > result;
	for (size_t h = 0; h < inp.size(); ++h)
	{
		bool found = false;
		for (size_t r = 0; r < result.size(); ++r)
			if (hitsTouching(result[r], *(inp[h])))
		{
			result[r].push_back(inp[h]); found = true; break;
		}
		if (!found)
		{
			result.push_back(std::vector< art::Ptr<recob::Hit> >());
			result.back().push_back(inp[h]);
		}
	}

	// merge...
}
// ------------------------------------------------------

void TrackShowerHits::produce(art::Event & evt)
{
	std::unique_ptr< std::vector< recob::Hit > > track_hits(new std::vector< recob::Hit >);

	if (sortHits(evt))
	{
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			for (const auto & v : c_t_v_hits[tpc_iter->Cryostat][tpc_iter->TPC])
			{
				simpleClustering(v.second);
			}
		}

	}
	else mf::LogWarning("TrackShowerHits") << "Hits not found in the event.";

	evt.put(std::move(track_hits));
}
// ------------------------------------------------------

DEFINE_ART_MODULE(TrackShowerHits)

}

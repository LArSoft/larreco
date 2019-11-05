/////////////////////////////////////////////////////////////////////////////////
// Class:       EmTrackClusterId2out
// Module Type: producer
// File:        EmTrackClusterId2out_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// Module applies CNN to 2D image made of deconvoluted wire waveforms in order
// to distinguish EM-like activity from track-like objects. New clusters of
// hits are produced to include also unclustered hits and tag everything in
// a common way.
// NOTE: This module uses 2-output CNN models, see EmTrackClusterId and
// EmTrackMichelClusterId for usage of 3 and 4-output models.
//
/////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"
#include "lardata/ArtDataHelper/MVAWriter.h"

#include <memory>

namespace nnet {

class EmTrackClusterId2out : public art::EDProducer {
public:

	// these types to be replaced with use of feature proposed in redmine #12602
	typedef std::unordered_map< unsigned int, std::vector< size_t > > view_keymap;
	typedef std::unordered_map< unsigned int, view_keymap > tpc_view_keymap;
	typedef std::unordered_map< unsigned int, tpc_view_keymap > cryo_tpc_view_keymap;

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<nnet::PointIdAlg::Config> PointIdAlg {
			Name("PointIdAlg")
		};
                fhicl::Atom<size_t> BatchSize {
                        Name("BatchSize"),
                        Comment("number of samples processed in one batch")
                };

		fhicl::Atom<art::InputTag> WireLabel {
			Name("WireLabel"),
			Comment("tag of deconvoluted ADC on wires (recob::Wire)")
		};

		fhicl::Atom<art::InputTag> HitModuleLabel {
			Name("HitModuleLabel"),
			Comment("tag of hits to be EM/track tagged")
		};

		fhicl::Atom<art::InputTag> ClusterModuleLabel {
			Name("ClusterModuleLabel"),
			Comment("tag of clusters to be used as a source of EM/track tagged new clusters (incl. single-hit clusters ) using accumulated results from hits")
		};

		fhicl::Atom<art::InputTag> TrackModuleLabel {
			Name("TrackModuleLabel"),
			Comment("tag of 3D tracks to be EM/track tagged using accumulated results from hits in the best 2D projection")
		};

		fhicl::Sequence<int> Views {
			Name("Views"),
			Comment("tag clusters in selected views only, or in all views if empty list")
		};
	};
	using Parameters = art::EDProducer::Table<Config>;
	explicit EmTrackClusterId2out(Parameters const & p);

	EmTrackClusterId2out(EmTrackClusterId2out const &) = delete;
	EmTrackClusterId2out(EmTrackClusterId2out &&) = delete;
	EmTrackClusterId2out & operator = (EmTrackClusterId2out const &) = delete;
	EmTrackClusterId2out & operator = (EmTrackClusterId2out &&) = delete;

private:
	void produce(art::Event & e) override;

	bool isViewSelected(int view) const;

        size_t fBatchSize;
	PointIdAlg fPointIdAlg;
	anab::MVAWriter<2> fMVAWriter; // <-------------- using 2-output CNN model

	art::InputTag fWireProducerLabel;
	art::InputTag fHitModuleLabel;
	art::InputTag fClusterModuleLabel;
	art::InputTag fTrackModuleLabel;
	bool fDoClusters, fDoTracks;

	std::vector< int > fViews;

        art::InputTag fNewClustersTag; // input tag for the clusters produced by this module
};
// ------------------------------------------------------

EmTrackClusterId2out::EmTrackClusterId2out(EmTrackClusterId2out::Parameters const& config) :
        EDProducer{config},
	fBatchSize(config().BatchSize()),
        fPointIdAlg(config().PointIdAlg()),
        fMVAWriter(producesCollector(), "emtrack"),
	fWireProducerLabel(config().WireLabel()),
	fHitModuleLabel(config().HitModuleLabel()),
	fClusterModuleLabel(config().ClusterModuleLabel()),
	fTrackModuleLabel(config().TrackModuleLabel()),
	fViews(config().Views()),

	fNewClustersTag(
	    config.get_PSet().get<std::string>("module_label"), "",
	    art::ServiceHandle<art::TriggerNamesService const>()->getProcessName())
{
    fMVAWriter.produces_using< recob::Hit >();

    if (!fClusterModuleLabel.label().empty())
    {
    	produces< std::vector<recob::Cluster> >();
	    produces< art::Assns<recob::Cluster, recob::Hit> >();

        fMVAWriter.produces_using< recob::Cluster >();
        fDoClusters = true;
    }
    else { fDoClusters = false; }

    if (!fTrackModuleLabel.label().empty())
    {
        fMVAWriter.produces_using< recob::Track >();
        fDoTracks = true;
    }
    else { fDoTracks = false; }
}
// ------------------------------------------------------

void EmTrackClusterId2out::produce(art::Event & evt)
{
    mf::LogVerbatim("EmTrackClusterId2out") << "next event: " << evt.run() << " / " << evt.id().event();

	auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);

	unsigned int cryo, tpc, view;

    // ******************* get and sort hits ********************
	auto hitListHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
	std::vector< art::Ptr<recob::Hit> > hitPtrList;
	art::fill_ptr_vector(hitPtrList, hitListHandle);

    EmTrackClusterId2out::cryo_tpc_view_keymap hitMap;
	for (auto const& h : hitPtrList)
	{
		view = h->WireID().Plane;
		if (!isViewSelected(view)) continue;

		cryo = h->WireID().Cryostat;
		tpc = h->WireID().TPC;

		hitMap[cryo][tpc][view].push_back(h.key());
	}

    // ********************* classify hits **********************
    auto hitID = fMVAWriter.initOutputs<recob::Hit>(fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());

    std::vector< char > hitInFA(hitPtrList.size(), 0); // tag hits in fid. area as 1, use 0 for hits close to the projectrion edges
    for (auto const & pcryo : hitMap)
    {
        cryo = pcryo.first;
        for (auto const & ptpc : pcryo.second)
        {
            tpc = ptpc.first;
            for (auto const & pview : ptpc.second)
            {
                view = pview.first;
                if (!isViewSelected(view)) continue; // should not happen, hits were selected

                fPointIdAlg.setWireDriftData(*wireHandle, view, tpc, cryo);

                // (1) do all hits in this plane ------------------------------------------------
                for (size_t idx = 0; idx < pview.second.size(); idx += fBatchSize)
                {
                    std::vector< std::pair<unsigned int, float> > points;
                    std::vector< size_t > keys;
                    for (size_t k = 0; k < fBatchSize; ++k)
                    {
                        if (idx + k >= pview.second.size()) { break; } // careful about the tail

                        size_t h = pview.second[idx+k]; // h is the Ptr< recob::Hit >::key()
                        const recob::Hit & hit = *(hitPtrList[h]);
                        points.emplace_back(hit.WireID().Wire, hit.PeakTime());
                        keys.push_back(h);
                    }

                    auto batch_out = fPointIdAlg.predictIdVectors(points);
                    if (points.size() != batch_out.size())
                    {
                        throw cet::exception("EmTrackClusterId") << "hits processing failed" << std::endl;
                    }

                    for (size_t k = 0; k < points.size(); ++k)
                    {
                        size_t h = keys[k];
                        fMVAWriter.setOutput(hitID, h, batch_out[k]);
                        if (fPointIdAlg.isInsideFiducialRegion(points[k].first, points[k].second))
                        { hitInFA[h] = 1; }
                    }
                } // hits done ------------------------------------------------------------------
            }
        }
    }

    // (2) do clusters when hits are ready in all planes ----------------------------------------
    if (fDoClusters)
    {
        // **************** prepare for new clusters ****************
	    auto clusters = std::make_unique< std::vector< recob::Cluster > >();
	    auto clu2hit = std::make_unique< art::Assns< recob::Cluster, recob::Hit > >();

        // ************** get and sort input clusters ***************
    	auto cluListHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);
	    std::vector< art::Ptr<recob::Cluster> > cluPtrList;
	    art::fill_ptr_vector(cluPtrList, cluListHandle);

        EmTrackClusterId2out::cryo_tpc_view_keymap cluMap;
	    for (auto const& c : cluPtrList)
	    {
	    	view = c->Plane().Plane;
	    	if (!isViewSelected(view)) continue;

	    	cryo = c->Plane().Cryostat;
	    	tpc = c->Plane().TPC;

	    	cluMap[cryo][tpc][view].push_back(c.key());
	    }

        auto cluID = fMVAWriter.initOutputs<recob::Cluster>(fNewClustersTag, fPointIdAlg.outputLabels());

        unsigned int cidx = 0; // new clusters index
        art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fClusterModuleLabel);
        std::vector< bool > hitUsed(hitPtrList.size(), false); // tag hits used in clusters
        for (auto const & pcryo : cluMap)
        {
            cryo = pcryo.first;
            for (auto const & ptpc : pcryo.second)
            {
                tpc = ptpc.first;
                for (auto const & pview : ptpc.second)
                {
                    view = pview.first;
                    if (!isViewSelected(view)) continue; // should not happen, clusters were pre-selected

                    for (size_t c : pview.second) // c is the Ptr< recob::Cluster >::key()
                    {
		                auto v = hitsFromClusters.at(c);
		                if (v.empty()) continue;

                        for (auto const & hit : v)
                        {
                            if (hitUsed[hit.key()]) { mf::LogWarning("EmTrackClusterId2out") << "hit already used in another cluster"; }
                            hitUsed[hit.key()] = true;
                        }

                        auto vout = fMVAWriter.getOutput<recob::Hit>(v,
                            [&](art::Ptr<recob::Hit> const & ptr) { return (float)hitInFA[ptr.key()]; });

	    	            mf::LogVerbatim("EmTrackClusterId2out") << "cluster in tpc:" << tpc << " view:" << view
                            << " size:" << v.size() << " p:" << vout[0];

                		clusters->emplace_back(
	    		            recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
	    			            v.size(), 0.0F, 0.0F, cidx, (geo::View_t)view, v.front()->WireID().planeID()));
	    	            util::CreateAssn(*this, evt, *clusters, v, *clu2hit);
	    	            cidx++;

	    	            fMVAWriter.addOutput(cluID, vout); // add copy of the input cluster
                    }

                    // (2b) make single-hit clusters --------------------------------------------
                    for (size_t h : hitMap[cryo][tpc][view]) // h is the Ptr< recob::Hit >::key()
                    {
                        if (hitUsed[h]) continue;

                        auto vout = fMVAWriter.getOutput<recob::Hit>(h);

		                mf::LogVerbatim("EmTrackClusterId2out") << "single hit in tpc:" << tpc << " view:" << view
			                << " wire:" << hitPtrList[h]->WireID().Wire << " drift:" << hitPtrList[h]->PeakTime() << " p:" << vout[0];

		                art::PtrVector< recob::Hit > cluster_hits;
		                cluster_hits.push_back(hitPtrList[h]);
		                clusters->emplace_back(
			                recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
			                    1, 0.0F, 0.0F, cidx, (geo::View_t)view, hitPtrList[h]->WireID().planeID()));
		                util::CreateAssn(*this, evt, *clusters, cluster_hits, *clu2hit);
		                cidx++;

		                fMVAWriter.addOutput(cluID, vout); // add single-hit cluster tagging unclutered hit
                    }
                    mf::LogVerbatim("EmTrackClusterId2out") << "...produced " << cidx - pview.second.size() << " single-hit clusters.";
                }
            }
        }

        evt.put(std::move(clusters));
	    evt.put(std::move(clu2hit));
	} // all clusters done ----------------------------------------------------------------------

    // (3) do tracks when all hits in all cryo/tpc/plane are done -------------------------------
    if (fDoTracks)
    {
        auto trkListHandle = evt.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);
        art::FindManyP< recob::Hit > hitsFromTracks(trkListHandle, evt, fTrackModuleLabel);
        std::vector< std::vector< art::Ptr<recob::Hit> > > trkHitPtrList(trkListHandle->size());
        for (size_t t = 0; t < trkListHandle->size(); ++t)
        {
            auto v = hitsFromTracks.at(t);
            size_t nh[3] = { 0, 0, 0 };
            for (auto const & hptr : v) { ++nh[hptr->View()]; }
            size_t best_view = 2; // collection
            if ((nh[0] >= nh[1]) && (nh[0] > 2 * nh[2])) best_view = 0; // ind1
            if ((nh[1] >= nh[0]) && (nh[1] > 2 * nh[2])) best_view = 1; // ind2

            size_t k = 0;
            while (!isViewSelected(best_view))
            {
                best_view = (best_view + 1) % 3;
                if (++k > 3) { throw cet::exception("EmTrackClusterId2out") << "No views selected at all?" << std::endl; }
            }

            for (auto const & hptr : v)
            {
                if (hptr->View() == best_view) trkHitPtrList[t].emplace_back(hptr);
            }
        }

        auto trkID = fMVAWriter.initOutputs<recob::Track>(fTrackModuleLabel, trkHitPtrList.size(), fPointIdAlg.outputLabels());
        for (size_t t = 0; t < trkHitPtrList.size(); ++t) // t is the Ptr< recob::Track >::key()
        {
            auto vout = fMVAWriter.getOutput<recob::Hit>(trkHitPtrList[t],
                [&](art::Ptr<recob::Hit> const & ptr) { return (float)hitInFA[ptr.key()]; });
            fMVAWriter.setOutput(trkID, t, vout);
        }
    }
    // tracks done ------------------------------------------------------------------------------

	fMVAWriter.saveOutputs(evt);
}
// ------------------------------------------------------

bool EmTrackClusterId2out::isViewSelected(int view) const
{
	if (fViews.empty()) return true;
	else
	{
		for (auto k : fViews) if (k == view) { return true; }
		return false;
	}
}
// ------------------------------------------------------

DEFINE_ART_MODULE(EmTrackClusterId2out)

}

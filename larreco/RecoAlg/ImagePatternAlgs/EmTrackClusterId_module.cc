/////////////////////////////////////////////////////////////////////////////////
// Class:       EmTrackClusterId
// Module Type: producer
// File:        EmTrackClusterId_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// Module applies neural net to 2D image made of deconvoluted wire waveforms in
// order to distinguish EM-like activity from track-like objects. Clusters of
// hits that were recognized as EM-like event parts are produced. Module uses
// clusters made with any algorithm as input; optionally (recommended) also
// single, unclustered hits are tested and if they look like EM parts then
// single-hit clusters are produced and added to the output collection.
//
/////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"
#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/MVAWriter.h"

#include <memory>

namespace nnet {

class EmTrackClusterId : public art::EDProducer {
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
			Comment("tag of clusters, to be EM/track tagged using accumulated results from hits")
		};

		fhicl::Atom<double> Threshold {
			Name("Threshold"),
			Comment("tag cluster as EM-like if net output accumulated over hits > threshold")
		};

		fhicl::Sequence<int> Views {
			Name("Views"),
			Comment("tag clusters in selected views only, or in all views if empty list")
		};
	};
	using Parameters = art::EDProducer::Table<Config>;
	explicit EmTrackClusterId(Parameters const & p);

	EmTrackClusterId(EmTrackClusterId const &) = delete;
	EmTrackClusterId(EmTrackClusterId &&) = delete;
	EmTrackClusterId & operator = (EmTrackClusterId const &) = delete;
	EmTrackClusterId & operator = (EmTrackClusterId &&) = delete;

	void produce(art::Event & e) override;

private:
	bool isViewSelected(int view) const;

	PointIdAlg fPointIdAlg;
	anab::MVAWriter<3> fMVAWriter;

	art::InputTag fWireProducerLabel;
	art::InputTag fHitModuleLabel;
	art::InputTag fClusterModuleLabel;

	double fThreshold;

	std::vector< int > fViews;
};
// ------------------------------------------------------

EmTrackClusterId::EmTrackClusterId(EmTrackClusterId::Parameters const& config) :
	fPointIdAlg(config().PointIdAlg()), fMVAWriter(this, "emtrack"),
	fWireProducerLabel(config().WireLabel()),
	fHitModuleLabel(config().HitModuleLabel()),
	fClusterModuleLabel(config().ClusterModuleLabel()),
	fThreshold(config().Threshold()),
	fViews(config().Views())
{
	produces< std::vector<recob::Cluster> >();
	produces< art::Assns<recob::Cluster, recob::Hit> >();

    std::cout << fMVAWriter << std::endl;

	fMVAWriter.produces< recob::Hit >();
	fMVAWriter.produces< recob::Cluster >();

	std::cout << fMVAWriter << std::endl;
}
// ------------------------------------------------------

void EmTrackClusterId::produce(art::Event & evt)
{
    mf::LogVerbatim("EmTrackClusterId") << "next event: " << evt.run() << " / " << evt.id().event();

	auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);

	auto clusters = std::make_unique< std::vector< recob::Cluster > >();
	auto clu2hit = std::make_unique< art::Assns< recob::Cluster, recob::Hit > >();

	unsigned int cidx = 0;
	unsigned int cryo, tpc, view;

    // ******************* get and sort hits ********************
	auto hitListHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
	std::vector< art::Ptr<recob::Hit> > hitPtrList;
	art::fill_ptr_vector(hitPtrList, hitListHandle);

    EmTrackClusterId::cryo_tpc_view_keymap hitMap;
	for (auto const& h : hitPtrList)
	{
		view = h->WireID().Plane;
		if (!isViewSelected(view)) continue;

		cryo = h->WireID().Cryostat;
		tpc = h->WireID().TPC;

		hitMap[cryo][tpc][view].push_back(h.key());
	}
	
    // ***************** get and sort clusters ******************
	auto cluListHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);
	std::vector< art::Ptr<recob::Cluster> > cluPtrList;
	art::fill_ptr_vector(cluPtrList, cluListHandle);

    EmTrackClusterId::cryo_tpc_view_keymap cluMap;
	for (auto const& c : cluPtrList)
	{
		view = c->Plane().Plane;
		if (!isViewSelected(view)) continue;

		cryo = c->Plane().Cryostat;
		tpc = c->Plane().TPC;

		cluMap[cryo][tpc][view].push_back(c.key());
	}

    // *************** classify hits and clusters ****************
    auto hitID = fMVAWriter.initOutputs<recob::Hit>(fHitModuleLabel, hitPtrList.size());
    auto cluID = fMVAWriter.initOutputs<recob::Cluster>(fHitModuleLabel); /// USE THIS MODULE LABEL FOR CLUSTERS PRODUCED HERE

    art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fClusterModuleLabel);
    std::vector< bool > hitUsed(hitPtrList.size(), false); // tag hits used in clusters
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

                // (1) do all hits in this plane --------------------------------------------
                for (size_t h : pview.second) // h is the Ptr::key()
                {
                    const recob::Hit & hit = *(hitPtrList[h]);
                    fMVAWriter.setOutput(hitID, h,
                        fPointIdAlg.predictIdVector(hit.WireID().Wire, hit.PeakTime()));
                } // hits done --------------------------------------------------------------

                // (2a) do clusters when hits are ready in this plane -----------------------
                for (size_t c : cluMap[cryo][tpc][view]) // c is the Ptr::key()
                {
		            auto v = hitsFromClusters.at(c);
		            if (v.empty()) continue;

                    auto vout = fMVAWriter.getOutput<recob::Hit>(v,
                        [&](auto const & hit)
                        {
                            if (fPointIdAlg.isInsideFiducialRegion(hit.WireID().Wire, hit.PeakTime())) return 1.0F;
                            else return 0.0F;
                        });

                    for (auto const & hit : v)
                    {
                        if (hitUsed[hit.key()]) { mf::LogWarning("EmTrackClusterId") << "hit already used in another cluster"; }
                        hitUsed[hit.key()] = true;
                    }

		            float pvalue = vout[0] / (vout[0] + vout[1]);
		            mf::LogVerbatim("EmTrackClusterId") << "cluster in tpc:" << tpc << " view:" << view
                        << " size:" << v.size() << " p:" << pvalue;

		            if (pvalue > fThreshold) continue; // identified as track-like, for the moment we save only em-like parts

            		clusters->emplace_back(
			            recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
				            v.size(), 0.0F, 0.0F, cidx, (geo::View_t)view, v.front()->WireID().planeID()));
		            util::CreateAssn(*this, evt, *clusters, v, *clu2hit);
		            cidx++;

		            fMVAWriter.addOutput(cluID, vout);
                }
                // (2b) make single-hit clusters --------------------------------------------
                for (size_t h : pview.second) // h is the Ptr::key()
                {
                    if (hitUsed[h]) continue;

                    auto vout = fMVAWriter.getOutput<recob::Hit>(h);
		            float pvalue = vout[0] / (vout[0] + vout[1]);

		            mf::LogVerbatim("EmTrackClusterId") << "single hit in tpc:" << tpc << " view:" << view
			            << " wire:" << hitPtrList[h]->WireID().Wire << " drift:" << hitPtrList[h]->PeakTime() << " p:" << pvalue;

		            if (pvalue > fThreshold) continue; // identified as track-like, for the moment we save only em-like parts

		            art::PtrVector< recob::Hit > cluster_hits;
		            cluster_hits.push_back(hitPtrList[h]);
		            clusters->emplace_back(
			            recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
			                1, 0.0F, 0.0F, cidx, (geo::View_t)view, hitPtrList[h]->WireID().planeID()));
		            util::CreateAssn(*this, evt, *clusters, cluster_hits, *clu2hit);
		            cidx++;

		            fMVAWriter.addOutput(cluID, vout);
                } // all clusters done ------------------------------------------------------

                // (3) do tracks which are best seen in this plane --------------------------
                // ....TO DO....
                // tracks done --------------------------------------------------------------
            }
        }
    }

	evt.put(std::move(clusters));
	evt.put(std::move(clu2hit));

	fMVAWriter.saveOutputs(evt);
	std::cout << "MVA saved" << std::endl;
}
// ------------------------------------------------------

bool EmTrackClusterId::isViewSelected(int view) const
{
	if (fViews.empty()) return true;
	else
	{
		bool selected = false;
		for (auto k : fViews) if (k == view) { selected = true; break; }
		return selected;
	}
}
// ------------------------------------------------------

DEFINE_ART_MODULE(EmTrackClusterId)

}


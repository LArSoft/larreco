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

#include <memory>

namespace nnet {

class EmTrackClusterId : public art::EDProducer {
public:

	// these types to be replaced with use of feature proposed in redmine #12602
	typedef std::unordered_map< unsigned int, std::vector< size_t > > view_hitmap;
	typedef std::unordered_map< unsigned int, view_hitmap > tpc_view_hitmap;
	typedef std::unordered_map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

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

    std::vector<float> pValuesMultiply(
        std::vector< art::Ptr<recob::Hit> > const & hits,
        const std::vector< std::vector<float> > & pValues,
        size_t nClasses) const;

	PointIdAlg fPointIdAlg;

	art::InputTag fWireProducerLabel;
	art::InputTag fHitModuleLabel;
	art::InputTag fClusterModuleLabel;

	double fThreshold;

	std::vector< int > fViews;
};
// ------------------------------------------------------

EmTrackClusterId::EmTrackClusterId(EmTrackClusterId::Parameters const& config) :
	fPointIdAlg(config().PointIdAlg()),
	fWireProducerLabel(config().WireLabel()),
	fHitModuleLabel(config().HitModuleLabel()),
	fClusterModuleLabel(config().ClusterModuleLabel()),
	fThreshold(config().Threshold()),
	fViews(config().Views())
{
	produces< std::vector<recob::Cluster> >();
	produces< art::Assns<recob::Cluster, recob::Hit> >();
}
// ------------------------------------------------------

std::vector<float> EmTrackClusterId::pValuesMultiply(
    std::vector< art::Ptr<recob::Hit> > const & hits,
    const std::vector< std::vector<float> > & pValues,
    size_t nClasses) const
{
    if (pValues.empty())
    {
        mf::LogError("EmTrackClusterId") << "no pValues!";
        return std::vector<float>();
    }

	std::vector<float> result(nClasses, 0);
	if (result.empty()) return result;

	double pmin = 1.0e-6, pmax = 1.0 - pmin;
	double log_pmin = log(pmin), log_pmax = log(pmax);
	double totarea = 0.0;
	size_t nhits = 0;

	for (auto const & h : hits)
	{
		unsigned int wire = h->WireID().Wire;
		float drift = h->PeakTime();

		if (!fPointIdAlg.isInsideFiducialRegion(wire, drift)) continue;

		double area = 1.0; // h->SummedADC();

		auto vout = pValues[h.key()];
		for (size_t i = 0; i < vout.size(); ++i)
		{
			if (vout[i] < pmin) vout[i] = log_pmin;
			else if (vout[i] > pmax) vout[i] = log_pmax;
			else vout[i] = log(vout[i]);

			result[i] += area * vout[i];
		}
		totarea += area;
		nhits++;
	}

	if (nhits)
	{
		double totp = 0.0;
		for (size_t i = 0; i < result.size(); ++i)
		{
			result[i] = exp(result[i] / totarea);
			totp += result[i];
		}
		for (size_t i = 0; i < result.size(); ++i)
		{
			result[i] /= totp;
		}
	}
	else std::fill(result.begin(), result.end(), 1.0 / result.size());

	return result;
}

void EmTrackClusterId::produce(art::Event & evt)
{
    mf::LogVerbatim("EmTrackClusterId") << "next event: " << evt.run() << " / " << evt.id().event();

	auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);
	auto cluListHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);

	auto clusters = std::make_unique< std::vector< recob::Cluster > >();
	auto clu2hit = std::make_unique< art::Assns< recob::Cluster, recob::Hit > >();

	unsigned int cidx = 0;
	const unsigned int emTag = 0x10000;

    // ******************* get and sort hits ********************
	auto hitListHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);

	std::vector< art::Ptr<recob::Hit> > allhitlist;
	art::fill_ptr_vector(allhitlist, hitListHandle);

    EmTrackClusterId::cryo_tpc_view_hitmap hitMap;
	unsigned int cryo, tpc, view;
	for (auto const& h : allhitlist)
	{
		view = h->WireID().Plane;
		if (!isViewSelected(view)) continue;

		cryo = h->WireID().Cryostat;
		tpc = h->WireID().TPC;

		hitMap[cryo][tpc][view].push_back(h.key());
	}

    // ******************** classify all hits *******************
    std::vector< std::vector<float> > pValues(allhitlist.size());
    std::vector< bool > hUsed(allhitlist.size(), false);
    for (auto & cryo : hitMap)
        for (auto & tpc : cryo.second)
            for (auto & view : tpc.second)
            {
                fPointIdAlg.setWireDriftData(*wireHandle, view.first, tpc.first, cryo.first);
                for (auto & h : view.second)
                {
                    const recob::Hit & hit = *(allhitlist[h]);
				    pValues[h] = fPointIdAlg.predictIdVector(hit.WireID().Wire, hit.PeakTime());
                }
            }

	// ******************* classify clusters ********************
	art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fClusterModuleLabel);
	for (size_t i = 0; i < hitsFromClusters.size(); ++i)
	{
		auto v = hitsFromClusters.at(i);
		if (v.empty()) continue;

        int tpc = v.front()->WireID().TPC;
		int view = v.front()->View();
		if (!isViewSelected(view)) continue;

		auto vout = pValuesMultiply(v, pValues, fPointIdAlg.NClasses());
		float pvalue = vout[0] / (vout[0] + vout[1]);

        for (auto const & h : v)
        {
            if (hUsed[h.key()])
            {
                mf::LogWarning("EmTrackClusterId") << "hit already used in another cluster";
            }
            hUsed[h.key()] = true;
        }

		mf::LogVerbatim("EmTrackClusterId") << "cluster in tpc:" << tpc << " view:" << view
			<< " size:" << v.size() << " p:" << pvalue;

		if (pvalue > fThreshold) continue; // identified as track-like, for the moment we save only em-like parts

		clusters->emplace_back(
			recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
				v.size(), 0.0F, 0.0F, cidx + emTag, (geo::View_t)v.front()->View(), v.front()->WireID().planeID()));
		util::CreateAssn(*this, evt, *clusters, v, *clu2hit);
		cidx++;
	}
	// **********************************************************

	// *************** add single hits as clusters **************
	for (auto const & hi : allhitlist)
	{
	    const auto & vout = pValues[hi.key()];

	    if (hUsed[hi.key()] || vout.empty()) continue;

    	int view = hi->View();
		int tpc = hi->WireID().TPC;

		float pvalue = vout[0] / (vout[0] + vout[1]);

		mf::LogVerbatim("EmTrackClusterId") << "hit in tpc:" << tpc << " view:" << view
			<< " wire:" << hi->WireID().Wire << " drift:" << hi->PeakTime() << " p:" << pvalue;

		if (pvalue > fThreshold) continue; // identified as track-like, for the moment we save only em-like parts

		art::PtrVector< recob::Hit > cluster_hits;
		cluster_hits.push_back(hi);
		clusters->emplace_back(
			recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
			1, 0.0F, 0.0F, cidx + emTag, (geo::View_t)hi->View(), hi->WireID().planeID()));
		util::CreateAssn(*this, evt, *clusters, cluster_hits, *clu2hit);
		cidx++;
	}
	// **********************************************************

	evt.put(std::move(clusters));
	evt.put(std::move(clu2hit));
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


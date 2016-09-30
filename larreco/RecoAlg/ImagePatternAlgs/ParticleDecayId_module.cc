/////////////////////////////////////////////////////////////////////////////////
// Class:       ParticleDecayId
// Module Type: producer
// File:        ParticleDecayId_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// THIS IS STIIL DEVELOPMENT NOW - PLEASE WAIT A FEW DAYS FOR THE CODE IS MAKING
// ACTUAL VERTICES
//
/////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"

#include <memory>
#include <TVector3.h>

namespace nnet {

class ParticleDecayId : public art::EDProducer {
public:

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

		fhicl::Atom<art::InputTag> TrackModuleLabel {
			Name("TrackModuleLabel"),
			Comment("tag of tracks where decays points are to be found")
		};

		fhicl::Atom<double> Threshold {
			Name("Threshold"),
			Comment("tag decay point if it is detected in at least two planes with net output > threshold")
		};

		fhicl::Atom<int> SkipView {
			Name("SkipView"),
			Comment("use all views to find decays if -1, or skip the view with provided index and use only the two other views")
		};
	};
	using Parameters = art::EDProducer::Table<Config>;
	explicit ParticleDecayId(Parameters const & p);

	ParticleDecayId(ParticleDecayId const &) = delete;
	ParticleDecayId(ParticleDecayId &&) = delete;
	ParticleDecayId & operator = (ParticleDecayId const &) = delete;
	ParticleDecayId & operator = (ParticleDecayId &&) = delete;

	void produce(art::Event & e) override;

private:
    bool DetectDecay(
        const std::vector<recob::Wire> & wires,
        const std::vector< art::Ptr<recob::Hit> > & hits,
        std::map< size_t, TVector3 > & spoints,
        double* xyz);

	PointIdAlg fPointIdAlg;

	art::InputTag fWireProducerLabel;
	art::InputTag fTrackModuleLabel;

	double fThreshold;

	int fSkipView;
};
// ------------------------------------------------------

ParticleDecayId::ParticleDecayId(ParticleDecayId::Parameters const& config) :
	fPointIdAlg(config().PointIdAlg()),
	fWireProducerLabel(config().WireLabel()),
	fTrackModuleLabel(config().TrackModuleLabel()),
	fThreshold(config().Threshold()),
	fSkipView(config().SkipView())
{
	produces< std::vector<recob::Vertex> >();
	produces< art::Assns<recob::Vertex, recob::Track> >();
}
// ------------------------------------------------------

void ParticleDecayId::produce(art::Event & evt)
{
    std::cout << "event " << evt.id().event() << std::endl << std::endl;

    auto vtxs = std::make_unique< std::vector< recob::Vertex > >();
    auto vtx2trk = std::make_unique< art::Assns< recob::Vertex, recob::Track > >();

    auto vid = getProductID< std::vector<recob::Vertex> >(evt);

    auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);
    auto trkListHandle = evt.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);
    auto spListHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fTrackModuleLabel);

    art::FindManyP< recob::Hit > hitsFromTracks(trkListHandle, evt, fTrackModuleLabel);
    art::FindManyP< recob::SpacePoint > spFromTracks(trkListHandle, evt, fTrackModuleLabel);
    art::FindManyP< recob::Hit > hitsFromSPoints(spListHandle, evt, fTrackModuleLabel);

	for (size_t i = 0; i < hitsFromTracks.size(); ++i)
	{
		auto hits = hitsFromTracks.at(i);
		auto spoints = spFromTracks.at(i);
		if (hits.empty()) continue;

        std::map< size_t, TVector3 > trkSpacePoints;
        for (const auto & p : spoints)
        {
            auto sp_hits = hitsFromSPoints.at(p.key());
            for (const auto & h : sp_hits)
            {
                trkSpacePoints[h.key()] = TVector3(p->XYZ()[0], p->XYZ()[1], p->XYZ()[2]);
            }
        }

        double xyz[3];
        if (DetectDecay(*wireHandle, hits, trkSpacePoints, xyz))
        {
            size_t vidx = vtxs->size();
            vtxs->push_back(recob::Vertex(xyz, vidx));

            art::Ptr<recob::Track> tptr(trkListHandle, i);
            art::Ptr<recob::Vertex> vptr(vid, vidx, evt.productGetter(vid));
            vtx2trk->addSingle(vptr, tptr);
        }
    }

    evt.put(std::move(vtxs));
    evt.put(std::move(vtx2trk));
}
// ------------------------------------------------------

bool ParticleDecayId::DetectDecay(
    const std::vector<recob::Wire> & wires,
    const std::vector< art::Ptr<recob::Hit> > & hits,
    std::map< size_t, TVector3 > & spoints,
    double* xyz)
{
    const size_t nviews = 3;

    std::vector< art::Ptr<recob::Hit> > wire_drift[nviews];
    for (size_t i = 0; i < hits.size(); ++i) // split hits between views
    {
        wire_drift[hits[i]->View()].push_back(hits[i]);
    }

    std::vector< float > outputs[nviews];
    for (size_t v = 0; v < nviews; ++v)      // calculate nn outputs for each view (hopefully not changing cryo/tpc many times)
    {
        outputs[v].resize(wire_drift[v].size(), 0);
        for (size_t i = 0; i < wire_drift[v].size(); ++i)
        {
    		int tpc = wire_drift[v][i]->WireID().TPC;
	    	int cryo = wire_drift[v][i]->WireID().Cryostat;

	    	fPointIdAlg.setWireDriftData(wires, v, tpc, cryo);

	        outputs[v][i] = fPointIdAlg.predictIdVector(wire_drift[v][i]->WireID().Wire, wire_drift[v][i]->PeakTime())[0]; // p(decay)
	    }
	    std::cout << std::endl;
    }

    std::vector< std::pair<size_t, float> > candidates[nviews];
    for (size_t v = 0; v < nviews; ++v)
    {
        size_t idx = 0;
        while (idx < outputs[v].size())
        {
            if (outputs[v][idx] > fThreshold)
            {
                size_t ci = idx;
                float max = outputs[v][idx];
                ++idx;

                while ((idx < outputs[v].size()) && (outputs[v][idx] > fThreshold))
                {
                    if (outputs[v][idx] > max) { max = outputs[v][idx]; ci = idx; }
                    ++idx;
                }
                candidates[v].emplace_back(ci, max);
            }
            else ++idx;
        }
    }

    for (size_t v = 0; v < nviews; ++v)
    {
        for (size_t i = 0; i < candidates[v].size(); ++i)
        {
            TVector3 sp( spoints[wire_drift[v][candidates[v][i].first].key()] );

            std::cout << "view:" << v << " idx:" << candidates[v][i].first
                << " w:" << wire_drift[v][candidates[v][i].first]->WireID().Wire
                << " d:" << wire_drift[v][candidates[v][i].first]->PeakTime()
                << " x:" << sp.X() << " y:" << sp.Y() << " z:" << sp.Z()
                << " p:" << candidates[v][i].second << std::endl;
        }
    }

    std::cout << std::endl;

    return false;
}

DEFINE_ART_MODULE(ParticleDecayId)

}


/////////////////////////////////////////////////////////////////////////////////
// Class:       ParticleDecayId
// Module Type: producer
// File:        ParticleDecayId_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// THIS IS STIIL DEVELOPMENT NOW - CODE MAKING VERTICES MAY STILL CHANGE STRATEGY
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

#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"

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

		fhicl::Atom<double> RoiThreshold {
			Name("RoiThreshold"),
			Comment("search for decay points where the net output > ROI threshold")
		};

		fhicl::Atom<double> PointThreshold {
			Name("PointThreshold"),
			Comment("tag decay point if it is detected in at least two planes with net outputs product > POINT threshold")
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
        std::vector< std::pair<TVector3, double> > & result);

	PointIdAlg fPointIdAlg;

	art::InputTag fWireProducerLabel;
	art::InputTag fTrackModuleLabel;

	double fRoiThreshold, fPointThreshold;

	int fSkipView;
};
// ------------------------------------------------------

ParticleDecayId::ParticleDecayId(ParticleDecayId::Parameters const& config) :
        EDProducer{config},
	fPointIdAlg(config().PointIdAlg()),
	fWireProducerLabel(config().WireLabel()),
	fTrackModuleLabel(config().TrackModuleLabel()),
	fRoiThreshold(config().RoiThreshold()),
	fPointThreshold(config().PointThreshold()),
	fSkipView(config().SkipView())
{
	produces< std::vector<recob::Vertex> >();
	produces< art::Assns<recob::Vertex, recob::Track> >();
}
// ------------------------------------------------------

void ParticleDecayId::produce(art::Event & evt)
{
    std::cout << std::endl << "event " << evt.id().event() << std::endl;

    auto vtxs = std::make_unique< std::vector< recob::Vertex > >();
    auto vtx2trk = std::make_unique< art::Assns< recob::Vertex, recob::Track > >();

    //auto vid = getProductID< std::vector<recob::Vertex> >(evt);

    auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);
    auto trkListHandle = evt.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);
    auto spListHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fTrackModuleLabel);

    art::FindManyP< recob::Hit > hitsFromTracks(trkListHandle, evt, fTrackModuleLabel);
    art::FindManyP< recob::SpacePoint > spFromTracks(trkListHandle, evt, fTrackModuleLabel);
    art::FindManyP< recob::Hit > hitsFromSPoints(spListHandle, evt, fTrackModuleLabel);

    std::vector< std::pair<TVector3, double> > decays;
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

        DetectDecay(*wireHandle, hits, trkSpacePoints, decays);
    }

    double xyz[3];
    for (const auto & p3d : decays)
    {

        xyz[0] = p3d.first.X(); xyz[1] = p3d.first.Y(); xyz[2] = p3d.first.Z();
        std::cout << "   detected: [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "] p:" << p3d.second << std::endl;

        size_t vidx = vtxs->size();
        vtxs->push_back(recob::Vertex(xyz, vidx));

        // to do: assn to eg. appropriate track
        // selected among set of connected tracks
        //art::Ptr<recob::Track> tptr(trkListHandle, i);
        //art::Ptr<recob::Vertex> vptr(vid, vidx, evt.productGetter(vid));
        //vtx2trk->addSingle(vptr, tptr);
    }

    evt.put(std::move(vtxs));
    evt.put(std::move(vtx2trk));
}
// ------------------------------------------------------

bool ParticleDecayId::DetectDecay(
    const std::vector<recob::Wire> & wires,
    const std::vector< art::Ptr<recob::Hit> > & hits,
    std::map< size_t, TVector3 > & spoints,
    std::vector< std::pair<TVector3, double> > & result)
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
    }

    std::vector< std::pair<size_t, float> > candidates2d[nviews];
    std::vector< std::pair<TVector3, float> > candidates3d[nviews];
    for (size_t v = 0; v < nviews; ++v)
    {
        size_t idx = 0;
        while (idx < outputs[v].size())
        {
            if (outputs[v][idx] > fRoiThreshold)
            {
                size_t ci = idx;
                float max = outputs[v][idx];
                ++idx;

                while ((idx < outputs[v].size()) && (outputs[v][idx] > fRoiThreshold))
                {
                    if (outputs[v][idx] > max) { max = outputs[v][idx]; ci = idx; }
                    ++idx;
                }
                candidates2d[v].emplace_back(ci, max);
                candidates3d[v].emplace_back(spoints[wire_drift[v][ci].key()], max);
            }
            else ++idx;
        }
    }

    double min_dist = 2.0; // [cm], threshold for today to distinguish between two different candidates,
                           // if belo threshold, then use 3D point corresponding to higher cnn output

    // need coincidence of high cnn out in two views, then look if there is another close candidate
    // and again select by cnn output value, would like to have few strong candidates
    bool found = false;
    for (size_t v = 0; v < nviews - 1; ++v)
    {
        for (size_t i = 0; i < candidates3d[v].size(); ++i)
        {
            TVector3 c0(candidates3d[v][i].first);
            float p0 = candidates3d[v][i].second;

            for (size_t u = v + 1; u < nviews; ++u)
            {
                for (size_t j = 0; j < candidates3d[v].size(); ++j)
                {
                    TVector3 c1(candidates3d[v][j].first);
                    float p1 = candidates3d[v][j].second;

                    if ((c0 - c1).Mag() < min_dist)
                    {
                        TVector3 c(c0);
                        if (p1 > p0) { c = c1; }
                        double p = p0 * p1;

                        if (p > fPointThreshold)
                        {
                            double d, dmin = min_dist;
                            size_t kmin = 0;
                            for (size_t k = 0; k < result.size(); ++k)
                            {
                                d = (result[k].first - c).Mag();
                                if (d < dmin)
                                {
                                    dmin = d; kmin = k;
                                }
                            }
                            if (dmin < min_dist)
                            {
                                if (result[kmin].second < p) // replace previously found point
                                {
                                    result[kmin].first = c;
                                    result[kmin].second = p;
                                    found = true;
                                }
                            }
                            else // nothing close in the list, add new point
                            {
                                result.emplace_back(c, p);
                                found = true;
                            }
                        } // if (p > fPointThreshold)
                    } // coincidence: points from views u and v are close
                } // loop over points in view u
            } // loop over views u
        } // loop over points in view v
    } // loop over views v
    return found;
}

DEFINE_ART_MODULE(ParticleDecayId)

}

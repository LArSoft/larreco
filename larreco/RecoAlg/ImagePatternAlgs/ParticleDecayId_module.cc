/////////////////////////////////////////////////////////////////////////////////
// Class:       ParticleDecayId
// Module Type: producer
// File:        ParticleDecayId_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// THIS IS JUST A PLACEHOLDER NOW - PLEASE WAIT A FEW DAYS FOR THE WORKING CODE
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
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"

#include <memory>

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
	//produces< std::vector<recob::Cluster> >();
	//produces< art::Assns<recob::Cluster, recob::Hit> >();
}

// ------------------------------------------------------

void ParticleDecayId::produce(art::Event & evt)
{

}
// ------------------------------------------------------

DEFINE_ART_MODULE(ParticleDecayId)

}


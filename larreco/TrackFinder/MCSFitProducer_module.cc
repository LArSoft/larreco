/// \class MCSFitProducer
///
/// \brief Producer 
///
/// \author G. Cerati, based on code from L. Kalousis and D. Kaleko
///

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

#include <memory>

namespace trkf {

  class MCSFitProducer : public art::EDProducer {
  public:

    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputTrajectoryLabel {
	Name("inputTrajectoryLabel"),
	Comment("Label of recob::TrackTrajectory Collection to be fit")
      };
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<MCSFitProducer::Inputs> inputs {
	Name("inputs"),
      };
      fhicl::Table<TrajectoryMCSFitter::Config> fitter {
	Name("fitter")
      };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit MCSFitProducer(Parameters const & p);
    ~MCSFitProducer();

    // Plugins should not be copied or assigned.
    MCSFitProducer(MCSFitProducer const &) = delete;
    MCSFitProducer(MCSFitProducer &&) = delete;
    MCSFitProducer & operator = (MCSFitProducer const &) = delete;
    MCSFitProducer & operator = (MCSFitProducer &&) = delete;

    void produce(art::Event & e) override;

  private:
    Parameters p_;
    art::InputTag trajectoryInputTag;
    TrajectoryMCSFitter mcsfitter;
  };
}

trkf::MCSFitProducer::MCSFitProducer(trkf::MCSFitProducer::Parameters const & p)
  : p_(p), mcsfitter(p_().fitter)
{
  produces<std::vector<recob::MCSFitResult> >();
}

trkf::MCSFitProducer::~MCSFitProducer() {}

void trkf::MCSFitProducer::produce(art::Event & e)
{
  //
  auto output  = std::make_unique<std::vector<recob::MCSFitResult> >();
  //
  art::Handle<std::vector<recob::TrackTrajectory> > inputTrajectoryH;
  bool ok = e.getByLabel(trajectoryInputTag,inputTrajectoryH);
  if (!ok) throw cet::exception("MCSFitProducer") << "Cannot find recob::TrackTrajectory art::Handle with inputTag " << trajectoryInputTag;
  const std::vector<recob::TrackTrajectory>& trajectoryVec = *(inputTrajectoryH.product());
  for (const auto& inTraj : trajectoryVec) {
    //fit
    recob::MCSFitResult result = mcsfitter.fitMcs(inTraj);
    output->emplace_back(std::move(result));
  }
  e.put(std::move(output));
}

DEFINE_ART_MODULE(trkf::MCSFitProducer)

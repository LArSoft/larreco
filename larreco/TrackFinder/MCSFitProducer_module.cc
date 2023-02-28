#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

#include <memory>

namespace trkf {
  /**
   * @file  MCSFitProducer_module.cc
   * @class trkf::MCSFitProducer
   *
   * @brief Producer for TrajectoryMCSFitter.
   *
   * Producer for TrajectoryMCSFitter, which performs a Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory.
   * It reads a recob::Track collection and produces a collection of recob::MCSFitResult where the elements are in the same order as the input collection (no explicit association is written).
   *
   * For configuration options see MCSFitProducer#Inputs and MCSFitProducer#Config
   *
   * @author  G. Cerati (FNAL, MicroBooNE), based on code from L. Kalousis and D. Kaleko
   * @date    2017
   * @version 1.0
   */
  class MCSFitProducer : public art::EDProducer {
  public:
    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputLabel{
        Name("inputLabel"),
        Comment("Label of recob::TrackTrajectory Collection to be fit")};
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<MCSFitProducer::Inputs> inputs{
        Name("inputs"),
      };
      fhicl::Table<TrajectoryMCSFitter::Config> fitter{Name("fitter")};
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit MCSFitProducer(Parameters const& p);
    ~MCSFitProducer();

    // Plugins should not be copied or assigned.
    MCSFitProducer(MCSFitProducer const&) = delete;
    MCSFitProducer(MCSFitProducer&&) = delete;
    MCSFitProducer& operator=(MCSFitProducer const&) = delete;
    MCSFitProducer& operator=(MCSFitProducer&&) = delete;

  private:
    void produce(art::Event& e) override;

    Parameters p_;
    art::InputTag inputTag;
    TrajectoryMCSFitter mcsfitter;
  };
}

trkf::MCSFitProducer::MCSFitProducer(trkf::MCSFitProducer::Parameters const& p)
  : EDProducer{p}, p_(p), mcsfitter(p_().fitter)
{
  inputTag = art::InputTag(p_().inputs().inputLabel());
  produces<std::vector<recob::MCSFitResult>>();
}

trkf::MCSFitProducer::~MCSFitProducer() {}

void trkf::MCSFitProducer::produce(art::Event& e)
{
  //
  auto output = std::make_unique<std::vector<recob::MCSFitResult>>();
  //
  art::Handle<std::vector<recob::Track>> inputH;
  bool ok = e.getByLabel(inputTag, inputH);
  if (!ok)
    throw cet::exception("MCSFitProducer")
      << "Cannot find input art::Handle with inputTag " << inputTag;
  const auto& inputVec = *(inputH.product());
  for (const auto& element : inputVec) {
    //fit
    recob::MCSFitResult result = mcsfitter.fitMcs(element);
    output->emplace_back(std::move(result));
  }
  e.put(std::move(output));
}

DEFINE_ART_MODULE(trkf::MCSFitProducer)

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

#include "larreco/RecoAlg/TauMCSFitter.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/OptionalTable.h"
#include <string>

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
   *
   * Supports selection between the standard TrajectoryMCSFitter
   * and TauMCSFitter through FHiCL configuration.
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

   // Configuration for selecting the MCS fitting algorithm.
   struct Config {
     using Name = fhicl::Name;
     using Comment = fhicl::Comment;

    fhicl::Atom<std::string> algorithm{
      Name("algorithm"),
      Comment("MCS fitting algorithm: 'standard' or 'tau'"),
      "standard"
    };

    fhicl::Table<MCSFitProducer::Inputs> inputs{
      Name("inputs")
    };

    fhicl::Table<TrajectoryMCSFitter::Config> fitter{
      Name("fitter")
    };

    fhicl::OptionalTable<TauMCSFitter::Config> tauFitter{
      Name("tauFitter"),
      Comment("Configuration for double gaussian MCS fitter")
    };
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
    std::string algorithm;
    TrajectoryMCSFitter mcsfitter;
    std::unique_ptr<TauMCSFitter> taumcsfitter;
  };
}

trkf::MCSFitProducer::MCSFitProducer(
  trkf::MCSFitProducer::Parameters const& p)
  : EDProducer{p}
  , p_(p)
  , inputTag(p_().inputs().inputLabel())
  , algorithm(p_().algorithm())
  , mcsfitter(p_().fitter)
{
  if (algorithm == "tau") {

    auto tauConfig = p_().tauFitter();

    if (!tauConfig) {
      throw cet::exception("MCSFitProducer")
        << "algorithm is set to 'tau', but no tauFitter "
        << "configuration was provided.";
    }

    taumcsfitter =
      std::make_unique<TauMCSFitter>(*tauConfig);
  }
  else if (algorithm != "standard") {
    throw cet::exception("MCSFitProducer")
      << "Unknown MCS fitting algorithm '" << algorithm
      << "'. Valid options are 'standard' and 'tau'.";
  }

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
     auto result =
       (algorithm == "standard")
       ? mcsfitter.fitMcs(element)
       : taumcsfitter->fitMcs(element);

     output->emplace_back(std::move(result));

  }
  e.put(std::move(output));
}

DEFINE_ART_MODULE(trkf::MCSFitProducer)

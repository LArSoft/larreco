////////////////////////////////////////////////////////////////////////
//
// \file HoughLineFinder_module.cc
//
// \author kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//
// \file HoughLineFinder.cxx
//
// \author joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters found by DBSCAN
//  after deconvolution and hit finding.
//
//  The algorithm is based on:
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", Machine
//  Vision and Applications 3, 87 (1990)
////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nurandom
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/HoughBaseAlg.h"

namespace cluster {

  class HoughLineFinder : public art::EDProducer {
  public:
    explicit HoughLineFinder(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    std::string fDBScanModuleLabel;
    unsigned int fHoughSeed;
    HoughBaseAlg fHLAlg; ///< object that does the Hough Transform
    CLHEP::HepRandomEngine& fEngine;
  };

  //------------------------------------------------------------------------------
  HoughLineFinder::HoughLineFinder(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fDBScanModuleLabel{pset.get<std::string>("DBScanModuleLabel")}
    , fHoughSeed{pset.get<unsigned int>("HoughSeed", 0)}
    , fHLAlg(pset.get<fhicl::ParameterSet>("HoughBaseAlg"))
    // Create random number engine needed for PPHT; obtain the random seed from
    // NuRandomService, unless overridden in configuration with key "Seed" remember that
    // HoughSeed will override this on each event if specified
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(createEngine(0),
                                                                                 pset,
                                                                                 "Seed"))
  {
    produces<std::vector<recob::Cluster>>();
    produces<art::Assns<recob::Cluster, recob::Hit>>();
  }

  //------------------------------------------------------------------------------
  void HoughLineFinder::produce(art::Event& evt)
  {
    //////////////////////////////////////////////////////
    // here is how to get a collection of objects out of the file
    // and connect it to a art::Handle
    //////////////////////////////////////////////////////
    auto const clusterListHandle =
      evt.getValidHandle<std::vector<recob::Cluster>>(fDBScanModuleLabel);

    std::vector<art::Ptr<recob::Cluster>> clusIn;
    clusIn.reserve(clusterListHandle->size());
    for (unsigned int ii = 0; ii < clusterListHandle->size(); ++ii) {
      clusIn.emplace_back(clusterListHandle, ii);
    }

    //Point to a collection of clusters to output.
    auto ccol = std::make_unique<std::vector<recob::Cluster>>();
    auto assn = std::make_unique<art::Assns<recob::Cluster, recob::Hit>>();

    // make a std::vector< art::PtrVector<recob::Hit> >
    // to hold the associated hits of the Hough Transform
    std::vector<art::PtrVector<recob::Hit>> clusHitsOut;

    // If a nonzero random number seed has been provided,
    // overwrite the seed already initialized
    if (fHoughSeed != 0) { fEngine.setSeed(fHoughSeed, 0); }

    size_t const numclus =
      fHLAlg.FastTransform(clusIn, *ccol, clusHitsOut, fEngine, evt, fDBScanModuleLabel);

    MF_LOG_DEBUG("HoughLineClusters") << "found " << numclus << "clusters with HoughBaseAlg";

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "HoughLineFinder Summary:";
    for (size_t i = 0; i < ccol->size(); ++i) {
      mf::LogVerbatim("Summary") << ccol->at(i);

      // associate the hits to this cluster
      util::CreateAssn(*this, evt, *ccol, clusHitsOut[i], *assn, i);
    }

    evt.put(move(ccol));
    evt.put(move(assn));
  }

} // end namespace

DEFINE_ART_MODULE(cluster::HoughLineFinder)

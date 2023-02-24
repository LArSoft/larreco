////////////////////////////////////////////////////////////////////////
//
// \file DBcluster_module.cc
//
// \author kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////

//Framework includes:
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/DBScanAlg.h"

#include "TH1.h"
#include <cstdlib>
#include <iomanip>

namespace cluster {

  //---------------------------------------------------------------
  class DBcluster : public art::EDProducer {
  public:
    explicit DBcluster(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt);
    void beginJob();

    TH1F* fhitwidth;
    TH1F* fhitwidth_ind_test;
    TH1F* fhitwidth_coll_test;

    std::string fhitsModuleLabel;

    DBScanAlg fDBScan; ///< object that implements the DB scan algorithm
  };

}

namespace cluster {

  //-------------------------------------------------
  DBcluster::DBcluster(fhicl::ParameterSet const& pset)
    : EDProducer{pset}, fDBScan(pset.get<fhicl::ParameterSet>("DBScanAlg"))
  {
    fhitsModuleLabel = pset.get<std::string>("HitsModuleLabel");

    produces<std::vector<recob::Cluster>>();
    produces<art::Assns<recob::Cluster, recob::Hit>>();
  }

  //-------------------------------------------------
  void DBcluster::beginJob()
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    fhitwidth = tfs->make<TH1F>(" fhitwidth", "width of hits in cm", 50000, 0, 5);
    fhitwidth_ind_test = tfs->make<TH1F>("fhitwidth_ind_test", "width of hits in cm", 50000, 0, 5);
    fhitwidth_coll_test =
      tfs->make<TH1F>("fhitwidth_coll_test", "width of hits in cm", 50000, 0, 5);
  }

  //-----------------------------------------------------------------
  void DBcluster::produce(art::Event& evt)
  {

    //get a collection of clusters
    std::unique_ptr<std::vector<recob::Cluster>> ccol(new std::vector<recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>> assn(
      new art::Assns<recob::Cluster, recob::Hit>);

    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here; configuration would happen here,
    // but we are using the default configuration for that algorithm
    ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;

    art::ServiceHandle<geo::Geometry const> geom;
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

    art::Handle<std::vector<recob::Hit>> hitcol;
    evt.getByLabel(fhitsModuleLabel, hitcol);

    // loop over all hits in the event and look for clusters (for each plane)
    std::vector<art::Ptr<recob::Hit>> allhits;

    // get channel quality service:
    lariov::ChannelStatusProvider const& channelStatus =
      art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    lariov::ChannelStatusProvider::ChannelSet_t const BadChannels = channelStatus.BadChannels();

    // make a map of the geo::PlaneID to vectors of art::Ptr<recob::Hit>
    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit>>> planeIDToHits;
    for (size_t i = 0; i < hitcol->size(); ++i)
      planeIDToHits[hitcol->at(i).WireID().planeID()].push_back(art::Ptr<recob::Hit>(hitcol, i));

    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);
    util::GeometryUtilities const gser{*geom, wireReadoutGeom, clock_data, det_prop};
    for (auto& itr : planeIDToHits) {

      geo::SigType_t sigType = wireReadoutGeom.SignalType(itr.first);
      allhits.resize(itr.second.size());
      allhits.swap(itr.second);

      fDBScan.InitScan(clock_data, det_prop, allhits, BadChannels);

      //----------------------------------------------------------------
      for (unsigned int j = 0; j < fDBScan.fps.size(); ++j) {

        if (allhits.size() != fDBScan.fps.size()) break;

        fhitwidth->Fill(fDBScan.fps[j][2]);

        if (sigType == geo::kInduction) fhitwidth_ind_test->Fill(fDBScan.fps[j][2]);
        if (sigType == geo::kCollection) fhitwidth_coll_test->Fill(fDBScan.fps[j][2]);
      }

      //*******************************************************************
      fDBScan.run_cluster();

      for (size_t i = 0; i < fDBScan.fclusters.size(); ++i) {
        art::PtrVector<recob::Hit> clusterHits;

        for (size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j) {
          if (fDBScan.fpointId_to_clusterId[j] == i) { clusterHits.push_back(allhits[j]); }
        }

        if (clusterHits.empty()) continue;

        /// \todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
        const geo::WireID& wireID = clusterHits.front()->WireID();
        unsigned int sw = wireID.Wire;
        unsigned int ew = clusterHits.back()->WireID().Wire;

        // feed the algorithm with all the cluster hits
        ClusterParamAlgo.ImportHits(gser, clusterHits);

        // create the recob::Cluster directly in the vector
        ClusterCreator cluster(gser,
                               ClusterParamAlgo,                     // algo
                               float(sw),                            // start_wire
                               0.,                                   // sigma_start_wire
                               clusterHits.front()->PeakTime(),      // start_tick
                               clusterHits.front()->SigmaPeakTime(), // sigma_start_tick
                               float(ew),                            // end_wire
                               0.,                                   // sigma_end_wire,
                               clusterHits.back()->PeakTime(),       // end_tick
                               clusterHits.back()->SigmaPeakTime(),  // sigma_end_tick
                               ccol->size(),                         // ID
                               clusterHits.front()->View(),          // view
                               wireID.planeID(),                     // plane
                               recob::Cluster::Sentry                // sentry
        );

        ccol->emplace_back(cluster.move());

        // associate the hits to this cluster
        util::CreateAssn(evt, *ccol, clusterHits, *assn);

      } //end loop over fclusters

      allhits.clear();
    } // end loop over PlaneID map

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "DBcluster Summary:";
    for (unsigned int i = 0; i < ccol->size(); ++i)
      mf::LogVerbatim("Summary") << ccol->at(i);

    evt.put(std::move(ccol));
    evt.put(std::move(assn));
  }

} // end namespace

DEFINE_ART_MODULE(cluster::DBcluster)

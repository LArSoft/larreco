////////////////////////////////////////////////////////////////////////
//
// ShowerCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#include <string>

// ROOT includes

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace shwf {
  class ShowerCheater : public art::EDProducer {
  public:
    explicit ShowerCheater(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt);

    std::string fCheatedClusterLabel; ///< label for module creating recob::Cluster objects
    std::string fG4ModuleLabel;       ///< label for module running G4 and making particles, etc
  };
}

namespace shwf {

  //--------------------------------------------------------------------
  ShowerCheater::ShowerCheater(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    fCheatedClusterLabel = pset.get<std::string>("CheatedClusterLabel", "cluster");
    fG4ModuleLabel = pset.get<std::string>("G4ModuleLabel", "largeant");

    produces<std::vector<recob::Shower>>();
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Shower, recob::Cluster>>();
    produces<art::Assns<recob::Shower, recob::SpacePoint>>();
    produces<art::Assns<recob::Shower, recob::Hit>>();
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();
  }

  //--------------------------------------------------------------------
  void
  ShowerCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    art::ServiceHandle<geo::Geometry const> geo;

    // grab the clusters that have been reconstructed
    art::Handle<std::vector<recob::Cluster>> clustercol;
    evt.getByLabel(fCheatedClusterLabel, clustercol);

    art::FindManyP<recob::Hit> fmh(clustercol, evt, fCheatedClusterLabel);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector<art::Ptr<recob::Cluster>> clusters;
    art::fill_ptr_vector(clusters, clustercol);

    // make a map of vectors of art::Ptrs keyed by eveID values
    std::map<int, std::vector<std::pair<size_t, art::Ptr<recob::Cluster>>>> eveClusterMap;

    // loop over all clusters and fill in the map
    for (size_t c = 0; c < clusters.size(); ++c) {

      // in the ClusterCheater module we set the cluster ID to be
      // the eve particle track ID*1000 + plane*100 + tpc*10 + cryostat number.  The
      // floor function on the cluster ID / 1000 will give us
      // the eve track ID
      int eveID = floor(clusters[c]->ID() / 1000.);

      std::pair<size_t, art::Ptr<recob::Cluster>> indexPtr(c, clusters[c]);

      eveClusterMap[eveID].push_back(indexPtr);

    } // end loop over clusters

    // loop over the map and make prongs
    std::unique_ptr<std::vector<recob::Shower>> showercol(new std::vector<recob::Shower>);
    std::unique_ptr<std::vector<recob::SpacePoint>> spcol(new std::vector<recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Shower, recob::Cluster>> scassn(
      new art::Assns<recob::Shower, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Shower, recob::Hit>> shassn(
      new art::Assns<recob::Shower, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Shower, recob::SpacePoint>> sspassn(
      new art::Assns<recob::Shower, recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Hit, recob::SpacePoint>> sphassn(
      new art::Assns<recob::Hit, recob::SpacePoint>);

    for (auto const& clusterMapItr : eveClusterMap) {

      // separate out the hits for each particle into the different views
      std::vector<std::pair<size_t, art::Ptr<recob::Cluster>>> const& eveClusters =
        clusterMapItr.second;

      size_t startSPIndx = spcol->size();

      double totalCharge = 0.;

      std::vector<art::Ptr<recob::Cluster>> ptrvs;
      std::vector<size_t> idxs;

      for (auto const& idxPtr : eveClusters) {
        idxs.push_back(idxPtr.first);
        ptrvs.push_back(idxPtr.second);

        // need to make the space points for this prong
        // loop over the hits for this cluster and make
        // a space point for each one
        // set the SpacePoint ID to be the cluster ID*10000
        // + the hit index in the cluster PtrVector of hits
        std::vector<art::Ptr<recob::Hit>> const& hits = fmh.at(idxPtr.first);

        for (size_t h = 0; h < hits.size(); ++h) {
          art::Ptr<recob::Hit> hit = hits[h];
          // add up the charge from the hits on the collection plane
          if (hit->SignalType() == geo::kCollection) totalCharge += hit->Integral();
          std::vector<double> xyz = bt_serv->HitToXYZ(hit);
          double sperr[6] = {0.01, 0.01, 0.1, 0.001, 0.001, 0.001};

          // make the space point and set its ID and XYZ
          // the errors and chi^2 are set to "good" values as we know the information perfectly
          recob::SpacePoint sp(&xyz[0], sperr, 0.9, idxPtr.second->ID() * 10000 + h);
          spcol->push_back(sp);

          // associate the space point to the hit
          util::CreateAssn(*this, evt, *spcol, hit, *sphassn);

        } // end loop over hits
      }   // end loop over pairs of index values and cluster Ptrs

      size_t endSPIndx = spcol->size();

      // is this prong electro-magnetic in nature or
      // hadronic/muonic?  EM --> shower, everything else is a track
      if (std::abs(pi_serv->ParticleList()[clusterMapItr.first]->PdgCode()) == 11 ||
          std::abs(pi_serv->ParticleList()[clusterMapItr.first]->PdgCode()) == 22 ||
          std::abs(pi_serv->ParticleList()[clusterMapItr.first]->PdgCode()) == 111) {

        mf::LogInfo("ShowerCheater")
          << "prong of " << clusterMapItr.first << " is a shower with pdg code "
          << pi_serv->ParticleList()[clusterMapItr.first]->PdgCode();

        // get the direction cosine for the eve ID particle
        // just use the same for both the start and end of the prong
        const TLorentzVector initmom = pi_serv->ParticleList()[clusterMapItr.first]->Momentum();
        TVector3 dcos(
          initmom.Px() / initmom.Mag(), initmom.Py() / initmom.Mag(), initmom.Pz() / initmom.Mag());
        TVector3 dcosErr(1.e-3, 1.e-3, 1.e-3);

        /// \todo figure out the max transverse width of the shower in the x and y directions
        //double maxTransWidth[2] = { util::kBogusD, util::kBogusD };
        //double distanceMaxWidth = util::kBogusD;

        // add a prong to the collection.  Make the prong
        // ID be the same as the track ID for the eve particle
        recob::Shower s;
        s.set_id(showercol->size());
        s.set_direction(dcos);
        s.set_direction_err(dcosErr);
        /*
	showercol->push_back(recob::Shower(dcos, dcosErr, maxTransWidth,
					   distanceMaxWidth, totalCharge, clusterMapItr.first));
	*/
        showercol->push_back(s);
        // associate the shower with its clusters
        util::CreateAssn(*this, evt, *showercol, ptrvs, *scassn);

        // get the hits associated with each cluster and associate those with the shower
        for (size_t i = 0; i < idxs.size(); ++i) {
          std::vector<art::Ptr<recob::Hit>> hits = fmh.at(i);
          util::CreateAssn(*this, evt, *showercol, hits, *shassn);
        }

        // associate the shower with the space points
        util::CreateAssn(*this, evt, *showercol, *spcol, *sspassn, startSPIndx, endSPIndx);

        mf::LogInfo("ShowerCheater") << "adding shower: \n"
                                     << showercol->back() << "\nto collection.";

      } // end if this is a shower
    }   // end loop over the map

    evt.put(std::move(showercol));
    evt.put(std::move(spcol));
    evt.put(std::move(scassn));
    evt.put(std::move(shassn));
    evt.put(std::move(sspassn));
    evt.put(std::move(sphassn));

    return;

  } // end produce

} // end namespace

namespace shwf {

  DEFINE_ART_MODULE(ShowerCheater)

}

////////////////////////////////////////////////////////////////////////
//
// DBSCANfinderAna class
//
// \author kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <string>

#include <TH1F.h>
#include <TH2F.h>

//Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nug4/ParticleNavigation/EmEveIdCalculator.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Core/EDAnalyzer.h"

///Cluster finding and building
namespace cluster {

  class DBclusterAna : public art::EDAnalyzer {

  public:
    explicit DBclusterAna(fhicl::ParameterSet const& pset);
    virtual ~DBclusterAna();

  private:
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

    TH1F* fNoParticles_pdg;
    TH1F* fNoParticles_trackid;
    TH1F* fNoParticles_trackid_mother;
    TH1F* fNoParticles_trackid_per_event;
    TH1F* fNoParticles_pdg_per_event;
    TH1F* fCl_for_Muon;
    TH1F* fNoClustersInEvent;
    TH1F* fPercentNoise;
    TH1F* fno_of_clusters_per_track;
    TH1F* fPercent_lost_muon_hits;
    TH1F* fPercent_lost_electron_hits;
    TH1F* fPercent_lost_positron_hits;
    TH1F* fPercent_lost_111_hits;
    TH1F* fPercent_lost_211_hits;
    TH1F* fPercent_lost_m211_hits;
    TH1F* fPercent_lost_2212_hits;
    TH1F* fPercent_lost_2112_hits;

    TH1F* fPercent_lost_muon_energy;
    TH1F* fPercent_lost_electron_energy;
    TH1F* fPercent_lost_positron_energy;
    TH1F* fPercent_lost_111_energy;
    TH1F* fPercent_lost_211_energy;
    TH1F* fPercent_lost_m211_energy;
    TH1F* fPercent_lost_2212_energy;
    TH1F* fPercent_lost_2112_energy;
    TH1F* fEnergy;
    TH2F* fbrian_in;
    TH2F* fbrian_coll;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fClusterFinderModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;

  }; // class DBclusterAna

}

namespace cluster {

  //--------------------------------------------------------------------
  DBclusterAna::DBclusterAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fDigitModuleLabel(pset.get<std::string>("DigitModuleLabel"))
    , fHitsModuleLabel(pset.get<std::string>("HitsModuleLabel"))
    , fLArG4ModuleLabel(pset.get<std::string>("LArGeantModuleLabel"))
    , fClusterFinderModuleLabel(pset.get<std::string>("ClusterFinderModuleLabel"))
    , fCalDataModuleLabel(pset.get<std::string>("CalDataModuleLabel"))
    , fGenieGenModuleLabel(pset.get<std::string>("GenieGenModuleLabel"))
  {}

  //------------------------------------------------------------------
  DBclusterAna::~DBclusterAna() {}

  void
  DBclusterAna::beginJob()
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    fNoParticles_pdg_per_event = tfs->make<TH1F>(
      "fNoParticles_pdg_per_event", "Average # of Particles per cluster for each event", 500, 0, 5);
    fNoParticles_pdg = tfs->make<TH1F>(
      "fNoParticles_pdg", "Number of Particles in a Cluster for each cluster", 500, 0, 5);
    fNoParticles_trackid = tfs->make<TH1F>(
      "fNoParticles_trackid", "Number of different TrackIDs in a Cluster", 300, 0, 30);

    fNoParticles_trackid_mother =
      tfs->make<TH1F>("fNoParticles_trackid_mother",
                      "Number of different TrackIDs in a Cluster(using mother)for each cluster",
                      300,
                      0,
                      30);

    fNoParticles_trackid_per_event =
      tfs->make<TH1F>("fNoParticles_trackid_per_event",
                      "Avg Number of different TrackIDs per Cluster per event",
                      300,
                      0,
                      30);
    fCl_for_Muon =
      tfs->make<TH1F>("fCl_for_Muon", "Number of Clusters for Muon per plane (pdg)", 1500, 0, 15);

    fNoClustersInEvent =
      tfs->make<TH1F>("fNoClustersInEvent", "Number of Clusters in an Event", 50, 0, 50);

    fPercentNoise =
      tfs->make<TH1F>("fPercentNoise", "% of hits that were marked as Noise by DBSCAN", 250, 0, 25);

    fno_of_clusters_per_track = tfs->make<TH1F>(
      "fno_of_clusters_per_track", "Number of Clusters per TrackID per plane", 1500, 0, 15);

    fPercent_lost_muon_hits =
      tfs->make<TH1F>("fPercent_lost_muon_hits",
                      "Number of muon hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_electron_hits =
      tfs->make<TH1F>("fPercent_lost_electron_hits",
                      "Number of electron hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_positron_hits =
      tfs->make<TH1F>("fPercent_lost_positron_hits",
                      "Number of positron hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_111_hits =
      tfs->make<TH1F>("fPercent_lost_111_hits",
                      "Number of pion(111) hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_211_hits =
      tfs->make<TH1F>("fPercent_lost_211_hits",
                      "Number of pion(211) hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_m211_hits =
      tfs->make<TH1F>("fPercent_lost_m211_hits",
                      "Number of pion(-211) hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_2212_hits =
      tfs->make<TH1F>("fPercent_lost_2212_hits",
                      "Number of proton hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_2112_hits =
      tfs->make<TH1F>("fPercent_lost_2112_hits",
                      "Number of neutron hits excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);

    fPercent_lost_muon_energy = tfs->make<TH1F>("fPercent_lost_muon_energy",
                                                " muon energy excluded by dbscan in % (per Event)",
                                                10000,
                                                0,
                                                100);
    fPercent_lost_electron_energy =
      tfs->make<TH1F>("fPercent_lost_electron_energy",
                      "electron energy excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_positron_energy =
      tfs->make<TH1F>("fPercent_lost_positron_energy",
                      " positron energy excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_111_energy =
      tfs->make<TH1F>("fPercent_lost_111_energy",
                      "pion(111) energy excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_211_energy =
      tfs->make<TH1F>("fPercent_lost_211_energy",
                      "pion(211) energy excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_m211_energy =
      tfs->make<TH1F>("fPercent_lost_m211_energy",
                      " pion(-211) energy excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);
    fPercent_lost_2212_energy = tfs->make<TH1F>("fPercent_lost_2212_energy",
                                                "proton energy excluded by dbscan in % (per Event)",
                                                10000,
                                                0,
                                                100);
    fPercent_lost_2112_energy =
      tfs->make<TH1F>("fPercent_lost_2112_energy",
                      "neutron energy excluded by dbscan in % (per Event)",
                      10000,
                      0,
                      100);

    fEnergy = tfs->make<TH1F>("fEnergy", "energy for each voxel", 100000, 0, 0.0005);

    fbrian_in = tfs->make<TH2F>("fbrian_in",
                                ";# Electrons deposited; # Electrons detected by hitfinder",
                                1000,
                                0,
                                10000000,
                                1000,
                                0,
                                10000000);
    fbrian_coll = tfs->make<TH2F>("fbrian_coll",
                                  ";# Electrons deposited; # Electrons detected by hitfinder",
                                  1000,
                                  0,
                                  10000000,
                                  1000,
                                  0,
                                  10000000);
  }

  void
  DBclusterAna::analyze(const art::Event& evt)
  {
    std::cout << "run    : " << evt.run() << std::endl;
    std::cout << "event  : " << evt.id().event() << std::endl;
    //----------------------------------------------------------------

    /* This is basically a module for studying MC efficiency/purity. Kick out now if not MC. EC, 8-Oct-2010 */
    if (evt.isRealData()) {
      std::cout << "**** DBclusterAna: Bailing. Don't call this module if you're not MC. "
                << std::endl;
      exit(1);
    }

    art::ServiceHandle<geo::Geometry const> geom;
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    art::Handle<std::vector<raw::RawDigit>> rdListHandle;
    evt.getByLabel(fDigitModuleLabel, rdListHandle);
    art::Handle<std::vector<recob::Hit>> hitListHandle;
    evt.getByLabel(fHitsModuleLabel, hitListHandle);
    art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
    evt.getByLabel(fGenieGenModuleLabel, mctruthListHandle);
    art::Handle<std::vector<recob::Cluster>> clusterListHandle;
    evt.getByLabel(fClusterFinderModuleLabel, clusterListHandle);
    art::Handle<std::vector<recob::Wire>> wireListHandle;
    evt.getByLabel(fCalDataModuleLabel, wireListHandle);

    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterFinderModuleLabel);

    //----------------------------------------------------------------

    //----------------------------------------------------------------

    art::PtrVector<raw::RawDigit> rawdigits;

    for (size_t ii = 0; ii < rdListHandle->size(); ++ii) {
      art::Ptr<raw::RawDigit> rawdigit(rdListHandle, ii);
      rawdigits.push_back(rawdigit);
    }

    //get the simb::MCParticle collection from the art::Event and then use the
    //Simulation/SimListUtils object to create a sim::ParticleList from the art::Event.
    pi_serv->SetEveIdCalculator(new sim::EmEveIdCalculator);

    sim::ParticleList const& _particleList = pi_serv->ParticleList();

    std::vector<int> mc_trackids;

    for (auto i = _particleList.begin(); i != _particleList.end(); ++i) {
      int trackID = (*i).first;
      mc_trackids.push_back(trackID);
    }

    //----------------------------------------------------------------
    art::PtrVector<simb::MCTruth> mclist;
    for (size_t ii = 0; ii < mctruthListHandle->size(); ++ii) {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle, ii);
      mclist.push_back(mctparticle);
    }

    //--------------------------------------------------
    std::vector<art::Ptr<recob::Hit>> hits_vec;

    //--------------------------------------------------
    std::vector<art::Ptr<recob::Hit>> hits;
    art::fill_ptr_vector(hits, hitListHandle);
    //---------------------------------------------------

    art::PtrVector<recob::Cluster> clusters;
    for (unsigned int ii = 0; ii < clusterListHandle->size(); ++ii) {
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle, ii);
      clusters.push_back(clusterHolder);
    }

    std::cout << "in Efficiency, clusters.size()= " << clusters.size() << std::endl;

    //---------------------------------------------------------------
    art::PtrVector<recob::Wire> wirelist;

    for (size_t ii = 0; ii < wireListHandle->size(); ++ii) {
      art::Ptr<recob::Wire> wireHolder(wireListHandle, ii);

      wirelist.push_back(wireHolder);
    }

    //...........................................................................
    // How many different particles do we have in a cluster???
    //  -- I will answer that by 3 methods
    //  1) count different TrackIDs in each cluster
    //  2)  count different TrackIDs in each cluster with the use of "Mother" to get rid of the ones that just randomly got their TrackIDs changed by Geant4
    //  3) count different pdg codes in each cluster <--probably not very good b/c if you have 2 "different" electrons in a cluster it's going to count them as one.
    //..........................................................................
    // How many clusters it takes to contain a particle???
    // take each TrackID and count in how many clusters it appears
    //.........................................................................

    double no_of_particles_in_cluster = 0;
    double sum_vec_trackid = 0;
    double no_of_clusters = 0;
    double total_no_hits_in_clusters = 0;
    std::vector<int> vec_pdg;
    std::vector<int> vec_trackid, vec_trackid_mother, vec_trackid_mother_en;
    std::vector<int> all_trackids;
    std::vector<int> ids;
    std::vector<int>::iterator it, it2, it3, it4, it5, it6, it7, it8;
    int no_cl_for_muon = 0;
    int no_cl_for_electron = 0;
    int no_cl_for_positron = 0;
    int no_cl_for_pion_111 = 0;
    int no_cl_for_pion_211 = 0;
    int no_cl_for_pion_m211 = 0;
    int no_cl_for_proton = 0;
    double noCluster = 0;
    double _hit_13 = 0, _hit_11 = 0, _hit_m_11 = 0, _hit_111 = 0, _hit_211 = 0, _hit_m211 = 0,
           _hit_2212 = 0, _hit_2112 = 0;
    double _en_13 = 0, _en_11 = 0, _en_m11 = 0, _en_111 = 0, _en_211 = 0, _en_m211 = 0,
           _en_2212 = 0, _en_2112 = 0;
    std::vector<double> diff_vec;

    double hit_energy = 0;
    double total_Q_cluster_hits = 0;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    if (clusters.size() != 0 && hits.size() != 0) {
      for (unsigned int plane = 0; plane < geom->Nplanes(); ++plane) {
        geo::View_t view = geom->Plane(plane).View();
        for (size_t j = 0; j < clusters.size(); ++j) {

          if (clusters[j]->View() == view) {

            std::vector<art::Ptr<recob::Hit>> _hits = fmh.at(j);
            art::Ptr<recob::Hit> _hits_ptr;

            for (size_t p = 0; p < _hits.size(); ++p) {
              _hits_ptr = _hits[p];
              hits_vec.push_back(_hits_ptr);
              total_Q_cluster_hits += _hits[p]->Integral();
            }

            std::vector<art::Ptr<recob::Hit>>::iterator itr = hits_vec.begin();

            while (itr != hits_vec.end()) {
              diff_vec.clear();

              hit_energy = _hits[itr - hits_vec.begin()]->Integral();

              std::vector<sim::TrackIDE> trackides = bt_serv->HitToTrackIDEs(clockData, *itr);
              std::vector<sim::TrackIDE> eveides = bt_serv->HitToEveTrackIDEs(clockData, *itr);

              std::vector<sim::TrackIDE>::iterator idesitr = trackides.begin();

              while (idesitr != trackides.end()) {

                vec_trackid.push_back((*idesitr).trackID);

                const simb::MCParticle* particle = _particleList.at((*idesitr).trackID);
                int pdg = particle->PdgCode();
                if (pdg == 13 || pdg == -13) {
                  _hit_13++;
                  _en_13 += hit_energy * ((*idesitr).energyFrac);
                }
                idesitr++;
              }

              itr++;

            } // loop thru hits

            it = find(vec_pdg.begin(), vec_pdg.end(), 13);
            if (it != vec_pdg.end()) { no_cl_for_muon++; }

            it2 = find(vec_pdg.begin(), vec_pdg.end(), 11);
            if (it2 != vec_pdg.end()) { no_cl_for_electron++; }

            it3 = find(vec_pdg.begin(), vec_pdg.end(), -11);
            if (it3 != vec_pdg.end()) { no_cl_for_positron++; }

            it4 = find(vec_pdg.begin(), vec_pdg.end(), 111);
            if (it4 != vec_pdg.end()) { no_cl_for_pion_111++; }
            it6 = find(vec_pdg.begin(), vec_pdg.end(), 211);
            if (it6 != vec_pdg.end()) { no_cl_for_pion_211++; }
            it7 = find(vec_pdg.begin(), vec_pdg.end(), -211);
            if (it7 != vec_pdg.end()) { no_cl_for_pion_m211++; }
            it8 = find(vec_pdg.begin(), vec_pdg.end(), 2212);
            if (it8 != vec_pdg.end()) { no_cl_for_proton++; }

            sort(vec_pdg.begin(), vec_pdg.end());
            vec_pdg.erase(unique(vec_pdg.begin(), vec_pdg.end()), vec_pdg.end());

            //same for vec_trackid:

            sort(vec_trackid.begin(), vec_trackid.end());
            vec_trackid.erase(unique(vec_trackid.begin(), vec_trackid.end()), vec_trackid.end());
            for (unsigned int ii = 0; ii < vec_trackid.size(); ++ii) {
              all_trackids.push_back(vec_trackid[ii]);
            }

            //...............................................................
            //Also Make vec_trackid_mother unique:

            sort(vec_trackid_mother.begin(), vec_trackid_mother.end());
            vec_trackid_mother.erase(unique(vec_trackid_mother.begin(), vec_trackid_mother.end()),
                                     vec_trackid_mother.end());

            // Also Make vec_trackid_mother_en unique:

            sort(vec_trackid_mother_en.begin(), vec_trackid_mother_en.end());
            vec_trackid_mother_en.erase(
              unique(vec_trackid_mother_en.begin(), vec_trackid_mother_en.end()),
              vec_trackid_mother_en.end());

            //........................................................................

            // Q: How many clusters it takes to contain a certain particle?
            for (unsigned int ii = 0; ii < mc_trackids.size(); ++ii) {
              it5 = find(vec_trackid.begin(), vec_trackid.end(), mc_trackids[ii]);
              if (it5 != vec_trackid.end()) {
                ids.push_back(mc_trackids[ii]); // then make it unique
                noCluster++;
              }
            }

            fNoParticles_pdg->Fill(vec_pdg.size());

            fNoParticles_trackid->Fill(vec_trackid.size());

            fNoParticles_trackid_mother->Fill(vec_trackid_mother.size());

            no_of_clusters++;

            no_of_particles_in_cluster += vec_pdg.size();
            sum_vec_trackid += vec_trackid_mother.size();
            total_no_hits_in_clusters += _hits.size();

            vec_pdg.clear();

            vec_trackid.clear();

            vec_trackid_mother.clear();

            vec_trackid_mother_en.clear();

          } // end if cluster is in correct view

          hits_vec.clear();
        } // for each cluster

        sort(all_trackids.begin(), all_trackids.end());
        all_trackids.erase(unique(all_trackids.begin(), all_trackids.end()), all_trackids.end());

        // now I have a vector(all_trackids) that only contains unique trackids
        // go to it and search for each trackid in every cluster

        sort(ids.begin(), ids.end());
        ids.erase(unique(ids.begin(), ids.end()), ids.end());
        double no_of_clusters_per_track = noCluster / ids.size();

        // Filling Histograms:

        if (no_cl_for_muon != 0) { fCl_for_Muon->Fill(no_cl_for_muon); }

        fno_of_clusters_per_track->Fill(no_of_clusters_per_track);

        no_cl_for_muon = 0;
        no_cl_for_electron = 0;
        no_cl_for_positron = 0;
        no_cl_for_pion_111 = 0;
        no_cl_for_pion_211 = 0;
        no_cl_for_pion_m211 = 0;
        no_cl_for_proton = 0;
        noCluster = 0;
        ids.clear();

      } // for each plane
      double result = no_of_particles_in_cluster / no_of_clusters;

      double no_noise_hits = hits.size() - total_no_hits_in_clusters;
      double percent_noise = double(no_noise_hits / hits.size()) * 100;

      double no_trackid_per_cl_per_event = sum_vec_trackid / no_of_clusters;
      fNoParticles_trackid_per_event->Fill(no_trackid_per_cl_per_event);
      fNoParticles_pdg_per_event->Fill(result);
      fNoClustersInEvent->Fill(clusters.size());

      fPercentNoise->Fill(percent_noise);

      //-----------------------------------------------------------------------
      for (unsigned int i = 0; i < mclist.size(); ++i) {
        art::Ptr<simb::MCTruth> mc(mclist[i]);

        for (int ii = 0; ii < mc->NParticles(); ++ii) {
          simb::MCParticle part(mc->GetParticle(ii));
          std::cout << "FROM MC TRUTH,the particle's pdg code is: " << part.PdgCode() << std::endl;
          std::cout << "with energy= " << part.E();
          if (abs(part.PdgCode()) == 13) { std::cout << " I have a muon!!!" << std::endl; }
          if (abs(part.PdgCode()) == 111) { std::cout << " I have a pi zero!!!" << std::endl; }
        }
      }
    }

    // Now I can load just hits (not clusters) and see how many of what kind are
    // lost due to dbscan marking them as noise:

    // Load hits and copy most of the code from above:

    double hit_13 = 0, hit_11 = 0, hit_m_11 = 0, hit_111 = 0, hit_211 = 0, hit_m211 = 0,
           hit_2212 = 0, hit_2112 = 0;
    double en_13 = 0, en_11 = 0, en_m11 = 0, en_111 = 0, en_211 = 0, en_m211 = 0, en_2212 = 0,
           en_2112 = 0;

    std::vector<art::Ptr<recob::Hit>>::iterator itr = hits.begin();
    while (itr != hits.end()) {

      std::vector<sim::TrackIDE> trackides = bt_serv->HitToTrackIDEs(clockData, *itr);
      std::vector<sim::TrackIDE> eveides = bt_serv->HitToEveTrackIDEs(clockData, *itr);

      std::vector<sim::TrackIDE>::iterator idesitr = trackides.begin();

      hit_energy = hits[itr - hits.begin()]->Integral();

      while (idesitr != trackides.end()) {

        const simb::MCParticle* particle = _particleList.at((*idesitr).trackID);

        int pdg = particle->PdgCode();
        diff_vec.clear();

        if (pdg == 13 || pdg == -13) {
          hit_13++;
          en_13 += hit_energy * ((*idesitr).energyFrac);
        }
        idesitr++;

      } // trackIDs

      itr++;
    } // hits

    if (hit_13 != 0) { fPercent_lost_muon_hits->Fill(100 - ((_hit_13 / hit_13) * 100)); }
    if (hit_11 != 0) { fPercent_lost_electron_hits->Fill(100 - ((_hit_11 / hit_11) * 100)); }
    if (hit_m_11 != 0) { fPercent_lost_positron_hits->Fill(100 - ((_hit_m_11 / hit_m_11) * 100)); }
    if (hit_111 != 0) { fPercent_lost_111_hits->Fill(100 - ((_hit_111 / hit_111) * 100)); }

    if (hit_211 != 0) { fPercent_lost_211_hits->Fill(100 - ((_hit_211 / hit_211) * 100)); }

    if (hit_m211 != 0) { fPercent_lost_m211_hits->Fill(100 - ((_hit_m211 / hit_m211) * 100)); }
    if (hit_2212 != 0) { fPercent_lost_2212_hits->Fill(100 - ((_hit_2212 / hit_2212) * 100)); }
    if (hit_2112 != 0) { fPercent_lost_2112_hits->Fill(100 - ((_hit_2112 / hit_2112) * 100)); }

    std::cout << "****** mu E from clusters = " << _en_13 << std::endl;
    std::cout << "****** mu E from hits = " << en_13 << std::endl;

    if (en_13 != 0) { fPercent_lost_muon_energy->Fill(100 - ((_en_13 / en_13) * 100)); }
    if (en_11 != 0) { fPercent_lost_electron_energy->Fill(100 - ((_en_11 / en_11) * 100)); }
    if (en_m11 != 0) { fPercent_lost_positron_energy->Fill(100 - ((_en_m11 / en_m11) * 100)); }
    if (en_111 != 0) { fPercent_lost_111_energy->Fill(100 - ((_en_111 / en_111) * 100)); }
    if (en_211 != 0) { fPercent_lost_211_energy->Fill(100 - ((_en_211 / en_211) * 100)); }
    if (en_m211 != 0) { fPercent_lost_m211_energy->Fill(100 - ((_en_m211 / en_m211) * 100)); }
    if (en_2212 != 0) { fPercent_lost_2212_energy->Fill(100 - ((_en_2212 / en_2212) * 100)); }
    if (en_2112 != 0) { fPercent_lost_2112_energy->Fill(100 - ((_en_2112 / en_2112) * 100)); }
  }

} // end namespace

DEFINE_ART_MODULE(cluster::DBclusterAna)

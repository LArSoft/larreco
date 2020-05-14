////////////////////////////////////////////////////////////////////////
//
//  ClusterMatchAlg source
//
////////////////////////////////////////////////////////////////////////

#include "ClusterMatchAlg.h"

#include "TString.h"
#include "TTree.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace cluster {

  //##################################################################
  ClusterMatchAlg::ClusterMatchAlg(fhicl::ParameterSet const& pset) : _ModName_MCTruth("")
  //##################################################################
  {

    _debug_mode = pset.get<bool>("DebugMode");
    _store_sps = pset.get<bool>("StoreSpacePoint");
    _num_sps_cut = pset.get<size_t>("CutParam_NumSpacePoint");
    _overlay_tratio_cut = pset.get<double>("CutParam_OverlayTimeFraction");
    _qratio_cut = pset.get<double>("CutParam_SumChargeRatio");
    std::vector<size_t> algo_list = pset.get<std::vector<size_t>>("MatchAlgoList");

    _sps_algo = new trkf::SpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
    _match_tree = 0;
    _cluster_tree = 0;
    _save_cluster_info = true;
    _det_params_prepared = false;

    for (size_t i = 0; i < (size_t)(kMATCH_METHOD_MAX); ++i)

      _match_methods[i] = false;

    for (auto const v : algo_list) {

      if (v >= (size_t)(kMATCH_METHOD_MAX))

        mf::LogError("ClusterMatchAlg") << Form("Invalid algorithm enum: %zu", v);

      else
        _match_methods[v] = true;
    }

    ReportConfig();

    ClearEventInfo();
  }

  //########################################
  void
  ClusterMatchAlg::ReportConfig() const
  //########################################
  {
    std::ostringstream msg;
    msg << std::endl
        << " ClusterMatchAlg Configuration:              " << std::endl
        << "---------------------------------------------" << std::endl;
    msg << " Debug Mode ... " << (_debug_mode ? "enabled!" : "disabled!") << std::endl;
    msg << " RoughZ ....... " << (_match_methods[kRoughZ] ? "enabled!" : "disabled!") << std::endl;
    msg << " RoughT ....... " << (_match_methods[kRoughT] ? "enabled!" : "disabled!") << std::endl;
    msg << " SpacePoint ... " << (_match_methods[kSpacePoint] ? "enabled!" : "disabled!")
        << std::endl;
    msg << " SumCharge .... " << (_match_methods[kSumCharge] ? "enabled!" : "disabled!")
        << std::endl;
    msg << std::endl;
    msg << " Overlay-Time Fraction Cut : " << _overlay_tratio_cut << std::endl
        << " Charge-Ratio Diff. Cut    : " << _qratio_cut << std::endl
        << " Minimum # of SpacePoint   : " << _num_sps_cut << std::endl
        << std::endl;
    msg << "---------------------------------------------" << std::endl;

    mf::LogWarning("ClusterMatchAlg") << msg.str();
  }

  void
  ClusterMatchAlg::ClearMatchInputInfo()
  {
    _ucluster_v.clear();
    _vcluster_v.clear();
    _wcluster_v.clear();

    _uhits_v.clear();
    _vhits_v.clear();
    _whits_v.clear();
  }

  void
  ClusterMatchAlg::ClearMatchOutputInfo()
  {
    _matched_uclusters_v.clear();
    _matched_vclusters_v.clear();
    _matched_wclusters_v.clear();
    _matched_sps_v.clear();
  }

  void
  ClusterMatchAlg::ClearTTreeInfo()
  {

    _mc_E = 0;
    _mc_Px = 0;
    _mc_Py = 0;
    _mc_Pz = 0;
    _mc_Vx = 0;
    _mc_Vy = 0;
    _mc_Vz = 0;
    _pdgid = 0;
    _tot_u = 0;
    _tot_v = 0;
    _tot_w = 0;
    _tot_pass_z = 0;
    _tot_pass_t = 0;
    _tot_pass_sps = 0;
    _tot_pass_qsum = 0;
    _qratio_v.clear();
    _uv_tratio_v.clear();
    _vw_tratio_v.clear();
    _wu_tratio_v.clear();
    _u_nhits_v.clear();
    _v_nhits_v.clear();
    _w_nhits_v.clear();
    _nsps.clear();

    _view_v.clear();
    _charge_v.clear();
    _nhits_v.clear();
    _tstart_min_v.clear();
    _tstart_max_v.clear();
    _tpeak_min_v.clear();
    _tpeak_max_v.clear();
    _tend_min_v.clear();
    _tend_max_v.clear();
  }

  //##################################################################
  void
  ClusterMatchAlg::ClearEventInfo()
  //##################################################################
  {
    // Clear input event data holders
    ClearMatchInputInfo();

    // Clear result data holders
    ClearMatchOutputInfo();

    /// Clear TTree variables
    ClearTTreeInfo();
  }

  //##################################################################
  void
  ClusterMatchAlg::PrepareTTree()
  //##################################################################
  {
    if (!_match_tree) {
      art::ServiceHandle<art::TFileService const> fileService;
      _match_tree = fileService->make<TTree>("match_tree", "");
      _match_tree->Branch("mc_E", &_mc_E, "mc_E/D");
      _match_tree->Branch("mc_Px", &_mc_Px, "mc_Px/D");
      _match_tree->Branch("mc_Py", &_mc_Py, "mc_Py/D");
      _match_tree->Branch("mc_Pz", &_mc_Pz, "mc_Pz/D");
      _match_tree->Branch("mc_Vx", &_mc_Vx, "mc_Vx/D");
      _match_tree->Branch("mc_Vy", &_mc_Vy, "mc_Vy/D");
      _match_tree->Branch("mc_Vz", &_mc_Vz, "mc_Vz/D");

      _match_tree->Branch("pdgid", &_pdgid, "pdgid/I");
      _match_tree->Branch("tot_u", &_tot_u, "tot_u/s");
      _match_tree->Branch("tot_v", &_tot_v, "tot_v/s");
      _match_tree->Branch("tot_w", &_tot_w, "tot_w/s");
      _match_tree->Branch("tot_pass_t", &_tot_pass_t, "tot_pass_t/s");
      _match_tree->Branch("tot_pass_z", &_tot_pass_z, "tot_pass_z/s");
      _match_tree->Branch("tot_pass_sps", &_tot_pass_sps, "tot_pass_sps/s");
      _match_tree->Branch("tot_pass_qsum", &_tot_pass_qsum, "tot_pass_qsum/s");

      _match_tree->Branch("uv_tratio_v", "std::vector<double>", &_uv_tratio_v);
      _match_tree->Branch("vw_tratio_v", "std::vector<double>", &_vw_tratio_v);
      _match_tree->Branch("wu_tratio_v", "std::vector<double>", &_wu_tratio_v);

      _match_tree->Branch("qratio_v", "std::vector<double>", &_qratio_v);
      _match_tree->Branch("u_nhits_v", "std::vector<UShort_t>", &_u_nhits_v);
      _match_tree->Branch("v_nhits_v", "std::vector<UShort_t>", &_v_nhits_v);
      _match_tree->Branch("w_nhits_v", "std::vector<UShort_t>", &_w_nhits_v);
      _match_tree->Branch("nsps", "std::vector<UShort_t>", &_nsps);
    }
    if (_save_cluster_info && !_cluster_tree) {
      art::ServiceHandle<art::TFileService const> fileService;
      _cluster_tree = fileService->make<TTree>("cluster_tree", "");
      _cluster_tree->Branch("view_v", "std::vector<uint16_t>", &_view_v);
      _cluster_tree->Branch("charge_v", "std::vector<double>", &_charge_v);
      _cluster_tree->Branch("nhits_v", "std::vector<uint16_t>", &_nhits_v);
      _cluster_tree->Branch("tstart_min_v", "std::vector<double>", &_tstart_min_v);
      _cluster_tree->Branch("tstart_max_v", "std::vector<double>", &_tstart_max_v);
      _cluster_tree->Branch("tpeak_min_v", "std::vector<double>", &_tpeak_min_v);
      _cluster_tree->Branch("tpeak_max_v", "std::vector<double>", &_tpeak_max_v);
      _cluster_tree->Branch("tend_min_v", "std::vector<double>", &_tend_min_v);
      _cluster_tree->Branch("tend_max_v", "std::vector<double>", &_tend_max_v);
    }
  }

  //##########################################################################################
  void
  ClusterMatchAlg::FillMCInfo(const art::Event& evt)
  //##########################################################################################
  {
    if (!_ModName_MCTruth.size()) return;

    std::vector<const simb::MCTruth*> mciArray;

    try {

      evt.getView(_ModName_MCTruth, mciArray);
    }
    catch (art::Exception const& e) {

      if (e.categoryCode() != art::errors::ProductNotFound) throw;
    }

    for (size_t i = 0; i < mciArray.size(); ++i) {

      if (i == 1) {
        mf::LogWarning("ClusterMatchAlg") << " Ignoring > 2nd MCTruth in MC generator...";
        break;
      }
      const simb::MCTruth* mci_ptr(mciArray.at(i));

      for (size_t j = 0; j < (size_t)(mci_ptr->NParticles()); ++j) {

        if (j == 1) {
          mf::LogWarning("ClusterMatchAlg") << " Ignoring > 2nd MCParticle in MC generator...";
          break;
        }

        const simb::MCParticle part(mci_ptr->GetParticle(j));

        _pdgid = part.PdgCode();
        _mc_E = part.E();
        _mc_Px = part.Px();
        _mc_Py = part.Py();
        _mc_Pz = part.Pz();
        _mc_Vx = part.Vx();
        _mc_Vy = part.Vy();
        _mc_Vz = part.Vz();
      }
    }
  }

  void
  ClusterMatchAlg::PrepareDetParams()
  {
    if (!_det_params_prepared) {
      // Total number of planes
      art::ServiceHandle<geo::Geometry const> geo;
      _tot_planes = geo->Nplanes();

      // Ask DetectorPrperties about time-offset among different wire planes ... used to correct timing
      // difference among different wire planes in the following loop.
      const detinfo::DetectorProperties* det_h =
        lar::providerFrom<detinfo::DetectorPropertiesService>();
      _time_offset_uplane = det_h->GetXTicksOffset(geo::kU, 0, 0);
      _time_offset_vplane = det_h->GetXTicksOffset(geo::kV, 0, 0);
      _time_offset_wplane = 0;
      if (_tot_planes > 2) _time_offset_wplane = det_h->GetXTicksOffset(geo::kW, 0, 0);
      _det_params_prepared = true;
    }
  }

  void
  ClusterMatchAlg::AppendClusterInfo(const recob::Cluster& in_cluster,
                                     const std::vector<art::Ptr<recob::Hit>>& in_hit_v)
  {

    PrepareDetParams();
    cluster_match_info ci(in_cluster.ID());
    ci.view = in_cluster.View();

    art::PtrVector<recob::Hit> hit_ptrv;
    FillHitInfo(ci, hit_ptrv, in_hit_v);

    // Save created art::PtrVector & cluster_match_info struct object
    switch (ci.view) {
    case geo::kU:
      _uhits_v.push_back(hit_ptrv);
      _ucluster_v.push_back(ci);
      AppendClusterTreeVariables(ci);
      break;
    case geo::kV:
      _vhits_v.push_back(hit_ptrv);
      _vcluster_v.push_back(ci);
      AppendClusterTreeVariables(ci);
      break;
    case geo::kW:
      _whits_v.push_back(hit_ptrv);
      _wcluster_v.push_back(ci);
      AppendClusterTreeVariables(ci);
      break;
    default:
      mf::LogError("ClusterMatchAlg") << Form("Found an invalid plane ID: %d", in_cluster.View());
    }
  }

  void
  ClusterMatchAlg::AppendClusterInfo(const art::Ptr<recob::Cluster> in_cluster,
                                     const std::vector<art::Ptr<recob::Hit>>& in_hit_v)
  {

    PrepareDetParams();
    cluster_match_info ci(in_cluster->ID());
    ci.view = in_cluster->View();

    art::PtrVector<recob::Hit> hit_ptrv;
    FillHitInfo(ci, hit_ptrv, in_hit_v);

    // Save created art::PtrVector & cluster_match_info struct object
    switch (ci.view) {
    case geo::kU:
      _uhits_v.push_back(hit_ptrv);
      _ucluster_v.push_back(ci);
      AppendClusterTreeVariables(ci);
      break;
    case geo::kV:
      _vhits_v.push_back(hit_ptrv);
      _vcluster_v.push_back(ci);
      AppendClusterTreeVariables(ci);
      break;
    case geo::kW:
      _whits_v.push_back(hit_ptrv);
      _wcluster_v.push_back(ci);
      AppendClusterTreeVariables(ci);
      break;
    default:
      mf::LogError("ClusterMatchAlg") << Form("Found an invalid plane ID: %d", in_cluster->View());
    }
  }

  void
  ClusterMatchAlg::FillHitInfo(cluster_match_info& ci,
                               art::PtrVector<recob::Hit>& out_hit_v,
                               const std::vector<art::Ptr<recob::Hit>>& in_hit_v)
  {

    out_hit_v.reserve(in_hit_v.size());

    double time_offset = 0;
    if (ci.view == geo::kU)
      time_offset = _time_offset_uplane;
    else if (ci.view == geo::kV)
      time_offset = _time_offset_vplane;
    else if (ci.view == geo::kW)
      time_offset = _time_offset_wplane;

    // Loop over hits in this cluster
    for (auto const hit : in_hit_v) {

      unsigned int wire = hit->WireID().Wire;
      double tstart = hit->PeakTimePlusRMS(-1.) - time_offset;
      double tpeak = hit->PeakTime() - time_offset;
      double tend = hit->PeakTimePlusRMS(+1.) - time_offset;

      ci.sum_charge += hit->Integral();

      ci.wire_max = (ci.wire_max < wire) ? wire : ci.wire_max;
      ci.wire_min = (ci.wire_min > wire) ? wire : ci.wire_min;

      ci.start_time_max = (ci.start_time_max < tstart) ? tstart : ci.start_time_max;
      ci.peak_time_max = (ci.peak_time_max < tpeak) ? tpeak : ci.peak_time_max;
      ci.end_time_max = (ci.end_time_max < tend) ? tend : ci.end_time_max;

      ci.start_time_min = (ci.start_time_min > tstart) ? tstart : ci.start_time_min;
      ci.peak_time_min = (ci.peak_time_min > tpeak) ? tpeak : ci.peak_time_min;
      ci.end_time_min = (ci.end_time_min > tend) ? tend : ci.end_time_min;

      out_hit_v.push_back(hit);
    }

    ci.nhits = in_hit_v.size();
  }

  void
  ClusterMatchAlg::AppendClusterTreeVariables(const cluster_match_info& ci)
  {

    if (_cluster_tree) {
      _view_v.push_back(ci.view);
      _charge_v.push_back(ci.sum_charge);
      _nhits_v.push_back(ci.nhits);
      _tstart_min_v.push_back(ci.start_time_min);
      _tstart_max_v.push_back(ci.start_time_max);
      _tpeak_min_v.push_back(ci.peak_time_min);
      _tpeak_max_v.push_back(ci.peak_time_max);
      _tend_min_v.push_back(ci.end_time_min);
      _tend_max_v.push_back(ci.end_time_max);
    }
  }

  //########################################################################################
  bool
  ClusterMatchAlg::Match_RoughZ(const cluster_match_info& ci1,
                                const cluster_match_info& ci2,
                                const geo::View_t v1,
                                const geo::View_t v2) const
  //########################################################################################
  {
    art::ServiceHandle<geo::Geometry const> geo_h;
    double y, z_min, z_max;
    y = z_min = z_max = -1;
    geo_h->IntersectionPoint(ci1.wire_min, ci2.wire_min, v1, v2, 0, 0, y, z_min);
    geo_h->IntersectionPoint(ci1.wire_max, ci2.wire_max, v1, v2, 0, 0, y, z_max);
    return (z_max > z_min);
  }

  //###########################################################################################
  bool
  ClusterMatchAlg::Match_RoughTime(const cluster_match_info& ci1, const cluster_match_info& ci2)
  //###########################################################################################
  {
    //return (!(ci1.end_time_max < ci2.start_time_min || ci2.end_time_max < ci1.start_time_min));
    double time_overlay = std::min(ci1.end_time_max, ci2.end_time_max) -
                          std::max(ci1.start_time_min, ci2.start_time_min);

    //if(time_overlay <= 0 && !_debug_mode) return false;

    double overlay_tratio =
      time_overlay /
      (ci1.end_time_max - ci1.start_time_min + ci2.end_time_max - ci2.start_time_min) * 2.;

    if ((ci1.view == geo::kU && ci2.view == geo::kV) ||
        (ci1.view == geo::kV && ci2.view == geo::kU))
      _uv_tratio_v.push_back(overlay_tratio);
    else if ((ci1.view == geo::kV && ci2.view == geo::kW) ||
             (ci1.view == geo::kW && ci2.view == geo::kV))
      _vw_tratio_v.push_back(overlay_tratio);
    else if ((ci1.view == geo::kW && ci2.view == geo::kU) ||
             (ci1.view == geo::kU && ci2.view == geo::kW))
      _wu_tratio_v.push_back(overlay_tratio);

    return (overlay_tratio > _overlay_tratio_cut);
  }

  //##############################################################################################
  bool
  ClusterMatchAlg::Match_SumCharge(const cluster_match_info& uc, const cluster_match_info& vc)
  //##############################################################################################
  {
    double qratio = (uc.sum_charge) / (vc.sum_charge);

    // For quality check log
    _qratio_v.push_back(qratio);

    return ((1 - _qratio_cut) < qratio && (qratio) < (1 + _qratio_cut));
  }

  //#####################################################################################################
  bool
  ClusterMatchAlg::Match_SpacePoint(const size_t uindex,
                                    const size_t vindex,
                                    const size_t windex,
                                    std::vector<recob::SpacePoint>& sps_v)
  //#####################################################################################################
  {
    bool use_wplane = _tot_planes > 2;

    if (uindex >= _ucluster_v.size() || vindex >= _vcluster_v.size() ||
        (use_wplane && (windex >= _wcluster_v.size()))) {

      mf::LogError("ClusterMatchAlg")
        << std::endl
        << Form(
             "Requested to cluster-index (U,V,W) = (%zu,%zu,%zu) where max-length is (%zu,%zu,%zu)",
             uindex,
             vindex,
             windex,
             _ucluster_v.size(),
             _vcluster_v.size(),
             _wcluster_v.size())
        << std::endl;
      return false;
    }

    // Define a time range in which hits are used for spacepoint finding ... here "peak time" is the relevant one
    double trange_min =
      std::min(_ucluster_v.at(uindex).peak_time_min, _vcluster_v.at(vindex).peak_time_min);
    if (use_wplane) trange_min = std::min(trange_min, _wcluster_v.at(windex).peak_time_min);

    double trange_max =
      std::max(_ucluster_v.at(uindex).peak_time_max, _vcluster_v.at(vindex).peak_time_max);
    if (use_wplane) trange_max = std::max(trange_max, _wcluster_v.at(windex).peak_time_max);

    // Space-point algorithm applies additional dT
    trange_min -= _sps_algo->maxDT();
    trange_max += _sps_algo->maxDT();

    // Make PtrVector<recob::Hit> for relevant Hits
    art::PtrVector<recob::Hit> hit_group;
    size_t max_size = _uhits_v.at(uindex).size() + _vhits_v.at(vindex).size();
    if (use_wplane) max_size += _whits_v.at(windex).size();
    hit_group.reserve(max_size);
    // Loop over hits in U-plane
    for (auto const hit : _uhits_v.at(uindex)) {
      if (hit->PeakTime() < trange_min) continue;
      if (hit->PeakTime() > trange_max) continue;
      hit_group.push_back(hit);
    }
    // Check if any hit found in this plane
    size_t u_nhits = hit_group.size();
    if (!u_nhits && !_debug_mode) return false;
    // Loop over hits in V-plane
    for (auto const hit : _vhits_v.at(vindex)) {
      if (hit->PeakTime() < trange_min) continue;
      if (hit->PeakTime() > trange_max) continue;
      hit_group.push_back(hit);
    }
    // Check if any hit found in this plane
    size_t v_nhits = hit_group.size() - u_nhits;
    if (!(v_nhits) && !_debug_mode) return false;

    // Loop over hits in W-plane
    if (use_wplane) {
      for (auto const hit : _whits_v.at(windex)) {
        if (hit->PeakTime() < trange_min) continue;
        if (hit->PeakTime() > trange_max) continue;
        hit_group.push_back(hit);
      }
    }
    // Check if any hit found in this plane
    size_t w_nhits = hit_group.size() - u_nhits - v_nhits;
    if (!(w_nhits) && use_wplane && !_debug_mode) return false;

    // Run SpacePoint finder algo
    if (u_nhits && v_nhits && (!use_wplane || (w_nhits && use_wplane))) {
      _sps_algo->clearHitMap();
      _sps_algo->makeSpacePoints(hit_group, sps_v);
    }

    size_t nsps = sps_v.size();
    _u_nhits_v.push_back(u_nhits);
    _v_nhits_v.push_back(v_nhits);
    if (use_wplane) _w_nhits_v.push_back(w_nhits);
    _nsps.push_back(nsps);

    if (nsps < _num_sps_cut) return false;
    return true;
  }

  //#################################################################################
  std::vector<std::vector<unsigned int>>
  ClusterMatchAlg::GetMatchedClusters() const
  //#################################################################################
  {
    std::vector<std::vector<unsigned int>> result;
    result.push_back(_matched_uclusters_v);
    result.push_back(_matched_vclusters_v);
    result.push_back(_matched_wclusters_v);
    return result;
  }

  //#######################################################################
  void
  ClusterMatchAlg::MatchTwoPlanes()
  //#######################################################################
  {
    std::ostringstream msg;
    msg << Form("Received (U,V,W) = (%zu,%zu,%zu) clusters...",
                _uhits_v.size(),
                _vhits_v.size(),
                _whits_v.size())
        << std::endl;
    _tot_u = _ucluster_v.size();
    _tot_v = _vcluster_v.size();
    _tot_w = _wcluster_v.size();

    if (!(_tot_u + _tot_v + _tot_w)) {

      mf::LogError(__PRETTY_FUNCTION__)
        << "No input cluster info found! Aborting the function call...";

      return;
    }

    // Initialization
    PrepareTTree();
    ClearMatchOutputInfo();

    bool overlay_2d = false;
    for (size_t uci_index = 0; uci_index < _ucluster_v.size(); ++uci_index) {

      for (size_t vci_index = 0; vci_index < _vcluster_v.size(); ++vci_index) {

        overlay_2d = true;

        // Apply cuts
        // Rough z-position overlay cut
        if (_match_methods[kRoughZ]) {

          if (Match_RoughZ(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index), geo::kU, geo::kV))
            _tot_pass_z++;
          else if (!_debug_mode)
            continue;
          else
            overlay_2d = false;
        }

        // Sum charge cut
        if (_match_methods[kSumCharge]) {

          if (Match_SumCharge(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index)))
            _tot_pass_qsum++;
          else if (!_debug_mode)
            continue;
          else
            overlay_2d = false;
        }

        // Rough time overlap cut
        if (_match_methods[kRoughT]) {

          if (Match_RoughTime(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index)))
            _tot_pass_t++;
          else if (!_debug_mode)
            continue;
          else
            overlay_2d = false;
        }

        // SpacePoint cut
        std::vector<recob::SpacePoint> sps_v;
        if (_match_methods[kSpacePoint]) {

          if (Match_SpacePoint(uci_index, vci_index, 0, sps_v))
            _tot_pass_sps++;
          else if (!_debug_mode)
            continue;
          else
            overlay_2d = false;
        }

        if (overlay_2d) {
          _matched_uclusters_v.push_back((unsigned int)(_ucluster_v[uci_index].cluster_index));
          _matched_vclusters_v.push_back((unsigned int)(_vcluster_v[vci_index].cluster_index));
          if (_store_sps) _matched_sps_v.push_back(sps_v);
        }
      } // end of ... _vcluster_v loop
    }   // end of ... _ucluster_v loop

    // Report
    msg << std::endl
        << Form("Found %zu matched cluster pairs...", _matched_uclusters_v.size()) << std::endl;
    for (size_t i = 0; i < _matched_uclusters_v.size(); ++i) {

      if (i == 0) msg << "Listing matched clusters (U,V)..." << std::endl;

      msg << Form("Pair %-2zu: (%-3d, %-3d)", i, _matched_uclusters_v[i], _matched_vclusters_v[i])
          << std::endl;
    }
    msg << std::endl;
    mf::LogWarning("ClusterMatchAlg") << msg.str();

    if (_match_tree) _match_tree->Fill();
    if (_cluster_tree) _cluster_tree->Fill();

    // Clear input event data holders
    ClearMatchInputInfo();
    /// Clear TTree variables
    ClearTTreeInfo();
  }

  //#######################################################################
  void
  ClusterMatchAlg::MatchThreePlanes()
  //#######################################################################
  {
    std::ostringstream msg;
    msg << Form("Received (U,V,W) = (%zu,%zu,%zu) clusters...",
                _uhits_v.size(),
                _vhits_v.size(),
                _whits_v.size())
        << std::endl;
    _tot_u = _ucluster_v.size();
    _tot_v = _vcluster_v.size();
    _tot_w = _wcluster_v.size();

    if (!(_tot_u + _tot_v + _tot_w)) {

      mf::LogError(__PRETTY_FUNCTION__)
        << "No input cluster info found! Aborting the function call...";

      return;
    }

    // Clear match information
    PrepareTTree();
    ClearMatchOutputInfo();

    bool overlay_2d = true;
    bool overlay_3d = true;
    // Loop over all possible u-v-w cluster combination
    for (size_t uci_index = 0; uci_index < _ucluster_v.size(); ++uci_index) {

      for (size_t vci_index = 0; vci_index < _vcluster_v.size(); ++vci_index) {

        // Apply cuts that can be done with U&V planes here
        overlay_2d = true;

        // Rough z-position overlay cut
        if (_match_methods[kRoughZ]) {

          if (Match_RoughZ(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index), geo::kU, geo::kV))
            _tot_pass_z++;
          else if (!_debug_mode)
            continue;
          else
            overlay_2d = false;
        }

        // Sum charge cut
        if (_match_methods[kSumCharge]) {

          if (Match_SumCharge(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index)))
            _tot_pass_qsum++;
          else if (!_debug_mode)
            continue;
          else
            overlay_2d = false;
        }

        for (size_t wci_index = 0; wci_index < _wcluster_v.size(); ++wci_index) {

          overlay_3d = overlay_2d;
          // Apply cuts that requires 3 planes here

          // Rough time overlap cut
          if (_match_methods[kRoughT]) {

            bool rough_time_match =
              Match_RoughTime(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index));
            if (!_debug_mode && !rough_time_match) continue;

            rough_time_match =
              (Match_RoughTime(_vcluster_v.at(vci_index), _wcluster_v.at(wci_index)) &&
               rough_time_match);
            if (!_debug_mode && !rough_time_match) continue;

            rough_time_match =
              (Match_RoughTime(_wcluster_v.at(wci_index), _ucluster_v.at(uci_index)) &&
               rough_time_match);

            overlay_3d = overlay_3d && rough_time_match;
            if (rough_time_match)
              _tot_pass_t++;
            else if (!_debug_mode)
              continue;
          }

          // SpacePoint cut
          std::vector<recob::SpacePoint> sps_v;
          if (_match_methods[kSpacePoint]) {

            if (Match_SpacePoint(uci_index, vci_index, wci_index, sps_v))
              _tot_pass_sps++;
            else if (!_debug_mode)
              continue;
            else
              overlay_3d = false;
          }

          if (overlay_3d) {
            _matched_uclusters_v.push_back((unsigned int)(_ucluster_v[uci_index].cluster_index));
            _matched_vclusters_v.push_back((unsigned int)(_vcluster_v[vci_index].cluster_index));
            _matched_wclusters_v.push_back((unsigned int)(_wcluster_v[wci_index].cluster_index));
            if (_store_sps) _matched_sps_v.push_back(sps_v);
          }
        } // end of ... _wcluster_v loop
      }   // end of ... _vcluster_v loop
    }     // end of ... _ucluster_v loop

    // Report
    msg << std::endl
        << Form("Found %zu matched cluster pairs...", _matched_uclusters_v.size()) << std::endl;
    for (size_t i = 0; i < _matched_uclusters_v.size(); ++i) {

      if (i == 0) msg << "Listing matched clusters (U,V,W)..." << std::endl;

      msg << Form("Pair %-2zu: (%-3d, %-3d, %-3d)",
                  i,
                  _matched_uclusters_v[i],
                  _matched_vclusters_v[i],
                  _matched_wclusters_v[i])
          << std::endl;
    }
    msg << std::endl;
    mf::LogWarning("ClusterMatchAlg") << msg.str();

    if (_match_tree) _match_tree->Fill();
    if (_cluster_tree) _cluster_tree->Fill();
    // Clear input event data holders
    ClearMatchInputInfo();
    /// Clear TTree variables
    ClearTTreeInfo();
  }

} // namespace match

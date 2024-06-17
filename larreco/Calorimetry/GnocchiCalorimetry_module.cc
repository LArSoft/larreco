// A gnew Calorimetry module
//
// Re-write of the Calorimetry module ported into art by E. Chruch
// from the original (ArgoNeuT) calorimetry module by Ornella Palamara and Maddalena Antonello.
//
// Gray Putnam
// grayputnam@uchicago.edu

#include <cmath>
#include <limits>  // std::numeric_limits<>
#include <numeric> // std::accumulate
#include <optional>
#include <string>

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/Calorimetry/INormalizeCharge.h"

#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// ROOT includes
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  constexpr unsigned int int_max_as_unsigned_int{std::numeric_limits<int>::max()};
  typedef std::pair<unsigned, std::vector<unsigned>> OrganizedHits;
}

namespace calo {

  class GnocchiCalorimetry : public art::EDProducer {
  public:
    struct Config {
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;

      enum ChargeMethod : unsigned {
        cmAmplitude = 0,
        cmIntegral = 1,
        cmSummedADC = 2,
        cmSummedIntegral = 3,
      };

      fhicl::Atom<std::string> TrackModuleLabel{Name("TrackModuleLabel"),
                                                Comment("Module label for track producer.")};

      fhicl::Atom<std::string> T0ModuleLabel{Name("T0ModuleLabel"),
                                             Comment("Module label for T0 time producer."),
                                             ""};
      fhicl::Atom<std::string> PFPModuleLabel{Name("PFPModuleLabel"),
                                             Comment("Module label for PFP producer. To be used to associate T0 with tracks."),
                                             ""};

      fhicl::Atom<std::string> AssocHitModuleLabel{
        Name("AssocHitModuleLabel"),
        Comment("Module label for association between tracks and hits. If not set, defaults to "
                "TrackModuleLabel."),
        ""};

      fhicl::Atom<unsigned> ChargeMethod{
        Name("ChargeMethod"),
        Comment("Method used to extract charge from a hit. Options: 0==Amplitude(), 1==Integral(), "
                "2==SummedADC(), 3==SummedIntegral(). See the ChargeMethod enum.")};

      fhicl::Atom<bool> FieldDistortion{
        Name("FieldDistortion"),
        Comment("True if field distortion (i.e. from space charge) is included in the input.")};

      fhicl::Atom<bool> FieldDistortionEfield{
        Name("FieldDistortionEfield"),
        Comment("True if field distortion (i.e. from space charge) is included in the input.")};

      fhicl::Atom<bool> TrackIsFieldDistortionCorrected{
        Name("TrackIsFieldDistortionCorrected"),
        Comment("Whether the space-points on the input tracks have their points corrected for the "
                "field distortions. "
                "I.e. whether the track trajectory points represent charge as seen by wires or the "
                "3D particle trajectory.")};

      fhicl::Atom<unsigned> Cryostat{Name("Cryostat"),
                                     Comment("Which cryostat number the input tracks occupy.")};

      fhicl::Atom<float> FieldDistortionCorrectionXSign{
        Name("FieldDistortionCorrectionXSign"),
        Comment("Sign of the field distortion correction to be applied in the X direction. "
                "Positive by default."),
        1.};

      fhicl::Table<calo::CalorimetryAlg::Config> CalorimetryAlgConfig{
        Name("CaloAlg"),
        Comment("Configuration for the calo::CalorimetryAlg")};

      fhicl::DelegatedParameter NormTools{
        Name("NormTools"),
        Comment("List of INormalizeCharge tool configurations to use.")};
    };

    using Parameters = art::EDProducer::Table<Config>;

    explicit GnocchiCalorimetry(Parameters const& config);

    void produce(art::Event& evt) override;

  private:
    Config fConfig;
    CalorimetryAlg fCaloAlg;
    std::vector<std::unique_ptr<INormalizeCharge>> fNormTools;

    // helper functions
    std::vector<std::vector<OrganizedHits>> OrganizeHits(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const std::vector<const recob::TrackHitMeta*>& thms,
      const recob::Track& track,
      unsigned nplanes);
    std::vector<std::vector<OrganizedHits>> OrganizeHitsIndividual(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const std::vector<const recob::TrackHitMeta*>& thms,
      const recob::Track& track,
      unsigned nplanes);
    std::vector<std::vector<OrganizedHits>> OrganizeHitsSnippets(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const std::vector<const recob::TrackHitMeta*>& thms,
      const recob::Track& track,
      unsigned nplanes);
    bool HitIsValid(const recob::TrackHitMeta* thm, const recob::Track& track);
    geo::Point_t GetLocation(const recob::Track& track,
                             const art::Ptr<recob::Hit> hit,
                             const recob::TrackHitMeta* meta);
    geo::Point_t GetLocationAtWires(const recob::Track& track,
                                    const art::Ptr<recob::Hit> hit,
                                    const recob::TrackHitMeta* meta);
    geo::Point_t WireToTrajectoryPosition(const geo::Point_t& loc, const geo::TPCID& tpc);
    geo::Point_t TrajectoryToWirePosition(const geo::Point_t& loc, const geo::TPCID& tpc);
    double GetPitch(const recob::Track& track,
                    const art::Ptr<recob::Hit> hit,
                    const recob::TrackHitMeta* meta);
    double GetCharge(const art::Ptr<recob::Hit> hit,
                     const std::vector<recob::Hit const*>& sharedHits);
    double GetEfield(const detinfo::DetectorPropertiesData& dprop,
                     const recob::Track& track,
                     const art::Ptr<recob::Hit> hit,
                     const recob::TrackHitMeta* meta);
    double Normalize(double dQdx,
                     const art::Event& e,
                     const recob::Hit& h,
                     const geo::Point_t& location,
                     const geo::Vector_t& direction,
                     double t0);
  };

} // end namespace calo

calo::GnocchiCalorimetry::GnocchiCalorimetry(Parameters const& param)
  : EDProducer{param}, fConfig(param()), fCaloAlg(fConfig.CalorimetryAlgConfig())
{
  produces<std::vector<anab::Calorimetry>>();
  produces<art::Assns<recob::Track, anab::Calorimetry>>();
  std::vector<fhicl::ParameterSet> norm_tool_configs =
    fConfig.NormTools.get<std::vector<fhicl::ParameterSet>>();
  for (const fhicl::ParameterSet& p : norm_tool_configs) {
    fNormTools.push_back(art::make_tool<INormalizeCharge>(p));
  }
}

void calo::GnocchiCalorimetry::produce(art::Event& evt)
{
  // Get services
  art::ServiceHandle<geo::Geometry const> geom;
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const det_prop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);

  size_t nplanes = geom->Nplanes();

  // Define output collections
  std::unique_ptr<std::vector<anab::Calorimetry>> outputCalo(new std::vector<anab::Calorimetry>);
  std::unique_ptr<art::Assns<recob::Track, anab::Calorimetry>> outputCaloAssn(
    new art::Assns<recob::Track, anab::Calorimetry>);

  // collect input
  art::Handle<std::vector<recob::Track>> trackListHandle;
  std::vector<art::Ptr<recob::Track>> tracklist;
  if (evt.getByLabel(fConfig.TrackModuleLabel(), trackListHandle)) {
    art::fill_ptr_vector(tracklist, trackListHandle);
  }

  // get the label to collect this hits
  const std::string& hitLabel = (fConfig.AssocHitModuleLabel() == "") ?
                                  fConfig.TrackModuleLabel() :
                                  fConfig.AssocHitModuleLabel();
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHits(trackListHandle, evt, hitLabel);


  // must be valid if the T0 module label is non-empty
  art::FindManyP<anab::T0> fmT0s(trackListHandle, evt, fConfig.T0ModuleLabel());
  art::FindManyP<recob::PFParticle> fmPFPs(trackListHandle, evt, fConfig.TrackModuleLabel());
  std::vector<art::Ptr<recob::PFParticle>> pfpList;
  if (fConfig.PFPModuleLabel().size()) {
    for (std::size_t itrk = 0, sz = fmPFPs.size(); itrk < sz; ++itrk) {
      // assuming every recob::Track has exactly one associated recob::PFParticle
      pfpList.push_back(fmPFPs.at(itrk).at(0));
    }
  }
  art::FindManyP<anab::T0> fmPFPT0s(pfpList, evt, fConfig.T0ModuleLabel());

  for (auto const& nt : fNormTools)
    nt->setup(evt);

  // iterate over all the tracks
  for (unsigned trk_i = 0; trk_i < tracklist.size(); trk_i++) {
    const recob::Track& track = *tracklist[trk_i];

    // collect input for this track
    const std::vector<art::Ptr<recob::Hit>>& hits = fmHits.at(trk_i);
    const std::vector<const recob::TrackHitMeta*>& thms = fmHits.data(trk_i);

    double T0 = 0;
    if (fConfig.T0ModuleLabel().size()) {
      if (fConfig.PFPModuleLabel().size()) {
	const std::vector<art::Ptr<anab::T0>>& this_t0s = fmPFPT0s.at(trk_i);
	if (this_t0s.size()) T0 = this_t0s.at(0)->Time();
      } else {
	const std::vector<art::Ptr<anab::T0>>& this_t0s = fmT0s.at(trk_i);
	if (this_t0s.size()) T0 = this_t0s.at(0)->Time();
      }
    }

    // organize the hits by plane
    std::vector<std::vector<OrganizedHits>> hit_indices = OrganizeHits(hits, thms, track, nplanes);

    for (unsigned plane_i = 0; plane_i < nplanes; plane_i++) {

      float kinetic_energy = 0.;
      std::vector<float> dEdxs;
      std::vector<float> dQdxs;
      std::vector<float> resranges;
      std::vector<float> deadwireresranges;
      float range = 0.;
      std::vector<float> pitches;
      std::vector<geo::Point_t> xyzs;
      std::vector<size_t> tp_indices;
      geo::PlaneID plane;

      // setup the plane ID
      plane.Plane = plane_i;
      plane.TPC = 0; // arbitrary -- tracks can cross TPC boundaries
      plane.Cryostat = fConfig.Cryostat();
      plane.isValid = true;

      std::vector<float> lengths;
      for (unsigned hit_i = 0; hit_i < hit_indices[plane_i].size(); hit_i++) {
        unsigned hit_index = hit_indices[plane_i][hit_i].first;

        std::vector<recob::Hit const*> sharedHits = {};
        sharedHits.reserve(hit_indices[plane_i][hit_i].second.size());
        for (const unsigned shared_hit_index : hit_indices[plane_i][hit_i].second)
          sharedHits.push_back(hits[shared_hit_index].get());

        // Get the location of this point
        geo::Point_t location = GetLocation(track, hits[hit_index], thms[hit_index]);

        // Get the pitch
        double pitch = GetPitch(track, hits[hit_index], thms[hit_index]);

        // And the charge
        double charge = GetCharge(hits[hit_index], sharedHits);

        // Get the EField
        double EField = GetEfield(det_prop, track, hits[hit_index], thms[hit_index]);

        // Angle to the drift electric field (in x direction), in units of degrees
        geo::Vector_t direction = track.DirectionAtPoint(thms[hit_index]->Index());
        double phi = acos(abs(direction.x())) * 180 / M_PI;

        double dQdx = charge / pitch;

        // Normalize out the detector response
        dQdx = Normalize(dQdx,
                         evt,
                         *hits[hit_index],
                         track.LocationAtPoint(thms[hit_index]->Index()),
                         track.DirectionAtPoint(thms[hit_index]->Index()),
                         T0);

        // turn into dEdx
        double dEdx = (fConfig.ChargeMethod() == calo::GnocchiCalorimetry::Config::cmAmplitude) ?
                        fCaloAlg.dEdx_AMP(clock_data,
                                          det_prop,
                                          dQdx,
                                          hits[hit_index]->PeakTime(),
                                          hits[hit_index]->WireID().Plane,
                                          T0,
                                          EField,
                                          phi) :
                        fCaloAlg.dEdx_AREA(clock_data,
                                           det_prop,
                                           dQdx,
                                           hits[hit_index]->PeakTime(),
                                           hits[hit_index]->WireID().Plane,
                                           T0,
                                           EField,
                                           phi);

        // save the length between each pair of hits
        if (xyzs.size() == 0) { lengths.push_back(0.); }
        else {
          lengths.push_back((location - xyzs.back()).r());
        }

        // save stuff
        dEdxs.push_back(dEdx);
        dQdxs.push_back(dQdx);
        pitches.push_back(pitch);
        xyzs.push_back(location);
        kinetic_energy += dEdx * pitch;

        // TODO: FIXME
        // It seems weird that the "trajectory-point-index" actually is the
        // index of the hit... is this a bug in the documentation
        // of anab::Calorimetry?
        //
        // i.e. -- I think this piece of code should actually be:
        // tp_indices.push_back(thms[hit_index]->Index());
        tp_indices.push_back(hits[hit_index].key());

      } // end iterate over hits

      // turn the lengths vector into a residual-range vector and total length
      if (lengths.size() > 1) {
        range = std::accumulate(lengths.begin(), lengths.end(), 0.);

        // check the direction that the hits are going in the track:
        // upstream (end-start) or downstream (start-end)
        bool is_downstream =
          (track.Trajectory().Start() - xyzs[0]).r() +
            (track.Trajectory().End() - xyzs.back()).r() <
          (track.Trajectory().End() - xyzs[0]).r() + (track.Trajectory().Start() - xyzs.back()).r();

        resranges.resize(lengths.size());
        if (is_downstream) {
          resranges[lengths.size() - 1] = lengths.back() / 2.;
          for (int i_len = lengths.size() - 2; i_len >= 0; i_len--) {
            resranges[i_len] = resranges[i_len + 1] + lengths[i_len + 1];
          }
        }
        else {
          resranges[0] = lengths[1] / 2.;
          for (unsigned i_len = 1; i_len < lengths.size(); i_len++) {
            resranges[i_len] = resranges[i_len - 1] + lengths[i_len];
          }
        }
      }

      // save the Calorimetry output
      //
      // Bogus if less than two hits on this plane
      if (lengths.size() > 1) {
        outputCalo->push_back(anab::Calorimetry(kinetic_energy,
                                                dEdxs,
                                                dQdxs,
                                                resranges,
                                                deadwireresranges,
                                                range,
                                                pitches,
                                                xyzs,
                                                tp_indices,
                                                plane));
      }
      else {
        outputCalo->push_back(
          anab::Calorimetry(util::kBogusD, {}, {}, {}, {}, util::kBogusD, {}, {}, {}, plane));
      }

      util::CreateAssn(evt, *outputCalo, tracklist[trk_i], *outputCaloAssn);

    } // end iterate over planes

  } // end iterate over tracks

  evt.put(std::move(outputCalo));
  evt.put(std::move(outputCaloAssn));

  return;
}

std::vector<std::vector<OrganizedHits>> calo::GnocchiCalorimetry::OrganizeHits(
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<const recob::TrackHitMeta*>& thms,
  const recob::Track& track,
  unsigned nplanes)
{
  // charge is computed per hit -- we organize hits indivudally
  if (fConfig.ChargeMethod() == calo::GnocchiCalorimetry::Config::cmIntegral ||
      fConfig.ChargeMethod() == calo::GnocchiCalorimetry::Config::cmAmplitude) {
    return OrganizeHitsIndividual(hits, thms, track, nplanes);
  }
  // charge is computed per snippet -- we organize hits by snippet
  else {
    return OrganizeHitsSnippets(hits, thms, track, nplanes);
  }
}

std::vector<std::vector<OrganizedHits>> calo::GnocchiCalorimetry::OrganizeHitsIndividual(
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<const recob::TrackHitMeta*>& thms,
  const recob::Track& track,
  unsigned nplanes)
{
  std::vector<std::vector<OrganizedHits>> ret(nplanes);
  for (unsigned i = 0; i < hits.size(); i++) {
    if (HitIsValid(thms[i], track)) { ret[hits[i]->WireID().Plane].push_back({i, {}}); }
  }

  return ret;
}

std::vector<std::vector<OrganizedHits>> calo::GnocchiCalorimetry::OrganizeHitsSnippets(
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<const recob::TrackHitMeta*>& thms,
  const recob::Track& track,
  unsigned nplanes)
{
  // In this case, we need to only accept one hit in each snippet
  // Snippets are counted by the Start, End, and Wire. If all these are the same for a hit, then they are on the same snippet.
  //
  // If there are multiple valid hits on the same snippet, we need a way to pick the best one.
  // (TODO: find a good way). The current method is to take the one with the highest charge integral.
  struct HitIdentifier {
    int startTick;
    int endTick;
    int wire;
    float integral;

    // construct
    explicit HitIdentifier(const recob::Hit& hit)
      : startTick(hit.StartTick())
      , endTick(hit.EndTick())
      , wire(hit.WireID().Wire)
      , integral(hit.Integral())
    {}

    // Defines whether two hits are on the same snippet
    inline bool operator==(const HitIdentifier& rhs) const
    {
      return startTick == rhs.startTick && endTick == rhs.endTick && wire == rhs.wire;
    }

    // Defines which hit to pick between two both on the same snippet
    inline bool operator>(const HitIdentifier& rhs) const { return integral > rhs.integral; }
  };

  std::vector<std::vector<OrganizedHits>> ret(nplanes);
  std::vector<std::vector<HitIdentifier>> hit_idents(nplanes);
  for (unsigned i = 0; i < hits.size(); i++) {
    if (HitIsValid(thms[i], track)) {
      HitIdentifier this_ident(*hits[i]);

      // check if we have found a hit on this snippet before
      bool found_snippet = false;
      auto const plane = hits[i]->WireID().Plane;
      for (unsigned j = 0; j < ret[plane].size(); j++) {
        if (this_ident == hit_idents[plane][j]) {
          found_snippet = true;
          if (this_ident > hit_idents[plane][j]) {
            ret[plane][j].second.push_back(ret[plane][j].first);
            ret[plane][j].first = i;
            hit_idents[plane][j] = this_ident;
          }
          else {
            ret[plane][j].second.push_back(i);
          }
          break;
        }
      }
      if (!found_snippet) {
        ret[plane].push_back({i, {}});
        hit_idents[plane].push_back(this_ident);
      }
    }
  }
  return ret;
}

bool calo::GnocchiCalorimetry::HitIsValid(const recob::TrackHitMeta* thm, const recob::Track& track)
{
  if (thm->Index() == int_max_as_unsigned_int) return false;
  return track.HasValidPoint(thm->Index());
}

geo::Point_t calo::GnocchiCalorimetry::GetLocation(const recob::Track& track,
                                                   const art::Ptr<recob::Hit> hit,
                                                   const recob::TrackHitMeta* meta)
{
  geo::Point_t loc = track.LocationAtPoint(meta->Index());
  return !fConfig.TrackIsFieldDistortionCorrected() ? WireToTrajectoryPosition(loc, hit->WireID()) :
                                                      loc;
}

geo::Point_t calo::GnocchiCalorimetry::WireToTrajectoryPosition(const geo::Point_t& loc,
                                                                const geo::TPCID& tpc)
{
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t ret = loc;

  if (sce->EnableCalSpatialSCE() && fConfig.FieldDistortion()) {
    geo::Vector_t offset = sce->GetCalPosOffsets(ret, tpc.TPC);

    ret.SetX(ret.X() + fConfig.FieldDistortionCorrectionXSign() * offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
}

geo::Point_t calo::GnocchiCalorimetry::GetLocationAtWires(const recob::Track& track,
                                                          const art::Ptr<recob::Hit> hit,
                                                          const recob::TrackHitMeta* meta)
{
  geo::Point_t loc = track.LocationAtPoint(meta->Index());
  return fConfig.TrackIsFieldDistortionCorrected() ? TrajectoryToWirePosition(loc, hit->WireID()) :
                                                     loc;
}

geo::Point_t calo::GnocchiCalorimetry::TrajectoryToWirePosition(const geo::Point_t& loc,
                                                                const geo::TPCID& tpc)
{
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  art::ServiceHandle<geo::Geometry const> geom;

  geo::Point_t ret = loc;

  if (sce->EnableCalSpatialSCE() && fConfig.FieldDistortion()) {
    // Returned X is the drift -- multiply by the drift direction to undo this
    int corr = geom->TPC(tpc).DriftDir().X();

    geo::Vector_t offset = sce->GetPosOffsets(ret);

    ret.SetX(ret.X() + corr * fConfig.FieldDistortionCorrectionXSign() * offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
}

double calo::GnocchiCalorimetry::GetPitch(const recob::Track& track,
                                          const art::Ptr<recob::Hit> hit,
                                          const recob::TrackHitMeta* meta)
{
  art::ServiceHandle<geo::Geometry const> geom;
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  double angleToVert =
    geom->WireAngleToVertical(hit->View(), hit->WireID().asPlaneID()) - 0.5 * ::util::pi<>();

  geo::Vector_t dir;

  // "dir" should be the direction that the wires see. If the track already has the field
  // distortion corrections applied, then we need to de-apply them to get the direction as
  // seen by the wire planes
  if (sce->EnableCalSpatialSCE() && fConfig.FieldDistortion() &&
      fConfig.TrackIsFieldDistortionCorrected()) {
    geo::Point_t loc = track.LocationAtPoint(meta->Index());

    // compute the dir of the track trajectory
    geo::Vector_t track_dir = track.DirectionAtPoint(meta->Index());
    geo::Point_t loc_mdx = loc - track_dir * (geom->WirePitch(hit->View()) / 2.);
    geo::Point_t loc_pdx = loc + track_dir * (geom->WirePitch(hit->View()) / 2.);

    loc_mdx = TrajectoryToWirePosition(loc_mdx, hit->WireID());
    loc_pdx = TrajectoryToWirePosition(loc_pdx, hit->WireID());

    // Direction at wires
    dir = (loc_pdx - loc_mdx) / (loc_mdx - loc_pdx).r();
  }
  // If there is no space charge or the track is not yet corrected, then the dir
  // is the track is what we want
  else {
    dir = track.DirectionAtPoint(meta->Index());
  }

  double cosgamma = std::abs(std::sin(angleToVert) * dir.Y() + std::cos(angleToVert) * dir.Z());
  double pitch;
  if (cosgamma) { pitch = geom->WirePitch(hit->View()) / cosgamma; }
  else {
    pitch = 0.;
  }

  // now take the pitch computed on the wires and correct it back to the particle trajectory
  geo::Point_t loc_w = GetLocationAtWires(track, hit, meta);

  geo::Point_t locw_pdx_traj = WireToTrajectoryPosition(loc_w + pitch * dir, hit->WireID());
  geo::Point_t loc = WireToTrajectoryPosition(loc_w, hit->WireID());

  pitch = (locw_pdx_traj - loc).R();

  return pitch;
}

double calo::GnocchiCalorimetry::GetCharge(const art::Ptr<recob::Hit> hit,
                                           const std::vector<recob::Hit const*>& sharedHits)
{
  switch (fConfig.ChargeMethod()) {
  case calo::GnocchiCalorimetry::Config::cmIntegral: return hit->Integral();
  case calo::GnocchiCalorimetry::Config::cmAmplitude: return hit->PeakAmplitude();
  case calo::GnocchiCalorimetry::Config::cmSummedADC: return hit->SummedADC();
  case calo::GnocchiCalorimetry::Config::cmSummedIntegral:
    return std::accumulate(
      sharedHits.cbegin(),
      sharedHits.cend(),
      hit->Integral(),
      [](double sum, recob::Hit const* sharedHit) { return sum + sharedHit->Integral(); });
  default: return 0.;
  }
  return 0.;
}

double calo::GnocchiCalorimetry::Normalize(double dQdx,
                                           const art::Event& e,
                                           const recob::Hit& h,
                                           const geo::Point_t& location,
                                           const geo::Vector_t& direction,
                                           double t0)
{
  double ret = dQdx;
  for (auto const& nt : fNormTools) {
    ret = nt->Normalize(ret, e, h, location, direction, t0);
  }

  return ret;
}

double calo::GnocchiCalorimetry::GetEfield(const detinfo::DetectorPropertiesData& dprop,
                                           const recob::Track& track,
                                           const art::Ptr<recob::Hit> hit,
                                           const recob::TrackHitMeta* meta)
{
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  double EField = dprop.Efield();
  if (sce->EnableSimEfieldSCE() && fConfig.FieldDistortionEfield()) {
    // Gets relative E field Distortions
    geo::Vector_t EFieldOffsets = sce->GetEfieldOffsets(GetLocation(track, hit, meta));
    // Add 1 in X direction as this is the direction of the drift field
    EFieldOffsets = EFieldOffsets + geo::Vector_t{1, 0, 0};
    // Convert to Absolute E Field from relative
    EFieldOffsets = EField * EFieldOffsets;
    // We only care about the magnitude for recombination
    EField = EFieldOffsets.r();
  }
  return EField;
}

DEFINE_ART_MODULE(calo::GnocchiCalorimetry)

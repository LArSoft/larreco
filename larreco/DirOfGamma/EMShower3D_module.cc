////////////////////////////////////////////////////////////////////////
// Class:       EMShower3D
// Module Type: producer
// File:        EMShower3D_module.cc
// Authors: dorota.stefan@cern.ch robert.sulej@cern.ch
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"

#include <memory>

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/DirOfGamma/DirOfGamma.h"

struct IniSeg {
  size_t idcl1;
  size_t idcl2;
  size_t idcl3;
  size_t view1;
  size_t view2;
  size_t view3;
  pma::Track3D* track;
  std::vector<art::Ptr<recob::Hit>> hits1;
  std::vector<art::Ptr<recob::Hit>> hits2;
  std::vector<art::Ptr<recob::Hit>> hits3;
};

namespace ems {
  class EMShower3D;
}

class ems::EMShower3D : public art::EDProducer {
public:
  explicit EMShower3D(fhicl::ParameterSet const& p);

  EMShower3D(EMShower3D const&) = delete;
  EMShower3D(EMShower3D&&) = delete;
  EMShower3D& operator=(EMShower3D const&) = delete;
  EMShower3D& operator=(EMShower3D&&) = delete;

private:
  void produce(art::Event& e) override;

  recob::Track ConvertFrom(pma::Track3D& src);
  recob::Track ConvertFrom2(pma::Track3D& src);
  recob::Cluster ConvertFrom(const std::vector<art::Ptr<recob::Hit>>& src);

  std::vector<ems::DirOfGamma*> CollectShower2D(art::Event const& e);

  void Link(art::Event const& e, std::vector<ems::DirOfGamma*> input);

  // Remove tracks which are too close to each other
  void Reoptimize();

  void Make3DSeg(art::Event const& e, std::vector<ems::DirOfGamma*> pair);

  bool Validate(art::Event const& e, const pma::Track3D& src, size_t plane);
  bool Validate(std::vector<ems::DirOfGamma*> input,
                size_t id1,
                size_t id2,
                size_t c1,
                size_t c2,
                size_t plane3);

  void FilterOutSmallParts(double r2d,
                           const std::vector<art::Ptr<recob::Hit>>& hits_in,
                           std::vector<art::Ptr<recob::Hit>>& hits_out);

  bool GetCloseHits(double r2d,
                    const std::vector<art::Ptr<recob::Hit>>& hits_in,
                    std::vector<size_t>& used,
                    std::vector<art::Ptr<recob::Hit>>& hits_out);

  bool Has(const std::vector<size_t>& v, size_t idx);

  size_t LinkCandidates(art::Event const& e, std::vector<ems::DirOfGamma*> input, size_t id);

  std::vector<IniSeg> fInisegs;
  std::vector<IniSeg> fSeltracks;
  std::vector<IniSeg> fPMA3D;

  std::vector<std::vector<art::Ptr<recob::Hit>>> fClusters;

  std::vector<size_t> fClustersNotUsed;
  std::vector<size_t> fTracksNotUsed;

  unsigned int fTrkIndex;
  unsigned int fClIndex;
  unsigned int fIniIndex;

  std::string fCluModuleLabel;
  std::string fTrk3DModuleLabel;

  pma::ProjectionMatchingAlg fProjectionMatchingAlg;
  calo::CalorimetryAlg fCalorimetryAlg;

  art::Handle<std::vector<recob::Cluster>> fCluListHandle;
};

ems::EMShower3D::EMShower3D(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fProjectionMatchingAlg(p.get<fhicl::ParameterSet>("ProjectionMatchingAlg"))
  , fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  fCluModuleLabel = p.get<std::string>("ClustersModuleLabel");
  fTrk3DModuleLabel = p.get<std::string>("Trk3DModuleLabel");

  produces<std::vector<recob::Track>>();
  produces<std::vector<recob::Vertex>>();
  produces<std::vector<recob::Cluster>>();
  produces<std::vector<recob::SpacePoint>>();
  produces<art::Assns<recob::Track, recob::Hit>>();
  produces<art::Assns<recob::Track, recob::Vertex>>();
  produces<art::Assns<recob::Cluster, recob::Hit>>();
  produces<art::Assns<recob::Track, recob::SpacePoint>>();
  produces<art::Assns<recob::SpacePoint, recob::Hit>>();
  produces<art::Assns<recob::Track, recob::Cluster>>();
}

recob::Cluster
ems::EMShower3D::ConvertFrom(const std::vector<art::Ptr<recob::Hit>>& src)
{

  return recob::Cluster(0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        0.0F,
                        src.size(),
                        0.0F,
                        0.0F,
                        fClIndex,
                        src[0]->View(),
                        src[0]->WireID().planeID());
}

recob::Track
ems::EMShower3D::ConvertFrom(pma::Track3D& src)
{
  auto const* geom = lar::providerFrom<geo::Geometry>();
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  double avdrift = (src.front()->Point3D().X() + src.back()->Point3D().X()) * 0.5;
  unsigned int nplanes = geom->Nplanes(src.front()->TPC(), src.front()->Cryo());
  size_t nusedhitsmax = 0;
  int bestplane = -1;
  for (unsigned int p = 0; p < nplanes; ++p) {
    unsigned int nusedP = 0;
    fProjectionMatchingAlg.selectInitialHits(src, p, &nusedP);

    if (nusedP > nusedhitsmax) {
      nusedhitsmax = nusedP;
      bestplane = int(p);
    }
  }

  std::vector<std::vector<double>> vdedx;
  std::vector<double> dedx;

  for (unsigned int p = 0; p < nplanes; ++p) {
    unsigned int nusedP = 0;
    double dqdxplane = fProjectionMatchingAlg.selectInitialHits(src, p, &nusedP);
    double timeP = detprop->ConvertXToTicks(avdrift, p, src.front()->TPC(), src.front()->Cryo());
    double dEdxplane = fCalorimetryAlg.dEdx_AREA(dqdxplane, timeP, p);
    dedx.push_back(dEdxplane);
    if (int(p) == bestplane)
      dedx.push_back(1);
    else
      dedx.push_back(0);
    vdedx.push_back(dedx);
  }

  std::vector<TVector3> xyz, dircos;

  for (size_t i = 0; i < src.size(); i++) {
    xyz.push_back(src[i]->Point3D());

    if (i < src.size() - 1) {
      TVector3 dc(src[i + 1]->Point3D());
      dc -= src[i]->Point3D();
      dc *= 1.0 / dc.Mag();
      dircos.push_back(dc);
    }
    else
      dircos.push_back(dircos.back());
  }

  return recob::Track(recob::TrackTrajectory(recob::tracking::convertCollToPoint(xyz),
                                             recob::tracking::convertCollToVector(dircos),
                                             recob::Track::Flags_t(xyz.size()),
                                             false),
                      0,
                      -1.,
                      0,
                      recob::tracking::SMatrixSym55(),
                      recob::tracking::SMatrixSym55(),
                      fIniIndex);
}

recob::Track
ems::EMShower3D::ConvertFrom2(pma::Track3D& src)
{

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* geom = lar::providerFrom<geo::Geometry>();

  double avdrift = (src.front()->Point3D().X() + src.back()->Point3D().X()) * 0.5;
  unsigned int nplanes = geom->Nplanes(src.front()->TPC(), src.front()->Cryo());
  size_t nusedhitsmax = 0;
  int bestplane = -1;
  for (unsigned int p = 0; p < nplanes; ++p) {
    unsigned int nusedP = 0;
    fProjectionMatchingAlg.selectInitialHits(src, p, &nusedP);

    if (nusedP > nusedhitsmax) {
      nusedhitsmax = nusedP;
      bestplane = int(p);
    }
  }

  std::vector<std::vector<double>> vdedx;
  std::vector<double> dedx;

  for (unsigned int p = 0; p < nplanes; ++p) {
    unsigned int nusedP = 0;
    double dqdxplane = fProjectionMatchingAlg.selectInitialHits(src, p, &nusedP);
    double timeP = detprop->ConvertXToTicks(avdrift, p, src.front()->TPC(), src.front()->Cryo());
    double dEdxplane = fCalorimetryAlg.dEdx_AREA(dqdxplane, timeP, p);
    dedx.push_back(dEdxplane);
    if (int(p) == bestplane)
      dedx.push_back(1);
    else
      dedx.push_back(0);
    vdedx.push_back(dedx);
  }

  std::vector<TVector3> xyz, dircos;

  for (size_t i = 0; i < src.size(); i++) {
    xyz.push_back(src[i]->Point3D());

    if (i < src.size() - 1) {
      TVector3 dc(src[i + 1]->Point3D());
      dc -= src[i]->Point3D();
      dc *= 1.0 / dc.Mag();
      dircos.push_back(dc);
    }
    else
      dircos.push_back(dircos.back());
  }

  return recob::Track(recob::TrackTrajectory(recob::tracking::convertCollToPoint(xyz),
                                             recob::tracking::convertCollToVector(dircos),
                                             recob::Track::Flags_t(xyz.size()),
                                             false),
                      0,
                      -1.,
                      0,
                      recob::tracking::SMatrixSym55(),
                      recob::tracking::SMatrixSym55(),
                      fIniIndex);
}

void
ems::EMShower3D::produce(art::Event& e)
{
  art::ServiceHandle<geo::Geometry const> geom;
  fSeltracks.clear();
  fInisegs.clear();
  fClusters.clear();
  fPMA3D.clear();
  fClustersNotUsed.clear();

  std::unique_ptr<std::vector<recob::Track>> tracks(new std::vector<recob::Track>);
  std::unique_ptr<std::vector<recob::Vertex>> vertices(new std::vector<recob::Vertex>);
  std::unique_ptr<std::vector<recob::Cluster>> clusters(new std::vector<recob::Cluster>);
  std::unique_ptr<std::vector<recob::SpacePoint>> allsp(new std::vector<recob::SpacePoint>);

  std::unique_ptr<art::Assns<recob::Track, recob::Hit>> trk2hit(
    new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Track, recob::Vertex>> trk2vtx(
    new art::Assns<recob::Track, recob::Vertex>);
  std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>> cl2hit(
    new art::Assns<recob::Cluster, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Track, recob::Cluster>> trk2cl(
    new art::Assns<recob::Track, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint>> trk2sp(
    new art::Assns<recob::Track, recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit>> sp2hit(
    new art::Assns<recob::SpacePoint, recob::Hit>);

  if (e.getByLabel(fCluModuleLabel, fCluListHandle)) {
    fClustersNotUsed.clear();
    fInisegs.clear();
    art::FindManyP<recob::Hit> fb(fCluListHandle, e, fCluModuleLabel);

    for (size_t id = 0; id < fCluListHandle->size(); id++) {
      std::vector<art::Ptr<recob::Hit>> hitlist;
      hitlist = fb.at(id);

      if (hitlist.size() > 5) fClustersNotUsed.push_back(id);
    }

    std::vector<ems::DirOfGamma*> showernviews = CollectShower2D(e);

    Link(e, showernviews);

    while (fInisegs.size()) {
      fSeltracks.push_back(fInisegs[0]);
      fInisegs.erase(fInisegs.begin() + 0);
    }

    Reoptimize();

    // conversion from pma track to recob::track

    size_t spStart = 0, spEnd = 0;
    double sp_pos[3], sp_err[6], vtx_pos[3];
    for (size_t i = 0; i < 6; i++)
      sp_err[i] = 1.0;

    fTrkIndex = 0;

    for (auto const trk : fSeltracks) {
      tracks->push_back(ConvertFrom(*(trk.track)));

      vtx_pos[0] = trk.track->front()->Point3D().X();
      vtx_pos[1] = trk.track->front()->Point3D().Y();
      vtx_pos[2] = trk.track->front()->Point3D().Z();
      vertices->push_back(recob::Vertex(vtx_pos, fTrkIndex));

      fTrkIndex++;

      std::vector<art::Ptr<recob::Cluster>> cl2d;
      cl2d.push_back(art::Ptr<recob::Cluster>(fCluListHandle, trk.idcl1));
      cl2d.push_back(art::Ptr<recob::Cluster>(fCluListHandle, trk.idcl2));

      std::vector<art::Ptr<recob::Hit>> hits2d;
      art::PtrVector<recob::Hit> sp_hits;

      spStart = allsp->size();
      for (int h = trk.track->size() - 1; h >= 0; h--) {
        pma::Hit3D* h3d = (*trk.track)[h];
        hits2d.push_back(h3d->Hit2DPtr());

        if ((h == 0) || (sp_pos[0] != h3d->Point3D().X()) || (sp_pos[1] != h3d->Point3D().Y()) ||
            (sp_pos[2] != h3d->Point3D().Z())) {
          if (sp_hits.size()) // hits assigned to the previous sp
          {
            util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
            sp_hits.clear();
          }
          sp_pos[0] = h3d->Point3D().X();
          sp_pos[1] = h3d->Point3D().Y();
          sp_pos[2] = h3d->Point3D().Z();
          allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
        }
        sp_hits.push_back(h3d->Hit2DPtr());
      }
      if (sp_hits.size()) // hits assigned to the last sp
      {
        util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
      }
      spEnd = allsp->size();

      if (vertices->size()) {
        size_t vtx_idx = (size_t)(vertices->size() - 1);
        util::CreateAssn(*this, e, *tracks, *vertices, *trk2vtx, vtx_idx, vtx_idx + 1);
      }

      if (cl2d.size()) { util::CreateAssn(*this, e, *tracks, cl2d, *trk2cl); }

      if (hits2d.size()) {
        util::CreateAssn(*this, e, *tracks, *allsp, *trk2sp, spStart, spEnd);
        util::CreateAssn(*this, e, *tracks, hits2d, *trk2hit);
      }
    }

    fIniIndex = fTrkIndex + 1;
    for (auto const trk : fPMA3D) {
      tracks->push_back(ConvertFrom2(*(trk.track)));

      fIniIndex++;

      std::vector<art::Ptr<recob::Cluster>> cl2d;
      cl2d.push_back(art::Ptr<recob::Cluster>(fCluListHandle, trk.idcl1));
      cl2d.push_back(art::Ptr<recob::Cluster>(fCluListHandle, trk.idcl2));

      std::vector<art::Ptr<recob::Hit>> hits2d;
      art::PtrVector<recob::Hit> sp_hits;

      spStart = allsp->size();
      for (int h = trk.track->size() - 1; h >= 0; h--) {
        pma::Hit3D* h3d = (*trk.track)[h];
        hits2d.push_back(h3d->Hit2DPtr());

        if ((h == 0) || (sp_pos[0] != h3d->Point3D().X()) || (sp_pos[1] != h3d->Point3D().Y()) ||
            (sp_pos[2] != h3d->Point3D().Z())) {
          if (sp_hits.size()) // hits assigned to the previous sp
          {
            util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
            sp_hits.clear();
          }
          sp_pos[0] = h3d->Point3D().X();
          sp_pos[1] = h3d->Point3D().Y();
          sp_pos[2] = h3d->Point3D().Z();
          allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
        }
        sp_hits.push_back(h3d->Hit2DPtr());
      }
      if (sp_hits.size()) // hits assigned to the last sp
      {
        util::CreateAssn(*this, e, *allsp, sp_hits, *sp2hit);
      }
      spEnd = allsp->size();

      if (cl2d.size()) { util::CreateAssn(*this, e, *tracks, cl2d, *trk2cl); }

      if (hits2d.size()) {
        util::CreateAssn(*this, e, *tracks, *allsp, *trk2sp, spStart, spEnd);
        util::CreateAssn(*this, e, *tracks, hits2d, *trk2hit);
      }
    }

    // create cluster from hits, which were an input to find initial part of the cascade.
    fClIndex = 0;
    for (auto const& cl : fClusters)
      if (cl.size()) {
        clusters->push_back(ConvertFrom(cl));
        fClIndex++;

        util::CreateAssn(*this, e, *clusters, cl, *cl2hit);
      }

    for (unsigned int i = 0; i < showernviews.size(); i++)
      delete showernviews[i];

    for (unsigned int i = 0; i < fSeltracks.size(); i++)
      delete fSeltracks[i].track;

    for (unsigned int i = 0; i < fInisegs.size(); i++)
      delete fInisegs[i].track;

    for (unsigned int i = 0; i < fPMA3D.size(); i++)
      delete fPMA3D[i].track;
  }

  e.put(std::move(tracks));
  e.put(std::move(vertices));
  e.put(std::move(clusters));
  e.put(std::move(allsp));

  e.put(std::move(trk2hit));
  e.put(std::move(trk2vtx));
  e.put(std::move(cl2hit));
  e.put(std::move(trk2cl));
  e.put(std::move(trk2sp));
  e.put(std::move(sp2hit));
}

void
ems::EMShower3D::Reoptimize()
{
  if (!fSeltracks.size()) return;
  const float min_dist = 3.0F;
  size_t ta = 0;
  while (ta < (fSeltracks.size() - 1)) {
    size_t tb = ta + 1;
    bool found = false;
    while (tb < fSeltracks.size()) {
      if (ta == tb) {
        tb++;
        continue;
      }

      TVector3 p1 = fSeltracks[ta].track->front()->Point3D();
      TVector3 p2 = fSeltracks[tb].track->front()->Point3D();
      float dist = std::sqrt(pma::Dist2(p1, p2));

      if (dist < min_dist)
        if ((fSeltracks[ta].idcl1 == fSeltracks[tb].idcl1) ||
            (fSeltracks[ta].idcl1 == fSeltracks[tb].idcl2) ||
            (fSeltracks[ta].idcl2 == fSeltracks[tb].idcl1) ||
            (fSeltracks[ta].idcl2 == fSeltracks[tb].idcl2)) {
          found = true;
          size_t view3 = fSeltracks[ta].view1;
          size_t idcl3 = fSeltracks[ta].idcl1;
          std::vector<art::Ptr<recob::Hit>> hits3 = fSeltracks[ta].hits1;
          std::vector<art::Ptr<recob::Hit>> hits = fSeltracks[ta].hits1;
          for (size_t h = 0; h < fSeltracks[ta].hits2.size(); ++h)
            hits.push_back(fSeltracks[ta].hits2[h]);

          if ((fSeltracks[tb].view1 != fSeltracks[ta].view1) &&
              (fSeltracks[tb].view1 != fSeltracks[ta].view2)) {
            view3 = fSeltracks[tb].view1;
            for (size_t h = 0; h < fSeltracks[tb].hits1.size(); ++h)
              hits.push_back(fSeltracks[tb].hits1[h]);
          }
          if ((fSeltracks[tb].view2 != fSeltracks[ta].view1) &&
              (fSeltracks[tb].view2 != fSeltracks[ta].view2)) {
            view3 = fSeltracks[tb].view2;
            for (size_t h = 0; h < fSeltracks[tb].hits2.size(); ++h)
              hits.push_back(fSeltracks[tb].hits2[h]);
          }

          if ((view3 == fSeltracks[ta].view1) || (view3 == fSeltracks[ta].view2)) {
            delete fSeltracks[ta].track;
            fSeltracks.erase(fSeltracks.begin() + ta);
          }
          else {
            pma::Track3D* track = fProjectionMatchingAlg.buildSegment(hits);

            if (pma::Dist2(track->back()->Point3D(), fSeltracks[ta].track->front()->Point3D()) <
                pma::Dist2(track->front()->Point3D(), fSeltracks[ta].track->front()->Point3D()))
              track->Flip();

            IniSeg initrack;
            initrack.idcl1 = fSeltracks[ta].idcl1;
            initrack.idcl3 = idcl3;
            initrack.view1 = fSeltracks[ta].view1;
            initrack.view3 = view3;
            initrack.hits1 = fSeltracks[ta].hits1;
            initrack.hits3 = hits3;
            initrack.idcl2 = fSeltracks[ta].idcl2;
            initrack.view2 = fSeltracks[ta].view2;
            initrack.hits2 = fSeltracks[ta].hits2;
            initrack.track = track;

            delete fSeltracks[tb].track;
            delete fSeltracks[ta].track;
            fSeltracks.erase(fSeltracks.begin() + tb);
            fSeltracks.erase(fSeltracks.begin() + ta);
            fSeltracks.push_back(initrack);
          }
        }

      if (found) break;
      tb++;
    }

    if (!found) ta++;
  }
}

std::vector<ems::DirOfGamma*>
ems::EMShower3D::CollectShower2D(art::Event const& e)
{
  std::vector<ems::DirOfGamma*> input;

  if (e.getByLabel(fCluModuleLabel, fCluListHandle)) {
    art::FindManyP<recob::Hit> fb(fCluListHandle, e, fCluModuleLabel);
    for (unsigned int c = 0; c < fCluListHandle->size(); c++) {
      std::vector<art::Ptr<recob::Hit>> hitlist;
      hitlist = fb.at(c);

      if (hitlist.size() > 5) {
        std::vector<art::Ptr<recob::Hit>> hits_out;
        FilterOutSmallParts(2.0, hitlist, hits_out);

        if (hits_out.size() > 5) {
          fClusters.push_back(hits_out);

          ems::DirOfGamma* sh = new ems::DirOfGamma(hits_out, 14, c);

          if (sh->GetHits2D().size()) input.push_back(sh);
        }
      }
    }
  }

  return input;
}

void
ems::EMShower3D::Link(art::Event const& e, std::vector<ems::DirOfGamma*> input)
{
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  std::vector<std::vector<size_t>> saveids;
  std::vector<size_t> saveidsnotusedcls;
  size_t i = 0;

  while (i < input.size()) {
    if (!input[i]->GetCandidates().size()) {
      i++;
      continue;
    }

    double mindist = 1.0; // cm
    std::vector<ems::DirOfGamma*> pairs;

    size_t startview = input[i]->GetFirstHit()->WireID().Plane;
    size_t tpc = input[i]->GetFirstHit()->WireID().TPC;
    size_t cryo = input[i]->GetFirstHit()->WireID().Cryostat;

    float t1 = detprop->ConvertTicksToX(input[i]->GetFirstHit()->PeakTime(), startview, tpc, cryo);

    unsigned int idsave = 0;
    for (unsigned int j = 0; j < input.size(); j++) {
      if (!input[j]->GetCandidates().size()) continue;

      size_t secondview = input[j]->GetFirstHit()->WireID().Plane;
      size_t tpc_j = input[j]->GetFirstHit()->WireID().TPC;
      size_t cryo_j = input[j]->GetFirstHit()->WireID().Cryostat;

      if ((i != j) && (secondview != startview) && (tpc == tpc_j) && (cryo == cryo_j)) {
        float t2 =
          detprop->ConvertTicksToX(input[j]->GetFirstHit()->PeakTime(), secondview, tpc_j, cryo_j);
        float dist = fabs(t2 - t1);

        if (dist < mindist) {
          mindist = dist;
          pairs.clear();
          pairs.push_back(input[i]);
          pairs.push_back(input[j]);
          idsave = j;
        }
      }
    }

    bool exist = false;
    for (unsigned int v = 0; v < saveids.size(); v++)
      if ((saveids[v][0] == i) || (saveids[v][0] == idsave))
        if ((saveids[v][1] == i) || (saveids[v][1] == idsave)) exist = true;

    if (pairs.size()) {
      if (!exist) Make3DSeg(e, pairs);

      std::vector<size_t> ids;
      ids.push_back(i);
      ids.push_back(idsave);
      saveids.push_back(ids);
    }
    else {
      saveidsnotusedcls.push_back(i);
    }

    i++;
  }

  i = 0;
  while (i < saveidsnotusedcls.size()) {
    LinkCandidates(e, input, i);
    i++;
  }
}

size_t
ems::EMShower3D::LinkCandidates(art::Event const& e, std::vector<ems::DirOfGamma*> input, size_t id)
{
  art::ServiceHandle<geo::Geometry const> geom;

  size_t index = id;
  bool found = false;

  if (input[id]->GetCandidates().size() < 2) { return index; }

  double mindist = 3.0; // cm
  std::vector<ems::DirOfGamma*> pairs;

  size_t idcsave = 0;
  size_t idcjsave = 0;
  size_t c = 0;
  size_t idsave = 0;
  while (c < input[id]->GetCandidates().size()) {

    size_t startview = input[id]->GetCandidates()[c].GetPlane();
    size_t tpc = input[id]->GetCandidates()[c].GetTPC();
    size_t cryo = input[id]->GetCandidates()[c].GetCryo();

    float t1 = input[id]->GetCandidates()[c].GetPosition().Y(); // y --> drift in 2D space.

    // loop over 2D showers
    for (size_t j = 0; j < input.size(); ++j) {
      if (!input[j]->GetCandidates().size()) continue;
      if (j == id) continue;

      // loop over candidates
      for (size_t cj = 0; cj < input[j]->GetCandidates().size(); ++cj) {
        size_t secondview = input[j]->GetCandidates()[cj].GetPlane();
        size_t tpc_j = input[j]->GetCandidates()[cj].GetTPC();
        size_t cryo_j = input[j]->GetCandidates()[cj].GetCryo();

        size_t thirdview = startview;

        const geo::CryostatGeo& cryostat = geom->Cryostat(cryo);
        for (size_t p = 0; p < cryostat.MaxPlanes(); p++)
          if ((p == startview) || (p == secondview)) { continue; }
          else {
            thirdview = p;
            break;
          }

        if ((startview != secondview) && (tpc == tpc_j) &&
            (cryo == cryo_j)) // && Validate(input, id, cj, thirdview))
        {
          float t2 = input[j]->GetCandidates()[cj].GetPosition().Y();
          float dist = fabs(t2 - t1);

          if ((dist < mindist) && Validate(input, id, j, c, cj, thirdview)) {
            mindist = dist;
            pairs.clear();
            pairs.push_back(input[id]);
            pairs.push_back(input[j]);
            idsave = j;
            index = j;
            idcsave = c;
            idcjsave = cj;
            found = true;
          }
        }
      }
    }

    c++;
  }

  if (found && pairs.size()) {
    input[id]->SetIdCandidate(idcsave);
    input[idsave]->SetIdCandidate(idcjsave);
    Make3DSeg(e, pairs);
  }

  return index;
}

void
ems::EMShower3D::Make3DSeg(art::Event const& e, std::vector<ems::DirOfGamma*> pair)
{
  if (pair.size() < 2) return;

  // to build a track correctly 2d hits must belong to the same tpc
  size_t tpc1 = pair[0]->GetFirstHit()->WireID().TPC;
  size_t tpc2 = pair[1]->GetFirstHit()->WireID().TPC;

  std::vector<art::Ptr<recob::Hit>> vec1 = pair[0]->GetIniHits();
  std::vector<art::Ptr<recob::Hit>> vec2 = pair[1]->GetIniHits();

  if ((vec1.size() < 3) && (vec2.size() < 3)) return;

  std::vector<art::Ptr<recob::Hit>> hitscl1uniquetpc;
  std::vector<art::Ptr<recob::Hit>> hitscl2uniquetpc;

  if (tpc1 == tpc2)
    for (size_t i = 0; i < vec1.size(); ++i)
      for (size_t j = 0; j < vec2.size(); ++j)
        if ((vec1[i]->WireID().TPC == vec2[j]->WireID().TPC) && (tpc1 == vec2[j]->WireID().TPC)) {
          hitscl1uniquetpc.push_back(vec1[i]);
          hitscl2uniquetpc.push_back(vec2[j]);
        }

  if ((hitscl1uniquetpc.size() > 2) && (hitscl2uniquetpc.size() > 2)) {
    pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(hitscl1uniquetpc, hitscl2uniquetpc);

    //pma::Track3D* trk = fProjectionMatchingAlg.buildSegment(vec1, vec2);

    // turn the track that front is at vertex - easier to handle associations.
    if ((trk->back()->Hit2DPtr() == pair[0]->GetFirstHit()) ||
        (trk->back()->Hit2DPtr() == pair[1]->GetFirstHit()))
      trk->Flip();

    IniSeg initrack;
    initrack.idcl1 = pair[0]->GetIdCl();
    initrack.view1 = pair[0]->GetFirstHit()->WireID().Plane;
    initrack.hits1 = hitscl1uniquetpc;
    initrack.idcl2 = pair[1]->GetIdCl();
    initrack.view2 = pair[1]->GetFirstHit()->WireID().Plane;
    initrack.hits2 = hitscl2uniquetpc;
    initrack.track = trk;

    fInisegs.push_back(initrack);
  }
}

bool
ems::EMShower3D::Validate(std::vector<ems::DirOfGamma*> input,
                          size_t id1,
                          size_t id2,
                          size_t c1,
                          size_t c2,
                          size_t plane3)
{
  bool result = false;
  if (id1 == id2) return false;

  std::vector<art::Ptr<recob::Hit>> vec1 = input[id1]->GetCandidates()[c1].GetIniHits();
  std::vector<art::Ptr<recob::Hit>> vec2 = input[id2]->GetCandidates()[c2].GetIniHits();

  if ((vec1.size() < 3) || (vec2.size() < 3)) return false;

  std::vector<art::Ptr<recob::Hit>> hitscl1uniquetpc;
  std::vector<art::Ptr<recob::Hit>> hitscl2uniquetpc;

  size_t tpc = vec1[0]->WireID().TPC;
  for (size_t i = 0; i < vec1.size(); ++i)
    for (size_t j = 0; j < vec2.size(); ++j)
      if ((vec1[i]->WireID().TPC == tpc) && (vec2[j]->WireID().TPC == tpc)) {
        hitscl1uniquetpc.push_back(vec1[i]);
        hitscl2uniquetpc.push_back(vec2[j]);
      }

  if ((hitscl1uniquetpc.size() < 3) || (hitscl2uniquetpc.size() < 3)) return false;

  pma::Track3D* track = fProjectionMatchingAlg.buildSegment(hitscl1uniquetpc, hitscl2uniquetpc);
  for (size_t i = 0; i < input.size(); ++i) {
    std::vector<Hit2D*> hits2dcl = input[i]->GetHits2D();
    for (size_t h = 0; h < hits2dcl.size(); ++h) {
      TVector2 pfront = pma::GetProjectionToPlane(
        track->front()->Point3D(), plane3, track->FrontTPC(), track->FrontCryo());
      TVector2 pback = pma::GetProjectionToPlane(
        track->back()->Point3D(), plane3, track->BackTPC(), track->BackCryo());
      if ((pma::Dist2(hits2dcl[h]->GetPointCm(), pfront) < 1.0F) &&
          (pma::Dist2(hits2dcl[h]->GetPointCm(), pback) < 1.0F)) {
        result = true;
        break;
      }
    }
  }
  delete track;
  return result;
}

bool
ems::EMShower3D::Validate(art::Event const& e, const pma::Track3D& src, size_t plane)
{
  bool result = false;

  art::FindManyP<recob::Hit> fbc(fCluListHandle, e, fCluModuleLabel);
  std::vector<art::Ptr<recob::Hit>> hitscl;
  for (size_t id = 0; id < fClustersNotUsed.size(); id++) {
    std::vector<art::Ptr<recob::Hit>> hits = fbc.at(fClustersNotUsed[id]);
    for (size_t i = 0; i < hits.size(); i++)
      hitscl.push_back(hits[i]);
  }

  if (fProjectionMatchingAlg.validate(src, hitscl) > 0.2) result = true;

  return result;
}

bool
ems::EMShower3D::Has(const std::vector<size_t>& v, size_t idx)
{
  for (auto c : v)
    if (c == idx) return true;
  return false;
}

bool
ems::EMShower3D::GetCloseHits(double r2d,
                              const std::vector<art::Ptr<recob::Hit>>& hits_in,
                              std::vector<size_t>& used,
                              std::vector<art::Ptr<recob::Hit>>& hits_out)
{

  hits_out.clear();

  const double gapMargin = 5.0; // can be changed to f(id_tpc1, id_tpc2)
  size_t idx = 0;

  while ((idx < hits_in.size()) && Has(used, idx))
    idx++;

  if (idx < hits_in.size()) {
    hits_out.push_back(hits_in[idx]);
    used.push_back(idx);

    double r2d2 = r2d * r2d;
    double gapMargin2 = sqrt(2 * gapMargin * gapMargin);
    gapMargin2 = (gapMargin2 + r2d) * (gapMargin2 + r2d);

    bool collect = true;
    while (collect) {
      collect = false;
      for (size_t i = 0; i < hits_in.size(); i++)
        if (!Has(used, i)) {
          art::Ptr<recob::Hit> hi = hits_in[i];
          TVector2 hi_cm = pma::WireDriftToCm(hi->WireID().Wire,
                                              hi->PeakTime(),
                                              hi->WireID().Plane,
                                              hi->WireID().TPC,
                                              hi->WireID().Cryostat);

          bool accept = false;
          //for (auto const& ho : hits_out)
          for (size_t idx_o = 0; idx_o < hits_out.size(); idx_o++) {
            art::Ptr<recob::Hit> ho = hits_out[idx_o];

            double d2 = pma::Dist2(hi_cm,
                                   pma::WireDriftToCm(ho->WireID().Wire,
                                                      ho->PeakTime(),
                                                      ho->WireID().Plane,
                                                      ho->WireID().TPC,
                                                      ho->WireID().Cryostat));

            if (hi->WireID().TPC == ho->WireID().TPC) {
              if (d2 < r2d2) {
                accept = true;
                break;
              }
            }
            else {
              if (d2 < gapMargin2) {
                accept = true;
                break;
              }
            }
          }
          if (accept) {
            collect = true;
            hits_out.push_back(hi);
            used.push_back(i);
          }
        }
    }
    return true;
  }
  else
    return false;
}

void
ems::EMShower3D::FilterOutSmallParts(double r2d,
                                     const std::vector<art::Ptr<recob::Hit>>& hits_in,
                                     std::vector<art::Ptr<recob::Hit>>& hits_out)
{
  size_t min_size = hits_in.size() / 5;
  if (min_size < 3) min_size = 3;

  std::vector<size_t> used;
  std::vector<art::Ptr<recob::Hit>> close_hits;

  while (GetCloseHits(r2d, hits_in, used, close_hits)) {
    if (close_hits.size() > min_size)
      for (auto h : close_hits)
        hits_out.push_back(h);
  }
}
DEFINE_ART_MODULE(ems::EMShower3D)

// Christopher Backhouse - bckhouse@fnal.gov

// Test file at Caltech: nfs/raid11/dunesam/prodgenie_nu_dune10kt_1x2x6_mcc7.0/prodgenie_nu_dune10kt_1x2x6_63_20160811T171439_merged.root

// C/C++ standard libraries
#include <string>
#include <vector>
#include <iostream>

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larsim/MCCheater/BackTracker.h"

#include "TGraph.h"
#include "TPad.h"

#include "Solver.h"
#include "HashTuple.h"

template<class T> T sqr(T x){return x*x;}

namespace reco3d
{

class Reco3D : public art::EDProducer
{
public:

  explicit Reco3D(const fhicl::ParameterSet& pset);
  virtual ~Reco3D();

  void produce(art::Event& evt);
  void beginJob();
  void endJob();

protected:
  void FillSystemToSpacePoints(const std::vector<CollectionWireHit*> cwires,
                               std::vector<recob::SpacePoint>& pts) const;

  bool fFit;

  double fAlpha;

  const detinfo::DetectorProperties* detprop;
};

DEFINE_ART_MODULE(Reco3D)

// ---------------------------------------------------------------------------
Reco3D::Reco3D(const fhicl::ParameterSet& pset)
  : fFit(pset.get<bool>("Fit")),
    fAlpha(pset.get<double>("Alpha"))
{
  produces<std::vector<recob::SpacePoint>>("pre");
  if(fFit){
    produces<std::vector<recob::SpacePoint>>();
    produces<std::vector<recob::SpacePoint>>("noreg");
  }
}

// ---------------------------------------------------------------------------
Reco3D::~Reco3D()
{
}

// ---------------------------------------------------------------------------
void Reco3D::beginJob()
{
  detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
}

// ---------------------------------------------------------------------------
void Reco3D::endJob()
{
}

// ---------------------------------------------------------------------------
// tpc makes sure the intersection is on the side we were expecting
bool ISect(const geo::Geometry* geom, int chanA, int chanB, geo::TPCID tpc,
           geo::WireIDIntersection& pt)
{
  // string has a default implementation for hash. Perhaps TPCID should
  // implement a hash function directly.
  typedef std::tuple<int, int, std::string/*geo::TPCID*/> Key_t;
  static std::unordered_map<Key_t, bool> bCache;
  static std::unordered_map<Key_t, geo::WireIDIntersection> pCache;

  // Prevent cache from growing without bound
  if(bCache.size() > 1e8){
    std::cout << "Clearing Intersection caches" << std::endl;
    bCache.clear();
    pCache.clear();
  }

  const Key_t key = std::make_tuple(chanA, chanB, std::string(tpc));
  if(bCache.count(key)){
    if(bCache[key]) pt = pCache[key];
    return bCache[key];
  }

  const std::vector<geo::WireID> awires = geom->ChannelToWire(chanA);
  const std::vector<geo::WireID> bwires = geom->ChannelToWire(chanB);

  for(geo::WireID awire: awires){
    if(geo::TPCID(awire) != tpc) continue;
    for(geo::WireID bwire: bwires){
      if(geo::TPCID(bwire) != tpc) continue;

      if(geom->WireIDsIntersect(awire, bwire, pt)){
        bCache[key] = true;
        pCache[key] = pt;
        return true;
      }
    }
  }

  bCache[key] = false;
  return false;
}

// ---------------------------------------------------------------------------
bool ISect(const geo::Geometry* geom, int chanA, int chanB, geo::TPCID tpc)
{
  geo::WireIDIntersection junk;
  return ISect(geom, chanA, chanB, tpc, junk);
}

// ---------------------------------------------------------------------------
bool CloseTime(double ta, double tb)
{
  // TODO - figure out cut value. Emprically 5 is a bit small
  return fabs(ta-tb) < 10;
}

// ---------------------------------------------------------------------------
bool CloseSpace(geo::WireIDIntersection ra, geo::WireIDIntersection rb)
{
  TVector3 pa(ra.y, ra.z, 0);
  TVector3 pb(rb.y, rb.z, 0);

  // TODO - figure out cut value. Empirically .25 is a bit small
  return (pa-pb).Mag() < .5;
}

// ---------------------------------------------------------------------------
void FastForward(std::vector<InductionWireHit*>::iterator& it,
                 double target,
                 const std::vector<InductionWireHit*>::const_iterator& end)
{
  while(it != end &&
        (*it)->fTime < target &&
        !CloseTime((*it)->fTime, target)) ++it;
}

// ---------------------------------------------------------------------------
void BuildSystem(std::vector<art::Ptr<recob::Hit>>& xhits,
                 std::vector<art::Ptr<recob::Hit>>& uhits,
                 std::vector<art::Ptr<recob::Hit>>& vhits,
                 std::vector<CollectionWireHit*>& cwires,
                 std::vector<InductionWireHit*>& iwires,
                 bool incNei)
{
  art::ServiceHandle<geo::Geometry> geom;
  const detinfo::DetectorProperties* detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
  // Maps from TPC to the induction wires. Normally want to access them this
  // way.
  std::map<geo::TPCID, std::vector<InductionWireHit*>> uwires, vwires;

  for(auto& ihits: {uhits, vhits}){
      for(auto& hit: ihits){
      const std::vector<geo::TPCID> tpcs = geom->ROPtoTPCs(geom->ChannelToROP(hit->Channel()));

      // TODO: Empirically, total collection charge is about 5% high of total
      // induction charge, which might cause "spare" charge to go where it's
      // not wanted.
      InductionWireHit* iwire = new InductionWireHit(hit->Channel(), hit->PeakTime(), hit->Integral() * .95);
      iwires.emplace_back(iwire);

      assert(tpcs.size() == 2);

      for(geo::TPCID tpc: tpcs){
        if(hit->View() == geo::kU) uwires[tpc].push_back(iwire);
        if(hit->View() == geo::kV) vwires[tpc].push_back(iwire);
      } // end for tpc
    } // end for hit
  } // end for U/V

  std::map<geo::TPCID, std::vector<art::Ptr<recob::Hit>>> xhits_by_tpc;
  for(auto& xhit: xhits){
    const std::vector<geo::TPCID> tpcs = geom->ROPtoTPCs(geom->ChannelToROP(xhit->Channel()));
    assert(tpcs.size() == 1);
    const geo::TPCID tpc = tpcs[0];
    xhits_by_tpc[tpc].push_back(xhit);
  }

  struct UVCrossing
  {
    geo::TPCID tpc;
    InductionWireHit *u, *v;

    bool operator<(const UVCrossing& x) const
    {
      return std::make_tuple(tpc, u, v) < std::make_tuple(x.tpc, x.u, x.v);
    }
  };

  // Build a table of UV crossers up front
  std::cout << "Building UV table..." << std::endl;
  std::map<UVCrossing, bool> isectUV;
  std::map<UVCrossing, geo::WireIDIntersection> ptsUV;

  for(auto it: uwires){
    const geo::TPCID tpc = it.first;

    auto vwires_begin = vwires[tpc].begin();

    for(InductionWireHit* uwire: uwires[tpc]){

      // Fast-forward up to the first vwire that could be relevant
      FastForward(vwires_begin, uwire->fTime, vwires[tpc].end());

      for(auto vit = vwires_begin; vit != vwires[tpc].end(); ++vit){
        InductionWireHit* vwire = *vit;

        // No more vwires can be relevant, bail out
        if(vwire->fTime > uwire->fTime &&
           !CloseTime(uwire->fTime, vwire->fTime)) break;

        const UVCrossing key = {tpc, uwire, vwire};

        isectUV[key] = ISect(geom.get(),
                             uwire->fChannel, vwire->fChannel, tpc,
                             ptsUV[key]);
      } // end for vwire
    } // end for uwire
  } // end for tpc

  std::cout << "Finding XUV coincidences..." << std::endl;
  std::vector<SpaceCharge*> spaceCharges;

  for(auto it: xhits_by_tpc){
    const geo::TPCID tpc = it.first;

    auto uwires_begin = uwires[tpc].begin();
    auto vwires_begin = vwires[tpc].begin();

    for(auto& hit: it.second){

      const std::vector<geo::WireID> ws = geom->ChannelToWire(hit->Channel());
      assert(ws.size() == 1);

      FastForward(uwires_begin, hit->PeakTime(), uwires[tpc].end());
      FastForward(vwires_begin, hit->PeakTime(), vwires[tpc].end());

      // Figure out which vwires intersect this xwire here so we don't do N^2
      // nesting inside the uwire loop below.
      std::vector<InductionWireHit*> vwires_cross;
      std::unordered_map<InductionWireHit*, geo::WireIDIntersection> ptsXV;
      vwires_cross.reserve(vwires[tpc].size()); // avoid reallocations
      for(auto vit = vwires_begin; vit != vwires[tpc].end(); ++vit){
        InductionWireHit* vwire = *vit;

        if(vwire->fTime > hit->PeakTime() &&
           !CloseTime(vwire->fTime, hit->PeakTime())) break;

        if(ISect(geom.get(), hit->Channel(), vwire->fChannel, tpc, ptsXV[vwire]))
          vwires_cross.push_back(vwire);
      } // end for vwire

      std::vector<SpaceCharge*> crossers;
      for(auto uit = uwires_begin; uit != uwires[tpc].end(); ++uit){
        InductionWireHit* uwire = *uit;

        if(uwire->fTime > hit->PeakTime() &&
           !CloseTime(uwire->fTime, hit->PeakTime())) break;

        geo::WireIDIntersection ptXU;
        if(!ISect(geom.get(), hit->Channel(), uwire->fChannel, tpc, ptXU)) continue;

        for(InductionWireHit* vwire: vwires_cross){

          const geo::WireIDIntersection ptXV = ptsXV[vwire];
          if(!CloseSpace(ptXU, ptXV)) continue;

          if(!isectUV[{tpc, uwire, vwire}]) continue;

          const geo::WireIDIntersection ptUV = ptsUV[{tpc, uwire, vwire}];

          if(!CloseSpace(ptXU, ptUV) ||
             !CloseSpace(ptXV, ptUV)) continue;

          // This average aleviates the problem with a single collection wire
          // matching multiple hits on an induction wire at different times.
          const double t = (hit->PeakTime()+uwire->fTime+vwire->fTime)/3;
          const double xpos = detprop->ConvertTicksToX(t, ws[0]);
          // TODO exactly which 3D position to set for this point?
          // Don't have a cwire object yet, set it later
          SpaceCharge* sc = new SpaceCharge(xpos, ptXU.y, ptXU.z,
                                            0, uwire, vwire);
          spaceCharges.push_back(sc);
          crossers.push_back(sc);
        } // end for vwire
      } // end for uwire

      CollectionWireHit* cwire = new CollectionWireHit(hit->Channel(), hit->PeakTime(), hit->Integral(), crossers);
      cwires.push_back(cwire);
      for(SpaceCharge* sc: crossers) sc->fCWire = cwire;

    } // end for hit
  } // end for it (tpc)


  if(incNei){
    static const double kCritDist = 5;

    // Could use a QuadTree or VPTree etc, but seems like overkill
    class IntCoord
    {
    public:
      IntCoord(const SpaceCharge& sc)
        : fX(sc.fX/kCritDist),
          fY(sc.fY/kCritDist),
          fZ(sc.fZ/kCritDist)
      {
      }

      bool operator<(const IntCoord& i) const
      {
        return std::make_tuple(fX, fY, fZ) < std::make_tuple(i.fX, i.fY, i.fZ);
      }

      std::vector<IntCoord> Neighbours() const
      {
        std::vector<IntCoord> ret;
        for(int dx = -1; dx <= +1; ++dx){
          for(int dy = -1; dy <= +1; ++dy){
            for(int dz = -1; dz <= +1; ++dz){
              ret.push_back(IntCoord(fX+dx, fY+dy, fZ+dz));
            }
          }
        }
        return ret;
      }
    protected:
      IntCoord(int x, int y, int z) : fX(x), fY(y), fZ(z) {}

      int fX, fY, fZ;
    };

    std::map<IntCoord, std::vector<SpaceCharge*>> scMap;
    for(SpaceCharge* sc: spaceCharges){
      scMap[IntCoord(*sc)].push_back(sc);
    }

    std::cout << "Neighbour search..." << std::endl;
    std::cout << spaceCharges.size() << std::endl;
    // Now that we know all the space charges, can go through and assign neighbours

    int Ntests = 0;
    int Nnei = 0;
    for(SpaceCharge* sc1: spaceCharges){
      IntCoord ic(*sc1);
      for(IntCoord icn: ic.Neighbours()){
        for(SpaceCharge* sc2: scMap[icn]){

          ++Ntests;

          if(sc1 == sc2) continue;
          /*const*/ double dist2 = sqr(sc1->fX-sc2->fX) + sqr(sc1->fY-sc2->fY) + sqr(sc1->fZ-sc2->fZ);

          if(dist2 > sqr(kCritDist)) continue;

          if(dist2 == 0){
            std::cout << "ZERO DISTANCE SOMEHOW?" << std::endl;
            std::cout << sc1->fCWire << " " << sc1->fWire1 << " " << sc1->fWire2 << std::endl;
            std::cout << sc2->fCWire << " " << sc2->fWire1 << " " << sc2->fWire2 << std::endl;
            std::cout << dist2 << " " << sc1->fX << " " << sc2->fX << " " << sc1->fY << " " << sc2->fY << " " << sc1->fZ << " " << sc2->fZ << std::endl;
            continue;
            dist2 = sqr(kCritDist);
          }

          ++Nnei;

          // This is a pretty random guess
          const double coupling = exp(-sqrt(dist2)/2);
          sc1->fNeighbours.emplace_back(sc2, coupling);

          if(isnan(1/sqrt(dist2)) || isinf(1/sqrt(dist2))){
            std::cout << dist2 << " " << sc1->fX << " " << sc2->fX << " " << sc1->fY << " " << sc2->fY << " " << sc1->fZ << " " << sc2->fZ << std::endl;
            abort();
          }
        } // end for sc2
      } // end for icn

      // The neighbours lists use the most memory, so be careful to trim
      sc1->fNeighbours.shrink_to_fit();
    } // end for sc1

    for(SpaceCharge* sc: spaceCharges){
      for(Neighbour& nei: sc->fNeighbours){
        sc->fNeiPotential += nei.fCoupling * nei.fSC->fPred;
      }
    }

    std::cout << Ntests << " tests to find " << Nnei << std::endl;
  }
}

// ---------------------------------------------------------------------------
void Reco3D::FillSystemToSpacePoints(const std::vector<CollectionWireHit*> cwires,
                                     std::vector<recob::SpacePoint>& pts) const
{
  const double err[6] = {0,};

  for(const CollectionWireHit* cwire: cwires){
    for(const SpaceCharge* sc: cwire->fCrossings){
      if(sc->fPred == 0) continue;

      // TODO find somewhere to save the charge too
      const double xyz[3] = {sc->fX, sc->fY, sc->fZ};
      pts.emplace_back(xyz, err, 0);
    }
  }
}

// ---------------------------------------------------------------------------
void Reco3D::produce(art::Event& evt)
{
  art::Handle<std::vector<recob::Hit>> hits;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel("gaushit", hits))
    art::fill_ptr_vector(hitlist, hits);

  // Skip very small events
  if(hits->size() < 20){
    auto spcol_pre = std::make_unique<std::vector<recob::SpacePoint>>();
    evt.put(std::move(spcol_pre), "pre");
    if(fFit){
      auto spcol = std::make_unique<std::vector<recob::SpacePoint>>();
      evt.put(std::move(spcol));
      auto spcol_noreg = std::make_unique<std::vector<recob::SpacePoint>>();
      evt.put(std::move(spcol_noreg), "noreg");
    }
    return;
  }

  art::ServiceHandle<geo::Geometry> geom;

  std::vector<art::Ptr<recob::Hit>> xhits, uhits, vhits;
  for(auto & hit: hitlist){
    if(hit->SignalType() == geo::kCollection){
      xhits.push_back(hit);
    }
    else{
      if(hit->View() == geo::kU) uhits.push_back(hit);
      if(hit->View() == geo::kV) vhits.push_back(hit);
    }
  } // end for hit

  // BuildSystem requires the hits to be sorted in time
  for(auto v: {&xhits, &uhits, &vhits}){
    std::sort(v->begin(), v->end(),
              [](art::Ptr<recob::Hit>& a, art::Ptr<recob::Hit>& b)
              {
                return a->PeakTime() < b->PeakTime();
              });
  }

  std::vector<CollectionWireHit*> cwires;
  // So we can find them all to free the memory
  std::vector<InductionWireHit*> iwires;

  BuildSystem(xhits, uhits, vhits, cwires, iwires, fAlpha != 0);

  auto spcol_pre = std::make_unique<std::vector<recob::SpacePoint>>();
  FillSystemToSpacePoints(cwires, *spcol_pre);
  evt.put(std::move(spcol_pre), "pre");

  if(fFit){
    std::cout << "Iterating..." << std::endl;
    double prevMetric = Metric(cwires, 0);//fAlpha);
    std::cout << "Begin: " << prevMetric << std::endl;
    for(int i = 0;; ++i){
      Iterate(cwires, 0);//fAlpha);
      const double metric = Metric(cwires, 0);//fAlpha);
      std::cout << i << " " << metric << std::endl;
      if(fabs(metric-prevMetric) < 1e-3*fabs(prevMetric)) break;
      if(i > 100) break;
      //    if(metric/prevMetric > .9999) break;
      prevMetric = metric;
    }

    auto spcol_noreg = std::make_unique<std::vector<recob::SpacePoint>>();
    FillSystemToSpacePoints(cwires, *spcol_noreg);
    evt.put(std::move(spcol_noreg), "noreg");

    prevMetric = Metric(cwires, fAlpha);
    std::cout << "Begin: " << prevMetric << std::endl;
    for(int i = 0;; ++i){
      Iterate(cwires, fAlpha);
      const double metric = Metric(cwires, fAlpha);
      std::cout << i << " " << metric << std::endl;
      if(fabs(metric-prevMetric) < 1e-3*fabs(prevMetric)) break;
      if(i > 100) break;
      //    if(metric/prevMetric > .9999) break;
      prevMetric = metric;
    }

    auto spcol = std::make_unique<std::vector<recob::SpacePoint>>();
    FillSystemToSpacePoints(cwires, *spcol);
    evt.put(std::move(spcol));
  } // end if fFit

  for(InductionWireHit* i: iwires) delete i;
  for(CollectionWireHit* c: cwires) delete c;
}

} // end namespace reco3d

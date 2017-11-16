// Christopher Backhouse - bckhouse@fnal.gov

// Test file at Caltech: /nfs/raid11/dunesam/prodgenie_nu_dune10kt_1x2x6_mcc7.0/prodgenie_nu_dune10kt_1x2x6_63_20160811T171439_merged.root

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
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larsim/MCCheater/BackTracker.h"

#include "TGraph.h"
#include "TH1.h"
#include "TPad.h"

#include "Solver.h"
#include "HashTuple.h"

template<class T> T sqr(T x){return x*x;}

namespace reco3d
{

/// The position in the drift direction that an induction wire's time implies
/// depends on which TPC we assume the hit originated in.
struct InductionWireWithXPos
{
  InductionWireWithXPos(InductionWireHit* w, double x) : iwire(w), xpos(x) {}

  bool operator<(const InductionWireWithXPos& w) const {return xpos < w.xpos;}

  InductionWireHit* iwire;
  double xpos;
};

class SpacePointSolver : public art::EDProducer
{
public:

  explicit SpacePointSolver(const fhicl::ParameterSet& pset);
  virtual ~SpacePointSolver();

  void produce(art::Event& evt);
  void beginJob();
  void endJob();

protected:
  double HitToXPos(const recob::Hit& hit, geo::TPCID tpc) const;

  bool CloseDrift(double xa, double xb) const;

  bool CloseSpace(geo::WireIDIntersection ra,
                  geo::WireIDIntersection rb) const;

  void FastForward(std::vector<InductionWireWithXPos>::iterator& it,
                   double target,
                   const std::vector<InductionWireWithXPos>::const_iterator& end) const;

  bool ISect(int chanA, int chanB, geo::TPCID tpc,
             geo::WireIDIntersection& pt) const;

  bool ISect(int chanA, int chanB, geo::TPCID tpc) const;

  void AddNeighbours(const std::vector<SpaceCharge*>& spaceCharges) const;

  typedef std::map<const WireHit*, art::Ptr<recob::Hit>> HitMap_t;

  void BuildSystemSP(const std::vector<art::Ptr<recob::Hit>>& xhits,
                     const std::vector<art::Ptr<recob::Hit>>& uhits,
                     const std::vector<art::Ptr<recob::Hit>>& vhits,
                     std::vector<CollectionWireHit*>& cwires,
                     std::vector<InductionWireHit*>& iwires,
                     bool incNei,
                     HitMap_t& hitmap) const;

  void BuildSystemDP(const std::vector<art::Ptr<recob::Hit>>& xhits,
                     const std::vector<art::Ptr<recob::Hit>>& uhits,
                     std::vector<CollectionWireHit*>& cwires,
                     std::vector<InductionWireHit*>& iwires,
                     bool incNei,
                     HitMap_t& hitmap) const;

  void FillSystemToSpacePoints(const std::vector<CollectionWireHit*> cwires,
                               std::vector<recob::SpacePoint>& pts) const;

  void FillAssns(art::Event& evt,
                 const std::vector<CollectionWireHit*> cwires,
                 const std::vector<recob::SpacePoint>& pts,
                 art::Assns<recob::Hit, recob::SpacePoint>& assn,
                 const HitMap_t& hitmap) const;

  std::string fHitLabel;

  bool fFit;

  double fAlpha;

  TH1* fDeltaX;

  const detinfo::DetectorProperties* detprop;
  const geo::GeometryCore* geom;
};

DEFINE_ART_MODULE(SpacePointSolver)

// ---------------------------------------------------------------------------
SpacePointSolver::SpacePointSolver(const fhicl::ParameterSet& pset)
  : fHitLabel(pset.get<std::string>("HitLabel")),
    fFit(pset.get<bool>("Fit")),
    fAlpha(pset.get<double>("Alpha"))
{
  produces<std::vector<recob::SpacePoint>>("pre");
  if(fFit){
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();
    produces<std::vector<recob::SpacePoint>>("noreg");
  }
}

// ---------------------------------------------------------------------------
SpacePointSolver::~SpacePointSolver()
{
}

// ---------------------------------------------------------------------------
void SpacePointSolver::beginJob()
{
  detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
  geom = art::ServiceHandle<geo::Geometry>()->provider();

  art::ServiceHandle<art::TFileService> tfs;
  fDeltaX = tfs->make<TH1F>("deltax", ";#Deltax (cm)", 100, -2, +2);
}

// ---------------------------------------------------------------------------
void SpacePointSolver::endJob()
{
}

// ---------------------------------------------------------------------------
// tpc makes sure the intersection is on the side we were expecting
bool SpacePointSolver::ISect(int chanA, int chanB, geo::TPCID tpc,
                             geo::WireIDIntersection& pt) const
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
bool SpacePointSolver::ISect(int chanA, int chanB, geo::TPCID tpc) const
{
  geo::WireIDIntersection junk;
  return ISect(chanA, chanB, tpc, junk);
}

// ---------------------------------------------------------------------------
bool SpacePointSolver::CloseDrift(double xa, double xb) const
{
  // Used to cut at 10 ticks (for basically empirical reasons). Reproduce that
  // in x.
  // Sampling rate is in ns/ticks
  // Drift velocity is in cm/us
  //  static const double k = 10*detprop->SamplingRate()*1e-3*detprop->DriftVelocity();
  const double k = 0.2;//0.4; // 1 sigma on deltax plot

  // TODO - figure out cut value
  return fabs(xa-xb) < k;
}

// ---------------------------------------------------------------------------
bool SpacePointSolver::CloseSpace(geo::WireIDIntersection ra,
                                  geo::WireIDIntersection rb) const
{
  TVector3 pa(ra.y, ra.z, 0);
  TVector3 pb(rb.y, rb.z, 0);

  // TODO - figure out cut value. Empirically .25 is a bit small
  //  return (pa-pb).Mag() < .5;

  return (pa-pb).Mag() < .35;
}

// ---------------------------------------------------------------------------
void SpacePointSolver::
FastForward(std::vector<InductionWireWithXPos>::iterator& it,
            double target,
            const std::vector<InductionWireWithXPos>::const_iterator& end) const
{
  while(it != end && it->xpos < target && !CloseDrift(it->xpos, target)) ++it;
}

// ---------------------------------------------------------------------------
double SpacePointSolver::HitToXPos(const recob::Hit& hit, geo::TPCID tpc) const
{
  const std::vector<geo::WireID> ws = geom->ChannelToWire(hit.Channel());
  for(geo::WireID w: ws){
    if(geo::TPCID(w) == tpc)
      return detprop->ConvertTicksToX(hit.PeakTime(), w);
  }

  std::cout << "Wire does not exist on given TPC!" << std::endl;
  abort();
}

// ---------------------------------------------------------------------------
void SpacePointSolver::
AddNeighbours(const std::vector<SpaceCharge*>& spaceCharges) const
{
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

// ---------------------------------------------------------------------------
void SpacePointSolver::
BuildSystemSP(const std::vector<art::Ptr<recob::Hit>>& xhits,
              const std::vector<art::Ptr<recob::Hit>>& uhits,
              const std::vector<art::Ptr<recob::Hit>>& vhits,
              std::vector<CollectionWireHit*>& cwires,
              std::vector<InductionWireHit*>& iwires,
              bool incNei,
              HitMap_t& hitmap) const
{
  std::map<geo::TPCID, std::vector<art::Ptr<recob::Hit>>> xhits_by_tpc;
  for(auto& xhit: xhits){
    const std::vector<geo::TPCID> tpcs = geom->ROPtoTPCs(geom->ChannelToROP(xhit->Channel()));
    assert(tpcs.size() == 1);
    const geo::TPCID tpc = tpcs[0];
    xhits_by_tpc[tpc].push_back(xhit);
  }

  // Maps from TPC to the induction wires. Normally want to access them this
  // way.
  std::map<geo::TPCID, std::vector<InductionWireWithXPos>> uwires, vwires;

  for(auto& ihits: {uhits, vhits}){
      for(auto& hit: ihits){
      const std::vector<geo::TPCID> tpcs = geom->ROPtoTPCs(geom->ChannelToROP(hit->Channel()));

      // TODO: Empirically, total collection charge is about 5% high of total
      // induction charge, which might cause "spare" charge to go where it's
      // not wanted.
      InductionWireHit* iwire = new InductionWireHit(hit->Channel(), hit->Integral() * .95);
      iwires.emplace_back(iwire);

      for(geo::TPCID tpc: tpcs){
        if(xhits_by_tpc.count(tpc) == 0) continue;

        const double xpos = HitToXPos(*hit, tpc);

        if(hit->View() == geo::kU) uwires[tpc].emplace_back(iwire, xpos);
        if(hit->View() == geo::kV) vwires[tpc].emplace_back(iwire, xpos);
      } // end for tpc
    } // end for hit
  } // end for U/V

  for(auto it = uwires.begin(); it != uwires.end(); ++it){
    std::sort(it->second.begin(), it->second.end());
  }
  for(auto it = vwires.begin(); it != vwires.end(); ++it){
    std::sort(it->second.begin(), it->second.end());
  }

  for(auto it = xhits_by_tpc.begin(); it != xhits_by_tpc.end(); ++it){
    const geo::TPCID tpc = it->first;
    std::sort(it->second.begin(), it->second.end(),
              [this, tpc](art::Ptr<recob::Hit>& a, art::Ptr<recob::Hit>& b)
              {
                return HitToXPos(*a, tpc) < HitToXPos(*b, tpc);
              });
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

    for(InductionWireWithXPos uwire: uwires[tpc]){

      // Fast-forward up to the first vwire that could be relevant
      FastForward(vwires_begin, uwire.xpos, vwires[tpc].end());

      for(auto vit = vwires_begin; vit != vwires[tpc].end(); ++vit){
        const InductionWireWithXPos vwire = *vit;

        // No more vwires can be relevant, bail out
        if(vwire.xpos > uwire.xpos &&
           !CloseDrift(uwire.xpos, vwire.xpos)) break;

        const UVCrossing key = {tpc, uwire.iwire, vwire.iwire};

        isectUV[key] = ISect(uwire.iwire->fChannel, vwire.iwire->fChannel, tpc,
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
      const double xpos = HitToXPos(*hit, tpc);

      FastForward(uwires_begin, xpos, uwires[tpc].end());
      FastForward(vwires_begin, xpos, vwires[tpc].end());

      // Figure out which vwires intersect this xwire here so we don't do N^2
      // nesting inside the uwire loop below.
      std::vector<InductionWireWithXPos> vwires_cross;
      std::unordered_map<InductionWireHit*, geo::WireIDIntersection> ptsXV;
      vwires_cross.reserve(vwires[tpc].size()); // avoid reallocations
      for(auto vit = vwires_begin; vit != vwires[tpc].end(); ++vit){
        InductionWireWithXPos vwire = *vit;

        if(vwire.xpos > xpos && !CloseDrift(vwire.xpos, xpos)) break;

        if(ISect(hit->Channel(), vwire.iwire->fChannel, tpc, ptsXV[vwire.iwire]))
          vwires_cross.push_back(vwire);
      } // end for vwire

      std::vector<SpaceCharge*> crossers;
      for(auto uit = uwires_begin; uit != uwires[tpc].end(); ++uit){
        const InductionWireWithXPos uwire = *uit;

        if(uwire.xpos > xpos && !CloseDrift(uwire.xpos, xpos)) break;

        geo::WireIDIntersection ptXU;
        if(!ISect(hit->Channel(), uwire.iwire->fChannel, tpc, ptXU)) continue;

        for(const InductionWireWithXPos& vwire: vwires_cross){

          const geo::WireIDIntersection ptXV = ptsXV[vwire.iwire];
          if(!CloseSpace(ptXU, ptXV)) continue;

          if(!isectUV[{tpc, uwire.iwire, vwire.iwire}]) continue;

          const geo::WireIDIntersection ptUV = ptsUV[{tpc, uwire.iwire, vwire.iwire}];

          if(!CloseSpace(ptXU, ptUV) ||
             !CloseSpace(ptXV, ptUV)) continue;

          fDeltaX->Fill(xpos-uwire.xpos);
          fDeltaX->Fill(xpos-vwire.xpos);
          fDeltaX->Fill(uwire.xpos-vwire.xpos);

          // TODO exactly which 3D position to set for this point? This average
          // aleviates the problem with a single collection wire matching
          // multiple hits on an induction wire at different times.

          // Don't have a cwire object yet, set it later
          SpaceCharge* sc = new SpaceCharge((xpos+uwire.xpos+vwire.xpos)/3,
                                            ptXU.y, ptXU.z,
                                            0, uwire.iwire, vwire.iwire);
          spaceCharges.push_back(sc);
          crossers.push_back(sc);
          hitmap[uwire.iwire] = hit;
          hitmap[vwire.iwire] = hit;
        } // end for vwire
      } // end for uwire

      CollectionWireHit* cwire = new CollectionWireHit(hit->Channel(), hit->Integral(), crossers);
      hitmap[cwire] = hit;
      cwires.push_back(cwire);
      for(SpaceCharge* sc: crossers) sc->fCWire = cwire;
    } // end for hit
  } // end for it (tpc)

  if(incNei) AddNeighbours(spaceCharges);
}

// ---------------------------------------------------------------------------
void SpacePointSolver::
BuildSystemDP(const std::vector<art::Ptr<recob::Hit>>& xhits,
              const std::vector<art::Ptr<recob::Hit>>& uhits,
              std::vector<CollectionWireHit*>& cwires,
              std::vector<InductionWireHit*>& iwires,
              bool incNei,
              HitMap_t& hitmap) const
{
  std::map<geo::TPCID, std::vector<art::Ptr<recob::Hit>>> xhits_by_tpc;
  for(auto& xhit: xhits){
    const std::vector<geo::TPCID> tpcs = geom->ROPtoTPCs(geom->ChannelToROP(xhit->Channel()));
    assert(tpcs.size() == 1);
    const geo::TPCID tpc = tpcs[0];
    xhits_by_tpc[tpc].push_back(xhit);
  }

  // Maps from TPC to the induction wires. Normally want to access them this
  // way.
  std::map<geo::TPCID, std::vector<InductionWireWithXPos>> uwires;

  for(auto& hit: uhits){
    const std::vector<geo::TPCID> tpcs = geom->ROPtoTPCs(geom->ChannelToROP(hit->Channel()));

    InductionWireHit* iwire = new InductionWireHit(hit->Channel(), hit->Integral());
    iwires.emplace_back(iwire);

    for(geo::TPCID tpc: tpcs){
      if(xhits_by_tpc.count(tpc) == 0) continue;

      const double xpos = HitToXPos(*hit, tpc);

      uwires[tpc].emplace_back(iwire, xpos);
    } // end for tpc
  } // end for hit

  for(auto it = uwires.begin(); it != uwires.end(); ++it){
    std::sort(it->second.begin(), it->second.end());
  }

  for(auto it = xhits_by_tpc.begin(); it != xhits_by_tpc.end(); ++it){
    const geo::TPCID tpc = it->first;
    std::sort(it->second.begin(), it->second.end(),
              [this, tpc](art::Ptr<recob::Hit>& a, art::Ptr<recob::Hit>& b)
              {
                return HitToXPos(*a, tpc) < HitToXPos(*b, tpc);
              });
  }


  std::cout << "Finding UV coincidences..." << std::endl;
  std::vector<SpaceCharge*> spaceCharges;

  for(auto it: xhits_by_tpc){
    const geo::TPCID tpc = it.first;

    auto uwires_begin = uwires[tpc].begin();

    for(auto& hit: it.second){
      const double xpos = HitToXPos(*hit, tpc);

      std::vector<SpaceCharge*> crossers;

      FastForward(uwires_begin, xpos, uwires[tpc].end());

      // Figure out which uwires intersect this xwire here.
      for(auto uit = uwires_begin; uit != uwires[tpc].end(); ++uit){
        InductionWireWithXPos uwire = *uit;

        if(uwire.xpos > xpos && !CloseDrift(uwire.xpos, xpos)) break;

        geo::WireIDIntersection ptXU;
        if(ISect(hit->Channel(), uwire.iwire->fChannel, tpc, ptXU)){

          // Don't have a cwire object yet, set it later
          SpaceCharge* sc = new SpaceCharge((xpos+uwire.xpos)/2, ptXU.y, ptXU.z,
                                            0, 0, uwire.iwire);
          spaceCharges.push_back(sc);
          crossers.push_back(sc);
          hitmap[uwire.iwire] = hit;
        } // end for uwire
      } // end for uit

      CollectionWireHit* cwire = new CollectionWireHit(hit->Channel(), hit->Integral(), crossers);
      hitmap[cwire] = hit;
      cwires.push_back(cwire);
      for(SpaceCharge* sc: crossers) sc->fCWire = cwire;
    } // end for hit
  } // end for it (tpc)

  if(incNei) AddNeighbours(spaceCharges);
}

// ---------------------------------------------------------------------------
void SpacePointSolver::
FillSystemToSpacePoints(const std::vector<CollectionWireHit*> cwires,
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
void SpacePointSolver::
FillAssns(art::Event& evt,
          const std::vector<CollectionWireHit*> cwires,
          const std::vector<recob::SpacePoint>& pts,
          art::Assns<recob::Hit, recob::SpacePoint>& assn,
          const HitMap_t& hitmap) const
{
  unsigned int ptidx = 0;

  // Must follow FillSystemToSpacePoints()'s looping order here
  for(const CollectionWireHit* cwire: cwires){
    for(const SpaceCharge* sc: cwire->fCrossings){
      if(sc->fPred == 0) continue;

      auto it = hitmap.find(cwire);
      assert(it != hitmap.end());
      util::CreateAssn(*this, evt, pts, it->second, assn, "", ptidx);

      if(sc->fWire1){
        auto it1 = hitmap.find(sc->fWire1);
        assert(it1);
        util::CreateAssn(*this, evt, pts, it1->second, assn, "", ptidx);
      }
      if(sc->fWire2){
        auto it2 = hitmap.find(sc->fWire2);
        assert(it2);
        util::CreateAssn(*this, evt, pts, it2->second, assn, "", ptidx);
      }

      ++ptidx;
    }
  }

  if(ptidx != pts.size()){
    std::cout << "Didn't manage to use up all the pts when making Assns!" << std::endl;
    abort();
  }
}

// ---------------------------------------------------------------------------
void SpacePointSolver::produce(art::Event& evt)
{
  art::Handle<std::vector<recob::Hit>> hits;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel(fHitLabel, hits))
    art::fill_ptr_vector(hitlist, hits);

  // Skip very small events
  if(hits->size() < 20){
    auto spcol_pre = std::make_unique<std::vector<recob::SpacePoint>>();
    evt.put(std::move(spcol_pre), "pre");
    if(fFit){
      auto spcol = std::make_unique<std::vector<recob::SpacePoint>>();
      evt.put(std::move(spcol));
      auto assns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
      evt.put(std::move(assns));
      auto spcol_noreg = std::make_unique<std::vector<recob::SpacePoint>>();
      evt.put(std::move(spcol_noreg), "noreg");
    }
    return;
  }

  art::ServiceHandle<geo::Geometry> geom;

  bool isDP = false;
  std::vector<art::Ptr<recob::Hit>> xhits, uhits, vhits;
  for(auto& hit: hitlist){
    if(isnan(hit->Integral()) || isinf(hit->Integral())){
      std::cout << "WARNING: bad recob::Hit::Integral() = "
                << hit->Integral()
                << ". Skipping." << std::endl;
      continue;
    }

    if(hit->SignalType() == geo::kCollection){
      // For DualPhase, both view are collection. Arbitrarily map V to the main
      // "X" view and keep U as-is. For Argoneut and Lariat, collection=V is
      // also the right convention.
      if(hit->View() == geo::kZ || hit->View() == geo::kV){
        xhits.push_back(hit);
      }
      else{
        uhits.push_back(hit);
        isDP = true;
      }
    }
    else{
      if(hit->View() == geo::kU) uhits.push_back(hit);
      if(hit->View() == geo::kV) vhits.push_back(hit);
    }
  } // end for hit

  std::vector<CollectionWireHit*> cwires;
  // So we can find them all to free the memory
  std::vector<InductionWireHit*> iwires;

  HitMap_t hitmap;
  if(isDP)
    BuildSystemDP(xhits, uhits, cwires, iwires, fAlpha != 0, hitmap);
  else
    BuildSystemSP(xhits, uhits, vhits, cwires, iwires, fAlpha != 0, hitmap);

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
    auto assns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
    FillSystemToSpacePoints(cwires, *spcol);
    FillAssns(evt, cwires, *spcol, *assns, hitmap);
    evt.put(std::move(spcol));
    evt.put(std::move(assns));
  } // end if fFit

  for(InductionWireHit* i: iwires) delete i;
  for(CollectionWireHit* c: cwires) delete c;
}

} // end namespace reco3d

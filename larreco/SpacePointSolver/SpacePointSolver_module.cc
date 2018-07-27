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
#include "lardata/ArtDataHelper/ChargedSpacePointCreator.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

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

  void BuildSystemXUV(const std::vector<art::Ptr<recob::Hit>>& xhits,
                      const std::vector<art::Ptr<recob::Hit>>& uhits,
                      const std::vector<art::Ptr<recob::Hit>>& vhits,
                      const std::vector<raw::ChannelID_t>& ubadchans,
                      const std::vector<raw::ChannelID_t>& vbadchans,
                      std::vector<CollectionWireHit*>& cwires,
                      std::vector<InductionWireHit*>& iwires,
                      bool incNei,
                      HitMap_t& hitmap) const;

  void BuildSystemXU(const std::vector<art::Ptr<recob::Hit>>& xhits,
                     const std::vector<art::Ptr<recob::Hit>>& uhits,
                     std::vector<CollectionWireHit*>& cwires,
                     std::vector<InductionWireHit*>& iwires,
                     bool incNei,
                     HitMap_t& hitmap) const;

  /// return whether the point was inserted (only happens when it has charge)
  bool AddSpacePoint(const SpaceCharge& sc,
                     int id,
                     recob::ChargedSpacePointCollectionCreator& points) const;

  void FillSystemToSpacePoints(const std::vector<CollectionWireHit*>& cwires,
                               recob::ChargedSpacePointCollectionCreator& pts) const;

  void FillSystemToSpacePointsAndAssns(const std::vector<CollectionWireHit*>& cwires,
                                       const HitMap_t& hitmap,
                                       recob::ChargedSpacePointCollectionCreator& points,
                                       art::Assns<recob::SpacePoint, recob::Hit>& assn) const;

  std::string fHitLabel;

  bool fFit;
  bool fAllowBadInductionHit;

  double fAlpha;

  double fDistThresh;
  double fDistThreshDrift;

  TH1* fDeltaX;

  const detinfo::DetectorProperties* detprop;
  const geo::GeometryCore* geom;
};

DEFINE_ART_MODULE(SpacePointSolver)

// ---------------------------------------------------------------------------
SpacePointSolver::SpacePointSolver(const fhicl::ParameterSet& pset)
  : fHitLabel(pset.get<std::string>("HitLabel")),
    fFit(pset.get<bool>("Fit")),
    fAllowBadInductionHit(pset.get<bool>("AllowBadInductionHit")),
    fAlpha(pset.get<double>("Alpha")),
    fDistThresh(pset.get<double>("WireIntersectThreshold")),
    fDistThreshDrift(pset.get<double>("WireIntersectThresholdDriftDir"))
{
  recob::ChargedSpacePointCollectionCreator::produces(*this, "pre");
  if(fFit){
    recob::ChargedSpacePointCollectionCreator::produces(*this);
    produces<art::Assns<recob::SpacePoint, recob::Hit>>();
    recob::ChargedSpacePointCollectionCreator::produces(*this, "noreg");
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
  return fabs(xa-xb) < fDistThreshDrift;
}

// ---------------------------------------------------------------------------
bool SpacePointSolver::CloseSpace(geo::WireIDIntersection ra,
                                  geo::WireIDIntersection rb) const
{
  const TVector3 pa(ra.y, ra.z, 0);
  const TVector3 pb(rb.y, rb.z, 0);

  return (pa-pb).Mag() < fDistThresh;
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
BuildSystemXUV(const std::vector<art::Ptr<recob::Hit>>& xhits,
               const std::vector<art::Ptr<recob::Hit>>& uhits,
               const std::vector<art::Ptr<recob::Hit>>& vhits,
               const std::vector<raw::ChannelID_t>& ubadchans,
               const std::vector<raw::ChannelID_t>& vbadchans,
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

  std::map<geo::TPCID, std::vector<raw::ChannelID_t>> ubad_by_tpc, vbad_by_tpc;
  for(raw::ChannelID_t chan: ubadchans){
    for(geo::TPCID tpc: geom->ROPtoTPCs(geom->ChannelToROP(chan))){
      ubad_by_tpc[tpc].push_back(chan);
    }
  }
  for(raw::ChannelID_t chan: vbadchans){
    for(geo::TPCID tpc: geom->ROPtoTPCs(geom->ChannelToROP(chan))){
      vbad_by_tpc[tpc].push_back(chan);
    }
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
      hitmap[iwire] = hit;

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
    //    raw::ChannelID_t uchan, vchan;
    int uchan, vchan;

    bool operator<(const UVCrossing& x) const
    {
      return std::make_tuple(tpc, uchan, vchan) < std::make_tuple(x.tpc, x.uchan, x.vchan);
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

        const UVCrossing key = {tpc, uwire.iwire->fChannel, vwire.iwire->fChannel};

        isectUV[key] = ISect(uwire.iwire->fChannel, vwire.iwire->fChannel, tpc,
                             ptsUV[key]);
      } // end for vwire

      for(raw::ChannelID_t vbad: vbad_by_tpc[tpc]){
        const UVCrossing key = {tpc, uwire.iwire->fChannel, int(vbad)};
        isectUV[key] = ISect(uwire.iwire->fChannel, vbad, tpc, ptsUV[key]);
      }
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
      std::unordered_map</*raw::ChannelID_t*/int, geo::WireIDIntersection> ptsXV;
      vwires_cross.reserve(vwires[tpc].size()); // avoid reallocations
      for(auto vit = vwires_begin; vit != vwires[tpc].end(); ++vit){
        InductionWireWithXPos vwire = *vit;

        if(vwire.xpos > xpos && !CloseDrift(vwire.xpos, xpos)) break;

        if(ISect(hit->Channel(), vwire.iwire->fChannel, tpc, ptsXV[vwire.iwire->fChannel]))
          vwires_cross.push_back(vwire);
      } // end for vwire

      // Figure out which bad V channels cross this X wire.
      std::vector<raw::ChannelID_t> vbad_cross;
      for(raw::ChannelID_t vbad: vbad_by_tpc[tpc]){
        if(ISect(hit->Channel(), vbad, tpc, ptsXV[vbad]))
          vbad_cross.push_back(vbad);
      }

      std::vector<SpaceCharge*> crossers;
      for(auto uit = uwires_begin; uit != uwires[tpc].end(); ++uit){
        const InductionWireWithXPos uwire = *uit;

        if(uwire.xpos > xpos && !CloseDrift(uwire.xpos, xpos)) break;

        geo::WireIDIntersection ptXU;
        if(!ISect(hit->Channel(), uwire.iwire->fChannel, tpc, ptXU)) continue;

        for(const InductionWireWithXPos& vwire: vwires_cross){

          const geo::WireIDIntersection ptXV = ptsXV[vwire.iwire->fChannel];
          if(!CloseSpace(ptXU, ptXV)) continue;

          const UVCrossing key = {tpc, uwire.iwire->fChannel, vwire.iwire->fChannel};
          if(!isectUV[key]) continue;

          const geo::WireIDIntersection ptUV = ptsUV[key];

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
                                            (ptXU.y+ptXV.y+ptUV.y)/3,
                                            (ptXU.z+ptXV.z+ptUV.z)/3,
                                            0, uwire.iwire, vwire.iwire);
          spaceCharges.push_back(sc);
          crossers.push_back(sc);
        } // end for vwire

        for(raw::ChannelID_t vbad: vbad_cross){
          const UVCrossing key = {tpc, uwire.iwire->fChannel, int(vbad)};
          if(!isectUV[key]) continue;
          const geo::WireIDIntersection ptUV = ptsUV[key];
          if(!CloseSpace(ptXU, ptUV)) continue;

          SpaceCharge* sc = new SpaceCharge((xpos+uwire.xpos)/2,
                                            (ptXU.y+ptUV.y)/2,
                                            (ptXU.z+ptUV.z)/2,
                                            0, uwire.iwire, 0);
          spaceCharges.push_back(sc);
          crossers.push_back(sc);
        }
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
BuildSystemXU(const std::vector<art::Ptr<recob::Hit>>& xhits,
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
    hitmap[iwire] = hit;

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
bool SpacePointSolver::
AddSpacePoint(const SpaceCharge& sc,
              int id,
              recob::ChargedSpacePointCollectionCreator& points) const
{
  static const double err[6] = {0,};

  const float charge = sc.fPred;
  if(charge == 0) return false;

  const double xyz[3] = {sc.fX, sc.fY, sc.fZ};
  points.add({ xyz, err, 0.0, id }, charge);

  return true;
}

// ---------------------------------------------------------------------------
void SpacePointSolver::
FillSystemToSpacePoints(const std::vector<CollectionWireHit*>& cwires,
                        recob::ChargedSpacePointCollectionCreator& points) const
{
  int iPoint = 0;
  for(const CollectionWireHit* cwire: cwires){
    for(const SpaceCharge* sc: cwire->fCrossings){
      AddSpacePoint(*sc, iPoint++, points);
    } // for sc
  } // for cwire
}


// ---------------------------------------------------------------------------
void SpacePointSolver::
FillSystemToSpacePointsAndAssns(const std::vector<CollectionWireHit*>& cwires,
                                const HitMap_t& hitmap,
                                recob::ChargedSpacePointCollectionCreator& points,
                                art::Assns<recob::SpacePoint, recob::Hit>& assn) const
{
  int iPoint = 0;
  for(const CollectionWireHit* cwire: cwires){
    for(const SpaceCharge* sc: cwire->fCrossings){
      // fill the space point and reconstructed charge information;
      // if the point is filtered out, it's not inserted (no association either)
      if(!AddSpacePoint(*sc, iPoint++, points)) continue;

      // now fill the associations to the last added space point
      const auto& spsPtr = points.lastSpacePointPtr();

      const auto& hit = hitmap.at(cwire);
      assn.addSingle(spsPtr, hit);

      if(sc->fWire1){
        assn.addSingle(spsPtr, hitmap.at(sc->fWire1));
      }
      if(sc->fWire2){
        assn.addSingle(spsPtr, hitmap.at(sc->fWire2));
      }
    } // for sc
  } // for cwire
}


// ---------------------------------------------------------------------------
void SpacePointSolver::produce(art::Event& evt)
{
  art::Handle<std::vector<recob::Hit>> hits;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel(fHitLabel, hits))
    art::fill_ptr_vector(hitlist, hits);

  recob::ChargedSpacePointCollectionCreator spcol_pre(evt, *this, "pre");
  recob::ChargedSpacePointCollectionCreator spcol_noreg(evt, *this, "noreg");
  recob::ChargedSpacePointCollectionCreator spcol(evt, *this);
  auto assns = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();

  // Skip very small events
  if(hits->size() < 20){
    spcol_pre.put();
    if(fFit){
      spcol.put();
      evt.put(std::move(assns));
      spcol_noreg.put();
    }
    return;
  }

  art::ServiceHandle<geo::Geometry> geom;

  bool is2view = false;
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
      // "X" view. For Argoneut and Lariat, collection=V is also the right
      // convention.
      if(hit->View() == geo::kZ){
        xhits.push_back(hit);
      }
      if(hit->View() == geo::kV){
        xhits.push_back(hit);
        is2view = true;
      }
      if(hit->View() == geo::kU || hit->View() == geo::kY){
        uhits.push_back(hit);
        is2view = true;
      }
    }
    else{
      if(hit->View() == geo::kU) uhits.push_back(hit);
      if(hit->View() == geo::kV) vhits.push_back(hit);
    }
  } // end for hit

  std::vector<raw::ChannelID_t> ubadchans, vbadchans;
  if(fAllowBadInductionHit){
    // NB - current implementation only allows intersection of vbadchan with
    // good hits from X and U. Would like to refactor before adding the
    // complexity of XV+U.
    for(raw::ChannelID_t cid: art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider().BadChannels()){
      if(geom->SignalType(cid) != geo::kCollection){
        if(geom->View(cid) == geo::kU) ubadchans.push_back(cid);
        if(geom->View(cid) == geo::kV) vbadchans.push_back(cid);
      }
    }
  }


  std::vector<CollectionWireHit*> cwires;
  // So we can find them all to free the memory
  std::vector<InductionWireHit*> iwires;

  HitMap_t hitmap;
  if(is2view)
    BuildSystemXU(xhits, uhits, cwires, iwires, fAlpha != 0, hitmap);
  else
    BuildSystemXUV(xhits, uhits, vhits, ubadchans, vbadchans,
                   cwires, iwires, fAlpha != 0, hitmap);

  FillSystemToSpacePoints(cwires, spcol_pre);
  spcol_pre.put();

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

    FillSystemToSpacePoints(cwires, spcol_noreg);
    spcol_noreg.put();

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

    FillSystemToSpacePointsAndAssns(cwires, hitmap, spcol, *assns);
    spcol.put();
    evt.put(std::move(assns));
  } // end if fFit

  for(InductionWireHit* i: iwires) delete i;
  for(CollectionWireHit* c: cwires) delete c;
}

} // end namespace reco3d

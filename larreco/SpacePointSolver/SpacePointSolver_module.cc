// Christopher Backhouse - bckhouse@fnal.gov

// Test file at Caltech: /nfs/raid11/dunesam/prodgenie_nu_dune10kt_1x2x6_mcc7.0/prodgenie_nu_dune10kt_1x2x6_63_20160811T171439_merged.root

// C/C++ standard libraries
#include <string>
#include <iostream>

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/ChargedSpacePointCreator.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "larreco/SpacePointSolver/HitReaders/IHitReader.h"

#include "Solver.h"
#include "TripletFinder.h"

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

private:
  void produce(art::Event& evt) override;
  void beginJob() override;

  void AddNeighbours(const std::vector<SpaceCharge*>& spaceCharges) const;

  typedef std::map<const WireHit*, const recob::Hit*> HitMap_t;

  void BuildSystem(const std::vector<HitTriplet>& triplets,
                   std::vector<CollectionWireHit*>& cwires,
                   std::vector<InductionWireHit*>& iwires,
                   std::vector<SpaceCharge*>& orphanSCs,
                   bool incNei,
                   HitMap_t& hitmap) const;

  void Minimize(const std::vector<CollectionWireHit*>& cwires,
                const std::vector<SpaceCharge*>& orphanSCs,
                double alpha,
                int maxiterations);

  /// return whether the point was inserted (only happens when it has charge)
  bool AddSpacePoint(const SpaceCharge& sc,
                     int id,
                     recob::ChargedSpacePointCollectionCreator& points) const;

  void FillSystemToSpacePoints(const std::vector<CollectionWireHit*>& cwires,
                               const std::vector<SpaceCharge*>& orphanSCs,
                               recob::ChargedSpacePointCollectionCreator& pts) const;

  void FillSystemToSpacePointsAndAssns(const std::vector<art::Ptr<recob::Hit>>& hitlist,
                                       const std::vector<CollectionWireHit*>& cwires,
                                       const std::vector<SpaceCharge*>& orphanSCs,
                                       const HitMap_t& hitmap,
                                       recob::ChargedSpacePointCollectionCreator& points,
                                       art::Assns<recob::SpacePoint, recob::Hit>& assn) const;

  std::string fHitLabel;

  bool fFit;
  bool fAllowBadInductionHit, fAllowBadCollectionHit;

  double fAlpha;

  double fDistThresh;
  double fDistThreshDrift;

  int fMaxIterationsNoReg;
  int fMaxIterationsReg;

  double fXHitOffset;

  const detinfo::DetectorProperties* detprop;
  const geo::GeometryCore* geom;
  std::unique_ptr<reco3d::IHitReader> fHitReader; ///<  Expt specific tool for reading hits
};

DEFINE_ART_MODULE(SpacePointSolver)

// ---------------------------------------------------------------------------
SpacePointSolver::SpacePointSolver(const fhicl::ParameterSet& pset) :
    EDProducer{pset},
    fHitLabel(pset.get<std::string>("HitLabel")),
    fFit(pset.get<bool>("Fit")),
    fAllowBadInductionHit(pset.get<bool>("AllowBadInductionHit")),
    fAllowBadCollectionHit(pset.get<bool>("AllowBadCollectionHit")),
    fAlpha(pset.get<double>("Alpha")),
    fDistThresh(pset.get<double>("WireIntersectThreshold")),
    fDistThreshDrift(pset.get<double>("WireIntersectThresholdDriftDir")),
    fMaxIterationsNoReg(pset.get<int>("MaxIterationsNoReg")),
    fMaxIterationsReg(pset.get<int>("MaxIterationsReg")),
    fXHitOffset(pset.get<double>("XHitOffset"))
{
  recob::ChargedSpacePointCollectionCreator::produces(producesCollector(), "pre");
  if(fFit){
    recob::ChargedSpacePointCollectionCreator::produces(producesCollector());
    produces<art::Assns<recob::SpacePoint, recob::Hit>>();
    recob::ChargedSpacePointCollectionCreator::produces(producesCollector(), "noreg");
  }

  fHitReader = art::make_tool<reco3d::IHitReader>(pset.get<fhicl::ParameterSet>("HitReaderTool"));
}

// ---------------------------------------------------------------------------
void SpacePointSolver::beginJob()
{
  detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();
  geom = art::ServiceHandle<geo::Geometry const>()->provider();
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

  std::cout << Ntests << " tests to find " << Nnei << " neighbours" << std::endl;
}

// ---------------------------------------------------------------------------
void SpacePointSolver::
BuildSystem(const std::vector<HitTriplet>& triplets,
            std::vector<CollectionWireHit*>& cwires,
            std::vector<InductionWireHit*>& iwires,
            std::vector<SpaceCharge*>& orphanSCs,
            bool incNei,
            HitMap_t& hitmap) const
{
  std::set<const recob::Hit*> ihits;
  std::set<const recob::Hit*> chits;
  for(const HitTriplet& trip: triplets){
    if(trip.x) chits.insert(trip.x);
    if(trip.u) ihits.insert(trip.u);
    if(trip.v) ihits.insert(trip.v);
  }

  std::map<const recob::Hit*, InductionWireHit*> inductionMap;
  for(const recob::Hit* hit: ihits){
    InductionWireHit* iwire = new InductionWireHit(hit->Channel(),
                                                   hit->Integral());
    inductionMap[hit] = iwire;
    iwires.emplace_back(iwire);
    hitmap[iwire] = hit;
  }

  std::map<const recob::Hit*, std::vector<SpaceCharge*>> collectionMap;
  std::map<const recob::Hit*, std::vector<SpaceCharge*>> collectionMapBad;

  std::set<InductionWireHit*> satisfiedInduction;

  for(const HitTriplet& trip: triplets){
    // Don't have a cwire object yet, set it later
    SpaceCharge* sc = new SpaceCharge(trip.pt.x,
                                      trip.pt.y,
                                      trip.pt.z,
                                      0,
                                      inductionMap[trip.u],
                                      inductionMap[trip.v]);

    if(trip.u && trip.v){
      collectionMap[trip.x].push_back(sc);
      if(trip.x){
        satisfiedInduction.insert(inductionMap[trip.u]);
        satisfiedInduction.insert(inductionMap[trip.v]);
      }
    }
    else{
      collectionMapBad[trip.x].push_back(sc);
    }
  }

  std::vector<SpaceCharge*> spaceCharges;

  for(const recob::Hit* hit: chits){
    // Find the space charges associated with this hit
    std::vector<SpaceCharge*>& scs = collectionMap[hit];
    if(scs.empty()){
      // If there are no full triplets try the triplets with one bad channel
      scs = collectionMapBad[hit];
    }
    else{
      // If there were good triplets, delete the bad hit ones
      for(SpaceCharge* sc: collectionMapBad[hit]) delete sc;
    }
    // Still no space points, don't bother making a wire
    if(scs.empty()) continue;

    CollectionWireHit* cwire = new CollectionWireHit(hit->Channel(),
                                                     hit->Integral(),
                                                     scs);
    hitmap[cwire] = hit;
    cwires.push_back(cwire);
    spaceCharges.insert(spaceCharges.end(), scs.begin(), scs.end());
    for(SpaceCharge* sc: scs) sc->fCWire = cwire;
  } // end for hit

  // Space charges whose collection wire is bad, which we have no other way of
  // addressing.
  for(SpaceCharge* sc: collectionMap[0]){
    // Only count orphans where an induction wire has no other explanation
    if(satisfiedInduction.count(sc->fWire1) == 0 ||
       satisfiedInduction.count(sc->fWire2) == 0){
      orphanSCs.push_back(sc);
    }
    else{
      delete sc;
    }
  }
  spaceCharges.insert(spaceCharges.end(), orphanSCs.begin(), orphanSCs.end());

  std::cout << cwires.size() << " collection wire objects" << std::endl;
  std::cout << spaceCharges.size() << " potential space points" << std::endl;

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
                        const std::vector<SpaceCharge*>& orphanSCs,
                        recob::ChargedSpacePointCollectionCreator& points) const
{
  int iPoint = 0;
  for(const CollectionWireHit* cwire: cwires){
    for(const SpaceCharge* sc: cwire->fCrossings){
      AddSpacePoint(*sc, iPoint++, points);
    } // for sc
  } // for cwire

  for(const SpaceCharge* sc: orphanSCs) AddSpacePoint(*sc, iPoint++, points);
}


// ---------------------------------------------------------------------------
void SpacePointSolver::
FillSystemToSpacePointsAndAssns(const std::vector<art::Ptr<recob::Hit>>& hitlist,
                                const std::vector<CollectionWireHit*>& cwires,
                                const std::vector<SpaceCharge*>& orphanSCs,
                                const HitMap_t& hitmap,
                                recob::ChargedSpacePointCollectionCreator& points,
                                art::Assns<recob::SpacePoint, recob::Hit>& assn) const
{
  std::map<const recob::Hit*, art::Ptr<recob::Hit>> ptrmap;
  for(art::Ptr<recob::Hit> hit: hitlist) ptrmap[hit.get()] = hit;

  std::vector<const SpaceCharge*> scs;
  for(const SpaceCharge* sc: orphanSCs) scs.push_back(sc);
  for(const CollectionWireHit* cwire: cwires)
    for(const SpaceCharge* sc: cwire->fCrossings)
      scs.push_back(sc);

  int iPoint = 0;

  for(const SpaceCharge* sc: scs){
    if(!AddSpacePoint(*sc, iPoint++, points)) continue;
    const auto& spsPtr = points.lastSpacePointPtr();

    if(sc->fCWire){
      assn.addSingle(spsPtr, ptrmap[hitmap.at(sc->fCWire)]);
    }
    if(sc->fWire1){
      assn.addSingle(spsPtr, ptrmap[hitmap.at(sc->fWire1)]);
    }
    if(sc->fWire2){
      assn.addSingle(spsPtr, ptrmap[hitmap.at(sc->fWire2)]);
    }
  }
}

// ---------------------------------------------------------------------------
void SpacePointSolver::Minimize(const std::vector<CollectionWireHit*>& cwires,
                                const std::vector<SpaceCharge*>& orphanSCs,
                                double alpha,
                                int maxiterations)
{
  double prevMetric = Metric(cwires, alpha);
  std::cout << "Begin: " << prevMetric << std::endl;
  for(int i = 0; i < maxiterations; ++i){
    Iterate(cwires, orphanSCs, alpha);
    const double metric = Metric(cwires, alpha);
    std::cout << i << " " << metric << std::endl;
    if(metric > prevMetric){
      std::cout << "Warning: metric increased" << std::endl;
      return;
    }
    if(fabs(metric-prevMetric) < 1e-3*fabs(prevMetric)) return;
    prevMetric = metric;
  }
}

// ---------------------------------------------------------------------------
void SpacePointSolver::produce(art::Event& evt)
{
  art::Handle<std::vector<recob::Hit>> hits;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if(evt.getByLabel(fHitLabel, hits))
    art::fill_ptr_vector(hitlist, hits);

  auto spcol_pre = recob::ChargedSpacePointCollectionCreator::forPtrs(evt, "pre");
  auto spcol_noreg = recob::ChargedSpacePointCollectionCreator::forPtrs(evt, "noreg");
  auto spcol = recob::ChargedSpacePointCollectionCreator::forPtrs(evt);
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

  art::ServiceHandle<geo::Geometry const> geom;

  std::vector<art::Ptr<recob::Hit>> xhits, uhits, vhits;
  bool is2view = fHitReader->readHits(hitlist, xhits, uhits, vhits);

  std::vector<raw::ChannelID_t> xbadchans, ubadchans, vbadchans;
  if(fAllowBadInductionHit || fAllowBadCollectionHit){
    for(raw::ChannelID_t cid: art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().BadChannels()){
      if(geom->SignalType(cid) == geo::kCollection){
        if(fAllowBadCollectionHit && geom->View(cid) == geo::kZ){
          xbadchans.push_back(cid);
        }
      }
      else{
        if(fAllowBadInductionHit){
          if(geom->View(cid) == geo::kU) ubadchans.push_back(cid);
          if(geom->View(cid) == geo::kV) vbadchans.push_back(cid);
        }
      }
    }
  }
  std::cout << xbadchans.size() << " X, "
            << ubadchans.size() << " U, "
            << vbadchans.size() << " V bad channels" << std::endl;


  std::vector<CollectionWireHit*> cwires;
  // So we can find them all to free the memory
  std::vector<InductionWireHit*> iwires;
  // Nodes with a bad collection wire that we otherwise can't address
  std::vector<SpaceCharge*> orphanSCs;

  HitMap_t hitmap;
  if(is2view){
    std::cout << "Finding 2-view coincidences..." << std::endl;
    TripletFinder tf(xhits, uhits, {},
                     xbadchans, ubadchans, {},
                     fDistThresh, fDistThreshDrift, fXHitOffset);
    BuildSystem(tf.TripletsTwoView(),
                cwires, iwires, orphanSCs,
                fAlpha != 0, hitmap);
  }
  else{
    std::cout << "Finding XUV coincidences..." << std::endl;
    TripletFinder tf(xhits, uhits, vhits,
                     xbadchans, ubadchans, vbadchans,
                     fDistThresh, fDistThreshDrift, fXHitOffset);
    BuildSystem(tf.Triplets(),
                cwires, iwires, orphanSCs,
                fAlpha != 0, hitmap);
  }

  FillSystemToSpacePoints(cwires, orphanSCs, spcol_pre);
  spcol_pre.put();

  if(fFit){
    std::cout << "Iterating with no regularization..." << std::endl;
    Minimize(cwires, orphanSCs, 0, fMaxIterationsNoReg);

    FillSystemToSpacePoints(cwires, orphanSCs, spcol_noreg);
    spcol_noreg.put();

    std::cout << "Now with regularization..." << std::endl;
    Minimize(cwires, orphanSCs, fAlpha, fMaxIterationsReg);

    FillSystemToSpacePointsAndAssns(hitlist, cwires, orphanSCs, hitmap, spcol, *assns);
    spcol.put();
    evt.put(std::move(assns));
  } // end if fFit

  for(InductionWireHit* i: iwires) delete i;
  for(CollectionWireHit* c: cwires) delete c;
  for(SpaceCharge* s: orphanSCs) delete s;
}

} // end namespace reco3d

//
// Name: FeatureTracker.h
//
// Purpose:  This module takes features found in 2D and uses them
//            to produce seeds for 3D tracking.
//
// Ben Jones, MIT
//

#include "TVector3.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larreco/RecoAlg/CornerFinderAlg.h"
#include "larreco/RecoAlg/SpacePointAlg.h"

namespace trkf {

  class FeatureTracker : public art::EDProducer {
  public:
    explicit FeatureTracker(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    void GetProjectedEnds(detinfo::DetectorPropertiesData const& detProp,
                          TVector3 xyz,
                          std::vector<double>& uvw,
                          std::vector<double>& t,
                          int tpc = 0,
                          int cryo = 0);

    std::map<int, std::vector<double>> ExtractEndPointTimes(
      std::vector<recob::EndPoint2D> EndPoints);

    std::vector<recob::SpacePoint> Get3DFeaturePoints(
      detinfo::DetectorClocksData const& clockData,
      detinfo::DetectorPropertiesData const& detProp,
      std::vector<recob::EndPoint2D> const& EndPoints,
      art::PtrVector<recob::Hit> const& Hits);

    std::vector<recob::Seed> GetValidLines(detinfo::DetectorPropertiesData const& detProp,
                                           std::vector<recob::SpacePoint>& sps,
                                           bool ApplyFilter = true);

    void RemoveDuplicatePaths(std::vector<recob::Seed>& Seeds,
                              std::vector<int>& ConnectionPoint1,
                              std::vector<int>& ConnectionPoint2);

    trkf::SpacePointAlg fSP;
    corner::CornerFinderAlg fCorner;

    std::string fHitModuleLabel;
    std::string fCalDataModuleLabel;

    double fLineIntFraction;
    double fLineIntThreshold;

    std::map<int, std::vector<double>> fEndPointTimes;
    art::ServiceHandle<geo::Geometry const> fGeometryHandle;
  };

  FeatureTracker::FeatureTracker(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    , fSP(pset.get<fhicl::ParameterSet>("SpacepointPset"))
    , fCorner(pset.get<fhicl::ParameterSet>("CornerPset"))
    , fHitModuleLabel{pset.get<std::string>("HitModuleLabel")}
    , fCalDataModuleLabel{pset.get<std::string>("CornerPset.CalDataModuleLabel")}
    , fLineIntFraction{pset.get<double>("LineIntFraction")}
    , fLineIntThreshold{pset.get<double>("LineIntThreshold")}
  {
    produces<std::vector<recob::Seed>>();
  }

  void
  FeatureTracker::produce(art::Event& evt)
  {

    // Extract hits PtrVector from event
    art::Handle<std::vector<recob::Hit>> hith;
    evt.getByLabel(fHitModuleLabel, hith);

    art::PtrVector<recob::Hit> hitvec;
    for (unsigned int i = 0; i < hith->size(); ++i) {
      art::Ptr<recob::Hit> prod(hith, i);
      hitvec.push_back(prod);
    }

    //We need to grab out the wires.
    art::Handle<std::vector<recob::Wire>> wireHandle;
    evt.getByLabel(fCalDataModuleLabel, wireHandle);
    std::vector<recob::Wire> const& wireVec(*wireHandle);

    //First, have it process the wires.
    fCorner.GrabWires(wireVec, *fGeometryHandle);

    std::vector<recob::EndPoint2D> EndPoints;
    fCorner.get_feature_points(EndPoints, *fGeometryHandle);

    fEndPointTimes = ExtractEndPointTimes(EndPoints);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    std::vector<recob::SpacePoint> sps = Get3DFeaturePoints(clockData, detProp, EndPoints, hitvec);
    std::vector<recob::Seed> SeedsToStore = GetValidLines(detProp, sps, true);

    std::unique_ptr<std::vector<recob::Seed>> seeds(new std::vector<recob::Seed>);

    for (size_t i = 0; i != SeedsToStore.size(); ++i)
      seeds->push_back(SeedsToStore.at(i));

    evt.put(std::move(seeds));
  }

  //---------------------------------------------------------------------

  std::map<int, std::vector<double>>
  FeatureTracker::ExtractEndPointTimes(std::vector<recob::EndPoint2D> EndPoints)
  {
    std::map<int, std::vector<double>> EndPointTimesInPlane;
    for (size_t i = 0; i != EndPoints.size(); ++i) {
      EndPointTimesInPlane[EndPoints.at(i).View()].push_back(EndPoints.at(i).DriftTime());
    }

    for (std::map<int, std::vector<double>>::iterator itEpTime = EndPointTimesInPlane.begin();
         itEpTime != EndPointTimesInPlane.end();
         ++itEpTime) {
      std::sort(itEpTime->second.begin(), itEpTime->second.end());
    }
    return EndPointTimesInPlane;
  }

  //---------------------------------------------------------------------

  std::vector<recob::Seed>
  FeatureTracker::GetValidLines(detinfo::DetectorPropertiesData const& detProp,
                                std::vector<recob::SpacePoint>& spts,
                                bool ApplyFilter)
  {
    std::vector<recob::Seed> SeedsToReturn;

    std::vector<int> ConnectionPoint1;
    std::vector<int> ConnectionPoint2;
    std::map<int, std::vector<int>> SeedConnections;

    for (size_t i = 0; i != spts.size(); ++i) {
      for (size_t j = 0; j != i; ++j) {

        TVector3 xyz_i;
        TVector3 xyz_j;

        std::vector<double> t_i, t_j;

        std::vector<double> uvw_i;
        std::vector<double> uvw_j;

        for (size_t d = 0; d != 3; ++d) {
          xyz_i[d] = spts.at(i).XYZ()[d];
          xyz_j[d] = spts.at(j).XYZ()[d];
        }

        GetProjectedEnds(detProp, xyz_i, uvw_i, t_i, 0, 0);
        GetProjectedEnds(detProp, xyz_j, uvw_j, t_j, 0, 0);

        bool ThisLineGood = true;

        for (size_t p = 0; p != uvw_i.size(); ++p) {
          TH2F const& RawHist = fCorner.GetWireDataHist(p);

          double lineint = fCorner.line_integral(
            RawHist, uvw_i.at(p), t_i.at(p), uvw_j.at(p), t_j.at(p), fLineIntThreshold);

          if (lineint < fLineIntFraction) { ThisLineGood = false; }
        }
        if (ThisLineGood) {
          double Err[3];
          double Pos[3];
          double Dir[3];

          for (size_t d = 0; d != 3; ++d) {
            Pos[d] = 0.5 * (xyz_i[d] + xyz_j[d]);
            Dir[d] = 0.5 * (xyz_i[d] - xyz_j[d]);
            Err[d] = 0;
          }

          ConnectionPoint1.push_back(i);
          ConnectionPoint2.push_back(j);

          SeedsToReturn.push_back(recob::Seed(Pos, Dir, Err, Err));
        }
      }
    }

    if (ApplyFilter) {
      RemoveDuplicatePaths(SeedsToReturn, ConnectionPoint1, ConnectionPoint2);
      mf::LogInfo("FeatureTracker")
        << "Seeds after filter " << SeedsToReturn.size() << " seeds" << std::endl;
    }

    return SeedsToReturn;
  }

  //--------------------------------------------------

  void
  FeatureTracker::RemoveDuplicatePaths(std::vector<recob::Seed>& Seeds,
                                       std::vector<int>& ConnectionPoint1,
                                       std::vector<int>& ConnectionPoint2)
  {

    std::map<int, bool> MarkedForRemoval;

    std::map<int, std::vector<int>> SeedsSharingPoint;
    for (size_t i = 0; i != Seeds.size(); ++i) {
      SeedsSharingPoint[ConnectionPoint1[i]].push_back(i);
      SeedsSharingPoint[ConnectionPoint2[i]].push_back(i);
    }

    for (size_t s = 0; s != Seeds.size(); ++s) {

      int StartPt = ConnectionPoint1.at(s);
      int EndPt = ConnectionPoint2.at(s);
      int MidPt = -1;

      for (size_t SeedsWithThisStart = 0; SeedsWithThisStart != SeedsSharingPoint[StartPt].size();
           SeedsWithThisStart++) {
        int i = SeedsSharingPoint[StartPt].at(SeedsWithThisStart);
        if (ConnectionPoint1.at(i) == StartPt)
          MidPt = ConnectionPoint2.at(i);
        else if (ConnectionPoint2.at(i) == StartPt)
          MidPt = ConnectionPoint1.at(i);

        for (size_t SeedsWithThisMid = 0; SeedsWithThisMid != SeedsSharingPoint[MidPt].size();
             SeedsWithThisMid++) {
          int j = SeedsSharingPoint[MidPt].at(SeedsWithThisMid);
          if ((ConnectionPoint1.at(j) == EndPt) || (ConnectionPoint2.at(j) == EndPt)) {

            double Lengthi = Seeds.at(i).GetLength();
            double Lengthj = Seeds.at(j).GetLength();
            double Lengths = Seeds.at(s).GetLength();

            if ((Lengths > Lengthi) && (Lengths > Lengthj)) {
              MarkedForRemoval[i] = true;
              MarkedForRemoval[j] = true;
            }

            if ((Lengthi > Lengths) && (Lengthi > Lengthj)) {
              MarkedForRemoval[s] = true;
              MarkedForRemoval[j] = true;
            }
            if ((Lengthj > Lengthi) && (Lengthj > Lengths)) {
              MarkedForRemoval[s] = true;
              MarkedForRemoval[i] = true;
            }
          }
        }
      }
    }
    for (std::map<int, bool>::reverse_iterator itrem = MarkedForRemoval.rbegin();
         itrem != MarkedForRemoval.rend();
         ++itrem) {
      if (itrem->second == true) {
        Seeds.erase(Seeds.begin() + itrem->first);
        ConnectionPoint1.erase(ConnectionPoint1.begin() + itrem->first);
        ConnectionPoint2.erase(ConnectionPoint2.begin() + itrem->first);
      }
    }
  }

  //---------------------------------------------------------------------

  void
  FeatureTracker::GetProjectedEnds(detinfo::DetectorPropertiesData const& detProp,
                                   TVector3 xyz,
                                   std::vector<double>& uvw,
                                   std::vector<double>& t,
                                   int tpc,
                                   int cryo)
  {
    art::ServiceHandle<geo::Geometry const> geo;

    int NPlanes = geo->Cryostat(cryo).TPC(tpc).Nplanes();

    uvw.resize(NPlanes);
    t.resize(NPlanes);

    for (int plane = 0; plane != NPlanes; plane++) {
      uvw[plane] = geo->NearestWire(xyz, plane, tpc, cryo);
      t[plane] = detProp.ConvertXToTicks(xyz[0], plane, tpc, cryo);
    }
  }

  //----------------------------------------------------------------------

  std::vector<recob::SpacePoint>
  FeatureTracker::Get3DFeaturePoints(detinfo::DetectorClocksData const& clockData,
                                     detinfo::DetectorPropertiesData const& detProp,
                                     std::vector<recob::EndPoint2D> const& EndPoints,
                                     art::PtrVector<recob::Hit> const& Hits)
  {

    art::PtrVector<recob::Hit> HitsForEndPoints;

    // Loop through the hits looking for the ones which match corners
    for (std::vector<recob::EndPoint2D>::const_iterator itEP = EndPoints.begin();
         itEP != EndPoints.end();
         ++itEP) {
      int EPMatchCount = 0;

      for (art::PtrVector<recob::Hit>::const_iterator itHit = Hits.begin(); itHit != Hits.end();
           ++itHit) {

        if ((itEP->View() == (*itHit)->View()) &&
            (itEP->WireID().Wire == (*itHit)->WireID().Wire)) {
          HitsForEndPoints.push_back(*itHit);
          EPMatchCount++;
        }
      }
    }
    std::vector<recob::SpacePoint> spts;
    fSP.makeSpacePoints(clockData, detProp, HitsForEndPoints, spts);

    for (size_t i = 0; i != spts.size(); ++i) {
      for (size_t p = 0; p != 3; ++p) {
        double Closest = 100000;
        double spt_t = detProp.ConvertXToTicks(spts.at(i).XYZ()[0], p, 0, 0);
        for (size_t epTime = 0; epTime != fEndPointTimes[p].size(); ++epTime) {
          if (fabs(fEndPointTimes[p].at(epTime) - spt_t) < Closest) {
            Closest = fabs(fEndPointTimes[p].at(epTime) - spt_t);
          }
        }
      }
    }
    return spts;
  }

  DEFINE_ART_MODULE(FeatureTracker)
}

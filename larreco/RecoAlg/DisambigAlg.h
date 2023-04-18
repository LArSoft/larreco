/////////////////////////////////////////////////////////////////
//  \fileDisambigAlg.h
//  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////
#ifndef DisambigAlg_H
#define DisambigAlg_H

#include <map>
#include <utility> // std::pair<>
#include <vector>

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace fhicl {
  class ParamterSet;
}

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/APAGeometryAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace apa {

  //---------------------------------------------------------------
  class DisambigAlg {
  public:
    explicit DisambigAlg(fhicl::ParameterSet const& pset);

    void RunDisambig(detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     art::Handle<std::vector<recob::Hit>> GausHits);

    void TrivialDisambig(detinfo::DetectorClocksData const& clockData,
                         detinfo::DetectorPropertiesData const& detProp,
                         unsigned int apa); ///< Make the easiest and safest disambiguations in apa
    void Crawl(unsigned int apa);           ///< Extend what we disambiguation we do have in apa
    unsigned int FindChanTimeEndPts(detinfo::DetectorPropertiesData const& detProp,
                                    unsigned int apa); ///< Basic endpoint-hit finder per apa
    void UseEndPts(detinfo::DetectorPropertiesData const& detProp,
                   unsigned int apa); ///< Try to associate endpoint hits and
                                      ///< crawl from there
    unsigned int CompareViews(
      detinfo::DetectorPropertiesData const& detProp,
      unsigned int apa); ///< Compare U and V to see if one says something about the other
    void AssessDisambigSoFar(
      unsigned int apa); ///< See how much disambiguation has been done in this apa so far

    std::map<unsigned int, double> fUeffSoFar;
    std::map<unsigned int, double> fVeffSoFar;
    std::map<unsigned int, unsigned int> fnUSoFar;
    std::map<unsigned int, unsigned int> fnVSoFar;
    std::map<unsigned int, unsigned int> fnDUSoFar;
    std::map<unsigned int, unsigned int> fnDVSoFar;

    std::vector<std::pair<art::Ptr<recob::Hit>, geo::WireID>> fDisambigHits;
    ///< The final list of hits to pass back to be made

  private:
    // other classes we will use
    apa::APAGeometryAlg fAPAGeo;
    art::ServiceHandle<geo::Geometry const> geom;
    // **temporarily** here to look at performance without noise hits
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

    // Hits organization
    std::map<raw::ChannelID_t, std::vector<art::Ptr<recob::Hit>>> fChannelToHits;
    std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>> fAPAToUVHits, fAPAToZHits;
    std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>> fAPAToHits;
    ///\ todo: Channel/APA to hits can be done in a unified way
    std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>> fAPAToEndPHits;
    std::map<unsigned int, std::vector<std::pair<art::Ptr<recob::Hit>, geo::WireID>>> fAPAToDHits;
    ///< Hold the disambiguations per APA

    // data/function to keep track of disambiguation along the way
    std::map<std::pair<double, double>, geo::WireID> fChanTimeToWid;
    ///< If a hit is disambiguated, map its chan and peak time to the chosen wireID
    std::map<unsigned int, std::map<std::pair<double, double>, bool>> fHasBeenDisambiged;
    ///< Convenient way to keep track of disambiguation so far
    void MakeDisambigHit(art::Ptr<recob::Hit> const& hit, geo::WireID, unsigned int apa);
    ///< Makes a disambiguated hit while keeping track of what has already been disambiguated

    // Functions that support disambiguation methods
    unsigned int MakeCloseHits(int ext, geo::WireID wid, double Dmin, double Dmax);
    ///< Having disambiguated a time range on a wireID, extend to neighboring channels
    bool HitsOverlapInTime(detinfo::DetectorPropertiesData const& detProp,
                           recob::Hit const& hitA,
                           recob::Hit const& hitB);
    bool HitsReasonablyMatch(art::Ptr<recob::Hit> hitA, art::Ptr<recob::Hit> hitB);
    ///\ todo: Write function that compares hits more detailedly

    // Configure the disambiguation
    bool fCrawl;
    bool fUseEndP;
    bool fCompareViews;
    unsigned int fNChanJumps; ///< Number of channels the crawl can jump over
    double fCloseHitsRadius;  ///< Distance (cm) away from a hit to look when
                              ///< checking if it's an endpoint
    double fMaxEndPDegRange;  ///< Within the close hits radius, how spread can
                              ///< the majority of the activity be around a
                              ///< possible endpoint

  }; // class DisambigAlg

} // namespace apa

#endif // ifndef DisambigAlg_H

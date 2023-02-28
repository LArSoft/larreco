////////////////////////////////////////////////////////////////////////
//
// \file TrackStitcher
//
// \author echurch@fnal.gov
//
//  This algorithm is designed to join tracks that point in roughly same direction
//  and whose endpoints are suitably close.
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <ostream>
#include <utility>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/RecoAlg/StitchAlg.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <string>
#include <vector>

namespace trkf {

  class TrackStitcher : public art::EDProducer {
  public:
    explicit TrackStitcher(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    art::PtrVector<recob::Hit> GetHitsFromComponentTracks(const art::PtrVector<recob::Track>&,
                                                          const art::Event& evt);
    art::PtrVector<recob::SpacePoint> GetSpacePointsFromComponentTracks(
      const art::PtrVector<recob::Track>&,
      const art::Event& evt);
    std::vector<art::Ptr<recob::Hit>> GetHitsFromAssdSpacePoints(
      const art::PtrVector<recob::SpacePoint>&,
      const art::Event& evt,
      std::vector<std::pair<std::vector<art::Ptr<recob::Hit>>::const_iterator,
                            std::vector<art::Ptr<recob::Hit>>::const_iterator>>& vpi);
    std::string fTrackModuleLabel; // label for input collection
    std::string fSpptModuleLabel;  // label for input collection
    bool fStizatch;                // CommonComponentStitch
    StitchAlg fStitchAlg;

  }; // class TrackStitcher

} // end namespace for declarations

namespace trkf {

  //-------------------------------------------------
  TrackStitcher::TrackStitcher(fhicl::ParameterSet const& pset)
    : EDProducer{pset}, fStitchAlg(pset.get<fhicl::ParameterSet>("StitchAlg"))
  {
    fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");
    fSpptModuleLabel = pset.get<std::string>("SpptModuleLabel");
    fStizatch = pset.get<bool>("CommonComponentStitch", true);

    produces<std::vector<recob::Track>>();
    produces<std::vector<art::PtrVector<recob::Track>>>();
    produces<art::Assns<recob::Track, recob::Hit>>();
    produces<art::Assns<recob::Track, recob::SpacePoint>>();
    produces<art::Assns<recob::SpacePoint, recob::Hit>>();
  }

  //------------------------------------------------------------------------------------//
  void TrackStitcher::produce(art::Event& evt)
  {

    // get services
    art::ServiceHandle<geo::Geometry const> geom;

    //////////////////////////////////////////////////////
    // Make a std::unique_ptr<> for the thing you want to put into the event
    //////////////////////////////////////////////////////
    // tcol is the collection of new tracks
    std::unique_ptr<std::vector<recob::Track>> tcol(new std::vector<recob::Track>);
    std::unique_ptr<art::PtrVector<recob::SpacePoint>> scol(new art::PtrVector<recob::SpacePoint>);
    // tvcol is the collection of vectors that comprise each tcol
    std::unique_ptr<std::vector<art::PtrVector<recob::Track>>> tvcol(
      new std::vector<art::PtrVector<recob::Track>>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>> thassn(
      new art::Assns<recob::Track, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint>> tsptassn(
      new art::Assns<recob::Track, recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit>> spthassn(
      new art::Assns<recob::SpacePoint, recob::Hit>);

    // Get the original Spacepoints. Trackers other than CosmicTracker wrote the
    // SpacePoints as a PtrVec of vecs. If they look like that, flatten into one vec.

    art::Handle<std::vector<recob::SpacePoint>> sppth;
    try {
      mf::LogWarning("TrackStitcher")
        << "Trying to read Track3DKalmanXYZ-style PtrVector of std::vector of SpacePoints"
        << std::endl;
      art::Handle<std::vector<art::PtrVector<recob::SpacePoint>>> sppth;
      evt.getByLabel(fSpptModuleLabel, sppth);
      for (size_t ii = 0; ii < sppth->size(); ii++)
        for (size_t jj = 0; jj < sppth->at(ii).size(); ii++) {
          art::Ptr<recob::SpacePoint> sptmp(sppth->at(ii).at(jj));
          scol->push_back(sptmp);
        }
    }
    catch (...) {
      mf::LogWarning("TrackStitcher")
        << "Trying instead to read CosmicTracker-style already-flattened vector of SpacePoints"
        << std::endl;
      art::Handle<std::vector<recob::SpacePoint>> sppthf;
      evt.getByLabel(fSpptModuleLabel, sppthf);
      for (size_t ii = 0; ii < sppthf->size(); ii++) {
        art::Ptr<recob::SpacePoint> sptmpf(sppthf, ii);
        scol->push_back(sptmpf);
      }
    }

    // Find the best match for each track's Head and Tail.
    fStitchAlg.FindHeadsAndTails(evt, fTrackModuleLabel);
    // walk through each vertex of one track to its match on another, and so on and stitch 'em.
    fStitchAlg.WalkStitch();
    // search composite tracks and stitch further if there are components in common. Do it until all are stitched.
    bool stizatch(fStizatch);
    while (stizatch) {
      stizatch = fStitchAlg.CommonComponentStitch();
    }
    mf::LogVerbatim("TrackStitcher.beginning") << "There are " << fStitchAlg.ftListHandle->size()
                                               << " Tracks in this event before stitching.";

    fStitchAlg.GetTracks(*tcol);
    fStitchAlg.GetTrackComposites(*tvcol);

    if (tcol->size() != tvcol->size())
      throw cet::exception("TrackStitcher")
        << "Tracks and TrackComposites do not match: " << tcol->size() << " vs " << tvcol->size()
        << "\n";

    std::vector<size_t> spIndices(scol->size());
    // create spIndices, index array for searching into original scol SpacePoints.
    for (size_t ii = 0; ii < scol->size(); ii++) {
      spIndices[ii] = ii;
    }

    for (size_t ii = 0; ii < tvcol->size(); ii++) {
      const art::PtrVector<recob::Hit>& hits(GetHitsFromComponentTracks(tvcol->at(ii), evt));
      // Now make the Assns of relevant Hits to stitched Track
      util::CreateAssn(*this, evt, *tcol, hits, *thassn, ii);
      const art::PtrVector<recob::SpacePoint>& sppts(
        GetSpacePointsFromComponentTracks(tvcol->at(ii), evt));
      // Now make the Assns of relevant Sppts to stitched Track
      util::CreateAssn(*this, evt, *tcol, sppts, *tsptassn, ii);

      // Now Assns of sppts to hits. For this Sppt
      // I call the function to bring back the vec of associated Hits and the vector of
      // pairs of iterators that allow to pull those Hits needed from each Sppt.
      std::vector<std::pair<std::vector<art::Ptr<recob::Hit>>::const_iterator,
                            std::vector<art::Ptr<recob::Hit>>::const_iterator>>
        pits;

      std::vector<art::Ptr<recob::Hit>> hitsFromSppts;
      art::FindManyP<recob::Hit> hitAssns(sppts, evt, fSpptModuleLabel);

      size_t start(0), finish(0);
      for (unsigned int ii = 0; ii < sppts.size(); ++ii) {
        hitsFromSppts.insert(hitsFromSppts.end(), hitAssns.at(ii).begin(), hitAssns.at(ii).end());
        finish = start + (size_t)(hitAssns.at(ii).end() - hitAssns.at(ii).begin());
        std::pair<std::vector<art::Ptr<recob::Hit>>::const_iterator,
                  std::vector<art::Ptr<recob::Hit>>::const_iterator>
          pithittmp(hitAssns.at(ii).begin(), hitAssns.at(ii).end());
        pits.push_back(pithittmp);
        start += (finish + 1);
      }
      //	std::cout << "TrackStitcher_module: scol->size() is " << scol->size() << std::endl;
      //	std::cout << "TrackStitcher_module: sppts.size() is " << sppts.size() << std::endl;
      for (size_t jj = 0; jj < sppts.size(); jj++) {
        // find jjth sppt in the list of scol. Meaning, find kkth element of sppth.
        size_t ll(scol->size());
        // this gives indices into the vector of original spacepoints in which
        // to look for our sppts.
        size_t off(0);
        for (auto& kk : spIndices) {
          const art::Ptr<recob::SpacePoint> spptnc(scol->at(kk));
          if (spptnc != sppts.at(jj)) {
            off++;
            continue;
          }
          ll = kk;
          //		std::cout << "TrackStitcher_module: index into spacepoints for which to write out sppt-hit Assns is " << ll << std::endl;
          // drop this one for future searches, since we've used it.

          break;
        }
        if (ll < scol->size()) {
          std::vector<art::Ptr<recob::Hit>> hitsThisSppt;
          hitsThisSppt.insert(hitsThisSppt.begin(), pits.at(jj).first, pits.at(jj).second);
          util::CreateAssn(*this, evt, scol->at(ll), hitsThisSppt, *spthassn);
        }
      }
    }

    mf::LogVerbatim("TrackStitcher.end")
      << "There are " << tvcol->size() << " Tracks in this event after stitching.";
    evt.put(std::move(tcol));
    evt.put(std::move(tvcol));
    // Add Hit-to-Track and Sppt-to-Track Assns.
    evt.put(std::move(thassn));
    evt.put(std::move(tsptassn));
    evt.put(std::move(spthassn));
  }

  art::PtrVector<recob::Hit> TrackStitcher::GetHitsFromComponentTracks(
    const art::PtrVector<recob::Track>& tcomp,
    const art::Event& evtGHFCT)
  {

    art::PtrVector<recob::Hit> hits;
    art::FindManyP<recob::Hit> hitAssns(tcomp, evtGHFCT, fTrackModuleLabel);

    for (unsigned int ii = 0; ii < tcomp.size(); ++ii) {
      hits.insert(hits.end(), hitAssns.at(ii).begin(), hitAssns.at(ii).end());
    }

    //    const art::PtrVector<recob::Hit> chits(hits);
    return hits;
  }

  art::PtrVector<recob::SpacePoint> TrackStitcher::GetSpacePointsFromComponentTracks(
    const art::PtrVector<recob::Track>& tcomp,
    const art::Event& evtGHFCT)
  {

    art::PtrVector<recob::SpacePoint> sppts;
    art::FindManyP<recob::SpacePoint> spptAssns(tcomp, evtGHFCT, fTrackModuleLabel);
    for (unsigned int ii = 0; ii < tcomp.size(); ++ii) {
      sppts.insert(sppts.end(), spptAssns.at(ii).begin(), spptAssns.at(ii).end());
    }

    //    const art::PtrVector<recob::Hit> chits(hits);
    return sppts;
  }

  std::vector<art::Ptr<recob::Hit>> TrackStitcher::GetHitsFromAssdSpacePoints(
    const art::PtrVector<recob::SpacePoint>& sppts,
    const art::Event& evtGHFCT,
    std::vector<std::pair<std::vector<art::Ptr<recob::Hit>>::const_iterator,
                          std::vector<art::Ptr<recob::Hit>>::const_iterator>>& pithit)
  {

    std::vector<art::Ptr<recob::Hit>> hits;
    art::FindManyP<recob::Hit> hitAssns(sppts, evtGHFCT, fSpptModuleLabel);

    size_t start(0), finish(0);
    for (unsigned int ii = 0; ii < sppts.size(); ++ii) {
      hits.insert(hits.end(), hitAssns.at(ii).begin(), hitAssns.at(ii).end());
      finish = start + (size_t)(hitAssns.at(ii).end() - hitAssns.at(ii).begin());
      std::pair<std::vector<art::Ptr<recob::Hit>>::const_iterator,
                std::vector<art::Ptr<recob::Hit>>::const_iterator>
        pithittmp(hitAssns.at(ii).begin(), hitAssns.at(ii).end());
      pithit.push_back(pithittmp);
      start += (finish + 1);
    }

    return hits;
  }

  DEFINE_ART_MODULE(TrackStitcher)

} // end namespace

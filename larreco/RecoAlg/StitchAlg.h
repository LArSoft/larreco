/////////////////////////////////////////////////////////////////
//  \StitchAlg.h
//  echurch@fnal.gov
////////////////////////////////////////////////////////////////////

#ifndef STITCHALG_H
#define STITCHALG_H

// art libraries
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/PtrVector.h"

// LArSoft libraries
#include "lardataobj/RecoBase/Track.h"

namespace fhicl {
  class ParameterSet;
}
namespace art {
  class Event;
}

// C/C++ standard libraries
#include <string>
#include <tuple>
#include <vector>

namespace trkf {

  class StitchAlg {

  public:
    StitchAlg(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& pset);

    void FindHeadsAndTails(const art::Event& e, const std::string& t);
    void FirstStitch(const std::vector<art::PtrVector<recob::Track>>::iterator itvvArg,
                     const std::vector<recob::Track>::iterator itvArg);
    void WalkStitch();
    bool CommonComponentStitch();

    void GetTrackComposites(std::vector<art::PtrVector<recob::Track>>& c) const
    {
      c = fTrackComposite;
    }
    void GetTracks(std::vector<recob::Track>& t) const { t = fTrackVec; }

    art::Handle<std::vector<recob::Track>> ftListHandle;

  private:
    std::vector<std::tuple<std::string, int, int, double, double>> fh;
    std::vector<std::tuple<std::string, int, int, double, double>> ft;
    int ftNo;
    double fCosAngTol;
    double fSepTol;

    std::vector<art::PtrVector<recob::Track>> fTrackComposite;
    std::vector<recob::Track> fTrackVec;
    std::vector<std::vector<std::string>> fHT;
  };

} // namespace

#endif // ifndef STITCHALG_H

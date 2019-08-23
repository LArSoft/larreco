////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgStitching
// Author:      L.Whitehead (leigh.howard.whitehead@cern.ch) January 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgStitching_h
#define PMAlgStitching_h

#include <map>

#include "fhiclcpp/types/Atom.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

class TVector3;

namespace pma{
  class PMAlgStitching;
  class TrkCandidateColl;
}

class pma::PMAlgStitching{

public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int> StitchingThreshold {
      Name("StitchingThreshold"),
      Comment("The maximum value allowed for the stitching score. Has dimensions of length, and 10.0(cm) is a reasonable value.")
    };

    fhicl::Atom<unsigned int> NodesFromEnd{
      Name("NodesFromEnd"),
      Comment("Number of nodes we step back from the ends of the tracks to perform the stitching extrapolation.")
    };

  };

  // Constructor
  PMAlgStitching(const pma::PMAlgStitching::Config &config);

  // CPA and APA stitching wrappers
  void StitchTracksCPA(pma::TrkCandidateColl &tracks);
  void StitchTracksAPA(pma::TrkCandidateColl &tracks);

private:
  // Main function of the algorithm
  void StitchTracks(pma::TrkCandidateColl &tracks, bool isCPA);

  double GetOptimalStitchShift(TVector3 &pos1, TVector3 &pos2, TVector3 &dir1, TVector3 &dir2, double &shift) const;
  double GetTrackPairDelta(TVector3 &pos1, TVector3 &pos2, TVector3 &dir1, TVector3 &dir2) const;

  void GetTPCXOffsets();
  double GetTPCOffset(unsigned int tpc, unsigned int cryo, bool isCPA);

  std::map<geo::TPCID,double> fTPCXOffsetsAPA;
  std::map<geo::TPCID,double> fTPCXOffsetsCPA;

  // Tuneable parameters
  double fStitchingThreshold; // The maximum stitching score allowed for a successful stitch.
  unsigned int fNodesFromEnd; // Number of nodes we step back to make the stitch extrapolation. Require to mitigate end effects on tracks.

};

#endif

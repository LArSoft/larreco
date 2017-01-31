////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgStitching
// Author:      L.Whitehead (leigh.howard.whitehead@cern.ch) January 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgStitching_h
#define PMAlgStitching_h

#include <map>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "TVector3.h"

namespace pma{
  class PMAlgStitching;
}

class pma::PMAlgStitching{

public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Sequence<int> TrackingOnlyPdg {
      Name("ExampleThing"),
      Comment("DescribeExampleThing")
    };
  };

  // Constructor
//  PMAlgStitching(pma::TrkCandidateColl &inputTracks, const pma::PMAlgStitching::Config &config);
  PMAlgStitching(pma::TrkCandidateColl &inputTracks);

  // Destructor
  ~PMAlgStitching();

  // CPA and APA stitching wrappers
  void StitchTracksCPA();
  void StitchTracksAPA();

private: 
  // Main function of the algorithm
  void StitchTracks(bool isCPA);
  
  double GetTrackPairDelta(TVector3 &pos1, TVector3 &pos2, TVector3 &dir1, TVector3 &dir2);

  void GetTPCXOffsets();
  double GetTPCOffset(unsigned int tpc, unsigned int cryo, bool isCPA);

  std::map<geo::TPCID,double> fTPCXOffsetsAPA;
  std::map<geo::TPCID,double> fTPCXOffsetsCPA;

  pma::TrkCandidateColl& fInputTracks;
};

#endif


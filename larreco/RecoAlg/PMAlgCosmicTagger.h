////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgCosmicTagger
// Author:      L. Whitehead (leigh.howard.whitehead@cern.ch),
//              R. Sulej (robert.sulej@cern.ch) March 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgCosmicTagger_h
#define PMAlgCosmicTagger_h

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

class TVector3;

namespace detinfo {
  class DetectorClocksData;
}

namespace pma {
  class PMAlgCosmicTagger;
  class TrkCandidateColl;
}

#include <vector>

class pma::PMAlgCosmicTagger {

public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<bool> TagOutOfDriftTracks{Name("TagOutOfDriftTracks"),
                                          Comment("Tag tracks sticking out of 1 drift window.")};
    fhicl::Atom<double> OutOfDriftMargin{
      Name("OutOfDriftMargin"),
      Comment("The minimum distance beyond 1 drift window required for tagging track as a cosmic "
              "background.")};

    fhicl::Atom<bool> TagFullHeightTracks{Name("TagFullHeightTracks"),
                                          Comment("Tag tracks crossing full detector height")};
    fhicl::Atom<bool> TagFullWidthTracks{Name("TagFullWidthTracks"),
                                         Comment("Tag tracks crossing full detector width")};
    fhicl::Atom<bool> TagFullLengthTracks{Name("TagFullLengthTracks"),
                                          Comment("Tag tracks crossing full detector length")};
    fhicl::Atom<double> FullCrossingMargin{
      Name("FullCrossingMargin"),
      Comment("The maximum distance between the track length and detector length for full detector "
              "crossing tracks.")};
    fhicl::Atom<bool> TagNonBeamT0Tracks{
      Name("TagNonBeamT0Tracks"),
      Comment("Tag particles with reconstructed T0 not consistent with the beam")};
    fhicl::Atom<double> NonBeamT0Margin{
      Name("NonBeamT0Margin"),
      Comment("Tag only those events at least <margin> from the beam time")};
    fhicl::Atom<bool> TagTopFrontBack{Name("TagTopFrontBack"),
                                      Comment("Tag tracks that enter through the top of the "
                                              "detector and exit through the front or back")};
    fhicl::Atom<double> TopFrontBackMargin{
      Name("TopFrontBackMargin"),
      Comment("Distance tolerence from the top, front and back of the detector")};
    fhicl::Atom<bool> TagApparentStopper{
      Name("TagApparentStopper"),
      Comment("Tag tracks that enter through the top of the detector appear to stop (without "
              "seeing evidence of stopping)")};
    fhicl::Atom<double> ApparentStopperMargin{
      Name("ApparentStopperMargin"),
      Comment(
        "Distance tolerence from the top of the detector to be considered coming in from the top")};
    fhicl::Atom<bool> VetoActualStopper{
      Name("VetoActualStopper"),
      Comment("If true: use de/dx information to identify stopping muons but do not tag them")};
    fhicl::Atom<double> StopperBuffer{
      Name("StopperBuffer"),
      Comment(
        "Should find no tracks starting within this distance from the end point of the track")};
  };

  PMAlgCosmicTagger(const pma::PMAlgCosmicTagger::Config& config)
    : fTagOutOfDriftTracks(config.TagOutOfDriftTracks())
    , fOutOfDriftMargin(config.OutOfDriftMargin())
    ,

    fTagFullHeightTracks(config.TagFullHeightTracks())
    , fTagFullWidthTracks(config.TagFullWidthTracks())
    , fTagFullLengthTracks(config.TagFullLengthTracks())
    , fFullCrossingMargin(config.FullCrossingMargin())
    ,

    fTagNonBeamT0Tracks(config.TagNonBeamT0Tracks())
    , fNonBeamT0Margin(config.NonBeamT0Margin())
    , fTagTopFrontBack(config.TagTopFrontBack())
    , fTopFrontBackMargin(config.TopFrontBackMargin())
    , fTagApparentStopper(config.TagApparentStopper())
    , fApparentStopperMargin(config.ApparentStopperMargin())
    , fVetoActualStopper(config.VetoActualStopper())
    , fStopperBuffer(config.StopperBuffer())
  {}

  bool
  tagAny() const
  {
    return (fTagOutOfDriftTracks || fTagFullHeightTracks || fTagFullWidthTracks ||
            fTagFullLengthTracks || fTagNonBeamT0Tracks || fTagApparentStopper || fTagTopFrontBack);
  }

  void tag(detinfo::DetectorClocksData const& clockData, pma::TrkCandidateColl& tracks);

private:
  size_t outOfDriftWindow(pma::TrkCandidateColl& tracks) const;
  size_t fullHeightCrossing(pma::TrkCandidateColl& tracks) const;
  size_t fullWidthCrossing(pma::TrkCandidateColl& tracks) const;
  size_t fullLengthCrossing(pma::TrkCandidateColl& tracks) const;
  size_t fullCrossingTagger(pma::TrkCandidateColl& tracks, int direction) const;
  size_t nonBeamT0Tag(detinfo::DetectorClocksData const& clockData,
                      pma::TrkCandidateColl& tracks) const;
  size_t tagTopFrontBack(pma::TrkCandidateColl& tracks) const;
  size_t tagApparentStopper(pma::TrkCandidateColl& tracks) const;

  // Convenience functions to see if we have a vertex at the top of the detector
  bool isTopVertex(const TVector3& pos, double tolerance, short int dirIndx) const;
  // or at the front / back walls
  bool isFrontBackVertex(const TVector3& pos, double tolerance, short int dirIndx) const;

  void GetDimensions(); // Use the geometry to get the extent of the detector in x, y and z.
  short int ConvertDirToInt(const TVector3& dir) const; // Is the direction along x, y or z?
  // Tagging parameters
  bool fTagOutOfDriftTracks; // Tag tracks sticking out of 1 drift window.
  double fOutOfDriftMargin;  // Min distance [cm] beyond 1 drift window required
                             // for tagging track as a cosmic background.

  bool fTagFullHeightTracks;  // Tag tracks crossing full height
  bool fTagFullWidthTracks;   // Tag tracks crossing full heightwidth
  bool fTagFullLengthTracks;  // Tag tracks crossing full heightlength
  double fFullCrossingMargin; // Max distance [cm] between track dimension and
                              // detector dimension for crossing tracks.

  bool fTagNonBeamT0Tracks; // Tag tracks that have a reconstructed T0 outside of the beam range.
  double fNonBeamT0Margin;  // Range outside which we should consider events not beam related.

  bool fTagTopFrontBack;
  double fTopFrontBackMargin;

  bool fTagApparentStopper;
  double fApparentStopperMargin;
  bool fVetoActualStopper;
  double fStopperBuffer; // A distance from the end of the track within which we
                         // should find no other track starting. Helps to
                         // prevent identifying broken tracks or interacting
                         // particles as stoppers.

  // The dimensions of the detector from the geometry
  std::vector<double> fDimensionsMin;
  std::vector<double> fDimensionsMax;
};

#endif

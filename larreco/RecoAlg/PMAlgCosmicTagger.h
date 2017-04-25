////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgCosmicTagger
// Author:      ..., R.Sulej (..., robert.sulej@cern.ch) March 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgCosmicTagger_h
#define PMAlgCosmicTagger_h

#include <map>

#include "fhiclcpp/types/Atom.h"
// #include "fhiclcpp/types/Sequence.h"

#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "TVector3.h"

namespace pma
{
    class PMAlgCosmicTagger;
}

class pma::PMAlgCosmicTagger
{

public:

    struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<bool> TagOutOfDriftTracks {
            Name("TagOutOfDriftTracks"),
            Comment("Tag tracks sticking out of 1 drift window.")
        };
        fhicl::Atom<double> OutOfDriftMargin {
            Name("OutOfDriftMargin"),
            Comment("The minimum distance beyond 1 drift window required for tagging track as a cosmic background.")
        };

        fhicl::Atom<bool> TagFullHeightTracks {
            Name("TagFullHeightTracks"),
            Comment("Tag tracks crossing full detector height")
        };
        fhicl::Atom<bool> TagFullWidthTracks {
            Name("TagFullWidthTracks"),
            Comment("Tag tracks crossing full detector width")
        };
        fhicl::Atom<bool> TagFullLengthTracks {
            Name("TagFullLengthTracks"),
            Comment("Tag tracks crossing full detector length")
        };
        fhicl::Atom<double> FullCrossingMargin {
            Name("FullCrossingMargin"),
            Comment("The maximum distance between the track length and detector length for full detector crossing tracks.")
        };
				fhicl::Atom<bool> TagNonBeamT0Tracks {
					Name("TagNonBeamT0Tracks"),
					Comment("Tag particles with reconstructed T0 not consistent with the beam")
				};
				fhicl::Atom<double> NonBeamT0Margin {
					Name("NonBeamT0Margin"),
					Comment("Tag only those events at least <margin> from the beam time")
				};
    };

    PMAlgCosmicTagger(const pma::PMAlgCosmicTagger::Config &config) :
        fTagOutOfDriftTracks(config.TagOutOfDriftTracks()),
        fOutOfDriftMargin(config.OutOfDriftMargin()),
        
        fTagFullHeightTracks(config.TagFullHeightTracks()),
        fTagFullWidthTracks(config.TagFullWidthTracks()),
        fTagFullLengthTracks(config.TagFullLengthTracks()),
        fFullCrossingMargin(config.FullCrossingMargin()),

				fTagNonBeamT0Tracks(config.TagNonBeamT0Tracks()),
				fNonBeamT0Margin(config.NonBeamT0Margin())
    { }

    bool tagAny() const { return (fTagOutOfDriftTracks || fTagFullHeightTracks || fTagFullWidthTracks || fTagFullLengthTracks || fTagNonBeamT0Tracks); }

    void tag(pma::TrkCandidateColl& tracks);

private:

    size_t outOfDriftWindow(pma::TrkCandidateColl& tracks);
		size_t fullHeightCrossing(pma::TrkCandidateColl& tracks);
		size_t fullWidthCrossing(pma::TrkCandidateColl& tracks);
		size_t fullLengthCrossing(pma::TrkCandidateColl& tracks);
		size_t fullCrossingTagger(pma::TrkCandidateColl& tracks, int direction);
		size_t nonBeamT0Tag(pma::TrkCandidateColl& tracks);

		void GetDimensions(); // Use the geometry to get the extent of the detector in x, y and z.
		int ConvertDirToInt(TVector3 &dir); // Is the direction along x, y or z?
    // Tagging parameters
    bool fTagOutOfDriftTracks; // Tag tracks sticking out of 1 drift window.
    double fOutOfDriftMargin;  // Min distance [cm] beyond 1 drift window required for tagging track as a cosmic background.

    bool fTagFullHeightTracks;    // Tag tracks crossing full height
    bool fTagFullWidthTracks;     // Tag tracks crossing full heightwidth
    bool fTagFullLengthTracks;    // Tag tracks crossing full heightlength
    double fFullCrossingMargin;     // Max distance [cm] between track dimension and detector dimension for crossing tracks.

		bool fTagNonBeamT0Tracks;  // Tag tracks that have a reconstructed T0 outside of the beam range.
		double fNonBeamT0Margin;   // Range outside which we should consider events not beam related.

		std::vector<double> fDimensions; // The size of the detector in the x, y and z dimensions.
};

#endif


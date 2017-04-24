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

        fhicl::Atom<bool> TagTopDownTracks {
            Name("TagTopDownTracks"),
            Comment("Tag tracks crossing full Y range (if it is not drift).")
        };
        fhicl::Atom<double> TopDownMargin {
            Name("TopDownMargin"),
            Comment("The maximum distance from top and down required for tagging track as crossing full Y.")
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
        
        fTagTopDownTracks(config.TagTopDownTracks()),
        fTopDownMargin(config.TopDownMargin()),

				fTagNonBeamT0Tracks(config.TagNonBeamT0Tracks()),
				fNonBeamT0Margin(config.NonBeamT0Margin())
    { }

    bool tagAny() const { return (fTagOutOfDriftTracks || fTagTopDownTracks /* || fTagSomething || fTagSomethingElse */); }

    void tag(pma::TrkCandidateColl& tracks);

private:

    size_t outOfDriftWindow(pma::TrkCandidateColl& tracks);
    size_t topDownCrossing(pma::TrkCandidateColl& tracks);
		size_t nonBeamT0Tag(pma::TrkCandidateColl& tracks);

    // Tagging parameters
    bool fTagOutOfDriftTracks; // Tag tracks sticking out of 1 drift window.
    double fOutOfDriftMargin;  // Min distance [cm] beyond 1 drift window required for tagging track as a cosmic background.

    bool fTagTopDownTracks;    // Tag tracks crossing full Y range (if it is not drift).
    double fTopDownMargin;     // Max distance [cm] from top and down required for tagging track as crossing full Y.

		bool fTagNonBeamT0Tracks;  // Tag tracks that have a reconstructed T0 outside of the beam range.
		double fNonBeamT0Margin;   // Range outside which we should consider events not beam related.
};

#endif


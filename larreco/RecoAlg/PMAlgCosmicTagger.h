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
    };

    PMAlgCosmicTagger(const pma::PMAlgCosmicTagger::Config &config) :
        fTagOutOfDriftTracks(config.TagOutOfDriftTracks()),
        fOutOfDriftMargin(config.OutOfDriftMargin()),
        
        fTagTopDownTracks(config.TagTopDownTracks()),
        fTopDownMargin(config.TopDownMargin())
    { }

    bool tagAny() const { return (fTagOutOfDriftTracks || fTagTopDownTracks /* || fTagSomething || fTagSomethingElse */); }

    void tag(pma::TrkCandidateColl& tracks);

private:

    size_t outOfDriftWindow(pma::TrkCandidateColl& tracks);
    size_t topDownCrossing(pma::TrkCandidateColl& tracks);

    // Tagging parameters
    bool fTagOutOfDriftTracks; // Tag tracks sticking out of 1 drift window.
    double fOutOfDriftMargin;  // Min distance [cm] beyond 1 drift window required for tagging track as a cosmic background.

    bool fTagTopDownTracks;    // Tag tracks crossing full Y range (if it is not drift).
    double fTopDownMargin;     // Max distance [cm] from top and down required for tagging track as crossing full Y.
};

#endif


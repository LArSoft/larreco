//////////////////////////////////////////////////////////////////////////////
// \version 
//
// \brief Wrappers for saving MVA results into art::Event
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAWRITER_H
#define ANAB_MVAWRITER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "lardataobj/AnalysisBase/MVAOutput.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/MVAWrapperBase.h"

namespace anab {

/// Helper for registering in the art::EDProducer all data products needed for
/// N-output MVA results: keep MVADescriptions<N> for all types T in one collection
/// while separate instance names are used for the MVA output value collections for
/// each type T.
/// Use one instance of this class per one MVA model, applied to one or more types.
template <size_t N>
class MVAWriter : public MVAWrapperBase {
public:

    /// Name provided to the constructor is used as an instance name for MVADescription<N>
    /// and MVAOutput<N> (combined with the processed data product names); good idea is
    /// to use this name as an indication what MVA was used on the data (like eg. "emtrack"
    /// for outputs from a model distinguishing EM from track-like hits and clusters).
    MVAWriter(const char* name = "") : fInstanceName(name) { }

    friend std::ostream& operator<< (std::ostream &o, MVAWriter const& a)
    {
        o << "MVAWriter for " << a.fInstanceName << std::endl;
        return o;
    }

private:

    std::string fInstanceName;
};

} // namespace anab

//----------------------------------------------------------------------------
// MVAWriter functions.
//

//----------------------------------------------------------------------------
#endif //ANAB_MVAREADER


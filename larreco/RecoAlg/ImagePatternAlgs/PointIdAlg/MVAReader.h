//////////////////////////////////////////////////////////////////////////////
// \version 
//
// \brief Wrappers for accessing MVA results and associated data products
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAREADER_H
#define ANAB_MVAREADER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/MVAWrapperBase.h"

namespace anab {

/// Wrapper for navigation through reconstructed objects of type T and associated
/// N-outputs MVA results with their metadata (this class is not a data product).
template <class T, size_t N>
class MVAReader : public MVAWrapperBase {
public:

    /// Creathe the wrapper for MVA data stored in the event evt with the provided input tag
    /// (the same tag which was used to save MVA results with MVAWriter class).
    MVAReader(const art::Event & evt, const art::InputTag & tag);

    /// Access data product at index "key".
    T const & item(size_t key) const { return (*fDataHandle)[key]; }
    std::vector<T> const & items() const { return *fDataHandle; }

    /// Access the vector of data products.
    std::vector< FeatureVector<N> > const & outputs() const { return *fOutputs; }

    /// Get copy of the MVA output vector at index "key".
    std::array<float, N> getOutput(size_t key) const
    {
        std::array<float, N> vout;
        for (size_t i = 0; i < N; ++i) vout[i] = (*fOutputs)[key][i];
        return vout;
    }

    /// Get copy of the MVA output vector idicated with art::Ptr::key().
    std::array<float, N> getOutput(art::Ptr<T> const & item) const
    { return getOutput(item.key()); }

    /// Get MVA results accumulated over the vector of items (eg. over hits associated to a cluster).
    std::array<float, N> getOutput(std::vector< art::Ptr<T> > const & items) const
    { return pAccumulate(items, *fOutputs); }

    /// Get MVA results accumulated with provided weights over the vector of items
    /// (eg. over clusters associated to a track, weighted by the cluster size; or
    /// over hits associated to a cluster, weighted by the hit area).
    std::array<float, N> getOutput(std::vector< art::Ptr<T> > const & items,
        std::vector<float> const & weights) const
    { return pAccumulate(items, weights, *fOutputs); }

    /// Get MVA results accumulated with provided weighting function over the vector
    /// of items (eg. over clusters associated to a track, weighted by the cluster size;
    /// or over hits associated to a cluster, weighted by the hit area).
    std::array<float, N> getOutput(std::vector< art::Ptr<T> > const & items,
        std::function<float (T const &)> fweight) const
    { return pAccumulate(items, fweight, *fOutputs); }

    /// Get the number of contained items (no. of data product objects equal to no. of MVA output vectors).
    size_t size() const { return fOutputs->size(); }

    /// Get the length of a single vector of MVA values.
    size_t length() const { return N; }

    /// Get the input tag (string representation) of data product used to calculate MVA values.
    const std::string & dataTag() const { return fDescription->dataTag(); }

    /// Access the data product handle.
    const art::Handle< std::vector<T> > & dataHandle() const { return fDataHandle; }

    /// Meaning/name of the index'th column in the collection of MVA output vectors.
    const std::string & outputName(size_t index) const { return fDescription->outputName(index); }

    friend std::ostream& operator<< (std::ostream &o, MVAReader const& a)
    {
        o << "MVAReader:" << std::endl << *(a.fDescription) << std::endl;
        return o;
    }

private:

    MVADescription<N> const * fDescription;
    std::vector< FeatureVector<N> > const * fOutputs;
    art::Handle< std::vector<T> > fDataHandle;
};

} // namespace anab

//----------------------------------------------------------------------------
// MVAReader functions.
//
template <class T, size_t N>
anab::MVAReader<T, N>::MVAReader(const art::Event & evt, const art::InputTag & tag) :
    fDescription(0)
{
    if (!N) { throw cet::exception("MVAReader") << "MVA size should be > 0." << std::endl; }

    auto descriptionHandle = evt.getValidHandle< std::vector< anab::MVADescription<N> > >(tag);

    // search for MVADescription<N> produced for the type T, with the instance name from the tag
    std::string outputInstanceName = tag.instance() + getProductName(typeid(T));
    for (auto const & dscr : *descriptionHandle)
    {
        if (dscr.outputInstance() == outputInstanceName)
        {
            fDescription = &dscr; break;
        }
    }
    if (!fDescription) { throw cet::exception("MVAReader") << "MVA description not found for " << outputInstanceName << std::endl; }

    fOutputs = &*(evt.getValidHandle< std::vector< FeatureVector<N> > >( art::InputTag(tag.label(), fDescription->outputInstance(), tag.process()) ));

    if (!evt.getByLabel( fDescription->dataTag(), fDataHandle ))
    {
        throw cet::exception("MVAReader") << "Associated data product handle failed: " << *(fDataHandle.whyFailed()) << std::endl;
    }

    if (fOutputs->size() != fDataHandle->size())
    {
        throw cet::exception("MVAReader") << "MVA outputs and data products sizes inconsistent: " << fOutputs->size() << "!=" << fDataHandle->size() << std::endl;
    }
}
//----------------------------------------------------------------------------

#endif //ANAB_MVAREADER


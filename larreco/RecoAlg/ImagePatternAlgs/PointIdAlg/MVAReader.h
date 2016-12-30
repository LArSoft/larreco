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
    T const & item(size_t key) const { return (*fDataProd)[key]; }

    /// Get copy of the MVA output vector at index "key".
    std::array<float, N> getOutput(size_t key) const
    {
        std::array<float, N> vout;
        for (size_t i = 0; i < N; ++i) vout[i] = (*fOutputs)[key][i];
        return vout;
    }

    /// Get copy of the MVA output vector idicated with art::Ptr::key().
    std::array<float, N> getOutput(art::Ptr<T> const & item) const
    {
        std::array<float, N> vout;
        for (size_t i = 0; i < N; ++i) vout[i] = (*fOutputs)[item.key()][i];
        return vout;
    }

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

    /// Get the number of contained items (no. of data product objjects equal to no. of MVA output vectors).
    size_t size() const { return fOutputs->size(); }

    /// Get the length of a single vector of MVA values.
    size_t length() const { return N; }

    const std::string & dataTag() const { return fDescription->dataTag(); }

    const std::string & outputName(size_t index) const { return fDescription->outputName(index); }

    friend std::ostream& operator<< (std::ostream &o, MVAReader const& a)
    {
        o << "MVAReader:" << std::endl << *(a.fDescription) << std::endl;
        return o;
    }

private:

    MVADescription<N> const * fDescription;
    std::vector< MVAOutput<N> > const * fOutputs;
    std::vector<T> const * fDataProd;
};

} // namespace anab

//----------------------------------------------------------------------------
// MVAReader functions.
//
template <class T, size_t N>
anab::MVAReader<T, N>::MVAReader(const art::Event & evt, const art::InputTag & tag) :
    fDescription(0),
    fDataProd(0)
{
    if (!N) { throw cet::exception("MVAReader") << "MVA size should be > 0."; }

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
    if (!fDescription) { throw cet::exception("MVAReader") << "MVA description not found."; }

    fOutputs = &*(evt.getValidHandle< std::vector< MVAOutput<N> > >( art::InputTag(tag.label(), fDescription->outputInstance(), tag.process()) ));

    fDataProd = &*(evt.getValidHandle< std::vector<T> >( fDescription->dataTag() ));

    if (fOutputs->size() != fDataProd->size()) { throw cet::exception("MVAReader") << "MVA outputs and data products sizes inconsistent."; }
}
//----------------------------------------------------------------------------

#endif //ANAB_MVAREADER


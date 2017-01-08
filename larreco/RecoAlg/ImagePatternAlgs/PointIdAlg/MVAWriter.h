//////////////////////////////////////////////////////////////////////////////
// \version 
//
// \brief Wrapper for saving MVA results into art::Event
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAWRITER_H
#define ANAB_MVAWRITER_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/MVAWrapperBase.h"

namespace anab {

/// Index to the MVA output collection, used when result vectors are added or set.
typedef size_t MVAOutput_ID;

/// Helper for registering in the art::EDProducer all data products needed for
/// N-output MVA results: keep MVADescriptions<N> for all types T in one collection
/// while separate instance names are used for the MVA output value collections for
/// each type T.
/// Use one instance of this class per one MVA model, applied to one or more types.
template <size_t N>
class MVAWriter : public MVAWrapperBase {
public:

    /// Name provided to the constructor is used as an instance name for MVADescription<N>
    /// and FeatureVector<N> (for which it is combined with the processed data product names).
    /// Good idea is to use the name as an indication of what MVA model was used on the data
    /// (like eg. "emtrack" for outputs from a model distinguishing EM from track-like hits
    /// and clusters). The name is used as an instance name for the MVADescription data product
    /// which lets you to save multiple MVA results from a single art module.
    MVAWriter(art::EDProducer* module, const char* name = "") :
        fProducer(module), fInstanceName(name),
        fIsDescriptionRegistered(false),
        fDescriptions(nullptr)
    { }

    /// Register the collection of metadata type MVADescription<N> (once for all data types
    /// for which MVA is saved) and the collection of FeatureVectors<N> (using data type name
    /// added to fInstanceName as instance name of the collection made for the type T).
    template <class T>
    void produces_using();

    /// Initialize container for MVA outputs and, if not yet done, the container for
    /// metadata, then creates metadata for data products of type T. MVA output container
    /// is initialized to hold dataSize vectors (if dataSize > 0): use setOutput() to
    /// store values.
    /// Returns index of outputs which should be used when saving actual output values.
    template <class T>
    MVAOutput_ID initOutputs(std::string const & dataTag, size_t dataSize,
        std::vector< std::string > const & names = std::vector< std::string >(N, ""));

    template <class T>
    MVAOutput_ID initOutputs(art::InputTag const & dataTag, size_t dataSize,
        std::vector< std::string > const & names = std::vector< std::string >(N, ""))
    { return initOutputs<T>(dataTag.encode(), dataSize, names); }

    void setOutput(MVAOutput_ID id, size_t key, std::array<float, N> const & values) { (*(fOutputs[id]))[key] = values; }
    void setOutput(MVAOutput_ID id, size_t key, std::array<double, N> const & values) { (*(fOutputs[id]))[key] = values; }
    void setOutput(MVAOutput_ID id, size_t key, std::vector<float> const & values) { (*(fOutputs[id]))[key] = values; }
    void setOutput(MVAOutput_ID id, size_t key, std::vector<double> const & values) { (*(fOutputs[id]))[key] = values; }


    /// Initialize container for MVA outputs and, if not yet done, the container for
    /// metadata, then creates metadata for data products of type T. MVA output container
    /// is initialized as EMPTY and MVA vectors should be added with addOutput() function.
    /// Returns index of outputs which should be used when adding actual output values.
    template <class T>
    MVAOutput_ID initOutputs(art::InputTag const & dataTag,
        std::vector< std::string > const & names = std::vector< std::string >(N, ""))
    { return initOutputs<T>(dataTag.encode(), 0, names); }

    template <class T>
    MVAOutput_ID initOutputs(
        std::vector< std::string > const & names = std::vector< std::string >(N, ""))
    { return initOutputs<T>(std::string(""), 0, names); }

    void addOutput(MVAOutput_ID id, std::array<float, N> const & values) { fOutputs[id]->emplace_back(values); }
    void addOutput(MVAOutput_ID id, std::array<double, N> const & values) { fOutputs[id]->emplace_back(values); }
    void addOutput(MVAOutput_ID id, std::vector<float> const & values) { fOutputs[id]->emplace_back(values); }
    void addOutput(MVAOutput_ID id, std::vector<double> const & values) { fOutputs[id]->emplace_back(values); }

    /// Set tag of associated data products in case it was not ready at the initialization time.
    void setDataTag(MVAOutput_ID id, art::InputTag const & dataTag) { (*fDescriptions)[id].setDataTag(dataTag.encode()); }

    /// Check consistency and save all the results in the event.
    void saveOutputs(art::Event & evt);

    /// Get MVA results accumulated over the vector of items (eg. over hits associated to a cluster).
    /// NOTE: MVA outputs for these items has to be added to the MVAWriter first!
    template <class T>
    std::array<float, N> getOutput(std::vector< art::Ptr<T> > const & items) const
    { return pAccumulate<T, N>(items, *(fOutputs[getProductID<T>()])); }

    /// Get MVA results accumulated with provided weights over the vector of items
    /// (eg. over clusters associated to a track, weighted by the cluster size; or
    /// over hits associated to a cluster, weighted by the hit area).
    /// NOTE: MVA outputs for these items has to be added to the MVAWriter first!
    template <class T>
    std::array<float, N> getOutput(std::vector< art::Ptr<T> > const & items,
        std::vector<float> const & weights) const
    { return pAccumulate<T, N>(items, weights, *(fOutputs[getProductID<T>()])); }

    /// Get MVA results accumulated with provided weighting function over the vector
    /// of items (eg. over clusters associated to a track, weighted by the cluster size;
    /// or over hits associated to a cluster, weighted by the hit area).
    /// NOTE: MVA outputs for these items has to be added to the MVAWriter first!
    template <class T>
    std::array<float, N> getOutput(std::vector< art::Ptr<T> > const & items,
        std::function<float (T const &)> fweight) const
    { return pAccumulate<T, N>(items, fweight, *(fOutputs[getProductID<T>()])); }

    template <class T>
    std::array<float, N> getOutput(std::vector< art::Ptr<T> > const & items,
        std::function<float (art::Ptr<T> const &)> fweight) const
    { return pAccumulate<T, N>(items, fweight, *(fOutputs[getProductID<T>()])); }

    /// Get copy of the MVA output vector for the type T, at index "key".
    template <class T>
    std::array<float, N> getOutput(size_t key) const
    {
        std::array<float, N> vout;
        auto const & src = ( *(fOutputs[getProductID<T>()]) )[key];
        for (size_t i = 0; i < N; ++i) vout[i] = src[i];
        return vout;
    }

    /// Get copy of the MVA output vector for the type T, idicated with art::Ptr::key().
    template <class T>
    std::array<float, N> getOutput(art::Ptr<T> const & item) const
    {
        std::array<float, N> vout;
        auto const & src = ( *(fOutputs[getProductID<T>()]) )[item.key()];
        for (size_t i = 0; i < N; ++i) vout[i] = src[i];
        return vout;
    }

    /// Get the number of contained MVA output vectors.
    size_t size(MVAOutput_ID id) const { return fOutputs[id]->size(); }

    /// Get the length of a single vector of MVA values.
    size_t length() const { return N; }

    friend std::ostream& operator<< (std::ostream &o, MVAWriter const& a)
    {
        o << "MVAWriter for " << a.fInstanceName << ", " << N << " outputs";
        if (!a.fRegisteredDataTypes.empty())
        {
            o << ", ready to write results made for:" << std::endl;
            for (auto const & n : a.fRegisteredDataTypes) { o << "\t" << n << std::endl; }
        }
        else { o << ", nothing registered for writing to the events" << std::endl; }
        return o;
    }

private:

    // Data initialized for the module life:
    art::EDProducer* fProducer;
    std::string fInstanceName;

    std::vector< std::string > fRegisteredDataTypes;
    bool fIsDescriptionRegistered;

    // Data collected for each event:
    std::unordered_map< size_t, MVAOutput_ID > fTypeHashToID;
    template <class T> MVAOutput_ID getProductID() const;

    std::vector< std::unique_ptr< std::vector< anab::FeatureVector<N> > > > fOutputs;
    std::unique_ptr< std::vector< anab::MVADescription<N> > > fDescriptions;
    void clearEventData()
    {
        fTypeHashToID.clear(); fOutputs.clear();
        fDescriptions.reset(nullptr);
    }

    /// Check if the the writer is configured to write results for data product type name.
    bool dataTypeRegistered(const std::string & dname) const;
    /// Check if the containers for results prepared for "tname" data type are ready.
    bool descriptionExists(const std::string & tname) const;
};

} // namespace anab

//----------------------------------------------------------------------------
// MVAWriter functions.
//
template <size_t N>
template <class T>
anab::MVAOutput_ID anab::MVAWriter<N>::getProductID() const
{
    auto const & ti = typeid(T);
    auto search = fTypeHashToID.find(getProductHash(ti));
    if (search != fTypeHashToID.end()) { return search->second; }
    else
    {
        throw cet::exception("MVAWriter") << "MVA not initialized for product " << getProductName(ti) << std::endl;
    }
}
//----------------------------------------------------------------------------

template <size_t N>
bool anab::MVAWriter<N>::dataTypeRegistered(const std::string & dname) const
{
    for (auto const & s : fRegisteredDataTypes)
    {
        if (s == dname) { return true; }
    }
    return false;
}
//----------------------------------------------------------------------------

template <size_t N>
template <class T>
void anab::MVAWriter<N>::produces_using()
{
    std::string dataName = getProductName(typeid(T));
    if (dataTypeRegistered(dataName))
    {
        throw cet::exception("MVAWriter") << "Type " << dataName << "was already registered." << std::endl;
    }

    if (!fIsDescriptionRegistered)
    {
        fProducer->produces< std::vector< anab::MVADescription<N> > >(fInstanceName);
        fIsDescriptionRegistered = true;
    }

    fProducer->produces< std::vector< anab::FeatureVector<N> > >(fInstanceName + dataName);
    fRegisteredDataTypes.push_back(dataName);
}
//----------------------------------------------------------------------------

template <size_t N>
bool anab::MVAWriter<N>::descriptionExists(const std::string & tname) const
{
    if (!fDescriptions) return false;

    std::string n = fInstanceName + tname;
    for (auto const & d : *fDescriptions)
    {
        if (d.outputInstance() == n) { return true; }
    }
    return false;
}
//----------------------------------------------------------------------------

template <size_t N>
template <class T>
anab::MVAOutput_ID anab::MVAWriter<N>::initOutputs(
    std::string const & dataTag, size_t dataSize,
    std::vector< std::string > const & names)
{
    size_t dataHash = getProductHash(typeid(T));
    std::string dataName = getProductName(typeid(T));

    if (!dataTypeRegistered(dataName))
    {
        throw cet::exception("MVAWriter") << "Type " << dataName << "not registered with produces_using() function." << std::endl;
    }

    if (!fDescriptions)
    {
        fDescriptions = std::make_unique< std::vector< anab::MVADescription<N> > >();
    }
    else if (descriptionExists(dataName))
    {
        throw cet::exception("MVAWriter") << "MVADescription<N> already initialized for " << dataName << std::endl;
    }
    fDescriptions->emplace_back(dataTag, fInstanceName + dataName);

    fOutputs.push_back( std::make_unique< std::vector< anab::FeatureVector<N> > >() );
    anab::MVAOutput_ID id = fOutputs.size() - 1;
    fTypeHashToID[dataHash] = id;

    if (dataSize) { fOutputs[id]->resize(dataSize, anab::FeatureVector<N>(0.0F)); }

    return id;
}
//----------------------------------------------------------------------------

template <size_t N>
void anab::MVAWriter<N>::saveOutputs(art::Event & evt)
{
    for (auto const & n : fRegisteredDataTypes)
    {
        if (!descriptionExists(n))
        {
            throw cet::exception("MVAWriter") << "No MVADescription<N> prepared for type " << n << std::endl;
        }
    }

    if (fOutputs.size() != fDescriptions->size())
    {
        throw cet::exception("MVAWriter") << "MVADescription<N> vector length not equal to the number of FeatureVector<N> vectors" << std::endl;
    }

    for (size_t i = 0; i < fOutputs.size(); ++i)
    {
        auto const & outInstName = (*fDescriptions)[i].outputInstance();
        if ((*fDescriptions)[i].dataTag().empty())
        {
            throw cet::exception("MVAWriter") << "MVADescription<N> reco data tag not set for " << outInstName << std::endl;
        }
        evt.put(std::move(fOutputs[i]), outInstName);
    }
    evt.put(std::move(fDescriptions), fInstanceName);
    clearEventData();
}
//----------------------------------------------------------------------------

#endif //ANAB_MVAREADER


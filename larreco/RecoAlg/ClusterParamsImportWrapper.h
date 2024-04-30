/** ****************************************************************************
 * @file   ClusterParamsImportWrapper.h
 * @brief  Wrapper for ClusterParamsAlgBase objects to accept arbitrary input
 * @author petrillo@fnal.gov
 * @date   January 22, 2015
 * @see    ClusterParamsAlgBase.h
 *
 * ****************************************************************************/

#ifndef CLUSTERPARAMSARTWRAPPER_H
#define CLUSTERPARAMSARTWRAPPER_H

// C/C++ standard libraries
#include <stdexcept> // std::logic_error
#include <utility>   // std::forward()

// LArSoft libraries
#include "lardata/Utilities/Dereference.h" // util::make_pointer()

namespace util {
  class GeometryUtilities;
}

namespace recob {
  class Hit;
}

/// Cluster reconstruction namespace
namespace cluster {

  /**
   * @brief Wrapper for ClusterParamsAlgBase objects to accept diverse input
   * @tparam Algo the ClusterParamsAlgBase-derived class to be wrapped
   * @see ClusterParamsAlgBase
   *
   * This simple wrapper class adds a non-virtual ImportHits() method that can
   * import Hits from different formats than std::vector<recob::Hit const*>.
   *
   * This also allows the algorithms derived from ClusterParamsAlgBase to stay
   * framework-agnostic.
   *
   */
  template <class Algo>
  class ClusterParamsImportWrapper : public Algo {
  public:
    using ClusterParamsAlg_t = Algo; ///< type of wrapped class

    /// Constructor: just forwards all the stuff to the wrapped class
    template <typename... Args>
    ClusterParamsImportWrapper(Args... args) : ClusterParamsAlg_t(std::forward<Args>(args)...)
    {}

    /// @{
    /// @name Hit import functions
    ///
    /// Methods to import hits int the algorithm.
    ///

    /**
     * @brief Calls SetHits() with the hits in the sequence
     * @tparam Iter type of iterator to source hits
     * @param begin iterator to the first hit source
     * @param end iterator to after-the-last hit source
     *
     * The type in the sequence should contain either recob::Hit or some sort
     * of pointer to it.
     */
    template <typename Iter>
    void ImportHits(util::GeometryUtilities const& gser, Iter begin, Iter end)
    {
      std::vector<recob::Hit const*> hits;
      std::transform(begin, end, std::back_inserter(hits), [](auto value) {
        return lar::util::make_pointer(value);
      });
      ClusterParamsAlg_t::SetHitsFromPointers(gser, hits);
    } // ImportHits()

    /**
     * @brief Calls SetHits() with the result of converted hits
     * @tparam Iter type of iterator to source hits
     * @tparam Convert type of operation to convert to recob::Hit const*
     * @param begin iterator to the first hit source
     * @param end iterator to after-the-last hit source
     * @param converter predicate to convert the pointed values to recob::Hit
     *
     * The converter should respect either of the forms:
     *
     *     recob::Hit converter(typename Iter::value_type)
     *     recob::Hit const* converter(typename Iter::value_type)
     *
     * The hit produced by the converter will be moved into a vector, and the
     * complete vector will be used to initialize the algorithm.
     */
    template <typename Iter, typename Convert>
    void ImportHits(Iter begin, Iter end, Convert converter)
    {
      std::vector<recob::Hit const*> hits;
      std::transform(begin, end, std::back_inserter(hits), [converter](auto value) {
        return lar::util::make_pointer(converter(value));
      });
      ClusterParamsAlg_t::SetHits(hits);
    } // ImportHits()

    /**
     * @brief Calls SetHits() with the hits in the sequence
     * @tparam Cont type of container to source hits
     * @param cont container of source hits
     *
     * The type in the container should contain either recob::Hit or some sort
     * of pointer to it.
     */
    template <typename Cont>
    void ImportHits(util::GeometryUtilities const& gser, Cont cont)
    {
      ImportHits(gser, std::begin(cont), std::end(cont));
    }

    /**
     * @brief Calls SetHits() with the result of converted hits
     * @tparam Cont type of container to source hits
     * @tparam Convert type of operation to convert to recob::Hit const*
     * @param cont container of source hits
     * @param converter predicate to convert the pointed values to recob::Hit
     *
     * The converter should respect either of the forms:
     *
     *     recob::Hit converter(typename Iter::value_type)
     *     recob::Hit const* converter(typename Iter::value_type)
     *
     * The hit produced by the converter will be moved into a vector, and the
     * complete vector will be used to initialize the algorithm.
     */
    template <typename Cont, typename Convert>
    void ImportHits(util::GeometryUtilities const& gser, Cont cont, Convert converter)
    {
      ImportHits(gser, std::begin(cont), std::end(cont), converter);
    }

    /// @}

  }; //class ClusterParamsImportWrapper

} //namespace cluster

#endif // CLUSTERPARAMSARTWRAPPER_H

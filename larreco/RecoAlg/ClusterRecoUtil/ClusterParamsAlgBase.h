/** ****************************************************************************
 * @file   ClusterParamsAlgBase.h
 * @brief  Interface for a algorithm class computing cluster parameters
 * @author petrillo@fnal.gov
 * @date   January 21, 2015
 * @see    ClusterParamsAlg.h StandardClusterParamsAlg.h
 *
 * Changes:
 * 20150121 Gianluca Petrillo (petrillo@fnal.gov)
 *   structure shaped on ClusterParamsAlg class by Andrzej Szelc
 *
 * ****************************************************************************/

#ifndef CLUSTERPARAMSALGBASE_H
#define CLUSTERPARAMSALGBASE_H

// C/C++ standard libraries
#include <algorithm> // std::transform()
#include <stdexcept> // std::logic_error
#include <vector>

// LArSoft libraries
#include "lardataobj/RecoBase/Hit.h"

/// Cluster reconstruction namespace
namespace cluster {

  /// Implementation details of cluster namespace
  namespace details {

    /// Type for a simple measurement: value and error
    template <typename T>
    class Measure_t : public std::pair<T, T> {
      using Base_t = std::pair<T, T>;

    public:
      using Data_t = T;

      /// Default constructor: initializes to 0
      Measure_t() : Base_t(Data_t(0), Data_t(0)) {}

      /// Constructor: initializes to the specified value, error is 0
      Measure_t(Data_t value) : Base_t(value, Data_t(0)) {}

      /// Constructor: initializes to the specified value and error
      Measure_t(Data_t value, Data_t error) : Base_t(value, error) {}

      Data_t
      value() const
      {
        return Base_t::first;
      }
      Data_t&
      value()
      {
        return Base_t::first;
      }

      Data_t
      error() const
      {
        return Base_t::second;
      }
      Data_t&
      error()
      {
        return Base_t::second;
      }
    }; // Measure_t

  } // namespace details

  /**
   * @brief Algorithm collection class computing cluster parameters
   * @see ClusterParamsAlg
   *
   * This class is an interface only.
   * The implementing classes should implement constructors able to read the
   * necessary information from hit lists, and any of the virtual algorithms.
   *
   * The interface allows for the single extraction of the information required
   * in recob::Cluster version 14.
   * The implementation can compute and cache different variables at once.
   * The accessor functions are non-const, allowing the caching of the quantity
   * after computation (without using mutable cache members).
   *
   * The default algorithm functions throw an exception.
   *
   * The algorithms require a list of hits as input. The structure chosen for
   * this interface is recob::Hit from LArSoft, since it is complete by
   * definition and it does not carry deep dependences (in fact, recob::Hit
   * version 14 has only larreco/SimpleTypesAndConstants.h as compile-time
   * dependency, and no external link-time dependency beside standard C++).
   *
   */
  class ClusterParamsAlgBase {
  public:
    /// Type used to return values with errors
    using Measure_t = details::Measure_t<float>;

    // no constructor provided for the abstract class

    /// Virtual destructor. Override freely.
    virtual ~ClusterParamsAlgBase() {}

    /**
     * @brief Restores the class to post-configuration, pre-initialization state
     *
     * This function is expected to be called before SetHits(), and the
     * implementation might call it in SetHits() itself.
     */
    virtual void
    Clear()
    {}

    /**
     * @brief Sets the list of input hits
     * @param hits list of pointers to hits
     * @throw undefined in case of error, this method can throw (anything)
     *
     * The hits are in the LArSoft format of recob::Hits, that should have
     * enough information for all the algorithms.
     * The hits are passed as constant pointers. The implementation is expected
     * to either copy the vector (not just to keep a reference to it, since
     * the vector might be temporary) or to translated the required information
     * from the hits into its own internal format.
     * The hits are expected to exist as long as this object is used, until
     * the next Clear() call.
     * It is recommended that this method call Clear() at the beginning.
     * This is left to the implementation, that might still implement a
     * different strategy.
     */
    virtual void SetHits(std::vector<recob::Hit const*> const& hits) = 0;

    /**
     * @brief Sets the list of input hits
     * @param hits list of hits (hits will not be modified)
     * @throw undefined in case of error, this method can throw (anything)
     *
     * The same general directions apply as for SetHits() version with pointers.
     * This version takes a list of recob::Hit, rather than their pointers.
     * It can simplify upstream handling when the original list is not
     * recob::Hit and creation of temporary recob::Hit is needed.
     * In that case, managing to obtain pointers to these temporary objects can
     * be inefficient.
     *
     * The default implementation provided here is not efficient either, since
     * it just creates an additional vector of recob::Hit pointers and uses it.
     * If an implementation is concerned with efficiency, it can reimplement
     * this to initialize the algorithm in a more direct way.
     */
    virtual void
    SetHits(std::vector<recob::Hit> const& hits)
    {
      std::vector<recob::Hit const*> hitptrs;
      hitptrs.reserve(hits.size());
      std::transform(hits.begin(),
                     hits.end(),
                     std::back_inserter(hitptrs),
                     [](recob::Hit const& hit) { return &hit; });
      SetHits(hitptrs);
    }

    /// Set the verbosity level
    virtual void
    SetVerbose(int level = 1)
    {
      verbose = level;
    }

    //@{
    /**
     * @brief Computes the charge on the first and last wire of the track
     * @return the charge in ADC counts, with uncertainty
     */
    virtual Measure_t
    StartCharge()
    {
      throw NotImplemented(__func__);
    }
    virtual Measure_t
    EndCharge()
    {
      throw NotImplemented(__func__);
    }
    //@}

    //@{
    /**
     * @brief Computes the angle at the start or end of the cluster
     * @return angle at the start of the cluster, in radians
     */
    virtual Measure_t
    StartAngle()
    {
      throw NotImplemented(__func__);
    }
    virtual Measure_t
    EndAngle()
    {
      throw NotImplemented(__func__);
    }
    //@}

    //@{
    /**
     * @brief Computes the opening angle at the start or end of the cluster
     * @return angle at the start of the cluster, in radians
     */
    virtual Measure_t
    StartOpeningAngle()
    {
      throw NotImplemented(__func__);
    }
    virtual Measure_t
    EndOpeningAngle()
    {
      throw NotImplemented(__func__);
    }
    //@}

    /// @name Cluster charge
    /// @{
    /**
     * @brief Computes the total charge of the cluster from Hit::Integral()
     * @return total charge of the cluster, in ADC count units
     * @see IntegralStdDev(), SummedADC()
     */
    virtual Measure_t
    Integral()
    {
      throw NotImplemented(__func__);
    }

    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see Integral()
     *
     * Hit charge is obtained by recob::Hit::Integral().
     */
    virtual Measure_t
    IntegralStdDev()
    {
      throw NotImplemented(__func__);
    }

    /**
     * @brief Computes the total charge of the cluster from Hit::SummedADC()
     * @return total charge of the cluster, in ADC count units
     * @see SummedADCStdDev(), Integral()
     */
    virtual Measure_t
    SummedADC()
    {
      throw NotImplemented(__func__);
    }

    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see SummedADC()
     *
     * Hit charge is obtained by recob::Hit::SummedADC().
     */
    virtual Measure_t
    SummedADCStdDev()
    {
      throw NotImplemented(__func__);
    }

    /// @}

    /// Returns the number of hits in the cluster
    virtual size_t
    NHits()
    {
      throw NotImplemented(__func__);
    }

    /**
     * @brief Fraction of wires in the cluster with more than one hit
     * @return fraction of wires with more than one hit, or 0 if no wires
     *
     * Returns a quantity defined as NMultiHitWires / NWires,
     * where NWires is the number of wires hosting at least one hit of this
     * cluster, and NMultiHitWires is the number of wires which have more
     * than just one hit.
     */
    virtual float
    MultipleHitDensity()
    {
      throw NotImplemented(__func__);
    }

    /**
     * @brief Computes the width of the cluster
     * @return width of the cluster
     *
     */
    virtual float
    Width()
    {
      throw NotImplemented(__func__);
    }

  protected:
    int verbose = 0; ///< verbosity level: 0 is normal, negative is even quieter

    static std::logic_error
    NotImplemented(std::string function_name)
    {
      return std::logic_error(function_name + "() not implemented.");
    }

  }; //class ClusterParamsAlgBase

} //namespace cluster

#endif // CLUSTERPARAMSALGBASE_H

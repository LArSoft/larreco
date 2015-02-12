/** ****************************************************************************
 * @file   OverriddenClusterParamsAlg.h
 * @brief  Overrides another ClusterParamsAlgBase class with selected constants
 * @author petrillo@fnal.gov
 * @date   February 3, 2015
 * @see    StandardClusterParamsAlg.cxx
 * 
 * ****************************************************************************/

#ifndef OVERRIDDENCLUSTERPARAMSALG_H
#define OVERRIDDENCLUSTERPARAMSALG_H

// C/C++ standard library
#include <vector>
#include <bitset>
#include <type_traits> // std::is_base_of<>
#include <utility> // std::forward<>

// LArSoft libraries
#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"
#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"



/// Cluster reconstruction namespace
namespace cluster {
  
  namespace details {
    /// Class holding a value of one among some selected types...
    class MultiValue {
        public:
      using Measure_t = cluster::details::Measure_t<float>;
      
      // I can't believe I am writing this shit
      union {
        Measure_t measure_value;
        float     float_value;
        size_t    size_t_value;
      };
      
      // set all the default stuff
      MultiValue(MultiValue const&) = default;
      MultiValue(MultiValue&&) = default;
      MultiValue& operator= (MultiValue const&) = default;
      MultiValue& operator= (MultiValue &&) = default;
      ~MultiValue() = default;
      
      /// Sets the value from a value of type T; undefined by default
      template <typename T>
      MultiValue& operator= (T);
      
      /// Converts the value to type T; undefined by default
      template <typename T>
      operator T() const;
      
    }; // MultiValue
    
  } // namespace details
  
  
  /** **************************************************************************
   * @brief Algorithm collection class computing cluster parameters
   * @tparam AlgoBase class of algorithms to be overridden
   * @see ClusterParamsAlg
   *
   * This class wraps a ClusterParamsAlgBase class, and overrides selected
   * methods with constant values.
   * The same effect can be obtained explicitly creating a new class with
   * the proper methods overridden. This one is a more convenient way to get
   * the same result, but it's slower and less flexible.
   */
  template <typename AlgoBase>
  class OverriddenClusterParamsAlg: public AlgoBase {
    static_assert(std::is_base_of<ClusterParamsAlgBase, AlgoBase>::value,
      "OverriddenClusterParamsAlg template parameter must derive"
      " from ClusterParamsAlgBase"
      );
      public:
    using Base_t = AlgoBase;
    using This_t = OverriddenClusterParamsAlg<AlgoBase>;
    using Measure_t = typename AlgoBase::Measure_t;
    
    typedef enum {
      cpStartAngle,          ///< StartAngle()
      cpEndAngle,            ///< EndAngle()
      cpStartCharge,         ///< StartCharge()
      cpEndCharge,           ///< EndCharge()
      cpStartOpeningAngle,   ///< StartOpeningAngle()
      cpEndOpeningAngle,     ///< EndOpeningAngle()
      cpIntegral,            ///< Integral()
      cpIntegralStdDev,      ///< IntegralStdDev()
      cpSummedADC,           ///< SummedADC()
      cpSummedADCStdDev,     ///< SummedADCStdDev()
      cpNHits,               ///< NHits()
      cpMultipleHitDensity,    ///< MultipleHitDensity()
      cpWidth,               ///< Width()
      NParameters            ///< total number of supported parameters
    } ParameterType_t; ///< type of cluster parameters
    
    /// Constructor; just forwards the arguments to the base class
    template <typename... Args>
    explicit OverriddenClusterParamsAlg(Args&&... args):
      Base_t(std::forward<Args>(args)...)
      {}
    
    /// Destructor
    virtual ~OverriddenClusterParamsAlg() = default;
    
    
    /**
     * @brief Overrides the specified cluster parameter
     * @param param which cluster parameter to override
     * @param value the value of the cluster parameter to be returned
     * @return this object
     * @see ReleaseParameter()
     *
     * For parameters without uncertainty, the uncertainty will be ignored.
     */
    This_t& OverrideParameter(ParameterType_t param, Measure_t value)
      {
        overridden_set.set((size_t) param);
        values[(size_t) param] = value;
        return *this;
      } // OverrideParameter()
    
    /**
     * @brief Cancels the override of the specified cluster parameter
     * @param param which cluster parameter not to override any more
     * @return this object
     * @see OverrideParameter()
     */
    This_t& ReleaseParameter(ParameterType_t param)
      { overridden_set.set((size_t) param); return *this; }
    
    
    /// Returns whether the specified parameter is currently overridden
    bool isOverridden(ParameterType_t param) const
      { return overridden_set.test((size_t) param); }
    
    
    /// @{
    /// @name Algorithm results
    
    //@{
    /**
     * @brief Computes the charge on the first and last wire of the track
     * @return the charge in ADC counts, with uncertainty
     * 
     * The implementation in ClusterParamsAlg provides an estimation of the
     * charge collected in the first or last 1 cm of the cluster, using a linear
     * fit on the deposited charge to reduce fluctuations.
     */
    virtual Measure_t StartCharge() override
      { return ReturnValue(cpStartCharge, &Base_t::StartCharge); }
    virtual Measure_t EndCharge() override
      { return ReturnValue(cpEndCharge, &Base_t::EndCharge); }
    //@}
    
    //@{
    /**
     * @brief Computes the angle of the cluster
     * @return angle of the cluster in the wire x time space, in radians
     * 
     * Uses the coordinates from the hits, weighted by charge (Hit::Integral())
     * to compute a slope in the homogenized wire x time space.
     * The homogenized space has both wires and ticks converted into a distance
     * (by using detector parameters: wire pitch and drift velocity).
     * 
     * The angle is in the @f$ [ -\pi, \pi ] @f$ range, with 0 corresponding to
     * a cluster parallel to the wire plane and @f$ \pi @f$ to a cluster
     * orthogonal to the wire plane, going farther from it.
     * 
     * @note Both the methods return the same value.
     */
    virtual Measure_t StartAngle() override
      { return ReturnValue(cpStartAngle, &Base_t::StartAngle); }
    
    virtual Measure_t EndAngle() override
      { return ReturnValue(cpEndAngle, &Base_t::EndAngle); }
    //@}
    
    //@{
    /**
     * @brief Computes the opening angle at the start or end of the cluster
     * @return angle at the start of the cluster, in radians
     * 
     * This algorithm returns an opening angle after weighting the hits by
     * their charge (as defined bu Hit::Integral());
     */
    virtual Measure_t StartOpeningAngle() override
      { return ReturnValue(cpStartOpeningAngle, &Base_t::StartOpeningAngle); }
    virtual Measure_t EndOpeningAngle() override
      { return ReturnValue(cpEndOpeningAngle, &Base_t::EndOpeningAngle); }
    //@}
    
    
    /// @name Cluster charge
    /// @{
    /**
     * @brief Computes the total charge of the cluster from Hit::Integral()
     * @return total charge of the cluster, in ADC count units
     * @see IntegralStdDev(), SummedADC()
     * 
     * ClusterParamsAlg computes the sum from all hits.
     */
    virtual Measure_t Integral() override
      { return ReturnValue(cpIntegral, &Base_t::Integral); }
    
    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see Integral()
     * 
     * ClusterParamsAlg computes the standard deviation of the sample of charges
     * from all hits.
     * Hit charge is obtained by recob::Hit::Integral().
     */
    virtual Measure_t IntegralStdDev() override
      { return ReturnValue(cpIntegralStdDev, &Base_t::IntegralStdDev); }
    
    /**
     * @brief Computes the total charge of the cluster from Hit::SummedADC()
     * @return total charge of the cluster, in ADC count units
     * @see SummedADCStdDev(), Integral()
     * 
     * ClusterParamsAlg computes the sum from all hits.
     */
    virtual Measure_t SummedADC() override
      { return ReturnValue(cpSummedADC, &Base_t::SummedADC); }
    
    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see SummedADC()
     * 
     * ClusterParamsAlg computes the standard deviation of the sample of charges
     * from all hits.
     * Hit charge is obtained by recob::Hit::SummedADC().
     */
    virtual Measure_t SummedADCStdDev() override
      { return ReturnValue(cpSummedADCStdDev, &Base_t::SummedADCStdDev); }
    
    /// @}
    
    
    /// Returns the number of hits in the cluster
    virtual size_t NHits() override
      { return ReturnValue(cpNHits, &Base_t::NHits); }
    
    
    /**
     * @brief Fraction of wires in the cluster with more than one hit
     * @return fraction of wires with more than one hit, or 0 if no wires
     *
     * Returns a quantity defined as NMultiHitWires / NWires,
     * where NWires is the number of wires hosting at least one hit of this
     * cluster, and NMultiHitWires is the number of wires which have more
     * than just one hit.
     */
    virtual float MultipleHitDensity() override
      { return ReturnValue(cpMultipleHitDensity, &Base_t::MultipleHitDensity); }
    
    /**
     * @brief Computes the width of the cluster
     * @return width of the cluster
     * 
     * @todo provide a description of the algorithm by words
     */
    virtual float Width() override
      { return ReturnValue(cpWidth, &Base_t::Width); }
    
    /// @}
    
    
      protected:
    
    using ValueFunction_t = float (Base_t::*) ();
    using MeasureFunction_t = Measure_t (Base_t::*) ();
    
    std::vector<details::MultiValue> values; ///< the overridden values
    std::bitset<NParameters> overridden_set; ///< bits for overriding
    
    template <typename Func>
    auto ReturnValue(ParameterType_t param, Func func)
      -> decltype((this->*func)())
      {
        if (isOverridden(param)) {
          // convert here to the return type of the function
          // (even if we are not using that function, it still defines the type)
          return values[(size_t) param];
        }
        else
          return (this->*func)();
      } // ReturnValue()
    
    
  }; // class OverriddenClusterParamsAlg
  
} // namespace cluster



//==============================================================================
//===  Template implementation
//==============================================================================

namespace cluster {
  
  namespace details {
/*
    /// Class holding a value of one among some selected types...
    class MultiValue {
        public:
      // I can't believe I am writing this shit...
      union {
        Measure_t measure;
        float     float_value;
        size_t    size_t_value;
      };
      
      // set all the default stuff
      MultiValue(MultiValue const&) = default;
      MultiValue(MultiValue&&) = default;
      MultiValue& operator= (MultiValue const&) = default;
      MultiValue& operator= (MultiValue &&) = default;
      ~MultiValue() = default;
      
      /// Sets the value from a value of type T; undefined by default
      template <typename T>
      MultiValue& operator= (T);
      
      /// Converts the value to type T; undefined by default
      template <typename T>
      operator T() const;
    };
*/
    
    // specialization: size_t
    template <>
    MultiValue& MultiValue::operator= (size_t value)
      { size_t_value = value; return *this; }
    
    template <>
    MultiValue::operator size_t () const { return size_t_value; }
    
    // specialization: float
    template <>
    MultiValue& MultiValue::operator= (float value)
      { float_value = value; return *this; }
    
    template <>
    MultiValue::operator float () const { return float_value; }
    
    // specialization: Measure_t
    template <>
    MultiValue& MultiValue::operator= (Measure_t value)
      { measure_value = value; return *this; }
    
    template <>
    MultiValue::operator MultiValue::Measure_t () const
      { return measure_value; }
    
  } // namespace details
  
} // namespace cluster

#endif // OVERRIDDENCLUSTERPARAMSALG_H



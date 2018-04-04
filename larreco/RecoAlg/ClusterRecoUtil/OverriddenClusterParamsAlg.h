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
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"



/// Cluster reconstruction namespace
namespace cluster {
  
  namespace details {
    /**
     * @brief Class holding a value of one among some selected types...
     * 
     * The result of default construction is not defined.
     * That's basically throwing away one of the pillars of C++.
     * 
     * This horrible class is supposed to keep a value that you give to it,
     * and give it back to you if you ask nicely.
     * Of course, if you ask something you did not give it first, it will become
     * naughty. In other words, it's caller's responsibility not to ask float
     * when it assigned size_t.
     */
    class MultiValue {
        public:
      using Measure_t = cluster::details::Measure_t<float>;
      
      // I can't believe I am writing this shit
      union {
        Measure_t measure_value;
        float     float_value;
        size_t    size_t_value;
      };
      
      /// Default constructor; it's here only to allow for vectors to be resized
      /// and its effect is undefined. This class is not to be considered valid
      /// until it's assigned a value with the operator= .
      MultiValue() {}
      
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
  class OverriddenClusterParamsAlg: public ClusterParamsAlgBase {
    static_assert(std::is_base_of<ClusterParamsAlgBase, AlgoBase>::value,
      "OverriddenClusterParamsAlg template parameter must derive"
      " from ClusterParamsAlgBase"
      );
      public:
    using Algo_t = AlgoBase;
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
      algo(std::forward<Args>(args)...),
      values(NParameters)
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
    /// @name Standard ClusterParamsAlgBase interface
    /// 
    /// The following methods replicate the ones of the templated Algo_t class.
    /// Except, of course, when they are overridden.
    
    /**
     * @brief Restores the class to post-configuration, pre-initialization state
     * @see Algo_t::Clear()
     */
    virtual void Clear() override { algo.Clear(); }
    
    
    /**
     * @brief Sets the list of input hits
     * @param hits list of pointers to hits
     * @throw undefined in case of error, this method can throw (anything)
     * @see Algo_t::SetHits().
     */
    virtual void SetHits(std::vector<recob::Hit const*> const& hits) override
      { algo.SetHits(hits); }
    
    /**
     * @brief Sets the list of input hits
     * @param hits list of hits (hits will not be modified)
     * @throw undefined in case of error, this method can throw (anything)
     * @see Algo_t::SetHits().
     */
    virtual void SetHits(std::vector<recob::Hit> const& hits) override
      { algo.SetHits(hits); }
    
    
    /// Set the verbosity level; @see Algo_t::SetVerbose().
    virtual void SetVerbose(int level = 1) override
      { ClusterParamsAlgBase::SetVerbose(level); algo.SetVerbose(level); }
    
    
    /// @{
    /// @name Algorithm results
    
    //@{
    /**
     * @brief Computes the charge on the first and last wire of the track
     * @return the charge in ADC counts, with uncertainty
     * @see Algo_t::StartCharge(), Algo_t::EndCharge()
     */
    virtual Measure_t StartCharge() override
      { return ReturnValue(cpStartCharge, &Algo_t::StartCharge); }
    virtual Measure_t EndCharge() override
      { return ReturnValue(cpEndCharge, &Algo_t::EndCharge); }
    //@}
    
    //@{
    /**
     * @brief Computes the angle of the cluster
     * @return angle of the cluster in the wire x time space, in radians
     * @see Algo_t::StartAngle(), Algo_t::EndAngle()
     * 
     * The angle is in the @f$ [ -\pi, \pi ] @f$ range, with 0 corresponding to
     * a cluster parallel to the wire plane and @f$ \pi @f$ to a cluster
     * orthogonal to the wire plane, going farther from it.
     */
    virtual Measure_t StartAngle() override
      { return ReturnValue(cpStartAngle, &Algo_t::StartAngle); }
    
    virtual Measure_t EndAngle() override
      { return ReturnValue(cpEndAngle, &Algo_t::EndAngle); }
    //@}
    
    //@{
    /**
     * @brief Computes the opening angle at the start or end of the cluster
     * @return angle at the start of the cluster, in radians
     * @see Algo_t::StartOpeningAngle(), Algo_t::EndOpeningAngle()
     */
    virtual Measure_t StartOpeningAngle() override
      { return ReturnValue(cpStartOpeningAngle, &Algo_t::StartOpeningAngle); }
    virtual Measure_t EndOpeningAngle() override
      { return ReturnValue(cpEndOpeningAngle, &Algo_t::EndOpeningAngle); }
    //@}
    
    
    /// @name Cluster charge
    /// @{
    /**
     * @brief Computes the total charge of the cluster from Hit::Integral()
     * @return total charge of the cluster, in ADC count units
     * @see IntegralStdDev(), SummedADC()
     * @see Algo_t::Integral()
     */
    virtual Measure_t Integral() override
      { return ReturnValue(cpIntegral, &Algo_t::Integral); }
    
    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see Integral()
     * @see Algo_t::IntegralStdDev()
     */
    virtual Measure_t IntegralStdDev() override
      { return ReturnValue(cpIntegralStdDev, &Algo_t::IntegralStdDev); }
    
    /**
     * @brief Computes the total charge of the cluster from Hit::SummedADC()
     * @return total charge of the cluster, in ADC count units
     * @see SummedADCStdDev(), Integral()
     * @see Algo_t::SummedADC()
     */
    virtual Measure_t SummedADC() override
      { return ReturnValue(cpSummedADC, &Algo_t::SummedADC); }
    
    /**
     * @brief Computes the standard deviation on the charge of the cluster hits
     * @return the standard deviation of charge of hits, in ADC count units
     * @see SummedADC()
     * @see Algo_t::SummedADCStdDev()
     */
    virtual Measure_t SummedADCStdDev() override
      { return ReturnValue(cpSummedADCStdDev, &Algo_t::SummedADCStdDev); }
    
    /// @}
    
    
    /// Returns the number of hits in the cluster
    virtual size_t NHits() override
      { return ReturnValue(cpNHits, &Algo_t::NHits); }
    
    
    /**
     * @brief Fraction of wires in the cluster with more than one hit
     * @return fraction of wires with more than one hit, or 0 if no wires
     * @see Algo_t::MultipleHitDensity()
     */
    virtual float MultipleHitDensity() override
      { return ReturnValue(cpMultipleHitDensity, &Algo_t::MultipleHitDensity); }
    
    /**
     * @brief Computes the width of the cluster
     * @return width of the cluster
     * @see Algo_t::Width()
     */
    virtual float Width() override
      { return ReturnValue(cpWidth, &Algo_t::Width); }
    
    /// @}
    
    /// @}
    
      protected:
    
    using ValueFunction_t = float (Algo_t::*) ();
    using MeasureFunction_t = Measure_t (Algo_t::*) ();
    
    Algo_t algo; ///< an instance of the wrapped algorithm class
    
    std::vector<details::MultiValue> values; ///< the overridden values
    std::bitset<NParameters> overridden_set; ///< bits for overriding
    
    template <typename Func>
    auto ReturnValue(ParameterType_t param, Func func)
      -> decltype((algo.*func)())
      {
        if (isOverridden(param)) {
          // convert here to the return type of the function
          // (even if we are not using that function, it still defines the type)
          return values[(size_t) param];
        }
        else
          return (algo.*func)();
      } // ReturnValue()
    
  }; // class OverriddenClusterParamsAlg
  
} // namespace cluster



//==============================================================================
//===  Template implementation
//==============================================================================

namespace cluster {
  
  namespace details {
    
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



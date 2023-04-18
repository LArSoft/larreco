/**
 * @file   GausFitCache.h
 * @brief  Provide caches for TF1 functions to be used with ROOT fitters
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 25th, 2015
 *
 *
 * This might be moved to another namespace than hit, and to a more generic
 * function cacher.
 */

#ifndef GAUSFITCACHE_H
#define GAUSFITCACHE_H 1

// C/C++ standard libraries
#include <string>
#include <vector>

// ROOT libraries
#include "RtypesCore.h" // Double_t
#include "TF1.h"

namespace hit {

  /** **************************************************************************
   * @brief A set of TF1 linear sum of base functions (Gaussians)
   *
   * This object holds and owns a number of functions that are linear sum of
   * "base" functions.
   * The original use is for Gaussian functions, and this class implements
   * Gaussians as base functions. This class is used as a base class for other
   * function caches that do not necessarily use Gaussian base functions.
   *
   * The idea is that if you need a temporary function to perform a fit and
   * extract the parameters, you can reuse the same function over and over
   * rather than creating a new one each time.
   *
   * The expected use for an algorithm is to have one of these caches as a data
   * member and when the algorithm code needs a function it obtains it with:
   *
   *     TF1* pFunc = fitcache.Get(nGaussians);
   *
   * Then `pFunc` is used for fitting but *never* destroyed.
   */
  class GausFitCache {
  public:
    /// Constructor; optionally set the name of the repository
    GausFitCache(std::string new_name = "GausFitCache") : name(new_name) {}

    /// Destructor
    virtual ~GausFitCache();

    /// Return the name of this cache
    std::string GetName() const { return name; }

    /**
     * @brief Returns a function sum of nFunc base functions
     * @param nGaus the number of base functions in the function
     * @return a pointer to a TF1 with the required characteristics
     *
     * The returned function must not be deleted at any time!
     *
     * This implementation returns a function sum of nFunc Gaussians.
     * The parameters are sorted by Gaussian: first normalization (not the area,
     * but the coefficient), first mean, first sigma, second normalization etc.
     */
    virtual TF1* Get(size_t nFunc);

    /// Returns a new function sum of nFunc base functions
    /// (caller needs to set limits and parameters)
    virtual TF1* GetClone(size_t nGaus);

    /// Returns a name for the function with nFunc base functions
    virtual std::string FunctionName(size_t nFunc) const;

  protected:
    std::string name; ///< name of the cache

    ///< Gaussian sum functions; n-th element is sum of n base functions
    std::vector<TF1*> funcs;

    /// Creates a new sum function
    virtual TF1* CreateFunction(size_t nFunc) const;

  }; // class GausFitCache

  namespace details {

    template <typename T>
    inline T sqr(T v)
    {
      return v * v;
    }

    /// Struct with member type corresponding to the NArg-th type from Args
    template <unsigned int NArg, typename FirstArg, typename... Args>
    struct TemplateArgumentHelper {
      using type = typename TemplateArgumentHelper<NArg - 1, Args...>::type;
    }; // struct TemplateArgumentHelper

    /// Struct with member type corresponding to the NArg-th type from Args
    template <unsigned int NArg, typename... Args>
    struct TemplateArgument {
      using type = typename TemplateArgumentHelper<NArg, Args...>::type;
    }; // struct TemplateArgument

    /**
     * @brief A sum of NFunc base functions Func
     * @tparam NFunc the number of base functions in the sum
     * @tparam Func the base function in the sum
     * @tparam NFuncParams the number of parameters required by Func
     *
     * This class provides in its `eval` member a compiled function suitable
     * to be wrapped into a ROOT's TF1 object.
     * The function is the sum of NFunc base functions.
     * Each base function is expected to use NFuncParams parameters. The
     * first function will use the first set of NFuncParams parameters,
     * the second one the next set of NFuncParams parameters, and so on.
     * If NFunc is 0, the value 0 is always returned.
     */
    template <unsigned int NFunc,
              Double_t Func(Double_t const*, Double_t const*),
              unsigned int NFuncParams>
    struct FuncSum {
      static Double_t eval(Double_t const*, Double_t const*);
      static constexpr unsigned int NParams = NFunc * NFuncParams;
    }; // struct FuncSum

    // partial specialization (declaration)
    template <Double_t Func(Double_t const*, Double_t const*), unsigned int NFuncParams>
    struct FuncSum<0U, Func, NFuncParams> {
      static Double_t eval(Double_t const*, Double_t const*);
      static constexpr unsigned int NParams = 0;
    }; // struct FuncSum<0>

    class CompiledGausFitCacheBaseStruct : public GausFitCache {
      /*
       * Since, as most of C++ metaprogramming code, this one is quite messy,
       * some explanations follow.
       *
       * The idea is to have a function sum of N base functions, and both N and
       * which base function can be specified as template parameters.
       *
       * The implementation of the sum of N base functions requires a
       * compile-time "for loop" that is implemented by recursion: a function
       * with a template parameter N=i calls the same with N=i-1. The loop is
       * interrupted by a template specialization for N=0 that does not call
       * N=-1. Unfortunately there is another template parameter that is carried
       * around, that is the base function. This means that the specialization
       * for N=0 is a partial specialization, since the function still needs
       * to be a template parameter. Partial specializations are only allowed
       * for classes (not for functions). Therefore we have here a class
       * with two template parameters and a static member eval() where the
       * total function is defined. Partial specializations of the full class
       * will implement the necessary N=0 exception to the recursion rule.
       *
       * Now, we want to put these in a vector of functions.
       * This is another loop, and since the function types are decided at
       * compile time, the loop must be done at compile time too.
       * Here it goes another class, moving around a vector.
       * Furthermore, now the actual type of function is a different one for
       * each element in the array: in this compile-time for-loop, the type
       * depends on the running parameter, and potentially on others.
       *
       *
       *
       *
       */
    public:
      using GausFitCache::GausFitCache;

      /// Throws an exception (ROOT does not support cloning compiled functions)
      virtual TF1* GetClone(size_t nGaus);

      /// Returns the maximum number of Gaussians in a function that we support
      virtual unsigned int MaxGaussians() const { return funcs.size() - 1; }

      /**
       * @brief Single Gaussian function
       * @param x vector; [0] variable value
       * @param params vector: [0] coefficient of the exponent [1] mean [2] sigma
       * @return the Gaussian evaluated at *x
       *
       * This function is equivalent to ROOT's "gaus".
       */
      static Double_t gaus(Double_t const* x, Double_t const* params);

      template <unsigned int CutOff>
      static Double_t gaus_trunc(Double_t const* x, Double_t const* params);

      template <unsigned int NGaus>
      static Double_t ngaus(Double_t const* x, Double_t const* params)
      {
        return gaus(x, params) + ngaus<NGaus - 1>(x, params + 3);
      }

      /// Sum of NGaus Gaussian functions truncated at CutOff sigmas
      template <unsigned int NGaus, unsigned int CutOff>
      static Double_t ngaus_trunc(Double_t const* x, Double_t const* params)
      {
        return FuncSum<NGaus, gaus_trunc<CutOff>, 3U>::eval(x, params);
      }

      /// Class around sum of NGaus Gaussian functions truncated at CutOff sigmas
      template <unsigned int NGaus, unsigned int CutOff>
      using NGaussTruncClass = FuncSum<NGaus, gaus_trunc<CutOff>, 3U>;

    protected:
      /**
       * @brief A helper class initializing the function vector
       * @tparam NFunc the maximum number of base functions in a function
       * @tparam Func the function encapsulating object
       * @tparam Others other template parameters of the function object
       *
       * This helper class provides a static function fill() that fills the
       * function vector of the specified CompiledGausFitCacheBaseStruct
       * with a sequence on NFunc+1 functions, of types `Func<0, Others...>`,
       * `Func<1, Others...>`, ... up to `Func<NFunc, Others...>` included.
       * Each function is actually wrapped by a ROOT's TF1.
       * Func type is a class, and the actual function must be found in the
       * member `Func::eval`. Also the Func class is required to have a
       * `NParams` static constexpr member whose value is the number of
       * `parameters that the `eval()` function takes.
       *
       * The template FuncSum is implementing all these requirements.
       */
      template <unsigned int NFunc, template <unsigned int> class Func>
      struct InitializeFuncSumVector {
        static void fill(CompiledGausFitCacheBaseStruct& cache);
      }; // struct InitializeFuncSumVector

      // partial specialization: 0 of any function
      template <template <unsigned int> class Func>
      struct InitializeFuncSumVector<0U, Func> {
        static void fill(CompiledGausFitCacheBaseStruct& cache);
      }; // struct InitializeFuncSumVector<0, Func>

      /// Returns a vector initialized with multigaussians
      template <unsigned int NGaus>
      void InitializeCompiledGausFitVector();

      /// Adds one function
      template <unsigned int NGaus>
      void AppendFunction();

      /// Throws an error asserting compiled functions can't be cretead run-time
      void CannotCreateFunction [[noreturn]] (size_t nGaus) const;

    }; // class CompiledGausFitCacheBaseStruct

  } // namespace details

  /** **************************************************************************
   * @brief A set of TF1 linear sum of Gaussians
   * @tparam MaxGaus the maximum number of Gaussians in the stored functions
   *
   * This class stores a predefined number MaxGaus of TF1
   * from pre-compiled functions.
   */
  template <unsigned int MaxGaus = 10>
  class CompiledGausFitCache : public details::CompiledGausFitCacheBaseStruct {
  public:
    /// Constructor: initializes all the functions
    CompiledGausFitCache(std::string new_name = "CompiledGausFitCache")
      : details::CompiledGausFitCacheBaseStruct(new_name)
    {
      InitializeCompiledGausFitVector<MaxGaus>();
    }

    virtual unsigned int MaxGaussians() const { return StoredMaxGaussians(); }

    /// Returns the maximum number of Gaussians in a function that we support
    constexpr unsigned int StoredMaxGaussians() const { return MaxGaus; }

  protected:
    /// Throws an error, since this class can't create functions run-time
    virtual TF1* CreateFunction [[noreturn]] (size_t nGaus) const { CannotCreateFunction(nGaus); }

  }; // class CompiledGausFitCache

  /** **************************************************************************
   * @brief A set of TF1 linear sum of truncated Gaussians
   * @tparam MaxGaus the maximum number of Gaussians in the stored functions
   * @tparam CutOff number of sigmas beyond which the function evaluates as 0
   *
   * This class stores a predefined number MaxGaus of TF1 from pre-compiled
   * functions. Each function is the linear sum of truncated Gaussians, at most
   * MaxGaus of them,
   */
  template <unsigned int MaxGaus = 10, unsigned int CutOff = 5>
  class CompiledTruncatedGausFitCache : public details::CompiledGausFitCacheBaseStruct {
    template <unsigned int NGaus>
    using CutOffNGaussianClass = NGaussTruncClass<NGaus, CutOff>;

  public:
    /// Constructor: initializes all the functions
    CompiledTruncatedGausFitCache(std::string new_name = "CompiledTruncatedGausFitCache")
      : details::CompiledGausFitCacheBaseStruct(new_name)
    {
      InitializeFuncSumVector<MaxGaus, CutOffNGaussianClass>::fill(*this);
    }

    virtual unsigned int MaxGaussians() const { return StoredMaxGaussians(); }

    /// Returns the maximum number of Gaussians in a function that we support
    constexpr unsigned int StoredMaxGaussians() const { return MaxGaus; }

  protected:
    /// Throws an error, since this class can't create functions run-time
    virtual TF1* CreateFunction [[noreturn]] (size_t nGaus) const { CannotCreateFunction(nGaus); }

  }; // class CompiledTruncatedGausFitCache

  //
  // template implementation
  //

  namespace details {

    // --- FuncSum -------------------------------------------------------------
    template <unsigned int NFunc,
              Double_t Func(Double_t const*, Double_t const*),
              unsigned int NFuncParams>
    constexpr unsigned int FuncSum<NFunc, Func, NFuncParams>::NParams;

    template <unsigned int NFunc,
              Double_t Func(Double_t const*, Double_t const*),
              unsigned int NFuncParams>
    Double_t FuncSum<NFunc, Func, NFuncParams>::eval(Double_t const* x, Double_t const* params)
    {
      return Func(x, params + NFuncParams * (NFunc - 1)) // use the last parameters
             + FuncSum<NFunc - 1, Func, NFuncParams>::eval(x, params);
    } // CompiledGausFitCacheBaseStruct::FuncSum<NFunc, Func>::eval()

    // partial specialization: 0 of any function
    template <Double_t Func(Double_t const*, Double_t const*), unsigned int NFuncParams>
    Double_t FuncSum<0U, Func, NFuncParams>::eval(Double_t const*, Double_t const*)
    {
      return 0.;
    }

    // --- CompiledGausFitCacheBaseStruct --------------------------------------

    template <unsigned int NGaus>
    void CompiledGausFitCacheBaseStruct::InitializeCompiledGausFitVector()
    {
      if (NGaus > 0) InitializeCompiledGausFitVector<NGaus - 1>();
      AppendFunction<NGaus>();
    } // CompiledGausFitCacheBaseStruct::InitializeCompiledGausFitVector()

    template <>
    inline void CompiledGausFitCacheBaseStruct::InitializeCompiledGausFitVector<0>()
    {
      AppendFunction<0>();
    }

    template <unsigned int NGaus>
    void CompiledGausFitCacheBaseStruct::AppendFunction()
    {
      // create a function in the ficticious range [ 0, 1 ]:
      funcs.push_back(new TF1(FunctionName(NGaus).c_str(), &ngaus<NGaus>, 0., 1., 3 * NGaus));
    } // CompiledGausFitCacheBaseStruct::AppendFunction()

    template <unsigned int CutOff>
    Double_t CompiledGausFitCacheBaseStruct::gaus_trunc(Double_t const* x, Double_t const* params)
    {
      const Double_t z = (x[0] - params[1]) / params[2];
      return ((z > -((Double_t)CutOff)) && (z < (Double_t)CutOff)) ?
               params[0] * std::exp(-0.5 * sqr(z)) :
               0.;
    } // CompiledGausFitCacheBaseStruct::gaus_trunc()

    template <>
    inline Double_t CompiledGausFitCacheBaseStruct::ngaus<0>(Double_t const* x,
                                                             Double_t const* params)
    {
      return 0.;
    }

    // --- CompiledGausFitCacheBaseStruct::InitializeFuncSumVector -------------

    template <unsigned int NFunc, template <unsigned int> class Func>
    void CompiledGausFitCacheBaseStruct::InitializeFuncSumVector<NFunc, Func>::fill(
      CompiledGausFitCacheBaseStruct& cache)
    {
      // first fill the lower functions
      InitializeFuncSumVector<NFunc - 1, Func>::fill(cache);
      // then add one
      cache.funcs.push_back(new TF1(
        cache.FunctionName(NFunc).c_str(), Func<NFunc>::eval, 0., 1., Func<NFunc>::NParams));
    } // InitializeFuncSumVector<NFunc, Func>::fill()

    template <template <unsigned int> class Func>
    void CompiledGausFitCacheBaseStruct::InitializeFuncSumVector<0U, Func>::fill(
      CompiledGausFitCacheBaseStruct& cache)
    {
      cache.funcs.push_back(new TF1(cache.FunctionName(0).c_str(), Func<0U>::eval, 0., 1., 0));
    } // InitializeFuncSumVector<0, Func>::fill()

    // -------------------------------------------------------------------------

  } // namespace details

  // ---------------------------------------------------------------------------

} // namespace hit

#endif // GAUSFITCACHE_H

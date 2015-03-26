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
#include "Rtypes.h" // Double_t
#include "TF1.h" // Double_t


namespace hit {
  
  /** **************************************************************************
   * @brief A set of TF1 linear sum of Gaussians
   * 
   * This object holds and owns a number of functions that are linear sum of
   * Gaussian functions.
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
    GausFitCache(std::string new_name = "GausFitCache"): name(new_name) {}
    
    /// Destructor
    virtual ~GausFitCache();
    
    /**
     * @brief Returns a function sum of nGaus Gaussians
     * @param nGaus the number of Gaussians in the function
     * @return a pointer to a TF1 with the required characteristics
     * 
     * The returned function must not be deleted at any time!
     * 
     * The parameters are sorted by Gaussian: first normalization (not the area,
     * but the coefficient), first mean, first sigma, second normalization etc.
     */
    virtual TF1* Get(size_t nGaus);
    
    /// Returns a new function sum of nGaus Gaussians
    /// (caller needs to set limits and parameters)
    TF1* GetClone(size_t nGaus);
    
    /// Returns a name for the function with nGaus Gaussians
    virtual std::string FunctionName(size_t nGaus) const;
    
      protected:
    
    std::string name; ///< name of the cache
    
    ///< Gaussian sum functions; n-th element is sum of n Gaussians
    std::vector<TF1*> gaussians;
    
    /// Creates a new Gaussian sum function
    virtual TF1* CreateFunction(size_t nGaus) const;
    
  }; // class GausFitCache
  
  
  
  namespace details {
    class CompiledGausFitCacheBaseStruct: public GausFitCache {
        public:
      using GausFitCache::GausFitCache;
      
        protected:
      /**
       * @brief Single Gaussian function
       * @param x vector; [0] variable value
       * @param params vector: [0] coefficient of the exponent [1] mean [2] sigma
       * @return the Gaussian evaluated at *x
       * 
       * This function is equivalent to ROOT's "gaus".
       */
      static Double_t gaus(Double_t const* x, Double_t const* params);
      
      template <unsigned int NGaus>
      static Double_t ngaus(Double_t const* x, Double_t const* params)
        { return gaus(x, params) + ngaus<NGaus-1>(x, params + 3); }
      
      /// Add the multi-Gaussian with NGaus Gaussians to the cache
      template <unsigned int NGaus>
      static void InitializeFunction();
      
      /// Returns a vector initialized with multigaussians
      template <unsigned int NGaus>
      void InitializeCompiledGausFitVector();
      
      /// Adds one function
      template <unsigned int NGaus>
      void AppendFunction();
      
    }; // class CompiledGausFitCacheBaseStruct
    
  } // namespace details
  
  
  /** **************************************************************************
   * @brief A set of TF1 linear sum of Gaussians
   * @tparam MAXGAUS the maximum number of Gaussians in the stored functions
   * 
   * This class stores a predefined number MAXGAUS of TF1
   * from pre-compiled functions.
   */
  template <unsigned int MAXGAUS = 10>
  class CompiledGausFitCache: public details::CompiledGausFitCacheBaseStruct
  {
      public:
    
    /// Constructor: initializes all the functions
    CompiledGausFitCache(std::string new_name = "CompiledGausFitCache"):
      GausFitCache(new_name)
      { InitializeCompiledGausFitVector<MAXGAUS>(); }
    
    
    /// Returns the maximum number of Gaussians in a function that we support
    unsigned int MaxGaus() const { return gaussians.size() - 1; }
    
      protected:
    
    /// Throws an error, since this class can't create functions run-time
    virtual TF1* CreateFunction(size_t nGaus) const;
    
  }; // class CompiledGausFitCache
  
  
  
  //
  // template implementation
  //
  namespace details {
    template <unsigned int NGaus>
    void CompiledGausFitCacheBaseStruct::InitializeCompiledGausFitVector() {
      if (NGaus > 0) InitializeCompiledGausFitVector<NGaus-1>();
      AppendFunction<NGaus>();
    } // CompiledGausFitCacheBaseStruct::InitializeCompiledGausFitVector()
    
    template <unsigned int NGaus>
    void CompiledGausFitCacheBaseStruct::AppendFunction()
      { gaussians.push_back(new TF1(FunctionName(NGaus).c_str(), &ngaus<NGaus>)); }
    
    
    template <>
    inline Double_t details::CompiledGausFitCacheBaseStruct::ngaus<0>
      (Double_t const* x, Double_t const* params)
      { return 0.; }
    
  } // namespace details
  
} // namespace hit


#endif // GAUSFITCACHE_H
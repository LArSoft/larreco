/**
 * @file   GausFitCache.cxx
 * @brief  Provide caches for TF1 functions to be used with ROOT fitters
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 25th, 2015
 * @see    GausFitCache.h
 */

// Library header
#include "GausFitCache.h"

// C/C++ standard libraries
#include <cmath> // std::sqrt(), std::abs()
#include <sstream>
#include <memory> // std::default_delete<>
#include <algorithm> // std::sort(), std::for_each()

// ROOT libraries
#include "TF1.h"

// Framework libraries
#include "art/Utilities/Exception.h"


namespace {
  
  template <typename T>
  inline T sqr(T v) { return v*v; }
  
} // local namespace


namespace hit {
  //----------------------------------------------------------------------------
  //--- GausFitCache
  //---
  GausFitCache::~GausFitCache() {
    std::for_each
      (gaussians.begin(), gaussians.end(), std::default_delete<TF1>());
  } // GausFitCache::~GausFitCache()
  
  
  //----------------------------------------------------------------------------
  TF1* GausFitCache::Get(size_t nGaus) {
    
    // expand the list if needed
    if (nGaus >= gaussians.size()) gaussians.resize(nGaus + 1, nullptr);
    
    // get the pointer we need...
    TF1*& pFunc = gaussians[nGaus];
    
    if (!pFunc) pFunc = CreateFunction(nGaus);
    
    return pFunc;
  } // GausFitCache::Get()
  
  //----------------------------------------------------------------------------
  TF1* GausFitCache::GetClone(size_t nGaus)
    { return static_cast<TF1*>(Get(nGaus)->Clone()); }
  
  
  //----------------------------------------------------------------------------
  TF1* GausFitCache::CreateFunction(size_t nGaus) const {
    
    std::string func_name = FunctionName(nGaus);
    
    // no gaussian, no nothing
    if (nGaus == 0) return new TF1(func_name.c_str(), "0");
    
    std::ostringstream sstr;
    sstr << "gaus(0)";
    for (size_t iGaus = 1; iGaus < nGaus; ++iGaus)
      sstr << " + gaus(" << (iGaus*3) << ")";
    
    // create and return the function
    return new TF1(func_name.c_str(), sstr.str().c_str());
  } // GausFitCache::CreateFunction()
  
  
  //----------------------------------------------------------------------------
  std::string GausFitCache::FunctionName(size_t nGaus) const {
    std::ostringstream sstr;
    sstr << name << "_" << nGaus;
    return sstr.str();
  } // GausFitCache::FunctionName()
  
  
  //----------------------------------------------------------------------------
  //--- CompiledGausFitCache
  //---
  template <unsigned int MAXGAUS>
  TF1* CompiledGausFitCache<MAXGAUS>::CreateFunction(size_t nGaus) const {
    throw art::Exception("CompiledGausFitCache")
      << "a " << nGaus
      << "-gaussian function was requested, but we have only up to "
      << (gaussians.size() - 1) << "-gaussian functions available";
  } // CompiledGausFitCache<>::CreateFunction()
  
  
  //----------------------------------------------------------------------------
  Double_t details::CompiledGausFitCacheBaseStruct::gaus
    (Double_t const* x, Double_t const* params)
  {
    return params[0] * std::exp(-0.5*sqr((x[0] - params[1])/params[2]));
  } // details::CompiledGausFitCacheBaseStruct::gaus()
  
  
  //----------------------------------------------------------------------------
  
} // namespace hit

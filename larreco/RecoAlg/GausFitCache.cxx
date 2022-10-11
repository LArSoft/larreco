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
#include <algorithm> // std::sort(), std::for_each()
#include <cmath>     // std::sqrt(), std::abs()
#include <memory>    // std::default_delete<>
#include <sstream>

// ROOT libraries
#include "TF1.h"

// Framework libraries
#include "canvas/Utilities/Exception.h" // for Exception, LogicError

namespace {

  template <typename T>
  inline T sqr(T v)
  {
    return v * v;
  }

} // local namespace

namespace hit {
  //----------------------------------------------------------------------------
  //--- GausFitCache
  //---
  GausFitCache::~GausFitCache()
  {
    std::for_each(funcs.begin(), funcs.end(), std::default_delete<TF1>());
  } // GausFitCache::~GausFitCache()

  //----------------------------------------------------------------------------
  TF1* GausFitCache::Get(size_t nFunc)
  {

    // expand the list if needed
    if (nFunc >= funcs.size()) funcs.resize(nFunc + 1, nullptr);

    // get the pointer we need...
    TF1*& pFunc = funcs[nFunc];

    if (!pFunc) pFunc = CreateFunction(nFunc);

    return pFunc;
  } // GausFitCache::Get()

  //----------------------------------------------------------------------------
  TF1* GausFitCache::GetClone(size_t nFunc) { return static_cast<TF1*>(Get(nFunc)->Clone()); }

  //----------------------------------------------------------------------------
  TF1* GausFitCache::CreateFunction(size_t nFunc) const
  {

    std::string func_name = FunctionName(nFunc);

    // no gaussian, no nothing
    if (nFunc == 0) return new TF1(func_name.c_str(), "0");

    std::ostringstream sstr;
    sstr << "gaus(0)";
    for (size_t iGaus = 1; iGaus < nFunc; ++iGaus)
      sstr << " + gaus(" << (iGaus * 3) << ")";

    // create and return the function
    return new TF1(func_name.c_str(), sstr.str().c_str());
  } // GausFitCache::CreateFunction()

  //----------------------------------------------------------------------------
  std::string GausFitCache::FunctionName(size_t nFunc) const
  {
    std::ostringstream sstr;
    sstr << name << "_" << nFunc;
    return sstr.str();
  } // GausFitCache::FunctionName()

  //----------------------------------------------------------------------------
  //--- CompiledGausFitCacheBaseStruct
  //---
  Double_t details::CompiledGausFitCacheBaseStruct::gaus(Double_t const* x, Double_t const* params)
  {
    return params[0] * std::exp(-0.5 * sqr((x[0] - params[1]) / params[2]));
  } // details::CompiledGausFitCacheBaseStruct::gaus()

  TF1* details::CompiledGausFitCacheBaseStruct::GetClone(size_t nFunc)
  {
    throw art::Exception(art::errors::LogicError)
      << "CompiledGausFitCacheBaseStruct: compiled functions can't be cloned";
  } // CompiledGausFitCacheBaseStruct::GetClone()

  void details::CompiledGausFitCacheBaseStruct::CannotCreateFunction(size_t nFunc) const
  {
    throw art::Exception(art::errors::LogicError)
      << name << " function cache can't create functions at run-time; " << nFunc
      << "-addend function was requested, but we have only functions with up to " << MaxGaussians()
      << " addends available\n";
  } // CompiledGausFitCache<>::CannotCreateFunction()

  //----------------------------------------------------------------------------

} // namespace hit

/**
 * @file   GausFitCache_test.cc
 * @brief  Test for classes in GausFitCache,h
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 27th, 2015
 * @see    GausFitCache,h
 */

// define the following symbol to produce plots of what's going on
// #define GAUSFITCACHE_TEST_DEBUG 1

// C/C++ standard libraries
#include <cassert>
#include <cmath>
#include <array>
#include <vector>
#include <limits> // std::numeric_limits<>
// #include <iostream>

// boost test libraries
#define BOOST_TEST_MODULE ( HitAnaAlg_test )
#include "cetlib/quiet_unit_test.hpp"
#include <cetlib/quiet_unit_test.hpp> // BOOST_CHECK_CLOSE

// ROOT libraries
#include "TH1D.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#ifdef GAUSFITCACHE_TEST_DEBUG
#  include "TCanvas.h"
#endif // GAUSFITCACHE_TEST_DEBUG

// LArSoft libraries
#include "larreco/RecoAlg/GausFitCache.h"


template <typename T>
inline T sqr(T v) { return v*v; }


double gaus(double x, double mean, double sigma, double amplitude) {
  register double z = (x - mean) / sigma;
  return amplitude * std::exp(-0.5*sqr(z));
} // gaus()


//******************************************************************************
BOOST_AUTO_TEST_SUITE( GausFitCacheSuite )


//******************************************************************************
// test that the test function behaves like a gaussian
BOOST_AUTO_TEST_CASE(TestGaussianTest)
{
  // static function to be tested:
  const auto gaus_func = ::gaus;
  
  const double exp1sigma = std::exp(-0.5*sqr(1.));
  const double exp2sigma = std::exp(-0.5*sqr(2.));
  
  // here there should be some more believable test...
  const double amplitude = 16.;
  for (double sigma = 0.5; sigma < 2.25; ++sigma) {
    for (double mean = -5; mean < 5.5; ++mean) {
      
    //  std::cout << "Testing local gaus(), mean=" << mean << " sigma=" << sigma
    //    << " amplitude=" << amplitude << std::endl;
      
      BOOST_CHECK      (gaus_func(-6.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func(-5.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func(-4.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      
      // we use tolerance of 10^-5 (0.001%)
      BOOST_CHECK_CLOSE(gaus_func(-2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func(-1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 0.0 * sigma + mean, mean, sigma, amplitude), amplitude, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      
      BOOST_CHECK      (gaus_func( 4.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func( 5.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func( 6.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      
    } // for mean
  } // for sigma
  
} // BOOST_AUTO_TEST_CASE(TestGaussianTest)


//******************************************************************************
// test that the library functions behaves like a gaussians

struct RootGausFuncWrapper {
  using Func_t = Double_t (*) (Double_t const*, Double_t const*);
  
  RootGausFuncWrapper(Func_t f): func(f) {}
  
  Double_t operator()
    (Double_t const x, Double_t mean, Double_t sigma, Double_t amplitude) const
    {
      // BUG the double brace syntax is required to work around clang bug 21629
      // (https://bugs.llvm.org/show_bug.cgi?id=21629)
      std::array<Double_t, 3> params = {{ amplitude, mean, sigma }};
      return func(&x, params.data());
    } // operator()
  
  Func_t func;
}; // RootGausFuncWrapper()


BOOST_AUTO_TEST_CASE(GaussianTest)
{
  // static function to be tested:
  const auto gaus_func = RootGausFuncWrapper
    (hit::details::CompiledGausFitCacheBaseStruct::gaus);
  
  const double exp1sigma = std::exp(-0.5*sqr(1.));
  const double exp2sigma = std::exp(-0.5*sqr(2.));
  
  const double amplitude = 16.;
  for (double sigma = 0.5; sigma < 2.25; ++sigma) {
    for (double mean = -5; mean < 5.5; ++mean) {
      
    //  std::cout << "Testing gaus, mean=" << mean << " sigma=" << sigma
    //    << " amplitude=" << amplitude << std::endl;
      
      BOOST_CHECK      (gaus_func(-6.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func(-5.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func(-4.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      
      // we use tolerance of 10^-5 (0.001%)
      BOOST_CHECK_CLOSE(gaus_func(-2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func(-1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 0.0 * sigma + mean, mean, sigma, amplitude), amplitude, 0.001);
      BOOST_CHECK_CLOSE(gaus_func(+1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func(+2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      
      BOOST_CHECK      (gaus_func(+4.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func(+5.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK      (gaus_func(+6.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      
    } // for mean
  } // for sigma
  
} // BOOST_AUTO_TEST_CASE(GaussianTest)


BOOST_AUTO_TEST_CASE(GaussianTrunc5Test)
{
  // static function to be tested:
  const auto gaus_func = RootGausFuncWrapper
    (hit::details::CompiledGausFitCacheBaseStruct::gaus_trunc<5>);
  
  const double exp1sigma = std::exp(-0.5*sqr(1.));
  const double exp2sigma = std::exp(-0.5*sqr(2.));
  
  const double amplitude = 16.;
  for (double sigma = 0.5; sigma < 2.25; ++sigma) {
    for (double mean = -5; mean < 5.5; ++mean) {
      
    //  std::cout << "Testing gaus_trunc<5>, mean=" << mean << " sigma=" << sigma
    //    << " amplitude=" << amplitude << std::endl;
      
      BOOST_CHECK_EQUAL(gaus_func(-6.5 * sigma + mean, mean, sigma, amplitude), 0.);
      BOOST_CHECK_EQUAL(gaus_func(-5.5 * sigma + mean, mean, sigma, amplitude), 0.);
      BOOST_CHECK      (gaus_func(-4.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      
      // we use tolerance of 10^-5 (0.001%)
      BOOST_CHECK_CLOSE(gaus_func(-2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func(-1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 0.0 * sigma + mean, mean, sigma, amplitude), amplitude, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      
      BOOST_CHECK      (gaus_func( 4.5 * sigma + mean, mean, sigma, amplitude) > 0.);
      BOOST_CHECK_EQUAL(gaus_func( 5.5 * sigma + mean, mean, sigma, amplitude), 0.);
      BOOST_CHECK_EQUAL(gaus_func( 6.5 * sigma + mean, mean, sigma, amplitude), 0.);
      
    } // for mean
  } // for sigma
  
} // BOOST_AUTO_TEST_CASE(GaussianTrunc5Test)


BOOST_AUTO_TEST_CASE(GaussianTrunc4Test)
{
  // static function to be tested:
  const auto gaus_func = RootGausFuncWrapper
    (hit::details::CompiledGausFitCacheBaseStruct::gaus_trunc<4>);
  
  const double exp1sigma = std::exp(-0.5*sqr(1.));
  const double exp2sigma = std::exp(-0.5*sqr(2.));
  
  const double amplitude = 16.;
  for (double sigma = 0.5; sigma < 2.25; ++sigma) {
    for (double mean = -5; mean < 5.5; ++mean) {
      
    //  std::cout << "Testing gaus_trunc<4>, mean=" << mean << " sigma=" << sigma
    //    << " amplitude=" << amplitude << std::endl;
      
      BOOST_CHECK_EQUAL(gaus_func(-6.5 * sigma + mean, mean, sigma, amplitude), 0.);
      BOOST_CHECK_EQUAL(gaus_func(-5.5 * sigma + mean, mean, sigma, amplitude), 0.);
      BOOST_CHECK_EQUAL(gaus_func(-4.5 * sigma + mean, mean, sigma, amplitude), 0.);
      
      // we use tolerance of 10^-5 (0.001%)
      BOOST_CHECK_CLOSE(gaus_func(-2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func(-1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 0.0 * sigma + mean, mean, sigma, amplitude), amplitude, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 1.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus_func( 2.0 * sigma + mean, mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      
      BOOST_CHECK_EQUAL(gaus_func( 4.5 * sigma + mean, mean, sigma, amplitude), 0.);
      BOOST_CHECK_EQUAL(gaus_func( 5.5 * sigma + mean, mean, sigma, amplitude), 0.);
      BOOST_CHECK_EQUAL(gaus_func( 6.5 * sigma + mean, mean, sigma, amplitude), 0.);
      
    } // for mean
  } // for sigma
  
} // BOOST_AUTO_TEST_CASE(GaussianTrunc4Test)


//******************************************************************************
// test the multi-gaussian functions

/// Expect for each Gaussian ROOT-like parameters: amplitude, mean, sigma
Double_t multi_gaus
  (Double_t x, const unsigned int nGaus, Double_t const* params)
{
  
  Double_t res = 0.;
  for (unsigned int iGaus = 0; iGaus < nGaus; ++iGaus) {
    // our gaus() takes the parameters as: mean, sigma, amplitude
    res += gaus(x, params[1], params[2], params[0]);
    params += 3;
  } // for
  return res;
} // multi_gaus()


// Returns the parameters in pFunc sorted to follow the expected ones in Params
std::vector<Double_t> SortGaussianResults
  (TF1 const* pFunc, Double_t const* Params)
{
  assert(pFunc->GetNpar() % 3 == 0); // Sorting non-Gaussian function!
  
  const size_t nGaus = pFunc->GetNpar() / 3;
  
  std::vector<size_t> BestMatch(nGaus, nGaus); // initialize with invalid number
  for (size_t iGaus = 0; iGaus < nGaus; ++iGaus) {
    double best_match_quality = std::numeric_limits<double>::max();
    for (size_t iFitGaus = 0; iFitGaus < nGaus; ++iFitGaus) {
      // so far, we compare all the parameters with the same weight
      const double match_quality
        = sqr(Params[iGaus*3 + 0] - pFunc->GetParameter(iFitGaus*3 + 0))
        + sqr(Params[iGaus*3 + 1] - pFunc->GetParameter(iFitGaus*3 + 1))
        + sqr(Params[iGaus*3 + 2] - pFunc->GetParameter(iFitGaus*3 + 2))
        ;
      if (match_quality < best_match_quality) {
        best_match_quality = match_quality;
        BestMatch[iGaus] = iFitGaus;
      } // if this is the best so far
    } // for fit
    
    assert(BestMatch[iGaus] < nGaus); // No reasonable solution matched!
    
    for (size_t i = 0; i < iGaus; ++i) {
      // Same fit Gaussian solution matches multiple solutions
      assert(BestMatch[i] != BestMatch[iGaus]);
    } // for i
    
  } // for solution
  
  std::vector<Double_t> sorted_params(nGaus * 3);
  for (size_t iGaus = 0; iGaus < nGaus; ++iGaus) {
    const size_t iBestFit = BestMatch[iGaus];
    for (size_t i = 0; i < 3; ++i)
      sorted_params[iGaus*3 + i] = pFunc->GetParameter(iBestFit*3 + i);
  }
  return sorted_params;
} // SortGaussianResults()


// Test a fit with a three-Gaussian function from the compiled cache
void ThreeGaussianFitTest(hit::GausFitCache& GausCache, float tol = 0.001) {
  std::string name = GausCache.GetName();
  
  constexpr size_t nGaus = 3; // we test three gaussians
  
  // sorted by increasing mean
  const Double_t Params[nGaus * 3] = { // amplitude, mean, sigma
    4.0, -2.0, 1.5,
    5.0,  0.1, 0.3,
    2.0,  1.0, 0.5
  }; // Params[]
  
  // fit test
  // - fill the input histogram
  TH1D Hist
    ("H3Gaus", ("Three-Gaussian test - " + name).c_str(), 200, -10., 10.);
  for (Int_t iBin = 1; iBin <= Hist.GetNbinsX(); ++iBin)
    Hist.SetBinContent(iBin, multi_gaus(Hist.GetBinCenter(iBin), nGaus, Params));
  
#ifdef GAUSFITCACHE_TEST_DEBUG
  Hist.Draw();
#endif // GAUSFITCACHE_TEST_DEBUG
  
  // - get the 3-Gaussian fitting function
  TF1* pFunc = GausCache.Get(nGaus);
  
  // - prepare the function with reasonable initial values
  const double MeanOffset = -((double(nGaus) - 1.) / 2.);
  const double MeanStep = 1.;
  for (size_t iGaus = 0; iGaus < nGaus; ++iGaus) {
    size_t BaseIndex = iGaus * 3;
    pFunc->SetParameter(BaseIndex + 0, 1.); // amplitude
    pFunc->SetParameter(BaseIndex + 1, MeanOffset + iGaus * MeanStep); // mean
    pFunc->SetParameter(BaseIndex + 2, 1.); // sigma
  } // for
  
  // function range is ignored in the fit, but used when drawing
  pFunc->SetRange(Hist.GetBinCenter(1), Hist.GetBinCenter(Hist.GetNbinsX()));
  
#ifdef GAUSFITCACHE_TEST_DEBUG
  TF1* pFCopy = pFunc->DrawCopy("L SAME");
  pFCopy->SetLineStyle(kDashed);
  pFCopy->SetLineColor(kCyan);
#endif // GAUSFITCACHE_TEST_DEBUG
  
  // - fit
  TFitResultPtr fit_res = Hist.Fit(pFunc, "WQ0NS");
  BOOST_CHECK(int(fit_res) == 0);
  
#ifdef GAUSFITCACHE_TEST_DEBUG
  pFCopy = pFunc->DrawCopy("L SAME");
  pFCopy->SetLineWidth(1);
  pFCopy->SetLineColor(kRed);
  
  gPad->SaveAs(("3GausTest_" + name + ".eps").c_str());
#endif // GAUSFITCACHE_TEST_DEBUG
  
  // - sort the results according to the expected ones
  std::vector<Double_t> results = SortGaussianResults(pFunc, Params);
  
  
  // - check them
  for (size_t iGaus = 0; iGaus < nGaus; ++iGaus) {
    const size_t iParam = iGaus * 3;
    
    // we use tolerance of 10^-3 (0.1%)
    // - check amplitude
    BOOST_CHECK_CLOSE(results[iParam + 0], Params[iParam + 0], tol * 100.);
    // - check mean
    BOOST_CHECK_CLOSE(results[iParam + 1], Params[iParam + 1], tol * 100.);
    // - check sigma
    BOOST_CHECK_CLOSE(results[iParam + 2], Params[iParam + 2], tol * 100.);
    
  } // for iGaus
  
  BOOST_WARN(pFunc->GetChisquare() < pFunc->GetNDF() / 2.);
  
#ifdef GAUSFITCACHE_TEST_DEBUG
  pFunc->Print();
  // always print the fit result
#else // !GAUSFITCACHE_TEST_DEBUG
  // in case of trouble, print the fit result
  if ((int(fit_res) != 0) || (pFunc->GetChisquare() >= pFunc->GetNDF() / 2.))
#endif // GAUSFITCACHE_TEST_DEBUG
  {
    fit_res->Print();
  }
  
} // ThreeGaussianFitTest()



// Test that the three-Gaussian function behaves as expected
BOOST_AUTO_TEST_CASE(ThreeGaussianTest)
{
  
  const Double_t Params[] = {
    4.0, -2.0, 1.5,
    5.0,  0.1, 0.3,
    2.0,  1.0, 0.5
  }; // Params[]
  
  // rounding errors on x on the look do not really matter
  for (double x = -10.; x <= +10.; x += 0.1) {
    const Double_t expected = multi_gaus(x, 3, Params);
    
    const Double_t tested
      = hit::details::CompiledGausFitCacheBaseStruct::ngaus<3>(&x, Params);
    
    // we use tolerance of 10^-5 (0.001%)
    BOOST_CHECK_CLOSE(tested, expected, 0.001);
  } // for x
  
} // BOOST_AUTO_TEST_CASE(ThreeGaussianTest)


// Test a fit with a three-Gaussian function from the run-time generated cache
BOOST_AUTO_TEST_CASE(RunTimeThreeGaussianFitTest)
{
  // max 20 Gaussians
  hit::GausFitCache GausCache("RunTimeGaussians");
  ThreeGaussianFitTest(GausCache, 2e-4); // 0.02% tolerance
} // BOOST_AUTO_TEST_CASE(RunTimeThreeGaussianFitTest)


// Test a fit with a three-Gaussian function from the compiled cache
BOOST_AUTO_TEST_CASE(CompiledThreeGaussianFitTest)
{
   // max 20 Gaussians
  hit::CompiledGausFitCache<20> GausCache("CompiledGaussians");
  ThreeGaussianFitTest(GausCache, 2e-4); // 0.02% tolerance
} // BOOST_AUTO_TEST_CASE(CompiledThreeGaussianFitTest)


// Test a fit with a three-Gaussian function (each truncated at 5 sigma)
// from the compiled cache
BOOST_AUTO_TEST_CASE(CompiledTruncated5ThreeGaussianFitTest)
{
   // max 20 Gaussians; truncate at 5 sigma
  hit::CompiledTruncatedGausFitCache<20, 5> GausCache
    ("CompiledTruncated5Gaussians");
  ThreeGaussianFitTest(GausCache, 1e-3); // 0.1% tolerance
} // BOOST_AUTO_TEST_CASE(CompiledTruncated5ThreeGaussianFitTest)


// Test a fit with a three-Gaussian function (each truncated at 4 sigma)
// from the compiled cache
BOOST_AUTO_TEST_CASE(CompiledTruncated4ThreeGaussianFitTest)
{
   // max 20 Gaussians; truncate at 4 sigma
  hit::CompiledTruncatedGausFitCache<20, 4> GausCache
    ("CompiledTruncated4Gaussians");
  ThreeGaussianFitTest(GausCache, 1e-3); // 0.1% tolerance
} // BOOST_AUTO_TEST_CASE(CompiledTruncated4ThreeGaussianFitTest)

// Test a fit with a three-Gaussian function (each truncated at 3 sigma)
// from the compiled cache
BOOST_AUTO_TEST_CASE(CompiledTruncated3ThreeGaussianFitTest)
{
   // max 20 Gaussians; truncate at 3 sigma
  hit::CompiledTruncatedGausFitCache<20, 3> GausCache
    ("CompiledTruncated3Gaussians");
  ThreeGaussianFitTest(GausCache, 1e-1); // 10% tolerance (seriously, it's that bad)
} // BOOST_AUTO_TEST_CASE(CompiledTruncated3ThreeGaussianFitTest)


BOOST_AUTO_TEST_SUITE_END()

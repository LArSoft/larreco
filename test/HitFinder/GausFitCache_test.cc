/**
 * @file   GausFitCache_test.cc
 * @brief  Test for classes in GausFitCache,h
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 27th, 2015
 * @see    GausFitCache,h
 */

// C/C++ standard libraries
#include <cmath>

// boost test libraries
#define BOOST_TEST_MODULE ( HitAnaAlg_test )
#include "cetlib/quiet_unit_test.hpp"

// LArSoft libraries
#include "larreco/HitFinder/GausFitCache.h"



double gaus(double x, double mean, double sigma, double amplitude) {
  double z = (x - mean) / sigma;
  return amplitude * std::exp(-0.5*z*z);
} // gaus()


// test that the test function behaves like a gaussian
BOOST_AUTO_TEST_CASE(GaussianTest)
{
  const exp1sigma = std::exp(-0.5);
  const exp2sigma = std::exp(-2.0);

  // here there should be some more believable test...
  const double amplitude = 16.;
  for (double sigma = 0.5; sigma < 2.25; ++sigma) {
    for (double mean = -5; mean < 5.5; ++mean) {

      const double x0 = mean / sigma;
      double value;

      // we use tolerance of 10^-5 (0.001%)
      BOOST_CHECK_CLOSE(gaus(-2., mean, sigma, amplitude), amplitude * exp2sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus(-1., mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus( 0., mean, sigma, amplitude), amplitude, 0.001);
      BOOST_CHECK_CLOSE(gaus( 1., mean, sigma, amplitude), amplitude * exp1sigma, 0.001);
      BOOST_CHECK_CLOSE(gaus( 2., mean, sigma, amplitude), amplitude * exp2sigma, 0.001);

    } // for mean
  } // for sigma

} // BOOST_AUTO_TEST_CASE(GaussianTest)


BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////
//  \file CalorimetryAlg.cxx
//
//  \brief Functions to calculate dE/dx. Based on code in Calorimetry.cxx
//
// andrzej.szelc@yale.edu
//
////////////////////////////////////////////////////////////////////////

// LArSoft includes
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larevt/CalibrationDBI/Interface/ElectronLifetimeProvider.h"
#include "larevt/CalibrationDBI/Interface/ElectronLifetimeService.h"

namespace calo {

  //--------------------------------------------------------------------
  CalorimetryAlg::CalorimetryAlg(const Config& config)
    : fCalAmpConstants{config.CalAmpConstants()}
    , fCalAreaConstants{config.CalAreaConstants()}
    , fUseModBox{config.CaloUseModBox()}
    , fLifeTimeForm{config.CaloLifeTimeForm()}
    , fDoLifeTimeCorrection{config.CaloDoLifeTimeCorrection()}
  {
    if (fLifeTimeForm != 0 and fLifeTimeForm != 1) {
      throw cet::exception("CalorimetryAlg")
        << "Unknow CaloLifeTimeForm " << fLifeTimeForm << '\n'
        << "Must select either '0' for exponential or '1' for exponential + "
           "constant.\n";
    }
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AMPLITUDE of the pulse
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                           detinfo::DetectorPropertiesData const& det_prop,
                           recob::Hit const& hit,
                           double const pitch,
                           double const T0) const
  {
    return dEdx_AMP(
      clock_data, det_prop, hit.PeakAmplitude() / pitch, hit.PeakTime(), hit.WireID().Plane, T0);
  }

  ///\todo The plane argument should really be for a view instead
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                           detinfo::DetectorPropertiesData const& det_prop,
                           double const dQ,
                           double const time,
                           double const pitch,
                           unsigned int const plane,
                           double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AMP(clock_data, det_prop, dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                           detinfo::DetectorPropertiesData const& det_prop,
                           double const dQdx,
                           double const time,
                           unsigned int const plane,
                           double const T0) const
  {
    double const fADCtoEl = fCalAmpConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0);
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AREA of the pulse
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            recob::Hit const& hit,
                            double const pitch,
                            double const T0) const
  {
    return dEdx_AREA(
      clock_data, det_prop, hit.Integral() / pitch, hit.PeakTime(), hit.WireID().Plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            double const dQ,
                            double const time,
                            double const pitch,
                            unsigned int const plane,
                            double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AREA(clock_data, det_prop, dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            double const dQdx,
                            double const time,
                            unsigned int const plane,
                            double const T0) const
  {
    double const fADCtoEl = fCalAreaConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0);
  }

  // Apply Lifetime and recombination correction.
  double
  CalorimetryAlg::dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clock_data,
                                   detinfo::DetectorPropertiesData const& det_prop,
                                   double dQdx_e,
                                   double const time,
                                   double const T0) const
  {
    if (fDoLifeTimeCorrection) {
      dQdx_e *= LifetimeCorrection(clock_data, det_prop, time, T0); // (dQdx_e in e/cm)
    }

    if (fUseModBox) { return det_prop.ModBoxCorrection(dQdx_e); }

    return det_prop.BirksCorrection(dQdx_e);
  }

  //------------------------------------------------------------------------------------//
  // for the time being copying from Calorimetry.cxx - should be decided where
  // to keep it.
  // ----------------------------------------------------------------------------------//
  double
  calo::CalorimetryAlg::LifetimeCorrection(detinfo::DetectorClocksData const& clock_data,
                                           detinfo::DetectorPropertiesData const& det_prop,
                                           double const time,
                                           double const T0) const
  {
    float const t = time - trigger_offset(clock_data);
    double const timetick = sampling_rate(clock_data) * 1.e-3; // time sample in microsec
    double const adjusted_time = t * timetick - T0 * 1e-3;     //  (in microsec)

    assert(fLifeTimeForm < 2);
    if (fLifeTimeForm == 0) {
      // Exponential form
      double const tau = det_prop.ElectronLifetime();
      return exp(adjusted_time / tau);
    }

    // Exponential+constant form
    auto const& elifetime_provider =
      art::ServiceHandle<lariov::ElectronLifetimeService const>()->GetProvider();
    return elifetime_provider.Lifetime(adjusted_time);
  }

} // namespace

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
    , fModBoxA{config.ModBoxA()}
    , fModBoxBF("ModBoxB", config.ModBoxBTF1().c_str(), 0, 90)
    , fBirksA{config.BirksA()}
    , fBirksKF("BirksK", config.BirksKTF1().c_str(), 0, 90)
  {
    if (fLifeTimeForm != 0 and fLifeTimeForm != 1) {
      throw cet::exception("CalorimetryAlg")
        << "Unknow CaloLifeTimeForm " << fLifeTimeForm << '\n'
        << "Must select either '0' for exponential or '1' for exponential + "
           "constant.\n";
    }

    // Set the parameters for the TF1s
    std::vector<double> modboxb_param;
    if (!config.ModBoxBParam(modboxb_param)) modboxb_param = {util::kModBoxB};

    for (unsigned i = 0; i < modboxb_param.size(); i++) {
      fModBoxBF.SetParameter(i, modboxb_param[i]);
    }

    std::vector<double> birksk_param;
    if (!config.BirksKParam(birksk_param)) birksk_param = {util::kRecombk};

    for (unsigned i = 0; i < birksk_param.size(); i++) {
      fBirksKF.SetParameter(i, birksk_param[i]);
    }
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AMPLITUDE of the pulse
  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                                  detinfo::DetectorPropertiesData const& det_prop,
                                  recob::Hit const& hit,
                                  double const pitch,
                                  double const T0) const
  {
    return dEdx_AMP(clock_data,
                    det_prop,
                    hit.PeakAmplitude() / pitch,
                    hit.PeakTime(),
                    hit.WireID().Plane,
                    T0,
                    det_prop.Efield());
  }

  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                                  detinfo::DetectorPropertiesData const& det_prop,
                                  double const dQdx,
                                  double const time,
                                  unsigned int const plane,
                                  double const T0) const
  {
    double const fADCtoEl = fCalAmpConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0, det_prop.Efield());
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AMPLITUDE of the pulse with EField
  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                                  detinfo::DetectorPropertiesData const& det_prop,
                                  recob::Hit const& hit,
                                  double const pitch,
                                  double const T0,
                                  double const EField,
                                  double const phi) const
  {
    return dEdx_AMP(clock_data,
                    det_prop,
                    hit.PeakAmplitude() / pitch,
                    hit.PeakTime(),
                    hit.WireID().Plane,
                    T0,
                    EField,
                    phi);
  }

  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                                  detinfo::DetectorPropertiesData const& det_prop,
                                  double const dQdx,
                                  double const time,
                                  unsigned int const plane,
                                  double const T0,
                                  double const EField,
                                  double const phi) const
  {
    double const fADCtoEl = fCalAmpConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0, EField, phi);
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AREA of the pulse
  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                                   detinfo::DetectorPropertiesData const& det_prop,
                                   recob::Hit const& hit,
                                   double const pitch,
                                   double const T0) const
  {
    return dEdx_AREA(clock_data,
                     det_prop,
                     hit.Integral() / pitch,
                     hit.PeakTime(),
                     hit.WireID().Plane,
                     T0,
                     det_prop.Efield());
  }

  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                                   detinfo::DetectorPropertiesData const& det_prop,
                                   double const dQdx,
                                   double const time,
                                   unsigned int const plane,
                                   double const T0) const
  {
    double const fADCtoEl = fCalAreaConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0, det_prop.Efield());
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AREA of the pulse with EField
  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                                   detinfo::DetectorPropertiesData const& det_prop,
                                   recob::Hit const& hit,
                                   double const pitch,
                                   double const T0,
                                   double const EField,
                                   double const phi) const
  {
    return dEdx_AREA(clock_data,
                     det_prop,
                     hit.Integral() / pitch,
                     hit.PeakTime(),
                     hit.WireID().Plane,
                     T0,
                     EField,
                     phi);
  }

  // ----------------------------------------------------------------------------------//
  double CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                                   detinfo::DetectorPropertiesData const& det_prop,
                                   double const dQdx,
                                   double const time,
                                   unsigned int const plane,
                                   double const T0,
                                   double const EField,
                                   double const phi) const
  {
    double const fADCtoEl = fCalAreaConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0, EField, phi);
  }

  // Apply Lifetime and recombination correction.
  double CalorimetryAlg::dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clock_data,
                                          detinfo::DetectorPropertiesData const& det_prop,
                                          double dQdx_e,
                                          double const time,
                                          double const T0) const
  {
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0, det_prop.Efield());
  }
  double CalorimetryAlg::dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clock_data,
                                          detinfo::DetectorPropertiesData const& det_prop,
                                          double dQdx_e,
                                          double const time,
                                          double const T0,
                                          double const EField,
                                          double const phi) const
  {
    if (fDoLifeTimeCorrection) {
      dQdx_e *= LifetimeCorrection(clock_data, det_prop, time, T0); // (dQdx_e in e/cm)
    }

    if (fUseModBox) { return ModBoxCorrection(dQdx_e, phi, det_prop.Density(), EField); }

    return BirksCorrection(dQdx_e, phi, det_prop.Density(), EField);
  }

  //------------------------------------------------------------------------------------//
  // for the time being copying from Calorimetry.cxx - should be decided where
  // to keep it.
  // ----------------------------------------------------------------------------------//
  double calo::CalorimetryAlg::LifetimeCorrection(detinfo::DetectorClocksData const& clock_data,
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

  // Reombination corrections

  // Modified box: allow for a general behavior of beta from phi, the angle between the track
  // and the electric field
  double calo::CalorimetryAlg::ModBoxCorrection(double dQdx,
                                                double phi,
                                                double rho,
                                                double E_field) const
  {
    // Modified Box model correction has better behavior than the Birks
    // correction at high values of dQ/dx.
    constexpr double Wion = 1000. / util::kGeVToElectrons; // 23.6 eV = 1e, Wion in MeV/e
    double const Beta = fModBoxBF.Eval(phi) / (rho * E_field);
    double const Alpha = fModBoxA;
    double const dEdx = (exp(Beta * Wion * dQdx) - Alpha) / Beta;

    return dEdx;
  }

  double calo::CalorimetryAlg::BirksCorrection(double dQdx,
                                               double phi,
                                               double rho,
                                               double E_field) const
  {
    // from: S.Amoruso et al., NIM A 523 (2004) 275

    double A = fBirksA;
    double K = fBirksKF.Eval(phi);                              // in KV/cm*(g/cm^2)/MeV
    constexpr double Wion = 1000. / util::kGeVToElectrons;      // 23.6 eV = 1e, Wion in MeV/e
    K /= rho;                                                   // KV/MeV
    double const dEdx = dQdx / (A / Wion - K / E_field * dQdx); // MeV/cm

    return dEdx;
  }

} // namespace

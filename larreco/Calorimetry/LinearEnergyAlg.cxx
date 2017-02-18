/**
 * @file   LinearEnergyAlg.cxx
 * @brief  Algorithm(s) calculating energy
 * @author Yun-Tse Tsai (yuntse@slac.stanford.edu) (adapting from the algorithm in LArLite)
 * @date   Feb 15, 2017
 * @see    LinearEnergyAlg.h
 *
 */

// LArSoft libraries
#include "larreco/Calorimetry/LinearEnergyAlg.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

// infrastructure and utilities
#include "cetlib/exception.h"

// C/C++ standard libraries
#include <stdexcept> // std::runtime_error()
#include <memory> // std::make_unique()

void calo::LinearEnergyAlg::initialize() {

  fElectronLifetime = detp->ElectronLifetime();
  fDeconNorm = 200;
  fRecombFactor = 1.;

}

void calo::LinearEnergyAlg::CalculateEnergy() {  // input clusters and hits, shower direction?

  for ( auto const &) {  //loop clusters
    for ( auto const& hit: hits ) {  // loop hits

      // lifetime correction
      double t = this->ConvertTickToDriftTime( hit->PeakTime(), hit->View() );
      double LifetimeCorr = std::exp( t / fElectronLifetime );

      // hit charge (ADC) -> Coulomb -> Number of electrons -> eV
      dE += hit->Integral() * kWion * fDeconNorm;
      // dE * lifetime correction
      dE *= LifetimeCorr;
      // dE / recombination factor
      double dEdx = 2.3;
      double RecombCorr = this->RecombinationCorrection( dEdx ) * kWion / dEdx;
      dE /= RecombCorr;
    }
  }
  return;
}

// Get the electron drift time; may move to somewhere else
double calo::LinearEnergyAlg::ConvertTickToDriftTime( double tick, geo::View_t plane ) const {
  double offTick = detp->GetXTicksOffset( plane, 0, 0 ) - dept->TriggerOffset();
  return detc->TPCTick2Time( tick - offTick );
}

double calo::LinearEnergyAlg::RecombinationCorrection( double dEdx ) {

  switch ( fRecombModel ) {
    case calo::kModBox:
      return this->ModBoxInverse( dEdx );
    case calo::kBirk:
      return this->BirksInverse( dEdx );
    case calo::kConstant:
      return fRecombFactor;
    default:
      return 1.;
  }
}

// Modified Box model correction, should be moved to somewhere else in the future
double calo::LinearEnergyAlg::ModBoxInverse( double dEdx ) const {
  // Modified Box model correction has better behavior than the Birks
  // correction at high values of dQ/dx.
  double rho     = detp->Density();                    // LAr density in g/cm^3
  double Efield  = detp->Efield();                     // Electric Field in the drift region in KV/cm
  double Beta    = util::kModBoxB / (rho * Efield);
  double Alpha   = util::kModBoxA;

  double dQdx = std::log ( Alpha + Beta * dEdx ) / ( Beta * kWion );

  return dQdx;
}

double calo::LinearEnergyAlg::BirksInverse( double dEdx ) const {
  // Correction for charge quenching using parameterization from
  // S.Amoruso et al., NIM A 523 (2004) 275

  double A3t     = util::kRecombA;
  double K3t     = util::kRecombk;                     // in KV/cm*(g/cm^2)/MeV
  double rho     = detp->Density();                    // LAr density in g/cm^3
  double Efield  = dept->Efield();                     // Electric Field in the drift region in KV/cm
  K3t           /= rho;                                // KV/MeV

  double dQdx = ( A3t/kWion ) / ( K3t / Efield * dEdx + 1);

  return dQdx;
}

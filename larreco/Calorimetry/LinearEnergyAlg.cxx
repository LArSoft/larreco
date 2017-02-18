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
#include "lardata/Utilities/ForEachAssociatedGroup.h" // util::associated_groups()
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"


// C/C++ standard libraries
#include <stdexcept> // std::runtime_error()
#include <cmath>


namespace {
  
  /**
   * @brief Returns the hits associated with the specified cluster.
   * @param cluster the cluster to retrieve the associated hits of
   * @param hitsPerCluster all the cluster-hit associations
   * @return a range of hits, that can be used in a range-for loop
   * @throw std::runtime_error if the cluster is not found
   */
  auto hitsAssociatedWith(
    art::Ptr<recob::Cluster> const& cluster,
    art::Assns<recob::Cluster, recob::Hit> const& hitsPerCluster
    )
  {
    //
    // Reminder: the cluster-hit association is a list of links:
    //     
    //     cluster #0  <--> hit #7
    //     cluster #0  <--> hit #8
    //     cluster #0  <--> hit #10
    //     cluster #1  <--> hit #4
    //     cluster #1  <--> hit #5
    //     cluster #1  <--> hit #9
    //     cluster #2  <--> hit #11
    //     cluster #2  <--> hit #12
    //     ...
    //     
    // associated_groups() presents this associations grouped by cluster:
    //     
    //     cluster #0  <--> hits { #7, #8, #10 }
    //     cluster #1  <--> hits { #4, #5, #9 }
    //     cluster #2  <--> hits { #11, #12 }
    //     ...
    //     
    // This function loops through this latter list, looking for the cluster
    // number that is stored in the `cluster` argument.
    //
    size_t iAssociatedCluster = 0;
    for (auto hits: util::associated_groups(hitsPerCluster)) {
      
      // if the iAssociatedCluster is the index we are looking for, we are done:
      if (iAssociatedCluster == cluster.key()) return hits;
      
      ++iAssociatedCluster; // else, try the next one!
      
    } // for all clusters
    throw std::runtime_error
      ("No hit associated with cluster " + std::to_string(cluster.key()) + " found!");
  } // calo::LinearEnergyAlg::hitsAssociatedWith()
  
} // local namespace



calo::LinearEnergyAlg::LinearEnergyAlg(Config const& config)
  : fUseArea( config.UseArea() )
  , fRecombFactor( 1. )
  , fElectronLifetime( 1e10 ) // needs to be read from service
  , fDeconNorm( 200 )
{
  if ( config.RecombinationModel() == "ModBox" ) fRecombModel = kModBox;
  else if ( config.RecombinationModel() == "Birk" ) fRecombModel = kBirk;
  else if ( config.RecombinationModel() == "Constant" ) fRecombModel = kConstant;
  else {
    throw std::runtime_error
      ( "Unsupported recombination mode: '" + config.RecombinationModel() + "'" );
  }
}


void calo::LinearEnergyAlg::initialize() {

  fElectronLifetime = detp->ElectronLifetime();

}



double calo::LinearEnergyAlg::CalculateHitEnergy(recob::Hit const& hit) const
{

  double const t = ConvertTickToDriftTime( hit.PeakTime(), hit.View() );
  double const LifetimeCorr = std::exp( t / fElectronLifetime );
  
  // hit charge (ADC) -> Coulomb -> Number of electrons -> eV
  double dE = hit.Integral() * kWion * fDeconNorm;
  
  // dE * lifetime correction
  dE *= LifetimeCorr;
  
  // dE / recombination factor
  double const dEdx = 2.3;
  double const RecombCorr = RecombinationCorrection( dEdx ) * kWion / dEdx;
  dE /= RecombCorr;
  
  return dE;
  
} // calo::LinearEnergyAlg::CalculateHitEnergy()



std::vector<double> calo::LinearEnergyAlg::CalculateEnergy(
  std::vector<art::Ptr<recob::Cluster>> const& clusters,
  art::Assns<recob::Cluster, recob::Hit> const& hitsPerCluster
  ) const
{  // input clusters and hits, shower direction?
  
  std::vector<double> clusterEnergies;
  
  for (art::Ptr<recob::Cluster> const& cluster: clusters) {
    
    // hitsAssociatedWith() searches for the right association links
    // to the cluster we are processing
    double const E = CalculateClusterEnergy
      (*cluster, hitsAssociatedWith(cluster, hitsPerCluster));
    
    clusterEnergies.push_back(E);
    
  } // for all clusters
  
  return clusterEnergies;
  
} // calo::LinearEnergyAlg::CalculateEnergy()

// Get the electron drift time; may move to somewhere else
double calo::LinearEnergyAlg::ConvertTickToDriftTime( double tick, geo::View_t plane ) const {
  double offTick = detp->GetXTicksOffset( plane, 0, 0 ) - detp->TriggerOffset();
  return detc->TPCTick2Time( tick - offTick );
}

double calo::LinearEnergyAlg::RecombinationCorrection( double dEdx ) const {

  switch ( fRecombModel ) {
    case kModBox:
      return this->ModBoxInverse( dEdx );
    case kBirk:
      return this->BirksInverse( dEdx );
    case kConstant:
      return fRecombFactor;
    default:
      throw std::logic_error
        ("Unexpected! recombination model in RecombinationCorrection()");
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
  double Efield  = detp->Efield();                     // Electric Field in the drift region in KV/cm
  K3t           /= rho;                                // KV/MeV

  double dQdx = ( A3t/kWion ) / ( K3t / Efield * dEdx + 1);

  return dQdx;
}

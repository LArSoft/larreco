/**
 * @file   LinearEnergyAlg.h
 * @brief  Algorithm(s) calculating energy
 * @author Yun-Tse Tsai (yuntse@slac.stanford.edu) (adapting from the algorithm in LArLite)
 * @date   Feb 15, 2017
 *
 */

#ifndef LARRECO_CALORIMETRY_LINEARENERGYALG_H
#define LARRECO_CALORIMETRY_LINEARENERGYALG_H


// LArSoft libraries
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// infrastructure and utilities
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <vector>
#include <memory> // std::unique_ptr<>

namespace calo {

  class LinearEnergyAlg {

    public:

    enum RecombinationModel_t {
      kModBox,
      kBirk,
      kConstant,
      nRecombinationModel
    }

    /// Algorithm configuration
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<bool> UseArea {
        Name("UseArea"),
        Comment("Whether to use the area of hit to count the deposited charges.")
      };

      fhicl::Atom<std::string> RecombinationModel {
        Name("RecombinationModel"),
        Comment("Which recombination model to use: ModBox, Birk, Constant.")
      };

    }; // Config

    /// @{
    /// @name Construction and configuration

    /**
     * @brief Constructor with configuration validation
     * @param config configuration parameter structure
     *
     * For the configuration, see `LinearEnergyAlg` documentation.
     */
    LinearEnergyAlg(Config const& config)
      : fUseArea( config.UseArea() ),
        fRecombModel( config.RecombinationModel() ) {
          if ( config.RecombinationModel() == "ModBox" ) fRecombModel = kModBox;
          else if ( config.RecombinationModel() == "Birk" ) fRecombModel = kBirk;
          else if ( config.RecombinationModel() == "Constant" ) fRecombModel = kConstant;
          else {
            throw std::runtime_error
            ( "Unsupported recombination mode: '" + config.RecombinationModel() + "'" );
          } {}

    /**
     * @brief Constructor with configuration validation
     * @param pset FHiCL configuration parameter set
     * @see SpacePointIsolationAlg(Config const&)
     *
     * Translates the parameter set into a configuration object and uses the
     * validating constructor to initialise the object.
     *
     * For the configuration, see `LinearEnergyAlg` documentation.
     */
    LinearEnergyAlg(fhicl::ParameterSet const& pset)
      : LinearEnergyAlg(fhicl::Table<Config>(pset, {})())
      {}

    /// @}

    /// @{
    /// @name Set up

    /**
     * @brief Sets up the algorithm
     * @param geometry the geometry service provider
     *
     * Acquires the geometry description.
     * This method must be called every time the geometry is changed.
     */
    void setup( detinfo::DetectorProperties const& detproperty, detinfo::DetectorClocks const& detclock, geo::GeometryCore const& geometry )
      { detp = &detproperty; detc = &detclock; geom = &geometry; initialize(); }

    /// @}

  private:

    /// Pointer to the geometry to be used
    geo::GeometryCore const* geom = nullptr;

    /// Pointer to the detector property
    detinfo::DetectorProperties const* detp = nullptr;

    /// Pointer to the detector clock
    detinfo::DetectorClocks const* detc = nullptr;

    bool   fUseArea;
    int    fRecombModel;
    // TODO
    // double fRecombA;
    // double fRecombk;
    // double fModBoxA;
    // double fModBoxB;
    double fRecombFactor;
    double fElectronLifetime;
    double fDeconNorm;

    static const double kWion           = 23.6e-6;  ///< ionization potenial in LAr, 23.6 eV = 1e, Wion in MeV/e
    // Conversion for energy deposited in GeV to number of ionization electrons produced
    static const double kRecombFactor   = 0.62;     ///< constant correction used in the current MicroBooNE shower reconstruction

    void initialize();

    void CalculateEnergy();

    double ConvertTickToDriftTime( double tick, geo::View_t plane ) const;

    double RecombinationCorrection( double dEdx );  ///< TODO: make it more flexible

    double ModBoxInverse( double dEdx ) const;

    double BirksInverse( double dEdx ) const;

  }  // class LinearEnergyAlg

}  // namespace calo
#endif  // LARRECO_CALORIMETRY_LINEARENERGYALG_H

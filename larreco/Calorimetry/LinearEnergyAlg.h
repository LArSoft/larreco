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
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::kModBoxA ...

// infrastructure and utilities
// #include "cetlib/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <vector>
#include <utility> // std::forward()
#include <type_traits> // std::is_same, std::decay_t

namespace calo {
  
  /**
   * @brief Calibrates the energy of the clusters.
   * 
   * Configuration
   * --------------
   * 
   * Example of configuration:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * EnergyAlgo: {
   *   UseArea: true
   *   Recombination: {
   *     Model: Birks
   *     A3t: 0.8
   *     k3t: 0.0486
   *   }
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * Parameters:
   * * *UseArea* (flag, _mandatory): whether to use hit integral
   *     (`recob::Hit::Integral()`) instead of hit peak amplitude for charge
   *     estimation
   * * *Recombination* (mandatory table): configures the recombination model
   *     to be used and its parameters:
   *     * *Model* (string, _mandatory_): specifies one of the
   *         supported recombination models: `Birks`, `ModBox` or `Constant`.
   *         Each requires its own specific configuration parameters.
   *         * `Birks`, Birks Law model
   *            (see S.Amoruso et al., NIM A 523 (2004) 275) requires
   *            configuration parameters:
   *             * *A3t* (real, default: `util::kRecombA`)
   *             * *k3t* (real, default: `util::kRecombk`) in kV/cm*(g/cm^2)/MeV
   *         * `ModBox`, box model requires configuration parameters:
   *             * *A* (real, default: `util::kModBoxA`)
   *             * *B* (real, default: `util::kModBoxB`) in kV/cm*(g/cm^2)/MeV
   *         * `Constant` recombination factor, requires configuration
   *             parameters:
   *             * *factor* (real, _mandatory_): the recombination factor,
   *                 uniformly applied to all hits
   * 
   */
  class LinearEnergyAlg {

    public:
    
    struct ModelName {
      static const std::string ModBox;
      static const std::string Birks;
      static const std::string Constant;
    };
    
    enum RecombinationModel_t {
      kModBox,
      kBirks,
      kConstant,
      nRecombinationModel
    };
    
    
    /// Configuration of parameters of the box model.
    struct RecombinationConfig {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      bool modelIsBirks() const
        { return Model() == ModelName::Birks; }
      bool modelIsModBox() const
        { return Model() == ModelName::ModBox; }
      bool modelIsConstant() const
        { return Model() == ModelName::Constant; }
      
      
      fhicl::Atom<std::string> Model {
        Name("Model"),
        Comment(std::string("Which recombination model to use: "
          + ModelName::ModBox + ", "
          + ModelName::Birks + ", "
          + ModelName::Constant + ".").c_str()
          )
        };
      
      fhicl::Atom<double> A {
        Name("A"),
        Comment("Parameter \"A\" of box model."),
        fhicl::use_if(this, &RecombinationConfig::modelIsModBox),
        util::kModBoxA
        };

      fhicl::Atom<double> B {
        Name("B"),
        Comment("Parameter \"B\" of box model [kV/cm*(g/cm^2)/MeV]."),
        fhicl::use_if(this, &RecombinationConfig::modelIsModBox),
        util::kModBoxB
        };
      
      fhicl::Atom<double> A3t {
        Name("A3t"),
        Comment("Recombination parameter \"A\" of Birks model."),
        fhicl::use_if(this, &RecombinationConfig::modelIsBirks),
        util::kRecombA
        };

      fhicl::Atom<double> k3t {
        Name("k3t"),
        Comment("Recombination parameter \"k\" of Birks model [kV/cm*(g/cm^2)/MeV]."),
        fhicl::use_if(this, &RecombinationConfig::modelIsBirks),
        util::kRecombk
        };
      
      fhicl::Atom<double> factor {
        Name("factor"),
        Comment("Constant recombination factor for \"constant\" model."),
        fhicl::use_if(this, &RecombinationConfig::modelIsConstant)
        };
      
    }; // RecombinationConfig
    
    
    /// Algorithm configuration
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<bool> UseArea {
        Name("UseArea"),
        Comment("Whether to use the area of hit to count the deposited charges.")
      };

      fhicl::Table<RecombinationConfig> Recombination {
        Name("Recombination"),
        Comment("Parameters of the recombination model")
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
    LinearEnergyAlg(Config const& config);
    
    
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
    
    /**
     * @brief Prints the current algorithm configuration.
     * @tparam Stream type of output stream 
     * @param out output stream where to write the information
     * @param indent indentation for all lines
     * @param firstIndent special indentation for the first line
     * 
     * The configuration parameters are printed in human-readable form.
     * The output starts at the current line of the stream, and it terminates
     * with a new line.
     * 
     * This method does not require the algorithm to have been set up (via
     * `Setup()`), since it just prints user configuration as given at
     * construction time.
     */
    template <typename Stream>
    void DumpConfiguration
      (Stream&& out, std::string indent, std::string firstIndent) const;
    
    /**
     * @brief Prints the current algorithm configuration.
     * @tparam Stream type of output stream 
     * @param out output stream where to write the information
     * @param indent (default: none) indentation for all lines (including the
     *               first one)
     * 
     * The configuration parameters are printed in human-readable form.
     * The output starts at the current line of the stream, and it terminates
     * with a new line.
     */
    template <typename Stream>
    void DumpConfiguration(Stream&& out, std::string indent = "") const
      { DumpConfiguration(std::forward<Stream>(out), indent, indent); }

    /// @}
    
    /// @{
    /// @name Operations
    
    /**
     * @brief Calculates the energy of a single cluster.
     * @tparam Hits a range of hits associated with the cluster
     * @param cluster list of clusters we want the energy of
     * @param hits pointers to all hits associated to the cluster
     * @return the calibrated energy of the specified cluster [GeV]
     * 
     * @todo Describe the algorithm
     * 
     * The `hits` are stored as (constant) pointers into a "range", that is any
     * object (e.g. a STL vector) supporting `begin()` and `end()` iterators.
     * 
     */
    template <typename Hits>
    double CalculateClusterEnergy
      (recob::Cluster const& cluster, Hits&& hits) const;
    
    /**
     * @brief Calculates the energy of the shower
     * @param clusters list of clusters we want the energy of
     * @param hitsPerCluster associations of all clusters to all hits
     * @return a vector with one energy per input cluster [GeV]
     * @see CalculateClusterEnergy()
     * 
     * A vector is returned with a entry, that is the cluster energy, for each
     * of the clusters in the input, in the same order.
     * See `CalculateClusterEnergy()` for the algorithm details.
     */
    std::vector<double> CalculateEnergy(
      std::vector<art::Ptr<recob::Cluster>> const& clusters,
      art::Assns<recob::Cluster, recob::Hit> const& hitsPerCluster
      ) const;
    
    /// @}
    
  private:
    
    struct ModBoxParameters {
      double A = util::kModBoxA;
      double B = util::kModBoxB;
    };
    
    struct BirksParameters {
      double A = util::kRecombA;
      double k = util::kRecombk;
    };
    
    struct ConstantRecombParameters {
      double factor = 1.0; // no sensible default value here
    };
    
    /// Pointer to the geometry to be used
    geo::GeometryCore const* geom = nullptr;

    /// Pointer to the detector property
    detinfo::DetectorProperties const* detp = nullptr;

    /// Pointer to the detector clock
    detinfo::DetectorClocks const* detc = nullptr;

    bool   fUseArea;
    int    fRecombModel;
    
    /// Parameters for recombination box model; filled only when this model is selected.
    ModBoxParameters recombModBoxParams;
    /// Parameters for recombination Birks model; filled only when this model is selected.
    BirksParameters recombBirksParams;
    /// Parameters for constant recombination factor; filled only when this model is selected.
    ConstantRecombParameters recombConstParams;
    
    // TODO
    // double fRecombA;
    // double fRecombk;
    // double fModBoxA;
    // double fModBoxB;
    double fRecombFactor;
    double fElectronLifetime;
    double fDeconNorm;

    static constexpr double kWion           = 23.6e-6;  ///< ionization potenial in LAr, 23.6 eV = 1e, Wion in MeV/e
    // Conversion for energy deposited in GeV to number of ionization electrons produced
    static constexpr double kRecombFactor   = 0.62;     ///< constant correction used in the current MicroBooNE shower reconstruction

    void initialize();
    
    /**
     * @brief Returns the corrected energy from the hit.
     * @param hit the hit to be corrected
     * @return the corrected hit energy [GeV]
     * 
     * @todo document the algorithm here
     * 
     */
    double CalculateHitEnergy(recob::Hit const& hit) const;

    double ConvertTickToDriftTime( double tick, geo::View_t plane ) const;

    double RecombinationCorrection( double dEdx ) const;  ///< TODO: make it more flexible

    double ModBoxInverse( double dEdx ) const;

    double BirksInverse( double dEdx ) const;

  };  // class LinearEnergyAlg

}  // namespace calo



//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Hits>
double calo::LinearEnergyAlg::CalculateClusterEnergy
  (recob::Cluster const& /* cluster */, Hits&& hits) const
{
  using std::begin;
  static_assert( // check the hits type
    std::is_same<std::decay_t<decltype(**begin(hits))>, recob::Hit>::value,
    "Hits must be a collection of recob::Hit pointers!" // any pointer will do
    );
  
  double E = 0.0; // total cluster energy
  
  for (auto hitPtr: hits) {
    
    E += CalculateHitEnergy(*hitPtr);
    
  } // for
    
  return E;
  
} // calo::LinearEnergyAlg::CalculateEnergy()

//------------------------------------------------------------------------------
template <typename Stream>
void calo::LinearEnergyAlg::DumpConfiguration
  (Stream&& out, std::string indent, std::string firstIndent) const
{
  
  out << firstIndent << "LinearEnergyAlg configuration:"
    << "\n" << indent << "  use hit " << (fUseArea? "area": "peak amplitude")
      << " for charge estimation"
    << "\n" << indent << "  recombination model: ";
  switch ( fRecombModel ) {
    case kModBox:
      out << ModelName::ModBox
        << "\n" << indent << "    A = " << recombModBoxParams.A
        << "\n" << indent << "    B = " << recombModBoxParams.B;
      break;
    case kBirks:
      out << ModelName::Birks
        << "\n" << indent << "    A = " << recombBirksParams.A
        << "\n" << indent << "    k = " << recombBirksParams.k;
      break;
    case kConstant:
      out << ModelName::Constant
        << "\n" << indent << "    k = " << recombConstParams.factor;
      break;
    default:
      out << "invalid (" << ((int) fRecombModel) << ")!!!";
  } // switch
  out << "\n";
  
} // calo::LinearEnergyAlg::DumpConfiguration()


//------------------------------------------------------------------------------





#endif  // LARRECO_CALORIMETRY_LINEARENERGYALG_H

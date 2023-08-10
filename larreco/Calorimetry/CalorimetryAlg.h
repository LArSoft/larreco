////////////////////////////////////////////////////////////////////////
// \file CalorimetryAlg.h
//
// \brief Functions to calculate dE/dx. Based on code in Calorimetry.cxx
//
// \author andrzej.szelc@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_CALORIMETRYALG_H
#define UTIL_CALORIMETRYALG_H

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TF1.h"

#include <vector>

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace recob {
  class Hit;
}

/// General LArSoft Utilities
namespace calo {
  class CalorimetryAlg {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Sequence<double> CalAmpConstants{
        Name("CalAmpConstants"),
        Comment("ADC to electrons constants for each plane.")};

      fhicl::Sequence<double> CalAreaConstants{
        Name("CalAreaConstants"),
        Comment("Area to electrons constants for each plane.")};

      fhicl::Atom<bool> CaloUseModBox{Name("CaloUseModBox"),
                                      Comment("Use modified box model if true, birks otherwise")};

      fhicl::Atom<int> CaloLifeTimeForm{Name("CaloLifeTimeForm"),
                                        Comment("0 = exponential, 1 = exponential + constant")};

      fhicl::Atom<bool> CaloDoLifeTimeCorrection{Name("CaloDoLifeTimeCorrection"),
                                                 Comment("Apply lifetime correction if true")};

      fhicl::Atom<double> ModBoxA{Name("ModBoxA"),
                                  Comment("Alpha value in modified box recombination."),
                                  util::kModBoxA};

      fhicl::Atom<std::string> ModBoxBTF1{
        Name("ModBoxBTF1"),
        Comment(
          "String compiled into a TF1. Should return the Mod-Box beta value as a function of phi."),
        "[0]"};

      fhicl::OptionalSequence<double> ModBoxBParam{
        Name("ModBoxBParam"),
        Comment("Parameters for the ModBoxBTF1 function.")};

      fhicl::Atom<double> BirksA{Name("BirksA"),
                                 Comment("Alpha value in modified box recombination."),
                                 util::kRecombA};

      fhicl::Atom<std::string> BirksKTF1{
        Name("BirksKTF1"),
        Comment(
          "String compiled into a TF1. Should return the Birks k value as a function of phi."),
        "[0]"};

      fhicl::OptionalSequence<double> BirksKParam{
        Name("BirksKParam"),
        Comment("Parameters for the BirksKTF1 function. List of doubles.")};
    };

    CalorimetryAlg(const fhicl::ParameterSet& pset)
      : CalorimetryAlg(fhicl::Table<Config>(pset, {})())
    {}

    CalorimetryAlg(const Config& config);

    double dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                    detinfo::DetectorPropertiesData const& det_prop,
                    recob::Hit const& hit,
                    double pitch,
                    double T0 = 0) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                    detinfo::DetectorPropertiesData const& det_prop,
                    double dQdx,
                    double time,
                    unsigned int plane,
                    double T0 = 0) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                    detinfo::DetectorPropertiesData const& det_prop,
                    recob::Hit const& hit,
                    double pitch,
                    double T0,
                    double EField,
                    double phi = 90) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                    detinfo::DetectorPropertiesData const& det_prop,
                    double dQdx,
                    double time,
                    unsigned int plane,
                    double T0,
                    double EField,
                    double phi = 90) const;

    // FIXME: How may of these are actually used?
    double dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     recob::Hit const& hit,
                     double pitch,
                     double T0 = 0) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     double dQdx,
                     double time,
                     unsigned int plane,
                     double T0 = 0) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     recob::Hit const& hit,
                     double pitch,
                     double T0,
                     double EField,
                     double phi = 90) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     double dQdx,
                     double time,
                     unsigned int plane,
                     double T0,
                     double EField,
                     double phi = 90) const;

    double ElectronsFromADCPeak(double adc, unsigned short plane) const
    {
      return adc / fCalAmpConstants[plane];
    }

    double ElectronsFromADCArea(double area, unsigned short plane) const
    {
      return area / fCalAreaConstants[plane];
    }

    double LifetimeCorrection(detinfo::DetectorClocksData const& clock_data,
                              detinfo::DetectorPropertiesData const& det_prop,
                              double time,
                              double T0 = 0) const;

    // Recombination corrections
    double BirksCorrection(double dQdx, double phi, double rho, double E_field) const;
    double ModBoxCorrection(double dQdx, double phi, double rho, double E_field) const;

  private:
    art::ServiceHandle<geo::Geometry const> geom;

    double dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            double dQdx_e,
                            double time,
                            double T0 = 0) const;
    double dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            double dQdx_e,
                            double time,
                            double T0,
                            double EField,
                            double phi = 90) const;

    std::vector<double> const fCalAmpConstants;
    std::vector<double> const fCalAreaConstants;
    bool const fUseModBox;
    int const fLifeTimeForm;
    bool const fDoLifeTimeCorrection;

    // Recombination parameters
    double fModBoxA; // Mod-Box alpha
    TF1 fModBoxBF;   // Function of phi to get the Mod-Box beta value
    double fBirksA;  // Birks A cosntant
    TF1 fBirksKF;    // Function of phi to get the Birks-k value

  }; // class CalorimetryAlg
} // namespace calo
#endif // UTIL_CALORIMETRYALG_H

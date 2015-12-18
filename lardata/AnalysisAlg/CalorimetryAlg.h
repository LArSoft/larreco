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

#include "fhiclcpp/ParameterSet.h"

#include "DetectorInfoServices/LArPropertiesService.h"
#include "Geometry/Geometry.h"
#include "DetectorInfoServices/DetectorPropertiesService.h"
#include <vector>

namespace recob { 
  class Hit; 
}


///General LArSoft Utilities
namespace calo{
    class CalorimetryAlg {
    public:


    CalorimetryAlg(fhicl::ParameterSet const& pset);
    
    ~CalorimetryAlg();
      
    void   reconfigure(fhicl::ParameterSet const& pset);
    
    double dEdx_AMP(art::Ptr< recob::Hit >  hit, double pitch, double T0=0) const;
    double dEdx_AMP(recob::Hit const&  hit, double pitch, double T0=0) const;
    double dEdx_AMP(double dQ, double time, double pitch, unsigned int plane, double T0=0) const;
    double dEdx_AMP(double dQdx,double time, unsigned int plane, double T0=0) const;
    
    double dEdx_AREA(art::Ptr< recob::Hit >  hit, double pitch, double T0=0) const;
    double dEdx_AREA(recob::Hit const&  hit, double pitch, double T0=0) const;
    double dEdx_AREA(double dQ,double time, double pitch, unsigned int plane, double T0=0) const;
    double dEdx_AREA(double dQdx,double time, unsigned int plane, double T0=0) const;
      
    double ElectronsFromADCPeak(double adc, unsigned short plane) const
    { return adc / fCalAmpConstants[plane]; }
      
    double ElectronsFromADCArea(double area, unsigned short plane) const
    { return area / fCalAreaConstants[plane]; }
    
    double LifetimeCorrection(double time, double T0=0) const;
    
  private:

    art::ServiceHandle<geo::Geometry> geom; 
    const detinfo::DetectorProperties* detprop;

    double dEdx_from_dQdx_e(double dQdx_e,double time, double T0=0) const;
   
    std::vector< double > fCalAmpConstants;
    std::vector< double > fCalAreaConstants;
    bool fUseModBox;
    
    
    }; // class CalorimetryAlg
} //namespace calo
#endif // UTIL_CALORIMETRYALG_H
